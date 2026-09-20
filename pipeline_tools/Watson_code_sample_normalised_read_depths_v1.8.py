#!/usr/bin/env python

'''''
Watson code for calculating and plotting normalised read depths across sample.
Version 1.8 (May 2021)

Input:
    1) mapped SSCS bam (overlapping reads clipped)
    2) sample name
    3) panel bed file
    4) chromosome ideogram

Outputs:
    Metrics files:
        1) metrics txt file containing mean read depths and normalised read depths for each probe region on the provided panel
    Plots:
        1) plot of normalised read depths across each individual chromosome
        2) plot of normalised read depths across each chromosome (all chromosomes on one plot)

Usage:
Watson_code_sample_normalised_read_depths_v1.8.py --infile mapped SSCS BAM --sample-name SAMPLENAME --panel_bed  panel bed file
                                --chromosome_ideogram txt file of band positions --out_directory output_file_directory

'''''
version = '1.8'

from argparse import ArgumentParser
import pysam
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MaxNLocator, MultipleLocator
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import timeit
import time
from datetime import date
from Bio.Seq import Seq
import csv

today = date.today()
date_today = today.strftime("%d/%m/%Y")

# Lists of colors for plots
c0 = (0.76, 0.76, 0.76)
c1 = (1.00, 0.18, 0.33);
c2 = (1.00, 0.23, 0.19);
c3 = (1.00, 0.58, 0.00);
c4 = (1.00, 0.80, 0.00);
c5 = (0.30, 0.85, 0.39);
c6 = (0.35, 0.78, 0.98);
c7 = (0.20, 0.67, 0.86);
c8 = (0.00, 0.48, 1.00);
c9 = (0.35, 0.34, 0.84);
c10 = (0.00, 0.31, 0.57);
c11 = (0.12, 0.29, 0.69);
c12 = (0.17, 0.17, 0.42);
c13 = (1.00, 1.00, 1.00);
c14 = (0.77, 0.04, 0.00);

grey1 = '#f7f7f7'
grey2 = '#cccccc'
grey3 = '#969696'
grey4 = '#636363'
grey5 = '#252525'


def ideograms(ideogram_file, chromosome): #function for plotting the chromosomal bands

    color_lookup = {'gneg': (1., 1., 1.),
                    'gpos25': (.6, .6, .6),
                    'gpos50': (.4, .4, .4),
                    'gpos75': (.2, .2, .2),
                   'gpos100': (0., 0., 0.),
                      'acen': (.8, .4, .4),
                      'gvar': (.8, .8, .8),
                     'stalk': (.9, .9, .9)}

    ideogram = open(ideogram_file)
    ideogram.readline()
    xranges = []
    colors = []
    mid_points = []
    labels = []

    for line in ideogram:
        chrom, start, stop, label, stain = line.strip().split('\t')
        start = int(start)
        stop = int(stop)
        width = stop - start
        mid_point = start + (width/2)
        if chrom == chromosome:
            xranges.append((start, width))
            colors.append(color_lookup[stain])
            mid_points.append(mid_point)
            labels.append(label)

    return xranges, [0, 0.9], colors, mid_points, labels

def plot_chromosome(chromosome_ideogram, chromosome, ax): #function for plotting the chromosomal bands on the chromosome

    xranges, yrange, colors, midpoints, labels = ideograms(chromosome_ideogram, chromosome)

    ax.broken_barh(xranges, yrange, facecolors= colors, edgecolor = 'black')

    ax.set_xticks(midpoints)
    ax.set_xticklabels(labels, rotation = 90, fontsize = 9)
    ax.set_yticks([])
    ax.text(-0.013, 0.35, chromosome, transform=ax.transAxes, fontsize = 15, ha = 'right')
    ax.xaxis.set_tick_params(width=0.8, color = grey3, length = 6)

    ax.minorticks_off()

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    return ax

def chromosome_sizes_bands(chromosome_ideogram):
    with open(chromosome_ideogram, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        chromosome_sizes = {}
        chromosome_bands = {}
        for row in read_reader:
            if row_count > 0:
                chromosome = row[0]
                start = int(row[1]) #start position of that part of the chromosome
                end = int(row[2]) #end position of that part of the chromosome
                band = row[3] #e.g. 'q33.1'
                arm = band[0] #e.g. 'q'
                part = row[4] #e.g. 'acen' if centromeric region of chromosome

                chromosome_sizes[chromosome]=end #the dictionary will keep getting overwritten until it gets to the end of the chromosome, so the last end position it adds will be the total length of the chromosome

                if chromosome in chromosome_bands.keys():
                    chromosome_bands[chromosome].append((start, end, band))
                else:
                    chromosome_bands[chromosome]=[(start, end, band)]

            row_count+=1

    return chromosome_sizes, chromosome_bands

def panel_positions_from_bed(panel_bed):
    with open(panel_bed, 'r') as bedfile:
        read_reader = csv.reader(bedfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        panel_positions = []

        for row in read_reader:
            if row_count > 2:
                chromosome = row[0] #e.g. 'chr2'
                start = row[1]
                stop = row[2]
                panel_positions.append((chromosome, int(start), int(stop)))
            row_count+=1

    return panel_positions

def panel_coverage_from_bed(panel_bed):
    with open(panel_bed, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        panel_coverage = {}
        for row in read_reader:
            if row_count > 2:
                chromosome = row[0] #e.g. chr2
                start = int(row[1]) #start position of that part of the chromosome
                end = int(row[2]) #end position of that part of the chromosome

                if chromosome in panel_coverage.keys():
                    panel_coverage[chromosome].append((start, end))
                else:
                    panel_coverage[chromosome]=([(start, end)])

            row_count+=1

    return panel_coverage

def normalised_read_depths(SSCS_bam, panel_positions, chromosome_bands, out_directory, sample_name):
    samfile = pysam.AlignmentFile(SSCS_bam, "rb")

    region_depths = {}
    all_positions_depths = 0
    total_positions = 0

    counter = 0
    for region in panel_positions: #iterate over each of the probe regions in the panel...
        start_time = time.time()
        chromosome = region[0].replace('chr', '') #the contigs in the samfile with the GATK reference genome are '1' rather than 'chr1'
        start = region[1]
        stop = region[2]

        position_depths = [] #create an empty list in which the depths at each of the positions in the region will be stored.
        for pileupcolumn in samfile.pileup(chromosome, start, stop, max_depth = 100000, stepper = 'all'): #create a pileup of all the reads that cover the probe
            if pileupcolumn.pos in range(start, stop+1): #just look at the pileups of the position in the probe regions (i.e. not the starts of regions that start outside of the probe region)
                depth = pileupcolumn.n #depth at that position
                position_depths.append(depth) #at the depths to the list for that probe region
                all_positions_depths+=depth #at the depth to the running list of all depths across the panel
                total_positions +=1 #add 1 to the total number of positions in the panel

        mean_region_depth = np.mean(position_depths) #calculate the mean depth across the probe region
        region_depths[(chromosome, start, stop)] = mean_region_depth #add the mean depth to the dictionary where the key is the probe region

        counter+=1
        if counter % 100 == 0:
            print('total probe regions processed: ', counter)
            print('time for last 100 probe regions to be writen = %s seconds' % int(time.time() - start_time))

    samfile.close()

    mean_depth_across_panel = all_positions_depths/total_positions
    metrics_dictionary = {}
    normalised_region_depths = {}
    for k, v in region_depths.items():
        chromosome = 'chr'+k[0] #e.g. chr2
        start = k[1]
        stop = k[2]
        normalised_depth = v/mean_depth_across_panel
        normalised_region_depths[(chromosome, start, stop)]=normalised_depth
        bands_on_chromosome = chromosome_bands[chromosome]
        bands = []
        for band in bands_on_chromosome:
            band_start = band[0]
            band_stop = band[1]
            band_name = band[2]
            if start < band_stop and stop > band_start:
                bands.append(band_name)
        metrics_dictionary[k]=(v, normalised_depth, bands[0])

    read_depth_metrics = open(out_directory+'/CNV_read_depths/'+sample_name+'_sample_normalised_read_depths.txt', 'w')
    read_depth_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    read_depth_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    read_depth_metrics.write('produced from code:\t' + 'Watson_code_sample_normalised_read_depths: version '+str(version) + '\n')
    read_depth_metrics.write('\n')
    read_depth_metrics.write('mean depth across panel:\t' + str(mean_depth_across_panel))
    read_depth_metrics.write('\n')
    read_depth_metrics.write('chromosome\t' + 'start\t' + 'stop\t' + 'band\t' + 'mean_region_depth\t' + 'normalised_region_depth\n')
    for k, v in metrics_dictionary.items():
        read_depth_metrics.write(k[0]+'\t' + str(k[1])+'\t'+ str(k[2])+'\t' + str(v[2])+'\t' + str(v[0])+'\t'+ str(v[1]) + '\n')
    read_depth_metrics.close()

    return normalised_region_depths

def mean_depth_across_bands(chromosome, chromosome_bands, normalised_region_depths, earliest_position, latest_position):
    chrom = 'chr'+chromosome
    band_depths = {}
    for a, b in chromosome_bands.items():
        if a == 'chr'+chromosome: #e.g. '2'
            for band in b: #iterate over the tuples of band start and stop positions in that chromosome
                band_start = int(band[0])
                band_stop = int(band[1])
                for k, v in normalised_region_depths.items(): #for each band, look to see what the mean coverage of the probes that covered that band were
                    if k[0] == 'chr'+chromosome: #e.g. 'chr2'
                        start = int(k[1])
                        stop = int(k[2])
                        depth = v
                        if start < band_stop and stop > band_start:
                            if band in band_depths.keys():
                                band_depths[band].append(depth)
                            else:
                                band_depths[band]=[depth]

    mean_band_depths = {}
    for k, v in band_depths.items():
        start = int(k[0])
        stop = int(k[1])
        if earliest_position <= stop and latest_position >= start: #only plot the mean depths for bands that are covered
            mean_band_depths[k]=np.mean(v)

    return mean_band_depths

def mean_depth_across_chromosome(chromosome, normalised_region_depths):
    chrom = 'chr'+chromosome

    depths = []
    for k, v in normalised_region_depths.items():
        if k[0] == chrom:
            depths.append(v)

    return np.mean(depths), min(depths), max(depths)

def plot_regions_targeted_by_panel(chromosome, panel_coverage, chromosome_sizes, ax, y_max):
    regions_covered = panel_coverage['chr'+chromosome]
    chromosome_size = chromosome_sizes['chr'+chromosome]

    for positions in regions_covered:
        start_position = positions[0]
        stop_position = positions[1]

        bottom = 0
        top = y_max

        x = [start_position, start_position, stop_position, stop_position]
        y = [bottom, top, top, bottom]
        ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with stars

        if stop_position == chromosome_size: #put a line at the end of the plot of final exon targeted
            ax.plot([chromosome_size-1, chromosome_size-1], [0, 800000], color = grey1, lw = 4, zorder = 6) #line at end of plot

    return ax

def individual_chromosome_coverage_depth_plot(chromosome, normalised_region_depths, chromosome_sizes, chromosome_bands, out_directory, sample_name, panel_coverage, chromosome_ideogram):
    plt.close('all')
    f, (ax1, ax2) = plt.subplots(2, 1, sharey=False, sharex = True, figsize=(20, 7))
    gs = matplotlib.gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[10,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    gs.update(hspace=0.05)

    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15
    exonfont=12

    chromosome_size = chromosome_sizes['chr'+chromosome]

    #READ DEPTHS
    all_region_depths = []
    all_start_positions = []
    all_end_positions = []
    for k, v in normalised_region_depths.items():
        if k[0] == 'chr'+chromosome:
            start = int(k[1])
            stop = int(k[2])
            region_depth = float(v)
            all_region_depths.append(region_depth)
            all_start_positions.append(start)
            all_end_positions.append(stop)
            ax1.scatter((stop+start)/2, region_depth, color = '#4292c6', s = 60, zorder = 100)

    if max(all_region_depths)>2:
        y_max = max(all_region_depths)*1.1
    else:
        y_max = 2

    earliest_position = min(all_start_positions)
    latest_position = max(all_start_positions)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel(chromosome, panel_coverage, chromosome_sizes, ax1, y_max)

    #CHROMOSOME IDEOGRAM
    plot_chromosome(chromosome_ideogram, 'chr'+chromosome, ax2)

    #PLOT MEAN DEPTH ACROSS BANDS
    band_depths = mean_depth_across_bands(chromosome, chromosome_bands, normalised_region_depths, earliest_position, latest_position)

    for k, v in band_depths.items():
        band_start = k[0]
        band_stop = k[1]
        band_mean_depth = v
        ax1.plot([band_start, band_stop], [band_mean_depth, band_mean_depth], color = c1, lw = 3, zorder = 500)

    #PLOT MEAN DEPTH ACROSS CHROMOSOME
    mean_depth, min_depth, max_depth = mean_depth_across_chromosome(chromosome, normalised_region_depths)

    ax1.plot([min(all_start_positions), chromosome_size],[1.0, 1.0], color = grey3, linestyle = ':', lw = 3, zorder = 500)

    mean_depth = round(mean_depth, 1)
    min_depth = round(min_depth, 1)
    max_depth = round(max_depth, 1)

    ax1.text(0.02, 0.95, 'mean normalised depth = '+str(mean_depth), transform=ax1.transAxes, fontsize = 13)
    ax1.text(0.02, 0.90, 'min normalised depth = '+str(min_depth), transform=ax1.transAxes, fontsize = 13)
    ax1.text(0.02, 0.85, 'max normalised depth = '+str(max_depth), transform=ax1.transAxes, fontsize = 13)


    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, chromosome_size)
    ax2.set_xlim(0, chromosome_size)

    # Changing the y axis to log scale
    ax1.set_ylim(0, y_max)

    #x-axis ticks
    x1_major_ticks = []
    x1_major_tick_labels = []
    ax1.set_xticks(x1_major_ticks)
    ax1.set_xticklabels(x1_major_tick_labels, fontsize = axisfont)

    ax1.tick_params(axis='y', which='major', labelsize=13)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(True)

    for axis in ['left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)

    legend_elements = [Line2D([0], [0], marker = 's', color='#deebf7', alpha=1.0, markersize = 12, \
                      lw=0, label='regions covered by TWIST CNV panel'),\
                  Line2D([0], [0], marker = '.', color=c1, alpha=1.0, markersize = 0, \
                      lw=3, label='mean normalised read depth across band')]

    ax1.legend(ncol=2, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13)

    #Title and axis labels
    ax1.set_title('Normalised probe region read depths across chr'+chromosome+' (from SSCS): '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('mean read depth in probe region \n(normalised by mean read depth across panel)', fontsize = axislabelfont)

    plt.tight_layout()

    return plt.savefig(out_directory+'/CNV_read_depths/'+sample_name+'_CNV_normalised_read_depths_SSCS_chr'+chromosome+'.pdf')


def plot_regions_targeted_by_panel_cumulative(ax, y_max, cumulative_chromosome_sizes, panel_coverage):

    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        cumulative_chromosome_start = cumulative_chromosome_sizes['chr'+i]
        regions_covered = panel_coverage['chr'+i]

        for positions in regions_covered:
            start_position = positions[0]+cumulative_chromosome_start
            stop_position = positions[1]+cumulative_chromosome_start

            bottom = 0
            top = y_max

            x = [start_position, start_position, stop_position, stop_position]
            y = [bottom, top, top, bottom]
            ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with stars

    return ax

def all_chromosome_depth_plot(normalised_region_depths, cumulative_chromosome_sizes, panel_coverage, out_directory, sample_name, chromosome_sizes, chromosome_bands):
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, sharey=False, sharex = True, figsize=(27, 5))

    axisfont=15
    titlefont=15
    axislabelfont=15

    #READ DEPTHS
    all_region_depths = []
    all_start_positions = []
    all_end_positions = []
    for k, v in normalised_region_depths.items():
        cumulative_chromosome_start = cumulative_chromosome_sizes[k[0]]
        start = int(k[1])+cumulative_chromosome_start
        stop = int(k[2])+cumulative_chromosome_start
        region_depth = float(v)
        all_region_depths.append(region_depth)
        all_start_positions.append(start)
        all_end_positions.append(stop)
        ax1.scatter((stop+start)/2, region_depth, color = '#4292c6', s = 25, zorder = 100)

    if max(all_region_depths)>2:
        y_max = max(all_region_depths)*1.1
    else:
        y_max = 2

    earliest_position = min(all_start_positions)
    latest_position = max(all_start_positions)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel_cumulative(ax1, y_max, cumulative_chromosome_sizes, panel_coverage)

    #PLOT MEAN DEPTH ACROSS BANDS
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        band_depths = mean_depth_across_bands(i, chromosome_bands, normalised_region_depths, earliest_position, latest_position)
        cumulative_chromosome_start = cumulative_chromosome_sizes['chr'+i]
        for k, v in band_depths.items():
            band_start = k[0]+cumulative_chromosome_start
            band_stop = k[1]+cumulative_chromosome_start
            band_mean_depth = v
            ax1.plot([band_start, band_stop], [band_mean_depth, band_mean_depth], color = c1, lw = 3, zorder = 5000)

    ax1.plot([0, cumulative_chromosome_sizes['chrX']+(chromosome_sizes['chrX'])],[1.0, 1.0], color = grey3, linestyle = ':', lw = 3, zorder = 5000)

    #PLOT A VERTICAL LINE TO DISTIGUISH THE CHROMOSOMES
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        if left_start !=0:
            ax1.plot([left_start, left_start], [0, y_max], color = 'k', lw = 1, zorder = 1000)

    end_X = cumulative_chromosome_sizes['chrX']+(chromosome_sizes['chrX'])
    ax1.plot([end_X, end_X], [0, y_max], color = 'k', lw = 1, zorder = 1000)

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, cumulative_chromosome_sizes['chrX']+(chromosome_sizes['chrX'])+100)

    # Changing the y axis to log scale
    ax1.set_ylim(0, y_max)

    #x-axis ticks
    x_axis_tick_positions = []
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        chromosome_size = chromosome_sizes[k]
        mid_point = left_start + (chromosome_size/2)
        x_axis_tick_positions.append(mid_point)

    x_major_ticks = x_axis_tick_positions
    x_major_tick_labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']
    ax1.set_xticks(x_major_ticks)
    ax1.set_xticklabels(x_major_tick_labels, fontsize = axisfont)

    ax1.tick_params(axis='y', which='major', labelsize=13)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(True)

    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)

    legend_elements = [Line2D([0], [0], marker = 's', color='#deebf7', alpha=1.0, markersize = 12, \
                      lw=0, label='regions covered by TWIST CNV panel'),\
                  Line2D([0], [0], marker = '.', color=c1, alpha=1.0, markersize = 0, \
                      lw=3, label='mean normalised read depth across chromosomal band')]

    ax1.legend(ncol=2, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13)

    #Title and axis labels
    ax1.set_title('Normalised probe region read depths (from SSCS): '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('normalised mean read depth)', fontsize = axislabelfont)
    ax1.set_xlabel('chromosome', fontsize = axislabelfont)

    plt.tight_layout()
    return plt.savefig(out_directory+'/CNV_read_depths/'+sample_name+'_CNV_normalised_read_depths_SSCS_all_chromosomes.pdf')


def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mapped SSCS bam (overlapping reads hard-clipped)", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--panel_bed", action="store", dest="panel_bed", help="input bed file for the panel", required=True)
    parser.add_argument("--chromosome-ideogram", action="store", dest="chromosome_ideogram", help="input chromosome ideogram file", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    SSCS_bam = o.infile
    sample_name = o.sample_name
    panel_bed = o.panel_bed
    chromosome_ideogram = o.chromosome_ideogram
    out_directory = o.out_directory

    #Get the chromosome sizes and chromosomal band sizes
    chromosome_sizes, chromosome_bands = chromosome_sizes_bands(chromosome_ideogram)
    #Get panel positions from BED file
    panel_positions = panel_positions_from_bed(panel_bed)
    #Get panel positions from BED file
    panel_coverage = panel_coverage_from_bed(panel_bed)

    #Calculate normalised read depths
    normalised_region_depths = normalised_read_depths(SSCS_bam, panel_positions, chromosome_bands, out_directory, sample_name)

    #PLOT EACH CHROMOSOME SEPARATELY
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        print('plotting read depths across chromosome'+str(i))
        individual_chromosome_coverage_depth_plot(i, normalised_region_depths, chromosome_sizes, chromosome_bands, out_directory, sample_name, panel_coverage, chromosome_ideogram)

    #PLOT ALL THE CHROMOSOMES ON ONE PLOT
    #creating a running start position for each chromosome (so can be plotted sequentially along x axis of plot)
    cumulative_chromosome_sizes = {}
    length = 0
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        chromosome_size = chromosome_sizes['chr'+i]
        cumulative_chromosome_sizes['chr'+i]=length
        length+=chromosome_size

    print('plotting read depths across all chromosomes on one plot')
    all_chromosome_depth_plot(normalised_region_depths, cumulative_chromosome_sizes, panel_coverage, out_directory, sample_name, chromosome_sizes, chromosome_bands)

    return print('Calculation of normalised read depths complete')

if __name__ == "__main__":
	main()
