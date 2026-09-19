#!/usr/bin/env python

'''''
Watson code for detecting KMT2A (MLL) partial tandem duplications.
Version 1.9 (May 2021)

Input:
    1) SSCS bam (overlapping reads hard-clipped)
    2) sample name

Outputs:
    Plots:
        1) normalised read depths across KMT2A (normalised to exon 27 mean read depth)
    Metrics files:
        1) metrics txt file containing info on read depths per position (and normalised read depth)
        2) metrics txt file containing info on mean read depths per intron or exon (and normalised read depth)

Usage:
Watson_code_KMT2A-PTD_calling_v1.9.py --infile SSCS BAM clipped --sample-name SAMPLENAME --KMT2A coordinates CSV file
                        -- panel-bed PANEL BED FILE --out_directory output_file_directory

'''''
version = '1.9'

from argparse import ArgumentParser
import csv
import pysam
import sys
import gzip
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MaxNLocator, MultipleLocator
import numpy as np
import timeit
import time
from datetime import date
from Bio.Seq import Seq

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


def KMT2A_coordinates_dict(KMT2A_coords):
    with open(KMT2A_coords, 'r') as csvfile:
        read_reader = csv.reader(csvfile)  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        KMT2A_coordinates = {}

        for row in read_reader:
            if row_count > 0:
                exon_intron = row[0]
                start = int(row[1])
                end = int(row[2])
                if exon_intron[0]=='E':
                    part = 'exon'
                else:
                    part = 'intron'
                KMT2A_coordinates[exon_intron]=(part, start, end)
            row_count+=1

    return KMT2A_coordinates

def KMT2A_panel_positions_list(KMT2A_coordinates, panel_bed, KMT2A_chromosome, KMT2A_start, KMT2A_stop):
    with open(panel_bed, 'r') as bedfile:
        read_reader = csv.reader(bedfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        KMT2A_panel_positions = []

        for row in read_reader:
            if row_count > 2:
                chromosome = row[0]
                start = int(row[1])
                stop = int(row[2])
                if chromosome == KMT2A_chromosome:
                    if start <  KMT2A_stop and stop > KMT2A_start:
                        for k, v in KMT2A_coordinates.items():
                            intron_exon_start = v[1]
                            intron_exon_stop = v[2]
                            region = k
                            if start < intron_exon_stop and stop > intron_exon_start:
                                KMT2A_panel_positions.append((chromosome, start, stop, region))

            row_count+=1

    return KMT2A_panel_positions

def KMT2A_panel_coverage(panel_bed):
    #panel coverage
    with open(panel_bed, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        panel_coverage = {}
        for row in read_reader:
            if row_count > 2:
                chromosome = row[0]
                start = int(row[1]) #start position of that part of the chromosome
                end = int(row[2]) #end position of that part of the chromosome

                if chromosome in panel_coverage.keys():
                    panel_coverage[chromosome].append((start, end))
                else:
                    panel_coverage[chromosome]=([(start, end)])

            row_count+=1

    return panel_coverage

def KMT2A_ideogram(KMT2A_coordinates):
    color_lookup = {'exon': '#cab2d6', 'intron': 'white'}

    xranges = []
    colors = []
    mid_points = []
    labels = []

    for k, v in KMT2A_coordinates.items():
        start = int(v[1])
        stop = int(v[2])
        label = k.split()[1]
        intron_or_exon = v[0]
        width = stop - start
        mid_point = start + (width/2)
        xranges.append((start, width))
        colors.append(color_lookup[intron_or_exon])
        if intron_or_exon == 'exon':

            mid_points.append(mid_point)
            labels.append(label)

    return xranges, [0, 1.0], colors, mid_points, labels

def plot_gene(KMT2A_coordinates, ax):

    xranges, yrange, colors, midpoints, labels = KMT2A_ideogram(KMT2A_coordinates)

    ax.broken_barh(xranges, yrange, facecolors= colors, edgecolor = 'black')

    ax.set_xticks(midpoints)
    ax.set_xticklabels(labels, rotation = 0, fontsize = 12)
    ax.set_yticks([])
    ax.text(-0.013, 0.35, 'KMT2A', transform=ax.transAxes, fontsize = 15, ha = 'right')
    ax.text(-0.013, -0.55, 'exons:', transform=ax.transAxes, fontsize = 12, ha = 'right')
    ax.xaxis.set_tick_params(width=0.8, color = grey3, length = 6)

    ax.minorticks_off()

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    return ax

def plot_regions_targeted_by_panel(panel_coverage, ax, y_max):

    regions_covered = panel_coverage['chr11']

    for positions in regions_covered:
        start_position = positions[0]
        stop_position = positions[1]

        bottom = 0
        top = y_max

        x = [start_position, start_position, stop_position, stop_position]
        y = [bottom, top, top, bottom]
        ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with stars

    return ax

def read_depths(SSCS_bam, KMT2A_panel_positions, KMT2A_start, KMT2A_stop, KMT2A_coordinates, out_directory, sample_name):
    samfile = pysam.AlignmentFile(SSCS_bam, "rb" )

    all_depths = {}
    all_positions_depths = 0
    total_positions = 0
    intron_exon_depths = {}

    counter = 0
    for region in KMT2A_panel_positions: #iterate over each of the probe regions in the panel...
        start_time = time.time()
        chromosome = region[0]
        start = region[1]
        stop = region[2]
        intron_exon = region[3]

        position_depths = [] #create an empty list in which the depths at each of the positions in the region will be stored.
        for pileupcolumn in samfile.pileup('11', KMT2A_start, KMT2A_stop, max_depth = 100000, stepper = 'all'): #create a pileup of all the reads that cover the probe
            if pileupcolumn.pos in range(start, stop+1): #just look at the pileups of the position in the probe regions (i.e. not the starts of regions that start outside of the probe region)
                depth = pileupcolumn.n #depth at that position
                all_positions_depths+=depth #at the depth to the running list of all depths across the KMT2A
                total_positions +=1 #add 1 to the total number of positions in KMT2A

                all_depths[(chromosome, pileupcolumn.pos, intron_exon)] = depth

                if intron_exon in intron_exon_depths.keys():
                    intron_exon_depths[intron_exon].append((chromosome, pileupcolumn.pos, depth))
                else:
                    intron_exon_depths[intron_exon]=[(chromosome, pileupcolumn.pos, depth)]

        counter+=1
        if counter % 10 == 0:
            print('total probe regions processed: ', counter)
            print('time for last 10 probe regions to be writen = %s seconds' % int(time.time() - start_time))

    samfile.close()

    mean_intron_exon_depths = {}
    for k, v in intron_exon_depths.items():
        depths = []
        start = KMT2A_coordinates[k][1]
        stop = KMT2A_coordinates[k][2]
        for i in v:
            depths.append(i[2])
        mean_depth = np.mean(depths)
        mean_intron_exon_depths[k, start, stop]=mean_depth

    exon_27_mean_depth = mean_intron_exon_depths[('Exon 27', 118373113, 118377361)]

    intron_exon_depths_normalised = {}
    for k, v in intron_exon_depths.items():
        normalised_depths = []
        for i in v:
            normalised_depths.append((i[0], i[1], i[2]/exon_27_mean_depth))
        intron_exon_depths_normalised[k] = normalised_depths

    mean_intron_exon_depths_normalised = {}
    for k, v in mean_intron_exon_depths.items():
        mean_intron_exon_depths_normalised[k]=v/exon_27_mean_depth

    exon3_27ratio = mean_intron_exon_depths_normalised['Exon 3', 118342377, 118345030]/1.0

    read_depth_metrics = open(out_directory+'/KMT2A-PTD/'+sample_name+'_KMT2A_read_depths.txt', 'w')
    read_depth_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    read_depth_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    read_depth_metrics.write('produced from code:\t' + 'Watson_code_KMT2A-PTD_read_depths: version '+str(version) + '\n')
    read_depth_metrics.write('\n')
    read_depth_metrics.write('exon3:27 ratio:\t' + str(exon3_27ratio))
    read_depth_metrics.write('\n')
    read_depth_metrics.write('KMT2A_intron_exon\t' +'chromosome\t' + 'position\t' + 'depth\t' + 'normalised_depth(to_mean_depth_exon_27)\n')
    for k, v in intron_exon_depths.items():
        for i in v:
            read_depth_metrics.write(k+'\t' + str(i[0])+'\t'+ str(i[1])+'\t' + str(i[2])+'\t' + str(i[2]/exon_27_mean_depth) + '\n')
    read_depth_metrics.close()

    mean_read_depth_metrics = open(out_directory+'/KMT2A-PTD/'+sample_name+'_mean_KMT2A_read_depths.txt', 'w')
    mean_read_depth_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    mean_read_depth_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    mean_read_depth_metrics.write('produced from code:\t' + 'Watson_code_KMT2A-PTD_read_depths: version '+str(version) + '\n')
    mean_read_depth_metrics.write('\n')
    mean_read_depth_metrics.write('exon3:27 ratio:\t' + str(exon3_27ratio))
    mean_read_depth_metrics.write('\n')
    mean_read_depth_metrics.write('KMT2A_intron_exon\t' + 'start\t' + 'stop\t' + 'mean_region_depth\t' + 'normalised_mean_region_depth\n')
    for k, v in mean_intron_exon_depths.items():
        normalised_mean_depths = mean_intron_exon_depths_normalised[k]
        mean_read_depth_metrics.write(k[0]+'\t' + str(k[1])+'\t'+ str(k[2])+'\t' + str(v)+'\t' + str(normalised_mean_depths) + '\n')
    mean_read_depth_metrics.close()

    return intron_exon_depths, intron_exon_depths_normalised, mean_intron_exon_depths, mean_intron_exon_depths_normalised, intron_exon_depths_normalised, exon3_27ratio

def read_depths_plot(intron_exon_depths_normalised, mean_intron_exon_depths_normalised, out_directory, sample_name, KMT2A_start, KMT2A_stop, KMT2A_coordinates, panel_coverage, exon3_27ratio):
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

    x = []
    y = []
    for k, v in intron_exon_depths_normalised.items():
        exon_or_intron = k
        for i in v:
            chromosome = i[0]
            position = i[1]
            normalised_depth = i[2]
            x.append(position)
            y.append(normalised_depth)

    ax1.scatter(x, y, s = 10, color = '#4292c6')

    y_max = max(y)*1.2

    ax1.plot([KMT2A_start-500, KMT2A_stop+500], [1, 1], color = grey3, lw = 2, linestyle = ':')

    #gene IDEOGRAM
    plot_gene(KMT2A_coordinates, ax2)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel(panel_coverage, ax1, y_max)

    #PLOT MEAN DEPTH ACROSS EACH EXON:
    for k, v in mean_intron_exon_depths_normalised.items():
        if k[0].split()[0] == 'Exon':
            ax1.plot([k[1], k[2]], [v, v], color = c1, lw = 3)


    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(KMT2A_start-500, KMT2A_stop+500)
    ax2.set_xlim(KMT2A_start-500, KMT2A_stop+500)
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

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    for axis in ['left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)

    legend_elements = [Line2D([0], [0], marker = 's', color='#deebf7', alpha=1.0, markersize = 12, \
                      lw=0, label='regions covered by TWIST CNV panel'),\
                  Line2D([0], [0], marker = '.', color=c1, alpha=1.0, markersize = 0, \
                      lw=3, label='mean normalised read depth across exon')]

    ax1.legend(ncol=2, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13)

    #Title and axis labels
    ax1.set_title('Read depths (normalised to exon 27) across KMT2A: '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('normalised read depth \n(normalised by exon 27 mean read depth)', fontsize = axislabelfont)

    ax1.text(0.01, 0.95, 'exon 3: exon 27 depth ratio = '+str(round(exon3_27ratio, 2)), transform=ax1.transAxes, fontsize = 14, ha = 'left')

    plt.tight_layout()

    return plt.savefig(out_directory+'/KMT2A-PTD/'+sample_name+'_watson_code_sample_normalised_read_depths_KMT2A.pdf')

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input clipped SSCS BAM file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    parser.add_argument("--KMT2A-coordinates", action="store", dest="KMT2A_coords", help="csv file containing coordinates of KMT2A introns and exons", required=True)
    parser.add_argument("--panel-bed", action="store", dest="panel_bed", help="input bed file for the panel", required=True)
    o = parser.parse_args()

    SSCS_bam = o.infile
    sample_name = o.sample_name
    out_directory = o.out_directory
    KMT2A_coords = o.KMT2A_coords
    panel_bed = o.panel_bed

    KMT2A_chromosome = 'chr11'
    KMT2A_start = 118339490
    KMT2A_stop = 118377361

    KMT2A_coordinates = KMT2A_coordinates_dict(KMT2A_coords)
    KMT2A_panel_positions = KMT2A_panel_positions_list(KMT2A_coordinates, panel_bed, KMT2A_chromosome, KMT2A_start, KMT2A_stop)
    panel_coverage = KMT2A_panel_coverage(panel_bed)
    KMT2A_intron_exon_positions = KMT2A_ideogram(KMT2A_coordinates)

    intron_exon_depths, intron_exon_depths_normalised, mean_intron_exon_depths, mean_intron_exon_depths_normalised, intron_exon_depths_normalised, exon3_27ratio = read_depths(SSCS_bam, KMT2A_panel_positions, KMT2A_start, KMT2A_stop, KMT2A_coordinates, out_directory, sample_name)

    read_depths_plot(intron_exon_depths_normalised, mean_intron_exon_depths_normalised, out_directory, sample_name, KMT2A_start, KMT2A_stop, KMT2A_coordinates, panel_coverage, exon3_27ratio)

    return print('KMT2A-PTD calling complete')

if __name__ == "__main__":
	main()
