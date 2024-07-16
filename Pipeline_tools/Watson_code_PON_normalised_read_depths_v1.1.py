#!/usr/bin/env python

'''''
Watson code for calculating and plotting read depths across sample, normalised to panel of normals (PON)(log2ratios).
Version 1.1 (December 202O)

Input:
    1) read depths text file
    2) panel of normals read depths text file
    2) sample name
    3) panel bed file
    4) chromosome ideogram

Outputs:
    Metrics files:
        1) metrics txt file containing mean read depths, normalised read depths, PON read depths, coefficient of variation, log2ratio and p-value for each probe region on the provided panel
    Plots:
        1) plot of log2ratios across each individual chromosome
        2) plot of log2ratios across each chromosome (all chromosomes on one plot)

Usage:
Watson_code_PON_normalised_read_depths.py --infile read depths text file --PON_file text file for PON --sample-name SAMPLENAME --panel_bed  panel bed file
                                --chromosome_ideogram txt file of band positions --out_directory output_file_directory

'''''
version = '1.1'

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
from Bio.Alphabet import IUPAC
import csv
import ast

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

#define the colors from colorbrewer2
orange1 = '#feedde'
orange2 = '#fdbe85'
orange3 = '#fd8d3c'
orange4 = '#e6550d'
orange5 = '#a63603'
purple1 = '#f2f0f7'
purple2 = '#cbc9e2'
purple3 = '#9e9ac8'
purple4 = '#756bb1'
purple5 = '#54278f'
grey1 = '#f7f7f7'
grey2 = '#cccccc'
grey3 = '#969696'
grey4 = '#636363'
grey5 = '#252525'

def chromosome_sizes_bands(chromosome_ideogram):
    with open(chromosome_ideogram, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        chromosome_sizes = {}
        chromosome_bands = {}
        for row in read_reader:
            if row_count > 0:
                chromosome = row[0].replace('chr', '')
                start = int(row[1]) #start position of that part of the chromosome
                end = int(row[2]) #end position of that part of the chromosome
                band = row[3] #e.g. 'q33.1'
                arm = band[0] #e.g. 'q'
                part = row[4] #e.g. 'acen' if centromeric region of chromosome

                chromosome_sizes[chromosome]=end #the dictionary will keep getting overwritten until it gets to the end of the chromosome, so the last end position it adds will be the total length of the chromosome

                if chromosome in chromosome_bands.keys():
                    chromosome_bands[chromosome].append((start, end))
                else:
                    chromosome_bands[chromosome]=[(start, end)]

            row_count+=1

    return chromosome_sizes, chromosome_bands

def panel_coverage_from_bed(panel_bed):
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

def ideograms(ideogram_file, chromosome):
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

def plot_chromosome(ideogram_file, chromosome, ax):

    xranges, yrange, colors, midpoints, labels = ideograms(ideogram_file, chromosome)

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

def plot_regions_targeted_by_panel(chromosome, panel_coverage, chromosome_sizes, ax, y_max, y_min):
    regions_covered = panel_coverage[chromosome]
    chromosome_size = chromosome_sizes[chromosome]

    for positions in regions_covered:
        start_position = positions[0]
        stop_position = positions[1]

        bottom = y_min
        top = y_max

        x = [start_position, start_position, stop_position, stop_position]
        y = [bottom, top, top, bottom]
        ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with stars

        if stop_position == chromosome_size: #put a line at the end of the plot of final exon targeted
            ax.plot([chromosome_size-1, chromosome_size-1], [0, 800000], color = grey1, lw = 4, zorder = 6) #line at end of plot

    return ax

def plot_regions_targeted_by_panel_cumulative(cumulative_chromosome_sizes, panel_coverage, ax, y_max, y_min):
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        cumulative_chromosome_start = cumulative_chromosome_sizes['chr'+i]
        regions_covered = panel_coverage['chr'+i]

        for positions in regions_covered:
            start_position = positions[0]+cumulative_chromosome_start
            stop_position = positions[1]+cumulative_chromosome_start

            bottom = y_min
            top = y_max

            x = [start_position, start_position, stop_position, stop_position]
            y = [bottom, top, top, bottom]
            ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with stars

    return ax

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
                        depth = v[5]
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
            depths.append(v[5])

    return np.mean(depths), min(depths), max(depths)

def chromosome_coverage_depth_plot(chromosome, chromosome_sizes, PON_normalised_read_depths, chromosome_bands, number_controls, sample_name, out_directory, panel_coverage, ideogram_file):
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
    for k, v in PON_normalised_read_depths.items():
        if k[0] == 'chr'+chromosome:
            start = int(k[1])
            stop = int(k[2])
            LRR = float(v[5])
            all_region_depths.append(LRR)
            all_start_positions.append(start)
            all_end_positions.append(stop)
            ax1.scatter((stop+start)/2, LRR, color = '#4292c6', s = 60, zorder = 100)

    if max(all_region_depths)>1 or min(all_region_depths)<-1:
        if max(all_region_depths)>abs(min(all_region_depths)):
            y_max = max(all_region_depths)*1.1
            y_min = -max(all_region_depths)*1.1
        else:
            y_max = abs(min(all_region_depths))*1.1
            y_min = min(all_region_depths)*1.1
    else:
        y_max = 1
        y_min = -1

    earliest_position = min(all_start_positions)
    latest_position = max(all_start_positions)

    #PLOT REGIONS WITH MINIMUM P-VALUE (I.E. VERY DIFFERENT TO CONTROLS)
    minimum_p_value = 1/((number_controls+1)/2)
    for k, v in PON_normalised_read_depths.items():
        if k[0] == 'chr'+chromosome:
            start = int(k[1])
            stop = int(k[2])
            LRR_p_value = float(v[6])
            if LRR_p_value == minimum_p_value:
                ax1.fill([start, start, stop, stop], [0, y_max, y_max, 0], color= orange2, fill = True, alpha = 1.0, linewidth = 1, zorder = 30)

    #PLOT VARIABLE REGIONS (CV > 20)
    minimum_p_value = 1/((number_controls+1)/2)
    for k, v in PON_normalised_read_depths.items():
        if k[0] == 'chr'+chromosome:
            start = int(k[1])
            stop = int(k[2])
            CV_value = float(v[3])
            if CV_value >=20:
                ax1.fill([start, start, stop, stop], [y_min, 0, 0, y_min], color= purple4, fill = True, alpha = 1.0, linewidth = 1, zorder = 40)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel('chr'+str(chromosome), panel_coverage, chromosome_sizes, ax1, y_max, y_min)

    #CHROMOSOME IDEOGRAM
    plot_chromosome(ideogram_file, 'chr'+str(chromosome), ax2)

    #PLOT MEAN DEPTH ACROSS BANDS
    band_depths = mean_depth_across_bands(chromosome, chromosome_bands, PON_normalised_read_depths, earliest_position, latest_position)

    for k, v in band_depths.items():
        band_start = k[0]
        band_stop = k[1]
        band_mean_depth = v
        ax1.plot([band_start, band_stop], [band_mean_depth, band_mean_depth], color = c1, lw = 3, zorder = 500)

    #PLOT MEAN DEPTH ACROSS CHROMOSOME
    mean_depth, min_depth, max_depth = mean_depth_across_chromosome(chromosome, PON_normalised_read_depths)

    ax1.plot([min(all_start_positions), chromosome_size],[0, 0], color = grey3, linestyle = ':', lw = 3, zorder = 500)

    mean_depth = round(mean_depth, 1)
    min_depth = round(min_depth, 1)
    max_depth = round(max_depth, 1)

    ax1.text(0.02, 0.95, 'mean log2ratio = '+str(mean_depth), transform=ax1.transAxes, fontsize = 13, zorder = 200)
    ax1.text(0.02, 0.90, 'min log2ratio = '+str(min_depth), transform=ax1.transAxes, fontsize = 13, zorder = 200)
    ax1.text(0.02, 0.85, 'max log2ratio = '+str(max_depth), transform=ax1.transAxes, fontsize = 13, zorder = 200)

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, chromosome_size)
    ax2.set_xlim(0, chromosome_size)

    # Changing the y axis to log scale
    ax1.set_ylim(y_min, y_max)

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
                      lw=3, label='mean log2ratio across band'),\
                      Line2D([0], [0], marker = '.', color=orange2, alpha=1.0, markersize = 0, \
                      lw=3, label='p-value for log2ratio = '+str(round(minimum_p_value, 2))+' (min possible)'),\
                      Line2D([0], [0], marker = '.', color=purple4, alpha=1.0, markersize = 0, \
                      lw=3, label='controls CV >=20')]

    ax1.legend(ncol=4, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13)

    #Title and axis labels
    ax1.set_title('Probe log2ratio read depths across chr'+chromosome+' (from SSCS): '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('log2ratio normalised read depths', fontsize = axislabelfont)

    plt.tight_layout()
    plt.savefig(out_directory+'/CNV_read_depths/PON_normalised_log2ratios/'+sample_name+'_CNV_logR_ratios_SSCS_chr'+chromosome+'.pdf')

    return plt.show()

def all_chromosome_depth_plot(PON_normalised_read_depths, cumulative_chromosome_sizes, panel_coverage, out_directory, sample_name, chromosome_sizes, chromosome_bands):
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, sharey=False, sharex = True, figsize=(27, 5))

    axisfont=15
    titlefont=15
    axislabelfont=15

    #READ DEPTHS
    all_region_depths = []
    all_start_positions = []
    all_end_positions = []
    for k, v in PON_normalised_read_depths.items():
        cumulative_chromosome_start = cumulative_chromosome_sizes[k[0]]
        start = int(k[1])+cumulative_chromosome_start
        stop = int(k[2])+cumulative_chromosome_start
        LRR = float(v[5])
        all_region_depths.append(LRR)
        all_start_positions.append(start)
        all_end_positions.append(stop)
        ax1.scatter((stop+start)/2, LRR, color = '#4292c6', s = 25, zorder = 100)

    if max(all_region_depths)>1 or min(all_region_depths)<-1:
        if max(all_region_depths)>abs(min(all_region_depths)):
            y_max = max(all_region_depths)*1.1
            y_min = -max(all_region_depths)*1.1
        else:
            y_max = abs(min(all_region_depths))*1.1
            y_min = min(all_region_depths)*1.1
    else:
        y_max = 1
        y_min = -1

    earliest_position = min(all_start_positions)
    latest_position = max(all_start_positions)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel_cumulative(cumulative_chromosome_sizes, panel_coverage, ax1, y_max, y_min)

    #PLOT MEAN DEPTH ACROSS BANDS
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        band_depths = mean_depth_across_bands(i, chromosome_bands, PON_normalised_read_depths, earliest_position, latest_position)
        cumulative_chromosome_start = cumulative_chromosome_sizes['chr'+i]
        for k, v in band_depths.items():
            band_start = k[0]+cumulative_chromosome_start
            band_stop = k[1]+cumulative_chromosome_start
            band_mean_depth = v
            ax1.plot([band_start, band_stop], [band_mean_depth, band_mean_depth], color = c1, lw = 3, zorder = 5000)

    ax1.plot([0, cumulative_chromosome_sizes['chrX']+(chromosome_sizes['chrX'])],[0, 0], color = grey3, linestyle = ':', lw = 3, zorder = 5000)

    #PLOT A VERTICAL LINE TO DISTIGUISH THE CHROMOSOMES
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        if left_start !=0:
            ax1.plot([left_start, left_start], [y_min, y_max], color = 'k', lw = 1, zorder = 1000)

    end_X = cumulative_chromosome_sizes['chrX']+(chromosome_sizes['chrX'])
    ax1.plot([end_X, end_X], [y_min, y_max], color = 'k', lw = 1, zorder = 1000)

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, cumulative_chromosome_sizes['chrX']+(chromosome_sizes['chrX'])+100)

    # Changing the y axis to log scale
    ax1.set_ylim(y_min, y_max)

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
                      lw=3, label='mean log2ratio across band')]

    ax1.legend(ncol=2, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13)

    #Title and axis labels
    ax1.set_title('Probe log2ratio read depths (from SSCS): '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('log2ratio normalised read depths', fontsize = axislabelfont)
    ax1.set_xlabel('chromosome', fontsize = axislabelfont)

    plt.tight_layout()
    return  plt.savefig(out_directory+'/CNV_read_depths/PON_normalised_log2ratios/'+sample_name+'_CNV_logR_ratios_SSCS_all_chromosomes.pdf')

def sample_read_depths(infile):
    with open(infile, 'r') as txtfile:
        read_reader = csv.reader(txtfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        sample_depths = {}

        for row in read_reader:
            if row_count > 5:
                chromosome = 'chr'+row[0]
                start = int(row[1])
                stop = int(row[2])
                band = row[3]
                mean_region_depth = float(row[4])
                normalised_region_depth = float(row[5])
                sample_depths[(chromosome, start, stop, band)] = (mean_region_depth, normalised_region_depth)
            row_count+=1

    return sample_depths

def PON_depths(PON_file):
    with open(PON_file, 'r') as txtfile:
        read_reader = csv.reader(txtfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        PON_read_depths = {}

        for row in read_reader:
            if row_count > 12:
                chromosome = row[0]
                start = int(row[1])
                stop = int(row[2])
                band = row[3]
                PON_mean_depth = float(row[4])
                PON_CV = float(row[5])
                all_controls = ast.literal_eval(row[6])
                number_controls = len(all_controls)
                PON_read_depths[(chromosome, start, stop, band)]=(PON_mean_depth, PON_CV, all_controls)
            row_count+=1

    return PON_read_depths, number_controls

def read_depths_normalised_by_PON(sample_depths, PON_read_depths):
    PON_normalised_read_depths = {}
    for k, v in sample_depths.items():
        mean_region_depth = v[0]
        normalised_region_depth = v[1]
        if k in PON_read_depths.keys():
            PON_mean = PON_read_depths[k][0]
            PON_CV = PON_read_depths[k][1]
            all_PON_read_depths = PON_read_depths[k][2]

            PON_sample_normalised_depth = normalised_region_depth/PON_mean
            log2ratio = np.log2(PON_sample_normalised_depth)
            p_val_LRR = p_value_LRR(all_PON_read_depths, normalised_region_depth)

        PON_normalised_read_depths[k]=(mean_region_depth, normalised_region_depth, PON_mean, PON_CV, PON_sample_normalised_depth, log2ratio, p_val_LRR)

    return PON_normalised_read_depths

def p_value_LRR(PON_read_depths, sample_read_depth): #p-value for log R r atio
    mid_position = len(PON_read_depths)/2

    position = 1
    for i in PON_read_depths:
        if sample_read_depth > i:
            position +=1

    if position > mid_position: #to take in to account the distribution is 2-sided
        position = (len(PON_read_depths)+2) - position

    p_value = position/((len(PON_read_depths)+1)/2)

    return p_value

def metrics_file(sample_name, out_directory, date_today, version, number_controls, PON_normalised_read_depths):
    LRR_metrics = open(out_directory+'/CNV_read_depths/PON_normalised_log2ratios/'+sample_name+'_PON_normalised_read_depths_and_LRR.txt', 'w')
    LRR_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    LRR_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    LRR_metrics.write('produced from code:\t' + 'Watson_code_PON_normalised_read_depths: version '+str(version) + '\n')
    LRR_metrics.write('total number of controls in panel of normals (PON) used for normalisation:\t' + str(number_controls) + '\n')
    LRR_metrics.write('\n')
    LRR_metrics.write('chromosome\t' + 'start\t' + 'stop\t' + 'band\t' + 'mean_region_depth\t' + 'normalised_region_depth\t' + 'mean_normalised_depth_PON\t' +\
                     'coefficient_of_variation_PON\t'+ 'sample_depth_normalised_by_PON\t' + 'log2ratio\t' + 'p-value\n')
    for k, v in PON_normalised_read_depths.items():
        LRR_metrics.write(k[0]+'\t' + str(k[1])+'\t'+ str(k[2])+'\t'+ str(k[3])+ '\t' + str(v[0])+'\t' + str(v[1])+'\t'+ str(v[2]) +\
                          '\t' + str(v[3]) + '\t' + str(v[4]) + '\t' + str(v[5]) + '\t' + str(v[6]) + '\n')
    LRR_metrics.close()

    return 'logRratio metrics file written for sample '+sample_name

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input sample normalised read depths text file", required=True)
    parser.add_argument("--PON_file", action="store", dest="PON_file", help="input panel of normals read depths text file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--panel_bed", action="store", dest="panel_bed", help="input bed file for the panel", required=True)
    parser.add_argument("--chromosome-ideogram", action="store", dest="chromosome_ideogram", help="input chromosome ideogram file", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    infile = o.infile
    PON_file = o.PON_file
    sample_name = o.sample_name
    panel_bed = o.panel_bed
    chromosome_ideogram = o.chromosome_ideogram
    out_directory = o.out_directory

    #Get the chromosome sizes and chromosomal band sizes
    chromosome_sizes, chromosome_bands = chromosome_sizes_bands(chromosome_ideogram)
    #Get panel positions from BED file
    panel_coverage = panel_coverage_from_bed(panel_bed)

    #Get sample-normalised read depths
    sample_depths = sample_read_depths(infile)

    #Get PON read depths
    PON_read_depths, number_controls = PON_depths(PON_file)

    #Calculate read depths normalised by PON
    PON_normalised_read_depths = read_depths_normalised_by_PON(sample_depths, PON_read_depths)

    #Create text file containing the logR ratios, p-values and normalised read depths
    metrics_file(sample_name, out_directory, date_today, version, number_controls, PON_normalised_read_depths)

    #PLOT EACH CHROMOSOME SEPARATELY
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        print('plotting log2ratios across chromosome'+str(i))
        chromosome_coverage_depth_plot(i, chromosome_sizes, PON_normalised_read_depths, chromosome_bands, number_controls, sample_name, out_directory, panel_coverage, chromosome_ideogram)

    #PLOT ALL THE CHROMOSOMES ON ONE PLOT
    #creating a running start position for each chromosome (so can be plotted sequentially along x axis of plot)
    cumulative_chromosome_sizes = {}
    length = 0
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        chromosome_size = chromosome_sizes['chr'+i]
        cumulative_chromosome_sizes['chr'+i]=length
        length+=chromosome_size

    print('plotting log2ratios across all chromosomes on one plot')
    all_chromosome_depth_plot(PON_normalised_read_depths, cumulative_chromosome_sizes, panel_coverage, out_directory, sample_name, chromosome_sizes, chromosome_bands)

    return print('Calculation of normalised read depths complete')

if __name__ == "__main__":
	main()
