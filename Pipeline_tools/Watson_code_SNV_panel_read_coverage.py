#!/usr/bin/env python

'''''
Watson code for plotting read depths across SNV panel.
Version 1.2 (June 2021)

Input:
    1) mapped bam
    2) sample name
    3) panel bed file

Outputs:
    Metrics files:
        1) metrics txt file containing read depths for each position on the SNV panel as well as mean panel depth
    Plots:
        1) plot of read depths across each targeted gene and 4 risk SNPs

Usage:
Watson_code_SNV_panel_read_coverage.py --infile mapped BAM --sample-name SAMPLENAME --panel_bed  panel bed file
                                    --BAM_type whether 'raw bam' 'SSCS' or 'duplex' --out_directory output_file_directory
                                    --min-family-size-SSCS minimum SSCS UMI family size used (for file savename)

'''''
version = '1.2'

from argparse import ArgumentParser
import pysam
import sys
import math
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

colors = [c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2, c4, c0, c3, c1, c5, \
          c6,c7, c8, c9, c10, c2,c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2,\
         c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2, c4, c0, c3, c1, c5, \
          c6,c7, c8, c9, c10, c2,c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2,\
         c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2, c4, c0, c3, c1, c5, \
          c6,c7, c8, c9, c10, c2,c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2]

def nucleotide_dictionary(exon_dictionary, strand): #create a dictionary where the key = coordinate and value = corresponding gene nucleotide number (rather than chromosomal position), so can plot the exons adjacent to each other in the plot
    nucleotide_dict = {}
    nucleotide = 1
    for k, v in exon_dictionary.items():
        if strand == '+':
            start = v[0]
            stop = v[1]
            positions = np.linspace(start-1, stop+1, (stop-start)+3)
            for i in positions:
                nucleotide_dict[i]=nucleotide
                nucleotide+=1
        if strand == '-':
            start = -v[0]
            stop = -v[1]
            positions = -np.linspace(start-1, stop+1, (stop-start)+3)
            for i in positions:
                nucleotide_dict[i]=nucleotide
                nucleotide+=1
    return nucleotide_dict

def exon_size_dictionary(exon_dictionary, strand): #create a dictionary where the key = exon number and the value = number of nucleotides
    exon_size_dict = {}
    for k, v in exon_dictionary.items():
        exon = k
        start = v[0]
        stop = v[1]

        if strand == '+':
            exon_size = (stop-start)+1
            exon_size_dict[str(exon)]=exon_size

        if strand == '-':
            exon_size = (start-stop)+1
            exon_size_dict[str(exon)]=exon_size

    return exon_size_dict

def exon_order(exon_size_dict): #create a list of the exons (e.g. 2, 3, 4, 5, 6, 7, 8)
    exon_orders_list = []
    for k, v in exon_size_dict.items():
        exon_orders_list.append(str(k))
    return exon_orders_list

def exon_nucleotide_start(exon_dictionary, nucleotide_dictionary): # create a list of the start positions of each exon in the gene (by gene nucleotide number)
    left_pos = []
    for k, v in exon_dictionary.items():
        position = v[0]
        nucleotide_number = nucleotide_dictionary[position]
        left_pos.append(nucleotide_number)

    return left_pos

def exons_targeted_by_panel(exon_dictionary, exons_targeted_list, nucleotide_dictionary, ax, y_max):

    final_exon_end = list(exon_dictionary.values())[-1][-1]
    final_exon_end_nucleotide = nucleotide_dictionary[final_exon_end]

    for i in exons_targeted_list:
        start_stop_positions = exon_dictionary[i]

        start_position = start_stop_positions[0]
        stop_position = start_stop_positions[1]

        start_nucleotide = nucleotide_dictionary[start_position]
        stop_nucleotide = nucleotide_dictionary[stop_position]

        bottom = 0
        top = y_max

        x = [start_nucleotide, start_nucleotide, stop_nucleotide, stop_nucleotide]
        y = [bottom, top, top, bottom]
        ax.fill(x, y, color= '#deebf7', fill = True, alpha = 0.5, linewidth = 1, zorder = 0) #fill in the box with stars
        ax.plot([start_nucleotide-1, start_nucleotide-1], [0, y_max], lw = 1, color = grey2)
        ax.plot([stop_nucleotide+1, stop_nucleotide+1], [0, y_max], lw = 1, color = grey2)

        if stop_nucleotide == final_exon_end_nucleotide: #put a line at the end of the plot of final exon targeted
            ax.plot([len(nucleotide_dictionary)-1, len(nucleotide_dictionary)-1], [0, y_max], color = grey2, lw = 4, zorder = 6) #line at end of plot

    return ax

def mean_exon_read_depths(read_depths_dict, exon_dictionary, strand):
    exon_read_depths = {}
    exon_position_read_depths = {}
    for k, v in read_depths_dict.items():
        chromosome = k[0]
        position = k[1]
        depth = v
        for a, b in exon_dictionary.items():
            if strand == '-':
                start = b[1]
                stop = b[0]
                if position >= start and position <=stop:
                    exon = a
                    if exon in exon_read_depths.keys():
                        exon_read_depths[exon].append(depth)
                    else:
                        exon_read_depths[exon]=[depth]

                    if exon in exon_position_read_depths.keys():
                        exon_position_read_depths[exon].append((position, depth))
                    else:
                        exon_position_read_depths[exon]=[(position, depth)]
            if strand == '+':
                start = b[0]
                stop = b[1]
                if position >= start and position <=stop:
                    exon = a
                    if exon in exon_read_depths.keys():
                        exon_read_depths[exon].append(depth)
                    else:
                        exon_read_depths[exon]=[depth]

                    if exon in exon_position_read_depths.keys():
                        exon_position_read_depths[exon].append((position, depth))
                    else:
                        exon_position_read_depths[exon]=[(position, depth)]

    mean_read_depths = {}
    for k, v in exon_read_depths.items():
        mean_read_depths[k]=np.mean(v)

    return mean_read_depths, exon_position_read_depths

def panel_positions_from_bed(panel_bed):
    with open(panel_bed, 'r') as bedfile:
        read_reader = csv.reader(bedfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        panel_positions = []

        for row in read_reader:
            if row_count > 1:
                chromosome = row[0] #e.g. 'chr2'
                start = row[1]
                stop = row[2]
                panel_positions.append((chromosome, int(start), int(stop)))
            row_count+=1

    return panel_positions

def gene_panel_positions(panel_positions, exon_dictionary, chromosome, strand):
    gene_regions_targeted = []
    for i in panel_positions:
        chrom = i[0].replace('chr', '')
        start = i[1]
        stop = i[2]
        if chrom == chromosome:
            for k, v in exon_dictionary.items():
                exon = k
                if strand == '-':
                    beginning = v[1]
                    end = v[0]
                    if start < end and stop > beginning:
                        gene_regions_targeted.append(i)
                if strand == '+':
                    beginning = v[0]
                    end = v[1]
                    if start < end and stop > beginning:
                        gene_regions_targeted.append(i)

    return gene_regions_targeted

def position_read_depths(bam, gene_panel_positions, gene_name, exon_dictionary, strand):
    samfile = pysam.AlignmentFile(bam, "rb")

    position_depths = {}
    all_positions_depths = 0
    total_positions = 0

    counter = 0
    for region in gene_panel_positions: #iterate over each of the probe regions in the panel...
        start_time = time.time()
        chromosome = region[0].replace('chr', '')
        start = region[1]
        stop = region[2]

        for pileupcolumn in samfile.pileup(chromosome, start, stop, max_depth = 1000000, stepper = 'all'): #create a pileup of all the reads that cover the probe
            if pileupcolumn.pos in range(start, stop+1): #just look at the pileups of the position in the probe regions (i.e. not the starts of regions that start outside of the probe region)
                depth = pileupcolumn.n #depth at that position
                position = pileupcolumn.pos
                exon = ''
                for k, v in exon_dictionary.items():
                    if strand == '-':
                        if v[1] <= position <= v[0]:
                            exon = 'exon '+str(k)
                    if strand == '+':
                        if v[0] <= position <= v[1]:
                            exon = 'exon '+str(k)
                position_depths[(chromosome, position, gene_name, exon)] = depth #record the depth at each position
                all_positions_depths+=depth #add the depth to the running list of all depths across the panel
                total_positions +=1 #add 1 to the total number of positions in the panel

        counter+=1
        if counter % 100 == 0:
            print('total probe regions processed: ', counter)
            print('time for last 100 probe regions to be writen = %s seconds' % int(time.time() - start_time))

    samfile.close()

    if total_positions >0:
        mean_depth = all_positions_depths/total_positions
    else:
        mean_depth = 0

    return position_depths, mean_depth

def gene_depths(in_bam, panel_positions, exon_dictionary, chromosome, strand, gene_name):
    panel_positions_gene = gene_panel_positions(panel_positions, exon_dictionary, chromosome, strand)
    read_depths_dict, mean_coverage = position_read_depths(in_bam, panel_positions_gene, gene_name, exon_dictionary, strand)

    return read_depths_dict, mean_coverage

def y_max_calculator(all_read_depths, rounding_number):

    def roundup(x, rounding_number):
        return int(math.ceil(x/rounding_number))*rounding_number #e.g. round up to nearest 50,000: rounding number = 50000

    max_depth = 0
    for depth in all_read_depths.values():
        if depth > max_depth:
            max_depth = depth

    return roundup(max_depth, rounding_number)

def gene_coverage_depth_plot(in_bam, read_depths_dict, mean_coverage, exon_dictionary, chromosome, strand, gene_name, bam_type, sample_name, exons_targeted_list, y_max, mean_depth_panel, out_directory, min_family_size):

    plt.close('all')
    f, (ax1, ax2) = plt.subplots(2, 1, sharey=True, sharex = True, figsize=(15, 7))
    gs = matplotlib.gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[10,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    gs.update(hspace=0.05)

    colors = [c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2, c4, c0, c3, c1, c5, \
              c6,c7, c8, c9, c10, c2,c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2,\
             c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2, c4, c0, c3, c1, c5, \
              c6,c7, c8, c9, c10, c2,c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2,\
             c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2, c4, c0, c3, c1, c5, \
              c6,c7, c8, c9, c10, c2,c4, c0, c3, c1, c5, c6, c7, c8, c9, c10, c2]

    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15
    exonfont=12

    nucleotides_dict = nucleotide_dictionary(exon_dictionary, strand)
    exons_sizes = exon_size_dictionary(exon_dictionary, strand)
    exons_orders = exon_order(exons_sizes)
    left_positions = exon_nucleotide_start(exon_dictionary, nucleotides_dict)

    #Plot mean exon read depths
    mean_depths, exon_read_depths = mean_exon_read_depths(read_depths_dict, exon_dictionary, strand)
    for k, v in exon_read_depths.items():
        x = []
        y = []
        exon = k
        positions = v
        for i in positions:
            x.append(nucleotides_dict[i[0]])
            y.append(i[1])
        ax1.plot(x, y, color = '#4292c6', lw = 3)

    for k, v in mean_depths.items():
        exon_start = nucleotides_dict[exon_dictionary[k][0]]
        exon_stop = nucleotides_dict[exon_dictionary[k][1]]
        mean = v
        ax1.plot([exon_start, exon_stop], [mean, mean], color = c1, lw = 3)

    ax1.text(0.98, 0.925, 'mean '+gene_name+' coverage = '+str(round(int(mean_coverage), 0)), ha = 'right', fontsize = 14, transform=ax1.transAxes, zorder = 50,
            bbox=dict(facecolor='white', edgecolor=grey2, boxstyle='round', pad=0.5))
    ax1.text(0.02, 0.925, 'chr '+chromosome, ha = 'left', fontsize = 14, transform=ax1.transAxes, zorder = 50,
            bbox=dict(facecolor='white', edgecolor=grey2, boxstyle='round', pad=0.5))

    #Plot mean depth across whole panel
    ax1.plot([0, len(nucleotides_dict)], [mean_depth_panel, mean_depth_panel], linestyle = ':', color = '#9e9ac8', lw = 3, zorder = 50)

    #PLOT EXONS
    xranges = []
    color_list = []
    mid_points = []
    labels = []

    n = 0
    for k, v in exon_dictionary.items():
        start = nucleotides_dict[int(v[0])]
        stop = nucleotides_dict[int(v[1])]
        width = stop - start
        mid_point = start + (width/2)
        xranges.append((start-1, width+2))
        color_list.append(colors[n])
        mid_points.append(mid_point)
        labels.append(k)
        n+=1

    ax2.broken_barh(xranges, [0, -10], facecolors= color_list, edgecolor = 'black')

    for label, position in zip(labels, mid_points):
        ax2.text(position, -5, label, ha='center', va='center', fontsize = exonfont)

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, len(nucleotides_dict))
    ax1.set_ylim(0, y_max)
    ax2.set_xlim(0, len(nucleotides_dict))

    #x-axis ticks
    x1_major_ticks = []
    x1_major_tick_labels = []
    x2_major_ticks = []
    x2_major_tick_labels = []
    ax1.set_xticks(x1_major_ticks)
    ax1.set_xticklabels(x1_major_tick_labels, fontsize = axisfont)
    ax2.set_xticks(x2_major_ticks)
    ax2.set_xticklabels(x2_major_tick_labels, fontsize = axisfont)

    #y-axis ticks
    y2_major_ticks = []
    y2_major_tick_labels = []
    ax2.set_yticks(y2_major_ticks)
    ax2.set_yticklabels(y2_major_tick_labels, fontsize = axisfont)

    ax1.tick_params(axis='y', which='major', labelsize=14)

    #Only show the required axis lines
    for ax in (ax1, ax2):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

    ax2.spines['left'].set_visible(False)

    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')

    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)

    #plot the regions targeted by the panel
    exons_targeted_by_panel(exon_dictionary, exons_targeted_list, nucleotides_dict, ax1, y_max)

    ax2.plot([len(nucleotides_dict)-1, len(nucleotides_dict)-1], [-10, -1], color = 'k', lw = 2, zorder = 6) #line at end of exons

    legend_elements = [Line2D([0], [0], marker = 's', color='#deebf7', alpha=1.0, markersize = 12, \
                      lw=0, label='regions covered by TWIST SNV panel'),\
                  Line2D([0], [0], marker = '.', color=c1, alpha=1.0, markersize = 0, \
                      lw=3, label='mean read depth across exon'),\
                      Line2D([0], [0], marker = '.', color='#9e9ac8', alpha=1.0, markersize = 0, \
                      lw=3, linestyle = ':', label='mean read depth across panel')]

    ax1.legend(ncol=3, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 14)

    #Title and axis labels
    ax1.set_title('Read depth across '+gene_name+': '+sample_name+' ('+bam_type+')', y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('read depth', fontsize = axislabelfont)
    ax2.set_xlabel('exons', fontsize = axislabelfont)

    save_type = bam_type.replace(' ', '_')

    if save_type in ['SSCS', 'duplex', 'DCS']:
        plt.savefig(out_directory+'/Gene_read_depths/'+sample_name+'_watson_code_read_depths_'+gene_name+'_'+save_type+'_MUFS'+str(min_family_size)+'.pdf')
    else:
        plt.savefig(out_directory+'/Gene_read_depths/'+sample_name+'_watson_code_read_depths_'+gene_name+'_'+save_type+'.pdf')

    return

def SNP_depths(in_bam, panel_positions, SNP_chromosome, SNP_position, SNP_name):
    SNP_probe = []
    for i in panel_positions:
        chromosome = i[0]
        start = i[1]
        stop = i[2]
        if chromosome == SNP_chromosome:
            if start <= SNP_position <= stop:
                SNP_probe.append(i)

    samfile = pysam.AlignmentFile(in_bam, "rb")
    position_depths = {}
    all_positions_depths = 0
    total_positions = 0

    counter = 0
    for region in SNP_probe: #iterate over each of the probe regions in the panel...
        start_time = time.time()
        chromosome = region[0]
        start = region[1]
        stop = region[2]

        for pileupcolumn in samfile.pileup(chromosome, start, stop, max_depth = 1000000, stepper = 'all'): #create a pileup of all the reads that cover the probe
            if pileupcolumn.pos in range(start, stop+1): #just look at the pileups of the position in the probe regions (i.e. not the starts of regions that start outside of the probe region)
                depth = pileupcolumn.n #depth at that position
                position_depths[(chromosome, pileupcolumn.pos, SNP_name+' probe region')] = depth #record the depth at each position
                all_positions_depths+=depth #add the depth to the running list of all depths across the panel
                total_positions +=1 #add 1 to the total number of positions in the panel

        counter+=1
    samfile.close()

    if total_positions > 0:
        mean_depth = all_positions_depths/total_positions
    else:
        mean_depth = 0

    return position_depths, mean_depth

def SNP_coverage_depth_plot(in_bam, position_depths, mean_depth, SNP_chromosome, SNP_position, SNP_name, bam_type, sample_name, y_max, mean_depth_panel, out_directory, min_family_size):

    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, sharey=True, sharex = False, figsize=(15, 6.3))

    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15
    exonfont=12

    #Plot mean exon read depths
    x = []
    y = []
    for k, v in position_depths.items():
        chromosome = k[0]
        position = k[1]
        depth = v
        x.append(position)
        y.append(depth)
    ax1.plot(x, y, color = '#4292c6', lw = 3)

    if len(x)>0:
        min_x = min(x)
        max_x = max(x)
    else:
        min_x = SNP_position-60
        max_x = SNP_position+60

    ax1.plot([min_x, max_x], [mean_depth, mean_depth], color = c1, lw = 3)

    #Plot mean depth across whole panel
    ax1.plot([min_x, max_x], [mean_depth_panel, mean_depth_panel], linestyle = ':', color = '#9e9ac8', lw = 3, zorder = 50)

    ax1.text(0.98, 0.925, 'mean coverage around '+SNP_name+' = '+str(round(int(mean_depth), 0)), ha = 'right', fontsize = 14, transform=ax1.transAxes, zorder = 50,
            bbox=dict(facecolor='white', edgecolor=grey2, boxstyle='round', pad=0.5))
    ax1.text(0.02, 0.925, 'chr '+SNP_chromosome, ha = 'left', fontsize = 14, transform=ax1.transAxes, zorder = 50,
            bbox=dict(facecolor='white', edgecolor=grey2, boxstyle='round', pad=0.5))

    #PLOT SNP
    ax1.plot([SNP_position, SNP_position], [0, y_max*0.83], color = grey4, lw = 3, linestyle = ':')
    ax1.text(SNP_position, y_max*0.85, SNP_name, fontsize = 14, ha = 'center')

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limit
    ax1.set_xlim(min_x, max_x)
    ax1.set_ylim(0, y_max)

    #x-axis ticks
    x1_major_ticks = [min_x, SNP_position, max_x]
    x1_major_tick_labels = [str(min_x), str(SNP_position), str(max_x)]
    ax1.set_xticks(x1_major_ticks)
    ax1.set_xticklabels(x1_major_tick_labels, fontsize = axisfont)

    #y-axis ticks
    ax1.tick_params(axis='y', which='major', labelsize=14)
    ax1.tick_params(axis='x', which='both', labelsize=14)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')

    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)

    #plot the regions targeted by the panel
    bottom = 0
    top = y_max
    x = [min_x, min_x, max_x, max_x]
    y = [bottom, top, top, bottom]

    ax1.fill(x, y, color= '#deebf7', fill = True, alpha = 0.5, linewidth = 1, zorder = 0) #fill in the box with stars
    ax1.plot([min_x-1, min_x-1], [0, y_max], lw = 1, color = grey2)
    ax1.plot([max_x+1, max_x+1], [0, y_max], lw = 1, color = grey2)

    legend_elements = [Line2D([0], [0], marker = 's', color='#deebf7', alpha=1.0, markersize = 12, \
                      lw=0, label='regions covered by TWIST SNV panel'),\
                  Line2D([0], [0], marker = '.', color=c1, alpha=1.0, markersize = 0, \
                      lw=3, label='mean read depth across exon'),\
                      Line2D([0], [0], marker = '.', color='#9e9ac8', alpha=1.0, markersize = 0, \
                      lw=3, linestyle = ':', label='mean read depth across panel')]

    ax1.legend(ncol=3, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 14)

    #Title and axis labels
    ax1.set_title('Read depth across '+SNP_name+' region: '+sample_name+' ('+bam_type+')', y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('read depth', fontsize = axislabelfont)
    ax1.set_xlabel('chromosomal position', fontsize = axislabelfont)

    save_type = bam_type.replace(' ', '_')

    if save_type in ['SSCS', 'duplex', 'DCS']:
        plt.savefig(out_directory+'/Gene_read_depths/'+sample_name+'_watson_code_read_depths_'+SNP_name+'_'+save_type+'_MUFS'+str(min_family_size)+'.pdf')
    else:
        plt.savefig(out_directory+'/Gene_read_depths/'+sample_name+'_watson_code_read_depths_'+SNP_name+'_'+save_type+'.pdf')

    return

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mapped SSCS bam (overlapping reads hard-clipped)", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--panel_bed", action="store", dest="panel_bed", help="input bed file for the panel", required=True)
    parser.add_argument("--min-family-size-SSCS", type=int, action="store", dest="min_family_size", help="min UMI family size used when calling SSCSS", required=True)
    parser.add_argument("--BAM_type", type=str, dest='bam_type', help="whether 'raw bam', 'SSCS' or 'duplex'", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    in_bam = o.infile
    sample_name = o.sample_name
    panel_bed = o.panel_bed
    bam_type = o.bam_type
    out_directory = o.out_directory
    min_family_size = o.min_family_size

    panel_positions = panel_positions_from_bed(panel_bed)

    ###### CALCULATE READ DEPTHS ########
    #ASXL1
    print('calculating read depths across ASXL1...')
    ASXL1_exons = {1: (30947550, 30947594), 2: (30954187, 30954269), 3: (30956818, 30956926), 4: (31015931, 31016051),
                   5: (31016128, 31016225), 6: (31017141, 31017234), 7: (31017704, 31017856), 8: (31019124, 31019287),
                   9: (31019386, 31019482), 10: (31020683, 31020788), 11: (31021087, 31021720), 12: (31022235, 31027121-1980)}
    ASXL1_exons_targeted = [11, 12]
    ASXL1_depths, ASXL1_mean_coverage = gene_depths(in_bam, panel_positions, ASXL1_exons, '20', '+', 'ASXL1')

    #BCOR
    print('calculating read depths across BCOR...')
    BCOR_exons = {1: (39956656, 39956468), 2: (39937222, 39937097), 3: (39935785, 39935707), 4: (39934433, 39931602), 5: (39930943, 39930890),
                  6: (39930412, 39930226), 7:  (39923852, 39923589), 8: (39923205, 39922861), 9: (39922324, 39921999), 10: (39921646, 39921392),
                  11: (39916574, 39916408), 12: (39914766, 39914621), 13: (39913586, 39913509), 14: (39913295, 39913139), 15: (39911653, 39910501)}
    BCOR_exons_targeted = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    BCOR_depths, BCOR_mean_coverage = gene_depths(in_bam, panel_positions, BCOR_exons, 'X', '-', 'BCOR')

    #BCORL1
    print('calculating read depths across BCORL1...')
    BCORL1_exons = {1: (129139208, 129139293), 2: (129146554, 129146644), 3: (129146926, 129150189), 4: (129154960, 129155125),\
                    5: (129156872, 129156952), 6: (129158965, 129159354), 7: (129162610, 129162836), 8: (129171342, 129171508),\
                    9: (129173112, 129173257), 10: (129184692, 129184769), 11: (129185835, 129185991), 12: (129189829, 129190111)}
    BCORL1_exons_targeted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    BCORL1_depths, BCORL1_mean_coverage = gene_depths(in_bam, panel_positions, BCORL1_exons, 'X', '+', 'BCORL1')

    #CBL
    print('calculating read depths across CBL...')
    CBL_exons= {1: (119077128, 119077322), 2: (119103158, 119103405), 3: (119142445, 119142591), 4: (119144578, 119144734), 5: (119145542, 119145663),\
                6: (119146707, 119146844), 7: (119148467, 119148554), 8: (119148876, 119149007), 9: (119149220, 119149423), 10: (119155679, 119155810),\
                11: (119155899, 119156276), 12: (119158562, 119158656), 13: (119167628, 119167744), 14: (119168094, 119168191), 15: (119169068, 119169250),\
                16: (119170205, 119170491)}
    CBL_exons_targeted = [8, 9, 16]
    CBL_depths, CBL_mean_coverage = gene_depths(in_bam, panel_positions, CBL_exons, '11', '+', 'CBL')

    #CEBPA
    print('calculating read depths across CEBPA...')
    CEBPA_exons = {1: (33793320, 33792244)}
    CEBPA_exons_targeted = [1]
    CEBPA_depths, CEBPA_mean_coverage = gene_depths(in_bam, panel_positions, CEBPA_exons, '19', '-', 'CEBPA')

    #CHEK2
    print('calculating read depths across CHEK2...')
    CHEK2_exons = {1: (29137832, 29137757), 2: (29130715, 29130391), 3: (29121355, 29121231), 4:  (29121112, 29120965),\
                   5: (29115473, 29115383), 6: (29108005, 29107897), 7: (29106047, 29105994), 8: (29099554, 29099493),\
                   9: (29095925, 29095826), 10: (29092975, 29092889), 11: (29091861, 29091698), 12: (29091230, 29091115),\
                   13: (29090105, 29090020), 14: (29085203, 29085123), 15: (29083974, 29083732)}
    CHEK2_exons_targeted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    CHEK2_depths, CHEK2_mean_coverage = gene_depths(in_bam, panel_positions, CHEK2_exons, '22', '-', 'CHEK2')

    #CSF3R
    print('calculating read depths across CSF3R...')
    CSF3R_exons = {3: (36945117, 36945034), 4: (36941274, 36940978), 5: (36939488, 36939365), 6: (36939223, 36939036), 7: (36938287, 36938118),\
                   8: (36937992, 36937839), 9: (36937740, 36937667), 10: (36937247, 36937034), 11: (36935441, 36935253),\
                   12: (36934858, 36934757), 13: (36933822, 36933676), 14: (36933563, 36933423), 15:  (36933252, 36933159),\
                   16: (36932912, 36932831), 17: (36932428, 36931644+314)}
    CSF3R_exons_targeted = [14, 15, 16, 17]
    CSF3R_depths, CSF3R_mean_coverage = gene_depths(in_bam, panel_positions, CSF3R_exons, '1', '-', 'CSF3R')

    #DDX41
    print('calculating read depths across DDX41...')
    DDX41_exons = {1: (176944470-524, 176943920), 2: (176943836, 176943726), 3: (176943448, 176943289), 4: (176943194, 176943120),\
                   5: (176942990, 176942930), 6: (176942822, 176942686), 7: (176942259, 176942187), 8: (176942070, 176941917),\
                   9: (176941838, 176941702), 10: (176940848, 176940686), 11: (176940485, 176940354), 12: (176940083, 176940012),\
                   13: (176939877, 176939781), 14:  (176939646, 176939497), 15: (176939394, 176939323), 16: (176939207, 176939097),\
                   17: (176938928, 176938578+213)}
    DDX41_exons_targeted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    DDX41_depths, DDX41_mean_coverage = gene_depths(in_bam, panel_positions, DDX41_exons, '5', '-', 'DDX41')

    #DNMT3A
    print('calculating read depths across DNMT3A...')
    DNMT3A_exons= {2: (25536853, 25536782), 3: (25523112, 25523008), 4: (25505580, 25505310), 5: (25498412, 25498369),\
              6: (25497956, 25497810), 7: (25471121, 25470906), 8: (25470618, 25470460), 9: (25470027, 25469920), 10: (25469645, 25469489),\
              11: (25469178, 25469029), 12: (25468933, 25468889), 13: (25468201, 25468122), 14: (25467521, 25467409), 15: (25467207, 25467024),\
              16: (25466851, 25466767), 17: (25464576, 25464431), 18: (25463599, 25463509), 19: (25463319, 25463171), 20: (25462084, 25461999),\
              21: (25459874, 25459805), 22: (25458694, 25458576), 23: (25457289, 25457148)}
    DNMT3A_exons_targeted = list(np.linspace(2, 23, (23-2)+1))
    DNMT3A_depths, DNMT3A_mean_coverage = gene_depths(in_bam, panel_positions, DNMT3A_exons, '2', '-', 'DNMT3A')

    #EZH2
    print('calculating read depths across EZH2...')
    EZH2_exons = {2: (148544404, 148544274), 3: (148543690, 148543562), 4: (148529842, 148529726), 5: (148526940, 148526820),\
                  6: (148525972, 148525832), 7:  (148524358, 148524256), 8: (148523724, 148523546), 9: (148516779, 148516688),\
                  10: (148515209, 148514969), 11: (148514483, 148514314), 12: (148513870, 148513776), 13: (148512638, 148512598),\
                  14: (148512131, 148512006), 15: (148511229, 148511051), 16: (148508812, 148508717), 17: (148507506, 148507425),\
                  18:  (148506482, 148506402), 19: (148506247, 148506163), 20: (148504798, 148504477+261)}
    EZH2_exons_targeted = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    EZH2_depths, EZH2_mean_coverage = gene_depths(in_bam, panel_positions, EZH2_exons, '7', '-', 'EZH2')

    #FLT3
    print('calculating read depths across FLT3...')
    FLT3_exons = {1: (28674647, 28674605), 2: (28644749, 28644628), 3: (28636206, 28636004), 4: (28631599, 28631484), 5: (28626811, 28626682), 6: (28624359, 28624232),\
                  7: (28623911, 28623772), 8: (28623674, 28623521), 9: (28622580, 28622412), 10: (28611425, 28611322), 11: (28610180, 28610072), 12: (28609810, 28609632),\
                  13: (28608544, 28608438), 14: (28608351, 28608219), 15: (28608128, 28608024), 16: (28602425, 28602315), 17: (28601378, 28601225), 18: (28599080, 28598998),\
                  19: (28597614, 28597487), 20: (28592726, 28592604), 21: (28589838, 28589727), 22: (28589393, 28589294), 23: (28588694, 28588589), 24: (28578311, 28578189)}
    FLT3_exons_targeted = [3, 6, 8, 11, 12, 13, 14, 15, 16, 17, 20]
    FLT3_depths, FLT3_mean_coverage = gene_depths(in_bam, panel_positions, FLT3_exons, '13', '-', 'FLT3')

    #GATA2
    print('calculating read depths across GATA2...')
    GATA2_exons = {3: (128205964, 128205646), 4: (128205211, 128204570), 5: (128202848, 128202703), 6: (128200787, 128200662), 7: (128200161, 128199862)}
    GATA2_exons_targeted = [5, 6, 7]
    GATA2_depths, GATA2_mean_coverage = gene_depths(in_bam, panel_positions, GATA2_exons, '3', '-', 'GATA2')

    #GNAS
    print('calculating read depths across GNAS...')
    GNAS_exons = {1: (57427769, 57430388), 2: (57470667, 57470739), 3: (57473996, 57474040), 4: (57478586, 57478640),\
                  5: (57478727, 57478846), 6: (57480438, 57480535), 7: (57484217, 57484271), 8: (57484405, 57484478),\
                  9: (57484576, 57484634), 10: (57484739, 57484859), 11: (57485006, 57485136), 12: ( 57485389, 57485456),\
                  13: (57485738, 57486247)}
    GNAS_exons_targeted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    GNAS_depths, GNAS_mean_coverage = gene_depths(in_bam, panel_positions, GNAS_exons, '20', '+', 'GNAS')

    #GNB1
    print('calculating read depths across GNB1...')
    GNB1_exons = {1: (1822495, 1822259), 2: (1770677, 1770629), 3: (1756938, 1756836), 4: (1749314, 1749276), 5: (1747301, 1747195), \
                  6: (1737977, 1737914), 7: (1736020, 1735858), 8: (1724750, 1724684), 9: (1722035, 1721834), 10: (1720708, 1720492),\
                  11: (1718876, 1718761), 12: (1718492, 1716729)}
    GNB1_exons_targeted = [5]
    GNB1_depths, GNB1_mean_coverage = gene_depths(in_bam, panel_positions, GNB1_exons, '1', '-', 'GNB1')

    #IDH1
    print('calculating read depths across IDH1...')
    IDH1_exons = {3: (209116275, 209116154), 4: (209113384, 209113093), 5: (209110148, 209110043), 6: (209108328, 209108151),\
                 7: (209106869, 209106718), 8: (209104727, 209104587), 9: (209103957, 209103795), 10: (209101893,209101803)}
    IDH1_exons_targeted = [4, 6]
    IDH1_depths, IDH1_mean_coverage = gene_depths(in_bam, panel_positions, IDH1_exons, '2', '-', 'IDH1')

    #IDH2
    print('calculating read depths across IDH2...')
    IDH2_exons = {1: (90645622, 90645508), 2: (90634876, 90634785), 3: (90633876, 90633711), 4: (90631979, 90631819), 5: (90631734, 90631591),\
                  6: (90630807, 90630671), 7: (90630495, 90630344), 8: (90628619, 90628507), 9: (90628330, 90628233), 10: (90628140, 90628048),\
                  11: (90627585, 90627494)}
    IDH2_exons_targeted = [4, 8]
    IDH2_depths, IDH2_mean_coverage = gene_depths(in_bam, panel_positions, IDH2_exons, '15', '-', 'IDH2')

    #JAK2
    print('calculating read depths across JAK2...')
    JAK2_exons = {3: (5021988, 5022213), 4: (5029783, 5029906), 5: (5044403, 5044520), 6: (5050686, 5050831), 7: (5054563, 5054884), 8: (5055669, 5055788),\
              9: (5064883, 5065040), 10: (5066678, 5066789), 11: (5069022, 5069208), 12: (5069925, 5070052), 13: (5072492, 5072626),\
              14: (5073698, 5073785), 15: (5077453, 5077580), 16: (5078306, 5078444), 17: (5080229, 5080380), 18: (5080533, 5080683), \
              19: (5081725, 5081861), 20: (5089674, 5089863), 21: (5090446, 5090570), 22: (5090739, 5090911), 23: (5123004, 5123121),\
              24: (5126333, 5126446), 25: (5126684, 5126793)}
    JAK2_exons_targeted = [6, 12, 14]
    JAK2_depths, JAK2_mean_coverage = gene_depths(in_bam, panel_positions, JAK2_exons, '9', '+', 'JAK2')

    #KIT
    print('calculating read depths across KIT...')
    KIT_exons = {1: (55524085+97, 55524248), 2: (55561678, 55561947), 3: (55564450, 55564731), 4: (55565796, 55565932), 5: (55569890, 55570058),\
                 6: (55573264, 55573453), 7: (55575590, 55575705), 8: (55589750, 55589864), 9: (55592023, 55592216), 10: (55593384, 55593490),\
                 11: (55593582, 55593708), 12: (55593989, 55594093), 13: (55594177, 55594287), 14: (55595501, 55595651),\
                 15: (55597494, 55597585), 16: (55598037, 55598164), 17: (55599236, 55599358), 18: (55602664, 55602775),\
                 19: (55602887, 55602986), 20: (55603341, 55603446), 21: (55604595, 55606881)}
    KIT_exons_targeted = [1, 2, 7, 8, 9, 10, 11, 12, 13, 16, 17]
    KIT_depths, KIT_mean_coverage = gene_depths(in_bam, panel_positions, KIT_exons, '4', '+', 'KIT')

    #KRAS
    print('calculating read depths across KRAS...')
    KRAS_exons= {2: (25398318, 25398208), 3: (25380346, 25380168), 4: (25378707, 25378548), 5: (25368494, 25368371)}
    KRAS_exons_targeted = [2, 3, 4, 5]
    KRAS_depths, KRAS_mean_coverage = gene_depths(in_bam, panel_positions, KRAS_exons, '12', '-', 'KRAS')

    #MPL
    print('calculating read depths across MPL...')
    MPL_exons = {1: (43803478, 43803598), 2: (43803770, 43803902), 3: (43804213, 43804391), 4: (43804942, 43805240), 5: (43805635, 43805797),\
                 6: (43806058, 43806184), 7: (43812116, 43812300), 8: (43812463, 43812605), 9: (43814514, 43814673), 10: (43814934, 43815030),\
                 11: (43817887, 43817974), 12: (43818189, 43818443)}
    MPL_exons_targeted = [9, 10, 11, 12]
    MPL_depths, MPL_mean_coverage = gene_depths(in_bam, panel_positions, MPL_exons, '1', '+', 'MPL')

    #NRAS
    print('calculating read depths across NRAS...')
    NRAS_exons = {2: (115258780, 115258671), 3: (115256599, 115256421), 4: (115252349, 115252190), 5: (115251275, 115251156)}
    NRAS_exons_targeted = [2, 3]
    NRAS_depths, NRAS_mean_coverage = gene_depths(in_bam, panel_positions, NRAS_exons, '1', '-', 'NRAS')

    #NPM1
    print('calculating read depths across NPM1...')
    NPM1_exons = {1: (170814953, 170815010), 2: (170817055, 170817134), 3: (170818309, 170818428), 4: (170818710, 170818803), 5: (170819714, 170819820),\
                  6: (170819918, 170819982), 7: (170827157, 170827214), 8: (170827843, 170827929), 9: (170832306, 170832407), 10: (170834704, 170834778),\
                  11: (170837531, 170837569)}
    NPM1_exons_targeted = [11]
    NPM1_depths, NPM1_mean_coverage = gene_depths(in_bam, panel_positions, NPM1_exons, '5', '+', 'NPM1')

    #PPM1D
    print('calculating read depths across PPM1D...')
    PPM1D_exons = {1: (58677776, 58678247), 2: (58700882, 58701110), 3: (58711214, 58711338), 4: (58725253, 58725443), 5: (58733960, 58734202), 6: (58740356, 58740914)}
    PPM1D_exons_targeted = [1, 5, 6]
    PPM1D_depths, PPM1D_mean_coverage = gene_depths(in_bam, panel_positions, PPM1D_exons, '17', '+', 'PPM1D')

    #PTPN11
    print('calculating read depths across PTPN11...')
    PTPN11_exons = {1: (112856916, 112856929), 2: (112884080, 112884202), 3: (112888122, 112888316), 4: (112890999, 112891191), 5: (112892368, 112892484),\
                    6: (112893754, 112893867), 7: (112910748, 112910844), 8: (112915455, 112915534), 9: (112915661, 112915819), 10: (112919878, 112920009),\
                    11: (112924279, 112924433), 12:  (112926247, 112926314), 13: (112926828, 112926979), 14: (112939948, 112940060), 15: (112942499, 112942568)}
    PTPN11_exons_targeted = [3, 7, 8, 13]
    PTPN11_depths, PTPN11_mean_coverage = gene_depths(in_bam, panel_positions, PTPN11_exons, '12', '+', 'PTPN11')

    #RUNX1
    print('calculating read depths across RUNX1...')
    RUNX1_exons = {1: (36259409, 36259140), 2: (36253010, 36252854), 3: (36231875, 36231771), 4: (36206898, 36206707), 5:  (36171759, 36171598), 6: (36164907, 36164432)}
    RUNX1_exons_targeted = [1, 2, 3, 4, 5, 6]
    RUNX1_depths, RUNX1_mean_coverage = gene_depths(in_bam, panel_positions, RUNX1_exons, '21', '-', 'RUNX1')

    #SF3B1
    print('calculating read depths across SF3B1...')
    SF3B1_exons = {1: (198299724, 198299696), 2: (198288698, 198288532), 3: (198285857, 198285753), 4: (198285266, 198285152), 5: (198283312, 198283233),\
                   6: (198281635, 198281465), 7: (198274731, 198274494), 8: (198273305, 198273093), 9: (198272843, 198272722), 10: (198270196, 198269999),\
                   11: (198269901, 198269800), 12: (198268488, 198268309), 13: (198267759, 198267673), 14: (198267550, 198267280), 15: (198266854, 198266709),\
                   16: (198266612, 198266466), 17: (198266249, 198266124), 18: (198265660, 198265439), 19: (198265158, 198264976), 20: (198264890, 198264779),\
                   21: (198263305, 198263185), 22: (198262840, 198262709), 23: (198261052, 198260780), 24: (198257912, 198257696), 25: (198257185, 198257027)}
    SF3B1_exons_targeted = [3, 4, 5, 6, 13, 14, 15, 16, 18, 24]
    SF3B1_depths, SF3B1_mean_coverage = gene_depths(in_bam, panel_positions, SF3B1_exons, '2', '-', 'SF3B1')

    #SRSF2
    print('calculating read depths across SRSF2...')
    SRSF2_exons = {1: (74733242, 74732881), 2: (74732546, 74732242)}
    SRSF2_exons_targeted = [1, 2]
    SRSF2_depths, SRSF2_mean_coverage = gene_depths(in_bam, panel_positions, SRSF2_exons, '17', '-', 'SRSF2')

    #STAG2
    print('calculating read depths across STAG2...')
    STAG2_exons = {3: (123156478, 123156521), 4: (123159690, 123159768), 5: (123164811, 123164975), 6: (123171377, 123171473), 7: (123176419, 123176495),\
                   8: (123179014, 123179218), 9: (123181204, 123181355), 10: (123182855, 123182928), 11: (123184036, 123184159), 12: (123184971, 123185069),\
                   13: (123185165, 123185244), 14: (123189978, 123190085), 15: (123191716, 123191827), 16:  (123195074, 123195191), 17: (123195621, 123195724),\
                   18: (123196752, 123196844), 19: (123196966, 123197055), 20: (123197698, 123197901), 21: (123199726, 123199796), 22: (123200025, 123200112),\
                   23: (123200206, 123200286), 24: (123202414, 123202506), 25: (123204999, 123205173), 26: (123210182, 123210321), 27:  (123211807, 123211908),\
                   28: (123215230, 123215378), 29: (123217271, 123217399), 30: (123220397, 123220620), 31: (123224425, 123224614), 32: (123227868, 123227994),\
                   33: (123229222, 123229299), 34: (123234424, 123234448)}
    STAG2_exons_targeted = list(np.linspace(3, 34, ((34-3)+1)))
    STAG2_depths, STAG2_mean_coverage = gene_depths(in_bam, panel_positions, STAG2_exons, 'X', '+', 'STAG2')

    #TET2
    print('calculating read depths across TET2...')
    TET2_exons= {3: (106155099, 106158508), 4: (106162495, 106162587), 5: (106163990, 106164084), 6: (106164726, 106164935),\
              7: (106180775, 106180926), 8: (106182915, 106183005), 9: (106190766, 106190904), 10: (106193720, 106194075), \
             11: (106196204, 106197676)}
    TET2_exons_targeted = list(np.linspace(3, 11, (11-3)+1))
    TET2_depths, TET2_mean_coverage = gene_depths(in_bam, panel_positions, TET2_exons, '4', '+', 'TET2')

    #TP53
    print('calculating read depths across TP53...')
    TP53_exons = {2: (7579912, 7579839), 3: (7579721, 7579700), 4: (7579590, 7579312), 5: (7578554, 7578371), 6: (7578289, 7578177),\
              7:  (7577608, 7577499), 8: (7577155, 7577019), 9: (7576926, 7576853), 10: (7574033, 7573927), 11: (7573008, 7572926)}
    TP53_exons_targeted = list(np.linspace(2, 11, (11-2)+1))
    TP53_depths, TP53_mean_coverage = gene_depths(in_bam, panel_positions, TP53_exons, '17', '-', 'TP53')

    #RAD21
    print('calculating read depths across RAD21...')
    RAD21_exons = {2: (117879000-32, 117878825), 3: (117875498, 117875369), 4: (117874179, 117874080), 5: (117870697, 117870591), 6: (117869712, 117869506),\
                   7: (117869010, 117868885), 8: (117868527, 117868405), 9: (117866707, 117866484), 10: (117864947, 117864788), 11: (117864335, 117864187),\
                   12: (117863006, 117862857), 13: (117861268, 117861185), 14: (117859930, 117858174+1565)}
    RAD21_exons_targeted = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    RAD21_depths, RAD21_mean_coverage = gene_depths(in_bam, panel_positions, RAD21_exons, '8', '-', 'RAD21')

    #U2AF1
    print('calculating read depths across U2AF1...')
    U2AF1_exons = {1: (44527604, 44527561), 2: (44524512, 44524425), 3: (44520629, 44520563), 4: (44515853, 44515804), 5: (44515646, 44515548),\
               6: (44514898, 44514765), 7: (44514673, 44514581), 8: (44513359, 44513211)}
    U2AF1_exons_targeted = [2, 6]
    U2AF1_depths, U2AF1_mean_coverage = gene_depths(in_bam, panel_positions, U2AF1_exons, '21', '-', 'U2AF1')

    #WT1
    print('calculating read depths across WT1...')
    WT1_exons = {1: (32457176-285, 32456246), 2: (32450165, 32450043), 3: (32449604, 32449502), 4: (32439200, 32439123), 5: (32438086, 32438036),\
                 6: (32421590, 32421494), 7: (32417953, 32417803), 8: (32414301, 32414212), 9:  (32413610, 32413518), 10: (32410725, 32409321+1278)}
    WT1_exons_targeted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    WT1_depths, WT1_mean_coverage = gene_depths(in_bam, panel_positions, WT1_exons, '11', '-', 'WT1')

    #ZRSR2
    print('calculating read depths across ZRSR2...')
    ZRSR2_exons = {1: (15808595+23, 15808659), 2: (15809057, 15809136), 3: (15817995, 15818076), 4: (15821811, 15821919),\
                   5: (15822234, 15822320), 6: (15826356, 15826394), 7: (15827323, 15827441), 8: (15833800, 15834013),\
                   9: (15836710, 15836765), 10: (15838330, 15838439), 11: (15840854, 15841383)}
    ZRSR2_exons_targeted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    ZRSR2_depths, ZRSR2_mean_coverage = gene_depths(in_bam, panel_positions, ZRSR2_exons, 'X', '+', 'ZRSR2')

    #SNPs
    print('calculating read depths across 4 SNPs...')
    rs10789158_depth, rs10789158_mean = SNP_depths(in_bam, panel_positions, '1', 63936248, 'rs10789158')
    rs3916765_depth, rs3916765_mean = SNP_depths(in_bam, panel_positions, '6', 32749936+60, 'rs3916765')
    rs1364429_depth, rs1364429_mean = SNP_depths(in_bam, panel_positions, '7', 134714538+60, 'rs1364429')
    rs2286510_depth, rs2286510_mean = SNP_depths(in_bam, panel_positions, '17', 9259404+60, 'rs2286510')

    ##### CREATE A MERGED DICTIONARY OF THE READ DEPTHS, WRITE TO FILE AND CALCULATE THE MAX READ DEPTH (TO SCALE THE AXES OF THE PLOTS) AND MEAN DEPTH #######
    all_read_depths = {**ASXL1_depths, **BCOR_depths, **BCORL1_depths, **CBL_depths, **CEBPA_depths,
                   **CHEK2_depths, **CSF3R_depths, **DDX41_depths, **DNMT3A_depths, **EZH2_depths,
                   **FLT3_depths, **GATA2_depths, **GNAS_depths, **GNB1_depths, **IDH1_depths, **IDH2_depths,
                   **JAK2_depths, **KIT_depths, **KRAS_depths, **MPL_depths, **NPM1_depths, **NRAS_depths,
                   **PPM1D_depths, **PTPN11_depths, **RAD21_depths, **RUNX1_depths, **SF3B1_depths, **SRSF2_depths,
                   **STAG2_depths, **TET2_depths,  **TP53_depths, **U2AF1_depths, **WT1_depths, **ZRSR2_depths,
                   **rs10789158_depth, **rs3916765_depth, **rs1364429_depth, **rs2286510_depth}

    ##### CALCULATE THE MAXIMUM Y TO PLOT ON THE Y AXIS OF ALL THE PLOTS (SO THEY HAVE THE SAME AXES) ######
    def roundup(x, rounding_number): #function to round up a number (to be used for y-axis maximum, e.g. round to nearest 50,000)
        return int(math.ceil(x/rounding_number))*rounding_number

    max_depth = 0
    total = 0
    for depth in all_read_depths.values():
        total+=depth
        if depth > max_depth:
            max_depth = depth #replacing the max depth with the new max depth

    number_positions = len(all_read_depths)
    mean_depth_panel = total/number_positions #calculate the mean depth across the whole panel (targeted regions)

    if max_depth >= 100000:
        roundup_number = 50000
    if 10000 <= max_depth < 100000:
        roundup_number = 15000
    if 1000 <= max_depth < 10000:
        roundup_number = 1500
    if 100 <= max_depth < 1000:
        roundup_number = 200

    y_max = roundup(max_depth, roundup_number) #calculate the maximum for the y-axis on the plots

    ###### CREATE TEXT FILE CONTAINING READ DEPTHS AND MEAN READ DEPTH ACROSS PANEL #######
    bam_type_name = bam_type.replace(' ', '_')
    read_depth_metrics = open(out_directory+'/Gene_read_depths/'+sample_name+'_watson_code_read_depth_metrics_'+bam_type_name+'.txt', 'w')
    read_depth_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    read_depth_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    read_depth_metrics.write('produced from code:\t' + 'Watson_code_SNV_panel_read_coverage: version '+str(version) + '\n')
    read_depth_metrics.write('type of BAM file analysed:\t' + bam_type + '\n')
    read_depth_metrics.write('\n')
    read_depth_metrics.write('mean read depth across panel:\t'+ str(mean_depth_panel)+ '\n')
    read_depth_metrics.write('\n')
    read_depth_metrics.write('chromosome\t'+ 'position\t' + 'gene\t' + 'exon\t'+ 'read depth'+'\n')

    for k, v in all_read_depths.items():
        chromosome = 'chr'+str(k[0])
        position = k[1]
        gene = k[2]
        exon = k[3]
        depth = v
        read_depth_metrics.write(chromosome+'\t' + str(position)+'\t' + gene+'\t' + exon+'\t' + str(depth)+'\n')

    read_depth_metrics.close()

    ###### PLOT READ DEPTHS ########
    print('plotting read depth plots...')
    ASXL1_plot = gene_coverage_depth_plot(in_bam, ASXL1_depths, ASXL1_mean_coverage, ASXL1_exons, '20', '+', 'ASXL1', bam_type, sample_name, ASXL1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    BCOR_plot = gene_coverage_depth_plot(in_bam, BCOR_depths, BCOR_mean_coverage, BCOR_exons, 'X', '-', 'BCOR', bam_type, sample_name, BCOR_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    BCORL1_plot = gene_coverage_depth_plot(in_bam, BCORL1_depths, BCORL1_mean_coverage, BCORL1_exons, 'X', '+', 'BCORL1', bam_type, sample_name, BCORL1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    CBL_plot = gene_coverage_depth_plot(in_bam, CBL_depths, CBL_mean_coverage, CBL_exons, '11', '+', 'CBL', bam_type, sample_name, CBL_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    CEBPA_plot = gene_coverage_depth_plot(in_bam, CEBPA_depths, CEBPA_mean_coverage, CEBPA_exons, '19', '-', 'CEBPA', bam_type, sample_name, CEBPA_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    CHEK2_plot = gene_coverage_depth_plot(in_bam, CHEK2_depths, CHEK2_mean_coverage, CHEK2_exons, '22', '-', 'CHEK2', bam_type, sample_name, CHEK2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    CSF3R_plot = gene_coverage_depth_plot(in_bam, CSF3R_depths, CSF3R_mean_coverage, CSF3R_exons, '1', '-', 'CSF3R', bam_type, sample_name, CSF3R_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    DDX41_plot = gene_coverage_depth_plot(in_bam, DDX41_depths, DDX41_mean_coverage, DDX41_exons, '5', '-', 'DDX41', bam_type, sample_name, DDX41_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    DNMT3A_plot = gene_coverage_depth_plot(in_bam, DNMT3A_depths, DNMT3A_mean_coverage, DNMT3A_exons, '2', '-', 'DNMT3A', bam_type, sample_name, DNMT3A_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    EZH2_plot = gene_coverage_depth_plot(in_bam, EZH2_depths, EZH2_mean_coverage, EZH2_exons, '7', '-', 'EZH2', bam_type, sample_name, EZH2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    FLT3_plot = gene_coverage_depth_plot(in_bam, FLT3_depths, FLT3_mean_coverage, FLT3_exons, '13', '-', 'FLT3', bam_type, sample_name, FLT3_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    GATA2_plot = gene_coverage_depth_plot(in_bam, GATA2_depths, GATA2_mean_coverage, GATA2_exons, '3', '-', 'GATA2', bam_type, sample_name, GATA2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    GNB1_plot = gene_coverage_depth_plot(in_bam, GNB1_depths, GNB1_mean_coverage, GNB1_exons, '1', '-', 'GNB1', bam_type, sample_name, GNB1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    GNAS_plot = gene_coverage_depth_plot(in_bam, GNAS_depths, GNAS_mean_coverage, GNAS_exons, '20', '+', 'GNAS', bam_type, sample_name, GNAS_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    IDH1_plot = gene_coverage_depth_plot(in_bam, IDH1_depths, IDH1_mean_coverage, IDH1_exons, '2', '-', 'IDH1', bam_type, sample_name, IDH1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    IDH2_plot = gene_coverage_depth_plot(in_bam, IDH2_depths, IDH2_mean_coverage, IDH2_exons, '15', '-', 'IDH2', bam_type, sample_name, IDH2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    JAK2_plot = gene_coverage_depth_plot(in_bam, JAK2_depths, JAK2_mean_coverage, JAK2_exons, '9', '+', 'JAK2', bam_type, sample_name, JAK2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    KRAS_plot = gene_coverage_depth_plot(in_bam, KRAS_depths, KRAS_mean_coverage, KRAS_exons, '12', '-', 'KRAS', bam_type, sample_name, KRAS_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    MPL_plot = gene_coverage_depth_plot(in_bam, MPL_depths, MPL_mean_coverage, MPL_exons, '1', '+', 'MPL', bam_type, sample_name, MPL_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    NRAS_plot = gene_coverage_depth_plot(in_bam, NRAS_depths, NRAS_mean_coverage, NRAS_exons, '1', '-', 'NRAS', bam_type, sample_name, NRAS_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    NPM1_plot = gene_coverage_depth_plot(in_bam, NPM1_depths, NPM1_mean_coverage, NPM1_exons, '5', '+', 'NPM1', bam_type, sample_name, NPM1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    RAD21_plot = gene_coverage_depth_plot(in_bam, RAD21_depths, RAD21_mean_coverage, RAD21_exons, '8', '-', 'RAD21', bam_type, sample_name, RAD21_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    RUNX1_plot = gene_coverage_depth_plot(in_bam, RUNX1_depths, RUNX1_mean_coverage, RUNX1_exons, '21', '-', 'RUNX1', bam_type, sample_name, RUNX1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    PPM1D_plot = gene_coverage_depth_plot(in_bam, PPM1D_depths, PPM1D_mean_coverage, PPM1D_exons, '17', '+', 'PPM1D', bam_type, sample_name, PPM1D_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    PTPN11_plot = gene_coverage_depth_plot(in_bam, PTPN11_depths, PTPN11_mean_coverage, PTPN11_exons, '12', '+', 'PTPN11', bam_type, sample_name, PTPN11_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    SF3B1_plot = gene_coverage_depth_plot(in_bam, SF3B1_depths, SF3B1_mean_coverage, SF3B1_exons, '2', '-', 'SF3B1', bam_type, sample_name, SF3B1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    SRSF2_plot = gene_coverage_depth_plot(in_bam, SRSF2_depths, SRSF2_mean_coverage, SRSF2_exons, '17', '-', 'SRSF2', bam_type, sample_name, SRSF2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    STAG2_plot = gene_coverage_depth_plot(in_bam, STAG2_depths, STAG2_mean_coverage, STAG2_exons, 'X', '+', 'STAG2', bam_type, sample_name, STAG2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    TET2_plot = gene_coverage_depth_plot(in_bam, TET2_depths, TET2_mean_coverage, TET2_exons, '4', '+', 'TET2', bam_type, sample_name, TET2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    TP53_plot = gene_coverage_depth_plot(in_bam, TP53_depths, TP53_mean_coverage, TP53_exons, '17', '-', 'TP53', bam_type, sample_name, TP53_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    U2AF1_plot = gene_coverage_depth_plot(in_bam, U2AF1_depths, U2AF1_mean_coverage, U2AF1_exons, '21', '-', 'U2AF1', bam_type, sample_name, U2AF1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    WT1_plot = gene_coverage_depth_plot(in_bam, WT1_depths, WT1_mean_coverage, WT1_exons, '11', '-', 'WT1', bam_type, sample_name, WT1_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    ZRSR2_plot = gene_coverage_depth_plot(in_bam, ZRSR2_depths, ZRSR2_mean_coverage, ZRSR2_exons, 'X', '+', 'ZRSR2', bam_type, sample_name, ZRSR2_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    KIT_plot = gene_coverage_depth_plot(in_bam, KIT_depths, KIT_mean_coverage, KIT_exons, '4', '+', 'KIT', bam_type, sample_name, KIT_exons_targeted, y_max, mean_depth_panel, out_directory, min_family_size)
    rs10789158_plot = SNP_coverage_depth_plot(in_bam, rs10789158_depth, rs10789158_mean, '1', 63936248, 'rs10789158', bam_type, sample_name, y_max, mean_depth_panel, out_directory, min_family_size)
    rs3916765_plot = SNP_coverage_depth_plot(in_bam, rs3916765_depth, rs3916765_mean, '6', 32749996, 'rs3916765', bam_type, sample_name, y_max, mean_depth_panel, out_directory, min_family_size)
    rs1364429_plot = SNP_coverage_depth_plot(in_bam, rs1364429_depth, rs1364429_mean, '7', 134714598, 'rs1364429', bam_type, sample_name, y_max, mean_depth_panel, out_directory, min_family_size)
    rs2286510_plot = SNP_coverage_depth_plot(in_bam, rs2286510_depth, rs2286510_mean, '17', 9259464, 'rs2286510', bam_type, sample_name, y_max, mean_depth_panel, out_directory, min_family_size)

    return print('Calculation and plotting of read depths complete')

if __name__ == "__main__":
	main()
