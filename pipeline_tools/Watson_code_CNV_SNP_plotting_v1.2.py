#!/usr/bin/env python

'''''
Watson code for creating a VCF file from samtools mpileup and for creating a txt file containing all variants and SNP sites in 1000 genomes covered by panel.
Version 1.2 (May 2021)

Input:
    1) samtools mpileup txt file
    2) sample name

Outputs:
    1) VCF file containing all positions in the targeted panel
    2) Text file containing all variant positions and all SNP positions from 1000 genomes covered by panel

Usage:
Watson_code_VCF_SNP_plotting_v1.2.py --targeted_SNPs file containing the SNPs specifically targeted by panel
                                    --sample-name sample_name
                                    --annovar_exonic annovar.exonic_variant_function file
                                    --annovar_variant annovar.variant_function file
                                    --annovar_invalid annovar.invalid_input file
                                    --panel_bed input bed file for the panel
                                    --chromosome_ideogram chromosome ideogram file
                                    --out-directory directory to save files in

'''''
version = '1.2'

from argparse import ArgumentParser
import pysam
import sys
import gzip
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker #plticker
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MaxNLocator, MultipleLocator
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
from array import array
import pandas as pd
import timeit
import time
import shelve
from datetime import date
from Bio.Seq import Seq
from ast import literal_eval
import csv

today = date.today()
date_today = today.strftime("%d/%m/%Y")
VCF_date_today = today.strftime("%Y%m%d")

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
blue1 = '#eff3ff'
blue2 = '#bdd7e7'
blue3 = '#6baed6'
blue4 = '#3182bd'
blue5 = '#08519c'
green1 = '#edf8e9'
green2 = '#bae4b3'
green3 = '#74c476'
green4 = '#31a354'
green5 = '#006d2c'
grey1 = '#f7f7f7'
grey2 = '#cccccc'
grey3 = '#969696'
grey4 = '#636363'
grey5 = '#252525'
purple1 = '#f2f0f7'
purple2 = '#cbc9e2'
purple3 = '#9e9ac8'
purple4 = '#756bb1'
purple5 = '#54278f'
red1 = '#fee5d9'
red2 = '#fcae91'
red3 = '#fb6a4a'
red4 = '#de2d26'
red5 = '#a50f15'


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

def plot_regions_targeted_by_panel(chromosome, panel_coverage, ax, y_max, chromosome_sizes):

    regions_covered = panel_coverage[chromosome]
    chromosome_size = chromosome_sizes[chromosome]

    for positions in regions_covered:
        start_position = positions[0]
        stop_position = positions[1]

        bottom = 0
        top = y_max

        x = [start_position, start_position, stop_position, stop_position]
        y = [bottom, top, top, bottom]
        ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with colour

        if stop_position == chromosome_size: #put a line at the end of the plot of final exon targeted
            ax.plot([chromosome_size-1, chromosome_size-1], [0, 800000], color = grey1, lw = 4, zorder = 6) #line at end of plot

    return ax

def plot_regions_targeted_by_panel_cumulative(cumulative_chromosome_sizes, panel_coverage, ax):
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        cumulative_chromosome_start = cumulative_chromosome_sizes['chr'+i]
        regions_covered = panel_coverage[i]

        for positions in regions_covered:
            start_position = positions[0]+cumulative_chromosome_start
            stop_position = positions[1]+cumulative_chromosome_start

            bottom = 0
            top = 1.02

            x = [start_position, start_position, stop_position, stop_position]
            y = [bottom, top, top, bottom]
            ax.fill(x, y, color= '#deebf7', fill = True, alpha = 1.0, linewidth = 1, zorder = 0) #fill in the box with stars

    return ax

def panel_coverage_dict(panel_bed):
    with open(panel_bed, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        panel_coverage = {}
        for row in read_reader:
            if row_count > 2:
                chromosome = row[0].replace('chr', '') #change 'chr2' to '2' etc.. (keep as string, because cannot convert X and Y to integers)
                start = int(row[1]) #start position of that part of the chromosome
                end = int(row[2]) #end position of that part of the chromosome
                if chromosome in panel_coverage.keys():
                    panel_coverage[chromosome].append((start, end))
                else:
                    panel_coverage[chromosome]=([(start, end)])

            row_count+=1
    return panel_coverage

def bases_covered_by_panel(chromosome, panel_coverage):
    regions_covered = panel_coverage[chromosome]
    bases = []
    for i in regions_covered:
        start = i[0]
        end = i[1]
        for base in np.linspace(start, end, (end-start)+1):
            bases.append(base)

    return bases

def targeted_SNPs(SNP_file):
    with open(SNP_file, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = ',')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        SNPs_targeted_dict = {}
        SNPs_targeted = []
        for row in read_reader:
            if row_count > 0:
                chromosome = row[1]
                position = row[2]
                ref = row[4]
                maj = row[5]
                minor = row[6]
                SNPs_targeted_dict[(chromosome, int(position), ref)] = (maj, minor)
                SNPs_targeted.append((chromosome, int(position), ref))

            row_count+=1
    return SNPs_targeted_dict

def merge_annovar_files(annovar_exonic, annovar_variant, annovar_invalid, sample_name, out_directory):
    exonic_df = pd.read_csv(annovar_exonic, delimiter = '\t',
                            names=["line", "syn_nonsyn", "exonic_annotation", "chromosome", "position", "end_position", \
                                   "REF", "ALT", "RSID", "MAF_1000_genomes", "total_depth", "variant_depth", "VAF", "variant_type", \
                                   "REF_mean_UMI_family_size", "variant_mean_UMI_family_size", "mean_position_in_read", "std_position_in_read",\
                                  "mean_base_quality", "std_base_quality", "mean_mapq", "REF_depth_by_strand",\
                                  "variant_depth_by_strand", "strand_bias_fisher_p_value", "strand_bias_fisher_p-value_phred"])


    #reorder the colums in the dataframe
    exonic_df2 = exonic_df[["chromosome", "position", "end_position", \
                                   "REF", "ALT", "RSID", "MAF_1000_genomes", "total_depth", "variant_depth", "VAF", "variant_type", \
                                   "REF_mean_UMI_family_size", "variant_mean_UMI_family_size", "mean_position_in_read", "std_position_in_read",\
                                  "mean_base_quality", "std_base_quality", "mean_mapq", "REF_depth_by_strand",\
                                  "variant_depth_by_strand", "strand_bias_fisher_p_value", "strand_bias_fisher_p-value_phred", "syn_nonsyn", "exonic_annotation"]]

    variant_df = pd.read_csv(annovar_variant, delimiter = '\t',
                            names=["intronic_exonic", "gene", "chromosome", "position", "end_position", \
                                   "REF", "ALT", "RSID", "MAF_1000_genomes", "total_depth", "variant_depth", "VAF", "variant_type", \
                                   "REF_mean_UMI_family_size", "variant_mean_UMI_family_size", "mean_position_in_read", "std_position_in_read",\
                                  "mean_base_quality", "std_base_quality", "mean_mapq", "REF_depth_by_strand",\
                                  "variant_depth_by_strand", "strand_bias_fisher_p_value", "strand_bias_fisher_p-value_phred"])

    variant_df2 = variant_df[["chromosome", "position", "end_position", \
                                   "REF", "ALT", "RSID", "MAF_1000_genomes", "total_depth", "variant_depth", "VAF", "variant_type", \
                                   "REF_mean_UMI_family_size", "variant_mean_UMI_family_size", "mean_position_in_read", "std_position_in_read",\
                                  "mean_base_quality", "std_base_quality", "mean_mapq", "REF_depth_by_strand",\
                                  "variant_depth_by_strand", "strand_bias_fisher_p_value", "strand_bias_fisher_p-value_phred", "intronic_exonic", "gene"]]

    result = pd.merge(variant_df2, exonic_df2, how = 'outer', on=["chromosome", "position", "end_position", \
                                   "REF", "ALT", "RSID", "MAF_1000_genomes", "total_depth", "variant_depth", "VAF", "variant_type", \
                                   "REF_mean_UMI_family_size", "variant_mean_UMI_family_size", "mean_position_in_read", "std_position_in_read",\
                                  "mean_base_quality", "std_base_quality", "mean_mapq", "REF_depth_by_strand",\
                                  "variant_depth_by_strand", "strand_bias_fisher_p_value", "strand_bias_fisher_p-value_phred"])

    result.to_csv(out_directory+'/'+sample_name+'_watson_code_just_SNPs_result_df.txt', index=False, sep = '\t')

    complex_df = pd.read_csv(annovar_invalid, delimiter = '\t',
                            names=["chromosome", "position", "end_position", \
                                   "REF", "ALT", "RSID", "MAF_1000_genomes", "total_depth", "variant_depth", "VAF", "variant_type", \
                                   "REF_mean_UMI_family_size", "variant_mean_UMI_family_size", "mean_position_in_read", "std_position_in_read",\
                                  "mean_base_quality", "std_base_quality", "mean_mapq", "REF_depth_by_strand",\
                                  "variant_depth_by_strand", "strand_bias_fisher_p_value", "strand_bias_fisher_p-value_phred"])

    complex_df = complex_df.iloc[6:]

    result2 = pd.concat([result, complex_df], axis = 0, join = 'outer', ignore_index = True, sort = False)
    result2.fillna('-', inplace = True)

    result2.to_csv(out_directory+'/'+sample_name+'_watson_code_SSCS_variant_calling_variants_and_SNPs_annovar_annotated.txt', index=False, sep = '\t')

    #also create a txt file of just the SNPs
    just_SNPs = result2[result2['RSID']!='-']

    return just_SNPs.to_csv(out_directory+'/'+sample_name+'_watson_code_SSCS_variant_calling_only_SNPs_annovar_annotated.txt', index=False, sep = '\t')

def variants_dictionary(sample_name, out_directory):
    with open(out_directory+'/'+sample_name+'_watson_code_SSCS_variant_calling_only_SNPs_annovar_annotated.txt', 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        variants_dict = {}
        for row in read_reader:
            if row_count > 5:
                chromosome = row[0]
                position = int(row[1])
                end_position = int(row[2])
                RSID = row[5]
                if RSID != '-':
                    MAF_1000_genomes = float(row[6])
                    reference = row[3]
                    alt = row[4]
                    total_depth = int(row[7])
                    variant_depth = int(row[8])
                    VAF = float(row[9])
                    variant_type = row[10]
                    if variant_type == 'SNV' or variant_type == 'REF':
                        variants_dict[(chromosome, position, reference)] = (alt, VAF, RSID, MAF_1000_genomes)

            row_count+=1

    return variants_dict

def chromosome_SNPs_plot(chromosome, variants_dict, panel_coverage, ideogram_file, SNPs_targeted_dict, chromosome_sizes, out_directory, sample_name):
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

    chromosome_size = chromosome_sizes[chromosome]
    bases_targeted = bases_covered_by_panel(chromosome, panel_coverage)

    #SNPs
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    for k, v in variants_dict.items():
        chrom = k[0]
        if chrom == 'chr'+chromosome:
            position = int(k[1])
            VAF = float(v[1])
            RSID = v[2]
            MAF_1000_genomes = v[3]
            if k in SNPs_targeted_dict.keys():
                x1.append(position)
                y1.append(VAF)
            else:
                x2.append(position)
                y2.append(VAF)

    ax1.scatter(x1, y1, color = '#4292c6', zorder = 200, label = 'target SNPs in panel')
    ax1.scatter(x2, y2, color = '#a1d99b', zorder = 100, label = 'other SNPs >1% MAF in 1000 genomes')

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel(str(chromosome), panel_coverage, ax1, 1.02, chromosome_sizes)

    #CHROMOSOME IDEOGRAM
    plot_chromosome(ideogram_file, 'chr'+str(chromosome), ax2)

    ax1.plot([0, chromosome_size], [0.5, 0.5], color = grey3, lw = 2, linestyle = ':')


    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, chromosome_size)
    ax2.set_xlim(0, chromosome_size)
    # ax2.set_ylim(-0.5, 0)

    # Changing the y axis to log scale
    ax1.set_ylim(-0.02, 1.02)

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

    ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

    legend_elements = [Line2D([0], [0], marker = 's', color='#deebf7', alpha=1.0, markersize = 12, \
                      lw=0, label='regions covered by TWIST CNV panel'),\
                      Line2D([0], [0], marker = 'o', color='#4292c6', alpha=1.0, markersize = 8, \
                      lw=0, label='target SNPs in panel'),
                      Line2D([0], [0], marker = 'o', color='#a1d99b', alpha=1.0, markersize = 8, \
                      lw=0, label='other SNPs >1% MAF in 1000 genomes')]

    ax1.legend(ncol=3, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13, frameon=True,\
              fancybox = True)

    #Title and axis labels
    ax1.set_title('BAFs across chr'+chromosome+': '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('BAF', fontsize = axislabelfont)

    plt.tight_layout()

    return plt.savefig(out_directory+'/SNPs_plots/'+sample_name+'_SSCS_SNPs_BAFs_chr'+str(chromosome)+'.pdf')

def all_chromosomes_SNPs_plot(SNPs_targeted_dict, variants_dict, cumulative_chromosome_sizes, panel_coverage, sample_name, chromosome_sizes, out_directory):
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, sharey=False, sharex = True, figsize=(27, 5))

    axisfont=15
    titlefont=15
    axislabelfont=15

    #SNPs
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    all_positions = []
    for k, v in variants_dict.items():
        cumulative_chromosome_start = cumulative_chromosome_sizes[k[0]]
        position = int(k[1])+cumulative_chromosome_start
        VAF = float(v[1])
        RSID = v[2]
        MAF_1000_genomes = v[3]
        all_positions.append(position)
        if k in SNPs_targeted_dict.keys():
            x1.append(position)
            y1.append(VAF)
        else:
            x2.append(position)
            y2.append(VAF)

    ax1.scatter(x1, y1, color = '#4292c6', zorder = 200, label = 'target SNPs in panel')
    ax1.scatter(x2, y2, color = '#a1d99b', zorder = 100, label = 'other SNPs >1% MAF in 1000 genomes')

    earliest_position = min(all_positions)
    latest_position = max(all_positions)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel_cumulative(cumulative_chromosome_sizes, panel_coverage, ax1)

    #PLOT A VERTICAL LINE TO DISTIGUISH THE CHROMOSOMES
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        if left_start !=0:
            ax1.plot([left_start, left_start], [-0.02, 1.02], color = 'k', lw = 1, zorder = 1000)

    end_X = cumulative_chromosome_sizes['chrX']+(chromosome_sizes['X'])
    ax1.plot([end_X, end_X], [-0.02, 1.02], color = 'k', lw = 1, zorder = 1000)

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, cumulative_chromosome_sizes['chrX']+(chromosome_sizes['X'])+100)

    # Changing the y axis limits
    ax1.set_ylim(0, 1.02)

    #x-axis ticks
    x_axis_tick_positions = []
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        chrom = k[3:]
        chromosome_size = chromosome_sizes[chrom]
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
                      Line2D([0], [0], marker = 'o', color='#4292c6', alpha=1.0, markersize = 8, \
                      lw=0, label='target SNPs in panel'),
                      Line2D([0], [0], marker = 'o', color='#a1d99b', alpha=1.0, markersize = 8, \
                      lw=0, label='other SNPs >1% MAF in 1000 genomes')]

    ax1.legend(ncol=3, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13, frameon=True,\
              fancybox = True)

    #Title and axis labels
    ax1.set_title('BAFs (from SSCS): '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('BAF', fontsize = axislabelfont)
    ax1.set_xlabel('chromosome', fontsize = axislabelfont)

    plt.tight_layout()

    return  plt.savefig(out_directory+'/SNPs_plots/'+sample_name+'_SSCS_SNPs_BAFs_all_chromosomes.pdf')

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--targeted_SNPs", action="store", dest="SNP_file", help="csv file detailing all the target SNPs in the panel", required=True)
    parser.add_argument("--annovar_exonic", action="store", dest="annovar_exonic", help="annovar.exonic_variant_function file", required=True)
    parser.add_argument("--annovar_variant", action="store", dest="annovar_variant", help="annovar.variant_function file", required=True)
    parser.add_argument("--annovar_invalid", action="store", dest="annovar_invalid", help="annovar.invalid_input file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--panel_bed", action="store", dest="panel_bed", help="input bed file for the panel", required=True)
    parser.add_argument("--chromosome-ideogram", action="store", dest="chromosome_ideogram", help="input chromosome ideogram file", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    SNP_file = o.SNP_file
    annovar_exonic = o.annovar_exonic
    annovar_variant = o.annovar_variant
    annovar_invalid = o.annovar_invalid
    sample_name = o.sample_name
    panel_bed = o.panel_bed
    ideogram_file = o.chromosome_ideogram
    out_directory = o.out_directory

    merge_annovar_files(annovar_exonic, annovar_variant, annovar_invalid, sample_name, out_directory)

    chromosome_sizes, chromosome_bands = chromosome_sizes_bands(ideogram_file)
    SNPs_targeted_dict = targeted_SNPs(SNP_file)
    panel_coverage = panel_coverage_dict(panel_bed)
    variants_dict = variants_dictionary(sample_name, out_directory)

    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        print('plotting BAFs across chromosome'+str(i))
        chromosome_SNPs_plot(i, variants_dict, panel_coverage, ideogram_file, SNPs_targeted_dict, chromosome_sizes, out_directory, sample_name)

    cumulative_chromosome_sizes = {}
    length = 0
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        chromosome_size = chromosome_sizes[i]
        cumulative_chromosome_sizes['chr'+i]=length
        length+=chromosome_size

    print('plotting BAFs across all chromosomes on one plot')
    all_chromosomes_SNPs_plot(SNPs_targeted_dict, variants_dict, cumulative_chromosome_sizes, panel_coverage, sample_name, chromosome_sizes, out_directory)

    return print('SNP plots complete')

if __name__ == "__main__":
	main()
