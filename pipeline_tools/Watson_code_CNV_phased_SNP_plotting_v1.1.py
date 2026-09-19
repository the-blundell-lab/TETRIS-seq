#!/usr/bin/env python

'''''
Watson code for plotting the phased SNPs (phased using Eagle).
Version 1.1 (December 202O)

Input:
    1) VCF file containing information on phased SNPs
    2) sample name

Outputs:
    1) BAF plot for each chromosome showing phased SNPs

Usage:
Watson_code_VCF_phased_SNP_plotting_v1.1.py --phased_VCF phased VCF file
                                    --sample-name sample_name
                                    --panel_bed input bed file for the panel
                                    --chromosome_ideogram chromosome ideogram file
                                    --out-directory directory to save files in

'''''
version = '1.1'

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
from Bio.Alphabet import IUPAC
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
        cumulative_chromosome_start = cumulative_chromosome_sizes[i]
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

def variants_dictionary(phased_VCF):
    with open(phased_VCF, 'r') as textfile:
        read_reader = csv.reader(textfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0
        variants_dict = {}
        for row in read_reader:
            if len(row)>3: #i.e. not the header
                if row[0]!='#CHROM': #i.e. not the header of the VCF columns
                    chromosome = row[0]
                    position = row[1]
                    RSID = row[2]
                    ref = row[3]
                    alt = row[4]
                    VAF = float(row[7].split(';')[5].split('=')[1])
                    GT = row[9].split(':')[0]
                    if GT[1]=='/':
                        phase = 'unphased'
                    if GT[1]=='|':
                        if GT[0]=='1' and GT[2]=='0':
                            phase = 'allele_1'
                        if GT[2]=='1' and GT[0]=='0':
                            phase = 'allele_2'
                        if GT[0]=='1' and GT[2]=='1':
                            phase = 'homozygous'
                    if RSID != '-':
                        variants_dict[(chromosome, position, ref, alt, RSID)]=(VAF, GT, phase)
            row_count+=1

    return variants_dict

def chromosome_SNPs_plot(chromosome, variants_dict, panel_coverage, ideogram_file, chromosome_sizes, out_directory, sample_name):
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
    x3 = []
    y3 = []
    x4 = []
    y4 = []

    for k, v in variants_dict.items():
        chrom = k[0]
        if chrom == chromosome: #e.g. '1'
            position = int(k[1])
            VAF = float(v[0])
            RSID = k[4]
            allele = v[2]
            if allele == 'allele_1':
                x1.append(position)
                y1.append(VAF)
            if allele == 'allele_2':
                x2.append(position)
                y2.append(VAF)
            if allele == 'homozygous':
                x3.append(position)
                y3.append(VAF)
            if allele == 'unphased':
                x4.append(position)
                y4.append(VAF)

    allele1 = '#4292c6'
    allele2 = '#a1d99b'
    hetero = '#fdae6b'
    unphased = grey2

    ax1.scatter(x4, y4, color = unphased, zorder = 100, label = 'unphased')
    ax1.scatter(x3, y3, color = hetero, zorder = 200, label = 'heterozygous')
    ax1.scatter(x1, y1, color = allele1, zorder = 200, label = 'allele 1')
    ax1.scatter(x2, y2, color = allele2, zorder = 100, label = 'allele 2')

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
                      Line2D([0], [0], marker = 'o', color=allele1, alpha=1.0, markersize = 8, \
                      lw=0, label='phased allele 1'),
                      Line2D([0], [0], marker = 'o', color=allele2, alpha=1.0, markersize = 8, \
                      lw=0, label='phased allele 2'),
                      Line2D([0], [0], marker = 'o', color=hetero, alpha=1.0, markersize = 8, \
                      lw=0, label='heterozygous'),
                      Line2D([0], [0], marker = 'o', color=unphased, alpha=1.0, markersize = 8, \
                      lw=0, label='unphased')]

    ax1.legend(ncol=5, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13, frameon=True,\
              fancybox = True)

    #Title and axis labels
    ax1.set_title('phased BAFs across chr'+chromosome+': '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('BAF', fontsize = axislabelfont)

    plt.tight_layout()

    return plt.savefig(out_directory+'/Phased_SNPs/Phased_SNPs_plots/'+sample_name+'_SSCS_phased_SNPs_BAFs_chr'+str(chromosome)+'.pdf')

def all_chromosomes_SNPs_plot(variants_dict, cumulative_chromosome_sizes, panel_coverage, sample_name, chromosome_sizes, out_directory):
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
    x3 = []
    y3 = []
    x4 = []
    y4 = []
    all_positions = []

    for k, v in variants_dict.items():
        cumulative_chromosome_start = cumulative_chromosome_sizes[k[0]]
        position = int(k[1])+cumulative_chromosome_start
        all_positions.append(position)
        VAF = float(v[0])
        RSID = k[4]
        allele = v[2]
        if allele == 'allele_1':
            x1.append(position)
            y1.append(VAF)
        if allele == 'allele_2':
            x2.append(position)
            y2.append(VAF)
        if allele == 'homozygous':
            x3.append(position)
            y3.append(VAF)
        if allele == 'unphased':
            x4.append(position)
            y4.append(VAF)

    allele1 = '#4292c6'
    allele2 = '#a1d99b'
    hetero = '#fdae6b'
    unphased = grey2

    ax1.scatter(x4, y4, color = unphased, zorder = 100, label = 'unphased')
    ax1.scatter(x3, y3, color = hetero, zorder = 200, label = 'heterozygous')
    ax1.scatter(x1, y1, color = allele1, zorder = 200, label = 'allele 1')
    ax1.scatter(x2, y2, color = allele2, zorder = 100, label = 'allele 2')

    earliest_position = min(all_positions)
    latest_position = max(all_positions)

    #REGIONS COVERED BY PANEL:
    plot_regions_targeted_by_panel_cumulative(cumulative_chromosome_sizes, panel_coverage, ax1)

    #PLOT A VERTICAL LINE TO DISTIGUISH THE CHROMOSOMES
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        if left_start !=0:
            ax1.plot([left_start, left_start], [-0.02, 1.02], color = 'k', lw = 1, zorder = 1000)

    end_X = cumulative_chromosome_sizes['X']+(chromosome_sizes['X'])
    ax1.plot([end_X, end_X], [-0.02, 1.02], color = 'k', lw = 1, zorder = 1000)

    # CONFIGURING THE GRAPH APPEARANCE
    #Set the x and y axis limits
    ax1.set_xlim(0, cumulative_chromosome_sizes['X']+(chromosome_sizes['X'])+100)

    # Changing the y axis limits
    ax1.set_ylim(0, 1.02)

    #x-axis ticks
    x_axis_tick_positions = []
    for k, v in cumulative_chromosome_sizes.items():
        left_start = v
        chrom = k
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
                      Line2D([0], [0], marker = 'o', color=allele1, alpha=1.0, markersize = 8, \
                      lw=0, label='phased allele 1'),
                      Line2D([0], [0], marker = 'o', color=allele2, alpha=1.0, markersize = 8, \
                      lw=0, label='phased allele 2'),
                      Line2D([0], [0], marker = 'o', color=hetero, alpha=1.0, markersize = 8, \
                      lw=0, label='heterozygous'),
                      Line2D([0], [0], marker = 'o', color=unphased, alpha=1.0, markersize = 8, \
                      lw=0, label='unphased')]

    ax1.legend(ncol=5, handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 1.07), fontsize = 13, frameon=True,\
              fancybox = True)

    #Title and axis labels
    ax1.set_title('phased BAFs (from SSCS): '+sample_name, y=1.125, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('BAF', fontsize = axislabelfont)
    ax1.set_xlabel('chromosome', fontsize = axislabelfont)

    plt.tight_layout()

    return  plt.savefig(out_directory+'/Phased_SNPs/Phased_SNPs_plots/'+sample_name+'_SSCS_SNPs_phased_BAFs_all_chromosomes.pdf')

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--phased_VCF", action="store", dest="phased_VCF", help="phased VCF file wildcard", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--panel_bed", action="store", dest="panel_bed", help="input bed file for the panel", required=True)
    parser.add_argument("--chromosome-ideogram", action="store", dest="chromosome_ideogram", help="input chromosome ideogram file", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    phased_VCF = o.phased_VCF
    sample_name = o.sample_name
    panel_bed = o.panel_bed
    ideogram_file = o.chromosome_ideogram
    out_directory = o.out_directory

    chromosome_sizes, chromosome_bands = chromosome_sizes_bands(ideogram_file)
    panel_coverage = panel_coverage_dict(panel_bed)
    variants_dict = variants_dictionary(phased_VCF)

    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        print('plotting BAFs across chromosome'+str(i))
        chromosome_SNPs_plot(i, variants_dict, panel_coverage, ideogram_file, chromosome_sizes, out_directory, sample_name)

    cumulative_chromosome_sizes = {}
    length = 0
    for i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']:
        chromosome_size = chromosome_sizes[i]
        cumulative_chromosome_sizes[i]=length
        length+=chromosome_size

    print('plotting BAFs across all chromosomes on one plot')
    all_chromosomes_SNPs_plot(variants_dict, cumulative_chromosome_sizes, panel_coverage, sample_name, chromosome_sizes, out_directory)

    return print('phased SNP plots complete')

if __name__ == "__main__":
	main()
