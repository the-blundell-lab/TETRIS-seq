#!/usr/bin/env python

'''''
Watson code for Grouping Reads by UMI and SSCS calling.
Version 2.1 (August 2021)

Input:
    1) mapped merged BAM (paired end), containing UMI information the RX tag
    2) sample name
    3) threshold for consensus calling
    4) max Ns permitted
    5) minimum mapping quality
    6) minimum per-base sequencing quality

Outputs:
    BAM files:
        1) paired-end BAM containing SSCs
        2) BAM file containing unpaired SSCSs (e.g. if the partner is filtered out due to too many Ns)
    Metrics files:
        1) metrics txt file (and histograms) containing info on distribution of read 1, read 2, forward and reverse reads from grouping by UMI
        2) metrics txt file (and histograms) containing info on distribution of read 1, read 2, forward and reverse reads for the SSCS
        3) metrics txt file (and plot) containing UMI read family size information and number of unique molecules
        4) UMI family size subplots broken down by read1, read2, forward and reverse.
        5) plot showing mapping quality distribution and plot showing number of reads filtered out due to low mapping quality

Usage:
Watson_code_SSCS_calling_2.1.py  --infile mapped_merged_BAM --sample-name sample_name --min-family-size NUMBER(default = 1, range 1-no max)"
                                 --threshold NUMBER(default=0.9, range 0-1) --max_N (default=1.0, range: 0-1) --outbam SSCS_bam_filename
                                 --min-mapping-quality NUMBER(default=20)
                                 --min-base-quality NUMBER(default=20) --unpaired-outbam unpaired_SSCS_bam_filename --out-directory directory to save files in

'''''
version = '2.1'

from argparse import ArgumentParser
import pysam
import sys
import gzip
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MaxNLocator, MultipleLocator
import numpy as np
from array import array
import timeit
import time
import shelve
from datetime import date
from Bio.Seq import Seq
from ast import literal_eval

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

#Group reads by UMI
def Group_reads_by_UMI(sample_name, merged_bam, out_directory, min_mapping_qual, min_family_size):

    infile = merged_bam #BAM file produced from MergeBamAlignment
    in_bam = pysam.Samfile(infile, "rb")
    bam_entry = in_bam.fetch(until_eof=True)

    read_number = 0
    good_flag_reads=0
    bad_flag_reads=0
    good_mapping_qual=0
    bad_mapping_qual=0
    read_dict = shelve.open('./TEMP/'+sample_name+'_grouped_reads_dict', flag = 'n')
    # read_dict = {}
    mapping_qualities = {} #dictionary storing the distribution of mapping qualities

    good_flag_list = [99, 147, 83, 163]
    # 99 = read paired (0x1), read mapped in proper pair (0x2), mate reverse strand (0x20), first in pair (0x40)
    # 147 = read paired (0x1), read mapped in proper pair (0x2), read reverse strand (0x10), second in pair (0x80)
    # 83 = read paired (0x1), read mapped in proper pair (0x2), read reverse strand (0x10), first in pair (0x40)
    # 163 = read paired (0x1), read mapped in proper pair (0x2), mate reverse strand (0x20), second in pair (0x80)

    read1_forward_inc_unmapped = 0
    read1_reverse_inc_unmapped = 0
    read2_forward_inc_unmapped = 0
    read2_reverse_inc_unmapped = 0

    read1_forward_raw = 0
    read1_reverse_raw = 0
    read2_forward_raw = 0
    read2_reverse_raw = 0

    read1_forward_raw_passed_flag = 0
    read1_reverse_raw_passed_flag = 0
    read2_forward_raw_passed_flag = 0
    read2_reverse_raw_passed_flag = 0

    read1_forward_good_mapping_qual = 0
    read1_reverse_good_mapping_qual = 0
    read2_forward_good_mapping_qual = 0
    read2_reverse_good_mapping_qual = 0

    start_time = time.time()
    for line in bam_entry:
        #### Collecting data on read strand direction and distribution (including unmapped) #####
        if line.is_read1 is True:
            if line.is_reverse is True:
                read1_reverse_inc_unmapped+=1
            if line.is_reverse is False:
                read1_forward_inc_unmapped+=1
        if line.is_read1 is False:
            if line.is_reverse is True:
                read2_reverse_inc_unmapped+=1
            if line.is_reverse is False:
                read2_forward_inc_unmapped+=1

        #### Grouping reads by UMI ######
        if line.is_unmapped is False: #check the read is mapped...
            sequence = line.seq #sequence
            flag = line.flag #e.g. 99 or 147
            cigar = line.cigartuples #cigars indicate the mapping of the read and indicate if there are insertions or deletions etc, e.g. 3M2I110M
            template_length = line.template_length #e.g. 179, -179 etc...(same as insert size (deprecated))
            chromosome = line.rname #chromosome (coded)
            coordinate = line.pos #the coordinate of the read (1 based coordinate)
            partner_chromosome = line.mrnm #next reference ID
            partner_coordinate = line.next_reference_start #the coordinate of the read pair (mpos) (1 based coordinate)
            qname = line.qname #query name
            UMI = line.get_tag('RX') #the UMI
            quality = list(line.query_qualities)
            mapq = line.mapping_quality

            if mapq in mapping_qualities.keys():
                mapping_qualities[mapq]+=1
            else:
                mapping_qualities[mapq]=1

            if line.is_read1 is True:
                read_pair = 'read1'
            if line.is_read2 is True:
                read_pair = 'read2'
            if line.is_reverse is True:
                strand = 'reverse'
            if line.is_reverse is False:
                strand = 'forward'

            read_info = (sequence, cigar, qname, quality, mapq)
            key = str((UMI, chromosome, coordinate, read_pair, strand, flag, partner_chromosome, partner_coordinate, template_length))

            ##### COLLECTING INFORMATION ON on read strand direction and distribution FOR MAPPED READS ########
            if read_pair == 'read1' and strand == 'reverse':
                read1_reverse_raw+=1
            if read_pair == 'read1' and strand == 'forward':
                read1_forward_raw+=1
            if read_pair == 'read2' and strand == 'reverse':
                read2_reverse_raw+=1
            if read_pair == 'read2' and strand == 'forward':
                read2_forward_raw+=1

            ##### CREATE DICTIONARY CONTAINING READS GROUPED BY UMIS ########
            if flag in good_flag_list:
                good_flag_reads+=1
                if mapq >= min_mapping_qual: #only add them to the dictionary if have mapping quality higher than threshold
                    good_mapping_qual+=1

                    if key in read_dict.keys(): #if the UMI-tagged molecule is already in the dictionary...
                        existing_read_info = read_dict[key]
                        existing_read_info.append(read_info)
                        read_dict[key] = existing_read_info #append the read to the family

                    if key not in read_dict.keys(): #look to see if the key is already in the dictionary..
                        read_dict[key] = [read_info] #if partner isn't there, then add this key

                else:
                    bad_mapping_qual+=1

            else:
                bad_flag_reads+=1

            ##### COLLECTING INFORMATION ON on read strand direction and distribution for MAPPED READS that passed filters ########
            if flag in good_flag_list:
                if read_pair == 'read1' and strand == 'reverse':
                    read1_reverse_raw_passed_flag+=1
                if read_pair == 'read1' and strand == 'forward':
                    read1_forward_raw_passed_flag+=1
                if read_pair == 'read2' and strand == 'reverse':
                    read2_reverse_raw_passed_flag+=1
                if read_pair == 'read2' and strand == 'forward':
                    read2_forward_raw_passed_flag+=1

                if mapq >= min_mapping_qual:
                    if read_pair == 'read1' and strand == 'reverse':
                        read1_reverse_good_mapping_qual+=1
                    if read_pair == 'read1' and strand == 'forward':
                        read1_forward_good_mapping_qual+=1
                    if read_pair == 'read2' and strand == 'reverse':
                        read2_reverse_good_mapping_qual+=1
                    if read_pair == 'read2' and strand == 'forward':
                        read2_forward_good_mapping_qual+=1

            read_number+=1

            if read_number % 1000000 == 0:
                print('Grouping reads by UMI for sample '+sample_name+': reads processed = ', read_number)
                print('time for last 1,000,000 reads to be processed = %s seconds' % int(time.time() - start_time))
                start_time = time.time() #reset the timer so it can calculate the time for the next 100,000 reads

    print('number of reads with good flags = ', good_flag_reads)
    print('number of reads with bad flags = ', bad_flag_reads)

    metrics_file = open(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_read_distributions_MUFS'+str(min_family_size)+'.txt', 'w')
    metrics_file.write('sample name :\t'+ str(sample_name)+ '\n')
    metrics_file.write('date of analysis :\t'+ str(date_today)+ '\n')
    metrics_file.write('produced from code:\t' + 'Watson_code_SSCS_calling: version '+str(version) + '\n')
    metrics_file.write('produced from function:\t' + 'Group_reads_by_UMI' + '\n')
    metrics_file.write('minimum mapping quality set to:\t' + str(min_mapping_qual) + '\n')
    metrics_file.write('\n')
    metrics_file.write('reads\t'+ 'read1_forward\t' + 'read1_reverse\t' + 'read2_forward\t' +  'read2_reverse\t' + '\n')
    metrics_file.write('including_unmapped\t' + str(read1_forward_inc_unmapped)+'\t' + str(read1_reverse_inc_unmapped)+'\t' + str(read2_forward_inc_unmapped)+'\t' +  str(read2_reverse_inc_unmapped)+'\t' + '\n')
    metrics_file.write('mapped\t' + str(read1_forward_raw)+'\t' + str(read1_reverse_raw)+'\t' + str(read2_forward_raw)+'\t' +  str(read2_reverse_raw)+'\t' + '\n')
    metrics_file.write('mapped_passed_flags\t' + str(read1_forward_raw_passed_flag)+'\t' + str(read1_reverse_raw_passed_flag)+'\t' + str(read2_forward_raw_passed_flag)+'\t' + str(read2_reverse_raw_passed_flag)+'\t' + '\n')
    metrics_file.write('mapped_passed_flags_good_mapping_quality (>= '+str(min_mapping_qual)+')\t' + str(read1_forward_good_mapping_qual)+'\t' + str(read2_forward_good_mapping_qual)+'\t' + str(read2_forward_good_mapping_qual)+'\t' + str(read2_reverse_good_mapping_qual)+'\t' + '\n')
    metrics_file.close()

    ######PLOT METRICS#####
    f, (axes) = plt.subplots(1, 3, figsize=(30, 8))

    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]

    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15
    exonfont=12

    #INCLUDING UNMAPPED
    x = ['read1_forward', 'read2_reverse', 'read2_forward', 'read1_reverse']
    y = [read1_forward_inc_unmapped, read2_reverse_inc_unmapped, read2_forward_inc_unmapped, read2_reverse_inc_unmapped]

    ax1.bar(x, y, color = [c1, c2, c3, c4], lw = 2, zorder = 5)
    ax1.set_title('Read strand/ direction distribution of raw reads (inc. unmapped)\n('+sample_name+')', y=1.05, fontsize = titlefont, fontweight='bold')

    #MAPPED
    x = ['read1_forward', 'read2_reverse', 'read2_forward', 'read1_reverse']
    y = [read1_forward_raw, read2_reverse_raw, read2_forward_raw, read2_reverse_raw]

    ax2.bar(x, y, color = [c1, c2, c3, c4], lw = 2, zorder = 5)
    ax2.set_title('Read strand/ direction distribution of raw reads (mapped)\n('+sample_name+')', y = 1.05, fontsize = titlefont, fontweight='bold')

    #MAPPED + GOOD FLAGS
    x = ['read1_forward', 'read2_reverse', 'read2_forward', 'read1_reverse']
    y = [read1_forward_raw_passed_flag, read2_reverse_raw_passed_flag, read2_forward_raw_passed_flag, read2_reverse_raw_passed_flag]

    ax3.bar(x, y, color = [c1, c2, c3, c4], lw = 2, zorder = 5)
    ax3.set_title('Read strand/ direction distribution of raw reads (mapped) \nthat passed flag filter\n('+sample_name+')', y = 1.05, fontsize = titlefont, fontweight='bold')

    for ax in axes.flatten():
        ax.set_xlabel('Read strand/ direction', fontsize = axislabelfont)
        ax.set_ylabel('Number of reads', fontsize = axislabelfont)
        ax.tick_params(axis='both',labelsize=axislabelfont)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color('#969696')
        ax.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
        ax.xaxis.set_tick_params(width=1, color = '#969696', length = 6)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_grouped_reads_UMI_read_strand_distributions_MUFS'+str(min_family_size)+'.pdf')

    #PLOT METRICS FOR MAPPING QUALITY##
    plt.close('all')
    f, (axes) = plt.subplots(1, 2, figsize=(20, 8))

    ax1 = axes[0]
    ax2 = axes[1]

    axisfont=15
    titlefont=15
    axislabelfont=15

    #Good mapping quality
    x = ['read1_forward', 'read2_reverse', 'read2_forward', 'read1_reverse']
    y = [read1_forward_good_mapping_qual, read1_reverse_good_mapping_qual, read2_forward_good_mapping_qual, read2_reverse_good_mapping_qual]

    ax1.bar(x, y, color = [c1, c2, c3, c4], lw = 2, zorder = 5)
    ax1.set_title('Read strand/ direction distribution of raw reads with mapping quality >= '+str(min_mapping_qual)+'\n('+sample_name+')', y=1.05, fontsize = titlefont, fontweight='bold')

    #Bad mapping quality
    x = ['read1_forward', 'read2_reverse', 'read2_forward', 'read1_reverse']
    y = [read1_forward_raw_passed_flag-read1_forward_good_mapping_qual,
         read1_reverse_raw_passed_flag-read1_reverse_good_mapping_qual,
         read2_forward_raw_passed_flag-read2_forward_good_mapping_qual,
         read2_reverse_raw_passed_flag-read2_reverse_good_mapping_qual]

    ax2.bar(x, y, color = [c1, c2, c3, c4], lw = 2, zorder = 5)
    ax2.set_title('Read strand/ direction distribution of raw reads filtered out due to mapping quality < '+str(min_mapping_qual)+'\n('+sample_name+')', y=1.05, fontsize = titlefont, fontweight='bold')

    for ax in axes.flatten():
        ax.set_xlabel('Read strand/ direction', fontsize = axislabelfont)
        ax.set_ylabel('Number of reads', fontsize = axislabelfont)
        ax.tick_params(axis='both',labelsize=axislabelfont)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color('#969696')
        ax.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
        ax.xaxis.set_tick_params(width=1, color = '#969696', length = 6)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_grouped_reads_strands_mapping_qualities_MUFS'+str(min_family_size)+'.pdf')

    #PLOT MAPPING QUALITY DISTRIBUTION##
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, figsize=(12, 8))

    axisfont=15
    titlefont=15
    axislabelfont=15

    x = []
    y = []
    for k, v in mapping_qualities.items():
        x.append(k)
        y.append(v)

    ax1.bar(x, y, color = blue2, zorder = 5)
    ax1.set_title('Mapping qualities pre-SSCS calling\n('+sample_name+')', y=1.05, fontsize = titlefont, fontweight='bold')

    ax1.set_xlabel('mapping quality', fontsize = axislabelfont)
    ax1.set_ylabel('number of reads', fontsize = axislabelfont)
    ax1.tick_params(axis='both',labelsize=axislabelfont)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(width=1, color = '#969696', length = 6)

    ax1.set_yscale('log')

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_grouped_reads_mapping_qualities_distribution_MUFS'+str(min_family_size)+'.pdf')

    return read_dict

def consensus_read_caller_with_qualities(grouped_reads, threshold, read_length, quality_threshold): #function for calling the consensus sequence from the UMI family
    #read length e.g. 146 for the IDT duplex adaptors kit
    #threshold is the proportion of reads at a position that should match to call a consensus, e.g. 0.9 (90%)
    #grouped_reads = the reads that share a UMI/ positions etc and have the same CIGAR (the most common CIGAR for that family)
    #grouped_reads = [(sequence, base qualities), (sequence, base qualities)...]

    consensus_read = '' #start the read as empty (and then will add to it...)
    quality_read = [] #the summed Phred scores for each position that make up the consensus (capped at 40)
    consensus_depth = []
    consensus_errors = []
    read_depths = []
    proportion_not_consensus = []
    prob_error = [] #the summed Phred scores for each position that make up the consensus (uncapped)

    for base_number in range(read_length): #goes through each base (with length e.g 146 bases)
        nucleotides = {'T': 0, 'C': 0, 'G': 0, 'A': 0, 'N': 0}
        base_qualities = {'T': 0, 'C': 0, 'G': 0, 'A': 0, 'N': 0}
        total_number_nucleotides = 0 #reset
        read_depth=0 #keep track of the number of reads that pass the min quality filter at this position

        for read in range(len(grouped_reads)): #go through each read in the list (where j is each read)
            try:
                base_quality = grouped_reads[read][1][base_number]
                if base_quality >= quality_threshold:
                    read_depth +=1
                    if grouped_reads[read][0][base_number] == 'T':
                        nucleotides['T'] +=1
                        base_qualities['T']+=base_quality #multiply the base qualtities for all reads at that position carrying that base
                    if grouped_reads[read][0][base_number] == 'C':
                        nucleotides['C'] +=1
                        base_qualities['C']+=base_quality
                    if grouped_reads[read][0][base_number] == 'G':
                        nucleotides['G'] +=1
                        base_qualities['G']+=base_quality
                    if grouped_reads[read][0][base_number] == 'A':
                        nucleotides['A'] +=1
                        base_qualities['A']+=base_quality
                    if grouped_reads[read][0][base_number] == 'N':
                        nucleotides['N'] +=1
                        base_qualities['N']+=base_quality
                    total_number_nucleotides +=1
            except:
                break

        try: #for each base, work out which nucleotides (if any) is >e.g. 90% of the nucleotides at that position
            for base_option in ['T', 'C', 'G', 'A', 'N']:
                if float(nucleotides[base_option]/total_number_nucleotides) > threshold: #if number of particular nucleotide divided by total is > e.g. 0.9 = nucleotide to call
                    consensus_read += base_option #add that nucleotide to the consensus read
                    quality_score = base_qualities[base_option]
                    prob_error.append(quality_score)
                    if quality_score > 40:
                        quality_read.append(40) #quality scores are unsigned char (1-255), but max seems to be 40 (was outputting No quality scores in BAM when too high (even if <==255))
                    if quality_score <= 40:
                        quality_read.append(quality_score)

                    consensus_depth.append(nucleotides[base_option])
                    consensus_errors.append(total_number_nucleotides-nucleotides[base_option])

                    proportion_supporting_consensus = float(nucleotides[base_option]/total_number_nucleotides)
                    proportion_not_supporting_consensus = 1-proportion_supporting_consensus

                    proportion_not_consensus.append(proportion_not_supporting_consensus)

                    break

                elif base_option == 'N': #if get to N and none of the nucleotides have been >e.g. 90% of the reads, then call an 'N'
                    consensus_read += 'N'
                    quality_read.append(0)
                    prob_error.append(0)
                    consensus_depth.append(0)
                    consensus_errors.append(total_number_nucleotides)

        except:
            consensus_read +='N'
            quality_read.append(0)
            prob_error.append(0)
            consensus_depth.append(0)
            consensus_errors.append(total_number_nucleotides)

        read_depths.append(read_depth) #append the read depth (for reads that pass the quality filter) at this position

    max_read_depth = max(read_depths) #max read depth for the consensus
    min_read_depth = min(read_depths) #min read depth for the consensus (i.e. when poor quality bases were filtered out)

    consensus_error_rate = np.mean(proportion_not_consensus) #fraction of bases in raw reads disagreeing with the final consensus calls

    quality_reads = array('B', quality_read)
    quality_reads_string = ''
    for i in quality_read:
        quality_reads_string+=chr(i+33)

    return consensus_read, quality_reads, quality_reads_string, len(grouped_reads), consensus_error_rate, consensus_depth, consensus_errors, max_read_depth, min_read_depth, prob_error

def SSCS_calling_BAM(sample_name, grouped_read_dict, merged_bam, threshold, min_family_size, max_N, SSCS_bam, unpaired_SSCS_bam, out_directory, quality_threshold, version, min_mapping_qual):
    #sample name, e.g. C92_024_s8
    #grouped_read_dict is the output of the Group_reads_by_UMI function
    #merged_bam acts as a template bam file
    #threshold is the proportion of reads at a position that should match to call a consensus, e.g. 0.9
    #min_family_size is the minimum number of reads required in a UMI famaily (with the same CIGAR) to call a consensus sequence

    #CREATE HEADER FOR NEW BAM FILE, THAT DOESN'T CONTAIN THE MERGED BAM'S PROGRAM GROUP ID (if included, causes problems with subsequent mergeBam)
    #Create a new PG that contains the command run for this code
    new_header = {}
    in_bam = pysam.Samfile(merged_bam, "rb")
    heading = in_bam.header
    new_header['HD'] = heading['HD']
    new_header['SQ'] = heading['SQ']
    new_header['RG'] = heading['RG']
    new_PG = [{'ID': 'Watson_code_SSCS_calling_v'+str(version),
    'PN': 'Watson_code_CNV_panel_v'+str(version),
    'VN': version,
    'CL': 'Watson_code_SSCS_calling_'+str(version)+'.py --infile '+merged_bam+' --sample-name '+sample_name+' --min-family-size '+str(min_family_size)+\
    ' --threshold '+str(threshold)+ ' --max_N '+str(max_N)+' --outbam '+ SSCS_bam+ \
    ' --min-mapping-quality '+str(min_mapping_qual)+' --min-base-quality '+str(quality_threshold)+' --unpaired-outbam '+unpaired_SSCS_bam+' --out-directory '+out_directory}]
    new_header['PG'] = new_PG

    #CREATE NEW BAM FILES
    out_bam = pysam.Samfile(SSCS_bam, "wb", header=new_header) #open a new BAM file to write the SSCS to
    unpaired_out_bam = pysam.Samfile(unpaired_SSCS_bam, "wb", header=new_header) #open a new BAM file to write the SSCS to
    # fastq = gzip.open(SSCS_fastq, "wt") #interleaved fastq file (r1 followed by r2 etc..)

    #N.B. the structure of the qname in the new files will be UMI : read1 UMI family size: read2 UMI family size: arbitrary counter number (specific to a pair of reads)

    def most_common(lst): #function for retrieving the most common element in a non-hashable type of list (i.e. to retrieve most common CIGAR from a list of CIGARS)
        return max((lst), key=lst.count)

    def reverse(quality_scores):
        return quality_scores[::-1]

    UMI_family_size_dict = {}
    # consensus_groups = shelve.open('./TEMP/'+sample_name+'_consensus_groups_dict.db')
    consensus_groups = {}
    #a dictionary where each key is a different UMI family (different key for each read pair) and the value is the consensus sequence, cigar, UMI family size, quality and new tags

    UMI_family = 0 #a running total of the number of UMI families used for consensus calling
    reads_included_in_consensus_calling = 0 #a running count of the number of reads that were included in consensus calling
    reads_filtered_out_due_to_cigar = 0 #a running count of the number of reads filtered out because they didn't have consensus CIGAR
    too_many_Ns_in_consensus = 0 #a running count of the number of reads filtered out because they had too many Ns in

    start_time = time.time() #start the timer (so can record how long the processes are taking)

    pair_counter = 0

    proportion_N_dict = {}

    for k, v in grouped_read_dict.items(): #reads grouped by UMI, with the following structure...
        #k = (UMI, chromosome, coordinates, read_pair, strand, flag, partner chrom, partner coord, template_length)
        #v = [(sequence, cigar, query name, quality, mapping quality), (sequence, cigar, query name), etc....]
        k = literal_eval(k) #convert the key (which is a string) pack to a tuple

        read_length = len(v[0][0]) #look at the first sequence and check how long it is (they are all the same)

        consensus_group = [] #a list to store the sequences that have the same UMI (if they have the most common CIGAR of the group)
        CIGAR_group = [] #a list to store all the different CIGARS in the UMI family

        #find the most common CIGAR in the UMI family group
        for i in v: #iterate over each read in the dictionary that has the same key and find the most common CIGAR string
            CIGAR = i[1]
            CIGAR_group.append(CIGAR)

        most_common_CIGAR = most_common(CIGAR_group) #if there are only 2 and they are different then it takes the first one

        #filter out the reads in the UMI family group that don't have the common CIGAR (create a list of the sequences that share the common CIGAR)
        for i in v: #iterate over each read in the dictionary that has the same key and find the sequences that have the most common CIGAR and put them in list
            sequence = i[0]
            CIGAR = i[1]
            quality = i[3] #list of qualities for each base
            mapq = i[4] #mapping quality of the read
            if CIGAR == most_common_CIGAR:
                consensus_group.append((sequence, quality))
                reads_included_in_consensus_calling +=1 #this read will be included in the consensus calling
            else:
                reads_filtered_out_due_to_cigar +=1 #the reads that don't have the common CIGAR will not be included in consensus calling

        #Store the UMI family size in the UMI family size dictionary (so can plot family size distribution)
        UMI_family_size = len(consensus_group)
        if UMI_family_size in UMI_family_size_dict.keys(): #Add the UMIs to a dictionary so can calculate UMI family size distribution
            UMI_family_size_dict[UMI_family_size]+=1
        else:
            UMI_family_size_dict[UMI_family_size]=1


        #make a consensus sequences from all the sequences that share the same CIGAR string (the most common CIGAR string)
        if UMI_family_size>=min_family_size: #check to make sure there are enough reads in the family (i.e. more than whatever minimum size was set)...
            UMI_family+=1 #add to the running counter of UMI families called
            consensus_read, quality_read, quality_reads_string, family_size, consensus_error_rate, consensus_depth, consensus_errors, max_read_depth, min_read_depth, prob_error = consensus_read_caller_with_qualities(consensus_group, threshold, read_length, quality_threshold) #call the consensus from the reads

            N_count = consensus_read.count('N') #count the number of Ns in the consensus read
            consensus_read_length = len(consensus_read)
            proportion_N = N_count/consensus_read_length

            if proportion_N in proportion_N_dict.keys(): #keep a tally of what proportion of the reads are N
                proportion_N_dict[proportion_N]+=1
            else:
                proportion_N_dict[proportion_N]=1

            if proportion_N < max_N: #if proportion of Ns in the consensus is less than the set limit:
                #k = str(UMI, chromosome, coordinates, read_pair, strand, flag, partner chrom, partner coord, template_length)
                UMI = k[0]
                chromosome = k[1]
                coordinate = k[2]
                read_pair = k[3]
                strand = k[4]
                flag = k[5]
                partner_chromosome = k[6]
                partner_coordinate = k[7]
                template_length = k[8]
                consensus = consensus_read
                quality = quality_read
                quality_string = quality_reads_string
                UMI_family_size = family_size
                CIGAR = most_common_CIGAR
                cD_tag = max_read_depth
                cM_tag = min_read_depth
                cE_tag = consensus_error_rate
                cd_tag = consensus_depth
                ce_tag = consensus_errors
                QS_tag = prob_error

                key = (UMI, chromosome, coordinate, read_pair, strand, flag, partner_chromosome, partner_coordinate, template_length)

                #create a dictionary where each key is a different UMI family (different key for each read pair) and the value is the consensus sequence and cigar
                consensus_groups[key]=(consensus, CIGAR, UMI_family_size, quality, cD_tag, cM_tag, cE_tag, cd_tag, ce_tag, QS_tag, quality_string)

                #look to see if the partner pair is already in the dictionary...
                if read_pair == 'read1':
                    partner_read_pair = 'read2'
                if read_pair == 'read2':
                    partner_read_pair = 'read1'
                if strand == 'forward':
                    partner_strand = 'reverse'
                if strand == 'reverse':
                    partner_strand = 'forward'

                if flag == 99:
                    partner_flag = 147
                if flag == 147:
                    partner_flag = 99
                if flag == 83:
                    partner_flag = 163
                if flag == 163:
                    partner_flag = 83

                if template_length <0: #i.e. negative number
                    partner_template_length = abs(template_length) #i.e. the positive version
                if template_length >0: #i.e. positive number
                    partner_template_length = -template_length
                if template_length == 0:
                    partner_template_length = 0

                partner_key = (UMI, partner_chromosome, partner_coordinate, partner_read_pair, partner_strand, partner_flag, chromosome,
                              coordinate, partner_template_length)

                #create the read entry for the BAM file...
                new_read = pysam.AlignedRead() #Create bam file reads for the new SSCS
                new_read.seq = consensus
                new_read.flag = flag #contains info on whether read 1, read2, forward or reverse
                new_read.cigartuples = CIGAR
                new_read.rname = chromosome
                new_read.pos = coordinate
                new_read.mapq = 225 #arbitrary mapping quality
                new_read.template_length = template_length
                new_read.query_qualities = quality
                new_read.mrnm = partner_chromosome
                new_read.next_reference_start = partner_coordinate
                new_read.set_tags([('cD', cD_tag, "i"), ('cM', cM_tag, "i"), ('cE', cE_tag, "f"), ('cd', cd_tag), ('ce', ce_tag), ('RG', sample_name), ('QS', QS_tag)])

                #look to see if the partner key is already in the dictionary...(it won't be for the first read looked at, but then it will find the first read when it's partner is later added)
                if partner_key in consensus_groups.keys(): #i.e. if the reads are paired (has the partner survived the N filter and min family size filter?)
                    #create the partner read entry for the BAM file...
                    partner = consensus_groups[partner_key] #retrieve the 'values' for the partner (i.e. (consensus, CIGAR, UMI family size))
                    partner_consensus = partner[0]
                    partner_CIGAR = partner[1]
                    partner_UMI_family_size = partner[2]
                    partner_quality = partner[3]
                    partner_cD_tag = partner[4]
                    partner_cM_tag = partner[5]
                    partner_cE_tag = partner[6]
                    partner_cd_tag = partner[7]
                    partner_ce_tag = partner[8]
                    partner_QS_tag = partner[9]
                    partner_quality_string = partner[10]

                    partner_read = pysam.AlignedRead() #Create bam file reads for the new SSCS
                    partner_read.seq = partner_consensus
                    partner_read.flag = partner_flag #contains info on whether read 1, read2, forward or reverse
                    partner_read.cigartuples = partner_CIGAR
                    partner_read.rname = partner_chromosome
                    partner_read.pos = partner_coordinate
                    partner_read.mapq = 225 #arbitrary mapping quality
                    partner_read.template_length = partner_template_length
                    partner_read.query_qualities = partner_quality
                    partner_read.mrnm = chromosome
                    partner_read.next_reference_start = coordinate
                    partner_read.set_tags([('cD', partner_cD_tag, "i"), ('cM', partner_cM_tag, "i"), ('cE', partner_cE_tag, "f"), ('cd', partner_cd_tag), ('ce', partner_ce_tag), ('RG', sample_name), ('QS', partner_QS_tag)])

                    if read_pair == 'read1':
                        new_read.qname = UMI+':'+str(UMI_family_size)+':'+str(partner_UMI_family_size)+':'+str(UMI_family) #read 1 UMI family size followed by read 2 UMI family size
                        partner_read.qname = UMI+':'+str(UMI_family_size)+':'+str(partner_UMI_family_size)+':'+str(UMI_family) #same qname as its partner

                        #write the read to the BAM file
                        out_bam.write(new_read) #write the new SSCS to new bam file (read 1 first)
                        #then write it's partner to the BAM file...
                        out_bam.write(partner_read) #write the new SSCS to new bam file

                        # #add a /1 or /2 on the end for the fastq qnames
                        # new_read.qname = UMI+':'+str(UMI_family_size)+':'+str(partner_UMI_family_size)+':'+str(UMI_family)+'/1' #read 1 UMI family size followed by read 2 UMI family size
                        # partner_read.qname = UMI+':'+str(UMI_family_size)+':'+str(partner_UMI_family_size)+':'+str(UMI_family)+'/2' #same qname as its partner
                        #
                        # if strand == 'forward':
                        #     #first write read 1...
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (new_read.qname, new_read.seq, quality_string))
                        #     #then write read 2...
                        #     #partner must be on reverse strand so needs to be reverse complemented...
                        #     sequence = Seq(partner_consensus)
                        #     partner_read.seq = str(sequence.reverse_complement())
                        #     partner_read.query_qualities = reverse(partner_quality)
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (partner_read.qname, partner_read.seq, partner_quality_string))
                        #
                        # if strand == 'reverse':
                        #     #first write read 1...
                        #     #take in to account it is on reverse strand so needs to be reverse complemented...
                        #     sequence = Seq(consensus)
                        #     new_read.seq = str(sequence.reverse_complement())
                        #     new_read.query_qualities = reverse(quality)
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (new_read.qname, new_read.seq, quality_string))
                        #     #then write read 2...
                        #     #partner is on the forward strand so does not need to be reverse complemented
                        #     partner_read.seq = partner_consensus
                        #     partner_read.query_qualities = partner_quality
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (partner_read.qname, partner_read.seq, partner_quality_string))

                    if read_pair == 'read2':
                        new_read.qname = UMI+':'+str(partner_UMI_family_size)+':'+str(UMI_family_size)+':'+str(UMI_family)
                        partner_read.qname = UMI+':'+str(partner_UMI_family_size)+':'+str(UMI_family_size)+':'+str(UMI_family)

                        #write the partner to the BAM file first...
                        out_bam.write(partner_read) #write the new SSCS to new bam file (read 1 first)
                        #then write the read to the BAM file...
                        out_bam.write(new_read) #write the new SSCS to new bam file

                        # #add a /1 or /2 on the end for the fastq qnames
                        # new_read.qname = UMI+':'+str(partner_UMI_family_size)+':'+str(UMI_family_size)+':'+str(UMI_family)+'/2'
                        # partner_read.qname = UMI+':'+str(partner_UMI_family_size)+':'+str(UMI_family_size)+':'+str(UMI_family)+'/1'
                        #
                        # if strand == 'forward':
                        #     #first write read 1..
                        #     #partner must be on reverse strand so needs to be reverse complemented...
                        #     sequence = Seq(partner_consensus)
                        #     partner_read.seq = str(sequence.reverse_complement())
                        #     partner_read.query_qualities = reverse(partner_quality)
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (partner_read.qname, partner_read.seq, partner_quality_string))
                        #     #then write read 2..
                        #     new_read.seq = consensus
                        #     new_read.query_qualities = quality
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (new_read.qname, new_read.seq, quality_string))
                        #
                        # if strand == 'reverse':
                        #     #first write read 1..
                        #     #partner is on the forward strand so does not need to be reverse complemented
                        #     partner_read.seq = partner_consensus
                        #     partner_read.query_qualities = partner_quality
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (partner_read.qname, partner_read.seq, partner_quality_string))
                        #     #then write read 2..
                        #     #take in to account it is on reverse strand so needs to be reverse complemented...
                        #     sequence = Seq(consensus)
                        #     new_read.seq = str(sequence.reverse_complement())
                        #     new_read.query_qualities = reverse(quality)
                        #     fastq.write('@:%s\n%s\n+\n%s\n' % (new_read.qname, new_read.seq, quality_string))

                    pair_counter+=1

                    del consensus_groups[partner_key] #delete the partner key from the dictionary now that it has been written to the BAM file/ fastq file
                    del consensus_groups[key]


            else: #if there are too many Ns in the consensus
                too_many_Ns_in_consensus+=1 #add to the running counter

            if UMI_family % 10000 == 0:
                print('Calling SSCSs: UMI families processed for sample '+sample_name+' = ', UMI_family)
                print('time for last 10,000 consensuses to be called = %s seconds' % int(time.time() - start_time))
                start_time = time.time() #reset the timer so it can calculate the time for the next 100 consensuses

    unpaired_counter = 0
    if len(consensus_groups)>0: #if there are any remaining UMI families in the dictinary (these will just be the remaining unpaired reads (i.e. those not already added to BAM/ fastq...)
        for k, v in consensus_groups.items():
            UMI = k[0]
            chromosome = k[1]
            coordinate = k[2]
            read_pair = k[3]
            strand = k[4]
            flag = k[5]
            partner_chromosome = k[6]
            partner_coordinate = k[7]
            template_length = k[8]
            consensus = v[0]
            CIGAR = v[1]
            UMI_family_size = v[2]
            quality_read = v[3]
            cD_tag = v[4]
            cM_tag = v[5]
            cE_tag = v[6]
            cd_tag = v[7]
            ce_tag = v[8]
            QS_tag = v[9]

            #create the read entry for the BAM file...
            new_read = pysam.AlignedRead() #Create bam file reads for the new SSCS
            new_read.seq = consensus
            new_read.cigartuples = CIGAR
            new_read.rname = chromosome
            new_read.pos = coordinate
            new_read.mapq = 225 #arbitrary mapping quality
            new_read.template_length = template_length
            new_read.query_qualities = quality_read
            new_read.mrnm = partner_chromosome
            new_read.next_reference_start = partner_coordinate
            new_read.set_tags([('cD', cD_tag, "i"), ('cM', cM_tag, "i"), ('cE', cE_tag, "f"), ('cd', cd_tag), ('ce', ce_tag), ('RG', sample_name), ('QS', QS_tag)])

            if read_pair == 'read1':
                new_read.qname = UMI+':'+str(UMI_family_size)+':'+str(0)+':'+str(unpaired_counter) #store that the partner (read2) had 0 UMI family size in the qname
            if read_pair == 'read2':
                new_read.qname = UMI+':'+str(0)+':'+str(UMI_family_size)+':'+str(unpaired_counter) #store that the partner (read1) had 0 UMI family size in the qname
            if strand == 'forward':
                new_read.flag = 0 #on forward strand, unpaired
            if strand == 'reverse':
                new_read.flag = 16 #on reverse strand, unpaired

            unpaired_out_bam.write(new_read) #write the unpaired reads to the unpaired BAM file
            unpaired_counter+=1 #add to the counter the number of unpaired reads


    out_bam.close() #close the new SSCS BAM file
    unpaired_out_bam.close() #close the new SSCS BAM file
    # fastq.close() #close the fastq file
    in_bam.close() #close the open in BAM file

    print('total number of reads included in consensus calling = ', reads_included_in_consensus_calling)
    print('total number of reads excluded due to non-common CIGAR = ', reads_filtered_out_due_to_cigar)
    print('total number of SSCS consensus reads (unique molecules) = ', UMI_family)
    print('total number of SSCS consensus reads filtered out due to too many Ns = ', too_many_Ns_in_consensus)
    print('total number of unpaired SSCS reads = ', unpaired_counter)
    print('total number of SSCS pairs = ', pair_counter)

    ##### CREATE METRICS FILE ########
    total_count = 0
    for k, v in UMI_family_size_dict.items():
        total_count+=v #count the number of UMI families
    print('total UMI family count = ', total_count)

    SSCS_UMI_metrics = open(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_SSCS_UMI_metrics_MUFS'+str(min_family_size)+'.txt', 'w')
    SSCS_UMI_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    SSCS_UMI_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    SSCS_UMI_metrics.write('produced from code:\t' + 'Watson_code_SSCS_calling: version '+str(version) + '\n')
    SSCS_UMI_metrics.write('produced from function:\t' + 'SSCS_calling_BAM' + '\n')
    SSCS_UMI_metrics.write('minimum family size required for SSCS consensus calling set to :\t'+ str(min_family_size)+ '\n')
    SSCS_UMI_metrics.write('minimum proportion of identical nucleotides required at each base set to :\t'+ str(threshold)+ '\n')
    SSCS_UMI_metrics.write('maximum proportion of Ns allowed :\t'+ str(max_N)+ '\n')
    SSCS_UMI_metrics.write('minimum base quality set to:\t' + str(quality_threshold) + '\n')
    SSCS_UMI_metrics.write('\n')
    SSCS_UMI_metrics.write('total number of reads included in SSCS consensus calling:\t'+ str(reads_included_in_consensus_calling)+ '\n')
    SSCS_UMI_metrics.write('total number of reads excluded due to non-common CIGAR:\t' + str(reads_filtered_out_due_to_cigar) + '\n')
    SSCS_UMI_metrics.write('total number of SSCS consensus reads (read 1 and read 2 separate):\t' + str(UMI_family) + '\n')
    SSCS_UMI_metrics.write('total number of unpaired SSCS reads:\t' + str(unpaired_counter) + '\n')
    SSCS_UMI_metrics.write('total number of SSCS pairs:\t' + str(pair_counter) + '\n')
    SSCS_UMI_metrics.write('\n')
    SSCS_UMI_metrics.write('UMI_family_size\t'+ 'count (read 1 and read 2 separate)\t' + 'fraction\t' + '\n')

    for k in sorted(UMI_family_size_dict.keys()):
        family_size = k
        number_molecules = UMI_family_size_dict[k]
        SSCS_UMI_metrics.write(str(family_size)+'\t' + str(number_molecules)+'\t' + str(number_molecules/total_count) + '\n')

    SSCS_UMI_metrics.close()

    SSCS_UMI_N_metrics = open(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_SSCS_UMI_N_metrics_MUFS'+str(min_family_size)+'.txt', 'w')
    SSCS_UMI_N_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    SSCS_UMI_N_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    SSCS_UMI_N_metrics.write('produced from code:\t' + 'Watson_code_SSCS_calling: version '+str(version) + '\n')
    SSCS_UMI_N_metrics.write('produced from function:\t' + 'SSCS_calling_BAM' + '\n')
    SSCS_UMI_N_metrics.write('minimum family size required for SSCS consensus calling set to :\t'+ str(min_family_size)+ '\n')
    SSCS_UMI_N_metrics.write('minimum proportion of identical nucleotides required at each base set to :\t'+ str(threshold)+ '\n')
    SSCS_UMI_N_metrics.write('maximum proportion of Ns allowed :\t'+ str(max_N)+ '\n')
    SSCS_UMI_N_metrics.write('\n')
    SSCS_UMI_N_metrics.write('total number of reads included in SSCS consensus calling:\t'+ str(reads_included_in_consensus_calling)+ '\n')
    SSCS_UMI_N_metrics.write('total number of reads excluded due to non-common CIGAR:\t' + str(reads_filtered_out_due_to_cigar) + '\n')
    SSCS_UMI_N_metrics.write('total number of SSCS consensus reads (read 1 and read 2 separate):\t' + str(UMI_family) + '\n')
    SSCS_UMI_N_metrics.write('total number of unpaired SSCS reads:\t' + str(unpaired_counter) + '\n')
    SSCS_UMI_N_metrics.write('total number of SSCS pairs:\t' + str(pair_counter) + '\n')
    SSCS_UMI_N_metrics.write('\n')
    SSCS_UMI_N_metrics.write('proportion of SSCS that is N\t'+ 'count (read 1 and read 2 separate)\t' + '\n')

    for k in sorted(proportion_N_dict.keys()):
        proportion_N = k
        number_molecules = proportion_N_dict[k]
        SSCS_UMI_N_metrics.write(str(proportion_N)+'\t' + str(number_molecules) + '\n')

    SSCS_UMI_N_metrics.close()

    ##### PLOT UMI FAMILY SIZE HISTOGRAM FOR ALL READS #####
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, figsize=(10, 8))
    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15

    x = []
    y = []
    total_reads = 0
    for k in sorted(UMI_family_size_dict.keys()):
        v = UMI_family_size_dict[k]
        x.append(int(k))
        y.append(int(v))
        total_reads+=(int(k)*int(v))

    mean_UMI_family_size = total_reads/total_count

    total_square_distance = 0
    for k in sorted(UMI_family_size_dict.keys()):
        v = UMI_family_size_dict[k]
        square_distance = v*((k-mean_UMI_family_size)**2)
        total_square_distance+=square_distance
    standard_deviation = np.sqrt(total_square_distance/total_count)

    ax1.plot(x, y, color = 'dodgerblue', lw = 2, zorder = 5)

    # CONFIGURING THE GRAPH APPEARANCE
    ax1.tick_params(axis='both',labelsize=axislabelfont)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_minor_locator(MultipleLocator(1))
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(which = 'minor', width=1, color = '#969696', length = 4)
    ax1.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    #Title and axis labels
    ax1.set_title('UMI family size distribution \n'+sample_name, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('number of read families', fontsize = axislabelfont)
    ax1.set_xlabel('UMI family size', fontsize = axislabelfont)

    ax1.text(0.95, 0.9, 'total reads = '+str(total_reads), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.85, 'total reads with family size of 1 = '+str(UMI_family_size_dict[1]), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.75, 'total unique molecules = '+str(UMI_family), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.7, 'total unique molecules (excl. family size 1) = '+str(UMI_family-UMI_family_size_dict[1]), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')

    ax1.text(0.95, 0.6, 'mean UMI family size = '+str(round(mean_UMI_family_size, 2))+' (std = '+str(round(standard_deviation, 2))+')', transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.55, 'max UMI family size = '+str(max(x)), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_SSCS_UMI_family_size_distribution_MUFS'+str(min_family_size)+'.pdf')


    ##### PLOT UMI FAMILY SIZE HISTOGRAM FOR ALL READS (JUST FAMILY SIZE 1-20) #####
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, figsize=(10, 8))
    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15

    x = []
    y = []
    total_reads = 0
    for k in sorted(UMI_family_size_dict.keys()):
        v = UMI_family_size_dict[k]
        x.append(int(k))
        y.append(int(v))
        total_reads+=(int(k)*int(v))

    mean_UMI_family_size = total_reads/total_count

    total_square_distance = 0
    for k in sorted(UMI_family_size_dict.keys()):
        v = UMI_family_size_dict[k]
        square_distance = v*((k-mean_UMI_family_size)**2)
        total_square_distance+=square_distance
    standard_deviation = np.sqrt(total_square_distance/total_count)

    ax1.plot(x, y, color = 'dodgerblue', lw = 2, zorder = 5)

    # CONFIGURING THE GRAPH APPEARANCE
    ax1.tick_params(axis='both',labelsize=axislabelfont)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_minor_locator(MultipleLocator(1))
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(which = 'minor', width=1, color = '#969696', length = 4)
    ax1.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    #Title and axis labels
    ax1.set_title('UMI family size distribution \n'+sample_name, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('number of read families', fontsize = axislabelfont)
    ax1.set_xlabel('UMI family size', fontsize = axislabelfont)
    ax1.set_xlim(0, 20)

    ax1.text(0.95, 0.9, 'total reads = '+str(total_reads), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.85, 'total reads with family size of 1 = '+str(UMI_family_size_dict[1]), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.75, 'total unique molecules = '+str(UMI_family), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.7, 'total unique molecules (excl. family size 1) = '+str(UMI_family-UMI_family_size_dict[1]), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')

    ax1.text(0.95, 0.6, 'mean UMI family size = '+str(round(mean_UMI_family_size, 2))+' (std = '+str(round(standard_deviation, 2))+')', transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')
    ax1.text(0.95, 0.55, 'max UMI family size = '+str(max(x)), transform=ax1.transAxes, fontsize = 14, horizontalalignment='right')

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_SSCS_UMI_family_size_distribution_1_20_MUFS'+str(min_family_size)+'.pdf')


    ##### PLOT PROPORTION OF SSCS THAT ARE N #####
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, figsize=(10, 8))
    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15

    total_reads = 0
    for k in sorted(proportion_N_dict.keys()):
        v = proportion_N_dict[k]
        total_reads+=v

    cumulative = total_reads
    cumulative_dict = {}
    for k in sorted(proportion_N_dict.keys()):
        v = proportion_N_dict[k]
        cumulative_dict[k] = cumulative/total_reads
        cumulative-=v

    x = []
    y = []
    for k, v in cumulative_dict.items():
        x.append(k)
        y.append(v)

    ax1.scatter(x, y, color = blue3, lw = 2, zorder = 5)
    ax1.plot(x, y, color = blue3, lw = 2, zorder = 5)

    # CONFIGURING THE GRAPH APPEARANCE
    ax1.tick_params(axis='both',labelsize=axislabelfont)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    #     ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_major_locator(MultipleLocator(0.1))
    ax1.yaxis.set_major_locator(MultipleLocator(0.1))
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(which = 'minor', width=1, color = '#969696', length = 4)
    ax1.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    #Title and axis labels
    ax1.set_title('Proportion of SSCS read called as N \n'+sample_name, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('cumulative fraction of SSCS reads \n(with N proportion >= x)', fontsize = axislabelfont)
    ax1.set_xlabel('proportion of SSCS read called as N', fontsize = axislabelfont)
    ax1.set_ylim(-0.01, 1.01)
    ax1.set_xlim(-0.01, 1.01)

    ax1.grid(which = 'both', lw = 0.8, linestyle = ':', color = grey3)

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_SSCS_proportion_N_distribution_MUFS'+str(min_family_size)+'.pdf')

    return print('SSCS calling (BAM) done')

def SSCS_metrics_plots(SSCS_bam, sample_name, out_directory, min_family_size):
    infile_SSCS = SSCS_bam
    in_bam_SSCS = pysam.Samfile(infile_SSCS, "rb")
    bam_entry_SSCS = in_bam_SSCS.fetch(until_eof=True)

    read_number = 0

    SSCS_read1_forward = 0
    SSCS_read1_reverse = 0
    SSCS_read2_forward = 0
    SSCS_read2_reverse = 0

    SSCS_read1_forward_dict = {}
    SSCS_read1_reverse_dict = {}
    SSCS_read2_forward_dict = {}
    SSCS_read2_reverse_dict = {}

    # for line in new_reads_list:
    for line in bam_entry_SSCS:
        header = line.qname
        UMI_family_size = int((header.split(':'))[1])

        if line.is_reverse is True:
            if line.is_read1 is True:
                SSCS_read1_reverse+=1
                if UMI_family_size in SSCS_read1_reverse_dict.keys():
                    SSCS_read1_reverse_dict[UMI_family_size]+=1
                if UMI_family_size not in SSCS_read1_reverse_dict.keys():
                    SSCS_read1_reverse_dict[UMI_family_size]=1
            if line.is_read1 is False:
                SSCS_read2_reverse+=1
                if UMI_family_size in SSCS_read2_reverse_dict.keys():
                    SSCS_read2_reverse_dict[UMI_family_size]+=1
                if UMI_family_size not in SSCS_read2_reverse_dict.keys():
                    SSCS_read2_reverse_dict[UMI_family_size]=1

        if line.is_reverse is False:
            if line.is_read1 is True:
                SSCS_read1_forward+=1
                if UMI_family_size in SSCS_read1_forward_dict.keys():
                    SSCS_read1_forward_dict[UMI_family_size]+=1
                if UMI_family_size not in SSCS_read1_forward_dict.keys():
                    SSCS_read1_forward_dict[UMI_family_size]=1
            if line.is_read1 is False:
                SSCS_read2_forward+=1
                if UMI_family_size in SSCS_read2_forward_dict.keys():
                    SSCS_read2_forward_dict[UMI_family_size]+=1
                if UMI_family_size not in SSCS_read2_forward_dict.keys():
                    SSCS_read2_forward_dict[UMI_family_size]=1

        read_number+=1

    ####### Plot the distribution of reads by strand ##########
    plt.close('all')
    f, (ax1) = plt.subplots(1, 1, figsize=(10, 8))

    m_size = 100
    axisfont=15
    titlefont=15
    axislabelfont=15

    x = ['read1_forward', 'read2_reverse', 'read2_forward', 'read1_reverse']
    y = [SSCS_read1_forward, SSCS_read2_reverse, SSCS_read2_forward, SSCS_read2_reverse]

    ax1.bar(x, y, color = [c1, c2, c3, c4], lw = 2, zorder = 5)

    # CONFIGURING THE GRAPH APPEARANCE
    ax1.tick_params(axis='both',labelsize=axislabelfont)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)
        ax1.spines[axis].set_color('#969696')
    ax1.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
    ax1.xaxis.set_tick_params(width=1, color = '#969696', length = 6)

    #Only show the required axis lines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    #Title and axis labels
    ax1.set_title('Read strand/ direction distribution of SSCSs \n('+sample_name+')', fontsize = titlefont, fontweight='bold')
    ax1.set_xlabel('Read strand/ direction', fontsize = axislabelfont)
    ax1.set_ylabel('Number of SSCSs', fontsize = axislabelfont)

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_read_strand_direction_SSCSs_MUFS'+str(min_family_size)+'.pdf')

    ####### Plot the UMI family size distribution by read/ strand ##########
    def UMI_family_size_plot(UMI_family_sizes, read_strand_direction, ax, x_lim):
        # plot
        m_size = 100
        axisfont=15
        titlefont=15
        axislabelfont=15

        x = []
        y = []

        total_number_reads = 0
        total_number_read_1 = 0
        total_number_reads_family_size_more_than_1 = 0
        total_number_unique_molecules = 0
        total_number_unique_molecules_excluding_1 = 0


        for k, v in sorted(UMI_family_sizes.items()):
            x.append(int(k))
            y.append(int(v))
            number_reads = int(k)*int(v)
            total_number_reads += number_reads
            if int(k) == 1:
                total_number_read_1 += int(v)
            if int(k) > 1:
                total_number_reads_family_size_more_than_1 += number_reads
                total_number_unique_molecules_excluding_1 += int(v)
            total_number_unique_molecules += int(v)

        mean_UMI_family_size = total_number_reads/total_number_unique_molecules

        total_square_distance = 0
        for k in sorted(UMI_family_sizes.keys()):
            v = UMI_family_sizes[k]
            square_distance = v*((k-mean_UMI_family_size)**2)
            total_square_distance+=square_distance
        standard_deviation = np.sqrt(total_square_distance/total_number_unique_molecules)

        ax.plot(x, y, color = 'dodgerblue', lw = 2, zorder = 5)

        # CONFIGURING THE GRAPH APPEARANCE
        ax.tick_params(axis='both',labelsize=axislabelfont)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color('#969696')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if x_lim<30:
            ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_tick_params(width=1, color = '#969696', length = 6)
        ax.xaxis.set_tick_params(which = 'minor', width=1, color = '#969696', length = 4)
        ax.xaxis.set_tick_params(which = 'major', width=1, color = '#969696', length = 6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(0, x_lim)

        #Title and axis labels
        ax.set_title('UMI family size distribution ('+sample_name+'):\n'+str(read_strand_direction), fontsize = titlefont, fontweight='bold')
        ax.set_ylabel('number of read families', fontsize = axislabelfont)
        ax.set_xlabel('UMI family size', fontsize = axislabelfont)

        ax.text(0.9, 0.9, 'total reads = '+str(total_number_reads), transform=ax.transAxes, fontsize = 14, horizontalalignment='right')
        ax.text(0.9, 0.85, 'total reads with family size of 1 = '+str(total_number_read_1), transform=ax.transAxes, fontsize = 14, horizontalalignment='right')
        ax.text(0.9, 0.75, 'total unique molecules = '+str(total_number_unique_molecules), transform=ax.transAxes, fontsize = 14, horizontalalignment='right')
        ax.text(0.9, 0.7, 'total unique molecules (excl. family size 1) = '+str(total_number_unique_molecules_excluding_1), transform=ax.transAxes, fontsize = 14, horizontalalignment='right')

        ax.text(0.9, 0.6, 'mean UMI family size = '+str(round(mean_UMI_family_size, 2))+' (std = '+str(round(standard_deviation, 2))+')', transform=ax.transAxes, fontsize = 14, horizontalalignment='right')
        ax.text(0.9, 0.55, 'max UMI family size = '+str(max(x)), transform=ax.transAxes, fontsize = 14, horizontalalignment='right')

        return ax

    plt.close('all')
    f, (axes) = plt.subplots(1, 4, figsize=(40, 8))

    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]

    family_sizes = []
    for k in SSCS_read1_forward_dict.keys():
        family_sizes.append(k)
    for k in SSCS_read2_forward_dict.keys():
        family_sizes.append(k)
    for k in SSCS_read1_reverse_dict.keys():
        family_sizes.append(k)
    for k in SSCS_read2_reverse_dict.keys():
        family_sizes.append(k)

    x_lim = max(family_sizes)

    UMI_family_size_plot(SSCS_read1_forward_dict, 'read 1 forward', ax1, x_lim)
    UMI_family_size_plot(SSCS_read2_forward_dict, 'read 2 forward', ax2, x_lim)
    UMI_family_size_plot(SSCS_read1_reverse_dict, 'read 1 reverse', ax3, x_lim)
    UMI_family_size_plot(SSCS_read2_reverse_dict, 'read 2 reverse', ax4, x_lim)

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_SSCS_UMI_family_size_different_reads_strand_direction_MUFS'+str(min_family_size)+'.pdf')

    return print('metrics plots saved')

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mapped merged BAM file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument('--min-family-size', type=int, default=1, dest='min_family_size', help="Minimum number of reads required to form a consensus. [default = 1]")
    parser.add_argument('--threshold', type=float, default=0.9, dest='threshold', help="Minimum number of matching nucleotides required to form a consensus. [default = 0.9]")
    parser.add_argument('--max_N', type=float, default=1, dest='max_N', help="Maximum fraction of Ns permitted in a consensus [default = 1.0]")
    parser.add_argument('--min-mapping-quality', type=int, default=20, dest='min_mapping_qual', help="Minimum mapping quality (reads with mapping quality below this will be excluded) [default = 20]")
    parser.add_argument('--min-base-quality', type=int, default=10, dest='quality_threshold', help="Minimum per-base quality (bases with quality scores below this will not be included in consensus calling) [default = 10]")
    parser.add_argument("--outbam", action="store", dest="outbam", help="output SSCS BAM file", required=True)
    parser.add_argument("--unpaired-outbam", action="store", dest="unpaired_outbam", help="output BAM file for unpaired SSCS", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    merged_bam = o.infile
    unpaired_SSCS_bam = o.unpaired_outbam
    SSCS_bam = o.outbam
    sample_name = o.sample_name
    min_family_size = o.min_family_size
    threshold = o.threshold
    out_directory = o.out_directory
    max_N = o.max_N
    min_mapping_qual=o.min_mapping_qual
    quality_threshold = o.quality_threshold

    #Group reads by UMI
    grouped_read_dict = Group_reads_by_UMI(sample_name, merged_bam, out_directory, min_mapping_qual, min_family_size)
    #Call SSCS
    SSCS_calling_BAM(sample_name, grouped_read_dict, merged_bam, threshold, min_family_size, max_N, SSCS_bam, unpaired_SSCS_bam, out_directory, quality_threshold, version, min_mapping_qual)
    SSCS_metrics_plots(SSCS_bam, sample_name, out_directory, min_family_size)

    return print('SSCS calling complete')

if __name__ == "__main__":
	main()
