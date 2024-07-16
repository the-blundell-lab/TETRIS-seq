#!/usr/bin/env python

'''''
Watson code for filtering SSCS BAM files by UMI family size.
Version 1.1 (December 202O)

Input:
    1) SSCS BAM file
    2) sample name
    3) minimum UMI family size

Outputs:
    BAM and fastq files:
        1) paired-end BAM containing SSCSs
        2) interleaved fastq file containing SSCSs
        3) BAM file containing unpaired SSCSs (e.g. if the partner is filtered out due to too small UMI family size)
    Metrics files:
        1) metrics txt file containing info on number of reads filtered out

Usage:
Watson_code_filter_SSCS_1.1.py  --infile SSCS_BAM --sample-name sample_name --min-family-size NUMBER(default = 1, range 1-no max)
                                --out-directory directory to save files in
'''''
version = '1.1'

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
from Bio.Alphabet import IUPAC
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

def filter_SSCS_bam(SSCS_bam, sample_name, min_family_size, out_directory, version):
    in_bam = pysam.Samfile(SSCS_bam, "rb")
    bam_entry = in_bam.fetch(until_eof=True)

    new_header = {}
    heading = in_bam.header
    new_header['HD'] = heading['HD']
    new_header['SQ'] = heading['SQ']
    new_header['RG'] = heading['RG']
    new_PG = [{'ID': 'Watson_code_filter_SSCS_v'+str(version),
    'PN': 'Watson_code_SNV_panel_v'+str(version),
    'VN': version,
    'CL': 'Watson_code_filter_SSCS_'+str(version)+'.py --infile '+SSCS_bam+' --sample-name '+sample_name+' --min-family-size '+str(min_family_size)+ '--out-directory '+out_directory}]
    new_header['PG'] = new_PG

    filtered_bam = SSCS_bam.replace('.bam', '_filtered_min_UMI_family_size_'+str(min_family_size)+'.bam')
    # unpaired_filtered_bam = SSCS_bam.replace('.bam', '_filtered_min_UMI_family_size_'+str(min_family_size)+'_unpaired.bam')
    filtered_fastq = SSCS_bam.replace('bam.bam', 'filtered_min_UMI_family_size_'+str(min_family_size)+'.fq.gz')

    #CREATE NEW BAM AND FASTQ FILES
    out_bam = pysam.Samfile(filtered_bam, "wb", header=new_header) #open a new BAM file to write the SSCS to
    # unpaired_out_bam = pysam.Samfile(unpaired_filtered_bam, "wb", header=new_header) #open a new BAM file to write the SSCS to
    fastq = gzip.open(filtered_fastq, "wt") #interleaved fastq file (r1 followed by r2 etc..

    paired_counter = 0
    unpaired_counter = 0
    filtered_out_counter = 0
    counter = 0

    start_time = time.time()
    for line in bam_entry:
        qname = line.qname #query name
        UMI = qname.split(':')[0]
        pair1_family_size = int(qname.split(':')[1])
        pair2_family_size = int(qname.split(':')[2])
        UMI_family = int(qname.split(':')[3])
        sequence = line.seq
        flag = line.flag
        cigar = line.cigartuples
        template_length = line.template_length
        chromosome = line.rname
        coordinate = line.pos
        partner_chromosome = line.mrnm
        partner_coordinate = line.next_reference_start
        mapq = line.mapping_quality
        quality = line.query_qualities
        cD_tag = line.get_tag('cD')
        cM_tag = line.get_tag('cM')
        cE_tag = line.get_tag('cE')
        cd_tag = line.get_tag('cd')
        ce_tag = line.get_tag('ce')
        RG_tag = line.get_tag('RG')
        QS_tag = line.get_tag('QS')
        if line.is_read1 is True:
            read_pair = 'read1'
        if line.is_read2 is True:
            read_pair = 'read2'
        if line.is_reverse is True:
            strand = 'reverse'
        if line.is_reverse is False:
            strand = 'forward'

        new_read = pysam.AlignedRead() #Create bam file reads for the new SSCS
        new_read.qname = qname
        new_read.seq = sequence
        new_read.flag = flag #contains info on whether read 1, read2, forward or reverse
        new_read.cigartuples = cigar
        new_read.rname = chromosome
        new_read.pos = coordinate
        new_read.mapq = mapq
        new_read.template_length = template_length
        new_read.query_qualities = quality
        new_read.mrnm = partner_chromosome
        new_read.next_reference_start = partner_coordinate
        new_read.set_tags([('cD', cD_tag), ('cM', cM_tag), ('cE', cE_tag), ('cd', cd_tag), ('ce', ce_tag), ('RG', sample_name), ('QS', QS_tag)])

        if (pair1_family_size >= min_family_size) and (pair2_family_size >= min_family_size):
            paired_counter+=1
            out_bam.write(new_read) #write the new SSCS to new bam file

            if read_pair == 'read1':
                qname = UMI+':'+str(pair1_family_size)+':'+str(pair2_family_size)+':'+str(UMI_family)+'/1'
            if read_pair == 'read2':
                qname = UMI+':'+str(pair1_family_size)+':'+str(pair2_family_size)+':'+str(UMI_family)+'/1'

            new_read.qname = qname
            fastq.write('@:%s\n%s\n+\n%s\n' % (new_read.qname, new_read.seq, new_read.query_qualities))

        if (pair1_family_size < min_family_size) and (pair2_family_size >= min_family_size):
            unpaired_counter+=1

        if (pair1_family_size >= min_family_size) and (pair2_family_size < min_family_size):
            unpaired_counter+=1

        if (pair1_family_size < min_family_size) and (pair2_family_size < min_family_size):
            filtered_out_counter+=1

        counter+=1
        if counter % 100000 == 0:
            print('SSCS reads processed by filter for sample '+sample_name+': ', counter)
            print('time for last 100,000 reads to be processed = %s seconds' % int(time.time() - start_time))
            start_time = time.time() #reset the timer so it can calculate the time for the next 100,000 reads

    out_bam.close() #close the new SSCS BAM file
    fastq.close() #close the fastq file
    in_bam.close() #close the open in BAM file

    metrics_file = open(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_filter_SSCS_metrics_min_UMI_family_size_'+str(min_family_size)+'.txt', 'w')
    metrics_file.write('sample name :\t'+ str(sample_name)+ '\n')
    metrics_file.write('date of analysis :\t'+ str(date_today)+ '\n')
    metrics_file.write('produced from code:\t' + 'Watson_code_filter_SSCS_calling: version '+str(version) + '\n')
    metrics_file.write('minimum UMI family size set to:\t' + str(min_family_size) + '\n')
    metrics_file.write('\n')
    metrics_file.write('reads passing filter with both pairs having UMI family size >= '+str(min_family_size)+': '+'\t' + str(paired_counter)+'\n')
    metrics_file.write('reads with one of read pair having UMI family size < '+str(min_family_size)+' (filtered out): '+'\t' + str(unpaired_counter)+'\n')
    metrics_file.write('reads with both pairs having UMI family size < '+str(min_family_size)+' (filtered out): '+'\t' + str(filtered_out_counter)+'\n')
    metrics_file.close()

    return 'filtered all reads'

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mapped merged BAM file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument('--min-family-size', type=int, default=1, dest='min_family_size', help="Minimum number of reads required to form a consensus. [default = 1]")
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    SSCS_bam = o.infile
    sample_name = o.sample_name
    min_family_size = o.min_family_size
    out_directory = o.out_directory

    filter_SSCS_bam(SSCS_bam, sample_name, min_family_size, out_directory, version)

    return print('SSCS filtering complete')

if __name__ == "__main__":
	main()
