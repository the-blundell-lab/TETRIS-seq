#!/usr/bin/env python

'''''
Watson code for converting mapped BAM to unmapped BAM
Version 1.0 (August 2021)

Input:
    1) mapped BAM (paired end)
    2) unmapped BAM name
    3) sample name

Outputs:
    BAM file:
        1) paired-end unmapped BAM (reverse reads are reverse complemented compared to in the mapped BAM)

Usage:
Watson_code_convert_mapped_to_unmapped_BAM_v1.py --infile mapped_BAM --unmapped_BAM unmapped_BAM --sample-name sample name

'''''
version = '1.0'

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

def create_unmapped_BAM(mapped_bam, unmapped_bam, sample_name):

    def reverse(quality_scores):
            return quality_scores[::-1]

    in_bam = pysam.Samfile(mapped_bam, "rb")
    mappedbam = in_bam.fetch(until_eof=True)

    new_header = {}
    heading = in_bam.header
    new_header['HD'] = heading['HD']
    new_header['SQ'] = heading['SQ']
    new_header['RG'] = heading['RG']
    new_PG = [{'ID': 'Watson_code_unmap_SSCS_bam'+str(version),
    'PN': 'Watson_code_convert_mapped_to_unmapped_BAM_v'+str(version)+'.py',
    'VN': version,
    'CL': 'Watson_code_convert_mapped_to_unmapped_BAM_v'+str(version)+'.py --infile '+mapped_bam+' --outbam '+ unmapped_bam+' --sample-name '+sample_name}]
    new_header['PG'] = new_PG

    out_bam = pysam.Samfile(unmapped_bam, "wb", header=new_header) #open a new BAM file to write to

    s_time = time.time() #start the timer (so can record how long the processes are taking)
    start_time = time.time() #start the timer (so can record how long the processes are taking)

    counter = 0
    for line in mappedbam:
        qname = line.qname
        sequence = line.seq
        cigar = line.cigartuples
        quality = line.query_qualities
        tags = line.get_tags()

        if line.is_read1 is True:
            flag = 77 #flag for unmapped read 1
        if line.is_read2 is True:
            flag = 141 #flag for unmapped read 2
        if line.is_reverse is True:
            strand = 'reverse'
            sequence = str(Seq(sequence).reverse_complement())
            quality = reverse(quality)

        new_read = pysam.AlignedRead() #Create bam file reads for the new SSCS
        new_read.qname = qname
        new_read.seq = sequence #reverse complemented if mapped read was on reverse strand
        new_read.flag = flag #contains info on whether read 1, read2
        new_read.cigartuples = None #shouldn't be present in an unmapped BAM
        new_read.rname = -1 #chromosome shouldn't be present in an unmapped BAM
        new_read.pos = -1 #coordinate shouldn't be present in an unmapped BAM
        new_read.mapq = 0 #mapping quality shouldn't be present in an unmapped BAM
        new_read.template_length = 0 #template length shouldn't be present in an unmapped BAM
        new_read.query_qualities = quality #reversed if the mapped read was on reverse strand
        new_read.mrnm = -1 #partner chromosome shouldn't be present if an unmapped BAM
        new_read.next_reference_start = -1 #partner coordinate shouldn't be present if an unmapped BAM
        new_read.set_tags(tags) #keep the tags from the unmapped BAM

        out_bam.write(new_read) #write the new BAM file

        counter+=1
        if counter % 10000 == 0:
            print('Converting mapped BAM to unmapped BAM for sample '+sample_name+': reads converted = ', counter)
            print('time for last 10,000 reads to be converted = %s seconds' % int(time.time() - start_time))
            print('total time elapsed = '+str(int(int(time.time() - s_time)/60))+' minutes')
            start_time = time.time() #reset the timer so it can calculate the time for the next 100 consensuses

    out_bam.close()

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mapped BAM file", required=True)
    parser.add_argument("--outbam", action="store", dest="outbam", help="output unmapped BAM file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample", required=True)
    o = parser.parse_args()

    mapped_bam = o.infile
    unmapped_bam = o.outbam
    sample_name = o.sample_name

    starting_time = time.time() #start the timer (so can record how long the processes are taking)
    create_unmapped_BAM(mapped_bam, unmapped_bam, sample_name)

    return print('Creation of unmapped BAM from mapped BAM completed in '+str(int(int(time.time() - starting_time)/60))+' minutes')

if __name__ == "__main__":
	main()
