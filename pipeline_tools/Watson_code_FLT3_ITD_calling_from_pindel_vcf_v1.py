#!/usr/bin/env python

'''''
Watson code for creating a simplified FLT-ITD output file from pindel VCF (that just looks for FLT3 ITD in the FLT3 ITD domain)
Version 1.0 (November 2021)

Input:
    1) pindel VCF file
    2) output file name
    3) sample name

Outputs:
    1) FLT3 ITD calls

Usage:
Watson_code_FLT3_ITD_calling_from_pindel_vcf_v1.py  --pindel_VCF pindel VCF file --outfile FLT3 ITD csv file --sample_name --min_depth

'''''
version = '1.0'

from argparse import ArgumentParser
import pysam
import sys
import gzip
import numpy as np
from array import array
import timeit
import time
from datetime import date
import os
import shlex
import subprocess
import csv
import pandas as pd

def create_dataframe(pindel_VCF, min_depth, DCS_or_SSCS):

    FLT3_ITD_results = {}

    FLT3_ITD_domain = (28608024, 28608351) #exons 14-15

    with open(pindel_VCF) as vcf:
        reader = csv.reader(vcf, delimiter = '\t')
        row_count=0
        for row in reader:
            if row[0][0]!='#':
                chromosome = row[0]
                position = row[1]
                ID = row[2]
                REF = row[3]
                ALT = row[4]
                QUAL = row[5]
                FILTER = row[6]
                INFO = row[7]
                END = int(INFO.split(';')[0].split('=')[1])
                HOMLEN = int(INFO.split(';')[1].split('=')[1])
                SVLEN = int(INFO.split(';')[2].split('=')[1])
                SVTYPE = INFO.split(';')[3]
                NTLEN = int(INFO.split(';')[4].split('=')[1])
                SVLEN_PLUS_NTLEN = SVLEN+NTLEN
                FORMAT = row[8]
                alt_depth = int(row[9].split(':')[1].split(',')[0])
                var_depth = int(row[9].split(':')[1].split(',')[1])
                total_depth = int(alt_depth+var_depth)
                if alt_depth==0:
                    if var_depth>0:
                        allelic_ratio = 1
                        VAF = 1
                    if var_depth ==0:
                        allelic_ratio = 0
                        VAF = 0
                if alt_depth!=0:
                    allelic_ratio = var_depth/alt_depth
                    VAF = var_depth/(var_depth+alt_depth)
                if int(total_depth) > int(min_depth):
                    if DCS_or_SSCS == 'SSCS':
                        if var_depth >1:
                            if FLT3_ITD_domain[0]<= int(position) <= FLT3_ITD_domain[1]:
                                FLT3_ITD_results[row_count]={'chromosome': 'chr'+chromosome, 'start': position, 'end': END, 'total_depth': total_depth,
                                                            'variant_depth': var_depth,
                                                             'VAF': VAF, 'FLT3 ITD allelic ratio': allelic_ratio,
                                                             'SVTYPE': SVTYPE, 'SVLEN': SVLEN, 'NTLEN': NTLEN, 'SVLEN+NTLEN': SVLEN_PLUS_NTLEN, 'REF': REF, 'ALT': ALT}
                    if DCS_or_SSCS == 'DCS':
                        if FLT3_ITD_domain[0]<= int(position) <= FLT3_ITD_domain[1]:
                            FLT3_ITD_results[row_count]={'chromosome': 'chr'+chromosome, 'start': position, 'end': END, 'total_depth': total_depth,
                                                        'variant_depth': var_depth,
                                                         'VAF': VAF, 'FLT3 ITD allelic ratio': allelic_ratio,
                                                         'SVTYPE': SVTYPE, 'SVLEN': SVLEN, 'NTLEN': NTLEN, 'SVLEN+NTLEN': SVLEN_PLUS_NTLEN, 'REF': REF, 'ALT': ALT}
        row_count+=1

    df = pd.DataFrame.from_dict(FLT3_ITD_results, orient = 'index')

    return df

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--pindel_VCF", action="store", dest="pindel_VCF", help="pindel VCF file", required=True)
    parser.add_argument("--outfile", type=str, dest='outfile', help="name of output config file", required=True)
    parser.add_argument("--min_depth", type=int, dest='min_depth', help="minimum total depth of positions to look at", required=True)
    parser.add_argument("--DCS_or_SSCS", type=str, dest='DCS_or_SSCS', help="whether DCS or SSCS file", required=True)
    parser.add_argument('--sample_name', type=str, dest='sample_name', help="sample_name")
    o = parser.parse_args()

    pindel_VCF = o.pindel_VCF
    outfile = o.outfile
    DCS_or_SSCS = o.DCS_or_SSCS
    sample_name=o.sample_name
    min_depth = o.min_depth

    FLT3_ITD_df = create_dataframe(pindel_VCF, min_depth, DCS_or_SSCS)
    FLT3_ITD_df.to_csv(outfile, index = False)

    return print('FLT3 ITD calling complete for '+sample_name)

if __name__ == "__main__":
	main()
