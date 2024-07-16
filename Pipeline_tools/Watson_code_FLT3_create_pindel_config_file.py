#!/usr/bin/env python

'''''
Watson code for creating a pindel configuration text file
Version 1.0 (November 2021)

Input:
    1) mapped merged BAM name (either SSCS or DCS)
    2) pindel output file name
    3) sequencing insert size

Outputs:
    1) pindel configuration text file

Usage:
Watson_code_FLT3_create_pindel_config_file.py  --infile mapped_merged_BAM --outfile pindel output file name
                                        --insert_size (e.g. 300) --SSCS_or_DCS whether SSCS or DCS

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


def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mapped merged BAM file (SSCS or DCS)", required=True)
    parser.add_argument("--outfile", type=str, dest='outfile', help="name of output config file", required=True)
    parser.add_argument('--insert_size', type=int, default=300, dest='insert_size', help="Sequencing insert size (default = 300)")
    parser.add_argument('--sample_name', type=str, dest='sample_name', help="sample_name")
    o = parser.parse_args()

    bam_file = o.infile
    pindel_config_file = o.outfile
    insert_size=o.insert_size
    sample_name=o.sample_name

    config_file = open(pindel_config_file, 'w')
    config_file.write(str(bam_file)+'\t'+str(insert_size)+'\t'+str(sample_name))
    config_file.close()

    return print('Pindel configuration file for '+sample_name+' created')

if __name__ == "__main__":
	main()
