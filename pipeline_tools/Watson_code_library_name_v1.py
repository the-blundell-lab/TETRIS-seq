#!/usr/bin/env python

'''''
Watson code for retrieving the library name from the fastq file name.
Version 1.1 (May 2021)

Input:
    1) fastq file name

Outputs:
    library name (e.g. SLX-20124)

Usage:
Watson_code_library_name_v1.py --fastq fastq_file_name

'''''
import csv
from argparse import ArgumentParser
import sys

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--fastq", type=str, dest='fastq_name', help="fastq file for processing (containing UDI index in filename)", required=True)
    o = parser.parse_args()

    fastqname = o.fastq_name

    library_name = fastqname.split('.')[0]

    return library_name

if __name__ == "__main__":
    sys.stdout.write(main())
	# main()
