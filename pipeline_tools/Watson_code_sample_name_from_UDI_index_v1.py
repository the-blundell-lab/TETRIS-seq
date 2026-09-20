#!/usr/bin/env python

'''''
Watson code for retrieving the sample name from a csv file linking sample index and sample name.
Version 1.1 (May 2021)

Input:
    1) csv file containing link between sample index and sample name
    2) fastq file name

Outputs:
    sample name

Usage:
Watson_code_sample_name_from_UDI_index_v1.py  --samples_csv_file csv_file --fastq fastq_file_name

'''''
import csv
from argparse import ArgumentParser
import sys

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--samples_csv_file", action="store", dest="infile", help="csv file containing sample names in 1st column and UDI indexes in 2nd column", required=True)
    parser.add_argument("--fastq", type=str, dest='fastq_name', help="fastq file for processing (containing UDI index in filename)", required=True)
    o = parser.parse_args()

    csv_file = o.infile
    fastqname = o.fastq_name

    UDI = fastqname.split('.')[1]
    # print(UDI)

    UDI_dictionary = {}
    with open(csv_file) as csvfile:
        readreader = csv.reader(csvfile)
        for row in readreader:
            UDI_dictionary[row[1]]=row[0]

    sample_name = UDI_dictionary[UDI]

    # print(sample_name)

    return sample_name

if __name__ == "__main__":
    sys.stdout.write(main())
	# main()
