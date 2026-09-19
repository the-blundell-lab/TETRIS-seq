#!/usr/bin/env python

'''''
Watson code for creating a text file from VarDict output
Version 1.2 (June 2021)

Input:
    1) VarDict java output file
    2) sample name

Outputs:
    1) Text file containing all variant positions where variant depth >0

Usage:
Watson_code_VarDict_to_TXT_1.2.py  --infile VarDictJava output file --sample-name sample_name
                                    --path-to-reference-genome path
                                    --variants_txt txt file to save the variants

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
from ast import literal_eval
import csv
import scipy
from scipy import stats

csv.field_size_limit(100000000)

today = date.today()
date_today = today.strftime("%d/%m/%Y")
VCF_date_today = today.strftime("%Y%m%d")

def fisher_exact(ref_F, ref_R, alt_F, alt_R):
    oddsratio, pvalue = stats.fisher_exact([[ref_F, ref_R], [alt_F, alt_R]])
    Phred_scale = -10*np.log10(pvalue)
    return pvalue, Phred_scale

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input VarDictJava file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--path-to-reference-genome", type=str, dest='reference_genome', help="path to reference genome", required=True)
    parser.add_argument("--variants_txt", action="store", dest="variants_txt", help="output file for text file of variants", required=True)
    o = parser.parse_args()

    vardictfile = o.infile
    sample_name = o.sample_name
    reference_genome = o.reference_genome
    variants_txt = o.variants_txt

    with open(vardictfile, 'r') as vcffile:
        read_reader = csv.reader(vcffile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0

        ####### CREATE TEXT FILE #######
        TEXT = open(variants_txt, 'w')
        TEXT.write('##fileDate='+str(VCF_date_today)+'\n')
        TEXT.write('##source=Watson_code_VarDict_to_TXT_v'+str(version)+'\n')
        TEXT.write('##reference='+str(reference_genome)+'\n')
        TEXT.write('##sample='+sample_name+'\n')
        TEXT.write('chromosome\t'+'position\t'+'end_position\t'+'REF\t'+'ALT\t'+'total_depth\t'+'variant_depth\t'+\
                  'VAF\t'+'variant_type\t'+'mean_position_in_read\t'+'std_position_in_read\t'+'mean_base_quality\t'+'std_base_quality\t'+'mean_mapq\t'+\
                   'REF_depth_by_strand\t'+'variant_depth_by_strand\t'+'strand_bias_fisher_p_value\t'+'strand_bias_fisher_p-value_phred\t'+\
                  'ratio_high_qual_to_low_qual_reads\t'+'VAF_for_high_qual_reads\t'+'adjusted_VAF_for_indels_due_to_local_realignment\t'+'number_bases_to_shift_3_prime_for_deletions_due_to_alterative_alignment\t'+\
                  'microsatellite\t'+'microsatellite_unit_length\t'+'average_number_mismatches_for_reads_containing_variant\t'+\
                   'high_qual_variant_reads\t'+'high_qual_coverage\t'+'5prime_ref_seq\t'+'3prime_ref_seg\t'+\
                  'duplication_rate_fraction\t'+'SV_split_reads_pairs_clusters\n')

        start_time = time.time()
        reset_time = time.time()

        for row in read_reader:
            chromosome = 'chr'+row[2]
            start = int(row[3])
            end = int(row[4])
            ref = row[5]
            alt = row[6]
            total_depth = row[7]
            variant_depth = int(row[8])

            ref_forward_depth = int(row[9])
            ref_rev_depth = int(row[10])
            REF_depth_by_strand = str(ref_forward_depth)+':'+str(ref_rev_depth)
            alt_forward_depth = int(row[11])
            alt_rev_depth = int(row[12])
            VAR_depth_by_strand = str(alt_forward_depth)+':'+str(alt_rev_depth)

            strand_bias_p_value, strand_bias_phred = fisher_exact(ref_forward_depth, ref_rev_depth, alt_forward_depth, alt_rev_depth)

            genotype = row[13]
            VAF = row[14]
            bias = row[15]
            mean_pos_read = row[16]
            mean_pos_read_std = row[17]
            mean_base_quality = row[18]
            mean_base_quality_std = row[19]
            mapq = row[20]
            qratio = row[21]
            hifreq = row[22]
            adj_VAF = row[23]
            shift3 = row[24]
            msi = row[25]
            msint = row[26]
            nm = row[27]
            hicnt = row[28]
            hicov = row[29]
            fivep_flank = row[30]
            threep_flank = row[31]
            position_description = row[32]
            vartype = row[33]
            duprate = row[34]
            sv_splits_pairs_clusters = row[35]

            if variant_depth !=0:
                TEXT.write(chromosome+'\t'+str(start)+'\t'+str(end)+'\t'+str(ref)+'\t'+str(alt)+'\t'+str(total_depth)+'\t'+str(variant_depth)+'\t'+
                str(VAF)+'\t'+str(vartype)+'\t'+str(mean_pos_read)+'\t'+str(mean_pos_read_std)+'\t'+str(mean_base_quality)+'\t'+str(mean_base_quality_std)+'\t'+
                str(mapq)+'\t'+str(REF_depth_by_strand)+'\t'+str(VAR_depth_by_strand)+'\t'+str(strand_bias_p_value)+'\t'+
                str(strand_bias_phred)+'\t'+str(qratio)+'\t'+str(hifreq)+'\t'+str(adj_VAF)+'\t'+str(shift3)+'\t'+str(msi)+'\t'+str(msint)+'\t'+
                str(nm)+'\t'+str(hicnt)+'\t'+str(hicov)+'\t'+str(fivep_flank)+'\t'+str(threep_flank)+'\t'+str(duprate)+'\t'+str(sv_splits_pairs_clusters)+'\n')

            row_count+=1
            if row_count % 10000 == 0:
                time_now = time.time()
                print('calling variants in '+chromosome+' in sample '+sample_name)
                print('10000 positions called in sample '+sample_name+' in %s seconds' % int(time_now - reset_time))
                print('total '+str(row_count)+' positions called in sample '+sample_name+' in '+str(round((time_now - start_time)/60, 2))+' minutes')
                reset_time = time.time()

    TEXT.close()

    return print('Text file creation complete')

if __name__ == "__main__":
	main()
