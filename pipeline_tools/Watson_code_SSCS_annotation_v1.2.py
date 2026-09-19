#!/usr/bin/env python

'''''
Watson code for creating a txt file combining the vcf information with the annovar annotations.
Version 1.2 (June 2021)

Input:
    1) SSCS variants and SNPs text file
    2) annovar annotated variants FILES

Outputs:
    1) txt file containing the annotated variants

Usage:
Watson_code_SSCS_annotation_v1.2.py  --variants_infile variants and SNPs text file --annotated_infile annovar annotated file
                                --min_family_size size minimum family size used for SSCS calling
                                --sample-name sample_name --outfile output txt file

'''''
version = '1.2'

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
import pandas as pd

csv.field_size_limit(100000000)

def genes_and_transcripts(gene_transcripts):
    with open(gene_transcripts) as csvfile:
        readreader = csv.reader(csvfile)
        gene_list = []
        transcripts_list = []
        for row in readreader:
            gene = row[0]
            transcript = row[1]
            NM_transcript = row[2]
            transcripts_list.append(NM_transcript)
            gene_list.append(gene)

    return gene_list, transcripts_list

def create_annotated_df(annotated_file, gene_list, transcripts_list):
    with open(annotated_file) as csvfile:
        readreader = csv.reader(csvfile)
        row_count = 0
        SSCS_ann_dict = {}

        for row in readreader:
            transcript_info = '.'
            exon = '.'
            cDNA = '.'
            protein_change = '.'
            cosmic_ID = '.'
            cosmic_occurence = '.'
            cosmic_haem_occurences = 0
            cosmic_total = 0

            if row_count > 2:
                chromosome = row[0]
                start = int(row[1])
                end = int(row[2])
                ref = row[3]
                alt = row[4]
                intronic_exonic = row[5]
                gene = row[6]
                if ';' in gene:
                    gene_options = gene.split(';')
                    for i in gene_options:
                        if i in gene_list:
                            gene = i
                            break
                if gene in gene_list:
                    non_exonic_changes = row[7]
                    exonic_function = row[8]

                    if exonic_function == '.':
                        change = non_exonic_changes.replace(':', ';')
                        changes = change.split(';')
                        n = 0
                        for i in changes:
                            if i in transcripts_list:
                                pos = n
                                transcript_info = changes[pos]
                                cDNA = changes[pos+1]
                                break
                            n+=1

                        cosmic_coding = row[10]
                        cosmic_noncoding = row[11]

                        if cosmic_noncoding != '.':
                            cosmic_haem_occurences = 0
                            cosmic_total = 0
                            cosmic_details = cosmic_noncoding.split(';')
                            cosmic_ID = cosmic_details[0].split('=')[1]
                            cosmic_occurence = cosmic_details[1].split('=')[1]
                            cosmic_sites = cosmic_occurence.split(',')
                            for site in cosmic_sites:
                                if site.split('(')[1]=='haematopoietic_and_lymphoid_tissue)':
                                    cosmic_haem_occurences = int(site.split('(')[0])
                                occurence = int(site.split('(')[0])
                                cosmic_total+=occurence

                        if cosmic_coding != '.':
                            cosmic_haem_occurences = 0
                            cosmic_total = 0
                            cosmic_details = cosmic_coding.split(';')
                            cosmic_ID = cosmic_details[0].split('=')[1]
                            cosmic_occurence = cosmic_details[1].split('=')[1]
                            cosmic_sites = cosmic_occurence.split(',')
                            for site in cosmic_sites:
                                if site.split('(')[1]=='haematopoietic_and_lymphoid_tissue)':
                                    cosmic_haem_occurences = int(site.split('(')[0])
                                occurence = int(site.split('(')[0])
                                cosmic_total+=occurence

                    else:
                        AA_change = row[9]
                        transcript = AA_change.split(':')
                        n = 0
                        for i in transcript:
                            if i in transcripts_list:
                                position = n
                                transcript_info = transcript[n]
                                exon = transcript[n+1]
                                cDNA = transcript[n+2]
                                if len(transcript)>4:
                                    if len(transcript)>(n+3):
                                        protein_change = transcript[n+3].split(',')[0]
                                    else:
                                        protein_change = '.'
                                else:
                                    protein_change = '.'
                                break
                            n+=1

                        cosmic_coding = row[10]
                        cosmic_noncoding = row[11]

                        if cosmic_noncoding != '.':
                            cosmic_haem_occurences = 0
                            cosmic_total = 0
                            cosmic_details = cosmic_noncoding.split(';')
                            cosmic_ID = cosmic_details[0].split('=')[1]
                            cosmic_occurence = cosmic_details[1].split('=')[1]
                            cosmic_sites = cosmic_occurence.split(',')
                            for site in cosmic_sites:
                                if site.split('(')[1]=='haematopoietic_and_lymphoid_tissue)':
                                    cosmic_haem_occurences = int(site.split('(')[0])
                                occurence = int(site.split('(')[0])
                                cosmic_total+=occurence

                        if cosmic_coding != '.':
                            cosmic_haem_occurences = 0
                            cosmic_total = 0
                            cosmic_details = cosmic_coding.split(';')
                            cosmic_ID = cosmic_details[0].split('=')[1]
                            cosmic_occurence = cosmic_details[1].split('=')[1]
                            cosmic_sites = cosmic_occurence.split(',')
                            for site in cosmic_sites:
                                if site.split('(')[1]=='haematopoietic_and_lymphoid_tissue)':
                                    cosmic_haem_occurences = int(site.split('(')[0])
                                occurence = int(site.split('(')[0])
                                cosmic_total+=occurence

                    exac_all = row[12]
                    exac_afr = row[13]
                    exac_amr = row[14]
                    exac_eas = row[15]
                    exac_fin = row[16]
                    exac_nfe = row[17]
                    exac_oth = row[18]
                    exac_sas = row[19]
                    gnomad_all = row[20]
                    gnomad_afr = row[21]
                    gnomad_amr = row[22]
                    gnomad_asj = row[23]
                    gnomad_eas = row[24]
                    gnomad_fin = row[25]
                    gnomad_nfe = row[26]
                    gnomad_oth = row[27]
                    clin_allele_id = row[28]
                    clin_dn = row[29]
                    clin_disdb = row[30]
                    clin_rev = row[31]
                    clin_sig = row[32]


                    SSCS_ann_dict[(chromosome, start, end, ref, alt)]=(gene, transcript_info, intronic_exonic, exon, cDNA, protein_change, exonic_function,
                                                                      cosmic_ID, cosmic_total, cosmic_haem_occurences, cosmic_occurence, exac_all, exac_afr,
                                                                      exac_amr, exac_eas, exac_fin, exac_nfe, exac_oth, exac_sas, gnomad_all,
                                                                      gnomad_afr, gnomad_amr, gnomad_asj, gnomad_eas, gnomad_fin, gnomad_nfe, gnomad_oth,
                                                                     clin_allele_id, clin_dn, clin_disdb, clin_rev, clin_sig)

            row_count+=1

    SSCS_ann_df = pd.DataFrame.from_dict(SSCS_ann_dict, orient = 'index', columns = ['gene', 'transcript', 'intronic_exonic', 'exon', 'cdNA', 'AA_change', 'exonic_function',
                                                                                  'cosmic_ID', 'cosmic_total', 'cosmic_haem_lymphoid', 'cosmic_sites',
                                                                                  'exac_all', 'exac_afr', 'exac_amr', 'exac_eas', 'eax_fin', 'exac_nfe',
                                                                                  'exac_oth', 'exac_sas', 'gnomad_all', 'gnomad_afr', 'gnomad_amr', 'gnomad_asj',
                                                                                  'gnomad_eas', 'gnomad_fin', 'gnomad_nfe', 'gnomad_oth', 'clinvar_allele_id',
                                                                                  'clinvar_dn', 'clinvar_disdb', 'clinvar_rev', 'clin_sig'])
    SSCS_ann_df = SSCS_ann_df.reset_index()

    return SSCS_ann_df

def create_variants_df(variants_file):
    with open(variants_file) as csvfile:
        readreader = csv.reader(csvfile, delimiter = '\t')
        row_count = 0
        SSCS_variants_dict = {}
        for row in readreader:
            if row_count > 5:
    #             print(row)
                chromosome = row[0]
                start = int(row[1])
                end = int(row[2])
                ref = row[3]
                alt = row[4]
                RSID = row[5]
                MAF_1000_genomes = row[6]
                total_depth = int(row[7])
                variant_depth = int(row[8])
                if variant_depth >0:
                    VAF = float(row[9])
                    variant_type = row[10]
                    REF_mean_UMI_family_size = float(row[11])
                    variant_mean_UMI_family_size = float(row[12])
                    mean_position_in_read = float(row[13])
                    std_position_in_read = float(row[14])
                    mean_base_quality = float(row[15])
                    std_base_quality = float(row[16])
                    mean_mapq = float(row[17])
                    REF_depth_by_strand = row[18]
                    variant_depth_by_strand = row[19]
                    strand_bias_fisher_p_value = float(row[20])
                    strand_bias_fisher_p_value_phred = float(row[21])

                    SSCS_variants_dict[((chromosome, start, end, ref, alt))] = (total_depth, variant_depth, VAF, variant_type, REF_mean_UMI_family_size,
                                                                              variant_mean_UMI_family_size, mean_position_in_read, std_position_in_read,
                                                                                mean_base_quality, std_base_quality,
                                                                              mean_mapq, REF_depth_by_strand, variant_depth_by_strand, strand_bias_fisher_p_value,
                                                                              strand_bias_fisher_p_value_phred, RSID, MAF_1000_genomes)

            row_count+=1

    SSCS_variants_df = pd.DataFrame.from_dict(SSCS_variants_dict, orient = 'index', columns = ['total_depth', 'variant_depth', 'VAF', 'variant_type', 'REF_mean_UMI_family_size',
                                                                                            'variant_mean_UMI_family_size', 'mean_position_in_read', 'std_position_in_read',
                                                                                               'mean_base_quality', 'std_base_quality',
                                                                                            'mean_mapq', 'REF_depth_by_strand', 'variant_depth_by_strand', 'strand_bias_fisher_p_value',
                                                                                            'strand_bias_fisher_p_value_phred', 'RSID', 'MAF_1000_genomes'])
    SSCS_variants_df = SSCS_variants_df.reset_index()

    return SSCS_variants_df

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--variants_infile", action="store", dest="variants_infile", help="input variants and SNPs txt file", required=True)
    parser.add_argument("--annotated_infile", action="store", dest="annotated_infile", help="input annovar annotated file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--gene_transcripts", action="store", dest="gene_transcripts", help="csv file containing gene transcripts used", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output annotated file", required=True)
    parser.add_argument('--min_family_size', type=int, default=1, dest='min_family_size', help="min UMI family size that was used for SSCS calling")
    o = parser.parse_args()

    variants_file = o.variants_infile
    annotated_file = o.annotated_infile
    sample_name = o.sample_name
    gene_transcripts = o.gene_transcripts
    outfile = o.outfile
    min_family_size = o.min_family_size

    gene_list, transcripts_list = genes_and_transcripts(gene_transcripts)

    #create dataframes for the annotated file and for the variants file
    annotated_df = create_annotated_df(annotated_file, gene_list, transcripts_list)
    variants_df = create_variants_df(variants_file)

    #merge the dataframes
    total_df = pd.merge(variants_df, annotated_df, how="left", on=["index"])
    total_df[['chromosome', 'start', 'end', 'REF', 'ALT']] = pd.DataFrame(total_df['index'].tolist(), index=total_df.index)

    #reorder the columns in the dataframe
    final_df = total_df[['chromosome', 'start', 'end', 'REF', 'ALT', 'total_depth', 'variant_depth', 'VAF',
                     'intronic_exonic', 'variant_type', 'gene', 'transcript', 'exon', 'cdNA',
       'AA_change', 'exonic_function', 'RSID', 'MAF_1000_genomes', 'cosmic_ID', 'cosmic_total',
       'cosmic_haem_lymphoid', 'cosmic_sites', 'REF_mean_UMI_family_size', 'variant_mean_UMI_family_size',
       'mean_position_in_read', 'std_position_in_read', 'mean_base_quality', 'std_base_quality', 'mean_mapq',
       'REF_depth_by_strand', 'variant_depth_by_strand',
       'strand_bias_fisher_p_value', 'strand_bias_fisher_p_value_phred',
       'exac_all', 'exac_afr',
       'exac_amr', 'exac_eas', 'eax_fin', 'exac_nfe', 'exac_oth', 'exac_sas',
       'gnomad_all', 'gnomad_afr', 'gnomad_amr', 'gnomad_asj', 'gnomad_eas',
       'gnomad_fin', 'gnomad_nfe', 'gnomad_oth', 'clinvar_allele_id',
       'clinvar_dn', 'clinvar_disdb', 'clinvar_rev', 'clin_sig']]

    #replace non-values with '.' and sort by VAF (descending)
    final_df = final_df.fillna('.')

    #SNPs in the panel
    SNPs = ['rs10789158', 'rs3916765', 'rs1364429', 'rs2286510']

    final_df = final_df[(final_df['RSID'].isin(SNPs)) | (final_df['gene']!='.')] #remove anything that isn't in a gene in the panel or one of the targeted SNPs
    final_df = final_df.sort_values(by=['VAF'], ascending = False)

    #save the dataframe to txt
    final_df.to_csv(outfile, index = False, sep = '\t')

    return print('annotated csv file creation complete')

if __name__ == "__main__":
	main()
