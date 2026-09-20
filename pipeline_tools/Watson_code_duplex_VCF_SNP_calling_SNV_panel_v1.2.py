#!/usr/bin/env python

'''''
Watson code for creating a VCF file from samtools mpileup and for creating a txt file containing all variants and SNP sites in 1000 genomes covered by panel.
Version 1.2 (June 2021)

Input:
    1) samtools mpileup txt file
    2) sample name

Outputs:
    1) VCF file containing all positions in the targeted panel
    2) Text file containing all variant positions and all SNP positions from 1000 genomes covered by panel

Usage:
Watson_code_duplex_VCF_SNP_calling_1.2.py  --infile mpileup --sample-name sample_name --path-to-reference-genome path
                                    --dbSNP_shelve_dictionary path to dbSNP shelve dictionary --ECS_type raw, SSCS or DCS
                                    --mapq_filter --read_pos_filter --bias_filter
                                    --min-family-size that was used for SSCS calling --out-directory directory to save files in
                                    --variants_txt txt file to save the variants --variants_vcf vcf file to save the variants

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

csv.field_size_limit(100000000)

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

def dbSNP_dictionary(dbSNP_directory, out_directory, sample_name):

    chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    dbSNP_dict = shelve.open(out_directory+'/dbSNP_dictionary.db', flag = 'n')

    for chro in chromosomes_list:
        print('adding '+chro+' dbSNP to dictionary for sample '+sample_name)
        with open(dbSNP_directory+'/UCSC_dbSNP153common_hg19_'+chro+'.interval', 'r') as txtfile:
            read_reader = csv.reader(txtfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
            row_count = 0
            for row in read_reader:
                if row_count > 0:
                    chrom = row[0]
                    start = int(row[1])
                    stop = int(row[2])
                    ref = row[4]
                    altcount = int(row[5])
                    alts = (list(row[6].split(',')))
                    alts_list = [x for x in alts if x]
                    RSID = row[3]
                    MAF_1000_genomes = float(row[9].split(',')[0])

                    key = str((chrom, start+1, ref))#the dbSNP interval file is 0-based coordinates, so add 1
                    alt_dict = {}
                    for alt in alts_list:
                        alt_dict[alt]=(RSID, MAF_1000_genomes)
                        dbSNP_dict[key]=alt_dict
                row_count+=1

    return dbSNP_dict

def fisher_exact(ref_F, ref_R, alt_F, alt_R):
    oddsratio, pvalue = stats.fisher_exact([[ref_F, ref_R], [alt_F, alt_R]])
    Phred_scale = -10*np.log10(pvalue)
    return pvalue, Phred_scale

def dbSNP_information(chromosome, position, reference, alternate, dbSNP_dict):
#     print('looking for dbSNP information')
    key = str((chromosome, int(position), reference))
    if key in dbSNP_dict.keys():
        dbSNP_results = dbSNP_dict[key]
        if alternate in dbSNP_results.keys():
            RSID = dbSNP_results[alternate][0]
            MAF = dbSNP_results[alternate][1]
            return RSID, MAF
        else:
            return '-', '-'
    else:
        return '-', '-'

def homozygous_SNP_position(chromosome, position, reference, dbSNP_dict):
    key = str((chromosome, int(position), reference))
    if key in dbSNP_dict.keys():
        dbSNP_results = dbSNP_dict[key]
        for k, v in dbSNP_results.items():
            RSID = v[0]
            MAF = v[1]
            return RSID, MAF
    else:
        return '-', '-'

#Functions for creating sections of VCF file
def basics(chromosome, position, ID, ref, alt, filt):
    return chromosome+'\t'+str(position)+'\t'+ID+'\t'+ref.upper()+'\t'+alt.upper()+'\t'+filt

def info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq):
    return str('SAMPLE='+sample_name)+';TYPE='+str(variant_type)+';DP='+str(depth)+';END='+str(end)+';VD='+\
                str(variant_depth)+';AF='+str(round(float(allele_frequency), 6))+';UMI_ref='+str(round(float(UMI_ref), 2))+\
                ';UMI_alt='+str(round(float(UMI_alt), 2))+';RSID='+ID+';MAF_1000_genomes='+str(MAF_1000_genomes)+\
                ';REFBIAS='+str(ref_strand_bias)+';VARBIAS='+var_strand_bias+';PMEAN='+str(mean_pos_read)+\
                ';PSTD='+str(round(std_pos_read))+\
                ';SBF='+str(round(strand_fisher, 5))+';FS='+str(FS)+';MQ='+str(abs(mean_mapq))#+';'+five_prime+';'+three_prime

def format_section(GT, DP, VD, AD, AF, RD, ALD):
    return 'GT:DP:VD:AD:AF:RD:ALD'+'\t'+str(GT)+':'+str(DP)+':'+str(VD)+':'+str(AD)+':'+str(round(float(AF), 6))+':'+str(RD)+':'+str(ALD)

#Functions for creating sections of TEXT file
def text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref, alt, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq):
    return chromosome+'\t'+str(position)+'\t'+str(end_position)+'\t'+ref.upper()+'\t'+alt.upper()+'\t'+ID+'\t'+str(MAF_1000_genomes)+'\t'+\
            str(depth)+'\t'+str(variant_depth)+'\t'+str(allele_frequency)+'\t'+variant_type+'\t'+\
            str(round(float(UMI_ref)))+'\t'+str(round(float(UMI_alt)))+'\t'+str(mean_pos_read)+'\t'+str(round(std_pos_read))+'\t'+\
            str(abs(mean_mapq))+'\t'+str(ref_strand_bias)+'\t'+str(var_strand_bias)+'\t'+str(round(strand_fisher, 5))+'\t'+str(FS)+'\n'

def position_analysis(reference, depth, bases, UMI_family_size, read_positions, mapq):

    base_counter = {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'X': {}, '+':{}, '-':{}, '*': 0, '<>': 0, '?': 0}

    strand_counter = {'A': {'F': 0, 'R': 0}, 'G': {'F': 0, 'R': 0}, 'C': {'F': 0, 'R': 0}, 'T': {'F': 0, 'R': 0},
                      '+': {}, '-': {}, '*': {'F': 0, 'R': 0}, '<>': {'F': 0, 'R': 0}, 'X': {}}

    UMI_family_size_counter = {'A': [], 'G': [], 'C': [], 'T': [], '+': {}, '-': {}, '*': [], '<>': [], 'X': {}}

    position_counter = {'A': [], 'G': [], 'C': [], 'T': [], '+': {}, '-': {}, '*': [], '<>': [],'X': {}}

    mapq_counter = {'A': [], 'G': [], 'C': [], 'T': [], '+': {}, '-': {}, '*': [], '<>': [],'X': {}}

    N_counter = 0

    while bases != '':
        if bases[0] in ['.', ',']: #i.e. if the position is the reference:
#             print("bases[0] in ['.', ',']")
            ##### BASE = REFERENCE SEQUENCE ########
            if (len(bases) == 1) or (bases[1] not in ['+', '-']): #i.e. if it isn't the end of the string and isn't followed by insertion or deletion
                #add the base to the counter and the strand counter
                base_counter[reference]+=1
                if bases[0] == '.':
                    strand_counter[reference]['F']+=1
                if bases[0] == ',':
                    strand_counter[reference]['R']+=1
                bases = bases[1:] #then remove that base from the base string (ready to look at the next base)

                #add the UMI family size of that base to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                UMI_family_size_counter[reference].append(UMI_family_size[0])
                UMI_family_size = UMI_family_size[1:]

                #add the position in read of that base to the position counter
                position_counter[reference].append(read_positions[0])
                read_positions = read_positions[1:]

                #add the mapping quality score of that base to the mapping quality score counter
                mapq_counter[reference].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                mapq = mapq[1:]

            ##### INSERTION OR DELETION FOLLOWING A REFERENCE BASE ########
            elif bases[1] in ['+', '-']: #if the reference is followed by insertion or deletion marker
                if bases[1] == '+':
                    indel_marker = '+'
                if bases[1] == '-':
                    indel_marker = '-'

                #add the insertion to the insertion counter and its strand to the strand counter
                if (bases[3].isdigit() == False): #i.e. if the length of the indel is <10 (i.e. only bases[2] is a number)
                    indel_length = int(bases[2])
                    indel_sequence = bases[3:(3+indel_length)] #the sequence that is inserted or deleted will be from 3rd position to length of insertion
                    indel_sequence_record = indel_sequence.upper()
                    bases = bases[(3+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[3].isdigit() == True) and (bases[4].isdigit() == False): #i.e. if the length of the indel is 10-100:
                    indel_length = int(str(bases[2])+str(bases[3]))
                    indel_sequence = bases[4:(4+indel_length)] #the sequence that is inserted or deleted will be from 4th position to length of insertion
                    indel_sequence_record = indel_sequence.upper()
                    bases = bases[(4+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[4].isdigit() == True):
                    indel_length = int(str(bases[2])+str(bases[3])+str(bases[4]))
                    indel_sequence = bases[5:(5+indel_length)] #the sequence that is inserted or deleted will be from 5th position to length of insertion
                    indel_sequence_record = indel_sequence.upper()
                    bases = bases[(5+indel_length):] #then remove that indel from the base string (ready to look at the next base)

                if 'N' not in indel_sequence_record: #i.e. if the indel does not contain an N
                    if indel_sequence_record in base_counter[indel_marker].keys():
                        base_counter[indel_marker][indel_sequence_record]+=1
                    else:
                        base_counter[indel_marker][indel_sequence_record]=1

                    if any(i for i in indel_sequence if i.isupper()): #if the indel is in uppercase it is on forward strand
                        if indel_sequence_record in strand_counter[indel_marker].keys():
                            strand_counter[indel_marker][indel_sequence_record]['F']+=1
                        else:
                            strand_counter[indel_marker][indel_sequence_record]={'F': 1, 'R':0}
                    if any(i for i in indel_sequence if i.islower()): #if the indel is in lowercase it is on forward strand
                        if indel_sequence_record in strand_counter[indel_marker].keys():
                            strand_counter[indel_marker][indel_sequence_record]['R']+=1
                        else:
                            strand_counter[indel_marker][indel_sequence_record]={'F': 0, 'R':1}

                    #add the UMI family size of the insertion to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                    if indel_sequence_record in UMI_family_size_counter[indel_marker].keys():
                        UMI_family_size_counter[indel_marker][indel_sequence_record].append(UMI_family_size[0])
                    else:
                        UMI_family_size_counter[indel_marker][indel_sequence_record]=[UMI_family_size[0]]
                    UMI_family_size = UMI_family_size[1:] #only 1 UMI family size for the whole insertion so don't need to add more than 1

                    #add the position in read of the insertion to the position counter
                    if indel_sequence_record in position_counter[indel_marker].keys():
                        position_counter[indel_marker][indel_sequence_record].append(read_positions[0])
                    else:
                        position_counter[indel_marker][indel_sequence_record]=[read_positions[0]]
                    read_positions = read_positions[1:]

                    #add the mapping quality score of the insertion to the mapping quality score counter
                    if indel_sequence_record in mapq_counter[indel_marker].keys():
                        mapq_counter[indel_marker][indel_sequence_record].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                    else:
                        mapq_counter[indel_marker][indel_sequence_record]=[ord(mapq[0])]
                    mapq = mapq[1:]

                if 'N' in indel_sequence_record: #if the indel on that read contains an N, ignore that read's position
                    N_counter+=1
                    UMI_family_size = UMI_family_size[1:]
                    read_positions = read_positions[1:]
                    mapq = mapq[1:]


        elif bases[0] in ['*', '#']: #deletion of reference base on forward ('*') or reverse ('#') strand
            ##### SIMPLE REFERENCE DELETION ########
#             print("bases[0] in ['*', '#']")
            if (len(bases) == 1) or (bases[1] not in ['+', '-']):
                #add the deletion to the deletion counter and its strand to the strand counter
                base_counter['*']+=1
                if bases[0] == '*':
                    strand_counter['*']['F']+=1
                if bases[0] == '#':
                    strand_counter['*']['R']+=1
                bases = bases[1:]

                #add the UMI family size of the deletion to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                UMI_family_size_counter['*'].append(UMI_family_size[0])
                UMI_family_size = UMI_family_size[1:]

                #add the position in read of the deletion to the position counter
                position_counter['*'].append(read_positions[0])
                read_positions = read_positions[1:]

                #add the mapping quality score of the deletion to the mapping quality score counter
                mapq_counter['*'].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                mapq = mapq[1:]

            ##### INSERTION OR DELETION FOLLOWING A REFERENCE DELETION ######## - RECORD AS COMPLEX
            elif bases[1] in ['+', '-']: #if the reference deletion is followed by insertion or deletion marker
                if bases[1] == '+':
                    indel_marker = '+'
                if bases[1] == '-':
                    indel_marker = '-'

                #add the insertion to the insertion counter and its strand to the strand counter
                if (bases[3].isdigit() == False): #i.e. if the length of the indel is <10 (i.e. only bases[2] is a number)
                    indel_length = int(bases[2])
                    indel_sequence = bases[3:(3+indel_length)] #the sequence that is inserted or deleted will be from 3rd position to length of insertion
                    indel_sequence_record = '*'+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(3+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[3].isdigit() == True) and (bases[4].isdigit() == False): #i.e. if the length of the indel is 10-100:
                    indel_length = int(str(bases[2])+str(bases[3]))
                    indel_sequence = bases[4:(4+indel_length)] #the sequence that is inserted or deleted will be from 4th position to length of insertion
                    indel_sequence_record = '*'+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(4+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[4].isdigit() == True):
                    indel_length = int(str(bases[2])+str(bases[3])+str(bases[4]))
                    indel_sequence = bases[5:(5+indel_length)] #the sequence that is inserted or deleted will be from 5th position to length of insertion
                    indel_sequence_record = '*'+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(5+indel_length):] #then remove that indel from the base string (ready to look at the next base)

                if 'N' not in indel_sequence_record: #i.e. if the indel does not contain an N
                    if indel_sequence_record in base_counter['X'].keys():
                        base_counter['X'][indel_sequence_record]+=1
                    else:
                        base_counter['X'][indel_sequence_record]=1

                    if any(i for i in indel_sequence if i.isupper()): #if the indel is in uppercase it is on forward strand
                        if indel_sequence_record in strand_counter['X'].keys():
                            strand_counter['X'][indel_sequence_record]['F']+=1
                        else:
                            strand_counter['X'][indel_sequence_record]={'F': 1, 'R':0}
                    if any(i for i in indel_sequence if i.islower()): #if the indel is in lowercase it is on forward strand
                        if indel_sequence_record in strand_counter['X'].keys():
                            strand_counter['X'][indel_sequence_record]['R']+=1
                        else:
                            strand_counter['X'][indel_sequence_record]={'F': 0, 'R':1}

                    #add the UMI family size of the insertion to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                    if indel_sequence_record in UMI_family_size_counter['X'].keys():
                        UMI_family_size_counter['X'][indel_sequence_record].append(UMI_family_size[0])
                    else:
                        UMI_family_size_counter['X'][indel_sequence_record]=[UMI_family_size[0]]
                    UMI_family_size = UMI_family_size[1:] #only 1 UMI family size for the whole insertion so don't need to add more than 1

                    #add the position in read of the insertion to the position counter
                    if indel_sequence_record in position_counter['X'].keys():
                        position_counter['X'][indel_sequence_record].append(read_positions[0])
                    else:
                        position_counter['X'][indel_sequence_record]=[read_positions[0]]
                    read_positions = read_positions[1:]

                    #add the mapping quality score of the insertion to the mapping quality score counter
                    if indel_sequence_record in mapq_counter['X'].keys():
                        mapq_counter['X'][indel_sequence_record].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                    else:
                        mapq_counter['X'][indel_sequence_record]=[ord(mapq[0])]
                    mapq = mapq[1:]


                if 'N' in indel_sequence_record: #if the indel on that read contains an N, ignore that read's position
                    N_counter+=1
                    UMI_family_size = UMI_family_size[1:]
                    read_positions = read_positions[1:]
                    mapq = mapq[1:]


        elif bases[0] in ['<', '>']: #reference skipped due to CIGAR 'N' ('alignment gap')
            ##### SIMPLE REFERENCE SKIP ########
#             print("bases[0] in ['<', '?']")
            if (len(bases) == 1) or (bases[1] not in ['+', '-']):
                base_counter['<>']+=1
                if bases[0] == '>':
                    strand_counter['<>']['F']+=1
                if bases[0] == '<':
                    strand_counter['<>']['R']+=1
                bases = bases[1:]

                #add the UMI family size of the deletion to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                UMI_family_size_counter['<>'].append(UMI_family_size[0])
                UMI_family_size = UMI_family_size[1:]

                #add the position in read of the deletion to the position counter
                position_counter['<>'].append(read_positions[0])
                read_positions = read_positions[1:]

                #add the mapping quality score of the deletion to the mapping quality score counter
                mapq_counter['<>'].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                mapq = mapq[1:]

            ##### INSERTION OR DELETION FOLLOWING A REFERENCE SKIP ######## - RECORD AS COMPLEX
            elif bases[1] in ['+', '-']: #if the reference deletion is followed by insertion or deletion marker
                if bases[1] == '+':
                    indel_marker = '+'
                if bases[1] == '-':
                    indel_marker = '-'

                #add the indel to the indel counter and its strand to the strand counter
                if (bases[3].isdigit() == False): #i.e. if the length of the indel is <10 (i.e. only bases[2] is a number)
                    indel_length = int(bases[2])
                    indel_sequence = bases[3:(3+indel_length)] #the sequence that is inserted or deleted will be from 3rd position to length of insertion
                    indel_sequence_record = '<>'+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(3+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[3].isdigit() == True) and (bases[4].isdigit() == False): #i.e. if the length of the indel is 10-100:
                    indel_length = int(str(bases[2])+str(bases[3]))
                    indel_sequence = bases[4:(4+indel_length)] #the sequence that is inserted or deleted will be from 4th position to length of insertion
                    indel_sequence_record = '<>'+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(4+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[4].isdigit() == True):
                    indel_length = int(str(bases[2])+str(bases[3])+str(bases[4]))
                    indel_sequence = bases[5:(5+indel_length)] #the sequence that is inserted or deleted will be from 5th position to length of insertion
                    indel_sequence_record = '<>'+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(5+indel_length):] #then remove that indel from the base string (ready to look at the next base)

                if 'N' not in indel_sequence_record: #i.e. if the indel does not contain an N
                    if indel_sequence_record in base_counter['X'].keys():
                        base_counter['X'][indel_sequence_record]+=1
                    else:
                        base_counter['X'][indel_sequence_record]=1

                    if any(i for i in indel_sequence if i.isupper()): #if the indel is in uppercase it is on forward strand
                        if indel_sequence_record in strand_counter['X'].keys():
                            strand_counter['X'][indel_sequence_record]['F']+=1
                        else:
                            strand_counter['X'][indel_sequence_record]={'F': 1, 'R':0}
                    if any(i for i in indel_sequence if i.islower()): #if the indel is in lowercase it is on forward strand
                        if indel_sequence_record in strand_counter['X'].keys():
                            strand_counter['X'][indel_sequence_record]['R']+=1
                        else:
                            strand_counter['X'][indel_sequence_record]={'F': 0, 'R':1}

                    #add the UMI family size of the insertion to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                    if indel_sequence_record in UMI_family_size_counter['X'].keys():
                        UMI_family_size_counter['X'][indel_sequence_record].append(UMI_family_size[0])
                    else:
                        UMI_family_size_counter['X'][indel_sequence_record]=[UMI_family_size[0]]
                    UMI_family_size = UMI_family_size[1:] #only 1 UMI family size for the whole insertion so don't need to add more than 1

                    #add the position in read of the insertion to the position counter
                    if indel_sequence_record in position_counter['X'].keys():
                        position_counter['X'][indel_sequence_record].append(read_positions[0])
                    else:
                        position_counter['X'][indel_sequence_record]=[read_positions[0]]
                    read_positions = read_positions[1:]

                    #add the mapping quality score of the insertion to the mapping quality score counter
                    if indel_sequence_record in mapq_counter['X'].keys():
                        mapq_counter['X'][indel_sequence_record].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                    else:
                        mapq_counter['X'][indel_sequence_record]=[ord(mapq[0])]
                    mapq = mapq[1:]

                if 'N' in indel_sequence_record: #if the indel on that read contains an N, ignore that read's position
                    N_counter+=1
                    UMI_family_size = UMI_family_size[1:]
                    read_positions = read_positions[1:]
                    mapq = mapq[1:]

        elif bases[0].upper() in ['A', 'C', 'T', 'G']:
            ##### SIMPLE SNV/SNP ########
            if (len(bases) == 1) or (bases[1] not in ['-', '+']): #i.e. if alternate base and not followed by insertion of deletion
                upper_base = bases[0].upper()
                base_counter[upper_base]+=1
                if bases[0].isupper(): #if the base is in uppercase it is on forward strand
                    strand_counter[upper_base]['F']+=1
                if bases[0].islower(): #if the base is in lowercase it is on forward strand
                    strand_counter[upper_base]['R']+=1
                UMI_family_size_counter[upper_base].append(UMI_family_size[0])
                UMI_family_size = UMI_family_size[1:]
                position_counter[upper_base].append(read_positions[0])
                read_positions = read_positions[1:]
                mapq_counter[upper_base].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                mapq = mapq[1:]
                bases = bases[1:]

            ##### INSERTION OR DELETION FOLLOWING A SNV/SNP ######## - RECORD AS COMPLEX
            elif bases[1] in ['+', '-']: #if the reference deletion is followed by insertion or deletion marker
                if bases[1] == '+':
                    indel_marker = '+'
                if bases[1] == '-':
                    indel_marker = '-'

                #add the insertion to the insertion counter and its strand to the strand counter
                if (bases[3].isdigit() == False): #i.e. if the length of the indel is <10 (i.e. only bases[2] is a number)
                    indel_length = int(bases[2])
                    indel_sequence = bases[3:(3+indel_length)] #the sequence that is inserted or deleted will be from 3rd position to length of insertion
                    indel_sequence_record = bases[0].upper()+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(3+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[3].isdigit() == True) and (bases[4].isdigit() == False): #i.e. if the length of the indel is 10-100:
                    indel_length = int(str(bases[2])+str(bases[3]))
                    indel_sequence = bases[4:(4+indel_length)] #the sequence that is inserted or deleted will be from 4th position to length of insertion
                    indel_sequence_record = bases[0].upper()+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(4+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[4].isdigit() == True):
                    indel_length = int(str(bases[2])+str(bases[3])+str(bases[4]))
                    indel_sequence = bases[5:(5+indel_length)] #the sequence that is inserted or deleted will be from 5th position to length of insertion
                    indel_sequence_record = bases[0].upper()+bases[1]+indel_sequence.upper() #i.e. record as e.g. 'T+GTCA'
                    bases = bases[(5+indel_length):] #then remove that indel from the base string (ready to look at the next base)

                if 'N' not in indel_sequence_record: #i.e. if the indel does not contain an N
                    if indel_sequence_record in base_counter['X'].keys():
                        base_counter['X'][indel_sequence_record]+=1
                    else:
                        base_counter['X'][indel_sequence_record]=1

                    if any(i for i in indel_sequence if i.isupper()): #if the indel is in uppercase it is on forward strand
                        if indel_sequence_record in strand_counter['X'].keys():
                            strand_counter['X'][indel_sequence_record]['F']+=1
                        else:
                            strand_counter['X'][indel_sequence_record]={'F': 1, 'R':0}
                    if any(i for i in indel_sequence if i.islower()): #if the indel is in lowercase it is on forward strand
                        if indel_sequence_record in strand_counter['X'].keys():
                            strand_counter['X'][indel_sequence_record]['R']+=1
                        else:
                            strand_counter['X'][indel_sequence_record]={'F': 0, 'R':1}

                    #add the UMI family size of the insertion to the UMI_family size counter and then remove that UMI family size from the UMI family size string
                    if indel_sequence_record in UMI_family_size_counter['X'].keys():
                        UMI_family_size_counter['X'][indel_sequence_record].append(UMI_family_size[0])
                    else:
                        UMI_family_size_counter['X'][indel_sequence_record]=[UMI_family_size[0]]
                    UMI_family_size = UMI_family_size[1:] #only 1 UMI family size for the whole insertion so don't need to add more than 1

                    #add the position in read of the insertion to the position counter
                    if indel_sequence_record in position_counter['X'].keys():
                        position_counter['X'][indel_sequence_record].append(read_positions[0])
                    else:
                        position_counter['X'][indel_sequence_record]=[read_positions[0]]
                    read_positions = read_positions[1:]

                    #add the mapping quality score of the insertion to the mapping quality score counter
                    if indel_sequence_record in mapq_counter['X'].keys():
                        mapq_counter['X'][indel_sequence_record].append(ord(mapq[0])) #ord converts the ASCII to (decimal) number
                    else:
                        mapq_counter['X'][indel_sequence_record]=[ord(mapq[0])]
                    mapq = mapq[1:]

                if 'N' in indel_sequence_record: #if the indel on that read contains an N, ignore that read's position
                    N_counter+=1
                    UMI_family_size = UMI_family_size[1:]
                    read_positions = read_positions[1:]
                    mapq = mapq[1:]

        elif bases[0].upper() in ['N']:
            N_counter+=1
            ##### SIMPLE SNV/SNP ########
            if (len(bases) == 1) or (bases[1] not in ['-', '+']): #i.e. if alternate base and not followed by insertion of deletion
                UMI_family_size = UMI_family_size[1:]
                read_positions = read_positions[1:]
                mapq = mapq[1:]
                bases = bases[1:]

            ##### INSERTION OR DELETION FOLLOWING A SNV/SNP ######## - RECORD AS COMPLEX
            elif bases[1] in ['+', '-']: #if the reference deletion is followed by insertion or deletion marker
                if (bases[3].isdigit() == False): #i.e. if the length of the indel is <10 (i.e. only bases[2] is a number)
                    indel_length = int(bases[2])
                    bases = bases[(3+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[3].isdigit() == True) and (bases[4].isdigit() == False): #i.e. if the length of the indel is 10-100:
                    indel_length = int(str(bases[2])+str(bases[3]))
                    bases = bases[(4+indel_length):] #then remove that indel from the base string (ready to look at the next base)
                elif (bases[4].isdigit() == True):
                    indel_length = int(str(bases[2])+str(bases[3])+str(bases[4]))
                    bases = bases[(5+indel_length):] #then remove that indel from the base string (ready to look at the next base)

                UMI_family_size = UMI_family_size[1:] #only 1 UMI family size for the whole insertion so don't need to add more than 1
                read_positions = read_positions[1:]
                mapq = mapq[1:]

    UMI_counter = {}
    for k, v in UMI_family_size_counter.items():
        if k in ['A', 'G', 'C', 'T', '*', '<>']:
            if len(v) > 0:
                UMI_counter[k]=(np.mean(v), np.std(v))
            else:
                UMI_counter[k]=(0, 0)
        if k in ['+', '-', 'X']:
            if len(v)>0:
                new_v = {}
                for indel, UMI_size in v.items():
                    new_v[indel]=(np.mean(UMI_size), np.std(UMI_size))
                UMI_counter[k]=new_v
            else:
                UMI_counter[k]=(0, 0)

    pos_in_read_counter = {}
    for k, v in position_counter.items():
        if k in ['A', 'G', 'C', 'T', '*', '<>']:
            if len(v) > 0:
                pos_in_read_counter[k]=(np.mean(v), np.std(v))
            else:
                pos_in_read_counter[k]=(0, 0)
        if k in ['+', '-', 'X']:
            if len(v)>0:
                new_v = {}
                for indel, pos_read in v.items():
                    new_v[indel]=(np.mean(pos_read), np.std(pos_read))
                pos_in_read_counter[k]=new_v
            else:
                pos_in_read_counter[k]=(0, 0)

    mapping_qualities_counter = {}
    for k, v in mapq_counter.items():
        if k in ['A', 'G', 'C', 'T', '*', '<>']:
            if len(v) > 0:
                mapping_qualities_counter[k]=(np.mean(v), np.std(v))
            else:
                mapping_qualities_counter[k]=(0, 0)
        if k in ['+', '-', 'X']:
            if len(v)>0:
                new_v = {}
                for indel, mapq in v.items():
                    new_v[indel]=(np.mean(mapq), np.std(mapq))
                mapping_qualities_counter[k]=new_v
            else:
                mapping_qualities_counter[k]=(0, 0)

    recal_depth = depth-N_counter

#     if depth!= recal_depth:
#         print('depth = ', depth)
#         print('recal_depth = ', recal_depth)

    return base_counter, strand_counter, UMI_counter, pos_in_read_counter, mapping_qualities_counter, recal_depth

def VCF_calling(base_counter, strand_counter, UMI_counter, pos_in_read_counter, mapping_qualities_counter, mapq_filter, read_pos_filter, bias_filter, reference, chromosome, position, dbSNP_dict, depth, sample_name, VCF, TEXT):
    alt_positions = 0

    #SNV calling
    bases = ['T', 'C', 'G', 'A'] #ignore Ns
    for i in bases: #iterate through each of the alternate bases so can create a VCF line for them
        #i.e. if one of the bases is called as an alternate allele:
        if i != reference:
            if base_counter[i] != 0:
                alt_positions+=1 #record that there is an ALT base called at this position

                RSID, MAF = dbSNP_information(chromosome, position, reference, i, dbSNP_dict)

                #FIRST FEW COLUMNS
                ID = RSID #INSERT RS NUMBER IF SNP
                ref = reference
                alt = i

                #INFO FIELDS
                variant_type = 'SNV'
                end = position #for SNVs that are only 1 base change
                variant_depth = base_counter[i]
                allele_frequency = variant_depth/depth
                UMI_ref = UMI_counter[ref][0] #mean UMI family size for the ref call
                UMI_alt = UMI_counter[alt][0] #mean UMI family size for the alt call
                ID = RSID
                MAF_1000_genomes = MAF

        #     strand_bias = str(strand_counter[i]['F'])+':'+str(strand_counter[i]['R']) ?WHAT IS THIS RATIO?
                ref_strand_bias = str(strand_counter[ref]['F'])+':'+str(strand_counter[ref]['R'])
                var_strand_bias = str(strand_counter[alt]['F'])+':'+str(strand_counter[alt]['R'])
                mean_pos_read = pos_in_read_counter[alt][0]
                std_pos_read = pos_in_read_counter[alt][1]
                strand_fisher, FS = fisher_exact(strand_counter[ref]['F'], strand_counter[ref]['R'],
                                                       strand_counter[alt]['F'], strand_counter[alt]['R'])
                FS = round(FS, 2)
                mean_mapq = mapping_qualities_counter[alt][0]

                #FILTER FIELD
                filt = 'pass' #set as default, but if any of the filter criteria are met it will be changed:
                filt_list = []
                if mean_mapq < mapq_filter:
                    filt_list.append('Q'+str(mapq_filter))
                if mean_pos_read < read_pos_filter:
                    filt_list.append('p'+str(read_pos_filter))
                if strand_fisher > bias_filter:
                    filt_list.append('Bias')
                if len(filt_list)>0:
                    filt = (str(filt_list)).replace('[', '').replace(']', '').replace(',',';').replace("'",'').replace(" ",'')

                #FORMAT FIELDS
                if allele_frequency <0.8:
                    GT = '1/0'
                else:
                    GT = '1/1'
                DP = depth
                VD = variant_depth
                AD = str(base_counter[ref])+','+str(variant_depth) #Allelic depths for the ref and alt alleles
                AF = allele_frequency
                RD = str(strand_counter[ref]['F'])+','+str(strand_counter[ref]['R'])
                ALD = str(strand_counter[alt]['F'])+','+str(strand_counter[alt]['R'])

                BASIC_INFO = basics(chromosome, position, ID, ref, alt, filt)
                INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
                FORMAT = format_section(GT, DP, VD, AD, AF, RD, ALD)

                VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

                #Write the text file (which only contains positions with variants)
                end_position = position+(len(ref)-1)
                TEXT_FILE_INFO = text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref, alt, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
                TEXT.write(TEXT_FILE_INFO)


    #INSERTION VARIANT calling (+)
    if len(base_counter['+']) > 0: #i.e. there is an insertion present
        for k, v in base_counter['+'].items(): #in case there is more than 1 insertion:
            ref = reference
            alt = k
            number_bases = v
            alt_positions+=1 #record that there is an ALT base called at this position

            RSID, MAF = dbSNP_information(chromosome, position, reference, ref+alt, dbSNP_dict)

            #INFO FIELDS
            info_name = sample_name
            variant_type = 'Insertion'
            total_depth = depth #combined depth
            end = position #this is what VarDictJava does
            variant_depth = number_bases
            allele_frequency = variant_depth/total_depth
            UMI_ref = UMI_counter[ref][0] #mean UMI family size for the ref call
            UMI_alt = UMI_counter['+'][alt][0] #mean UMI family size for the alt call
            ID = RSID
            MAF_1000_genomes = MAF

        #     strand_bias = str(strand_counter[i]['F'])+':'+str(strand_counter[i]['R']) ?WHAT IS THIS RATIO?
            ref_strand_bias = str(strand_counter[ref]['F'])+':'+str(strand_counter[ref]['R'])
            var_strand_bias = str(strand_counter['+'][alt]['F'])+':'+str(strand_counter['+'][alt]['R'])
            mean_pos_read = pos_in_read_counter['+'][alt][0]
            std_pos_read = pos_in_read_counter['+'][alt][1]
            strand_fisher, FS = fisher_exact(strand_counter[ref]['F'], strand_counter[ref]['R'],
                                                   strand_counter['+'][alt]['F'], strand_counter['+'][alt]['R'])
            FS = round(FS, 2)
            mean_mapq = mapping_qualities_counter['+'][alt][0]

            #FILTER FIELD
            filt = 'pass' #set as default, but if any of the filter criteria are met it will be changed:
            filt_list = []
            if mean_mapq < mapq_filter:
                filt_list.append('Q'+str(mapq_filter))
            if mean_pos_read < read_pos_filter:
                filt_list.append('p'+str(read_pos_filter))
            if strand_fisher > bias_filter:
                filt_list.append('Bias')
        #     if variant adjacent to insertion variant....:
        #         filt_list.append('InIns')
            if len(filt_list)>0:
                filt = (str(filt_list)).replace('[', '').replace(']', '').replace(',',';').replace("'",'').replace(" ",'')

            #FORMAT FIELDS
            if 0.5 <= allele_frequency <0.8:
                GT = '1/0'
            if allele_frequency >= 0.8:
                GT = '1/1'
            if 0 < allele_frequency < 0.5:
                GT = '0/1'
            DP = depth
            VD = variant_depth
            AD = str(base_counter[ref])+','+str(variant_depth) #Allelic depths for the ref and alt alleles
            AF = allele_frequency
            RD = str(strand_counter[ref]['F'])+','+str(strand_counter[ref]['R'])
            ALD = str(strand_counter['+'][alt]['F'])+','+str(strand_counter['+'][alt]['R'])

            #Write the VCF file (which contains all positions)
            BASIC_INFO = basics(chromosome, position, ID, ref, ref+alt, filt)
            INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
    #             INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, ref_strand_bias, strand_fisher, FS, mean_mapq, five_prime, three_prime)
            FORMAT = format_section(GT, DP, VD, AD, AF, RD, ALD)

            VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

            #Write the text file (which only contains positions with variants)
            end_position = position+(len(ref)-1)
            TEXT_FILE_INFO = text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref, ref+alt, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
            TEXT.write(TEXT_FILE_INFO)


    #DELETION VARIANT calling ('-')
    if len(base_counter['-']) > 0: #i.e. there is an insertion present
        for k, v in base_counter['-'].items(): #in case there is more than 1 insertion:
            ref = reference
            alt = k
            number_bases = v
            alt_positions+=1 #record that there is an ALT base called at this position

            RSID, MAF = dbSNP_information(chromosome, position, reference+alt, ref, dbSNP_dict)

            #INFO FIELDS
            info_name = sample_name
            variant_type = 'Deletion'
            total_depth = depth #combined depth
            end = str(abs(int(position)+len(alt)))
            variant_depth = v
            allele_frequency = variant_depth/total_depth
            UMI_ref = UMI_counter[ref][0] #mean UMI family size for the ref call
            UMI_alt = UMI_counter['-'][alt][0] #mean UMI family size for the alt call
            ID = RSID
            MAF_1000_genomes = MAF

        #     strand_bias = str(strand_counter[i]['F'])+':'+str(strand_counter[i]['R']) ?WHAT IS THIS RATIO?
            ref_strand_bias = str(strand_counter[ref]['F'])+':'+str(strand_counter[ref]['R'])
            var_strand_bias = str(strand_counter['-'][alt]['F'])+':'+str(strand_counter['-'][alt]['R'])
            mean_pos_read = pos_in_read_counter['-'][alt][0]
            std_pos_read = pos_in_read_counter['-'][alt][1]
            strand_fisher, FS = fisher_exact(strand_counter[ref]['F'], strand_counter[ref]['R'],
                                                   strand_counter['-'][alt]['F'], strand_counter['-'][alt]['R'])
            FS = round(FS, 2)
            mean_mapq = mapping_qualities_counter['-'][alt][0]

            #FILTER FIELD
            filt = 'pass' #set as default, but if any of the filter criteria are met it will be changed:
            filt_list = []
            if mean_mapq < mapq_filter:
                filt_list.append('Q'+str(mapq_filter))
            if mean_pos_read < read_pos_filter:
                filt_list.append('p'+str(read_pos_filter))
            if strand_fisher > bias_filter:
                filt_list.append('Bias')
        #     if variant adjacent to insertion variant....:
        #         filt_list.append('InIns')
            if len(filt_list)>0:
                filt = (str(filt_list)).replace('[', '').replace(']', '').replace(',',';').replace("'",'').replace(" ",'')

            #FORMAT FIELDS
            if 0.5 <= allele_frequency <0.8:
                GT = '1/0'
            if allele_frequency >= 0.8:
                GT = '1/1'
            if 0 < allele_frequency < 0.5:
                GT = '0/1'
            DP = depth
            VD = variant_depth
            AD = str(base_counter[ref])+','+str(variant_depth) #Allelic depths for the ref and alt alleles
            AF = allele_frequency
            RD = str(strand_counter[ref]['F'])+','+str(strand_counter[ref]['R'])
            ALD = str(strand_counter['-'][alt]['F'])+','+str(strand_counter['-'][alt]['R'])

            BASIC_INFO = basics(chromosome, position, ID, ref+alt, ref, filt)
            INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
    #             INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, ref_strand_bias, strand_fisher, FS, mean_mapq, five_prime, three_prime)
            FORMAT = format_section(GT, DP, VD, AD, AF, RD, ALD)

            VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

            #Write the text file (which only contains positions with variants)
            end_position = position+(len(ref+alt)-1)
            TEXT_FILE_INFO = text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref+alt, ref, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
            TEXT.write(TEXT_FILE_INFO)

    #COMPLEX VARIANT calling ('X')
    if len(base_counter['X']) > 0: #i.e. there is a complex variant present
        ref = reference
        for k, v in base_counter['X'].items(): #in case there is more than 1 insertion:
            alt = k
            number_bases = v
            alt_positions+=1 #record that there is an ALT base called at this position

            if alt[0] in ['A', 'C', 'T', 'G', 'N']: #i.e. SNV followed by indel
                if alt[1] == '+':
                    variant_type = 'Complex insertion (SNV followed by insertion)'
                    alt_sequence = alt[0]+alt[2:] #alt is e.g. N+ATG (ref = N, followed by ATG insertion), so changed sequence is alt[2:]
                if alt[1] == '-':
                    variant_type = 'Complex deletion (SNV followed by deletion)'
                    alt_sequence = alt
                RSID, MAF = dbSNP_information(chromosome, position, reference, alt_sequence, dbSNP_dict)

                #INFO FIELDS
                info_name = sample_name
                total_depth = depth #combined depth
                end = str(abs(int(position)+len(alt)-2)) #subtract the first character and the '+' or '-'
                variant_depth = v
                allele_frequency = variant_depth/total_depth
                UMI_ref = UMI_counter[ref][0] #mean UMI family size for the ref call
                UMI_alt = UMI_counter['X'][alt][0] #mean UMI family size for the alt call
                ID = RSID
                MAF_1000_genomes = MAF

            #     strand_bias = str(strand_counter[i]['F'])+':'+str(strand_counter[i]['R']) ?WHAT IS THIS RATIO?
                ref_strand_bias = str(strand_counter[ref]['F'])+':'+str(strand_counter[ref]['R'])
                var_strand_bias = str(strand_counter['X'][alt]['F'])+':'+str(strand_counter['X'][alt]['R'])
                mean_pos_read = pos_in_read_counter['X'][alt][0]
                std_pos_read = pos_in_read_counter['X'][alt][1]
                strand_fisher, FS = fisher_exact(strand_counter[ref]['F'], strand_counter[ref]['R'],
                                                       strand_counter['X'][alt]['F'], strand_counter['X'][alt]['R'])
                FS = round(FS, 2)
                mean_mapq = mapping_qualities_counter['X'][alt][0]

                #FILTER FIELD
                filt = 'pass' #set as default, but if any of the filter criteria are met it will be changed:
                filt_list = []
                if mean_mapq < mapq_filter:
                    filt_list.append('Q'+str(mapq_filter))
                if mean_pos_read < read_pos_filter:
                    filt_list.append('p'+str(read_pos_filter))
                if strand_fisher > bias_filter:
                    filt_list.append('Bias')
                if len(filt_list)>0:
                    filt = (str(filt_list)).replace('[', '').replace(']', '').replace(',',';').replace("'",'').replace(" ",'')

                #FORMAT FIELDS
                if 0.5 <= allele_frequency <0.8:
                    GT = '1/0'
                if allele_frequency >= 0.8:
                    GT = '1/1'
                if 0 < allele_frequency < 0.5:
                    GT = '0/1'
                DP = depth
                VD = variant_depth
                AD = str(base_counter[ref])+','+str(variant_depth) #Allelic depths for the ref and alt alleles
                AF = allele_frequency
                RD = str(strand_counter[ref]['F'])+','+str(strand_counter[ref]['R'])
                ALD = str(strand_counter['X'][alt]['F'])+','+str(strand_counter['X'][alt]['R'])

                BASIC_INFO = basics(chromosome, position, ID, ref, alt_sequence, filt)
                INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
        #             INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, ref_strand_bias, strand_fisher, FS, mean_mapq, five_prime, three_prime)
                FORMAT = format_section(GT, DP, VD, AD, AF, RD, ALD)

                VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

                #Write the text file (which only contains positions with variants)
                end_position = position+(len(ref)-1)
                TEXT_FILE_INFO = text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref, alt_sequence, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
                TEXT.write(TEXT_FILE_INFO)

            else:
                RSID, MAF = dbSNP_information(chromosome, position, reference, alt, dbSNP_dict)

                #INFO FIELDS
                info_name = sample_name
                variant_type = 'Complex'
                total_depth = depth #combined depth
                end = str(abs(int(position)+len(alt)-2)) #subtract the first character and the '+' or '-'
                variant_depth = v
                allele_frequency = variant_depth/total_depth
                UMI_ref = UMI_counter[ref][0] #mean UMI family size for the ref call
                UMI_alt = UMI_counter['X'][alt][0] #mean UMI family size for the alt call
                ID = RSID
                MAF_1000_genomes = MAF

            #     strand_bias = str(strand_counter[i]['F'])+':'+str(strand_counter[i]['R']) ?WHAT IS THIS RATIO?
                ref_strand_bias = str(strand_counter[ref]['F'])+':'+str(strand_counter[ref]['R'])
                var_strand_bias = str(strand_counter['X'][alt]['F'])+':'+str(strand_counter['X'][alt]['R'])
                mean_pos_read = pos_in_read_counter['X'][alt][0]
                std_pos_read = pos_in_read_counter['X'][alt][1]
                strand_fisher, FS = fisher_exact(strand_counter[ref]['F'], strand_counter[ref]['R'],
                                                       strand_counter['X'][alt]['F'], strand_counter['X'][alt]['R'])
                FS = round(FS, 2)
                mean_mapq = mapping_qualities_counter['X'][alt][0]

                #FILTER FIELD
                filt = 'pass' #set as default, but if any of the filter criteria are met it will be changed:
                filt_list = []
                if mean_mapq < mapq_filter:
                    filt_list.append('Q'+str(mapq_filter))
                if mean_pos_read < read_pos_filter:
                    filt_list.append('p'+str(read_pos_filter))
                if strand_fisher > bias_filter:
                    filt_list.append('Bias')
                if len(filt_list)>0:
                    filt = (str(filt_list)).replace('[', '').replace(']', '').replace(',',';').replace("'",'').replace(" ",'')

                #FORMAT FIELDS
                if 0.5 <= allele_frequency <0.8:
                    GT = '1/0'
                if allele_frequency >= 0.8:
                    GT = '1/1'
                if 0 < allele_frequency < 0.5:
                    GT = '0/1'
                DP = depth
                VD = variant_depth
                AD = str(base_counter[ref])+','+str(variant_depth) #Allelic depths for the ref and alt alleles
                AF = allele_frequency
                RD = str(strand_counter[ref]['F'])+','+str(strand_counter[ref]['R'])
                ALD = str(strand_counter['X'][alt]['F'])+','+str(strand_counter['X'][alt]['R'])

                BASIC_INFO = basics(chromosome, position, ID, ref, alt, filt)
                INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
                FORMAT = format_section(GT, DP, VD, AD, AF, RD, ALD)

                VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

                #Write the text file (which only contains positions with variants)
                end_position = position+(len(ref)-1)
                TEXT_FILE_INFO = text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref, alt, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
                TEXT.write(TEXT_FILE_INFO)

    #IF NO SNV, SNP OR INDEL AT THIS POSITION:
    if alt_positions == 0:  #i.e. if there are no ALT alleles at this position

        ref = reference

        RSID, MAF_1000_genomes = homozygous_SNP_position(chromosome, position, reference, dbSNP_dict)

        #FIRST FEW COLUMNS
        alt = '-'
        ID = RSID

        #INFO FIELDS
        info_name = sample_name
        variant_type = 'REF'
        total_depth = depth #combined depth
        end = position #for SNVs that are only 1 base change
        variant_depth = 0
        allele_frequency = 0
        UMI_ref = UMI_counter[ref][0] #mean UMI family size for the ref call
        UMI_alt = 0 #mean UMI family size for the alt call

    #     strand_bias = str(strand_counter[i]['F'])+':'+str(strand_counter[i]['R']) ?WHAT IS THIS RATIO?
        ref_strand_bias = str(strand_counter[ref]['F'])+':'+str(strand_counter[ref]['R'])
        var_strand_bias = '0:0'
        mean_pos_read = pos_in_read_counter[ref][0]
        std_pos_read = pos_in_read_counter[ref][1]
        strand_fisher = 0
        FS = 0
        mean_mapq = mapping_qualities_counter[ref][0]

        #FILTER FIELD
        filt = 'pass' #set as default, but if any of the filter criteria are met it will be changed:
        filt_list = []
        if mean_mapq < mapq_filter:
            filt_list.append('Q'+str(mapq_filter))
        if mean_pos_read < read_pos_filter:
            filt_list.append('p'+str(read_pos_filter))
        if strand_fisher > bias_filter:
            filt_list.append('Bias')
        if len(filt_list)>0:
            filt = (str(filt_list)).replace('[', '').replace(']', '').replace(',',';').replace("'",'').replace(" ",'')

        #FORMAT FIELDS
        GT = '0/0'
        DP = depth
        VD = 0
        AD = str(base_counter[ref])+','+str(0) #Allelic depths for the ref and alt alleles
        AF = 0
        RD = str(strand_counter[ref]['F'])+','+str(strand_counter[ref]['R'])
        ALD = str(0)+','+str(0)

        BASIC_INFO = basics(chromosome, position, ID, ref, alt, filt)
        INFO = info_section(sample_name, variant_type, depth, end, variant_depth, allele_frequency, UMI_ref, UMI_alt, ID, MAF_1000_genomes, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
        FORMAT = format_section(GT, DP, VD, AD, AF, RD, ALD)

        VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

        if ID != '-':
            #Write the text file (which only contains positions with variants and 1000 genome variants)
            end_position = position+(len(ref)-1)
            TEXT_FILE_INFO = text_output(chromosome, position, end_position, ID, MAF_1000_genomes, ref, alt, depth, variant_depth, allele_frequency, variant_type, UMI_ref, UMI_alt, ref_strand_bias, var_strand_bias, mean_pos_read, std_pos_read, strand_fisher, FS, mean_mapq)
            TEXT.write(TEXT_FILE_INFO)

    return

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input mpileup txt file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--path-to-reference-genome", type=str, dest='reference_genome', help="path to reference genome", required=True)
    parser.add_argument("--ECS_type", type=str, dest='ECS_type', help="raw, SSCS or DCS", required=True)
    parser.add_argument("--dbSNP_directory", action="store", dest="dbSNP_directory", help="directory where dbSNP interval files are stored", required=True)
    parser.add_argument('--mapq_filter', type=float, default=10, dest='mapq_filter', help="mapping quality filter (no filtering, just highlights in filter column) [default = 10]")
    parser.add_argument('--read_pos_filter', type=int, default=8, dest='read_pos_filter', help="mean position in read filter (no filtering, just highlights in filter column) [default = 8]")
    parser.add_argument('--bias_filter', type=int, default=100, dest='bias_filter', help="Strand Bias Fisher p-value (Phred) filter (no filtering, just highlights in filter column) [default = 100]")
    parser.add_argument('--min_family_size', type=int, default=1, dest='min_family_size', help="min UMI family size that was used for SSCS calling")
    parser.add_argument("--variants_txt", action="store", dest="variants_txt", help="output file for text file of variants", required=True)
    parser.add_argument("--variants_vcf", action="store", dest="variants_vcf", help="output file for VCF file of variants", required=True)
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    mpileup = o.infile
    sample_name = o.sample_name
    reference_genome = o.reference_genome
    dbSNP_directory = o.dbSNP_directory
    mapq_filter = o.mapq_filter
    read_pos_filter = o.read_pos_filter
    bias_filter = o.bias_filter
    out_directory = o.out_directory
    min_family_size = o.min_family_size
    variants_txt = o.variants_txt
    variants_vcf = o.variants_vcf
    ECS_type = o.ECS_type

    dbSNP_dict = dbSNP_dictionary(dbSNP_directory, out_directory, sample_name)

    with open(mpileup, 'r') as pileupfile:
        read_reader = csv.reader(pileupfile, delimiter = '\t')  #csv.reader returns a reader object which will iterate over lines in the csvfile
        row_count = 0

        ####### CREATE TEXT FILE #######
        TEXT = open(variants_txt, 'w')
        TEXT.write('##fileDate='+str(VCF_date_today)+'\n')
        TEXT.write('##source=Watson_code_VCF_v'+str(version)+'\n')
        TEXT.write('##reference='+str(reference_genome)+'\n')
        TEXT.write('##sample='+sample_name+'\n')
        TEXT.write('##mean_position_in_read, std_position_in_read, mean_mapq are all calculated from just the variant reads, except where only REF is present'+'\n')
        TEXT.write('chromosome\t'+'position\t'+'end_position\t'+'REF\t'+'ALT\t'+'RSID\t'+'MAF_1000_genomes\t'+'total_depth\t'+'variant_depth\t'+\
                  'VAF\t'+'variant_type\t'+'REF_mean_UMI_family_size\t'+'variant_mean_UMI_family_size\t'+'mean_position_in_read\t'+\
                   'std_position_in_read\t'+'mean_mapq\t'+\
                   'REF_depth_by_strand\t'+'variant_depth_by_strand\t'+'strand_bias_fisher_p_value\t'+'strand_bias_fisher_p-value_phred\n')

        ####### CREATE VCF FILE #########
        VCF = open(variants_vcf, 'w')
        VCF.write('##fileformat=VCFv4.2'+'\n')
        VCF.write('##fileDate='+str(VCF_date_today)+'\n')
        VCF.write('##source=Watson_code_VCF_v'+str(version)+'\n')
        VCF.write('##reference='+str(reference_genome)+'\n')

        #INFO FIELDS
        VCF.write('##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">\n')
        VCF.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">\n')
        VCF.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        VCF.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End Position of variant">\n')
        VCF.write('##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">\n')
        VCF.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n')
        VCF.write('##INFO=<ID=UMI_ref,Number=1,Type=Float,Description="Mean of the sum of the SSCS UMI family sizes for the pair that made up the duplex for ref (e.g. 4 in one SSCS, 6 in the other SSCS = total 10)">\n')
        VCF.write('##INFO=<ID=UMI_alt,Number=1,Type=Float,Description="MMean of the sum of the SSCS UMI family sizes for the pair that made up the duplex for alt (e.g. 4 in one SSCS, 6 in the other SSCS = total 10)">\n')
        VCF.write('##INFO=<ID=RSID,Number=1,Type=String,Description="RS ID number">\n')
        VCF.write('##INFO=<ID=MAF_1000_genomes,Number=1,Type=String,Description="MAF in 1000 genomes">\n')
        VCF.write('##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Reference depth by strand">\n')
        VCF.write('##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">\n')
        VCF.write('##INFO=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">\n')
        VCF.write('##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">\n') #Standard deviation of mean position in read
        VCF.write('##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">\n')
        VCF.write('##INFO=<ID=FS,Number=1,Type=Float,Description="Strand Bias Fisher p-value as Phred (-10log10(p-value))">\n')
        VCF.write('##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">\n')

        #FILTER FIELDS
        VCF.write('##FILTER=<ID=pass,Description="Passed all filters">\n')
        VCF.write('##FILTER=<ID=Q'+str(mapq_filter)+',Description="Mean Mapping Quality Below '+str(mapq_filter)+'">\n')
        VCF.write('##FILTER=<ID=p'+str(read_pos_filter)+',Description="Mean Position in Reads Less than '+str(read_pos_filter)+'">\n')
        VCF.write('##FILTER=<ID=Bias,Description="Strand Bias">\n')

        #FORMAT FIELDS
        VCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        VCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        VCF.write('##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">\n')
        VCF.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
        VCF.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n')
        VCF.write('##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">\n')
        VCF.write('##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">\n')

        VCF.write('#CHROM\t'+'POS\t'+'ID\t'+'REF\t'+'ALT\t'+'FILTER\t'+'INFO\t'+'FORMAT\t'+sample_name+'\n')

        start_time = time.time()
        reset_time = time.time()
        for row in read_reader:
            if row_count > 0:
                chromosome = 'chr'+row[0]
                position = int(row[1])
                reference = row[2].upper() #make sure reference is upper case
                depth = int(row[3])

                #Delete unnescesary base information
                bases = row[4]
                orig_bases = row[4]

                if depth !=0:
                    bases = bases.replace('$', '') #get rid of last position covered by read symbols
                    for i in bases:  #get rid of first position covered by read symbols ('^' and then it's assoc quality)
                        if i == '^':
                            first_pos = bases.find('^')
                            assoc_quality = bases[first_pos+1]
                            bases = bases.replace('^'+assoc_quality, '')

                    #Get quality scores and mapping qualities
                    UMI_family_size = row[5]
                    mapq = row[6]
                    UMI_family_sizes = []
                    for i in UMI_family_size:
                        UMI_family_sizes.append(ord(i)-33)

                    #Get positions in read
                    read_positions = []
                    position_in_read = row[7].split(',')

                    for i in position_in_read:
                        read_positions.append(int(i))

                    position_details = position_analysis(reference, depth, bases, UMI_family_sizes, read_positions, mapq)
                    base_counter = position_details[0]
                    strand_counter = position_details[1]
                    UMI_counter = position_details[2]
                    pos_in_read_counter = position_details[3]
                    mapping_qualities_counter = position_details[4]
                    recal_depth = position_details[5]

                    VCF_calling(base_counter, strand_counter, UMI_counter, pos_in_read_counter, mapping_qualities_counter,
                                   mapq_filter, read_pos_filter, bias_filter, reference, chromosome, position, dbSNP_dict, recal_depth, sample_name, VCF, TEXT)

                else: #if no coverage at that position
                    BASIC_INFO = basics(chromosome, position, '.', reference, '.', '.')
                    INFO = 'SAMPLE='+sample_name+';TYPE=.;DP=0;END=.;VD=0;AF=0;UMI_ref=0;UMI_alt=0;RSID=.;MAF_1000_genomes=.;REFBIAS=.;VARBIAS=.;PMEAN=.;PSTD=.;QUAL=.;QSTD=.;SBF=.;FS=.;MQ=.'#+';'+five_prime+';'+three_prime
                    FORMAT = format_section('.', 0, 0, '0,0', 0, '0,0', '0,0')
                    VCF.write(BASIC_INFO+'\t'+INFO +'\t'+FORMAT+'\n')

            row_count+=1
            if row_count % 10000 == 0:
                time_now = time.time()
                print('calling variants in '+chromosome+' in sample '+sample_name)
                print('10000 positions called in sample '+sample_name+' in %s seconds' % int(time_now - reset_time))
                print('total '+str(row_count)+' positions called in sample '+sample_name+' in '+str(round((time_now - start_time)/60, 2))+' minutes')
                reset_time = time.time()

    VCF.close()
    TEXT.close()

    return print('VCF and text file creation complete')

if __name__ == "__main__":
	main()
