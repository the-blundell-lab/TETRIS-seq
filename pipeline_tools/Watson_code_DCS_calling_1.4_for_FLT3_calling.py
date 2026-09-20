#!/usr/bin/env python

'''''
Watson code for Duplex consensus calling for FLT3 calling.
Version 1.4 (October 2021)

Input:
    1) SSCS BAM file (paired end)
    2) sample name

Outputs:
    BAM files:
        1) paired-end BAM containing DCSs
            - the 'quality scores' are actually the sum of the SSCS family sizes for the 2 SSCS reads forming that DCS.  This value is capped at 40 though.
            - the qname is UMIA-UMIB: SSCS family size for R1 forward SSCS: SSCS family size for R2 forward SSCS: SSCS family size for R1 reverse SSCS, SSCS family size for R2 reverse SSCS
    Metrics files and Figures:
        1) metrics txt file containing info on number of SSCS and number of DCS called and number SSCS filtered
        2) metrics txt file containing info on proportion of Ns in DCS reads
        3) plot showing proportion of Ns in DCS reads

Usage:
Watson_code_DCS_calling_1.4_for_FLT3_calling.py  --infile SSCS_BAM --sample-name sample_name --out-directory directory to save files in

'''''
version = '1.4'

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


def find_duplex_pairs(input_SSCS_bam, sample_name):
    in_bam_SSCS = pysam.Samfile(input_SSCS_bam, "rb")
    SSCS_bam = in_bam_SSCS.fetch(until_eof=True)
    n = 0
    duplex_dictionary = {}

    for line in SSCS_bam:
        qname = line.qname #query name
        UMI = qname.split(':')[0]
        UMI_1 = UMI.split('-')[0]
        UMI_2 = UMI.split('-')[1]
        UMI_partner = UMI_2+'-'+UMI_1

        pair1_family_size = int(qname.split(':')[1])
        pair2_family_size = int(qname.split(':')[2])
        UMI_family_number = int(qname.split(':')[3]) #arbitrary number

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
#        quality_reads_string = ''
#        for i in quality_read:
#            quality_reads_string+=chr(i+33)

        cD_tag = line.get_tag('cD') #max read depth for SSCS, e.g. 7
        cM_tag = line.get_tag('cM') #min read depth for SSCS, e.g. 3
        cE_tag = line.get_tag('cE') #consensus error rate, e.g. 0
        cd_tag = line.get_tag('cd') #consensus depth, e.g. array('B', [7, 7, 7, 7, 7, 7, 7, 7, 7,...] = number of bases supporting consensus
        ce_tag = line.get_tag('ce') #consensus errors, e.g. array('B', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0...] = number of bases not supporting consensus
        RG_tag = line.get_tag('RG') #read group (sample name), e.g. C92_022_s7
        QS_tag = line.get_tag('QS') #prob errors, e.g. array('H', [259, 259, 259, 259, 259, 259, 259, 259,....] = summed quality scores at each position
        if line.is_read1 is True:
            read_pair = 'read1'
            partner_read_pair = 'read2'
            UMI_family_size = pair1_family_size
        if line.is_read2 is True:
            read_pair = 'read2'
            partner_read_pair = 'read1'
            UMI_family_size = pair2_family_size
        if line.is_reverse is True:
            strand = 'reverse'
            partner_strand = 'reverse'
        if line.is_reverse is False:
            strand = 'forward'
            partner_strand = 'forward'

        key = (UMI, chromosome, coordinate, read_pair, strand, partner_chromosome, partner_coordinate, abs(template_length)) #there will only be one of these because this was the key used for SSCS calling
        partner_key = (UMI_partner, chromosome, coordinate, partner_read_pair, partner_strand, partner_chromosome, partner_coordinate, abs(template_length))

        read_info = (sequence, flag, cigar, qname, mapq, quality, cD_tag, cM_tag, cE_tag, cd_tag, ce_tag, RG_tag, QS_tag, UMI_family_size, read_pair, strand, abs(template_length), partner_chromosome, partner_coordinate, UMI)

        if key in duplex_dictionary.keys():
            duplex_dictionary[key].append(read_info)

        if key not in duplex_dictionary.keys():
            if partner_key in duplex_dictionary.keys():
                duplex_dictionary[partner_key].append(read_info)
            else:
                duplex_dictionary[key]=[read_info]

        n+=1
        if n % 100000 == 0:
            print(str(n)+' SSCS reads processed for sample '+sample_name)

    return duplex_dictionary

def duplex_consensus_read_caller(duplex_reads, read_length, quality_threshold): #function for calling the consensus sequence from the duplex UMI pairs
    #duplex_reads = [(sequence, base qualities (uncapped), UMI family size), (sequence, base qualities (uncapped), UMI family size)]

    consensus_read = '' #start the read as empty (and then will add to it...)
    consensus_depth = []
    total_consensus_reads = [] #the summed SSCS family sizes for each position that make up the consensus (capped at 40)
    read_depths = []
    prob_error = [] #the summed Phred scores for each position that make up the consensus (uncapped)

    for base_number in range(read_length): #goes through each base (with length e.g 146 bases)
        nucleotides = {'T': 0, 'C': 0, 'G': 0, 'A': 0, 'N': 0}
        base_qualities = {'T': 0, 'C': 0, 'G': 0, 'A': 0, 'N': 0}
        base_family_sizes = {'T': 0, 'C': 0, 'G': 0, 'A': 0, 'N': 0}
        total_number_nucleotides = 0 #reset
        read_depth=0 #keep track of the number of reads that pass the min quality filter at this position

        for read in range(len(duplex_reads)): #go through each read in the list (where j is each read) (there will be 2 reads)
            try:
                base_quality = duplex_reads[read][1][base_number] #this will have been capped at 40 - most will be 40 or above
                family_size = duplex_reads[read][2][base_number] #SSCS family size contributing to that SSCS consensus base
                if base_quality >= quality_threshold: #ignore that base if the quality is less than the set threshold
                    read_depth +=1
                    if duplex_reads[read][0][base_number] == 'T':
                        nucleotides['T'] +=1
                        base_qualities['T']+=base_quality #add the base qualtities for all reads at that position carrying that base
                        base_family_sizes['T']+=family_size #add the family size that created that consensus base
                    if duplex_reads[read][0][base_number] == 'C':
                        nucleotides['C'] +=1
                        base_qualities['C']+=base_quality
                        base_family_sizes['C']+=family_size #add the family size that created that consensus base
                    if duplex_reads[read][0][base_number] == 'G':
                        nucleotides['G'] +=1
                        base_qualities['G']+=base_quality
                        base_family_sizes['G']+=family_size #add the family size that created that consensus base
                    if duplex_reads[read][0][base_number] == 'A':
                        nucleotides['A'] +=1
                        base_qualities['A']+=base_quality
                        base_family_sizes['A']+=family_size #add the family size that created that consensus base
                    if duplex_reads[read][0][base_number] == 'N':
                        nucleotides['N'] +=1
                        base_qualities['N']+=base_quality
                        base_family_sizes['N']+=family_size #add the family size that created that consensus base
                    total_number_nucleotides +=1
            except:
                break

        try: #for each base, work out which nucleotides (if any) is >e.g. 90% of the nucleotides at that position
            for base_option in ['T', 'C', 'G', 'A', 'N']:
                if float(nucleotides[base_option]) == 2.0: #if only 1 read passed quality score, then s consensus is not called.  If both passed, but disagree, then this value will be <2
                    consensus_read += base_option #add that nucleotide to the consensus read

                    quality_score = base_qualities[base_option]
                    prob_error.append(quality_score)

                    total_family_size = base_family_sizes[base_option]

                    if total_family_size > 40:
                        total_consensus_reads.append(40) #will store in quality score position in bam and quality scores are unsigned char (1-255), but max seems to be 40 (was outputting No quality scores in BAM when too high (even if <==255))
                    if total_family_size <= 40:
                        total_consensus_reads.append(total_family_size)

                    consensus_depth.append(nucleotides[base_option]) #will usually be 2 unless 1 base was not excluded e.g. due to poor quality

                    break

                elif base_option == 'N': #if get to N and none of the nucleotides have been >e.g. 100% of the reads, then call an 'N'
                    consensus_read += 'N'
                    prob_error.append(0)
                    consensus_depth.append(0)
                    total_consensus_reads.append(0)

        except:
            consensus_read +='N'
            prob_error.append(0)
            consensus_depth.append(0)
            total_consensus_reads.append(0)

        read_depths.append(read_depth) #append the read depth (for reads that pass the quality filter) at this position

    max_read_depth = max(read_depths) #max read depth for the consensus
    min_read_depth = min(read_depths) #min read depth for the consensus (i.e. when poor quality bases were filtered out)

    total_SSCS_family_size = array('B', total_consensus_reads)

    return consensus_read, consensus_depth, max_read_depth, min_read_depth, prob_error, total_SSCS_family_size

def call_duplex_consensus_reads(input_SSCS_bam, DCS_bam, sample_name, out_directory, min_UMI_family_size, max_N, quality_threshold, duplex_dictionary, version):

    def reverse(quality_scores):
        return quality_scores[::-1]

    new_header = {}
    in_bam_SSCS = pysam.Samfile(input_SSCS_bam, "rb")
    heading = in_bam_SSCS.header
    new_header['HD'] = heading['HD']
    new_header['SQ'] = heading['SQ']
    new_header['RG'] = heading['RG']
    new_PG = [{'ID': 'Watson_code_DCS_calling_v'+str(version),
    'PN': 'Watson_code_SNV_panel_v'+str(version),
    'VN': version,
    'CL': 'Watson_code_DCS_calling_'+str(version)+'.py --infile '+input_SSCS_bam+' --sample-name '+sample_name+' --min-family-size '+str(min_UMI_family_size)+\
    ' --max_N '+str(max_N)+' --outbam '+ DCS_bam+'--min-base-quality '+str(quality_threshold)+' --out-directory '+out_directory}]
    new_header['PG'] = new_PG

    #CREATE NEW BAM AND FASTQ FILES
    out_bam = pysam.Samfile(DCS_bam, "wb", header=new_header) #open a new BAM file to write the SSCS tov
    # unpaired_out_bam = pysam.Samfile(unpaired_SSCS_bam, "wb", header=new_header) #open a new BAM file to write the SSCS to
    # fastq = gzip.open(DCS_fastq, "wt") #interleaved fastq file (r1 followed by r2 etc..)

    unpartnered_SSCS = 0
    DCS_created = 0
    partnered_DCS = 0
    filtered_SSCS = 0

    counter = 0
    duplex_counter = 0
    too_many_Ns_in_consensus = 0

    proportion_N_dict = {}
    consensus_dict = {}

    start_time = time.time()

    for k, v in duplex_dictionary.items():
        consensus_read = ''

        if len(v) == 1: #SSCS without a duplex partner
            UMI_family_size = v[0][13]
            if UMI_family_size>=min_UMI_family_size:
                unpartnered_SSCS+=1
            else:
                filtered_SSCS+=1

        if len(v) == 2: #SSCS with a duplex partner
            strand = k[4] #forward or reverse

            if v[0][14] == 'read1': #read a is whichever is read 1 and read b is whichever is read 2
                read_a = v[0]
                read_b = v[1]
            else:
                read_a = v[1]
                read_b = v[0]

            read_a_sequence = read_a[0]
            read_a_flag = read_a[1]
            read_a_CIGAR = read_a[2]
            read_a_qualities = read_a[5]
            read_a_cD = read_a[6] #max read depth for SSCS, e.g. 7
            read_a_cM = read_a[7] #min read depth for SSCS, e.g. 3
            read_a_cE = read_a[8] #consensus error rate, e.g. 0
            read_a_cd = read_a[9] #consensus depth, e.g. array('B', [7, 7, 7, 7, 7, 7, 7, 7, 7,...] = number of bases supporting consensus
            read_a_ce = read_a[10] #consensus errors, e.g. array('B', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0...] = number of bases not supporting consensus
            read_a_RG = read_a[11] #read group (sample name), e.g. C92_022_s7
            read_a_QS = read_a[12] #prob errors, e.g. array('H', [259, 259, 259, 259, 259, 259, 259, 259,....] = summed quality scores at each position
            read_a_UMI_family_size = read_a[13]
            read_a_pair = read_a[14]
            read_a_UMI = read_a[19] #read1 UMI

            read_b_sequence = read_b[0]
            read_b_flag = read_b[1]
            read_b_CIGAR = read_b[2]
            read_b_qualities = read_b[5]
            read_b_cD = read_b[6]
            read_b_cM = read_b[7]
            read_b_cE = read_b[8]
            read_b_cd = read_b[9]
            read_b_ce = read_b[10]
            read_b_RG = read_b[11]
            read_b_QS = read_b[12]
            read_b_UMI_family_size = read_b[13]
            read_b_pair = read_b[14]
            read_b_UMI = read_b[19] #read2 UMI

            read_a_info = (read_a_sequence, read_a_QS, read_a_cd)
            read_b_info = (read_b_sequence, read_b_QS, read_b_cd)

            duplex_reads = [read_a_info, read_b_info]
            read_length = len(read_a_sequence) #will be the same as read_b

            ###### CALL THE DUPLEX CONSENSUS FOR READS WHOSE PARTNER STRAND IS PRESENT ########
            if read_a_pair == 'read1' and read_b_pair == 'read2': #make sure they are duplex partners
                if read_a_UMI_family_size >=min_UMI_family_size:
                    if read_b_UMI_family_size >=min_UMI_family_size: #if both reads have minimum UMI family size
                        consensus_read, consensus_depth, max_depth, min_depth, uncapped_quality, total_SSCS_family_size = duplex_consensus_read_caller(duplex_reads, read_length, quality_threshold)
                        DCS_created+=1
                    else: #only read A has the minimum UMI family size
                        unpartnered_SSCS+=1
                        filtered_SSCS+=1

                if read_a_UMI_family_size <min_UMI_family_size: #if read A doesn't have minimum UMI family size
                    if read_b_UMI_family_size >=min_UMI_family_size: #if only read B meets minimum UMI family size
                        unpartnered_SSCS+=1
                        filtered_SSCS+=1
                    else: #if neither meet the minimum UMI family size
                        filtered_SSCS+=2

            elif read_a_pair == 'read2' and read_b_pair == 'read1': #make sure they are duplex partners
                if read_a_UMI_family_size >=min_UMI_family_size:
                    if read_b_UMI_family_size >=min_UMI_family_size: #if both reads have minimum UMI family size
                        consensus_read, consensus_depth, max_depth, min_depth, uncapped_quality, total_SSCS_family_size = duplex_consensus_read_caller(duplex_reads, read_length, quality_threshold)
                        DCS_created+=1
                    else: #only read A has the minimum UMI family size
                        unpartnered_SSCS+=1
                        filtered_SSCS+=1

                if read_a_UMI_family_size <min_UMI_family_size: #if read A doesn't have minimum UMI family size
                    if read_b_UMI_family_size >=min_UMI_family_size: #if only read B meets minimum UMI family size
                        unpartnered_SSCS+=1
                        filtered_SSCS+=1
                    else: #if neither meet the minimum UMI family size
                        filtered_SSCS+=2

            else: #both from the same read
                unpartnered_SSCS+=2

            ###### FILTER OUT DUPLEX CONSENSUS READS THAT HAVE TOO MANY N'S #######
            if len(consensus_read)>1: #if a duplex consensus was made
                N_count = consensus_read.count('N')
                consensus_read_length = len(consensus_read)
                proportion_N = N_count/consensus_read_length

                if proportion_N in proportion_N_dict.keys(): #keep a tally of what proportion of the reads are N
                    proportion_N_dict[proportion_N]+=1
                else:
                    proportion_N_dict[proportion_N]=1

                ####### CREATE THE READ1 AND READ2 ###########
                if proportion_N < max_N: #if proportion of Ns in the consensus is less than the set limit:
                    duplex_counter+=1

                    chromosome = k[1]
                    coordinate = k[2]
                    strand = k[4] #forward or reverse
                    if strand == 'forward':
                        read_pair = 'read1'
                        flag = 99
                        template_length = k[7]
                        UMI_name = read_a_UMI+':'+read_b_UMI #put the UMI of the top strand first then UMI of bottom strand next
                    if strand == 'reverse':
                        read_pair = 'read2'
                        flag = 147
                        template_length = -k[7] #template length is negative if on reverse strand
                        UMI_name = read_b_UMI+':'+read_a_UMI #put the UMI of the top strand first then UMI of bottom strand next
                    partner_chromosome = k[5]
                    partner_coordinate = k[6]
                    consensus = consensus_read
    #                 quality = quality_scores
                    quality = total_SSCS_family_size

                    quality_string = '' #for the output fastq files
                    for i in quality:
                        quality_string+=chr(i+33)

                    total_UMI_family_size = read_a_UMI_family_size + read_b_UMI_family_size
                    CIGAR = read_a_CIGAR

                    aD = read_a_cD #max depth for each SSCS, e.g. 5
                    bD = read_b_cD #max depth for each SSCS, e.g. 7
                    aM = read_a_cM #min depth for each SSCS, e.g. 5
                    bM = read_b_cM #min depth for each SSCS, e.g. 5
                    aE = read_a_cE #SSCS error rate
                    bE = read_b_cE #SSCS error rate
                    ad = read_a_cd #SSCS depth, e.g. array('B', [7, 7, 7, 7, 7, 7, 7, 7, 7,...] (number of bases supporting consensus)
                    bd = read_b_cd #SSCS depth, e.g. array('B', [7, 7, 7, 7, 7, 7, 7, 7, 7,...]
                    ae = read_a_ce #SSCS errors, e.g. array('B', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0...] (number of bases not supporting consensus)
                    be = read_b_ce #SSCS errors, e.g. array('B', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0...]
                    aq = read_a_QS #SSCS qualities (uncapped)
                    bq = read_b_QS #SSCS qualities (uncapped)
                    mq = uncapped_quality #DCS uncapped qualities

                    RG = read_a_RG

                    key = (UMI_name, chromosome, coordinate, read_pair, strand, flag, partner_chromosome, partner_coordinate, template_length)
                    value = (consensus, CIGAR, read_a_UMI_family_size, read_b_UMI_family_size, total_UMI_family_size, quality, aD, bD, aM, bM,
                             aE, bE, ad, bd, ae, be, aq, bq, mq, RG, quality_string)

                    consensus_dict[key]=value

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
                    new_read.set_tags([('aD', aD), ('bD', bD), ('aM', aM), ('bM', bM), ('aE', aE), ('bE', bE), ('ad', ad), ('bd', bd),
                                       ('ae', ae), ('be', be), ('aq', aq), ('bq', bq), ('mq', mq), ('RG', read_a_RG)])


                    #partner read details...
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

                    if template_length <0: #i.e. negative number
                        partner_template_length = abs(template_length) #i.e. the positive version
                    if template_length >0: #i.e. positive number
                        partner_template_length = -template_length
                    if template_length == 0:
                        partner_template_length = 0

                    partner_key = (UMI_name, partner_chromosome, partner_coordinate, partner_read_pair, partner_strand, partner_flag, chromosome,
                                  coordinate, partner_template_length)

                    #look to see if the partner key is already in the dictionary...(it won't be for the first read looked at, but then it will find the first read when it's partner is later added)
                    if partner_key in consensus_dict.keys(): #i.e. if the reads are paired (has the partner survived the N filter and min family size filter?)
                        #create the partner read entry for the BAM file...
                        partner = consensus_dict[partner_key] #retrieve the 'values' for the partner (i.e. (consensus, CIGAR, read_a_UMI family size))
                        partner_consensus = partner[0]
                        partner_CIGAR = partner[1]
                        partner_read_a_UMI_family_size = partner[2]
                        partner_read_b_UMI_family_size = partner[3]
                        partner_total_UMI_family_size = partner[4]
                        partner_quality = partner[5]
                        partner_aD = partner[6]
                        partner_bD = partner[7]
                        partner_aM = partner[8]
                        partner_bM = partner[9]
                        partner_aE = partner[10]
                        partner_bE = partner[11]
                        partner_ad = partner[12]
                        partner_bd = partner[13]
                        partner_ae = partner[14]
                        partner_be = partner[15]
                        partner_aq = partner[16]
                        partner_bq = partner[17]
                        partner_mq = partner[18]
                        partner_RG = partner[19]
                        partner_quality_string = partner[20]


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
                        partner_read.set_tags([('aD', partner_aD), ('bD', partner_bD), ('aM', partner_aM), ('bM', partner_bM), ('aE', partner_aE),
                                               ('bE', partner_bE), ('ad', partner_ad), ('bd', partner_bd), ('ae', partner_ae), ('be', partner_be),
                                               ('aq', partner_aq), ('bq', partner_bq), ('mq', partner_mq), ('RG', partner_RG)])


                        if read_pair == 'read1':
                            new_read.qname = UMI_name+':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+\
                            ':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+':'+str(duplex_counter)
                            partner_read.qname = UMI_name+':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+\
                            ':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+':'+str(duplex_counter)

                            #write the read to the BAM file
                            out_bam.write(new_read) #write the new DCS to new bam file (read 1 first)
                            #then write it's partner to the BAM file...
                            out_bam.write(partner_read) #write the new DCS to new bam file

                            # #add a /1 or /2 on the end for the fastq qnames
                            # new_read.qname = UMI_name+':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+\
                            # ':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+':'+str(duplex_counter)+'/1' #read 1 UMI family size followed by read 2 UMI family size
                            # partner_read.qname = UMI_name+':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+\
                            # ':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+':'+str(duplex_counter)+'/2' #same qname as its partner
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
                            new_read.qname = UMI_name+':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+\
                            ':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+':'+str(duplex_counter)
                            partner_read.qname = UMI_name+':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+\
                            ':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+':'+str(duplex_counter)

                            #write the partner to the BAM file first...
                            out_bam.write(partner_read) #write the new DCS to new bam file (read 1 first)
                            #then write the read to the BAM file...
                            out_bam.write(new_read) #write the new DCS to new bam file

                            # #add a /1 or /2 on the end for the fastq qnames
                            # new_read.qname = UMI_name+':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+\
                            # ':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+':'+str(duplex_counter)+'/2'
                            # partner_read.qname = UMI_name+':'+str(partner_read_a_UMI_family_size)+':'+str(partner_read_b_UMI_family_size)+\
                            # ':'+str(read_a_UMI_family_size)+':'+str(read_b_UMI_family_size)+':'+str(duplex_counter)+'/1'
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

                        del consensus_dict[partner_key] #delete the partner key from the dictionary now that it has been written to the BAM file/ fastq file
                        del consensus_dict[key]

                else: #if there are too many Ns in the consensus
                    too_many_Ns_in_consensus+=1 #add to the running counter

        counter+=1

        if counter % 10000 == 0:
            print('Calling DCSs: SSCS processed for sample '+sample_name+' = ', counter)
            print('time for last 10,000 SSCS to be processed = %s seconds' % int(time.time() - start_time))
            start_time = time.time() #reset the timer so it can calculate the time for the next 100 consensuses

    remaining_DCS = len(consensus_dict)


    ##### CREATE METRICS FILE ########
    DCS_metrics = open(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_DCS_metrics.txt', 'w')
    DCS_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    DCS_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    DCS_metrics.write('produced from code:\t' + 'Watson_code_DCS_calling: version '+str(version) + '\n')
    DCS_metrics.write('minimum family size required for SSCS consensus set to :\t'+ str(min_UMI_family_size)+ '\n')
    DCS_metrics.write('maximum proportion of Ns allowed in DCS:\t'+ str(max_N)+ '\n')
    DCS_metrics.write('minimum base quality set to:\t' + str(quality_threshold) + '\n')
    DCS_metrics.write('\n')
    DCS_metrics.write('total number of SSCS included in DCS consensus calling:\t'+ str(unpartnered_SSCS+filtered_SSCS+(DCS_created*2)+(too_many_Ns_in_consensus*2))+ '\n')
    DCS_metrics.write('total number of SSCS without duplex partner:\t' + str(unpartnered_SSCS) + '\n')
    DCS_metrics.write('total number of SSCS filtered out due to UMI family size:\t' + str(filtered_SSCS) + '\n')
    DCS_metrics.write('total number of DCS created (read 1 and read 2 separate):\t' + str(DCS_created+too_many_Ns_in_consensus) + '\n')
    DCS_metrics.write('total number of DCS filtered due to too many Ns:\t' + str(too_many_Ns_in_consensus) + '\n')
    DCS_metrics.write('total number of DCS after N filtering:\t' + str(DCS_created) + '\n')
    DCS_metrics.write('total number of unpartnered DCS:\t' + str(remaining_DCS)+ '\n')
    DCS_metrics.write('SSCS: DCS ratio:\t' + str((unpartnered_SSCS+filtered_SSCS+(DCS_created*2)+(too_many_Ns_in_consensus*2))/DCS_created) + '\n')
    DCS_metrics.write('\n')
    DCS_metrics.close()

    ####### CREATE PROPORTION N METRICS FILE ######
    DCS_UMI_N_metrics = open(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_DCS_UMI_N_metrics.txt', 'w')
    DCS_UMI_N_metrics.write('sample name :\t'+ str(sample_name)+ '\n')
    DCS_UMI_N_metrics.write('date of analysis :\t'+ str(date_today)+ '\n')
    DCS_UMI_N_metrics.write('produced from code:\t' + 'Watson_code_DCS_calling: version '+str(version) + '\n')
    DCS_UMI_N_metrics.write('minimum family size required for SSCS consensus set to :\t'+ str(min_UMI_family_size)+ '\n')
    DCS_UMI_N_metrics.write('maximum proportion of Ns allowed in DCS:\t'+ str(max_N)+ '\n')
    DCS_UMI_N_metrics.write('minimum base quality set to:\t' + str(quality_threshold) + '\n')
    DCS_UMI_N_metrics.write('\n')
    DCS_UMI_N_metrics.write('proportion of DCS that is N\t'+ 'count (read 1 and read 2 separate)\t' + '\n')

    for k in sorted(proportion_N_dict.keys()):
        proportion_N = k
        number_molecules = proportion_N_dict[k]
        DCS_UMI_N_metrics.write(str(proportion_N)+'\t' + str(number_molecules) + '\n')

    DCS_UMI_N_metrics.close()


    ##### PLOT PROPORTION OF DCS THAT ARE N #####
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
    ax1.set_title('Proportion of DCS read called as N \n'+sample_name, fontsize = titlefont, fontweight='bold')
    ax1.set_ylabel('cumulative fraction of DCS reads \n(with N proportion >= x)', fontsize = axislabelfont)
    ax1.set_xlabel('proportion of DCS read called as N', fontsize = axislabelfont)
    ax1.set_ylim(-0.01, 1.01)
    ax1.set_xlim(-0.01, 1.01)

    ax1.grid(which = 'both', lw = 0.8, linestyle = ':', color = grey3)

    plt.tight_layout()
    plt.savefig(out_directory+'/Metrics_and_images/'+sample_name+'_watson_code_DCS_proportion_N_distribution.pdf')

    ###### CLOSE THE FILES #######
    out_bam.close() #close the new DCS BAM file
    # fastq.close() #close the fastq file
    in_bam_SSCS.close() #close the open in BAM file

    print('total number of SSCS processed = ', unpartnered_SSCS+filtered_SSCS+(DCS_created*2)+(too_many_Ns_in_consensus*2))
    print('total number of SSCS without duplex partner = ', unpartnered_SSCS)
    print('total number of SSCS filtered out due to UMI family size = ', filtered_SSCS)
    print('total number of DCS created = ', DCS_created+too_many_Ns_in_consensus)
    print('total number of DCS filtered due to too many Ns = ', too_many_Ns_in_consensus)
    print('total number of DCS after N filtering = ', DCS_created)
    print('total number of unpartnered DCS = ', remaining_DCS)

    return print('DCS BAM created')

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="SSCS BAM file", required=True)
    parser.add_argument("--sample-name", type=str, dest='sample_name', help="name of sample to prefix file names with", required=True)
    parser.add_argument("--outbam", action="store", dest="outbam", help="output DCS BAM file", required=True)
    parser.add_argument("--min-family-size-SSCS", type=int, action="store", dest="min_UMI_family_size", help="min UMI family size used when calling SSCSS", required=True)
    parser.add_argument('--max_N', type=float, default=1, dest='max_N', help="Maximum fraction of Ns permitted in a consensus [default = 1.0]")
    parser.add_argument('--min-base-quality', type=int, default=10, dest='quality_threshold', help="Minimum per-base quality (bases with quality scores below this will not be included in consensus calling) [default = 10]")
    parser.add_argument("--out-directory", action="store", dest="out_directory", help="output directory where output files will be stored", required=True)
    o = parser.parse_args()

    input_SSCS_bam = o.infile
    DCS_bam = o.outbam
    sample_name = o.sample_name
    out_directory = o.out_directory
    min_UMI_family_size = o.min_UMI_family_size
    max_N = o.max_N
    quality_threshold = o.quality_threshold

    duplex_dictionary = find_duplex_pairs(input_SSCS_bam, sample_name)
    call_duplex_consensus_reads(input_SSCS_bam, DCS_bam, sample_name, out_directory, min_UMI_family_size, max_N, quality_threshold, duplex_dictionary, version)

    return print('DCS calling complete')

if __name__ == "__main__":
	main()
