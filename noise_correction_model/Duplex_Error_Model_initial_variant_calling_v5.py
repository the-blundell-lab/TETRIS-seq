#!/usr/bin/env python

'''''
Watson code for duplex position-specific error model, calling variants from across each sequencing lane
Version 5.0

Input:
    1) Sequencing library ID (assuming all the samples to run the error model are from the same sequencing library/ sequencing library sample sheet contains all the samples to be included in the error model)
    2) DCS or SSCS (whether running in DCS or SSCS variant calls)
    3) list of samples to exclude (if any)

    Also requires the following files:
    - CSV file containing sample info for all samples to be included in the error model: 'SNV_panel/'+library+'/'+library+'_sample_sheet.csv' (example shown in 'Data_files/Sample_indexes_example.csv'
    - VCF files (containing all positions in the panel) produced via Watson_code_SNV_panel_v1.7.sh, for all the samples to be included in the error model: 'UKCTOCS_sample_level_results_V2/'+sample_abbrev_name+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_all_positions.vcf'
    - Annotated variant call files produced via Watson_code_SNV_panel_v1.7.sh, for all the samples to be included in the error model: 'UKCTOCS_sample_level_results_V2/'+sample_abbrev+'/SNV/DCS/'+sample+'_SNV_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_annotated.txt'
    - Fasta file of reference genome (e.g. 'Homo_sapiens_assembly19.fasta')
    - CSV file containing every possible base change at each position in the custom panel (e.g. 'Data_files/TWIST_SNV_panel_all_possible_changes_panel_sites.csv')
    - CSV file containing every position in the custom panel that is also observed as mutated in COSMIC (+ COSMIC info) (e.g. 'Data_files/All_SNV_panel_positions_seen_in_COSMIC_haem_lymphoid.csv')

Outputs:
    1) Annotated variant call file for each sample, with information on whether called as errors or real variants ('UKCTOCS_sample_level_results_V2/'+sample_abbrev+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_watson_code_'+SSCS_or_DCS+'_beta_binomial_SNV_all_variant_calls.txt')
    2) File containing p-values for all positions in the panel called as errors ('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_p_values_for_variants_called_as_errors.csv')
    3) File containing p-values for all positions in the panel called as real ('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_p_values_for_variants_called_as_real.csv')
    4) MLE outcomes file for each position ('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_model_MLE_outcome_per_position.csv')

Usage:
Duplex_Error_Model_initial_variant_calling_v5.py  --library --SSCS_or_DCS --samples_to_exclude

'''''
version = '5.0'

# imported packages
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm
import scipy.special
import scipy.integrate as it
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.stats import kde
from scipy.stats import betabinom
from scipy.stats import binom
import copy
import glob, os
import re
from sklearn import datasets, linear_model
import pandas as pd
from decimal import *
from operator import itemgetter
from collections import OrderedDict
import timeit
import time
import csv
from pyfaidx import Fasta
from argparse import ArgumentParser
import pysam
import sys
import matplotlib
from array import array
from ast import literal_eval

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

def expand_info_depth(info): #extract total variant depth from the variant call file
    info_split = info.split(';')
    depth = info_split[2].split('=')[1]
    var_depth = info_split[4].split('=')[1]
    return depth

def expand_info_variant_depth(info): #extract variant depth from the variant call file
    info_split = info.split(';')
    depth = info_split[2].split('=')[1]
    var_depth = info_split[4].split('=')[1]
    return var_depth

def five_prime(chromosome, position, ref_genome):
    chrom = chromosome.split('chr')[1]
    return str(ref_genome[chrom][(position-1)-1].seq)

def three_prime(chromosome, position, ref_genome):
    chrom = chromosome.split('chr')[1]
    return str(ref_genome[chrom][(position-1)+1].seq)

def create_list_of_cases_and_controls(library, samples_to_exclude):
    cases = []
    controls = []

    sample_names_csv = 'SNV_panel/'+library+'/'+library+'_sample_sheet.csv'
    with open(sample_names_csv) as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if library == 'SLX_19284': #SLX_19284 sample names do not have _SNV on the end
                sample_name = row[0]
                if sample_name not in samples_to_exclude:
                    if sample_name+'_SNV' not in samples_to_exclude:
                        if 'C92' in sample_name:
                            cases.append(sample_name)
                        if 'CNTRL' in sample_name:
                            controls.append(sample_name)
            else:
                sample_name = row[0].split('_SNV')[0]
                if sample_name not in samples_to_exclude:
                    if 'C92' in sample_name:
                        cases.append(sample_name)
                    if 'CNTRL' in sample_name:
                        controls.append(sample_name)

    cases_and_controls = cases+controls
    print('number of cases included = ', len(cases))
    print('number of controls included = ', len(controls))

    return cases_and_controls

def list_of_samples_with_multiple_timepoints(cases_and_controls):
    #Create a list of the samples that are duplicates
    sample_abbreviation_names = {}
    for i in cases_and_controls:
        sample_abbrev = i.split('_')[0]+'_'+i.split('_')[1]
        if sample_abbrev in sample_abbreviation_names.keys():
            sample_abbreviation_names[sample_abbrev]+=1
        else:
            sample_abbreviation_names[sample_abbrev]=1
    sample_abbreviation_names

    samples_with_multi_timepoints = []
    for k, v in sample_abbreviation_names.items():
        if v>1:
            samples_with_multi_timepoints.append(k)

    print(samples_with_multi_timepoints)

    return samples_with_multi_timepoints

def create_sample_df_with_names(sample_name, overall_df, SSCS_or_DCS, library, ref_genome):

    bases = ['A', 'C', 'G', 'T']
    complementary_bases = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}

    #Import the variant calls
    sample_abbrev_name = sample_name.split('_')[0]+'_'+sample_name.split('_')[1] #e.g. C92_070_s2 -> C92_070

    vcf_path = 'UKCTOCS_sample_level_results_V2/'+sample_abbrev_name+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_all_positions.vcf'
    vcf_path2 = 'UKCTOCS_sample_level_results_V2/'+sample_abbrev_name+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_all_positions.vcf'

    if SSCS_or_DCS == 'SSCS':
        colnames = ['chromosome', 'position', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO']
    if SSCS_or_DCS == 'DCS':
        colnames = ['chromosome', 'position', 'ID', 'ref', 'alt', 'FILTER', 'INFO']

    try:
        df = pd.read_csv(vcf_path, index_col=False, delimiter = '\t', comment='#', names=colnames)
    except FileNotFoundError:
        df = pd.read_csv(vcf_path2, index_col=False, delimiter = '\t', comment='#', names=colnames)
    # df = df[df['alt']!='-']
    df[sample_name+'_depth']=df['INFO'].apply(expand_info_depth) #extract the total depth
    df[sample_name+'_var_depth']=df['INFO'].apply(expand_info_variant_depth) #extract the variant depth
    df['change']=df['ref']+'>'+df['alt'] #e.g. C>T
    df['position'] = pd.to_numeric(df['position']) #variant position on the chromosome
    df['three_prime'] = df.apply(lambda x: three_prime(x.chromosome, x.position, ref_genome), axis=1) #3' base of the variant
    df['five_prime'] = df.apply(lambda x: five_prime(x.chromosome, x.position, ref_genome), axis=1) #5' base of the variant
    df['trinucleotide']=df['five_prime']+'['+df['change']+']'+df['three_prime'] #e.g. A[C>T]G
    df = df[['chromosome', 'position', 'ref', 'alt', 'change', 'trinucleotide', sample_name+'_depth', sample_name+'_var_depth']]
    df = df[df['ref'].str.len() ==1] #i.e. only include variants were the ref = 1 (will exclude deletions, but not insertions?)

    df[sample_name]=+df[sample_name+'_var_depth']+','+df[sample_name+'_depth']+','+sample_name

    data_list = df.values.tolist()

    list_all_changes = []
    all_changes = {}

    for position in data_list:
        chromosome = position[0]
        pos = position[1]
        ref = position[2]
        alt = position[3]
        change = position[4]
        trinucleotide = position[5]
        fiveprime = trinucleotide[0]
        threeprime = trinucleotide[-1]
        depth = position[6]
        var_depth = position[7]
        info = position[8]

        for base in bases: #C, T, G, A
            if base != ref:
                if base != alt:
                    base_change = ref+'>'+base
                    trinucleotide_change = fiveprime+'['+base_change+']'+threeprime
                    if (chromosome, pos, ref, base, base_change, trinucleotide_change) not in all_changes.keys(): #don't overwrite if alt already there
                        all_changes[(chromosome, pos, ref, base, base_change, trinucleotide_change)]= ('0,'+depth+','+sample_name)
                if base == alt: #if there is a variant there, record it's variant depth
                    all_changes[(chromosome, pos, ref, alt, change, trinucleotide)]=(info)

    for k, v in all_changes.items():
        list_all_changes.append([k[0], k[1], k[2], k[3], k[4], k[5], v])

    df2 = pd.DataFrame(list_all_changes, columns =['chromosome', 'position', 'ref', 'alt', 'change', 'trinucleotide', sample_name])

    merged_df2 = overall_df.merge(df2, how="left", on=["chromosome", "position", "ref", "alt", "change", "trinucleotide"])
    merged_df2 = merged_df2.fillna(0)

    list_all_changes.clear()
    data_list.clear()

    return merged_df2

def create_dataframe(all_positions_df, SSCS_or_DCS, samples, library, ref_genome):
    combined_df = all_positions_df.copy()
    for sample_name in samples:
        print(sample_name)
        new_df = create_sample_df_with_names(sample_name, combined_df, SSCS_or_DCS, library, ref_genome)
        combined_df = new_df

    return combined_df

def annotated_positions(sample_names, SSCS_or_DCS):
    concat_df = pd.DataFrame(columns = ['position ID', 'intronic_exonic', 'variant_type', 'gene', 'AA_change'])
    for sample in sample_names:
        sample_abbrev = sample.split('_')[0]+'_'+sample.split('_')[1]
        annotated = 'UKCTOCS_sample_level_results_V2/'+sample_abbrev+'/SNV/DCS/'+sample+'_SNV_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_annotated.txt'

        df = pd.read_csv(annotated, delimiter = '\t')
        df['position ID'] = df[['chromosome', 'start', 'REF', 'ALT']].astype(str).agg(','.join, axis=1)
        df2 = df[['position ID', 'intronic_exonic', 'variant_type', 'gene', 'AA_change']]
        concat_df = pd.concat([concat_df, df2])

    concat_df = concat_df.sort_values(by = 'position ID')
    concat_df['variant ID'] = concat_df[['gene', 'intronic_exonic', 'variant_type', 'AA_change']].astype(str).agg(' '.join, axis=1)

    all_annotated = pd.DataFrame.to_dict(concat_df.reset_index(), orient = 'index')
    positions_annotated_dict = {}
    for k, v in all_annotated.items():
        pos = v['position ID']
        variant = v['variant ID']
        positions_annotated_dict[pos]=variant

    return positions_annotated_dict

def BB_likelihood(params, var_depths, total_depths):
    epsilon = params[0]
    delta = params[1]
    a = 1/delta
    b = (1-epsilon)/(delta*epsilon)

    log_sample_BBs = []
    sample_BBs = []

    for sample_variant_depth, sample_total_depth in zip(var_depths, total_depths):
        BB = betabinom.pmf(sample_variant_depth, sample_total_depth, a, b) #calculate the beta-binomial likelihood of measuring that number of variant reads, given the sample depth, position error rate and delta
        log_sample_BBs.append(np.log(BB))
        sample_BBs.append(BB)

    model_likelihood_log = np.sum(log_sample_BBs)

    if delta<1e-08/epsilon:
        model_likelihood_log = -100000000000

    if epsilon<0:
        model_likelihood_log = -100000000000

    if epsilon>1:
        model_likelihood_log = -100000000000

    return -model_likelihood_log

def a_b(epsilon, delta):
    a = 1/delta
    b = (1-epsilon)/(delta*epsilon)
    return a, b

def estimate_delta_method_of_moments(variant_reads, total_depth_reads):
    number_samples = len(variant_reads)
    mean_depth = np.mean(total_depth_reads)
    error_rate_estimate = (np.sum(variant_reads))/np.sum(total_depth_reads)

    #calculate the 1st raw sample moment
    k = 1
    m_1 = (sum([var_read**k for var_read in variant_reads]))/number_samples

    #calculate the 2nd raw sample moment
    k = 2
    m_2 = (sum([var_read**k for var_read in variant_reads]))/number_samples

    try:
        alpha = (mean_depth*m_1-m_2)/(mean_depth*((m_2/m_1)-m_1-1)+m_1)
    except ZeroDivisionError:
        alpha = 0

    try:
        delta = 1/alpha
    except ZeroDivisionError:
        delta = 0

    return delta

def beta_binomial_MLE(variant_reads_list, sample_depths, position_error_rate_initial_guess, position_delta_initial_guess):

    initial_guess = [position_error_rate_initial_guess, position_delta_initial_guess]

    #Calculate the most likely epsilon and delta
    optimization = scipy.optimize.minimize(BB_likelihood, initial_guess, args=(variant_reads_list, sample_depths, ),
                                            method='Nelder-Mead', options = {'maxiter': 10000, 'maxfev': 10000})
    optimization_minima = optimization['fun']
    optimization_outcome = optimization['success']
    optimization_message = optimization['message']
    epsilon = optimization['x'][0]
    delta = optimization['x'][1]

    position_MLE_results = (epsilon, delta, optimization_outcome, optimization_message, optimization_minima)

    return position_MLE_results

def remove_high_VAFs(sample_list, high_VAF_threshold):
    real_variants = {}
    remaining_samples = []

    real_variants_names = []

    var_depths = []
    tot_depths = []
    for sample in sample_list:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])

    position_error_rate = np.sum(var_depths)/np.sum(tot_depths)

    for sample in sample_list: #iterate through the list of samples [(sample_name, variant depth, total depth), (sample_name, variant depth, total depth).....]
        sample_name = sample[0]
        sample_variant_depth = sample[1]
        sample_total_depth = sample[2]
        sample_VAF = sample_variant_depth/sample_total_depth
        sample_specific_p_value_treshold = sample[5]

        if sample_VAF >= high_VAF_threshold:
            real_variants[sample_name] = {}
            real_variants[sample_name] = {'variant_depth': sample_variant_depth, 'total_depth': sample_total_depth, 'iteration called at': 'VAF >'+str(high_VAF_threshold),
                                          'p_value': np.nan, 'fitting_method': 'VAF >'+str(high_VAF_threshold), 'p-value treshold': 'not applicable'}
            real_variants_names.append(sample_name)
        else:
            remaining_samples.append((sample_name, sample_variant_depth, sample_total_depth, sample_VAF, 'VAF >'+str(high_VAF_threshold), np.nan, 'VAF >'+str(high_VAF_threshold), sample_specific_p_value_treshold))

    #calculate the error_rate of the samples that weren't called as real:
    var_depths = []
    tot_depths = []
    number_non_zero = 0
    for sample in remaining_samples:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])
        if sample[1]>0:
            number_non_zero+=1

    final_position_error_rate = np.sum(var_depths)/np.sum(tot_depths)

    return remaining_samples, real_variants, position_error_rate, final_position_error_rate, number_non_zero, real_variants_names

def binomial_calling(sample_list, iteration):
    real_variants = {}
    remaining_samples = []

    real_variants_names = []

    var_depths = []
    tot_depths = []
    for sample in sample_list:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])

    position_error_rate = np.sum(var_depths)/np.sum(tot_depths)

    for sample in sample_list: #iterate through the list of samples [(sample_name, variant depth, total depth), (sample_name, variant depth, total depth).....]
        sample_name = sample[0]
        sample_variant_depth = sample[1]
        sample_total_depth = sample[2]
        sample_VAF = sample_variant_depth/sample_total_depth
        sample_p_value_threshold = sample[7]
        #calculate the p-value for each sample
        p_value = binom.sf(sample_variant_depth-1, sample_total_depth, position_error_rate) #subtract 1 from variant depth because should be p-value for >= that number of variant reads
        if (sample_variant_depth >0) and p_value < sample_p_value_threshold: #if p-value < threshold, call as real
            real_variants[sample_name] = {}
            real_variants[sample_name] = {'variant_depth': sample_variant_depth, 'total_depth': sample_total_depth, 'iteration called at': iteration,
                                          'p_value': p_value, 'fitting_method': 'binomial', 'p-value treshold': sample_p_value_threshold}
            real_variants_names.append(sample_name)
        else: #if p-value not <threshold, then call as an error
            remaining_samples.append((sample_name, sample_variant_depth, sample_total_depth, sample_VAF, iteration, p_value, 'binomial', sample_p_value_threshold))

    #calculate the error_rate of the samples that weren't called as real:
    var_depths = []
    tot_depths = []
    number_non_zero = 0
    for sample in remaining_samples:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])
        if sample[1]>0:
            number_non_zero+=1

    if number_non_zero>0:
        final_position_error_rate = np.sum(var_depths)/np.sum(tot_depths)
    else:
        final_position_error_rate = '<'+str(1/np.sum(tot_depths))

    return remaining_samples, real_variants, position_error_rate, final_position_error_rate, number_non_zero, real_variants_names

def beta_binomial_no_iteration_no_exclusion(sample_list):
    real_variants = {}
    error_samples = []

    real_variants_names = []

    sample_list_to_fit = sample_list # do not exclude the most extreme value
#     sample_list_to_fit = sorted(sample_list, key=lambda x: x[3], reverse = True)[1:] #[1:] will exclude the first sample from the list

    #calculate an error rate estimate to use for the MLE
    var_depths = []
    tot_depths = []
    for sample in sample_list_to_fit:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])

    mean_depth = np.mean(tot_depths)

    error_rate_estimate = np.sum(var_depths)/np.sum(tot_depths)
    delta_estimate = estimate_delta_method_of_moments(var_depths, tot_depths)
    if delta_estimate < 0:
        delta_estimate = 1e-04/error_rate_estimate

    #Calculate the MLE for epsilon and delta
    MLE = beta_binomial_MLE(var_depths, tot_depths, error_rate_estimate, delta_estimate)
#     print(MLE)
    epsilon = MLE[0]
    delta = MLE[1]
    outcome = MLE[2]
    comment = MLE[3]
    minima = MLE[4]
    a, b = a_b(epsilon, delta)

    #Using the inferred value for epsilon and delta, calculate the p-value for each sample...
    for sample in sample_list: #sample list includes the excluded sample
        sample_name = sample[0]
        sample_variant_depth = sample[1]
        sample_total_depth = sample[2]
        sample_VAF = sample_variant_depth/sample_total_depth
        sample_p_value_threshold = sample[7]
        p_value = betabinom.sf(sample_variant_depth-1, sample_total_depth, a, b) #subtract 1 from variant depth because should be p-value for >= that number of variant reads

        if p_value == 0: #i.e. if the p-value is below the ability of the betabinom.sf's ability to determine it...do manually...
            for i in range(sample_variant_depth-1, sample_total_depth):
                pmf = betabinom.pmf(i, sample_total_depth, a, b)
                p_value+=pmf
                if pmf < p_value/1000:
                    break

        if (sample_variant_depth >0) and (p_value < sample_p_value_threshold):
            real_variants[sample_name] = {'variant_depth': sample_variant_depth, 'total_depth': sample_total_depth, 'iteration called at': 'no iteration',
                                          'p_value': p_value, 'fitting_method': 'beta-binomial', 'p-value treshold': sample_p_value_threshold}
            real_variants_names.append(sample_name)
        #if the p-value isn't less than the threshold, add the sample to a 'remaining variants' list
        else:
            error_samples.append((sample_name, sample_variant_depth, sample_total_depth, sample_VAF, 'no iteration', p_value, 'beta-binomial', sample_p_value_threshold))


    #calculate the error_rate and delta of the samples that weren't called as real:
    var_depths = []
    tot_depths = []
    number_non_zero = 0
    for sample in error_samples:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])
        if sample[1]>0:
            number_non_zero+=1

    mean_depth = np.mean(tot_depths)

    if number_non_zero>3: #fit a beta-binomial
        final_error_rate_estimate = np.sum(var_depths)/np.sum(tot_depths)
        final_delta_estimate = estimate_delta_method_of_moments(var_depths, tot_depths)
        if final_delta_estimate < 0:
            final_delta_estimate = 1e-04/final_error_rate_estimate

        MLE = beta_binomial_MLE(var_depths, tot_depths, final_error_rate_estimate, final_delta_estimate)
        final_epsilon = MLE[0]
        final_delta = MLE[1]
        final_outcome = MLE[2]
        final_comment = MLE[3]
        final_minima = MLE[4]

    if number_non_zero<=3: #use a binomial
        if number_non_zero == 0:
            final_error_rate_estimate = '<'+str(1/np.sum(tot_depths))
            final_epsilon = '<'+str(1/np.sum(tot_depths))
            final_delta_estimate = np.nan
            final_delta = np.nan
            final_outcome = np.nan
            final_comment = np.nan
            final_minima = np.nan
        else:
            final_error_rate_estimate = np.sum(var_depths)/np.sum(tot_depths)
            final_epsilon = np.sum(var_depths)/np.sum(tot_depths)
            final_delta_estimate = np.nan
            final_delta = np.nan
            final_outcome = np.nan
            final_comment = np.nan
            final_minima = np.nan

    return error_samples, real_variants, outcome, comment, epsilon, delta, error_rate_estimate, delta_estimate, minima, final_epsilon, final_delta, final_error_rate_estimate, final_delta_estimate, number_non_zero, real_variants_names

def beta_binomial_iteration(sample_list, iteration_number):
    real_variants = {}
    error_samples = []
    real_variants_names = []

    if iteration_number == 0:
        #for the first iteration (iteration 0), do not exclude the most extreme value
        sample_list_to_fit = sample_list
    if iteration_number > 0:
        #if not the first iteration, sort the sample list from highest to lowest VAF and then exclude the highest samples
        sample_list_to_fit = sorted(sample_list, key=lambda x: x[3], reverse = True)[1:] #[1:] will exclude the first sample from the list

    #calculate an error rate estimate to use for the MLE
    var_depths = []
    tot_depths = []
    for sample in sample_list_to_fit:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])

    mean_depth = np.mean(tot_depths)

    error_rate_estimate = np.sum(var_depths)/np.sum(tot_depths)
    delta_estimate = estimate_delta_method_of_moments(var_depths, tot_depths)
    if delta_estimate < 0:
        delta_estimate = 1e-04/error_rate_estimate

    #Calculate the MLE for epsilon and delta
    MLE = beta_binomial_MLE(var_depths, tot_depths, error_rate_estimate, delta_estimate)
#     print(MLE)
    epsilon = MLE[0]
    delta = MLE[1]
    outcome = MLE[2]
    comment = MLE[3]
    minima = MLE[4]
    a, b = a_b(epsilon, delta)

    #Using the inferred value for epsilon and delta, calculate the p-value for each sample (including the extreme excluded sample, if there was one)...
    for sample in sample_list: #sample list includes the excluded sample
        sample_name = sample[0]
        sample_variant_depth = sample[1]
        sample_total_depth = sample[2]
        sample_VAF = sample_variant_depth/sample_total_depth
        sample_p_value_threshold = sample[7] #sample position-specifc p-value (if the variant is present in the final timepoint sample then the p-value threshold used is 0.1)
        # print('sample_p_value_threshold = ', sample_p_value_threshold)
        p_value = betabinom.sf(sample_variant_depth-1, sample_total_depth, a, b) #subtract 1 from variant depth because should be p-value for >= that number of variant reads

        if p_value == 0: #i.e. if the p-value is below the ability of the betabinom.sf's ability to determine it...do manually...
            for i in range(sample_variant_depth-1, sample_total_depth):
                pmf = betabinom.pmf(i, sample_total_depth, a, b)
                p_value+=pmf
                if pmf < p_value/1000:
                    break

        if iteration_number==0: #don't call any variants as real with this first iteration
            error_samples.append((sample_name, sample_variant_depth, sample_total_depth, sample_VAF, iteration_number, p_value, 'beta-binomial (no calls at this iteration)', sample_p_value_threshold))
        if iteration_number >0:
            #if the p-value is less than the threshold, record that sample as a real variant
            if (sample_variant_depth >0) and (p_value < sample_p_value_threshold):
                real_variants[sample_name] = {'variant_depth': sample_variant_depth, 'total_depth': sample_total_depth, 'iteration called at': iteration_number,
                                              'p_value': p_value, 'fitting_method': 'beta-binomial', 'p-value threshold': sample_p_value_threshold}
                real_variants_names.append(sample_name)
            #if the p-value isn't less than the threshold, add the sample to a 'remaining variants' list
            else:
                error_samples.append((sample_name, sample_variant_depth, sample_total_depth, sample_VAF, iteration_number, p_value, 'beta-binomial', sample_p_value_threshold))

    #calculate the error_rate and delta of the samples that weren't called as real:
    var_depths = []
    tot_depths = []
    number_non_zero = 0
    for sample in error_samples:
        var_depths.append(sample[1])
        tot_depths.append(sample[2])
        if sample[1]>0:
            number_non_zero+=1

    mean_depth = np.mean(tot_depths)

    if number_non_zero>3:
        final_error_rate_estimate = np.sum(var_depths)/np.sum(tot_depths)
        final_delta_estimate = estimate_delta_method_of_moments(var_depths, tot_depths)
        if final_delta_estimate < 0:
            final_delta_estimate = 1e-04/final_error_rate_estimate

        MLE = beta_binomial_MLE(var_depths, tot_depths, final_error_rate_estimate, final_delta_estimate)
        final_epsilon = MLE[0]
        final_delta = MLE[1]
        final_outcome = MLE[2]
        final_comment = MLE[3]
        final_minima = MLE[4]

    if number_non_zero<=3:
        if number_non_zero == 0:
            final_error_rate_estimate = '<'+str(1/np.sum(tot_depths))
            final_delta_estimate = np.nan
            final_epsilon = np.sum(var_depths)/np.sum(tot_depths)
            final_delta = np.nan
            final_outcome = np.nan
            final_comment = np.nan
            final_minima = np.nan
        else:
            final_error_rate_estimate = np.sum(var_depths)/np.sum(tot_depths)
            final_delta_estimate = np.nan
            final_epsilon = np.sum(var_depths)/np.sum(tot_depths)
            final_delta = np.nan
            final_outcome = np.nan
            final_comment = np.nan
            final_minima = np.nan

    return error_samples, real_variants, outcome, comment, epsilon, delta, error_rate_estimate, delta_estimate, minima, final_epsilon, final_delta, final_error_rate_estimate, final_delta_estimate, number_non_zero, real_variants_names

def iterative_variant_calling(sample_list, number_non_zero, high_VAF_threshold):

    positions_MLEs = {}
    positions_real_variants = []

    positions_p_values_errors = {}
    positions_p_values_real_variants = {}

    total_variants_called = 0

    iterative_epsilons = []
    iterative_deltas = []
    iterative_outcomes = []
    iterative_comments = []
    iterative_variants_called = []
    iterative_fitting_methods = []
    initialisations = []

    #first exclude high VAF variants...
    remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = remove_high_VAFs(sample_list, high_VAF_threshold)

    if len(real_variants)>0: #i.e. if high VAF variants have been detected
        total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
        positions_real_variants.append(real_variants) #keep a list of the real variants
        final_delta = np.nan
        fitting_method = 'VAF >'+str(high_VAF_threshold)
        sample_list = remaining_samples #exclude the high VAF variants

        positions_p_values_errors['VAF >'+str(high_VAF_threshold)]=remaining_samples

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(np.nan)
        iterative_outcomes.append(np.nan)
        iterative_comments.append(np.nan)
        iterative_variants_called.append(len(real_variants))
        iterative_fitting_methods.append(fitting_method)
        initialisations.append((np.nan, np.nan))

    #if, after high VAF exclusion, <= 3 are non-zero variant reads, then fit a binomial (error rate = sum variant reads/ total variant reads)
    if number_non_zero<=3: #fit a binomial
        iteration_number = 1 #don't need to start at zero because not doing any exclusions here

        #don't exclude any samples when fit the binomial because it's not like the beta-binomial that overfits dispersion
        remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(remaining_samples, iteration_number)
        final_delta = np.nan #because binomial
        fitting_method = 'binomial'

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(np.nan)
        iterative_outcomes.append(np.nan)
        iterative_comments.append(np.nan)
        iterative_variants_called.append(len(real_variants))
        iterative_fitting_methods.append(fitting_method)
        initialisations.append((epsilon, np.nan))

        positions_p_values_errors['iteration 1']=remaining_samples
        positions_p_values_real_variants['iteration 1']=real_variants #samples and their p-values at this iteration

        while len(real_variants)>0: #i.e. if a real variant has been detected, keep iterating...
            total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
            positions_real_variants.append(real_variants) #keep a list of the real variants
            sample_list = remaining_samples #for the next iteration, only look at the samples that haven't been called as real
            iteration_number+=1

            if number_non_zero>0: #if there is more than 0 non-zero sample remaining that hasn't been caled as real
                remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(sample_list, iteration_number)
                iterative_fitting_methods.append('binomial')

            else: #if only non-zero sample remaining...
                real_variants = [] #this will stop the iteration, because no real variants being called
                epsilon = final_epsilon
                delta = np.nan
                final_delta = np.nan
                outcome = np.nan
                comment = np.nan
                iterative_fitting_methods.append('binomial')
                epsilon_initialisation = epsilon
                delta_initialisation = np.nan
                final_epsilon_initialisation = final_epsilon
                final_delta_initialisation = np.nan

            iterative_epsilons.append(epsilon)
            iterative_deltas.append(np.nan)
            iterative_outcomes.append(np.nan)
            iterative_comments.append(np.nan)
            iterative_variants_called.append(len(real_variants))
            initialisations.append((epsilon, np.nan))

            positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
            positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration

        final_epsilon_initialisation = epsilon
        final_delta_initialisation = np.nan


    #if more than one non-zero variant read, exclude the most extreme and fit a binomial (calculate p-values for both)
    if number_non_zero>3:
        final_fitting_method = 'beta-binomial'

        #include all samples at the position to start with...(but don't call real variants (i.e. real variants list will be empty))
        iteration_number = 0
        remaining_samples, real_variants, outcome, comment, epsilon, delta, epsilon_initialisation, delta_initialisation, minima, final_epsilon, final_delta, final_epsilon_initialisation, final_delta_initialisation, number_non_zero, real_variants_names = beta_binomial_iteration(remaining_samples, iteration_number)

        positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
        positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(delta)
        iterative_outcomes.append(outcome)
        iterative_comments.append(comment)
        iterative_variants_called.append(np.nan) #this will be 0
        iterative_fitting_methods.append('beta-binomial')
        initialisations.append((epsilon_initialisation, delta_initialisation))

        #after the 0th iteration, exclude the most extreme position and check for real variants again... (even if variant was detected at 0th iteration)
        iteration_number = 1
        remaining_samples, real_variants, outcome, comment, epsilon, delta, epsilon_initialisation, delta_initialisation, minima, final_epsilon, final_delta, final_epsilon_initialisation, final_delta_initialisation, number_non_zero, real_variants_names = beta_binomial_iteration(remaining_samples, iteration_number)

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(delta)
        iterative_outcomes.append(outcome)
        iterative_comments.append(comment)
        iterative_variants_called.append(len(real_variants))
        initialisations.append((epsilon_initialisation, delta_initialisation))

        positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
        positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration


        #if the extreme position was deemed to be real, then take the new list (without the real variant) and exclude the next most extreme and see if another real variant is called...
        while len(real_variants)>0: #i.e. if a real variant has been detected, keep iterating...
            total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
            positions_real_variants.append(real_variants) #keep a list of the real variants
            sample_list = remaining_samples #for the next iteration, only look at the samples that haven't been called as real
            iteration_number+=1

            if number_non_zero>3: #if there is more than 1 non-zero sample remaining that hasn't been caled as real
                remaining_samples, real_variants, outcome, comment, epsilon, delta, epsilon_initialisation, delta_initialisation, minima, final_epsilon, final_delta, final_epsilon_initialisation, final_delta_initialisation, number_non_zero, real_variants_names = beta_binomial_iteration(remaining_samples, iteration_number)
                iterative_fitting_methods.append('beta-binomial')

            else: #if <=3 non-zero sample remaining, then revert to binomial fitting
                if number_non_zero>0:
                    remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(remaining_samples, iteration_number)
                    delta = np.nan
                    final_delta = np.nan
                    outcome = np.nan
                    comment = np.nan
                    iterative_fitting_methods.append('binomial')
                    epsilon_initialisation = epsilon
                    delta_initialisation = np.nan
                    final_epsilon_initialisation = epsilon
                    final_delta_initialisation = np.nan

                else: # if all zero variant reads:
                    real_variants = [] #this will stop the iteration, because no real variants being called
                    tot_depths = []
                    for sample in remaining_samples:
                        tot_depths.append(sample[2])

                    final_epsilon = '<'+str(1/np.sum(tot_depths))

                    real_variants = []
                    epsilon = final_epsilon
                    delta = np.nan
                    final_delta = np.nan
                    outcome = np.nan
                    comment = np.nan
                    iterative_fitting_methods.append('binomial')
                    epsilon_initialisation = epsilon
                    delta_initialisation = np.nan
                    final_epsilon_initialisation = final_epsilon
                    final_delta_initialisation = np.nan

            iterative_epsilons.append(epsilon)
            iterative_deltas.append(delta)
            iterative_outcomes.append(outcome)
            iterative_comments.append(comment)
            iterative_variants_called.append(len(real_variants))
            initialisations.append((epsilon_initialisation, delta_initialisation))

            positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
            positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration

    positions_MLEs= {'iterative fitting methods': iterative_fitting_methods, 'final position error rate': final_epsilon, 'final delta': final_delta,
                                              'total variants called': total_variants_called,
                                             'epsilons': iterative_epsilons, 'deltas': iterative_deltas, 'MLE outcomes': iterative_outcomes,
                                             'MLE comments': iterative_comments, 'iterative variants called': iterative_variants_called,
                    'initialisation epsilons and deltas': initialisations, 'final initialisation epsilon and delta': (final_epsilon_initialisation,
                     final_delta_initialisation)}

    return positions_MLEs, positions_real_variants, total_variants_called, positions_p_values_errors, positions_p_values_real_variants

def iterative_variant_calling_for_multi_timepoints(sample_list, number_non_zero, high_VAF_threshold, samples_to_call_on_iterations):

    #samples_to_call_iteration_on = a list of the abbreviated name of the samples (e.g. C92_007) that are on this lane from multiple timepoints (for which one of them has been called real) and which should call iteratively

    positions_MLEs = {}
    positions_real_variants = []

    positions_p_values_errors = {}
    positions_p_values_real_variants = {}

    total_variants_called = 0

    iterative_epsilons = []
    iterative_deltas = []
    iterative_outcomes = []
    iterative_comments = []
    iterative_variants_called = []
    iterative_fitting_methods = []
    initialisations = []

    #first exclude high VAF variants...
    remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = remove_high_VAFs(sample_list, high_VAF_threshold)

    if len(real_variants)>0: #i.e. if high VAF variants have been detected
        total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
        positions_real_variants.append(real_variants) #keep a list of the real variants
        final_delta = np.nan
        fitting_method = 'VAF >'+str(high_VAF_threshold)
        sample_list = remaining_samples #exclude the high VAF variants

        positions_p_values_errors['VAF >'+str(high_VAF_threshold)]=remaining_samples

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(np.nan)
        iterative_outcomes.append(np.nan)
        iterative_comments.append(np.nan)
        iterative_variants_called.append(len(real_variants))
        iterative_fitting_methods.append(fitting_method)
        initialisations.append((np.nan, np.nan))

    #if, after high VAF exclusion, all non-zero variant reads, except 1, then fit a binomial (error rate = sum variant reads/ total variant reads)
    if number_non_zero<=3: #fit a binomial
        iteration_number = 1
        remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(remaining_samples, iteration_number)

        final_delta = np.nan #because binomial
        fitting_method = 'binomial'

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(np.nan)
        iterative_outcomes.append(np.nan)
        iterative_comments.append(np.nan)
        iterative_variants_called.append(len(real_variants))
        iterative_fitting_methods.append(fitting_method)
        initialisations.append((epsilon, np.nan))
        final_epsilon_initialisation = epsilon
        final_delta_initialisation = np.nan

        positions_p_values_errors['iteration 1']=remaining_samples
        positions_p_values_real_variants['iteration 1']=real_variants #samples and their p-values at this iteration

        samples_to_iterate_on_called_as_real = 0
        if len(real_variants)>0: #i.e. if a real variant has been detected, keep iterating, but only if it is in the list of multi-timepoint samples
            total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
            positions_real_variants.append(real_variants) #keep a list of the real variants
            sample_list = remaining_samples #for the next iteration, only look at the samples that haven't been called as real
            iteration_number+=1

            #check if the samples that have been called as real are the multi-timepoint samples
            real_variants_abbrev_names = []
            for i in real_variants_names:
                real_variants_abbrev_names.append(i.split('_')[0]+'_'+i.split('_')[1]) #the abbreviated names of the samples that have been called as real (e.g. C92_007)

            for i in samples_to_call_on_iterations:
                if i in real_variants_abbrev_names:
                    samples_to_iterate_on_called_as_real+=1

        #only carry on iterating if one of the multi-timepoint samples has been called as real...
        while samples_to_iterate_on_called_as_real >0:
            sample_list = remaining_samples #for the next iteration, only look at the samples that haven't been called as real
            iteration_number+=1

            if number_non_zero>0: #if there is more than 1 non-zero sample remaining that hasn't been caled as real
                remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(remaining_samples, iteration_number)
                iterative_fitting_methods.append('binomial')

                if len(real_variants)>0:
                    total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
                    positions_real_variants.append(real_variants) #keep a list of the real variants

                    #check if the samples that have been called as real are the multi-timepoint samples
                    real_variants_abbrev_names = []
                    for i in real_variants_names:
                        real_variants_abbrev_names.append(i.split('_')[0]+'_'+i.split('_')[1]) #the abbreviated names of the samples that have been called as real (e.g. C92_007)

                    samples_to_iterate_on_called_as_real = 0
                    for i in samples_to_call_on_iterations:
                        if i in real_variants_abbrev_names:
                            samples_to_iterate_on_called_as_real+=1

                else:
                    samples_to_iterate_on_called_as_real = 0


            else:#if only non-zero sample remaining...
                samples_to_iterate_on_called_as_real = 0 #this will stop the iteration, because no real variants being called
                real_variants = []
                epsilon = final_epsilon
                delta = np.nan
                final_delta = np.nan
                outcome = np.nan
                comment = np.nan
                iterative_fitting_methods.append('binomial')
                epsilon_initialisation = epsilon
                delta_initialisation = np.nan
                final_epsilon_initialisation = final_epsilon
                final_delta_initialisation = np.nan

            iterative_epsilons.append(epsilon)
            iterative_deltas.append(np.nan)
            iterative_outcomes.append(np.nan)
            iterative_comments.append(np.nan)
            iterative_variants_called.append(len(real_variants))
            initialisations.append((epsilon, np.nan))

            positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
            positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration

        final_epsilon_initialisation = epsilon
        final_delta_initialisation = np.nan


    #if more than one non-zero variant read, fit a beta-binomial (calculate p-values for both)
    if number_non_zero>3:
        final_fitting_method = 'beta-binomial'

        #check for real variants...
        iteration_number = 1
        remaining_samples, real_variants, outcome, comment, epsilon, delta, epsilon_initialisation, delta_initialisation, minima, final_epsilon, final_delta, final_epsilon_initialisation, final_delta_initialisation, number_non_zero, real_variants_names = beta_binomial_iteration(remaining_samples, iteration_number)

        iterative_epsilons.append(epsilon)
        iterative_deltas.append(delta)
        iterative_outcomes.append(outcome)
        iterative_comments.append(comment)
        iterative_variants_called.append(len(real_variants))
        initialisations.append((epsilon_initialisation, delta_initialisation))

        positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
        positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration

        #if a real variant was detected from one of the samples that has multiple timepoints at this position, then iterate again to see if any of the other timepoint are called as real...
        samples_to_iterate_on_called_as_real = 0
        if len(real_variants)>0: #i.e. if a real variant has been detected, keep iterating...
            total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
            positions_real_variants.append(real_variants) #keep a list of the real variants

            #check if the samples that have been called as real are the multi-timepoint samples
            real_variants_abbrev_names = []
            for i in real_variants_names:
                real_variants_abbrev_names.append(i.split('_')[0]+'_'+i.split('_')[1]) #the abbreviated names of the samples that have been called as real (e.g. C92_007)

            for i in samples_to_call_on_iterations:
                if i in real_variants_abbrev_names:
                    samples_to_iterate_on_called_as_real+=1


        #only carry on iterating if one of the multi-timepoint samples has been called as real...
        while samples_to_iterate_on_called_as_real >0:
            sample_list = remaining_samples #for the next iteration, only look at the samples that haven't been called as real
            iteration_number+=1

            if number_non_zero<=3: #if 3 or less non-zero sample remaining, then revert to binomial fitting
                if number_non_zero>0: #fit a binomial
                    remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(remaining_samples, iteration_number)
                    if len(real_variants)>0:
                        total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
                        positions_real_variants.append(real_variants) #keep a list of the real variants

                    delta = np.nan
                    final_delta = np.nan
                    outcome = np.nan
                    comment = np.nan
                    iterative_fitting_methods.append('binomial')
                    epsilon_initialisation = epsilon
                    delta_initialisation = np.nan
                    final_epsilon_initialisation = epsilon
                    final_delta_initialisation = np.nan

                    #check if the samples that have been called as real are the multi-timepoint samples
                    real_variants_abbrev_names = []
                    for i in real_variants_names:
                        real_variants_abbrev_names.append(i.split('_')[0]+'_'+i.split('_')[1]) #the abbreviated names of the samples that have been called as real (e.g. C92_007)

                    samples_to_iterate_on_called_as_real = 0
                    for i in samples_to_call_on_iterations:
                        if i in real_variants_abbrev_names:
                            samples_to_iterate_on_called_as_real+=1

                else:
                    samples_to_iterate_on_called_as_real = 0 #this will stop the iteration, because no real variants being called
                    real_variants = []
                    epsilon = final_epsilon
                    delta = np.nan
                    final_delta = np.nan
                    outcome = np.nan
                    comment = np.nan
                    iterative_fitting_methods.append('binomial')
                    epsilon_initialisation = epsilon
                    delta_initialisation = np.nan
                    final_epsilon_initialisation = final_epsilon
                    final_delta_initialisation = np.nan

            if number_non_zero>3: #if there is more than 1 non-zero sample remaining that hasn't been caled as real
                remaining_samples, real_variants, outcome, comment, epsilon, delta, epsilon_initialisation, delta_initialisation, minima, final_epsilon, final_delta, final_epsilon_initialisation, final_delta_initialisation, number_non_zero, real_variants_names = beta_binomial_iteration(remaining_samples, iteration_number)
                if len(real_variants)>0:
                    total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
                    positions_real_variants.append(real_variants) #keep a list of the real variants
                iterative_fitting_methods.append('beta-binomial')

                #check if the samples that have been called as real are the multi-timepoint samples
                real_variants_abbrev_names = []
                for i in real_variants_names:
                    real_variants_abbrev_names.append(i.split('_')[0]+'_'+i.split('_')[1]) #the abbreviated names of the samples that have been called as real (e.g. C92_007)

                samples_to_iterate_on_called_as_real = 0
                for i in samples_to_call_on_iterations:
                    if i in real_variants_abbrev_names:
                        samples_to_iterate_on_called_as_real+=1

            iterative_epsilons.append(epsilon)
            iterative_deltas.append(delta)
            iterative_outcomes.append(outcome)
            iterative_comments.append(comment)
            iterative_variants_called.append(len(real_variants))
            initialisations.append((epsilon_initialisation, delta_initialisation))

            positions_p_values_errors['iteration '+str(iteration_number)]=remaining_samples #samples and their p-values at this iteration
            positions_p_values_real_variants['iteration '+str(iteration_number)]=real_variants #samples and their p-values at this iteration

    positions_MLEs= {'iterative fitting methods': iterative_fitting_methods, 'final position error rate': final_epsilon, 'final delta': final_delta,
                                              'total variants called': total_variants_called,
                                             'epsilons': iterative_epsilons, 'deltas': iterative_deltas, 'MLE outcomes': iterative_outcomes,
                                             'MLE comments': iterative_comments, 'iterative variants called': iterative_variants_called,
                    'initialisation epsilons and deltas': initialisations, 'final initialisation epsilon and delta': (final_epsilon_initialisation,
                     final_delta_initialisation)}

    return positions_MLEs, positions_real_variants, total_variants_called, positions_p_values_errors, positions_p_values_real_variants

def non_iterative_variant_calling_no_exclusion_for_multi_timepoints(sample_list, number_non_zero, high_VAF_threshold):

    positions_MLEs = {}
    positions_real_variants = []

    positions_p_values_errors = {}
    positions_p_values_real_variants = {}

    total_variants_called = 0

    list_of_names_of_real_variants = [] #keep a simple list of the names of the samples that are called as real (i.e. just the names, not the p-value info etc.)

    #first exclude high VAF variants...
    remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = remove_high_VAFs(sample_list, high_VAF_threshold)

    if len(real_variants)>0: #i.e. if high VAF variants have been detected
        total_variants_called+=len(real_variants) #keep a running tally of how many real variants called at this position
        positions_real_variants.append(real_variants) #keep a list of the real variants
        fitting_method = 'VAF >'+str(high_VAF_threshold)
        sample_list = remaining_samples #exclude the high VAF variants

        positions_p_values_errors['VAF >'+str(high_VAF_threshold)]=remaining_samples

        for i in real_variants_names:
            list_of_names_of_real_variants.append(i)


    if number_non_zero<=3: #fit a binomial
        remaining_samples, real_variants, epsilon, final_epsilon, number_non_zero, real_variants_names = binomial_calling(remaining_samples, 1)
        final_delta = np.nan #because binomial
        outcome = np.nan
        comment = np.nan

        total_variants_called+=len(real_variants)
        positions_real_variants.append(real_variants)
        fitting_method = 'binomial'

        positions_p_values_errors['no iteration']=remaining_samples
        positions_p_values_real_variants['no iteration']=real_variants #samples and their p-values at this iteration

        for i in real_variants_names:
            list_of_names_of_real_variants.append(i)

    #if more than one non-zero variant read...
    if number_non_zero>3:
        remaining_samples, real_variants, outcome, comment, epsilon, delta, epsilon_initialisation, delta_initialisation, minima, final_epsilon, final_delta, final_epsilon_initialisation, final_delta_initialisation, number_non_zero, real_variants_names = beta_binomial_no_iteration_no_exclusion(remaining_samples)

        total_variants_called+=len(real_variants)
        positions_real_variants.append(real_variants) #keep a list of the real variants
        fitting_method = 'beta-binomial'

        positions_p_values_errors['no iteration']=remaining_samples
        positions_p_values_real_variants['no iteration']=real_variants #samples and their p-values at this iteration

        for i in real_variants_names:
            list_of_names_of_real_variants.append(i)

    positions_MLEs= {'fitting methods': fitting_method, 'final position error rate': final_epsilon, 'final delta': final_delta,
                                              'total variants called': total_variants_called,
                                             'MLE outcomes': outcome,
                                             'MLE comments': comment}

    return positions_MLEs, positions_real_variants, total_variants_called, positions_p_values_errors, positions_p_values_real_variants, list_of_names_of_real_variants

def noniterative_or_iterative_beta_binomial_binomial_variant_calling_for_multi_timepoints(df, high_VAF_threshold, p_value_threshold, COSMIC_positions, multi_timepoint_samples):

    data_values = df.values.tolist()

    n = 0
    start_time = time.time()
    reset_time = time.time()

    position_MLEs = {}
    position_real_variants = {}

    position_error_sample_p_values = {}
    position_real_sample_p_values = {}

    position_iterative_or_not = {}

    total_variants_called = 0
    total_positions_analysed = 0

    for position in data_values:
        chromosome = position[0]
        pos = position[1]
        ref = position[2]
        alt = position[3]

        position_label = (chromosome, pos, ref, alt)

        number_samples = 0
        number_non_zero = 0
        position_sample_list = []

        m=0
        for sample in position[6:]:
            if sample not in ['0', 0]: #i.e. if that sample had coverage at the position
                number_samples+=1
                var_depth = int(sample.split(',')[0])
                tot_depth = int(sample.split(',')[1])
                sample_name = sample.split(',')[2]
                sample_name_abbrev = sample_name.split('_')[0]+'_'+sample_name.split('_')[1]
                VAF = var_depth/tot_depth
#                 if sample_name_abbrev in final_timepoint_variants.keys():
#                     if position_label in final_timepoint_variants[sample_name_abbrev].keys():
#                         p_value_to_use = 0.1
#                     else:
#                         p_value_to_use = p_value_threshold
#                 else:
#                     p_value_to_use = p_value_threshold
                position_sample_list.append((sample_name, var_depth, tot_depth, VAF, sample_name_abbrev, p_value_threshold)) #sample-specific p-value threshold is stored in the position sample list
                if var_depth!=0:
                    number_non_zero+=1
            m+=1

        iterative_or_not = ''

        if number_samples>0: #only look at the positions that have coverage in at least 1 sample
            total_positions_analysed+=1
            if position_label in COSMIC_positions.keys(): #do iterative method if it is a position seen in COSMIC haem/ lymphoid (don't need to worry about samples with multi-timepoints here because already iterating)
                print('COSMIC position: '+str(position_label))
                iterative_or_not = 'iterative'
                MLE_outcomes, real_variants, number_variants_called, p_values_errors, p_values_real = iterative_variant_calling(position_sample_list, number_non_zero, high_VAF_threshold)
                position_MLEs[position_label]= MLE_outcomes
                position_real_variants[position_label] = real_variants
                position_error_sample_p_values[position_label]=p_values_errors
                position_real_sample_p_values[position_label]=p_values_real
                total_variants_called+=number_variants_called

            if position_label not in COSMIC_positions.keys(): # don't do iterative approach if not a position seen in COSMIC haem/ lymphoid
                iterative_or_not = 'non-iterative'
                MLE_outcomes, real_variants, number_variants_called, p_values_errors, p_values_real, list_of_names_of_real_variants = non_iterative_variant_calling_no_exclusion_for_multi_timepoints(position_sample_list, number_non_zero, high_VAF_threshold)

                #create a list of the samples that are from multi-timepoints and for which one of them has been called as real
                #first, get the abbreviated name of the samples that were called as real
                real_variants_abbrev_names = []
                for i in list_of_names_of_real_variants:
                    real_variants_abbrev_names.append(i.split('_')[0]+'_'+i.split('_')[1]) #i.e. C92_007_s2 -> C92_007

                samples_to_call_iteration_on = []
                for i in multi_timepoint_samples: #this is the abbreviated sample name, e.g. 'C92_007'
                    if i in real_variants_abbrev_names:
                        samples_to_call_iteration_on.append(i) # a list of the names of the samples that are from multi-timepoints and for which one of them has been called as real, e.g. ['C92_007']

                if len(samples_to_call_iteration_on)>0:
                    MLE_outcomes, real_variants, number_variants_called, p_values_errors, p_values_real = iterative_variant_calling_for_multi_timepoints(position_sample_list, number_non_zero, high_VAF_threshold, samples_to_call_iteration_on)

                #if none of the multi-timepoint samples have called real variants, then don't need to iterate

                position_MLEs[position_label]= MLE_outcomes
                position_real_variants[position_label] = real_variants
                position_error_sample_p_values[position_label]=p_values_errors
                position_real_sample_p_values[position_label]=p_values_real
                total_variants_called+=number_variants_called

        position_iterative_or_not[position_label]=iterative_or_not

        n+=1

        if n%10==0:
            time_now = time.time()
            print('total positions fitted = '+str(n)+' (last 10 in '+str(int(time_now - reset_time))+' seconds).  Total time elapsed = '+str(round((time_now - start_time)/60, 2))+' minutes')
            reset_time = time.time()

    return position_MLEs, position_real_variants, total_variants_called, position_iterative_or_not, total_positions_analysed, position_error_sample_p_values, position_real_sample_p_values

def save_MLE_outcomes(position_MLEs, merged_BB, library, SSCS_or_DCS):
    position_MLEs_rekeyed = {}
    for k, v in position_MLEs.items():
        position_ID = k[0]+','+str(k[1])+','+str(k[2])+','+str(k[3])
        position_MLEs_rekeyed[position_ID]=v

    MLE_df = pd.DataFrame.from_dict(position_MLEs_rekeyed, orient = 'index')
    MLE_df = MLE_df.replace('np.nan', np.nan)

    merged_BB_MLE = pd.merge(merged_BB, MLE_df, how = 'left', left_index = True, right_index = True)

    return merged_BB_MLE.to_csv('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_model_MLE_outcome_per_position.csv')

def save_error_sample_p_values(position_error_sample_p_values, library, SSCS_or_DCS): #p-values for the samples that were called as errors
    position_error_sample_p_values_rekeyed = {}
    for k, v in position_error_sample_p_values.items():
        position_ID = k[0]+','+str(k[1])+','+str(k[2])+','+str(k[3])
        position_error_sample_p_values_rekeyed[position_ID]=v

    errors_df = pd.DataFrame.from_dict(position_error_sample_p_values_rekeyed, orient = 'index')
    errors_df = errors_df.replace('np.nan', np.nan)

    return errors_df.to_csv('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_p_values_for_variants_called_as_errors.csv')

def save_real_sample_p_values(position_real_sample_p_values, library, SSCS_or_DCS): #p-values for the samples that were called as real
    position_real_sample_p_values_rekeyed = {}
    for k, v in position_real_sample_p_values.items():
        position_ID = k[0]+','+str(k[1])+','+str(k[2])+','+str(k[3])
        position_real_sample_p_values_rekeyed[position_ID]=v

    reals_df = pd.DataFrame.from_dict(position_real_sample_p_values_rekeyed, orient = 'index')
    reals_df = reals_df.replace('np.nan', np.nan)

    return reals_df.to_csv('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_p_values_for_variants_called_as_real.csv')

def save_all_variant_calls(position_real_variants, position_error_p_values, position_final_error_rate, SSCS_or_DCS, library):

    real_variants = {}
    for k, v in position_real_variants.items():
        if len(v)>0: #i.e. if there are variants
            real_variants[k]=v

    #create a sample based dictionary of DCS real variants
    sample_variants = {}

    for position, samples in real_variants.items():
        for sample_dictionary in samples:
            for sample_name, variants in sample_dictionary.items():
                position_ID = position[0]+','+str(position[1])+','+position[2]+','+position[3]
                final_position_error_rate = position_final_error_rate[position_ID]['final position error rate']
                final_delta = position_final_error_rate[position_ID]['final delta']
                total_variants_called = position_final_error_rate[position_ID]['total variants called']
                MLE_outcomes = position_final_error_rate[position_ID]['MLE outcomes']
                MLE_comments = position_final_error_rate[position_ID]['MLE comments']

                if 'iteration called at' in variants.keys():
                    iteration_called_at = variants['iteration called at']
                else:
                    iteration_called_at = ''
                if sample_name in sample_variants.keys():
                    sample_variants[sample_name][position_ID] = {'variant depth': variants['variant_depth'],
                                                                  'total depth': variants['total_depth'],
                                                                  'VAF': variants['variant_depth']/variants['total_depth'],
                                                                     'fitting method': variants['fitting_method'],
                                                                  'iteration called at':  iteration_called_at,
                                                                  'p-value': variants['p_value'], 'call': 'REAL VARIANT',
                                                                'position final error rate': final_position_error_rate,
                                                                'position final delta': final_delta,
                                                                'total variants called at position': total_variants_called,
                                                                'MLE outcomes': MLE_outcomes,
                                                                'MLE comments': MLE_comments}
                else:
                    sample_variants[sample_name] = {}
                    sample_variants[sample_name][position_ID] = {'variant depth': variants['variant_depth'],
                                                                  'total depth': variants['total_depth'],
                                                                  'VAF': variants['variant_depth']/variants['total_depth'],
                                                                     'fitting method': variants['fitting_method'],
                                                                  'iteration called at':  iteration_called_at,
                                                                  'p-value': variants['p_value'], 'call': 'REAL VARIANT',
                                                                'position final error rate': final_position_error_rate,
                                                                'position final delta': final_delta,
                                                                'total variants called at position': total_variants_called,
                                                                'MLE outcomes': MLE_outcomes,
                                                                'MLE comments': MLE_comments}

    # create a dataframe for each sample, which contains all the called real variants
    BB_sample_variants_dataframes = {}

    for sample_name, variants in sample_variants.items():
        variants_df = pd.DataFrame.from_dict(variants, orient = 'index')
        BB_sample_variants_dataframes[sample_name]=variants_df

    #create a sample based dictionary of DCS errors
    errors_dictionary = {}
    sample_errors = {}

    for k, v in position_error_p_values.items():
        sample_dictionary = {}
        position_ID = str(k[0])+','+str(k[1])+','+str(k[2])+','+str(k[3])
        final_position_error_rate = position_final_error_rate[position_ID]['final position error rate']
        final_delta = position_final_error_rate[position_ID]['final delta']
        total_variants_called = position_final_error_rate[position_ID]['total variants called']
        MLE_outcomes = position_final_error_rate[position_ID]['MLE outcomes']
        MLE_comments = position_final_error_rate[position_ID]['MLE comments']
        for iteration, samples in v.items():
            for i in samples:
                sample_ID = i[0]
                variant_depth = i[1]
                total_depth = i[2]
                VAF = i[3]
                fitting_method = i[6]
                iteration_call = i[4]
                p_value = i[5]
                sample_dictionary[sample_ID] = {'variant depth': variant_depth, 'total depth': total_depth, 'VAF': VAF, 'fitting method': fitting_method,
                                               'iteration called at': iteration_call, 'p-value': p_value, 'call': 'ERROR',
                                                'position final error rate': final_position_error_rate,
                                                'position final delta': final_delta, 'total variants called at position': total_variants_called,
                                                'MLE outcomes': MLE_outcomes,'MLE comments': MLE_comments}

        errors_dictionary[position_ID]=sample_dictionary


    for k, v in errors_dictionary.items():
        position_ID = k
        for sample_name, variants in v.items():
            if sample_name in sample_errors.keys():
                sample_errors[sample_name][position_ID] = variants
            else:
                sample_errors[sample_name] = {}
                sample_errors[sample_name][position_ID] = variants

    #exclude variants that were called as real from the errors_dictionary:
    sample_errors_excluding_real_dataframes = {}
    for sample_name, error_variants in sample_errors.items():
        error_variants_dict = {}
        positions_with_real_variants = list(sample_variants[sample_name].keys())
        for position, details in error_variants.items():
            if position not in positions_with_real_variants:
                error_variants_dict[position]=details

        error_variants_df = pd.DataFrame.from_dict(error_variants_dict, orient = 'index')
        sample_errors_excluding_real_dataframes[sample_name]=error_variants_df

    # create a dataframe for each sample, which contains all the error variants and all the real variants
    BB_variants_dataframes = {}
    for sample_name in sample_errors.keys():
        errors_df = sample_errors_excluding_real_dataframes[sample_name]
        variants_df = BB_sample_variants_dataframes[sample_name]
        all_variants = pd.concat([errors_df, variants_df])
        BB_variants_dataframes[sample_name]=all_variants

    #merge the called variants dataframes with the annotated variants file (just keeping the variants that were called)
    BB_sample_variants_dataframes_annotated = {}
    for sample_name in BB_variants_dataframes.keys():
        sample_abbrev = sample_name.split('_')[0]+'_'+sample_name.split('_')[1]

        annotated1 = 'UKCTOCS_sample_level_results_V2/'+sample_abbrev+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_annotated.txt'
        annotated2 = 'UKCTOCS_sample_level_results_V2/'+sample_abbrev+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_SNV_watson_code_'+SSCS_or_DCS+'_variants_MUFs_3_annotated.txt'

        try:
            df = pd.read_csv(annotated1, delimiter = '\t')
        except FileNotFoundError:
            df = pd.read_csv(annotated2, delimiter = '\t')

        df['position ID'] = df[['chromosome', 'start', 'REF', 'ALT']].astype(str).agg(','.join, axis=1)
        df = df.set_index('position ID')
        called_variants = BB_variants_dataframes[sample_name]
        annotated_variants = pd.merge(called_variants, df, how = "left", left_index=True, right_index = True)
        annotated_variants = annotated_variants.drop(columns=['variant depth', 'total depth', 'VAF_x'])
        annotated_variants = annotated_variants.reset_index()
        annotated_variants = annotated_variants.rename(columns={'index': 'position ID', 'VAF_y': 'VAF', 'variant_depth': 'variant depth', 'total_depth': 'total depth'})
        annotated_variants = annotated_variants.dropna(thresh=30) #drop row if 17 or more values are NaN (i.e is outside the panel target regions)
        annotated_variants['sample name']=sample_name #create a column that has the sample name in it

        # move the sample name column to head of list using index, pop and insert
        cols = list(annotated_variants)
        cols.insert(0, cols.pop(cols.index('sample name')))
        annotated_variants = annotated_variants.loc[:, cols]
        annotated_variants = annotated_variants.sort_values(by=['VAF'], ascending=False)

        BB_sample_variants_dataframes_annotated[sample_name]=annotated_variants

        #save a CSV for each sample containing the variant calls
        annotated_variants.to_csv('UKCTOCS_sample_level_results_V2/'+sample_abbrev+'/SNV/'+SSCS_or_DCS+'/'+sample_name+'_SNV_watson_code_'+SSCS_or_DCS+'_beta_binomial_SNV_all_variant_calls.txt', sep = '\t', index = False)

    BB_all_samples_df = pd.concat(BB_sample_variants_dataframes_annotated.values())
    BB_all_samples_df = BB_all_samples_df.sort_values(['sample name', 'VAF'], ascending=[True, False])

    #apply other filters
    BB_all_samples_df = BB_all_samples_df[BB_all_samples_df['strand_bias_fisher_p_value_phred']<=60] #gets rid of ~6
    BB_all_samples_df = BB_all_samples_df[BB_all_samples_df['mean_position_in_read']>=8] #gets rid of ~6
    BB_all_samples_df = BB_all_samples_df[BB_all_samples_df['variant_mean_UMI_family_size']>=6] #doesn't get rid of any

    #save a CSV of all the variant calls across all samples for this library
    BB_all_samples_df.to_csv('Data_files/'+library+'_'+SSCS_or_DCS+'_beta_binomial_all_variant_calls.csv', index = False)

    return print('sample specific variant call csv file saved AND all samples variant calls csv file saved')

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--library", action="store", dest="library", help="name of the sequencing library, e.g. SLX_20124", required=True)
    parser.add_argument("--SSCS_or_DCS", action="store", type=str, dest="SSCS_or_DCS", help="whether SSCS or DCS variant calls", required=True)
    parser.add_argument("--samples_to_exclude", action="store", type=list, dest="samples_to_exclude", help="list of samples to exclude", required=True)

    o = parser.parse_args()

    library = o.library
    SSCS_or_DCS = o.SSCS_or_DCS
    samples_to_exclude = o.samples_to_exclude

    ref_genome = Fasta('Homo_sapiens_assembly19.fasta')

    changes = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'] #group e.g. C>T and G>A togethe
    complementary_changes = {'C>T':'G>A', 'C>A':'G>T', 'C>G':'G>C', 'T>C':'A>G', 'T>A':'A>T', 'T>G':'A>C',
                            'G>A':'C>T','G>T':'C>A', 'G>C':'C>G', 'A>G':'T>C','A>T':'T>A', 'A>C':'T>G'}

    bases = ['A', 'C', 'G', 'T']
    complementary_bases = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}

    #Create a dataframe of all positions in the TWIST SNV panel
    all_positions_df = pd.read_csv('Data_files/TWIST_SNV_panel_all_possible_changes_panel_sites.csv', usecols=[1, 2, 3, 4, 5, 6])

    #create a dictionary of the COSMIC positions
    COSMIC_positions = {}
    cosmic_file = 'Data_files/All_SNV_panel_positions_seen_in_COSMIC_haem_lymphoid.csv'
    with open(cosmic_file) as csvfile:
        row_count=0
        reader = csv.reader(csvfile)
        for row in reader:
            if row_count>0:
                chromosome = row[0]
                position = int(row[1])
                ref = row[3]
                alt = row[4]
                COSMIC = int(row[5])
                COSMIC_positions[(chromosome, position, ref, alt)]=COSMIC
            row_count+=1

    # Step 1 = Create a list of all the cases and controls
    cases_and_controls = create_list_of_cases_and_controls(library, samples_to_exclude)
    #if you want to exclude any samples, include their names (e.g. 'C92_017_s2_SNV') in a list, otherwise put an empty list

    # Step 2 = create a dataframe for the library, containing information on the variant depths and total depths for all samples at each position
    library_df = create_dataframe(all_positions_df, SSCS_or_DCS, cases_and_controls, library, ref_genome)

    # Step 3 = get a list of the multi-timepoint samples on this lane
    multi_timepoint_samples = list_of_samples_with_multiple_timepoints(cases_and_controls)

    # Step 4 = run the beta-binomial model
    high_VAF_threshold = 0.1
    p_value_threshold = 6e-6
    MLEs, real_variants, total_variants_called, iterative_or_not, total_positions_analysed, errors_p_values, real_sample_p_values = noniterative_or_iterative_beta_binomial_binomial_variant_calling_for_multi_timepoints(library_df, high_VAF_threshold, p_value_threshold, COSMIC_positions, multi_timepoint_samples)

    # Step 5 = save output final_epsilon
    position_final_error_rate = {}
    for k, v in MLEs.items():
        position_ID = str(k[0])+','+str(k[1])+','+str(k[2])+','+str(k[3])
        position_final_error_rate[position_ID]=v

    save_all_variant_calls(real_variants, errors_p_values, position_final_error_rate, SSCS_or_DCS, library)

    save_real_sample_p_values(errors_p_values, library, SSCS_or_DCS)

    save_error_sample_p_values(errors_p_values, library, SSCS_or_DCS)

    save_MLE_outcomes(MLEs, library_df, library, SSCS_or_DCS)

    return print('duplex error model variant calling complete on library '+library)

if __name__ == "__main__":
	main()
