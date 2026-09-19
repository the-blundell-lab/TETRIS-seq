#!/usr/bin/env python

'''''
Watson code for grouping and filtering chromosomal rearrangements from output of Watson chromosomal rearrangement caller
Version 1.0 (July 2025)

Input:
    1) CSV file produced by chromosomal rearrangement caller, e.g. sample_name+'_translocations_found_just_those_specifically_targeted_both_sides_panel.csv'

Outputs:
    1) csv file containing chromosomal rearrangements grouped (according to type, notation, left and right gene and breakpoint within 500bp of each other), intra-gene rearrangements removed.
    2) xls file of the above, with highlighting for ease of reading in excel

    output files are saved in the same location as the input csv file

Usage:
Watson_code_filter_and_group_chromosomal_rearrangements.py  --infile CSV file from caller

'''''
version = '1.0'

from argparse import ArgumentParser
import os
import glob
import re
import csv
import gzip
import math
import time
import timeit
import random
import copy
from datetime import date
import numpy as np
import pandas as pd
from scipy.stats import binom
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator, LinearLocator
)
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import pysam
from pyfaidx import Fasta
import networkx as nx
from itertools import combinations

today = date.today()
date_today = today.strftime("%d/%m/%Y")

def extract_position(position_info):
    coordinate = position_info.split(' ')[2]
    coordinate_number = coordinate.strip('><')
    return int(coordinate_number)

def extract_position2(position_info):
    try:
        parts = position_info.split(' ')
        if len(parts) < 3:
            return None
        coordinate = parts[2].strip('><')
        return int(coordinate)
    except (ValueError, AttributeError, IndexError):
        return None

def maximum_of_two_columns(depth1, depth2):
    num1 = pd.to_numeric(depth1, errors='coerce')
    num2 = pd.to_numeric(depth2, errors='coerce')
    maxdepth = np.fmax(num1.fillna(-np.inf), num2.fillna(-np.inf))
    maxdepth = maxdepth.replace(-np.inf, np.nan)
    return maxdepth

def assign_proximity_groups(subdf): #for grouping together translocations whose coordinates are within 500bp of each other
    positions = subdf['position'].values
    subgroup = [0]
    current_group = 0
    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] > 500:
            current_group += 1
        subgroup.append(current_group)
    subdf['subgroup'] = current_group_base + np.array(subgroup)
    return subdf

def extract_chromosome(value):
    parts = value.split(' ')
    return parts[1]

def fusion_genes(left_gene, right_gene):
    return left_gene.astype(str)+'::'+right_gene.astype(str)

def filter_and_group_chromosomal_rearrangements(csv_input):
    df = pd.read_csv(csv_input)
    df['LEFT COORDINATE'] = df["LEFT SIDE (5') OF NON-INVERTED BREAKPOINT"].apply(extract_position)
    df['RIGHT COORDINATE'] = df["RIGHT SIDE (3') OF NON-INVERTED BREAKPOINT"].apply(extract_position)
    df = df[df['LEFT GENE'] != df['RIGHT GENE']] #remove rows where the rearrangement supposedly occurs within the same gene
    df['max_normal_depth'] = maximum_of_two_columns(df['normal_depth_left'], df['normal_depth_right'])
    df['max_VAF'] = maximum_of_two_columns(df['VAF_left'], df['VAF_right'])
    
    # STEP 1: Reset and label rows
    df = df.reset_index(drop=True)  # start fresh
    df['row_id'] = df.index         # assign a unique integer row ID
    
    df['GENE_PAIR'] = fusion_genes(df['LEFT GENE'], df['RIGHT GENE'])
    
    # STEP 2: Sort
    df = df.sort_values(by=['NOTATION', 'GENE_PAIR', 'LEFT GENE', 'RIGHT GENE', 'TYPE', 'LEFT COORDINATE', 'RIGHT COORDINATE'])
    
    # STEP 3: Build groups
    grouped_blocks = []
    group_id_base = 0
    
    for _, group in df.groupby(['NOTATION', 'LEFT GENE', 'RIGHT GENE', 'TYPE']):
        g = group.copy()
        G = nx.Graph()    
        row_ids = g['row_id'].tolist()
        G.add_nodes_from(row_ids)
        g['subgroup'] = np.nan
    
        for i, j in combinations(g.index, 2):
            row_i = g.loc[i]
            row_j = g.loc[j]
            if (
                abs(row_i['LEFT COORDINATE'] - row_j['LEFT COORDINATE']) <= 500
                and abs(row_i['RIGHT COORDINATE'] - row_j['RIGHT COORDINATE']) <= 500
            ):
                G.add_edge(row_i['row_id'], row_j['row_id'])
    
        for group_num, component in enumerate(nx.connected_components(G), start=group_id_base):
            idx = g['row_id'].isin(component)
            g.loc[idx, 'subgroup'] = group_num
    
        group_id_base = int(g['subgroup'].max()) + 1
        grouped_blocks.append(g)
    
    # STEP 4: Combine and annotate group sizes
    df_with_groups = pd.concat(grouped_blocks).reset_index(drop=True)
    df_with_groups['subgroup'] = df_with_groups['subgroup'].astype(int)
    group_sizes = df_with_groups['subgroup'].value_counts()
    df_with_groups['subgroup_size'] = df_with_groups['subgroup'].map(group_sizes)
    
    # STEP 5: Rank gene pairs by frequency
    gene_pair_sizes = df_with_groups['GENE_PAIR'].value_counts()
    gene_pair_rank = {gp: i for i, gp in enumerate(gene_pair_sizes.index, start=1)}
    df_with_groups['SORT_ORDER'] = df_with_groups['GENE_PAIR'].map(gene_pair_rank)
    
    summary_cols = [
        'D_NO_OVERLAP',
        ('D_NPP', 'C_SUPP_NPP'),
        ('D_NPP',),
        ('D_NPP', 'softclip_mapping'),
        ('D_SUPP_PP', 'C_PP'),
        ('D_SUPP_NPP', 'D_NPP'),
        ('concordant_1_end_mapping',),
        ('D_NPP_MPP', 'softclip_mapping'),
        ('C_SUPP_NPP', 'softclip_mapping'),
        ('D_SUPP_NPP', 'D_NPP_MPP'),
        ('C_PP',),
        ('D_SUPP_NPP', 'D_NPP', 'C_SUPP_NPP'),
        ('D_SUPP_NPP', 'C_NPP'),
        ('D_SUPP_NPP', 'C_SUPP_NPP'),
        ('D_SUPP_NPP', 'softclip_mapping'),
        ('D_NPP_MPP', 'C_SUPP_NPP'),
        ('C_PP', 'softclip_mapping'),
        ('C_NPP', 'softclip_mapping'),
        ('D_NPP_LEFT_INV', 'C_SUPP_NPP_LEFT'),
        ('D_NPP_LEFT_INV',),
        ('D_NPP_LEFT_INV', 'inv_softclip_mapping'),
        ('D_SUPP_NPP_LEFT', 'D_NPP_LEFT_INV'),
        ('D_SUPP_PP_INV', 'C_PP'),
        'D_NO_OVERLAP_INV',
        ('D_SUPP_NPP_INV', 'D_NPP'),
        ('D_SUPP_NPP_INV', 'D_NPP_LEFT_INV'),
        ('D_SUPP_NPP_INV', 'inv_softclip_mapping'),
        ('C_PP', 'inv_softclip_mapping'),
        ('D_SUPP_PP_INV', 'inv_softclip_mapping'),
        ('D_NPP', 'inv_softclip_mapping'),
        ('D_SUPP_NPP_LEFT', 'inv_softclip_mapping'),
        ('D_NPP_RIGHT_INV', 'C_SUPP_NPP_RIGHT'),
        ('C_SUPP_NPP_RIGHT', 'inv_softclip_mapping'),
        ('D_NPP_RIGHT_INV',),
        ('D_SUPP_NPP_RIGHT', 'D_NPP_RIGHT_INV'),
        ('D_NPP_RIGHT_INV', 'inv_softclip_mapping'),
        ('D_SUPP_NPP_INV', 'D_NPP_RIGHT_INV', 'C_SUPP_NPP_RIGHT'),
        ('D_SUPP_NPP_INV', 'D_NPP_RIGHT_INV'),
        ('C_SUPP_NPP_RIGHT',),
        ('D_SUPP_NPP_RIGHT', 'inv_softclip_mapping'),
        ('translocation_depth'),
        ('normal_depth_left'),
        ('normal_depth_right'),
        ('VAF_left'),
        ('VAF_right'),
        ('max_VAF')
    ]
    
    # STEP 6: Final sort for grouping
    df_sorted = df_with_groups.sort_values(
        by=['SORT_ORDER', 'GENE_PAIR', 'subgroup_size', 'subgroup', 'LEFT COORDINATE'],
        ascending=[True, True, False, True, True]
    ).reset_index(drop=True)
    
    # STEP 7: Build summary and output
    summary_blocks = []
    ungrouped_blocks = []
    
    # Define summary columns to sum
    summary_cols = [col for col in df.columns if col not in [
        'GROUP TAG', 'row_id', 'subgroup', 'subgroup_size',
        'LEFT COORDINATE', 'RIGHT COORDINATE', 'FIRST CHROMOSOME', 'SECOND CHROMOSOME',
        'LEFT GENE', 'RIGHT GENE', 'GENE_PAIR', 'TYPE', 'NOTATION'
    ]]
    
    for subgroup_id, group in df_sorted.groupby('subgroup', sort=False):
        # group = group.copy()
        group = group.sort_values(by='translocation_depth', ascending=False, na_position='last').copy()
        group_size = len(group)
        # summary_label = f"{group.iloc[0]['NOTATION']} {group.iloc[0]['LEFT GENE']}::{group.iloc[0]['RIGHT GENE']} {group.iloc[0]['TYPE']}"
        summary_label = f"GROUP: {group.iloc[0]['NOTATION']} {group.iloc[0]['LEFT GENE']}::{group.iloc[0]['RIGHT GENE']}"
    
        left_chr = extract_chromosome(group.iloc[0]["LEFT SIDE (5') OF NON-INVERTED BREAKPOINT"])
        right_chr = extract_chromosome(group.iloc[0]["RIGHT SIDE (3') OF NON-INVERTED BREAKPOINT"])
        
        left_range = f"chr{left_chr} {int(group['LEFT COORDINATE'].min())}-{int(group['LEFT COORDINATE'].max())}"
        right_range = f"chr{right_chr} {int(group['RIGHT COORDINATE'].min())}-{int(group['RIGHT COORDINATE'].max())}"
    
        if group_size > 1:
            summary_row = {
                'GROUP TAG': summary_label,
                'NOTATION': '',
                'LEFT GENE': '',
                'RIGHT GENE': '',
                'TYPE': '',
                "LEFT SIDE (5') OF NON-INVERTED BREAKPOINT": left_range,
                "RIGHT SIDE (3') OF NON-INVERTED BREAKPOINT": right_range,
                'row_id': '',
                'subgroup': subgroup_id,
                'subgroup_size': group_size,
            }
            
            for col in summary_cols:
                if col in group.columns:
                    if col not in ["LEFT SIDE (5') OF NON-INVERTED BREAKPOINT", "RIGHT SIDE (3') OF NON-INVERTED BREAKPOINT"]:
                        summary_row[col] = pd.to_numeric(group[col], errors='coerce').sum()
    
            group['GROUP TAG'] = ''
            summary_blocks.append(pd.DataFrame([summary_row]))
            summary_blocks.append(group)
    
            blank_row = pd.DataFrame([{col: '' for col in group.columns}])
            summary_blocks.append(blank_row)
    
        else:
            group['GROUP TAG'] = ''
            ungrouped_blocks.append(group)
    
    # print(summary_blocks)
    
    # Add ungrouped rows if present
    if ungrouped_blocks:
        ungrouped = pd.concat(ungrouped_blocks, ignore_index=True)
        ungrouped = ungrouped.sort_values(by='translocation_depth', ascending=False, na_position='last')
    
        ungrouped_summary = pd.DataFrame([{
            'GROUP TAG': 'UNGROUPED',
            **{col: '' for col in ungrouped.columns if col != 'GROUP TAG'}
        }])
    
        summary_blocks.append(ungrouped_summary)
        summary_blocks.append(ungrouped)
    
    # STEP 8: Finalize output
    final_df = pd.concat(summary_blocks, ignore_index=True)
    
    # Make GROUP TAG the first column
    cols = ['GROUP TAG'] + [col for col in final_df.columns if col != 'GROUP TAG']
    final_df = final_df[cols]
    
    # Column order preference
    fixed_columns = ['GROUP TAG', 'NOTATION', "LEFT SIDE (5') OF NON-INVERTED BREAKPOINT", "RIGHT SIDE (3') OF NON-INVERTED BREAKPOINT",
                     'TYPE', 'translocation_depth', 'normal_depth_left', 'normal_depth_right', 'max_normal_depth',
                     'VAF_left', 'VAF_right', 'max_VAF']
    other_columns = [col for col in final_df.columns if col not in fixed_columns]
    new_order = fixed_columns + other_columns
    final_df = final_df[new_order]
    final_df = final_df.drop(columns=['row_id', 'subgroup', 'SORT_ORDER', 'subgroup_size'], errors='ignore')
    
    output_csv_file = str(csv_input).replace('.csv', '_grouped_and_filtered.csv')
    output_xls_file = str(csv_input).replace('.csv', '_grouped_and_filtered.xlsx')

    final_df.to_csv(output_csv_file, index=False)
    
    # Save final_df to Excel with formatting
    output_file = output_xls_file
    
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        final_df.to_excel(writer, index=False, sheet_name='Summary')
    
        workbook  = writer.book
        worksheet = writer.sheets['Summary']
    
        # Freeze the top row
        worksheet.freeze_panes(1, 0)
    
        # Autosize all columns
        for i, col in enumerate(final_df.columns):
            # Find length of longest entry in the column (including header)
            series = final_df[col].astype(str)
            max_len = max(series.map(len).max(), len(str(col)))  # +2 for padding
            worksheet.set_column(i, i, min(max_len, 35))
    
        # Narrower widths for P to BC
        for col_idx in range(15, 55):
            worksheet.set_column(col_idx, col_idx, 8)
    
        # Define a bold and shaded format
        summary_format = workbook.add_format({
            'bold': True,
            'bg_color': '#D9E1F2',  # light blue shade
            'border': 0
        })
    
        # Loop through rows and apply format to summary rows
        for idx, value in enumerate(final_df['GROUP TAG']):
            if pd.notna(value) and value != '':
                worksheet.set_row(idx + 1, None, summary_format)  # +1 to account for header row
    
    return final_df

def main():
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="chromosomal rearrangement caller csv file", required=True)
    o = parser.parse_args()

    csv_input = o.infile

    filter_and_group_chromosomal_rearrangements(csv_input)

    return print('Grouped and filtered chromosomal rearrangement file produced')

if __name__ == "__main__":
	main()
