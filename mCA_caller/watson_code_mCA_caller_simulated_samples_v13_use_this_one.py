#!/usr/bin/env python3
"""
Run mCA Caller on Simulated Samples and Compare to Ground Truth

Processes all simulated samples, detects mCAs, and generates:
1. all_detected_mCAs.csv - All calls made by the algorithm
2. truth_comparison.csv - Ground truth vs detected calls with matching
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
from matplotlib import gridspec
import ruptures as rpt

import seaborn as sns
import warnings
import time
from scipy import stats
from scipy.stats import binomtest
warnings.filterwarnings('ignore')

# ==========================================
# CONFIGURATION
# ==========================================
SIMULATED_DIR = None  # set in main
MANIFEST_FILE = None
OUTPUT_DIR = None
PLOTS_DIR = None
PLOTS_UNDETECTED_DIR = None

ideogram_file = 'chromosome_ideogram_hg19.txt'

#mCA calling variables (based on allelic imbalance of heterozygous SNPs)
het_lo = 0.02 
het_hi = 0.98

#larger windows reduce standard error, allowing you to trust smaller deviations
window_hets=15 #hw many heterozygous SNPs go into each sliding window
step_hets=3 #how far you slide the window each time (in het SNPs)

loss_lrr = -0.12 #only call a loss if LRR <= this number (i.e. sensitivity cut-off) (0 = diploid, -0.58 = 100% loss, 0.12 = -15-20% cell fraction)
gain_lrr = 0.12 #only call a gain if LRR >= this number (i.e. sensitivity cut-off) (0 = diploit, +0.58 = 100% gain, +0.12 = ~20% cell fraction)
neutral_lrr = 0.04 #how close to zero the LRR must be to call something "copy neutral"

min_supporting_fraction = 0.25 #proportion of het SNPs in a window that need to support the mCA to call it
min_consecutive_hets = 5 #require allelic imbalance in at least 5 consecutive SNPs to call an mCA
min_spanning_hets = 3.0 #distance between 1st and last supporting het for the mCA needs to be at least this distance (Mb)
min_span_mb = 3 #minimum length (Mb) of CN-LOH you're willing to call (i.e. ignore anything below this length)
gap_mb=8 #maximum allowed distance (Mb) between consecutive supporting het SNPs before you stop (i.e. end at the last SNP before the gap)

#allelic imbalance thresholds (most important factor affecting sensitivity)
enter_thresh = 0.03 #the allelic imbalance (AI) threshold needed to start a BAF allelic imbalance region (i.e. enter when window mean AI >= enter_thresh) (higher threshold = more specific, less sensitive)
exit_thresh = 0.03 #the AI threshold needed to end a BAF allelic imbalance region (i.e. exit from window mean AI < exit_thresh) (if exit too high, call ends early; if exit is too low, call extends too far)
ai_point_thresh=0.03 #SNP level AI threshold used to refine the end boundary
support_ai_thresh = 0.03

multiscale_configs = [
    # (window_hets, step_hets, min_span_mb, min_consecutive_hets, min_spanning_hets)
    (40,  8,  8,  3,  5.0),   # Medium: catches 10-20 Mb events at 15-25% CF
    (80,  15, 15, 3,  10.0),  # Large:  catches >20 Mb / whole-arm at 10-20% CF
]


# Matching criteria - to decide if true/ false positive
LENGTH_TOLERANCE = 0.30  # ±30%
POSITION_TOLERANCE_MB = 5.0  # ±5 Mb for start/end
OVERLAP_THRESHOLD = 0.5  # 50% reciprocal overlap

# Centromere dictionary
centromere_dict = {
    'chr1': {'start': 121500000, 'end': 128900000}, 'chr2': {'start': 91800000, 'end': 96000000},
    'chr3': {'start': 87800000, 'end': 93900000}, 'chr4': {'start': 48200000, 'end': 52700000},
    'chr5': {'start': 46100000, 'end': 51400000}, 'chr6': {'start': 58500000, 'end': 63300000},
    'chr7': {'start': 58100000, 'end': 62100000}, 'chr8': {'start': 43200000, 'end': 47200000},
    'chr9': {'start': 42200000, 'end': 45500000}, 'chr10': {'start': 38000000, 'end': 42300000},
    'chr11': {'start': 51000000, 'end': 55800000}, 'chr12': {'start': 33200000, 'end': 37800000},
    'chr13': {'start': 16300000, 'end': 19500000}, 'chr14': {'start': 15600000, 'end': 19100000},
    'chr15': {'start': 17500000, 'end': 20700000}, 'chr16': {'start': 34600000, 'end': 38600000},
    'chr17': {'start': 22700000, 'end': 27400000}, 'chr18': {'start': 15400000, 'end': 19000000},
    'chr19': {'start': 24200000, 'end': 28100000}, 'chr20': {'start': 25700000, 'end': 30400000},
    'chr21': {'start': 10900000, 'end': 13000000}, 'chr22': {'start': 13700000, 'end': 17400000},
    'chrX': {'start': 58100000, 'end': 63800000}
}

chromosome_sizes = {
    'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
    'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
    'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
    'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560
}

all_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX']

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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

# ==========================================
# FUNCTIONS FOR CALLING MCAS
# ==========================================

def ideograms(ideogram_file, chromosome):
    
    color_lookup = {'gneg': (1., 1., 1.),
                    'gpos25': (.6, .6, .6),
                    'gpos50': (.4, .4, .4),
                    'gpos75': (.2, .2, .2),
                   'gpos100': (0., 0., 0.),
                      'acen': (.8, .4, .4),
                      'gvar': (.8, .8, .8),
                     'stalk': (.9, .9, .9)}
    
    ideogram = open(ideogram_file)
    ideogram.readline()
    xranges = []
    colors = []
    mid_points = []
    labels = []

    for line in ideogram:
        chrom, start, stop, label, stain = line.strip().split('\t')
        start = int(start)
        stop = int(stop)
        width = stop - start
        mid_point = start + (width/2)
        if chrom == chromosome:
            xranges.append((start, width))
            colors.append(color_lookup[stain])
            mid_points.append(mid_point)
            labels.append(label)
        
    return xranges, [0, 0.9], colors, mid_points, labels

def plot_chromosome(ideogram_file, chromosome, ax):

    xranges, yrange, colors, midpoints, labels = ideograms(ideogram_file, chromosome)

    ax.broken_barh(xranges, yrange, facecolors= colors, edgecolor = 'black')

    ax.set_xticks(midpoints)
    ax.set_xticklabels(labels, rotation = 90, fontsize = 9)
    ax.set_yticks([])
    ax.text(-0.013, 0.35, chromosome, transform=ax.transAxes, fontsize = 15, ha = 'right')
    ax.xaxis.set_tick_params(width=0.8, color = grey3, length = 6)

    ax.minorticks_off()

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    return ax

def load_centromeres(filepath):
    """
    Parses the UCSC cytoband file to find the start/end of centromeres ('acen').
    Returns a dictionary: {'1': (start, end), '2': (start, end), ...}
    """
    import pandas as pd
    
    df = pd.read_csv(filepath, sep="\t", comment='#', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'type'])
    
    # Filter for centromeric regions ('acen')
    acen = df[df['type'] == 'acen'].copy()
    
    centromeres = {}
    for chrom, grp in acen.groupby('chrom'):
        # Centromeres usually span two bands (p-arm end, q-arm start)
        # We take the overall min start and max end.
        start = grp['start'].min()
        end = grp['end'].max()
        
        centromeres[chrom] = {'start': start, 'end': end}
        
    return centromeres

def load_and_merge_data(sample_name, library):
    """
    Fast vectorized loading. Replaces create_joint_BAF_LRR_file.
    Returns a merged DataFrame directly, skipping the intermediate text file.
    """
    # print(f"Loading data for {sample_name}...")
    
    # 1. Load LRR File
    lrr_path = library+'/'+sample_name+'/PON_normalised_log2ratios_Feb2026/'+sample_name+'_PON_normalised_read_depths_and_LRR.txt'

    # Skip header rows (assuming 5 lines based on your row_count>5 logic)
    df_lrr = pd.read_csv(lrr_path, sep='\t', skiprows=5, header=None, 
                         names=['chromosome', 'start', 'stop', 'band', 'mean_region_depth', 'normalised_region_depth', 'log2ratio', 'coefficient_of_variation_PON', 'p-value', 'mean_normalised_depth_PON', 'sample_depth_normalised_by_PON'])
    
    # Keep only necessary columns and ensure numeric
    df_lrr = df_lrr[['chromosome', 'start', 'stop', 'log2ratio', 'p-value', 'coefficient_of_variation_PON']].copy()
    df_lrr['start'] = pd.to_numeric(df_lrr['start'], errors='coerce')
    df_lrr['stop'] = pd.to_numeric(df_lrr['stop'], errors='coerce')
    df_lrr = df_lrr.sort_values(['chromosome', 'start']) # Essential for merge_asof

    # 2. Load SNP File
    snp_path = library+'/'+sample_name+'/'+sample_name+'_watson_code_SSCS_variant_calling_only_SNPs_annovar_annotated.txt'
    if not os.path.exists(snp_path):
        snp_path = library+'/'+sample_name+'/'+sample_name+'_CNV_watson_code_SSCS_variant_calling_only_SNPs_annovar_annotated.txt'

    df_snp = pd.read_csv(snp_path, sep='\t')
    
    # Rename columns to match desired output
    df_snp = df_snp.rename(columns={
        'REF': 'ref', 'ALT': 'alt', 'VAF': 'BAF',
    })
    
    # Depending on your exact file header, you might need to adjust column names above.
    # Assuming standard columns: chromosome, position, ref, alt, RSID, BAF...
    
    # Standardize Chromosomes (remove 'chr' if present to match)
    df_lrr['chromosome'] = df_lrr['chromosome'].astype(str)
    df_snp['chromosome'] = df_snp['chromosome'].astype(str)
    df_snp['position'] = pd.to_numeric(df_snp['position'], errors='coerce')

    # 3. FAST MERGE (The Magic Step)
    # We split by chromosome first to make merge_asof happy (it requires sorted unique keys usually)
    merged_dfs = []
    
    for chrom in df_snp['chromosome'].unique():
        # Get subset for this chrom
        snps_c = df_snp[df_snp['chromosome'] == chrom].sort_values('position')
        lrr_c = df_lrr[df_lrr['chromosome'] == chrom]
        
        if lrr_c.empty or snps_c.empty:
            continue

        # 1. Ensure SNP positions are integers
        snps_c['position'] = snps_c['position'].astype('int64')

        # 2. Ensure LRR start/stop positions are integers
        # (Floats often happen if there was a NaN somewhere or if pandas inferred it loosely)
        lrr_c = lrr_c.dropna(subset=['start', 'stop']) # Safety: remove rows with no coordinates
        lrr_c['start'] = lrr_c['start'].astype('int64')
        lrr_c['stop']  = lrr_c['stop'].astype('int64')
            
        # merge_asof looks for the nearest match on 'start'
        # direction='backward' means we find the LRR probe where LRR_start <= SNP_pos
        m = pd.merge_asof(snps_c, lrr_c, left_on='position', right_on='start', 
                          direction='backward', suffixes=('', '_lrr'))
        
        # Filter: The match is only valid if SNP_pos <= LRR_stop
        m = m[m['position'] <= m['stop']]
        
        # Rename for compatibility with your cleaning function
        m = m.rename(columns={'coefficient_of_variation_PON': 'probe_CV', 'p-value': 'LRR_p_value', 'start': 'probe_start', 'stop': 'probe_stop', 'log2ratio': 'LRR'})
        
        # Calculate BAF deviation here to save time later
        m['BAF'] = pd.to_numeric(m['BAF'], errors='coerce')
        m['BAF_deviation'] = abs(0.5 - m['BAF'])
        
        merged_dfs.append(m)

    if not merged_dfs:
        return pd.DataFrame()
        
    final_df = pd.concat(merged_dfs, ignore_index=True)
    # print(f"Data loaded. Rows: {len(final_df)}")
    return final_df

def clean_dataframe(df, het_lo, het_hi): 
    """ Collapse multiple rows at same (chrom,pos) to a single representative. 
    Preference: pick BAF closest to 0.5 (most informative for allelic imbalance). 
    Isolates the informative SNPs for imbalance (based on MAF and allelic imbalance (AI)). 
    Prevents noisy depth probes from forcing false CNV calls by blanking LRR where probe_CV is too high. """ 
    
    df = df.copy() 
    df["BAF"] = pd.to_numeric(df["BAF"], errors="coerce") 
    df["pos"] = pd.to_numeric(df["position"], errors="coerce") 
    df["chrom"] = df["chromosome"].astype(str) 
    
    # score: distance to 0.5; smaller is better 
    df["dist05"] = (df["BAF"] - 0.5).abs() 
    
    # prefer het-like (in case there is both 0 and 0.49 etc) 
    df["is_het_like"] = df["BAF"].between(het_lo, het_hi) #flags which SNPs are informative for imbalance/LOH (i.e. if BAF is 0 o 1 then it is homozygous and carries no allelic-balance information) 
    
    # sort so het-like + closest to 0.5 comes first 
    df = df.sort_values(["chrom", "pos", "is_het_like", "dist05"], ascending=[True, True, False, True]) 

    out = df.drop_duplicates(["chrom", "pos"]).copy() #because we sorted, the first row kept per (chrom, pos) is the best candidate 
    out = out.drop(columns=["dist05"]) #remove the helper columns 
    
    # tidy column types 
    out["chrom"] = out["chromosome"].astype(str) 
    out["pos"] = pd.to_numeric(out["position"], errors="coerce") 
    out["BAF"] = pd.to_numeric(out["BAF"], errors="coerce") 
    out["LRR"] = pd.to_numeric(out["LRR"], errors="coerce") 
    out["probe_CV"] = pd.to_numeric(out["probe_CV"], errors="coerce") 
    out["LRR_p_value"] = pd.to_numeric(out["LRR_p_value"], errors="coerce") 
    
    #drop rows that are missing essential fields, then sort by genomic order 
    out = out.dropna(subset=["chrom", "pos", "BAF"]).sort_values(["chrom", "pos"]).reset_index(drop=True) 
    
    #calculate the MAF (minor allele fraction) at each locus, e.g. if BAF = 0.48, MAF = min(0.48, 0.52) = 0.48, if BAF = 0.8 (MAF = min(0.8, 0.20) = 0.2 
    #MAF is always between 0 and 0.5 (if normal heterozygous SNP: MAF = 0.5). For homozygotes (BAF ~0 or 1, set MAF to NaN) 
    out["maf"] = np.where(out["is_het_like"], np.minimum(out["BAF"], 1 - out["BAF"]), np.nan) 
    
    #calculate allelic imbalance magnitude, e.g. if MAF = 0.5, AI = 0 (balanced), if MAF = 0.25, AI = 0.25 (strong imbalance) 
    out["ai"] = 0.5 - out["maf"] # only defined for het-like 
    
    # LRR quality gate: if probe_CV too high, blank LRR (but still keep SNP for BAF track) 
    # Robust CV cutoff based on distribution
    cv_vals = out["probe_CV"].to_numpy(dtype=float)
    cv_vals = cv_vals[np.isfinite(cv_vals)]

    if cv_vals.size == 0:
        cv_med = cv_mad = cv_cut = np.nan
        out["lrr_filt"] = np.nan
    else:
        cv_med = np.median(cv_vals)
        cv_mad = np.median(np.abs(cv_vals - cv_med))
        cv_cut = cv_med + 3 * cv_mad
        out["lrr_filt"] = out["LRR"].where(out["probe_CV"] <= cv_cut, np.nan)

    # print("LRR kept fraction:", np.isfinite(out["lrr_filt"]).mean())
    # print("CV median:", cv_med, "CV cutoff:", cv_cut)
    
    return out

def baf_peaks_mu1_mu2(seg, min_per_side, het_lo=0.15, het_hi=0.85):

    b = pd.to_numeric(seg["BAF"], errors="coerce")
    b = b[b.between(het_lo, het_hi)].dropna()
    if b.empty:
        return np.nan, np.nan, np.nan
    lower = b[b < 0.5]
    upper = b[b > 0.5]
    
    # Standard two-sided estimation
    if (len(lower) >= min_per_side) and (len(upper) >= min_per_side):
        mu1 = float(np.median(lower))
        mu2 = float(np.median(upper))
        delta = mu2 - mu1
        return mu1, mu2, delta
    
    # Asymmetric: use dominant side and mirror around 0.5
    min_mirrored_delta = 0.20
    if len(lower) >= min_per_side and len(upper) < min_per_side:
        mu1 = float(np.median(lower))
        delta = 2 * (0.5 - mu1)
        mu2 = mu1 + delta
        if delta < min_mirrored_delta:
            return np.nan, np.nan, np.nan
        return mu1, mu2, delta
    
    if len(upper) >= min_per_side and len(lower) < min_per_side:
        mu2 = float(np.median(upper))
        delta = 2 * (mu2 - 0.5)
        mu1 = mu2 - delta
        if delta < min_mirrored_delta:
            return np.nan, np.nan, np.nan
        return mu1, mu2, delta
    
    return np.nan, np.nan, np.nan

def p_from_delta(delta, event):
    if np.isnan(delta):
        return np.nan
    delta = np.clip(delta, 0, 0.999999)

    if event == "CN-LOH":
        p = delta
    elif event == "GAIN":
        p = (2 * delta) / (1 - delta)
    elif event == "LOSS":
        p = (2 * delta) / (1 + delta)
    else:
        return np.nan

    return float(np.clip(p, 0, 1))

def expected_lrr_from_p(p, event):
    if np.isnan(p): 
        return np.nan
    p = float(np.clip(p, 0, 1))
    if event == "GAIN":
        return float(np.log2((2 + p) / 2))
    if event == "LOSS":
        return float(np.log2((2 - p) / 2))
    if event == "CN-LOH":
        return 0.0
    return np.nan

def p_from_lrr(mean_lrr, event):
    if np.isnan(mean_lrr):
        return np.nan
    l = float(mean_lrr)
    if event == "GAIN":
        return float(np.clip(2 * (2**l - 1), 0, 1))
    if event == "LOSS":
        return float(np.clip(2 * (1 - 2**l), 0, 1))
    return np.nan  # CN-LOH not identifiable from LRR alone

def classify_event_by_baf_lrr(delta, mean_lrr, lrr_sigma):
    """
    Choose event by how well the observed median LRR matches the LRR expected
    from the BAF-derived p under each hypothesis.
    lrr_sigma is your expected SD of segment-median LRR noise (tune ~0.05–0.12).
    """
    # If delta is NaN, fall back to LRR-only
    if np.isnan(delta):
        if np.isnan(mean_lrr):
            return {"event": "AI_event_unclear", "p_baf": np.nan, 
                    "lrr_expected_from_p_baf": np.nan, "p_lrr": np.nan,
                    "lrr_fit_score": np.inf, "flag_lrr_baf_mismatch": False}
        
        if abs(mean_lrr) <= 0.04:
            ev = "CN-LOH"
        elif mean_lrr >= 0.12:
            ev = "GAIN"
        elif mean_lrr <= -0.12:
            ev = "LOSS"
        else:
            ev = "AI_event_unclear"
        
        p_lrr = p_from_lrr(mean_lrr, ev)
        return {"event": ev, "p_baf": np.nan,
                "lrr_expected_from_p_baf": np.nan, "p_lrr": p_lrr,
                "lrr_fit_score": np.inf, "flag_lrr_baf_mismatch": False}
    
    # Original logic when delta is valid
    candidates = ["CN-LOH", "LOSS", "GAIN"]
    rows = []
    for ev in candidates:
        p_baf = p_from_delta(delta, ev)
        lrr_exp = expected_lrr_from_p(p_baf, ev)
        # score: squared z residual in LRR space
        if np.isnan(mean_lrr) or np.isnan(lrr_exp):
            score = np.inf
        else:
            score = ((mean_lrr - lrr_exp) / lrr_sigma) ** 2

        # also compute p from LRR for gain/loss to report discordance
        p_lrr = p_from_lrr(mean_lrr, ev)
        rows.append((ev, p_baf, lrr_exp, p_lrr, score))

    # pick smallest score
    best = min(rows, key=lambda x: x[-1])
    best_event, p_baf, lrr_exp, p_lrr, score = best

    # optional: a simple mismatch flag when both p estimates exist
    flag_mismatch = False
    if not np.isnan(p_lrr) and not np.isnan(p_baf):
        if abs(p_lrr - p_baf) > 0.25:
            flag_mismatch = True

    return {
        "event": best_event,
        "p_baf": p_baf,
        "lrr_expected_from_p_baf": lrr_exp,
        "p_lrr": p_lrr,
        "lrr_fit_score": float(score),
        "flag_lrr_baf_mismatch": flag_mismatch,
    }

def mad_sigma(x, min_points=50): #estimate standard deviation of noise using median absolute deviation
    """
    Robust sigma estimate using MAD. Returns NaN if too few points.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    
    if x.size < min_points:
        return np.nan
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return 1.4826 * mad  # MAD -> sigma for approx-Gaussian noise

def lrr_sigma_outside_intervals(d_chrom, bp_intervals, min_points=50):
    """
    Estimate per-chromosome LRR noise from markers OUTSIDE AI-called intervals.
    d_chrom: per-chrom df (sorted), must have columns: position, lrr_filt
    bp_intervals: list of (start_pos, end_pos) in bp
    """
    pos = d_chrom["position"].to_numpy(dtype=int)
    y   = d_chrom["lrr_filt"].to_numpy(dtype=float)

    mask = np.isfinite(y)
    for s, e in bp_intervals:
        mask &= ~((pos >= int(s)) & (pos <= int(e)))

    return mad_sigma(y[mask], min_points=min_points)

def max_consecutive_true(series):
    """Returns the maximum length of consecutive True values in a boolean series."""
    # Convert series to int (0/1), find changes, group by changes, count size
    return series.groupby((series != series.shift()).cumsum()).transform('size').where(series).max()

def add_mca_size_arrows(ax, calls_df, y=1.04, text_y=1.06, fontsize=9, color=grey4):
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

    for _, r in calls_df.iterrows():
        x0, x1 = float(r["start_pos"]), float(r["end_pos"])
        size_mb = (x1 - x0) / 1e6
        event = r["event"]

        ax.annotate("", xy=(x1, y), xytext=(x0, y), xycoords=trans, textcoords=trans, arrowprops=dict(arrowstyle="<->", lw=1.2, color=color),annotation_clip=False)
        ax.text((x0 + x1) / 2, text_y, f"{event}: {size_mb:.1f} Mb", transform=trans, ha="center", va="bottom", fontsize=fontsize, color='black', clip_on=False)

def plot_BAF_AI_windows(df, chrom, calls_df, wdf, enter_thresh, chromosome_sizes):
    event_colors = {
        "GAIN": red4,
        "LOSS": blue4,
        "CN-LOH": orange3,
        "AI_event_unclear": grey4
    }

    # FIX: Handle ALL None inputs
    if calls_df is None:
        calls_df = pd.DataFrame()
    
    if wdf is None:
        wdf = pd.DataFrame()
    
    if df is None:  # Just in case
        raise ValueError("df cannot be None")

    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)
    het = d[d["is_het_like"] & d["ai"].notna()][["position", "ai"]]
    
    # 5 panels: LRR, All BAFs, AI (het-only), Window mean AI, Ideogram
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(
        5, 1,
        figsize=(14, 8),
        sharex=False,
        gridspec_kw={"height_ratios": [8, 8, 8, 4, 1.5], "hspace": 0.08},
        constrained_layout=True
    )
    ax2.sharex(ax1)
    ax3.sharex(ax1)
    ax4.sharex(ax1)
    
    # --- PANEL 1: LRR ---
    ax1.scatter(d["position"], d["lrr_filt"], s=8, color=grey4)
    ax1.axhline(0, linestyle="--", color=grey2)
    ax1.set_ylabel("LRR (filtered)")
    
    # --- PANEL 2: ALL BAFs (NEW) ---
    # Separate het and homozygous SNPs
    het_data = d[d["is_het_like"]].copy()
    hom_data = d[~d["is_het_like"]].copy()
    
    # Plot homozygous SNPs first (background, light gray)
    if not hom_data.empty:
        ax2.scatter(hom_data["position"], hom_data["BAF"], s=10, color='lightgray', alpha=0.8, label='Homozygous SNPs', zorder = 10)
    
    # Plot heterozygous SNPs (gray background)
    ax2.scatter(het_data["position"], het_data["BAF"], s=10, color=grey4, alpha=1.0, zorder = 10)
    
    ax2.axhline(0.5, linestyle="--", color=grey2, zorder = 0)
    ax2.set_ylabel("BAF (all SNPs)")
    ax2.set_ylim(-0.05, 1.05)
    
    # --- PANEL 3: AI (het-only) ---
    ax3.scatter(het["position"], het["ai"], s=10, color=grey4)
    ax3.axhline(enter_thresh, linestyle="--", color=grey2)
    ax3.set_ylabel("AI (het-only)")
    
    # --- Overlay colored points inside called mCAs ---
    if len(calls_df) > 0:
        for _, r in calls_df.iterrows():
            col = event_colors.get(r["event"], "grey")
            mask = (d["position"] >= r["start_pos"]) & (d["position"] <= r["end_pos"])
            het_mask = mask & d["is_het_like"] & d["ai"].notna()
            
            # Panel 1: LRR points
            ax1.scatter(d.loc[mask, "position"], d.loc[mask, "lrr_filt"],
                        s=10, color=col, zorder=100)
            
            # Panel 2: All BAF points (het only - hom stays gray)
            ax2.scatter(d.loc[het_mask, "position"], d.loc[het_mask, "BAF"],
                        s=12, color=col, zorder=100)
            
            # Panel 3: AI points
            ax3.scatter(d.loc[het_mask, "position"], d.loc[het_mask, "ai"],
                        s=14, color=col, zorder=100)
    
    # --- PANEL 4: Window mean AI track ---
    if not wdf.empty:
        ax4.plot((wdf["w_start_pos"] + wdf["w_end_pos"]) / 2, wdf["w_mean_AI"], color=grey4)
        ax4.axhline(enter_thresh, linestyle="--", color=grey2)
        ax4.set_ylabel("Window \nmean AI")
    
    # --- PANEL 5: Ideogram ---
    plot_chromosome(ideogram_file, chrom, ax5)
    
    # Only show required axis lines
    for ax in (ax1, ax2, ax3, ax4, ax5):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.5)
        ax.spines["bottom"].set_linewidth(1.5)
        ax.spines["left"].set_color(grey3)
        ax.spines["bottom"].set_color(grey3)
    
    # X-limits
    chromosome_size = chromosome_sizes[chrom]
    chrom_size_mb = chromosome_size / 1e6
    ax5.text(1.01, 0.5, f"{chrom_size_mb:.0f} Mb", transform=ax5.transAxes, 
             ha="left", va="center", fontsize=12, color='black')
    
    for ax in (ax1, ax2, ax3, ax4, ax5):
        ax.set_xlim(0, chromosome_size)
    
    add_mca_size_arrows(ax1, calls_df)
    
    # Hide x tick labels on upper panels
    for ax in (ax1, ax2, ax3, ax4):
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False, 
                       top=False, labeltop=False)
    
    # Shade merged calls by event type
    if len(calls_df) > 0:
        for _, r in calls_df.iterrows():
            col = event_colors.get(r["event"], "grey")
            for ax in [ax1, ax2, ax3, ax4]:
                ax.axvspan(r["start_pos"], r["end_pos"], color=col, alpha=0.15)
    
    ax5.tick_params(axis="x", bottom=True, labelbottom=True, top=False, labeltop=False)
    
    plt.suptitle(f"{chrom}: BAF_AI windows + merged calls")
    
    return fig

def trim_event_boundaries(het_seg, window_size=5, threshold=0.045):
    """
    Two-Stage Refinement with Adaptive Logic:
    1. THE CHAINSAW: Checks for massive genomic gaps (Centromeres).
    2. THE SCALPEL: Trims edges, but adapts based on signal characteristics.
    """
    if len(het_seg) < window_size:
        return het_seg.iloc[0]['position'], het_seg.iloc[-1]['position']

    vals = het_seg['ai'].to_numpy()
    original_start = het_seg.iloc[0]['position']
    original_end = het_seg.iloc[-1]['position']
    span_mb = (original_end - original_start) / 1e6
    
    # --- ADAPTIVE THRESHOLD: Scale based on probe density ---
    n_hets = len(vals)
    het_density = n_hets / span_mb if span_mb > 0 else 99
    
    if het_density < 3.0:
        threshold = 0.03  # Sparse: gentler trim
    # else: keep the passed-in threshold (default 0.045)

    # --- ADAPTIVE DECISION: Should we trim? ---
    median_ai = np.median(vals)
    q25_ai = np.percentile(vals, 25)
    q75_ai = np.percentile(vals, 75)
    
    # Strong signal: skip trimming entirely
    if span_mb > 15.0 and q25_ai > 0.10 and (q75_ai - q25_ai) < 0.10:
        return original_start, original_end
    
    # Large event with subtle but statistically significant signal: skip trimming
    # Relaxed: lowered span threshold and z-score requirement for low-fraction events
    if n_hets >= 30 and span_mb > 10.0:
        ai_se = np.std(vals) / np.sqrt(n_hets)
        if ai_se > 0 and (median_ai - 0.03) / ai_se > 3.0:
            if median_ai < 0.15:
                return original_start, original_end
    
    # Very large events (>20 Mb) with any consistent signal: skip trimming
    if span_mb > 20.0 and n_hets >= 20 and median_ai > 0.04:
        return original_start, original_end

    # --- STAGE 1: THE CHAINSAW (Split at Centromeres) ---
    positions = het_seg['position'].to_numpy()
    dist_diffs = np.diff(positions)
    gap_indices = np.where(dist_diffs > 5_000_000)[0]

    if len(gap_indices) > 0:
        chunks = []
        start_idx = 0
        
        for gap_idx in gap_indices:
            end_idx = gap_idx + 1 
            chunks.append(het_seg.iloc[start_idx:end_idx])
            start_idx = end_idx
        chunks.append(het_seg.iloc[start_idx:])

        best_chunk = None
        max_chunk_signal = -1

        for chunk in chunks:
            if len(chunk) < 3: continue
            chunk_signal = chunk['ai'].median()
            if chunk_signal > max_chunk_signal:
                max_chunk_signal = chunk_signal
                best_chunk = chunk
        
        if best_chunk is not None:
            het_seg = best_chunk
            vals = het_seg['ai'].to_numpy()

    # --- STAGE 2: THE SCALPEL (Precise Edge Trimming) ---
    
    smooth_vals = pd.Series(vals).rolling(window=5, center=True, min_periods=1).median().to_numpy()
    peak_idx = np.argmax(smooth_vals)
    peak_val = smooth_vals[peak_idx]

    if peak_val < threshold:
        return het_seg.iloc[0]['position'], het_seg.iloc[-1]['position']

    robust_peak = np.percentile(vals, 90)
    cutoff = max(threshold, robust_peak * 0.3)  # Reduced from 0.4 to 0.3 — less aggressive trimming at edges

    # Walk outwards - require 4 consecutive low SNPs (was 3 — more conservative trimming)
    
    # RIGHT SIDE
    right_idx = len(vals) - 1
    consecutive_low = 0
    for i in range(peak_idx, len(vals) - 1):
        if vals[i] < cutoff:
            consecutive_low += 1
            if consecutive_low >= 4:
                right_idx = max(i - 3, peak_idx)
                break
        else:
            consecutive_low = 0
            
    # LEFT SIDE
    left_idx = 0
    consecutive_low = 0
    for i in range(peak_idx, 0, -1):
        if vals[i] < cutoff:
            consecutive_low += 1
            if consecutive_low >= 4:
                left_idx = min(i + 3, peak_idx)
                break
        else:
            consecutive_low = 0

    trimmed_start = int(het_seg.iloc[left_idx]['position'])
    trimmed_end = int(het_seg.iloc[right_idx]['position'])
    
    # --- SAFETY CHECK: Don't let trim reduce span by more than 60% ---
    original_span = (original_end - original_start) / 1e6
    trimmed_span = (trimmed_end - trimmed_start) / 1e6
    
    if original_span > 0 and trimmed_span < original_span * 0.4:
        # Trim is too aggressive — keep original bounds
        return original_start, original_end

    return trimmed_start, trimmed_end

def calculate_overlap(start1, end1, start2, end2):
    """Calculate reciprocal overlap between two regions"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start >= overlap_end:
        return 0.0
    
    overlap_len = overlap_end - overlap_start
    len1 = end1 - start1
    len2 = end2 - start2
    
    overlap1 = overlap_len / len1 if len1 > 0 else 0
    overlap2 = overlap_len / len2 if len2 > 0 else 0
    
    return min(overlap1, overlap2)

def snap_telomeric_boundaries(calls_df, df_by_chrom, chromosome_sizes, centromere_dict=None, max_gap_mb=3.0):
    if calls_df.empty:
        return calls_df
    
    calls_df = calls_df.copy()
    
    for idx, call in calls_df.iterrows():
        chrom = call['chromosome']
        
        if chrom in df_by_chrom:
            chrom_probes = df_by_chrom[chrom]
            first_probe = chrom_probes['position'].min()
            last_probe = chrom_probes['position'].max()
        else:
            continue
        
        chr_end = chromosome_sizes.get(chrom, last_probe)
        
        centromere_start = None
        centromere_end = None
        if centromere_dict and chrom in centromere_dict:
            centromere_start = centromere_dict[chrom]['start']
            centromere_end = centromere_dict[chrom]['end']
        
        call_center = (call['start_pos'] + call['end_pos']) / 2
        
        # Snap left — only if call is on p-arm and no normal hets in the gap
        if (call['start_pos'] - first_probe) / 1e6 < max_gap_mb:
            if centromere_start is None or call_center < centromere_start:
                # Check for normal hets between chromosome start and call start
                gap_snps = chrom_probes[(chrom_probes['position'] < call['start_pos']) & 
                                        (chrom_probes['BAF'] >= 0.35) & 
                                        (chrom_probes['BAF'] <= 0.65)]
                if len(gap_snps) < 3:
                    calls_df.at[idx, 'start_pos'] = 1
        
        # Snap right — only if call is on q-arm and no normal hets in the gap
        if (last_probe - call['end_pos']) / 1e6 < max_gap_mb:
            if centromere_end is None or call_center > centromere_end:
                gap_snps = chrom_probes[(chrom_probes['position'] > call['end_pos']) & 
                                        (chrom_probes['BAF'] >= 0.35) & 
                                        (chrom_probes['BAF'] <= 0.65)]
                if len(gap_snps) < 3:
                    calls_df.at[idx, 'end_pos'] = int(chr_end)
        
        calls_df.at[idx, 'span_mb'] = (calls_df.at[idx, 'end_pos'] - calls_df.at[idx, 'start_pos']) / 1e6
    
    return calls_df

def refine_boundaries_by_lrr(d, start_pos, end_pos, event_type,
                              search_margin_mb=3.0, window_snps=10):
    """Find sharpest LRR transition to tighten GAIN/LOSS boundaries."""
    search_start = max(0, start_pos - search_margin_mb * 1e6)
    search_end = end_pos + search_margin_mb * 1e6

    region = d[(d["position"] >= search_start) &
               (d["position"] <= search_end)].sort_values("position")

    if len(region) < window_snps * 3:
        return start_pos, end_pos

    lrr_vals = region["lrr_filt"].to_numpy()
    positions = region["position"].to_numpy()

    rolling_lrr = pd.Series(lrr_vals).rolling(
        window=window_snps, center=True, min_periods=5
    ).median().to_numpy()

    lrr_diff = np.diff(rolling_lrr)
    mid = (start_pos + end_pos) / 2

    if event_type == "LOSS":
        # Left: steepest drop
        left_mask = positions[1:] < mid
        if left_mask.sum() > 0:
            left_diffs = np.where(left_mask, lrr_diff, 0)
            idx = np.argmin(left_diffs)
            if left_diffs[idx] < -0.03 and positions[idx] > start_pos:
                start_pos = int(positions[idx])
        # Right: steepest rise
        right_mask = positions[1:] > mid
        if right_mask.sum() > 0:
            right_diffs = np.where(right_mask, lrr_diff, 0)
            idx = np.argmax(right_diffs)
            if right_diffs[idx] > 0.03 and positions[idx] < end_pos:
                end_pos = int(positions[idx])

    elif event_type == "GAIN":
        left_mask = positions[1:] < mid
        if left_mask.sum() > 0:
            left_diffs = np.where(left_mask, lrr_diff, 0)
            idx = np.argmax(left_diffs)
            if left_diffs[idx] > 0.03 and positions[idx] > start_pos:
                start_pos = int(positions[idx])
        right_mask = positions[1:] > mid
        if right_mask.sum() > 0:
            right_diffs = np.where(right_mask, lrr_diff, 0)
            idx = np.argmin(right_diffs)
            if right_diffs[idx] < -0.03 and positions[idx] < end_pos:
                end_pos = int(positions[idx])

    return int(start_pos), int(end_pos)

def refine_boundaries_changepoint(
    d,                      # Full chromosome DataFrame (sorted by position)
    rough_start,            # State machine start_pos (bp)
    rough_end,              # State machine end_pos (bp)
    enter_thresh=0.03,      # AI threshold used by the state machine
    search_margin_mb=5.0,   # How far beyond rough boundaries to search (Mb)
    min_segment_hets=5,     # Minimum hets in a segment to consider it "event"
    pen_multiplier=3.0,     # Penalty multiplier (higher = fewer changepoints)
    use_lrr=False,          # Also use LRR signal (better for GAIN/LOSS)
    lrr_weight=0.5,         # Weight of LRR in combined signal (0-1)
    centromere_dict=None,
    chrom=None,
    DEBUG=False,
):
    """
    Refine mCA boundaries using changepoint detection on per-SNP AI values.
    """

    rough_span = rough_end - rough_start
    search_margin = max(search_margin_mb * 1e6, rough_span * 0.3)

    # --- 1. Define search region ---
    search_start = max(0, rough_start - search_margin)
    search_end = rough_end + search_margin

    # Respect centromere: don't search across it
    if centromere_dict and chrom and chrom in centromere_dict:
        cen_start = centromere_dict[chrom]['start']
        cen_end = centromere_dict[chrom]['end']
        rough_mid = (rough_start + rough_end) / 2

        if rough_mid < cen_start:
            search_end = min(search_end, cen_start)
        elif rough_mid > cen_end:
            search_start = max(search_start, cen_end)

    # --- 2. Extract het SNP data in search region ---
    mask = (
        (d["position"] >= search_start)
        & (d["position"] <= search_end)
        & d["is_het_like"]
        & d["maf"].notna()
    )
    hets = d[mask].sort_values("position").reset_index(drop=True)

    if len(hets) < 15:
        if DEBUG:
            print(f"    [CPD] Too few hets in search region ({len(hets)}), keeping rough boundaries")
        return rough_start, rough_end

    positions = hets["position"].to_numpy(dtype=np.int64)
    ai_vals = hets["ai"].to_numpy(dtype=np.float64)

    # --- 3. Build signal for changepoint detection ---
    signal = ai_vals.copy().reshape(-1, 1)

    if use_lrr and "lrr_filt" in hets.columns:
        lrr_vals = hets["lrr_filt"].to_numpy(dtype=np.float64)
        lrr_valid = np.isfinite(lrr_vals)
        if lrr_valid.sum() > len(lrr_vals) * 0.5:
            lrr_clean = np.where(lrr_valid, lrr_vals, 0.0)
            lrr_abs = np.abs(lrr_clean)
            ai_range = np.percentile(ai_vals, 95) - np.percentile(ai_vals, 5)
            lrr_range = np.percentile(lrr_abs, 95) - np.percentile(lrr_abs, 5)
            if lrr_range > 0:
                lrr_scaled = lrr_abs * (ai_range / lrr_range)
            else:
                lrr_scaled = lrr_abs
            signal = np.column_stack([
                ai_vals * (1 - lrr_weight),
                lrr_scaled * lrr_weight
            ])

    # --- 4. Estimate penalty from background noise ---
    outside_mask = (positions < rough_start) | (positions > rough_end)
    bg_ai = ai_vals[outside_mask]

    if len(bg_ai) >= 20:
        bg_var = np.var(bg_ai)
    else:
        med = np.median(ai_vals)
        mad = np.median(np.abs(ai_vals - med))
        bg_var = (1.4826 * mad) ** 2

    # n must be defined BEFORE adaptive penalty and penalty calculation
    n = len(ai_vals)

    # Adaptive: lower penalty for sparse data to find subtle boundaries
    hets_per_mb = n / ((search_end - search_start) / 1e6) if search_end > search_start else n
    if hets_per_mb < 5:
        pen_multiplier = max(1.5, pen_multiplier * 0.5)

    # Penalty scales with variance and log(n) — BIC-like
    penalty = pen_multiplier * bg_var * np.log(n)

    if DEBUG:
        print(f"    [CPD] Search region: {search_start/1e6:.1f}-{search_end/1e6:.1f} Mb, "
              f"{len(hets)} hets, bg_var={bg_var:.5f}, penalty={penalty:.5f}, "
              f"hets_per_mb={hets_per_mb:.1f}, pen_mult={pen_multiplier:.1f}")

    # --- 5. Run changepoint detection ---
    try:
        algo = rpt.Pelt(model="l2", min_size=max(5, n // 20), jump=1).fit(signal)
        changepoints = algo.predict(pen=penalty)
    except Exception as e:
        if DEBUG:
            print(f"    [CPD] ruptures failed ({e}), keeping rough boundaries")
        return rough_start, rough_end

    bkps = [cp for cp in changepoints if cp < n]

    if DEBUG:
        bkp_positions_mb = [positions[min(bp, n-1)] / 1e6 for bp in bkps]
        print(f"    [CPD] Found {len(bkps)} changepoints at: {bkp_positions_mb}")

    if len(bkps) == 0:
        if np.median(ai_vals) > enter_thresh * 1.5:
            return int(positions[0]), int(positions[-1])
        return rough_start, rough_end

    # --- 6. Select the best pair of changepoints ---
    all_bkps = sorted([0] + bkps + [n])
    segments = []
    for i in range(len(all_bkps) - 1):
        seg_start_idx = all_bkps[i]
        seg_end_idx = all_bkps[i + 1]
        seg_ai = ai_vals[seg_start_idx:seg_end_idx]
        seg_positions = positions[seg_start_idx:seg_end_idx]

        segments.append({
            "start_idx": seg_start_idx,
            "end_idx": seg_end_idx,
            "start_pos": int(seg_positions[0]) if len(seg_positions) > 0 else 0,
            "end_pos": int(seg_positions[-1]) if len(seg_positions) > 0 else 0,
            "median_ai": float(np.median(seg_ai)),
            "mean_ai": float(np.mean(seg_ai)),
            "n_hets": len(seg_ai),
        })

    if DEBUG:
        for s in segments:
            print(f"    [CPD]   Segment {s['start_pos']/1e6:.1f}-{s['end_pos']/1e6:.1f} Mb: "
                  f"median_ai={s['median_ai']:.4f}, n={s['n_hets']}")

    # Classify segments as "event" or "background"
    bg_median = float(np.median(bg_ai)) if len(bg_ai) >= 10 else 0.015
    bg_sigma = float(1.4826 * np.median(np.abs(bg_ai - bg_median))) if len(bg_ai) >= 10 else 0.015

    event_threshold = bg_median + max(2.0 * bg_sigma, 0.015)

    for seg in segments:
        seg_overlaps_rough = (seg["end_pos"] >= rough_start) and (seg["start_pos"] <= rough_end)
        seg["is_event"] = (
            seg["median_ai"] > event_threshold
            and seg["n_hets"] >= min_segment_hets
            and seg_overlaps_rough
        )

    if DEBUG:
        print(f"    [CPD] Event threshold: {event_threshold:.4f} "
              f"(bg_median={bg_median:.4f}, bg_sigma={bg_sigma:.4f})")
        for s in segments:
            if s["is_event"]:
                print(f"    [CPD]   EVENT segment: {s['start_pos']/1e6:.1f}-{s['end_pos']/1e6:.1f} Mb")

    event_segments = [s for s in segments if s["is_event"]]

    if len(event_segments) == 0:
        lenient_threshold = bg_median + bg_sigma
        event_segments = [
            s for s in segments
            if s["median_ai"] > lenient_threshold
            and s["n_hets"] >= 5
            and (s["end_pos"] >= rough_start) and (s["start_pos"] <= rough_end)
        ]
        if len(event_segments) == 0:
            if DEBUG:
                print(f"    [CPD] No event segments even at lenient threshold, keeping rough")
            return rough_start, rough_end

    # Merge adjacent event segments (allow small gaps between them)
    event_segments.sort(key=lambda s: s["start_pos"])
    merged = [event_segments[0].copy()]
    for seg in event_segments[1:]:
        prev = merged[-1]
        gap_mb = (seg["start_pos"] - prev["end_pos"]) / 1e6
        if gap_mb < 3.0:
            prev["end_pos"] = seg["end_pos"]
            prev["end_idx"] = seg["end_idx"]
            prev["n_hets"] += seg["n_hets"]
        else:
            merged.append(seg.copy())

    # Pick the merged event region that best overlaps the rough interval
    best_region = None
    best_overlap = -1
    for region in merged:
        overlap_start = max(region["start_pos"], rough_start)
        overlap_end = min(region["end_pos"], rough_end)
        overlap = max(0, overlap_end - overlap_start)
        if overlap > best_overlap:
            best_overlap = overlap
            best_region = region

    if best_region is None:
        return rough_start, rough_end

    refined_start = best_region["start_pos"]
    refined_end = best_region["end_pos"]

    # --- 7. Fine-tune: place boundary at midpoint between last event het
    #         and first background het ---
    start_idx = best_region["start_idx"]
    end_idx = best_region["end_idx"]

    if start_idx > 0:
        prev_het_pos = int(positions[start_idx - 1])
        first_event_pos = int(positions[start_idx])
        refined_start = (prev_het_pos + first_event_pos) // 2

    if end_idx < n:
        last_event_pos = int(positions[end_idx - 1])
        next_het_pos = int(positions[min(end_idx, n - 1)])
        if end_idx < n:
            refined_end = (last_event_pos + next_het_pos) // 2
        else:
            refined_end = last_event_pos

    # CPD should only tighten boundaries, never widen beyond the state machine
    refined_start = max(refined_start, rough_start)
    refined_end = min(refined_end, rough_end)

    if DEBUG:
        print(f"    [CPD] Refined: {refined_start/1e6:.1f}-{refined_end/1e6:.1f} Mb "
              f"(rough was {rough_start/1e6:.1f}-{rough_end/1e6:.1f} Mb)")

    return int(refined_start), int(refined_end)

def add_truth_region_overlay(fig, truth_start, truth_end, truth_type, truth_fraction=None):
    """
    Overlay the expected truth region on a plot with hatched shading and dashed borders.
    Uses event-appropriate colors with white hatching for a distinct 'expected' appearance.
    """
    event_hatch_colors = {
        "CN-LOH": orange3,
        "CNLOH": orange3,
        "GAIN": red4,
        "LOSS": blue4,
    }
    
    # Normalise truth_type to match event_hatch_colors
    truth_type_upper = truth_type.upper().replace("-", "").replace("_", "")
    color = orange3  # default
    for key, val in event_hatch_colors.items():
        if key.upper().replace("-", "").replace("_", "") == truth_type_upper:
            color = val
            break
    
    # Get data axes (first 4 panels: LRR, BAF, AI, window AI)
    data_axes = [ax for ax in fig.axes[:4]]
    
    for ax in data_axes:
        # Hatched region — semi-transparent fill with diagonal hatching
        ax.axvspan(truth_start, truth_end, color=color, alpha=0.06, zorder=0)
        rect = mpatches.Rectangle(
            (truth_start, ax.get_ylim()[0]),
            truth_end - truth_start,
            ax.get_ylim()[1] - ax.get_ylim()[0],
            linewidth=0, fill=True,
            facecolor=color, alpha=0.10,
            hatch='///', edgecolor=color,
            zorder=0
        )
        ax.add_patch(rect)
        
        # Dashed vertical borders
        ax.axvline(truth_start, color=color, linestyle='--', linewidth=1.5, alpha=0.7, zorder=5)
        ax.axvline(truth_end, color=color, linestyle='--', linewidth=1.5, alpha=0.7, zorder=5)
    
    # Label above top panel
    top_ax = fig.axes[0]
    trans = mtransforms.blended_transform_factory(top_ax.transData, top_ax.transAxes)
    mid = (truth_start + truth_end) / 2
    size_mb = (truth_end - truth_start) / 1e6
    
    label = f"EXPECTED: {truth_type} {size_mb:.0f} Mb"
    if truth_fraction is not None:
        label += f" ({truth_fraction*100:.0f}%)"
    
    top_ax.text(mid, 1.08, label, transform=trans, ha='center', va='bottom',
                fontsize=9, color=color, fontstyle='italic', clip_on=False)

#Detect small (>3 Mb) ~100% CN-LOH and LOSS via het depletion (i.e. gaps in hets)
def detect_small_homozygous(df, chrom, min_span_mb=3, max_het_density=0.001, window_size_mb=2, strict_het_lo=0.2, strict_het_hi=0.8, centromere_dict=centromere_dict):
    """
    Detects small 100% CN-LOH and LOSS events (3-10 Mb) using fine-grained sliding windows.
    """
    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)
    
    if len(d) < 50:
        return pd.DataFrame()
    
    d['strict_het'] = (d['BAF'] >= strict_het_lo) & (d['BAF'] <= strict_het_hi)

    # d['is_het'] = d['is_het_like']
    positions = d['position'].to_numpy()
    
    chr_start = positions.min()
    chr_end = positions.max()
    chr_length_mb = (chr_end - chr_start) / 1e6
    
    # Small sliding windows
    window_size = window_size_mb * 1e6
    step = 1e6  # 1 Mb steps

    # Calculate chromosome-wide SNP density baseline
    total_snps = len(d)
    # Calculate baseline from median window density (robust to outlier clusters)
    window_densities = []
    pos = chr_start
    while pos < chr_end:
        window_end = pos + window_size
        mask = (positions >= pos) & (positions < window_end)
        n = mask.sum()
        if n >= 5:
            wmb = (min(window_end, chr_end) - pos) / 1e6
            window_densities.append(n / wmb)
        pos += step
    
    baseline_snp_density = np.median(window_densities) if window_densities else total_snps / chr_length_mb
    
    windows = []
    pos = chr_start
    
    while pos < chr_end:
        window_end = pos + window_size
        
        mask = (positions >= pos) & (positions < window_end)
        window_snps = d[mask]

        if len(window_snps) >= 5:  # Lower minimum for small windows
            try:
                n_total = len(window_snps)
                n_hets = window_snps['strict_het'].sum()
                het_density = n_hets / n_total if n_total > 0 else 0
                
                actual_window_size_mb = (min(window_end, chr_end) - pos) / 1e6
                snp_density = n_total / actual_window_size_mb
                
                median_lrr = np.nanmedian(window_snps['lrr_filt'])
                
                snp_ratio = snp_density / baseline_snp_density if baseline_snp_density > 0 else 0
                is_low_het = (het_density <= max_het_density) and (snp_ratio > 0.30 or median_lrr < -0.3) #snp_ratio was >0.6

                windows.append({
                    'start': int(pos),
                    'end': int(window_end),
                    'het_density': het_density,
                    'snp_ratio': snp_ratio,
                    'median_lrr': median_lrr,
                    'is_low_het': is_low_het
                })

            except Exception as e:
                print(f"    ERROR: {e}")
        
        pos += step
    
    if not windows:
        return pd.DataFrame()
    
    wdf = pd.DataFrame(windows)
    
    # Find contiguous low-het regions
    wdf['region_id'] = (wdf['is_low_het'] != wdf['is_low_het'].shift()).cumsum()
    
    calls = []
    for region_id, grp in wdf.groupby('region_id'):

        if not grp['is_low_het'].iloc[0]:
            continue
        
        start_pos = int(grp['start'].min())
        end_pos = int(grp['end'].max())
    
        span_mb = (end_pos - start_pos) / 1e6
        
        if span_mb < min_span_mb:
            continue
        
        # === EXTEND TO NORMAL HETS ===
        # Find where normal hets (BAF ~0.5) reappear
        normal_het_baf_low = 0.35
        normal_het_baf_high = 0.65
        
        # LEFT: Find where normal hets start
        left_hets = d[(d['position'] < start_pos) & d['strict_het'] & d['BAF'].notna()]
        if len(left_hets) > 0:
            left_hets = left_hets.sort_values('position', ascending=False)
            for idx, row in left_hets.head(100).iterrows():
                if normal_het_baf_low <= row['BAF'] <= normal_het_baf_high:
                    # Found normal het - set boundary after it
                    next_snp = d[d['position'] > row['position']]
                    if len(next_snp) > 0:
                        start_pos = int(next_snp.iloc[0]['position'])
                    break
        
        # RIGHT: Find where normal hets end
        right_hets = d[(d['position'] > end_pos) & d['strict_het'] & d['BAF'].notna()]
        if len(right_hets) > 0:
            right_hets = right_hets.sort_values('position', ascending=True)
            for idx, row in right_hets.head(100).iterrows():
                if normal_het_baf_low <= row['BAF'] <= normal_het_baf_high:
                    # Found normal het - set boundary before it
                    prev_snp = d[d['position'] < row['position']]
                    if len(prev_snp) > 0:
                        end_pos = int(prev_snp.iloc[-1]['position'])
                    break
        # ==============================
        
        span_mb = (end_pos - start_pos) / 1e6
        if span_mb < min_span_mb:
            continue

        # Split at centromere only if signal is one-sided
        if centromere_dict and chrom in centromere_dict:
            cen_start = centromere_dict[chrom]['start']
            cen_end = centromere_dict[chrom]['end']

            if start_pos < cen_start and end_pos > cen_end:
                # Check if both arms have low-het signal
                left_snps = d[(d['position'] >= start_pos) & (d['position'] < cen_start)]
                right_snps = d[(d['position'] > cen_end) & (d['position'] <= end_pos)]
                
                left_het_density = left_snps['strict_het'].mean() if len(left_snps) > 5 else 1.0
                right_het_density = right_snps['strict_het'].mean() if len(right_snps) > 5 else 1.0
                
                # If an arm has too few SNPs, treat it as NOT having signal
                left_low = left_het_density < 0.10 if len(left_snps) > 5 else False
                right_low = right_het_density < 0.10 if len(right_snps) > 5 else False
                
                if not left_low and not right_low:
                    continue  # Neither arm has real signal — it's just the centromere gap

                if left_low and not right_low:
                    end_pos = cen_start
                    span_mb = (end_pos - start_pos) / 1e6
                    if span_mb < min_span_mb:
                        continue
                elif right_low and not left_low:
                    start_pos = cen_end
                    span_mb = (end_pos - start_pos) / 1e6
                    if span_mb < min_span_mb:
                        continue
                # If both arms are low-het, keep the full span (genuine cross-centromere event)

            # If event overlaps centromere and has low probe count, reject
            if start_pos < cen_end and end_pos > cen_start:
                # Event spans centromere — require higher evidence
                # Recalculate with final boundaries
                cen_region_snps = d[(d['position'] >= start_pos) & (d['position'] <= end_pos)]
                cen_n_total = len(cen_region_snps)
                cen_density = cen_n_total / span_mb if span_mb > 0 else 0
                
                if cen_density < 5.0:
                    continue
                if cen_n_total < 30:
                    continue

        # Final validation
        region_snps = d[(d['position'] >= start_pos) & (d['position'] <= end_pos)]
        n_hets_actual = region_snps['strict_het'].sum()
        n_total_actual = len(region_snps)
        final_het_density = n_hets_actual / n_total_actual if n_total_actual > 0 else 1
        
        # Reject probe-sparse regions (centromeric/gap artefacts)
        snp_density = n_total_actual / span_mb if span_mb > 0 else 0
        if snp_density < 3.0:
            continue

        # Require truly zero hets when probe count is low
        if n_total_actual < 30:
            if n_hets_actual > 0:
                continue
        else:
            if final_het_density > 0.05:
                continue
        
        if n_total_actual < 20:  # Lower for small events (previously 30)
            continue
        
        # Check homozygous ratio
        hom_snps = region_snps[~region_snps['strict_het']]
        if len(hom_snps) < 0.90 * n_total_actual:
            continue
        
        median_lrr = np.nanmedian(region_snps['lrr_filt'])
        
        if np.isnan(median_lrr):
            continue
        
        # Calculate het depletion for cell fraction estimate
        bg_snps = d[~((d['position'] >= start_pos) & (d['position'] <= end_pos))]
        bg_hets = bg_snps[(bg_snps['BAF'] >= strict_het_lo) & (bg_snps['BAF'] <= strict_het_hi)]
        bg_het_fraction = len(bg_hets) / len(bg_snps) if len(bg_snps) > 0 else 0.3
        region_het_fraction = n_hets_actual / n_total_actual if n_total_actual > 0 else 0
        het_depletion = 1 - (region_het_fraction / bg_het_fraction) if bg_het_fraction > 0 else 1.0
        het_depletion = float(np.clip(het_depletion, 0, 1))
        
        # Classify
        # Compare region LRR to chromosome background
        bg_median_lrr = np.nanmedian(bg_snps['lrr_filt'])
        bg_lrr_sd = np.nanstd(bg_snps['lrr_filt'])

        lrr_contrast = median_lrr - bg_median_lrr

        # Only call LOSS if region LRR is significantly below background
        if lrr_contrast < -0.2 and median_lrr < -0.3 and (bg_lrr_sd == 0 or abs(lrr_contrast) / bg_lrr_sd > 2.0):
            event = 'LOSS'
            cf_lrr = 2 * (1 - 2**median_lrr)
            cf_lrr = float(np.clip(cf_lrr, 0, 1))
            cf_baf = het_depletion
        else:
            event = 'CN-LOH'
            cf_lrr = np.nan
            cf_baf = het_depletion
        
        cf_combined = (cf_lrr + cf_baf) / 2 if not np.isnan(cf_lrr) else cf_baf

        # Refine cell fraction using AI if hets are present in wider definition
        region_hets_wide = region_snps[region_snps['is_het_like'] & region_snps['maf'].notna()]
        if len(region_hets_wide) >= 10:
            mean_ai = region_hets_wide['ai'].mean()
            if event == 'CN-LOH':
                cf_from_ai = min(2 * mean_ai, 1.0)
            elif event == 'LOSS':
                cf_from_ai = min(4 * mean_ai / (1 + 2 * mean_ai), 1.0)
            else:
                cf_from_ai = cf_baf
            
            # Use AI-based estimate if it's lower than het depletion (more conservative)
            if cf_from_ai < cf_baf:
                cf_baf = cf_from_ai

        # Cell fraction agreement check
        if not np.isnan(cf_lrr):
            cf_diff = abs(cf_lrr - cf_baf)
            if cf_diff > 0.4:
                continue
        
        calls.append({
            'chromosome': chrom,
            'start_pos': start_pos,
            'end_pos': end_pos,
            'span_mb': span_mb,
            'event': event,
            'cell_fraction_BAF': cf_baf,
            'cell_fraction_LRR': cf_lrr,
            'cell_fraction_combined': cf_combined,
            'n_hets': int(n_hets_actual),
            'n_total_snps': int(n_total_actual),
            'het_density': float(final_het_density),
            'median_LRR': float(median_lrr),
            'mean_LRR': float(np.nanmean(region_snps['lrr_filt'])),
            'confidence': 'HIGH',
            'detector': 'small_homozygous',
            'note': f'homozygous small event (het density={final_het_density:.3f})'
        })
        
    return pd.DataFrame(calls)

#Detect large (>10 Mb) ~100% CN-LOH and LOSS via het depletion (i.e. gaps in hets)
def detect_large_homozygous(df, chrom, min_span_mb=10, max_het_density=0.01, window_size_mb=5, strict_het_lo=0.2, strict_het_hi=0.8, centromere_dict=centromere_dict):
    """
    Detects 100% cell fraction events with smart boundary extension.
    Extends until finding normal hets, or to chromosome arm ends if no normal hets found.
    Uses STRICT het definition to avoid counting borderline homozygous SNPs.
    """
    
    # Create strict het column
    df = df.copy()
    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)
    
    if len(d) < 50:
        return pd.DataFrame()
    
    d['strict_het'] = (d['BAF'] >= strict_het_lo) & (d['BAF'] <= strict_het_hi)
    positions = d['position'].to_numpy()
    
    # Get chromosome boundaries
    chr_start = positions.min()
    chr_end = positions.max()
    
    # Get centromere position if available
    centromere_start = None
    centromere_end = None
    if centromere_dict and chrom in centromere_dict:
        centromere_start = centromere_dict[chrom]['start']
        centromere_end = centromere_dict[chrom]['end']
    
    window_size = window_size_mb * 1e6
    
    # Create windows
    windows = []
    step = window_size / 2
    
    pos = chr_start
    while pos < chr_end:
        window_end = pos + window_size
        
        mask = (positions >= pos) & (positions < window_end)
        window_snps = d[mask]
        
        if len(window_snps) >= 20:
            n_total = len(window_snps)
            n_hets = window_snps['strict_het'].sum()
            het_density = n_hets / n_total
            
            median_lrr = np.nanmedian(window_snps['lrr_filt'])
            
            windows.append({
                'start': int(pos),
                'end': int(window_end),
                'het_density': het_density,
                'median_lrr': median_lrr,
                'is_low_het': het_density <= max_het_density
            })
        
        pos += step
    
    if not windows:
        return pd.DataFrame()
    
    wdf = pd.DataFrame(windows)

    # Find contiguous regions
    wdf['region_id'] = (wdf['is_low_het'] != wdf['is_low_het'].shift()).cumsum()
    
    preliminary_calls = []
    for region_id, grp in wdf.groupby('region_id'):
        if not grp['is_low_het'].iloc[0]:
            continue
        
        # Initial boundaries from windows
        start_pos = int(grp['start'].min())
        end_pos = int(grp['end'].max())
        span_mb = (end_pos - start_pos) / 1e6
    
        # === NEW: SPLIT HETEROGENEOUS LRR REGIONS ===
        region_windows = grp

        if len(region_windows) > 3:
            lrr_values = region_windows['median_lrr'].values
            lrr_std = np.nanstd(lrr_values)
            lrr_range = np.nanmax(lrr_values) - np.nanmin(lrr_values)
            
            if lrr_std > 0.25 or lrr_range > 0.6:
                split_indices = []
                for i in range(len(region_windows) - 1):
                    lrr_current = region_windows.iloc[i]['median_lrr']
                    lrr_next = region_windows.iloc[i+1]['median_lrr']
                    
                    if not np.isnan(lrr_current) and not np.isnan(lrr_next):
                        if abs(lrr_current - lrr_next) > 0.4:
                            split_indices.append(i)
                
                if split_indices:
                    sub_regions = []
                    start_idx = 0
                    
                    for split_idx in split_indices:
                        sub_windows = region_windows.iloc[start_idx:split_idx+1]
                        sub_start = int(sub_windows['start'].min())
                        sub_end = int(sub_windows['end'].max())
                        sub_span = (sub_end - sub_start) / 1e6
                        
                        if sub_span >= min_span_mb:
                            sub_regions.append((sub_start, sub_end, sub_span))
                        
                        start_idx = split_idx + 1
                    
                    # Add final region
                    sub_windows = region_windows.iloc[start_idx:]
                    sub_start = int(sub_windows['start'].min())
                    sub_end = int(sub_windows['end'].max())
                    sub_span = (sub_end - sub_start) / 1e6
                    
                    if sub_span >= min_span_mb:
                        sub_regions.append((sub_start, sub_end, sub_span))
                    
                    for sub_start, sub_end, sub_span in sub_regions:
                        preliminary_calls.append({
                            'start_pos': sub_start,
                            'end_pos': sub_end,
                            'span_mb': sub_span
                        })
                    
                    continue  # Skip normal processing
        # ========================================

        if span_mb < min_span_mb:
            continue

        # Save core boundaries from windows
        core_start = start_pos
        core_end = end_pos
        max_extension_mb = 5.0  # Never extend more than 5 Mb beyond core windows
        
        # === EXTEND TO NORMAL HETS OR CHROMOSOME ARM ENDS ===
        region_center = (start_pos + end_pos) / 2
        event_arm = None
        if centromere_start and centromere_end:
            if region_center < centromere_start:
                event_arm = 'p'
            elif region_center > centromere_end:
                event_arm = 'q'
            else:
                event_arm = 'both'
        
        all_hets = d[d['strict_het'] & d['BAF'].notna()]
        
        if len(all_hets) > 100:
            normal_het_baf_low = 0.35
            normal_het_baf_high = 0.65

            core_region = d[(d['position'] >= start_pos) & (d['position'] <= end_pos)]
            core_median_lrr = np.nanmedian(core_region['lrr_filt'])
            
            # LEFT EXTENSION
            found_left_boundary = False

            if event_arm == 'q' and centromere_end:
                if start_pos <= centromere_end + 2e6:
                    start_pos = max(start_pos, int(centromere_end))
                    found_left_boundary = True

            if not found_left_boundary:
                left_window_start = max(chr_start, start_pos - 8e6)
                left_window = d[(d['position'] >= left_window_start) & (d['position'] < start_pos)]
                if len(left_window) > 5:
                    left_median_lrr = np.nanmedian(left_window['lrr_filt'])
                    if not np.isnan(left_median_lrr) and not np.isnan(core_median_lrr):
                        lrr_diff = abs(left_median_lrr - core_median_lrr)
                        lrr_stop_thresh = min(0.4, abs(core_median_lrr) * 0.5) if abs(core_median_lrr) > 0.2 else 0.4
                        if lrr_diff > lrr_stop_thresh:
                            found_left_boundary = True
            
            if not found_left_boundary:
                left_search = d[d['position'] < start_pos].copy()
                left_hets = left_search[left_search['strict_het'] & left_search['BAF'].notna()]
                
                if len(left_hets) > 0:
                    left_hets = left_hets.sort_values('position', ascending=False)
                    
                    for idx, row in left_hets.head(200).iterrows():
                        if normal_het_baf_low <= row['BAF'] <= normal_het_baf_high:
                            next_snps = d[d['position'] > row['position']]
                            if len(next_snps) > 0:
                                start_pos = int(next_snps.iloc[0]['position'])
                                found_left_boundary = True
                            break
            
            if not found_left_boundary:
                if centromere_start and centromere_end:
                    if event_arm == 'p':
                        start_pos = int(chr_start)
                    elif event_arm == 'q':
                        start_pos = int(centromere_end)
                    else:
                        region_center = (start_pos + end_pos) / 2
                        if region_center < centromere_start:
                            start_pos = int(chr_start)
                        else:
                            start_pos = int(centromere_end)
                else:
                    start_pos = int(chr_start)
            
            # RIGHT EXTENSION
            found_right_boundary = False

            if event_arm == 'p' and centromere_start:
                if end_pos >= centromere_start - 2e6:
                    end_pos = min(end_pos, int(centromere_start))
                    found_right_boundary = True

            if not found_right_boundary:
                right_window_end = min(chr_end, end_pos + 8e6)
                right_window = d[(d['position'] > end_pos) & (d['position'] <= right_window_end)]
                if len(right_window) > 5:
                    right_median_lrr = np.nanmedian(right_window['lrr_filt'])
                    if not np.isnan(right_median_lrr) and not np.isnan(core_median_lrr):
                        lrr_diff = abs(right_median_lrr - core_median_lrr)
                        lrr_stop_thresh = min(0.4, abs(core_median_lrr) * 0.5) if abs(core_median_lrr) > 0.2 else 0.4
                        if lrr_diff > lrr_stop_thresh:
                            found_right_boundary = True
            
            if not found_right_boundary:
                right_search = d[d['position'] > end_pos].copy()
                right_hets = right_search[right_search['strict_het'] & right_search['BAF'].notna()]
                
                if len(right_hets) > 0:
                    right_hets = right_hets.sort_values('position', ascending=True)
                    
                    for idx, row in right_hets.head(200).iterrows():
                        if normal_het_baf_low <= row['BAF'] <= normal_het_baf_high:
                            prev_snps = d[d['position'] < row['position']]
                            if len(prev_snps) > 0:
                                end_pos = int(prev_snps.iloc[-1]['position'])
                                found_right_boundary = True
                            break
            
            if not found_right_boundary:
                if centromere_start and centromere_end:
                    if event_arm == 'q':
                        end_pos = int(chr_end)
                    elif event_arm == 'p':
                        end_pos = int(centromere_start)
                    else:
                        region_center = (start_pos + end_pos) / 2
                        if region_center > centromere_end:
                            end_pos = int(chr_end)
                        else:
                            end_pos = int(centromere_start)
                else:
                    end_pos = int(chr_end)
        # ==============================

        # Clamp extension to max distance from core windows
        start_pos = max(start_pos, int(core_start - max_extension_mb * 1e6))
        end_pos = min(end_pos, int(core_end + max_extension_mb * 1e6))

        # Trim to actual probe positions — don't extend into empty space
        region_snps = d[(d['position'] >= start_pos) & (d['position'] <= end_pos)]
        if len(region_snps) > 0:
            start_pos = int(region_snps['position'].min())
            end_pos = int(region_snps['position'].max())
        
        span_mb = (end_pos - start_pos) / 1e6
        
        if span_mb < min_span_mb:
            continue
        
        preliminary_calls.append({
            'start_pos': start_pos,
            'end_pos': end_pos,
            'span_mb': span_mb
        })
    
    if not preliminary_calls:
        return pd.DataFrame()
    
    # === MERGE OVERLAPPING CALLS ===
    preliminary_calls = sorted(preliminary_calls, key=lambda x: x['start_pos'])
    merged_calls = []
    current = preliminary_calls[0].copy()
    
    for call in preliminary_calls[1:]:
        gap_mb = (call['start_pos'] - current['end_pos']) / 1e6

        crosses_centromere = False
        if centromere_dict and chrom in centromere_dict:
            cen_s = centromere_dict[chrom]['start']
            cen_e = centromere_dict[chrom]['end']
            if current['end_pos'] < cen_s and call['start_pos'] > cen_e:
                crosses_centromere = True
            elif current['end_pos'] < cen_e and call['start_pos'] > cen_s:
                crosses_centromere = True

        if gap_mb < 3 and not crosses_centromere:
            current['end_pos'] = max(current['end_pos'], call['end_pos'])
            current['span_mb'] = (current['end_pos'] - current['start_pos']) / 1e6
        else:
            merged_calls.append(current)
            current = call.copy()
    
    merged_calls.append(current)
    # ================================
    
    # Final processing
    final_calls = []
    for call in merged_calls:
        start_pos = call['start_pos']
        end_pos = call['end_pos']
        span_mb = call['span_mb']
        
        # Get region stats
        region_snps = d[(d['position'] >= start_pos) & (d['position'] <= end_pos)]
        n_hets_actual = region_snps['strict_het'].sum()
        n_total_actual = len(region_snps)
        final_het_density = n_hets_actual / n_total_actual if n_total_actual > 0 else 1
        
        if final_het_density > 0.1:
            continue
        
        # Density check: require minimum SNPs per Mb
        snp_density = n_total_actual / span_mb if span_mb > 0 else 0
        if n_total_actual < 30 or snp_density < 1.5:
            continue

        # Also check that we have enough homozygous SNPs
        hom_snps = region_snps[~region_snps['strict_het']]
        hom_fraction = len(hom_snps) / n_total_actual
        
        if len(hom_snps) < 0.90 * n_total_actual:
            continue
        
        median_lrr = np.nanmedian(region_snps['lrr_filt'])
        
        if np.isnan(median_lrr):
            continue
        
        # Calculate het depletion for cell fraction estimate
        bg_snps = d[~((d['position'] >= start_pos) & (d['position'] <= end_pos))]
        bg_hets = bg_snps[(bg_snps['BAF'] >= strict_het_lo) & (bg_snps['BAF'] <= strict_het_hi)]
        bg_het_fraction = len(bg_hets) / len(bg_snps) if len(bg_snps) > 0 else 0.3
        region_het_fraction = n_hets_actual / n_total_actual if n_total_actual > 0 else 0
        het_depletion = 1 - (region_het_fraction / bg_het_fraction) if bg_het_fraction > 0 else 1.0
        het_depletion = float(np.clip(het_depletion, 0, 1))
        
        # Classify
        if median_lrr < -0.3:
            event = 'LOSS'
            cf_lrr = 2 * (1 - 2**median_lrr)
            cf_lrr = float(np.clip(cf_lrr, 0, 1))
            cf_baf = het_depletion
        else:
            event = 'CN-LOH'
            cf_lrr = np.nan
            cf_baf = het_depletion
        
        cf_combined = (cf_lrr + cf_baf) / 2 if not np.isnan(cf_lrr) else cf_baf

        # Refine cell fraction using AI if hets are present in wider definition
        region_hets_wide = region_snps[region_snps['is_het_like'] & region_snps['maf'].notna()]
        if len(region_hets_wide) >= 10:
            mean_ai = region_hets_wide['ai'].mean()
            if event == 'CN-LOH':
                cf_from_ai = min(2 * mean_ai, 1.0)
            elif event == 'LOSS':
                cf_from_ai = min(4 * mean_ai / (1 + 2 * mean_ai), 1.0)
            else:
                cf_from_ai = cf_baf
            
            # Use AI-based estimate if it's lower than het depletion (more conservative)
            if cf_from_ai < cf_baf:
                cf_baf = cf_from_ai
        
        final_calls.append({
            'chromosome': chrom,
            'start_pos': start_pos,
            'end_pos': end_pos,
            'span_mb': span_mb,
            'event': event,
            'cell_fraction_BAF': cf_baf,
            'cell_fraction_LRR': cf_lrr,
            'cell_fraction_combined': cf_combined,
            'n_hets': int(n_hets_actual),
            'n_total_snps': int(n_total_actual),
            'het_density': float(final_het_density),
            'median_LRR': float(median_lrr),
            'mean_LRR': float(np.nanmean(region_snps['lrr_filt'])),
            'confidence': 'HIGH',
            'detector': 'large_homozygous',
            'note': f'large homozygous (het density={final_het_density:.3f})'
        })
    
    return pd.DataFrame(final_calls)

#Detect small (>3 Mb) GAINS at >30% cell fraction via BAF split + LRR
def detect_small_gains(df, chrom, min_span_mb=3, window_size_mb=2, strict_het_lo=0.2, strict_het_hi=0.8):
    """
    Detects small GAIN events (2-10 Mb) using BAF split pattern.
    100% GAINs show elevated LRR + BAF split to 0.33/0.67 (AAB/ABB).
    """
    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)
    
    if len(d) < 50:
        return pd.DataFrame()
    
    # For GAIN, look for BAF splits at 0.33/0.67
    d['is_gain_het'] = ((d['BAF'] >= 0.25) & (d['BAF'] <= 0.42)) | \
                       ((d['BAF'] >= 0.58) & (d['BAF'] <= 0.75))
    
    window_size = window_size_mb * 1e6
    positions = d['position'].to_numpy()
    
    windows = []
    step = window_size / 2
    pos = positions.min()
    
    while pos < positions.max():
        window_end = pos + window_size
        mask = (positions >= pos) & (positions < window_end)
        window_snps = d[mask]
        
        if len(window_snps) >= 10:
            n_total = len(window_snps)
            n_gain_hets = window_snps['is_gain_het'].sum()
            gain_het_fraction = n_gain_hets / n_total
            
            median_lrr = np.nanmedian(window_snps['lrr_filt'])
            het_mask = (window_snps['BAF'] >= strict_het_lo) & (window_snps['BAF'] <= strict_het_hi)
            median_ai = np.nanmedian(window_snps[het_mask]['ai'])
            
            # GAIN signature: LRR > 0.15, >20% gain-hets, AI > 0.05
            is_gain = (median_lrr > 0.15 and 
                      gain_het_fraction > 0.15 and 
                      median_ai > 0.05)
            
            windows.append({
                'start': int(pos),
                'end': int(window_end),
                'median_lrr': median_lrr,
                'gain_het_fraction': gain_het_fraction,
                'median_ai': median_ai,
                'is_gain': is_gain
            })
        
        pos += step
    
    if not windows:
        return pd.DataFrame()
    
    wdf = pd.DataFrame(windows)
    wdf['region_id'] = (wdf['is_gain'] != wdf['is_gain'].shift()).cumsum()

    # After building windows, before region merging:
    # Calculate AI noise floor from non-gain windows
    non_gain_windows = wdf[~wdf['is_gain']]
    if len(non_gain_windows) > 5:
        bg_ai = non_gain_windows['median_ai'].dropna()
        ai_min_threshold = max(0.08, bg_ai.median() + 3 * bg_ai.std())
    else:
        ai_min_threshold = 0.15

    # if chrom == 'chr9':
    #     print(f"ai_min_threshold: {ai_min_threshold:.3f}")
    #     for w in windows:
    #         if w['start'] >= 20e6 and w['end'] <= 35e6:
    #             print(f"  Window {w['start']/1e6:.0f}-{w['end']/1e6:.0f}: "
    #                 f"lrr={w['median_lrr']:.3f} gain_het_frac={w['gain_het_fraction']:.3f} "
    #                 f"ai={w['median_ai']:.3f} is_gain={w['is_gain']}")
    
    calls = []
    for region_id, grp in wdf.groupby('region_id'):
        if not grp['is_gain'].iloc[0]:
            continue
        
        start_pos = int(grp['start'].min())
        end_pos = int(grp['end'].max())
        span_mb = (end_pos - start_pos) / 1e6
        
        if span_mb < min_span_mb:
            continue

        # # Require minimum SNP count
        # if len(region_snps) < 30:
        #     continue
        
        region_snps = d[(d['position'] >= start_pos) & (d['position'] <= end_pos)]
        median_lrr = np.nanmedian(region_snps['lrr_filt'])
        gain_het_fraction = region_snps['is_gain_het'].sum() / len(region_snps)
        
        het_mask = (region_snps['BAF'] >= strict_het_lo) & (region_snps['BAF'] <= strict_het_hi)
        median_ai = np.nanmedian(region_snps[het_mask]['ai'])

        # Median AI must be above noise floor
        if median_ai < ai_min_threshold:
            continue
        
        # Cell fraction estimates
        cf_lrr = 2 * (2**median_lrr - 1)
        cf_lrr = float(np.clip(cf_lrr, 0, 1))
        cf_baf = 4 * median_ai / (1 - 2 * median_ai) if median_ai < 0.5 else 1.0
        cf_baf = float(np.clip(cf_baf, 0, 1))
        cf_combined = (cf_lrr + cf_baf) / 2
        cf_diff = abs(cf_lrr - cf_baf)
        
        # Reject if LRR and BAF disagree
        if cf_diff > 0.2:
            continue
        
        # Reject if both too low
        if cf_lrr < 0.1 and cf_baf < 0.1:
            continue

        # Reject if average cell fraction <0.5:
        if cf_combined <0.3: #changed from 0.5 to 0.3
            continue

        # Require enough hets for reliable BAF assessment
        het_mask = (region_snps['BAF'] >= strict_het_lo) & (region_snps['BAF'] <= strict_het_hi)
        n_hets = het_mask.sum()
        if n_hets < 8:
            continue
        
        # NEW: Mean and median LRR must agree (catches outlier-driven false positives)
        mean_lrr = np.nanmean(region_snps['lrr_filt'])
        if mean_lrr < 0.05:  # Mean LRR must actually be elevated
            continue
        
        # Confidence
        if cf_diff < 0.15:
            confidence = 'HIGH'
        elif cf_diff < 0.25:
            confidence = 'MEDIUM'
        else:
            confidence = 'LOW'
        
        calls.append({
            'chromosome': chrom,
            'start_pos': start_pos,
            'end_pos': end_pos,
            'span_mb': span_mb,
            'event': 'GAIN',
            'confidence': confidence,
            'cell_fraction_LRR': cf_lrr,
            'cell_fraction_BAF': cf_baf,
            'cell_fraction_combined': cf_combined,
            'n_hets': int(region_snps[(region_snps['BAF'] >= strict_het_lo) & (region_snps['BAF'] <= strict_het_hi)].shape[0]),
            'n_total_snps': len(region_snps),
            'gain_het_fraction': gain_het_fraction,
            'median_LRR': float(median_lrr),
            'mean_LRR': float(np.nanmean(region_snps['lrr_filt'])),
            'median_AI': float(median_ai),
            'detector': 'small_gains',
            'note': f'gain (gain-het fraction={gain_het_fraction:.2f})'
        })
    
    return pd.DataFrame(calls)

#Detect any GAIN or LOSS at >30% cell fraction via LRR+BAF confirmation (debuggable version)
def detect_LRR_BAF(df, chrom, chromosome_sizes,
                                              centromere_dict=None,
                                              lrr_gain_thresh=0.12,
                                              lrr_loss_thresh=-0.12,
                                              min_consecutive_bins=5,
                                              min_hets_in_region=3,
                                              min_span_mb=3,
                                              baf_ai_thresh=0.03,
                                              min_baf_support_fraction=0.7,
                                              min_valid_lrr_fraction=0.7,
                                              true_het_baf_lo=0.25,
                                              true_het_baf_hi=0.75, debug_chrom=None, genome_baseline_lrr=None):
    """
    LRR-first detector with BAF confirmation.
    Collects ALL non-neutral bins, merges nearby same-state bins,
    then filters by total bin count and span.
    """
    
    # DEBUG = (chrom==debug_chrom) # <<< SET TO FALSE TO DISABLE DEBUG OUTPUT
    DEBUG = False # <<< SET TO FALSE TO DISABLE DEBUG OUTPUT
    # DEBUG = True

    if len(df) < 50:
        if DEBUG: print(f"  [LRR_BAF {chrom}] Too few SNPs: {len(df)} < 50")
        return pd.DataFrame()
    
    df = df.copy().sort_values('position').reset_index(drop=True)
    
    # Step 1: Classify LRR state per bin (first pass)
    df['lrr_smoothed'] = df['lrr_filt'].rolling(window=11, center=True, min_periods=3).median()

    all_lrr = df['lrr_filt'].dropna()
    if len(all_lrr) > 50:
        chr_baseline_lrr = np.nanmedian(all_lrr)
    else:
        chr_baseline_lrr = 0.0
    
    if DEBUG: print(f"  [LRR_BAF {chrom}] Initial baseline LRR: {chr_baseline_lrr:.4f}")
    
    df['lrr_state'] = 'NEUTRAL'
    valid_lrr_mask = df['lrr_filt'].notna()
    df.loc[valid_lrr_mask & (df['lrr_filt'] - chr_baseline_lrr > lrr_gain_thresh), 'lrr_state'] = 'GAIN'
    df.loc[valid_lrr_mask & (df['lrr_filt'] - chr_baseline_lrr < lrr_loss_thresh), 'lrr_state'] = 'LOSS'    
    
    n_gain_1 = (df['lrr_state'] == 'GAIN').sum()
    n_loss_1 = (df['lrr_state'] == 'LOSS').sum()
    if DEBUG: print(f"  [LRR_BAF {chrom}] First pass: {n_gain_1} GAIN bins, {n_loss_1} LOSS bins")
    
    # Step 1.5: Recompute baseline excluding candidate regions
    candidate_mask = df['lrr_state'] != 'NEUTRAL'
    background_lrr = df.loc[~candidate_mask, 'lrr_filt'].dropna()
    if len(background_lrr) > 50:
        chr_baseline_lrr = np.nanmedian(background_lrr)
        if DEBUG: print(f"  [LRR_BAF {chrom}] Refined baseline LRR: {chr_baseline_lrr:.4f} (from {len(background_lrr)} neutral bins)")
    
    # Step 1.6: Genome-wide baseline fallback for large events
    used_genome_baseline = False
    n_total_snps = len(df)
    n_neutral = (~candidate_mask).sum()
    neutral_fraction = n_neutral / n_total_snps if n_total_snps > 0 else 1.0
    
    if abs(chr_baseline_lrr) > 0.20 or neutral_fraction < 0.30:
        if genome_baseline_lrr is not None:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ⚠ Baseline contaminated (|{chr_baseline_lrr:.3f}| > 0.20 or neutral_frac={neutral_fraction:.2f})")
            if DEBUG: print(f"  [LRR_BAF {chrom}] → Using genome-wide baseline: {genome_baseline_lrr:.4f}")
            chr_baseline_lrr = genome_baseline_lrr
            used_genome_baseline = True
        else:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ⚠ Baseline may be contaminated but genome-wide baseline not provided")
    
    # Reclassify with final baseline
    if len(background_lrr) > 50 or used_genome_baseline:
        df['lrr_state'] = 'NEUTRAL'
        df.loc[valid_lrr_mask & (df['lrr_filt'] - chr_baseline_lrr > lrr_gain_thresh), 'lrr_state'] = 'GAIN'
        df.loc[valid_lrr_mask & (df['lrr_filt'] - chr_baseline_lrr < lrr_loss_thresh), 'lrr_state'] = 'LOSS'   
        
        n_gain_2 = (df['lrr_state'] == 'GAIN').sum()
        n_loss_2 = (df['lrr_state'] == 'LOSS').sum()
        if DEBUG: print(f"  [LRR_BAF {chrom}] Final pass: {n_gain_2} GAIN bins, {n_loss_2} LOSS bins")
    
    df['state_change'] = (df['lrr_state'] != df['lrr_state'].shift()).cumsum()
    
    # Step 2: Collect all preliminary LRR segments (no span check yet)
    preliminary_segments = []
    for state_id, state_group in df.groupby('state_change'):
        state = state_group['lrr_state'].iloc[0]
        if state == 'NEUTRAL':
            continue
        
        valid_lrr_count = state_group['lrr_filt'].notna().sum()
        valid_lrr_fraction = valid_lrr_count / len(state_group)
        
        if valid_lrr_fraction < min_valid_lrr_fraction:
            continue
        
        start = state_group['position'].min()
        end = state_group['position'].max()
        
        preliminary_segments.append({
            'state': state,
            'start': start,
            'end': end,
            'n_state_bins': len(state_group)
        })
    
    if DEBUG: 
        print(f"  [LRR_BAF {chrom}] Preliminary segments: {len(preliminary_segments)}")
        for seg in preliminary_segments:
            print(f"    {seg['state']}: {seg['start']/1e6:.1f}-{seg['end']/1e6:.1f} Mb ({seg['n_state_bins']} bins)")

    # Step 3: Merge nearby segments of same type (SEPARATELY)
    centromere_start = None
    centromere_end = None
    if centromere_dict and chrom in centromere_dict:
        centromere_start = centromere_dict[chrom]['start']
        centromere_end = centromere_dict[chrom]['end']
    
    def crosses_centromere(seg_start, seg_end):
        if centromere_start is None:
            return False
        return seg_start < centromere_start and seg_end > centromere_end
    
    def on_different_arms(seg1_end, seg2_start):
        if centromere_start is None:
            return False
        return seg1_end < centromere_start and seg2_start > centromere_end
    
    merged_segments = []
    for target_state in ['GAIN', 'LOSS']:
        state_segs = sorted(
            [s for s in preliminary_segments if s['state'] == target_state],
            key=lambda x: x['start']
        )
        if not state_segs:
            continue
        current = state_segs[0].copy()
        for next_seg in state_segs[1:]:
            gap_mb = (next_seg['start'] - current['end']) / 1e6
            if gap_mb < 3 and not on_different_arms(current['end'], next_seg['start']):
                if DEBUG: print(f"    Merging {target_state}: {current['start']/1e6:.1f}-{current['end']/1e6:.1f} + {next_seg['start']/1e6:.1f}-{next_seg['end']/1e6:.1f} (gap={gap_mb:.1f} Mb)")
                current['end'] = next_seg['end']
                current['n_state_bins'] += next_seg['n_state_bins']
            else:
                merged_segments.append(current)
                current = next_seg.copy()
        merged_segments.append(current)
    
    if DEBUG:
        print(f"  [LRR_BAF {chrom}] Merged segments: {len(merged_segments)}")
        for seg in merged_segments:
            print(f"    {seg['state']}: {seg['start']/1e6:.1f}-{seg['end']/1e6:.1f} Mb ({seg['n_state_bins']} bins, span={(seg['end']-seg['start'])/1e6:.1f} Mb)")

    # Step 3.1: Split any segment that still spans the centromere
    if centromere_start is not None:
        split_segments = []
        for seg in merged_segments:
            if crosses_centromere(seg['start'], seg['end']):
                if DEBUG: print(f"    Splitting centromere-spanning segment: {seg['start']/1e6:.1f}-{seg['end']/1e6:.1f} Mb")
                region_snps = df[(df['position'] >= seg['start']) & (df['position'] <= seg['end'])]
                p_snps = region_snps[region_snps['position'] < centromere_start]
                q_snps = region_snps[region_snps['position'] > centromere_end]
                
                if len(p_snps) >= min_consecutive_bins:
                    p_state_bins = len(p_snps[p_snps['lrr_state'] == seg['state']])
                    split_segments.append({
                        'state': seg['state'],
                        'start': seg['start'],
                        'end': int(p_snps['position'].max()),
                        'n_state_bins': p_state_bins
                    })
                if len(q_snps) >= min_consecutive_bins:
                    q_state_bins = len(q_snps[q_snps['lrr_state'] == seg['state']])
                    split_segments.append({
                        'state': seg['state'],
                        'start': int(q_snps['position'].min()),
                        'end': seg['end'],
                        'n_state_bins': q_state_bins
                    })
            else:
                split_segments.append(seg)
        merged_segments = split_segments

    # Step 3.5: Calculate AI noise floor EXCLUDING candidate regions
    candidate_mask = pd.Series(False, index=df.index)
    for segment in merged_segments:
        candidate_mask |= (df['position'] >= segment['start']) & (df['position'] <= segment['end'])
    
    background_snps = df[~candidate_mask]
    background_hets = background_snps[
        (background_snps['BAF'] >= true_het_baf_lo) & 
        (background_snps['BAF'] <= true_het_baf_hi)
    ]
    background_ai = background_hets['ai'].dropna()
    
    if len(background_ai) > 10:
        ai_noise_median = background_ai.median()
        ai_noise_std = background_ai.std()
        ai_min_threshold = ai_noise_median + 3 * ai_noise_std
    else:
        ai_noise_median = np.nan
        ai_noise_std = np.nan
        ai_min_threshold = 0.1
    
    if DEBUG: print(f"  [LRR_BAF {chrom}] AI noise floor: median={ai_noise_median:.3f}, std={ai_noise_std:.3f}, threshold={ai_min_threshold:.3f}")

    # Step 4: Process merged segments
    calls = []
    
    for seg_idx, segment in enumerate(merged_segments):
        state = segment['state']
        start = segment['start']
        end = segment['end']
        span_mb = (end - start) / 1e6
        number_bins = segment['n_state_bins']

        if DEBUG: print(f"\n  [LRR_BAF {chrom}] === Segment {seg_idx}: {state} {start/1e6:.1f}-{end/1e6:.1f} Mb ({number_bins} state bins) ===")

        if segment['n_state_bins'] < min_consecutive_bins:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: n_state_bins {number_bins} < {min_consecutive_bins}")
            continue
        
        region_snps = df[(df['position'] >= start) & (df['position'] <= end)]
        if DEBUG: print(f"  [LRR_BAF {chrom}] Region has {len(region_snps)} total SNPs")

        # Trim edges
        trim_window = max(5, len(region_snps) // 20)
        median_lrr_pretrim = region_snps['lrr_filt'].median()
        
        if state == 'GAIN':
            trim_gain_thresh = max(lrr_gain_thresh, median_lrr_pretrim * 0.5)
            trim_loss_thresh = lrr_loss_thresh * 2
        elif state == 'LOSS':
            trim_loss_thresh = min(lrr_loss_thresh, median_lrr_pretrim * 0.5)
            trim_gain_thresh = lrr_gain_thresh * 2

        if DEBUG: print(f"  [LRR_BAF {chrom}] Pre-trim median LRR: {median_lrr_pretrim:.3f}, trim_thresh: gain={trim_gain_thresh:.3f}, loss={trim_loss_thresh:.3f}")

        # Trim left
        for k in range(len(region_snps) - trim_window):
            window = region_snps.iloc[k:k+trim_window]
            window_lrr = window['lrr_filt'].median()
            if state == 'GAIN' and window_lrr > trim_gain_thresh:
                start = window['position'].iloc[0]
                break
            elif state == 'LOSS' and window_lrr < trim_loss_thresh:
                start = window['position'].iloc[0]
                break
        
        # Trim right
        for k in range(len(region_snps) - 1, trim_window - 1, -1):
            window = region_snps.iloc[k-trim_window+1:k+1]
            window_lrr = window['lrr_filt'].median()
            if state == 'GAIN' and window_lrr > trim_gain_thresh:
                end = window['position'].iloc[-1]
                break
            elif state == 'LOSS' and window_lrr < trim_loss_thresh:
                end = window['position'].iloc[-1]
                break
        
        if DEBUG: print(f"  [LRR_BAF {chrom}] After trim: {start/1e6:.1f}-{end/1e6:.1f} Mb")

        region_snps = df[(df['position'] >= start) & (df['position'] <= end)]

        state_bin_fraction = segment['n_state_bins'] / len(region_snps) if len(region_snps) > 0 else 0
        if DEBUG: print(f"  [LRR_BAF {chrom}] State bin fraction: {state_bin_fraction:.2f} ({segment['n_state_bins']}/{len(region_snps)})")
        if state_bin_fraction < 0.15:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: state_bin_fraction {state_bin_fraction:.2f} < 0.15")
            continue
        
        median_lrr = region_snps['lrr_filt'].median()
        if DEBUG: print(f"  [LRR_BAF {chrom}] Post-trim median LRR: {median_lrr:.3f}")
        if pd.isna(median_lrr):
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: median_lrr is NaN")
            continue

        lrr_sd = region_snps['lrr_filt'].std()
        if DEBUG: print(f"  [LRR_BAF {chrom}] LRR sd: {lrr_sd:.3f}")
        if lrr_sd > 0.35:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: LRR too noisy (sd={lrr_sd:.3f} > 0.35)")
            continue

        # BAF confirmation — adaptive het range
        if state == 'LOSS' and median_lrr <= lrr_loss_thresh:
            p_est = min(2 * (1 - 2**median_lrr), 1.0)
            baf_shift = p_est / (2 - p_est) if p_est < 2 else 0.45
            adaptive_lo = max(0.05, true_het_baf_lo - baf_shift)
            adaptive_hi = min(0.95, true_het_baf_hi + baf_shift)
        elif state == 'GAIN' and median_lrr >= lrr_gain_thresh:
            p_est = min(2 * (2**median_lrr - 1), 1.0)
            baf_shift = p_est / (2 + p_est) if p_est > 0 else 0
            adaptive_lo = max(0.05, true_het_baf_lo - baf_shift)
            adaptive_hi = min(0.95, true_het_baf_hi + baf_shift)
        else:
            adaptive_lo = true_het_baf_lo
            adaptive_hi = true_het_baf_hi

        if DEBUG: print(f"  [LRR_BAF {chrom}] Het range: [{adaptive_lo:.2f}, {adaptive_hi:.2f}]")

        true_hets_in_region = region_snps[
            (region_snps['BAF'] >= adaptive_lo) & 
            (region_snps['BAF'] <= adaptive_hi)]
        
        if DEBUG: print(f"  [LRR_BAF {chrom}] True hets in region: {len(true_hets_in_region)}")

        loss_bins = region_snps[region_snps['lrr_filt'] < lrr_loss_thresh]
        loss_bin_median = loss_bins['lrr_filt'].median() if len(loss_bins) > 0 else 0
        if DEBUG: print(f"  [LRR_BAF {chrom}] Loss bin median: {loss_bin_median:.3f} ({len(loss_bins)} loss bins)")

        # ALTERNATIVE PATH: Deep deletion
        if state == 'LOSS' and loss_bin_median < -0.7:
            if DEBUG: print(f"  [LRR_BAF {chrom}] Entering deep deletion path (loss_bin_median={loss_bin_median:.3f})")
            total_snps = len(region_snps)
            n_hets = len(true_hets_in_region)
            het_fraction = n_hets / total_snps if total_snps > 0 else 0
            
            bg_snps = df[~((df['position'] >= start) & (df['position'] <= end))]
            bg_hets = bg_snps[(bg_snps['BAF'] >= true_het_baf_lo) & 
                              (bg_snps['BAF'] <= true_het_baf_hi)]
            bg_het_fraction = len(bg_hets) / len(bg_snps) if len(bg_snps) > 0 else 0.3

            if loss_bin_median < -0.9:
                het_thresh = 0.5
            else:
                het_thresh = 0.25
            
            if DEBUG: print(f"  [LRR_BAF {chrom}] Deep del: het_frac={het_fraction:.3f}, bg_het_frac={bg_het_fraction:.3f}, thresh={het_thresh}")

            if bg_het_fraction > 0 and het_fraction < bg_het_fraction * het_thresh:
                cf_lrr = 2 * (1 - 2**median_lrr)
                cf_lrr = float(np.clip(cf_lrr, 0, 1))
                het_depletion = 1 - (het_fraction / bg_het_fraction) if bg_het_fraction > 0 else 1.0
                cf_baf = min(1.0, het_depletion)
                
                cf_combined = (cf_lrr + cf_baf) / 2
                cf_diff = abs(cf_lrr - cf_baf)
                
                if DEBUG: print(f"  [LRR_BAF {chrom}] Deep del: cf_lrr={cf_lrr:.3f}, cf_baf={cf_baf:.3f}, cf_diff={cf_diff:.3f}")
                
                if cf_diff > 0.3:
                    if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: deep del cf_diff {cf_diff:.3f} > 0.3")
                    continue

                span_mb = (end - start) / 1e6
                if span_mb < 2:
                    if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: deep del span {span_mb:.1f} < 2 Mb")
                    continue
                
                confidence = 'HIGH' if cf_diff < 0.15 else 'MEDIUM'
                
                valid_lrr_count = region_snps['lrr_filt'].notna().sum()
                valid_lrr_fraction = valid_lrr_count / len(region_snps)
                
                if DEBUG: print(f"  [LRR_BAF {chrom}] ✅ Deep deletion PASSED: {start/1e6:.1f}-{end/1e6:.1f} Mb, cf_combined={cf_combined:.3f}")
                
                calls.append({
                    'chromosome': chrom,
                    'start_pos': int(start),
                    'end_pos': int(end),
                    'span_mb': span_mb,
                    'event': state,
                    'confidence': confidence,
                    'cell_fraction_LRR': cf_lrr,
                    'cell_fraction_BAF': cf_baf,
                    'cell_fraction_combined': cf_combined,
                    'median_LRR': median_lrr,
                    'median_AI': np.nan,
                    'baf_support_fraction': het_depletion,
                    'n_bins': len(region_snps),
                    'n_total_snps': len(region_snps),
                    'n_true_hets': n_hets,
                    'valid_lrr_fraction': valid_lrr_fraction,
                    'detector': 'lrr_baf_combined'
                })
                continue
            else:
                if DEBUG: print(f"  [LRR_BAF {chrom}] Deep del path: het depletion insufficient, falling through to normal path")

        # Recalculate span after trimming
        span_mb = (end - start) / 1e6
        if state == 'GAIN' and median_lrr > 0.4:
            effective_min_span = min_span_mb * 0.7
        elif state == 'LOSS' and loss_bin_median < -0.7:
            effective_min_span = min_span_mb * 0.7
        else:
            effective_min_span = min_span_mb
        
        if DEBUG: print(f"  [LRR_BAF {chrom}] Span: {span_mb:.1f} Mb, effective_min_span: {effective_min_span:.1f} Mb")
            
        if span_mb < effective_min_span:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: span {span_mb:.1f} < {effective_min_span:.1f} Mb")
            continue

        # Median LRR must agree with call direction
        if state == 'GAIN' and median_lrr <= lrr_gain_thresh:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: GAIN but median_lrr {median_lrr:.3f} < {lrr_gain_thresh}")
            continue
        if state == 'LOSS' and median_lrr >= lrr_loss_thresh:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: LOSS but median_lrr {median_lrr:.3f} > {lrr_loss_thresh}")
            continue

        if len(true_hets_in_region) < min_hets_in_region:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: true_hets {len(true_hets_in_region)} < {min_hets_in_region}")
            continue
        
        ai_values = true_hets_in_region['ai'].dropna()
        
        if len(ai_values) < min_hets_in_region:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: ai_values {len(ai_values)} < {min_hets_in_region}")
            continue
        
        hets_with_ai = (ai_values > baf_ai_thresh).sum()
        ai_support_fraction = hets_with_ai / len(ai_values)
        median_ai = ai_values.median()
        
        if DEBUG: print(f"  [LRR_BAF {chrom}] AI: median={median_ai:.3f}, support_frac={ai_support_fraction:.3f} ({hets_with_ai}/{len(ai_values)}), ai_min_threshold={ai_min_threshold:.3f}")
        
        if ai_support_fraction >= min_baf_support_fraction:
            if DEBUG: print(f"  [LRR_BAF {chrom}] BAF confirms LRR (support {ai_support_fraction:.2f} >= {min_baf_support_fraction})")

            lrr_strength = abs(median_lrr - chr_baseline_lrr)
            if lrr_strength > 0.20:
                # Strong LRR — accept minimal BAF confirmation
                effective_ai_thresh = ai_min_threshold * 0.4
            elif lrr_strength > 0.08:
                effective_ai_thresh = ai_min_threshold * 0.5
            else:
                effective_ai_thresh = ai_min_threshold

            if median_ai < effective_ai_thresh:
                if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: median_ai {median_ai:.3f} < effective_ai_threshold {effective_ai_thresh:.3f}")
                continue

            if span_mb < 5:
                if lrr_strength > 0.20:
                    small_ai_thresh = 0.04
                elif lrr_strength > 0.08:
                    small_ai_thresh = 0.04
                else:
                    small_ai_thresh = 0.10
                if median_ai < small_ai_thresh: 
                    if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: small span ({span_mb:.1f} Mb) + weak AI ({small_ai_thresh:.3f} < 0.1)")
                    continue
            
            # Cell fraction estimates
            if state == 'GAIN':
                cf_lrr = 2 * (2**median_lrr - 1)
                cf_lrr = float(np.clip(cf_lrr, 0, 1))
                cf_baf = 4 * median_ai / (1 - 2 * median_ai) if median_ai < 0.5 else 1.0
                cf_baf = float(np.clip(cf_baf, 0, 1))
            else:  # LOSS
                cf_lrr = 2 * (1 - 2**median_lrr)
                cf_lrr = float(np.clip(cf_lrr, 0, 1))
                cf_baf = 4 * median_ai / (1 + 2 * median_ai)
                cf_baf = float(np.clip(cf_baf, 0, 1))
            
            cf_combined = (cf_lrr + cf_baf) / 2
            cf_diff = abs(cf_lrr - cf_baf)

            if DEBUG: print(f"  [LRR_BAF {chrom}] Cell fractions: cf_lrr={cf_lrr:.3f}, cf_baf={cf_baf:.3f}, cf_combined={cf_combined:.3f}, cf_diff={cf_diff:.3f}")

            if cf_combined < 0.3:
                if (state == 'GAIN' and 
                    median_lrr > 0.04 and
                    median_ai > 0.04 and
                    ai_support_fraction > 0.8 and
                    span_mb > 5 and
                    cf_diff < 0.10 and
                    len(true_hets_in_region) >= 15):
                    confidence = 'LOW'
                    if DEBUG: print(f"  [LRR_BAF {chrom}] cf_combined < 0.3 but RESCUED (weak GAIN)")
                elif (cf_diff < 0.15 and
                      ai_support_fraction >= 0.8 and
                      span_mb > 3 and
                      len(true_hets_in_region) >= 15):
                    confidence = 'LOW'
                    if DEBUG: print(f"  [LRR_BAF {chrom}] cf_combined < 0.3 but RESCUED (strong agreement: cf_diff={cf_diff:.3f})")
                else:
                    if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: cf_combined {cf_combined:.3f} < 0.3")
                    continue

            if cf_diff > 0.2:
                if lrr_strength > 0.3 and cf_diff < 0.4:
                    confidence = 'MEDIUM'
                    if DEBUG: print(f"  [LRR_BAF {chrom}] cf_diff {cf_diff:.3f} > 0.2 but RESCUED (strong LRR)")
                else:
                    if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: cf_diff {cf_diff:.3f} > 0.2")
                    continue
            if cf_lrr < 0.1 and cf_baf < 0.1:
                if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: both cf_lrr ({cf_lrr:.3f}) and cf_baf ({cf_baf:.3f}) < 0.1")
                continue
            
            if cf_diff < 0.15 and ai_support_fraction > 0.5:
                confidence = 'HIGH'
            elif cf_diff < 0.25:
                confidence = 'MEDIUM'
            else:
                confidence = 'LOW'
            
            valid_lrr_count = region_snps['lrr_filt'].notna().sum()
            valid_lrr_fraction = valid_lrr_count / len(region_snps)
            
            if DEBUG: print(f"  [LRR_BAF {chrom}] ✅ PASSED: {state} {start/1e6:.1f}-{end/1e6:.1f} Mb, cf_combined={cf_combined:.3f}, confidence={confidence}")
            
            calls.append({
                'chromosome': chrom,
                'start_pos': int(start),
                'end_pos': int(end),
                'span_mb': span_mb,
                'event': state,
                'confidence': confidence,
                'cell_fraction_LRR': cf_lrr,
                'cell_fraction_BAF': cf_baf,
                'cell_fraction_combined': cf_combined,
                'median_LRR': median_lrr,
                'median_AI': median_ai,
                'baf_support_fraction': ai_support_fraction,
                'n_bins': len(region_snps),
                'n_total_snps': len(region_snps),
                'n_hets': len(true_hets_in_region),
                'n_true_hets': len(true_hets_in_region),
                'valid_lrr_fraction': valid_lrr_fraction,
                'detector': 'lrr_baf'
            })
        else:
            if DEBUG: print(f"  [LRR_BAF {chrom}] ❌ KILLED: ai_support_fraction {ai_support_fraction:.3f} < {min_baf_support_fraction}")
    
    if DEBUG: print(f"  [LRR_BAF {chrom}] Final calls: {len(calls)}")
    return pd.DataFrame(calls) if calls else pd.DataFrame()

#Detect any event type any any cell fraction via AI sliding windows
def BAF_AI_windows_fast(df, chrom, window_hets, step_hets, support_ai_thresh, min_event_snps=20, centromere_dict = centromere_dict):

    expected_cols = [
        "chromosome", "w_start_pos", "w_end_pos", 
        "w_mean_AI", "w_mean_MAF", "w_mean_LRR", 
        "w_support_frac"
    ]

    # --- 2. DATA SETUP ---
    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)
    # Note: 'h' here includes the noisy 0.02 SNPs because of your upstream filter
    h = d[d["is_het_like"] & d["maf"].notna()].copy().reset_index(drop=True)

    if len(h) < 1: return pd.DataFrame(columns=expected_cols)


# --- 3. SPLIT BY CENTROMERE ---
    if chrom in centromere_dict:
        c_start, c_end = centromere_dict[chrom]['start'], centromere_dict[chrom]['end']
        p_arm = h[h['position'] < c_start].copy()
        q_arm = h[h['position'] > c_end].copy()
        arms_to_process = [p_arm, q_arm]
    else:
        arms_to_process = [h]

    # --- 4. COMPUTE WINDOWS PER ARM ---
    results = []

    for arm_hets in arms_to_process:
        if len(arm_hets) < 10:
            continue

        segment_data = arm_hets.copy()

        # Check if this is a telomeric arm edge
        is_p_arm_start = segment_data.iloc[0]['position'] < 5e6
        is_q_arm_end = segment_data.iloc[-1]['position'] > (d['position'].max() - 5e6)
        
        if is_p_arm_start:
            # Left-aligned: window looks backward, good for telomeric p-arm
            seg_roll_mean = lambda x: x.rolling(window=window_hets, min_periods=3).mean()
        elif is_q_arm_end:
            # Right-aligned: window looks forward, good for telomeric q-arm
            seg_roll_mean = lambda x: x.iloc[::-1].rolling(window=window_hets, min_periods=3).mean().iloc[::-1]
        else:
            seg_roll_mean = lambda x: x.rolling(window=window_hets, center=True, min_periods=3).mean()

        segment_data['w_mean_AI'] = seg_roll_mean(segment_data['ai'])
        segment_data['w_mean_MAF'] = seg_roll_mean(segment_data['maf'])
        segment_data['w_mean_LRR'] = seg_roll_mean(segment_data['lrr_filt'])
        
        segment_data['is_supported'] = (segment_data['ai'] >= support_ai_thresh).astype(float)
        segment_data['w_support_frac'] = seg_roll_mean(segment_data['is_supported'])

        # seg_roll_mean = lambda x: x.rolling(window=window_hets, center=True, min_periods=1).mean()
        # segment_data['w_mean_AI'] = seg_roll_mean(segment_data['ai'])
        # segment_data['w_mean_MAF'] = seg_roll_mean(segment_data['maf'])
        # segment_data['w_mean_LRR'] = seg_roll_mean(segment_data['lrr_filt'])
        
        # segment_data['is_supported'] = (segment_data['ai'] >= support_ai_thresh).astype(float)
        # segment_data['w_support_frac'] = seg_roll_mean(segment_data['is_supported'])

        # Subsample
        stride_idx = segment_data.index[::step_hets]
        first_idx = segment_data.index[:1]
        last_idx = segment_data.index[-1:]
        final_idx = stride_idx.union(first_idx).union(last_idx).sort_values()

        res = segment_data.loc[final_idx].copy()
        res['w_start_pos'] = res['position']
        res['w_end_pos'] = res['position']
        res['chromosome'] = chrom

        results.append(res)

    # --- 6. COMBINE ---
    if not results: return pd.DataFrame(columns=expected_cols)

    final_result = pd.concat(results, ignore_index=True)
    return final_result[expected_cols]

def plot_changepoint_refinement(
    d, rough_start, rough_end, refined_start, refined_end,
    chrom, enter_thresh=0.03, centromere_dict=None,
    search_margin_mb=5.0, truth_start=None, truth_end=None
):
    """
    Diagnostic plot showing the changepoint refinement.
    Useful for validating on simulated data.
    """
    import matplotlib.pyplot as plt

    rough_span = rough_end - rough_start
    search_margin = max(search_margin_mb * 1e6, rough_span * 0.3)
    search_start = max(0, rough_start - search_margin)
    search_end = rough_end + search_margin

    mask = (
        (d["position"] >= search_start)
        & (d["position"] <= search_end)
        & d["is_het_like"]
        & d["maf"].notna()
    )
    hets = d[mask].sort_values("position")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    # Panel 1: Per-SNP AI
    ax1.scatter(hets["position"] / 1e6, hets["ai"], s=8, c="grey", alpha=0.6)
    ax1.axhline(enter_thresh, ls="--", c="grey", lw=0.8)

    # Shade regions
    ax1.axvspan(rough_start / 1e6, rough_end / 1e6, alpha=0.15, color="red", label="Rough (state machine)")
    ax1.axvspan(refined_start / 1e6, refined_end / 1e6, alpha=0.20, color="blue", label="Refined (changepoint)")

    if truth_start is not None and truth_end is not None:
        ax1.axvline(truth_start / 1e6, c="green", ls="-", lw=2, label="Truth")
        ax1.axvline(truth_end / 1e6, c="green", ls="-", lw=2)

    ax1.set_ylabel("AI (per het SNP)")
    ax1.legend(fontsize=8)
    ax1.set_title(f"{chrom}: Boundary refinement comparison")

    # Panel 2: LRR
    all_snps = d[(d["position"] >= search_start) & (d["position"] <= search_end)]
    ax2.scatter(all_snps["position"] / 1e6, all_snps["lrr_filt"], s=6, c="grey", alpha=0.4)
    ax2.axhline(0, ls="--", c="grey", lw=0.8)
    ax2.axvspan(rough_start / 1e6, rough_end / 1e6, alpha=0.15, color="red")
    ax2.axvspan(refined_start / 1e6, refined_end / 1e6, alpha=0.20, color="blue")

    if truth_start is not None and truth_end is not None:
        ax2.axvline(truth_start / 1e6, c="green", ls="-", lw=2)
        ax2.axvline(truth_end / 1e6, c="green", ls="-", lw=2)

    ax2.set_ylabel("LRR (filtered)")
    ax2.set_xlabel("Position (Mb)")

    plt.tight_layout()
    return fig

def detect_BAF_AI(df, chrom, window_hets, step_hets, enter_thresh, 
                                  exit_thresh, gain_lrr, loss_lrr, neutral_lrr, min_span_mb, 
                                  ai_point_thresh, gap_mb, het_lo, het_hi, min_consecutive_hets, min_spanning_hets, min_supporting_fraction, support_ai_thresh, centromere_dict, 
                                  lrr_sigma_fallback=0.08, lrr_sigma_min_points=50, chromosome_sizes = chromosome_sizes, debug_chrom=None):

    # DEBUG = (chrom==debug_chrom) # <<< SET TO FALSE TO DISABLE DEBUG OUTPUT
    DEBUG = False # <<< SET TO FALSE TO DISABLE DEBUG OUTPUT
    # DEBUG = True


    # 1. Per chromosome view (Full Data)
    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)

    # 2. FAST SCAN: Generate windows
    support_ai_thresh = enter_thresh  # Use consistent threshold
    wdf = BAF_AI_windows_fast(df, chrom, window_hets, step_hets, support_ai_thresh=support_ai_thresh, min_event_snps=20, centromere_dict=centromere_dict)
    if wdf.empty:
        if DEBUG: print(f"  [{chrom}] No windows generated - empty wdf")
        return pd.DataFrame(), wdf

    wdf = wdf.copy()
    
    # Extract tracks for state machine
    ai_track = wdf["w_mean_AI"].to_numpy(dtype=float)
    support_track = wdf["w_support_frac"].to_numpy(dtype=float)

    # 3. State Machine (Identifies Intervals)
    intervals = []
    in_call = False
    L = None
    exit_count = 0
    EXIT_PERSISTENCE = 3

    for i, val in enumerate(ai_track):
        ok_support = (support_track[i] >= min_supporting_fraction)
        
        if in_call and chrom in centromere_dict:
            c_start = centromere_dict[chrom]['start']
            call_start_pos = wdf.iloc[L]['w_start_pos']
            current_pos = wdf.iloc[i]['w_start_pos']
            if call_start_pos < c_start and current_pos > c_start:
                intervals.append((L, i - 1)); in_call = False; L = None

        if not in_call:
            if (val >= enter_thresh) and ok_support:
                in_call = True; L = i
        else:
            if val < exit_thresh: 
                exit_count += 1
                if exit_count >= EXIT_PERSISTENCE:
                    intervals.append((L, i - EXIT_PERSISTENCE))
                    in_call = False; L = None; exit_count = 0
            else:
                exit_count = 0

    if in_call: intervals.append((L, len(ai_track) - 1))

    if not intervals: 
        if DEBUG: print(f"  [{chrom}] State machine found NO intervals")
        return pd.DataFrame(), wdf

    if DEBUG:
        for idx, (Lw, Rw) in enumerate(intervals):
            s = int(wdf["w_start_pos"].iloc[Lw])
            e = int(wdf["w_end_pos"].iloc[Rw])
            print(f"  [{chrom}] Interval {idx}: {s/1e6:.1f}-{e/1e6:.1f} Mb")

    # 4. Interval Refinement
    bp_intervals = []
    for Lw, Rw in intervals:
        s = int(wdf["w_start_pos"].iloc[Lw])
        e = int(wdf["w_end_pos"].iloc[Rw])
        bp_intervals.append((s, e))

    # 5. Sigma Estimation
    lrr_sigma = lrr_sigma_outside_intervals(d_chrom=d, bp_intervals=bp_intervals, min_points=lrr_sigma_min_points)
    if np.isnan(lrr_sigma) or lrr_sigma <= 0: lrr_sigma = float(lrr_sigma_fallback)

    # 6. Process Each Candidate
    calls = []
    
    for interval_idx, (Lw, Rw) in enumerate(intervals):
        start_pos = int(wdf["w_start_pos"].iloc[Lw])
        end_pos   = int(wdf["w_end_pos"].iloc[Rw])

        if DEBUG: print(f"\n  [{chrom}] === Processing interval {interval_idx}: "
                        f"{start_pos/1e6:.1f}-{end_pos/1e6:.1f} Mb ===")

        # --- v7 CHANGEPOINT BOUNDARY REFINEMENT ---
        rough_seg = d[(d["position"] >= start_pos) & (d["position"] <= end_pos)]
        rough_lrr = rough_seg["lrr_filt"].dropna()
        rough_median_lrr = float(np.nanmedian(rough_lrr)) if len(rough_lrr) > 5 else 0.0
        inner_start = start_pos + (end_pos - start_pos) * 0.25
        inner_end = end_pos - (end_pos - start_pos) * 0.25
        inner_seg = d[(d["position"] >= inner_start) & (d["position"] <= inner_end)]
        inner_lrr = inner_seg["lrr_filt"].dropna()
        inner_median_lrr = float(np.nanmedian(inner_lrr)) if len(inner_lrr) > 5 else 0.0
        use_lrr_for_cpd = abs(inner_median_lrr) > 0.06

        start_pos, end_pos = refine_boundaries_changepoint(
            d,
            rough_start=start_pos,
            rough_end=end_pos,
            enter_thresh=enter_thresh,
            search_margin_mb=5.0,
            min_segment_hets=5,
            pen_multiplier=3.0,
            use_lrr=use_lrr_for_cpd,
            lrr_weight=0.6,
            centromere_dict=centromere_dict,
            chrom=chrom,
            DEBUG=DEBUG,
        )

        if DEBUG: print(f"  [{chrom}] After CPD refinement: "
                        f"{start_pos/1e6:.1f}-{end_pos/1e6:.1f} Mb")

        # --- LRR BOUNDARY SHARPENING (for low-CF GAIN/LOSS only) ---
        if abs(inner_median_lrr) > 0.20:
            bg_lrr = d[(d["position"] < start_pos) | (d["position"] > end_pos)]["lrr_filt"].dropna()
            bg_median_lrr = float(np.nanmedian(bg_lrr)) if len(bg_lrr) > 20 else 0.0
            lrr_contrast_sharp = abs(inner_median_lrr - bg_median_lrr)
            
            inner_hets = d[(d["position"] >= inner_start) & (d["position"] <= inner_end) 
                          & d["is_het_like"] & d["maf"].notna()]
            inner_median_ai = float(np.nanmedian(inner_hets["ai"])) if len(inner_hets) > 5 else 0.0
            
            if lrr_contrast_sharp > 0.12 and inner_median_ai < 0.10:
                lrr_start, lrr_end = refine_boundaries_by_lrr(
                    d, start_pos, end_pos,
                    event_type="LOSS" if inner_median_lrr < 0 else "GAIN",
                    search_margin_mb=3.0, window_snps=10
                )
                if lrr_start > start_pos:
                    start_pos = lrr_start
                if lrr_end < end_pos:
                    end_pos = lrr_end
                if DEBUG: print(f"  [{chrom}] After LRR sharpening: "
                                f"{start_pos/1e6:.1f}-{end_pos/1e6:.1f} Mb")

        # === CHECK FOR LARGE GAPS ===
        seg_all = d[(d["position"] >= start_pos) & (d["position"] <= end_pos)]
        het_seg = seg_all[seg_all["is_het_like"] & seg_all["maf"].notna()].copy()

        if len(het_seg) >= 2:
            het_seg_sorted = het_seg.sort_values('position')
            het_seg_sorted['gap_to_next_mb'] = het_seg_sorted['position'].diff() / 1e6
            
            large_gaps = het_seg_sorted[het_seg_sorted['gap_to_next_mb'] > gap_mb]
            
            if len(large_gaps) > 0:
                gap_position = large_gaps.iloc[0]['position']
                last_snp_before_gap = het_seg_sorted[het_seg_sorted['position'] < gap_position].iloc[-1]['position']
                
                if DEBUG: print(f"  [{chrom}] Large gap at {gap_position/1e6:.1f} Mb, trimming end to {last_snp_before_gap/1e6:.1f} Mb")
                
                end_pos = int(last_snp_before_gap)
                span_mb = (end_pos - start_pos) / 1e6
                
                if span_mb < min_span_mb:
                    if DEBUG: print(f"  [{chrom}] ❌ KILLED at GAP CHECK: span after gap trim {span_mb:.1f} < {min_span_mb}")
                    continue

        span_mb = (end_pos - start_pos) / 1e6
        
        if span_mb < min_span_mb:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at MIN SPAN: {span_mb:.1f} < {min_span_mb}")
            continue

        # --- RE-SLICE & QC ---
        seg_all = d[(d["position"] >= start_pos) & (d["position"] <= end_pos)]
        het_seg = seg_all[seg_all["is_het_like"] & seg_all["maf"].notna()].copy()
        
        if len(het_seg) == 0:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED: 0 hets after re-slice")
            continue

        support_ai_thresh = enter_thresh
        het_seg["supports_event"] = het_seg["ai"] >= support_ai_thresh
        max_run = max_consecutive_true(het_seg["supports_event"])
        if np.isnan(max_run): max_run = 0

        if span_mb > 10.0:
            required_consecutive = 3
        else:
            required_consecutive = min_consecutive_hets
        
        if max_run < required_consecutive:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at CONSECUTIVE (2nd): {max_run} < {required_consecutive}")
            continue

        supp_pos = het_seg.loc[het_seg["ai"] >= support_ai_thresh, "position"].to_numpy()

        min_support_for_span = max(8, span_mb * 1.0)
        if len(supp_pos) < min_support_for_span:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at DENSITY CHECK: {len(supp_pos)} supporting hets < {min_support_for_span:.0f} required (span={span_mb:.1f} Mb)")
            continue

        support_span_mb = (supp_pos.max() - supp_pos.min()) / 1e6
        
        if support_span_mb < min_spanning_hets:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at SUPPORT SPAN (2nd): {support_span_mb:.1f} < {min_spanning_hets}")
            continue

        # --- CALCULATE FINAL STATS ---
        het_lrr_vals = seg_all[seg_all['is_het_like']]['lrr_filt'].dropna()
        all_lrr_vals = seg_all['lrr_filt'].dropna()

        mean_maf = float(np.nanmean(het_seg["maf"]))
        mean_ai  = float(np.nanmean(het_seg["ai"]))
        mean_lrr = float(np.nanmean(all_lrr_vals))
        median_lrr = float(np.nanmedian(all_lrr_vals))
        het_median_lrr = float(np.nanmedian(het_lrr_vals)) if len(het_lrr_vals) >= 5 else median_lrr

        if DEBUG: print(f"  [{chrom}] Stats: mean_ai={mean_ai:.3f}, median_lrr={median_lrr:.3f}, het_median_lrr={het_median_lrr:.3f}, n_hets={len(het_seg)}, span={span_mb:.1f} Mb")

        # --- BACKGROUND LRR CORRECTION ---
        is_background = (d['position'] < start_pos) | (d['position'] > end_pos)
        background_lrr_data = d.loc[is_background, 'lrr_filt'].dropna()
        
        chr_size_mb_bg = chromosome_sizes.get(chrom, 250e6) / 1e6
        event_fraction_bg = span_mb / chr_size_mb_bg
        bg_is_contaminated = False
        
        if len(background_lrr_data) > 50:
            background_median_lrr_for_class = np.nanmedian(background_lrr_data)
            if abs(background_median_lrr_for_class) > 0.20 or event_fraction_bg > 0.40:
                bg_is_contaminated = True
        else:
            background_median_lrr_for_class = 0.0
            if event_fraction_bg > 0.40:
                bg_is_contaminated = True
        
        if bg_is_contaminated:
            genome_bg_lrr = df[(df["chromosome"] != chrom)]['lrr_filt'].dropna()
            if len(genome_bg_lrr) > 100:
                genome_baseline_for_class = np.nanmedian(genome_bg_lrr)
                if DEBUG: print(f"  [{chrom}] ⚠ Background contaminated (bg_median={background_median_lrr_for_class:.3f}, event_frac={event_fraction_bg:.2f})")
                if DEBUG: print(f"  [{chrom}] → Using genome-wide baseline for classification: {genome_baseline_for_class:.4f}")
                background_median_lrr_for_class = genome_baseline_for_class
        
        het_median_lrr_corrected = het_median_lrr - background_median_lrr_for_class
        median_lrr_corrected = median_lrr - background_median_lrr_for_class
        
        if DEBUG: print(f"  [{chrom}] LRR correction: bg_median={background_median_lrr_for_class:.3f}, het_median_lrr_corrected={het_median_lrr_corrected:.3f}, bg_contaminated={bg_is_contaminated}")



        # --- CLASSIFICATION ---
        if abs(het_median_lrr_corrected) <= neutral_lrr: event = "CN-LOH"
        elif het_median_lrr_corrected >= gain_lrr: event = "GAIN"
        elif het_median_lrr_corrected <= loss_lrr: event = "LOSS"
        else: event = "AI_event_unclear"

        if DEBUG: print(f"  [{chrom}] Initial classification: {event} (het_median_lrr_corrected={het_median_lrr_corrected:.3f}, gain_lrr={gain_lrr}, loss_lrr={loss_lrr}, neutral_lrr={neutral_lrr})")

        min_per_side = max(3, min(8, len(het_seg) // 3))
        mu1, mu2, delta = baf_peaks_mu1_mu2(het_seg, min_per_side, het_lo, het_hi)

        lrr_confirms_event = (
            (event == "GAIN" and het_median_lrr_corrected >= gain_lrr) or
            (event == "LOSS" and het_median_lrr_corrected <= loss_lrr)
        )

        # Moderate LRR: clearly in the right direction but below strict threshold
        lrr_directional = (
            not lrr_confirms_event and (
                (het_median_lrr_corrected > 0.06 and event in ("GAIN", "AI_event_unclear")) or
                (het_median_lrr_corrected < -0.06 and event in ("LOSS", "AI_event_unclear"))
            )
        )
        
        cnloh_consistent = (
            event == "CN-LOH" and abs(het_median_lrr_corrected) <= 0.06
        )

        if len(het_seg) > 100:
            min_delta_floor = 0.05
        elif len(het_seg) > 50:
            min_delta_floor = 0.07
        else:
            min_delta_floor = 0.10

        if lrr_confirms_event:
            min_delta_floor *= 0.6
            scaling = 0.40
        elif lrr_directional:
            min_delta_floor *= 0.8
            scaling = 0.55
        elif cnloh_consistent:
            min_delta_floor *= 0.8
            scaling = 0.55
        else:
            scaling = 0.70

        min_delta_for_hets = max(min_delta_floor, scaling / np.sqrt(len(het_seg)))

        if DEBUG: print(f"  [{chrom}] Delta check: delta={(f'{delta:.4f}' if not np.isnan(delta) else 'NaN')}, min_delta={min_delta_for_hets:.4f}, lrr_confirms={lrr_confirms_event}, mu1={(f'{mu1:.3f}' if not np.isnan(mu1) else 'NaN')}, mu2={(f'{mu2:.3f}' if not np.isnan(mu2) else 'NaN')}") 
        
        if not np.isnan(delta) and delta < min_delta_for_hets:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at DELTA CHECK: {delta:.4f} < {min_delta_for_hets:.4f}")
            continue

        cls = classify_event_by_baf_lrr(delta, het_median_lrr_corrected, lrr_sigma)
        
        event  = cls["event"]
        cf_baf = cls["p_baf"]

        # Minimum cell fraction
        min_cf_baf = 0.12
        if np.isfinite(cf_baf) and cf_baf < min_cf_baf:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at MIN CF: cf_baf={cf_baf:.3f} < {min_cf_baf}")
            continue

        # CN-LOH requires stronger evidence
        # CN-LOH: adaptive het count check — scales with span and delta strength
        # v14: replaced hard density override (het_density<4 → 18) with adaptive formula
        # to handle low-density SNP arrays (~1-3 hets/Mb) where 7-10 Mb events yield 10-17 hets
        if event == "CN-LOH":
            het_density = len(het_seg) / span_mb if span_mb > 0 else 0

            # Adaptive floor: scale minimum hets by span, modulated by delta confidence
            if not np.isnan(delta) and delta > 0.30:
                # Strong BAF split — high confidence, lower requirement
                min_hets_cnloh = max(8, int(np.ceil(span_mb * 1.0)))
            elif not np.isnan(delta) and delta > 0.15:
                # Moderate BAF split — standard requirement
                min_hets_cnloh = max(10, int(np.ceil(span_mb * 1.2)))
            else:
                # Weak/missing BAF split — require more evidence
                min_hets_cnloh = max(12, int(np.ceil(span_mb * 1.5)))

            # Small events (<5 Mb) need extra caution — tighter absolute floor
            if span_mb < 5:
                min_hets_cnloh = max(min_hets_cnloh, 15)

            if len(het_seg) < min_hets_cnloh:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at CN-LOH HET COUNT: {len(het_seg)} < {min_hets_cnloh} (adaptive, density={het_density:.1f} hets/Mb, span={span_mb:.1f} Mb, delta={delta:.3f})")
                continue


        if event == "CN-LOH" and abs(het_median_lrr_corrected) < 0.06:
            bg_hets = d[~((d['position'] >= start_pos) & (d['position'] <= end_pos))]
            bg_hets = bg_hets[bg_hets['is_het_like'] & bg_hets['maf'].notna()]
            if len(bg_hets) > 20:
                bg_ai_median = bg_hets['ai'].median()
                bg_ai_mad = np.median(np.abs(bg_hets['ai'] - bg_ai_median))
                bg_ai_sigma = 1.4826 * bg_ai_mad
                event_ai_median = het_seg['ai'].median()
                ai_snr = (event_ai_median - bg_ai_median) / bg_ai_sigma if bg_ai_sigma > 0 else 0
                if DEBUG: print(f"  [{chrom}] CN-LOH noise check: event_ai={event_ai_median:.3f}, bg_ai={bg_ai_median:.3f}, bg_sigma={bg_ai_sigma:.3f}, SNR={ai_snr:.1f}")
                if ai_snr < 3.0:
                    if DEBUG: print(f"  [{chrom}] ❌ KILLED at CN-LOH SNR: {ai_snr:.1f} < 3.0")
                    continue

            hets_above = (het_seg['BAF'] > 0.5).sum()
            hets_below = (het_seg['BAF'] < 0.5).sum()
            minor_side = min(hets_above, hets_below)
            minor_fraction = minor_side / len(het_seg) if len(het_seg) > 0 else 0
            
            if DEBUG: print(f"  [{chrom}] CN-LOH symmetry: above={hets_above}, below={hets_below}, minor_frac={minor_fraction:.2f}")
            
            if len(het_seg) >= 50 and minor_fraction < 0.15:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at CN-LOH SYMMETRY: minor_fraction {minor_fraction:.2f} < 0.20")
                continue
            elif len(het_seg) >= 30 and minor_fraction < 0.08:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at CN-LOH SYMMETRY: minor_fraction {minor_fraction:.2f} < 0.10")
                continue
            elif len(het_seg) >= 20 and minor_side == 0:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at CN-LOH SYMMETRY: completely one-sided")
                continue

            region_lrr_count = seg_all['lrr_filt'].notna().sum()
            if DEBUG: print(f"  [{chrom}] CN-LOH LRR probes: {region_lrr_count}")
            
            if region_lrr_count < 10:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at CN-LOH LRR COUNT: {region_lrr_count} < 10")
                continue

        cf_baf_str = f"{cf_baf:.3f}" if np.isfinite(cf_baf) else "NaN"
        if DEBUG: print(f"  [{chrom}] classify_event_by_baf_lrr: event={event}, cf_baf={cf_baf_str}, lrr_fit_score={cls.get('lrr_fit_score', 'N/A')}")

        if event in ('GAIN', 'LOSS'):
            lrr_vals = seg_all['lrr_filt'].dropna()
            if len(lrr_vals) >= 10:
                lrr_check_baseline = background_median_lrr_for_class
                if event == 'GAIN':
                    agreeing_frac = (lrr_vals > lrr_check_baseline).sum() / len(lrr_vals)
                else:
                    agreeing_frac = (lrr_vals < lrr_check_baseline).sum() / len(lrr_vals)

                if DEBUG: print(f"  [{chrom}] LRR consistency: {agreeing_frac:.3f} (need >= 0.75)")

                if agreeing_frac < 0.75:
                    if abs(median_lrr) > 0.15 and agreeing_frac > 0.70:
                        confidence = "MEDIUM"
                        if DEBUG: print(f"  [{chrom}] LRR consistency {agreeing_frac:.3f} borderline but strong LRR ({median_lrr:.3f}), downgraded to MEDIUM")
                    else:
                        if DEBUG: print(f"  [{chrom}] ❌ KILLED at LRR CONSISTENCY: {agreeing_frac:.3f} < 0.75")
                        continue

        if np.isnan(delta) or event == "AI_event_unclear":
            if DEBUG: print(f"  [{chrom}] ❌ KILLED: delta NaN or AI_event_unclear")
            continue

        if event == "CN-LOH":
            if not np.isfinite(cf_baf):
                if DEBUG: print(f"  [{chrom}] ❌ KILLED: CN-LOH with non-finite cf_baf")
                continue

        cf_lrr = cls["p_lrr"]
        if event == "CN-LOH":
            cf_ai = 2 * mean_ai
        elif event == "GAIN":
            cf_ai = 4 * mean_ai / (1 - 2 * mean_ai) if mean_ai < 0.5 else 1.0
        elif event == "LOSS":
            cf_ai = 4 * mean_ai / (1 + 2 * mean_ai)
        else:
            cf_ai = 2 * mean_ai
        cf_ai = float(np.clip(cf_ai, 0, 1))
        
        # Flags & CONFIDENCE CALLS
        flag_small_delta = (not np.isnan(delta)) and (delta < 0.05)

        min_het_density = 1.0
        flag_low_hets = (len(het_seg) / span_mb) < min_het_density

        flag_lrr_baf_mismatch_cls = bool(cls.get("flag_lrr_baf_mismatch", False))

        flag_lrr_baf_mismatch = False
        if event in ("GAIN", "LOSS") and (not np.isnan(cf_lrr)) and (not np.isnan(cf_baf)):
            cf_mismatch = abs(cf_lrr - cf_baf)
            # Tighter tolerance when AI signal is weak — a genuine LOSS/GAIN at
            # the CF that LRR claims should produce proportionally stronger AI.
            # Weak AI + strong LRR disagreement = noisy LRR dip, not a real event.
            cf_mismatch_thresh = 0.15 if mean_ai < 0.06 else 0.30
            if cf_mismatch > cf_mismatch_thresh:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at CF MISMATCH: |cf_lrr - cf_baf| = {cf_mismatch:.3f} > {cf_mismatch_thresh} (mean_ai={mean_ai:.3f})")
                continue

        # GAIN/LOSS require LRR confirmation (pre-override check)
        if event in ("GAIN", "LOSS") and not np.isfinite(cf_lrr):
            if DEBUG: print(f"  [{chrom}] ❌ KILLED: {event} with no LRR cell fraction (cf_lrr=NaN)")
            continue

        flag_neutral_lrr_large_ai = (np.isfinite(median_lrr) and (abs(median_lrr) <= neutral_lrr) and (not np.isnan(delta)) and (delta > 0.20))

        # CONFIDENCE ASSIGNMENT - SIZE AWARE
        if span_mb > 15.0:
            if flag_small_delta or flag_low_hets: 
                confidence = "LOW"
            elif flag_lrr_baf_mismatch or flag_lrr_baf_mismatch_cls:
                confidence = "MEDIUM"
            else:
                confidence = "HIGH"
        if span_mb <= 15.0:
            if flag_small_delta:
                confidence = "LOW"
            elif flag_low_hets or flag_lrr_baf_mismatch or flag_lrr_baf_mismatch_cls:
                confidence = "MEDIUM"
            else:
                confidence = "HIGH"

        if DEBUG: print(f"  [{chrom}] Confidence: {confidence}, flags: small_delta={flag_small_delta}, low_hets={flag_low_hets}, lrr_baf_mm={flag_lrr_baf_mismatch}")

        # --- REVISED DESTROYER FILTERS ---
        final_supporting_snps = het_seg[het_seg["ai"] >= support_ai_thresh]
        
        # Filter 0: Window Consistency Check
        interval_windows = wdf.iloc[Lw:Rw+1]
        n_windows = len(interval_windows)
        n_passing_windows = ((interval_windows["w_mean_AI"] >= enter_thresh) & 
                             (interval_windows["w_support_frac"] >= min_supporting_fraction)).sum()
        window_pass_rate = n_passing_windows / n_windows if n_windows > 0 else 0
        
        if DEBUG: print(f"  [{chrom}] Filter 0 (Window consistency): {n_passing_windows}/{n_windows} = {window_pass_rate:.2%}")

        if span_mb < 10.0:
            if window_pass_rate < 0.70:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER 0: window_pass_rate {window_pass_rate:.2%} < 70% (small event)")
                continue
        else:
            if window_pass_rate < 0.50:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER 0: window_pass_rate {window_pass_rate:.2%} < 50% (large event)")
                continue
        
        # Filter A: LRR Noise
        lrr_sd = seg_all["lrr_filt"].std()
        
        if event == "CN-LOH":
            max_lrr_sd = 0.40
        else:
            max_lrr_sd = 0.35
        
        if DEBUG: print(f"  [{chrom}] Filter A (LRR noise): lrr_sd={lrr_sd:.3f}, max={max_lrr_sd}")

        if lrr_sd > max_lrr_sd:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER A: lrr_sd {lrr_sd:.3f} > {max_lrr_sd}")
            continue

        # Filter B: Weak Waves
        if len(final_supporting_snps) > 0: 
            current_ai_strength = final_supporting_snps["ai"].mean()
        else: 
            current_ai_strength = 0

        median_ai_event = het_seg['ai'].median()

        if event == "CN-LOH":
            if span_mb > 20.0: min_ai_strength = 0.04; min_median_ai = 0.035
            elif span_mb > 10.0: min_ai_strength = 0.04; min_median_ai = 0.035
            elif span_mb > 5.0: min_ai_strength = 0.04; min_median_ai = 0.035
            else: min_ai_strength = 0.08; min_median_ai = 0.06
        else:  # GAIN/LOSS
            if span_mb > 20.0: min_ai_strength = 0.035; min_median_ai = 0.03
            elif span_mb > 10.0: min_ai_strength = 0.04; min_median_ai = 0.035
            else: min_ai_strength = 0.04; min_median_ai = 0.035

        if lrr_sd < 0.001 and len(seg_all) > 5:
            if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER A: lrr_sd suspiciously low ({lrr_sd:.4f})")
            continue

        if DEBUG: print(f"  [{chrom}] Filter B: ai_strength={current_ai_strength:.3f} (min={min_ai_strength}), median_ai={median_ai_event:.3f} (min={min_median_ai})")

        if current_ai_strength < min_ai_strength:
            if event in ("GAIN", "LOSS") and abs(median_lrr) > 0.15 and current_ai_strength > min_ai_strength * 0.8:
                confidence = "MEDIUM"
                if DEBUG: print(f"  [{chrom}] Filter B: ai_strength {current_ai_strength:.3f} borderline but strong LRR ({median_lrr:.3f}), downgraded to MEDIUM")
            elif event == "CN-LOH" and not np.isnan(delta) and delta > 0.08 and current_ai_strength > min_ai_strength * 0.75:
                confidence = "MEDIUM"
                if DEBUG: print(f"  [{chrom}] Filter B: CN-LOH rescue — ai_strength {current_ai_strength:.3f} borderline but delta={delta:.3f}, downgraded to MEDIUM")
            else:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER B: ai_strength {current_ai_strength:.3f} < {min_ai_strength}")
                continue
        if median_ai_event < min_median_ai:
            if event in ("GAIN", "LOSS") and abs(median_lrr) > 0.12 and median_ai_event > min_median_ai * 0.7:
                confidence = "MEDIUM"  # rescue
            elif event == "CN-LOH" and not np.isnan(delta) and delta > 0.08 and median_ai_event > min_median_ai * 0.70:
                confidence = "MEDIUM"
                if DEBUG: print(f"  [{chrom}] Filter B (median): CN-LOH rescue — median_ai {median_ai_event:.3f} borderline but delta={delta:.3f}")
            else:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER B (median): median_ai {median_ai_event:.3f} < {min_median_ai}")
                continue

        # Filter C: Neighbor Contrast
        lrr_contrast = None
        is_background = (d['position'] < start_pos) | (d['position'] > end_pos)
        background_data = d.loc[is_background, 'lrr_filt']

        # === BLOCK 1: Forward override — CN-LOH → GAIN/LOSS if CFs agree ===
        lrr_override_fired = False
        if event == "CN-LOH" and len(background_data) > 50:
            bg_med = np.nanmedian(background_data)
            lrr_contrast_check = median_lrr - bg_med

            if abs(median_lrr) > 0.3:
                cf_agree_thresh = 0.40
            else:
                cf_agree_thresh = 0.20
            
            if lrr_contrast_check <= -0.12 or median_lrr <= -0.12:
                cf_lrr_check = 2 * (1 - 2**median_lrr)
                cf_lrr_check = float(np.clip(cf_lrr_check, 0, 1))
                cf_baf_loss = 4 * mean_ai / (1 + 2 * mean_ai)
                cf_baf_loss = float(np.clip(cf_baf_loss, 0, 1))
                if DEBUG: print(f"  [{chrom}] LRR override check LOSS: cf_lrr={cf_lrr_check:.3f}, cf_baf_loss={cf_baf_loss:.3f}, diff={abs(cf_lrr_check - cf_baf_loss):.3f}")
                if abs(cf_lrr_check - cf_baf_loss) < cf_agree_thresh:
                    event = "LOSS"
                    cf_baf = cf_baf_loss
                    lrr_override_fired = True
                    if DEBUG: print(f"  [{chrom}] LRR override: CN-LOH → LOSS")
                else:
                    if DEBUG: print(f"  [{chrom}] LRR override BLOCKED: CF mismatch — keeping CN-LOH")
                    
            elif lrr_contrast_check >= 0.12 or median_lrr >= 0.12:
                cf_lrr_check = 2 * (2**median_lrr - 1)
                cf_lrr_check = float(np.clip(cf_lrr_check, 0, 1))
                cf_baf_gain = 4 * mean_ai / (1 - 2 * mean_ai) if mean_ai < 0.5 else 1.0
                cf_baf_gain = float(np.clip(cf_baf_gain, 0, 1))
                if DEBUG: print(f"  [{chrom}] LRR override check GAIN: cf_lrr={cf_lrr_check:.3f}, cf_baf_gain={cf_baf_gain:.3f}, diff={abs(cf_lrr_check - cf_baf_gain):.3f}")
                if abs(cf_lrr_check - cf_baf_gain) < cf_agree_thresh:
                    event = "GAIN"
                    cf_baf = cf_baf_gain
                    lrr_override_fired = True
                    if DEBUG: print(f"  [{chrom}] LRR override: CN-LOH → GAIN")
                else:
                    if DEBUG: print(f"  [{chrom}] LRR override BLOCKED: CF mismatch — keeping CN-LOH")

        elif event == "CN-LOH" and len(background_data) <= 50:
            if median_lrr <= -0.15:
                cf_lrr_check = 2 * (1 - 2**median_lrr)
                cf_lrr_check = float(np.clip(cf_lrr_check, 0, 1))
                cf_baf_loss = 4 * mean_ai / (1 + 2 * mean_ai)
                cf_baf_loss = float(np.clip(cf_baf_loss, 0, 1))
                if abs(cf_lrr_check - cf_baf_loss) < 0.20:
                    event = "LOSS"
                    cf_baf = cf_baf_loss
                    lrr_override_fired = True
                    if DEBUG: print(f"  [{chrom}] Fallback override: CN-LOH → LOSS (raw median_lrr={median_lrr:.3f})")
            elif median_lrr >= 0.15:
                cf_lrr_check = 2 * (2**median_lrr - 1)
                cf_lrr_check = float(np.clip(cf_lrr_check, 0, 1))
                cf_baf_gain = 4 * mean_ai / (1 - 2 * mean_ai) if mean_ai < 0.5 else 1.0
                cf_baf_gain = float(np.clip(cf_baf_gain, 0, 1))
                if abs(cf_lrr_check - cf_baf_gain) < 0.20:
                    event = "GAIN"
                    cf_baf = cf_baf_gain
                    lrr_override_fired = True
                    if DEBUG: print(f"  [{chrom}] Fallback override: CN-LOH → GAIN (raw median_lrr={median_lrr:.3f})")

        # === BLOCK 2: Reverse override — GAIN/LOSS → CN-LOH if CFs disagree ===
        if event in ("GAIN", "LOSS") and not lrr_override_fired:
            if abs(median_lrr) > 0.15:
                if DEBUG: print(f"  [{chrom}] Reverse override SKIPPED: |median_lrr|={abs(median_lrr):.3f} > 0.15, keeping {event}")
            else:
                if event == "GAIN":
                    cf_lrr_check = 2 * (2**median_lrr - 1)
                    cf_baf_check = 4 * mean_ai / (1 - 2 * mean_ai) if mean_ai < 0.5 else 1.0
                else:
                    cf_lrr_check = 2 * (1 - 2**median_lrr)
                    cf_baf_check = 4 * mean_ai / (1 + 2 * mean_ai)
                cf_lrr_check = float(np.clip(cf_lrr_check, 0, 1))
                cf_baf_check = float(np.clip(cf_baf_check, 0, 1))
                cf_diff_check = abs(cf_lrr_check - cf_baf_check)

                if cf_diff_check > 0.25:
                    if DEBUG: print(f"  [{chrom}] Reverse override: {event} → CN-LOH (cf_lrr={cf_lrr_check:.3f}, cf_baf={cf_baf_check:.3f}, diff={cf_diff_check:.3f})")
                    event = "CN-LOH"
                    cf_baf = min(2 * mean_ai, 1.0)
                    cf_baf = float(np.clip(cf_baf, 0, 1))

        # === BLOCK 3: Post-override LRR boundary trim for reclassified GAIN/LOSS ===
        if lrr_override_fired and event in ("LOSS", "GAIN"):
            seg_for_refine = d[(d["position"] >= start_pos) & (d["position"] <= end_pos)].copy()
            refine_window = max(5, len(seg_for_refine) // 20)
            for k in range(len(seg_for_refine) - refine_window):
                window = seg_for_refine.iloc[k:k+refine_window]
                w_lrr = window['lrr_filt'].median()
                if event == "LOSS" and w_lrr < -0.10:
                    start_pos = int(window['position'].iloc[0])
                    break
                elif event == "GAIN" and w_lrr > 0.10:
                    start_pos = int(window['position'].iloc[0])
                    break
            for k in range(len(seg_for_refine) - 1, refine_window - 1, -1):
                window = seg_for_refine.iloc[k-refine_window+1:k+1]
                w_lrr = window['lrr_filt'].median()
                if event == "LOSS" and w_lrr < -0.10:
                    end_pos = int(window['position'].iloc[-1])
                    break
                elif event == "GAIN" and w_lrr > 0.10:
                    end_pos = int(window['position'].iloc[-1])
                    break
            seg_all = d[(d["position"] >= start_pos) & (d["position"] <= end_pos)]
            het_seg = seg_all[seg_all["is_het_like"] & seg_all["maf"].notna()].copy()
            all_lrr_vals = seg_all['lrr_filt'].dropna()
            median_lrr = float(np.nanmedian(all_lrr_vals))
            mean_ai = float(np.nanmean(het_seg["ai"]))
            span_mb = (end_pos - start_pos) / 1e6
            if DEBUG: print(f"  [{chrom}] Post-override LRR trim: {start_pos/1e6:.1f}-{end_pos/1e6:.1f} Mb, median_lrr={median_lrr:.3f}, span={span_mb:.1f} Mb")

            # RE-VALIDATE after trim
            if span_mb < min_span_mb:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED: post-override trim reduced span to {span_mb:.1f} Mb < {min_span_mb}")
                continue
            
            if len(het_seg) < 8:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED: post-override trim left only {len(het_seg)} hets < 8")
                continue
            
            min_support_post = max(8, span_mb * 1.0)
            n_supporting_post = (het_seg["ai"] >= enter_thresh).sum()
            if n_supporting_post < min_support_post:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED: post-override trim: {n_supporting_post} supporting hets < {min_support_post:.0f}")
                continue

        # === POST-OVERRIDE cf_lrr CHECK ===
        # After BLOCKs 1-3, event may have changed to GAIN/LOSS but cf_lrr is still from
        # the original classification. Re-check that GAIN/LOSS has valid LRR confirmation.
        if event in ("GAIN", "LOSS") and not np.isfinite(cf_lrr):
            if DEBUG: print(f"  [{chrom}] ❌ KILLED: {event} with no LRR cell fraction after override (cf_lrr=NaN)")
            continue

        # === EXISTING Filter C ===
        # Compute BAF statistical strength for adaptive LRR threshold
        _bg_hets_fc = d[~((d['position'] >= start_pos) & (d['position'] <= end_pos))]
        _bg_hets_fc = _bg_hets_fc[_bg_hets_fc['is_het_like'] & _bg_hets_fc['maf'].notna()]
        if len(_bg_hets_fc) > 20 and len(het_seg) > 0:
            _bg_ai_med = _bg_hets_fc['ai'].median()
            _bg_ai_mad = np.median(np.abs(_bg_hets_fc['ai'] - _bg_ai_med))
            _bg_ai_sigma = 1.4826 * _bg_ai_mad
            _event_ai_med = het_seg['ai'].median()
            _se = _bg_ai_sigma / np.sqrt(len(het_seg))
            baf_z = (_event_ai_med - _bg_ai_med) / _se if _se > 0 else 0
        else:
            baf_z = 0

        # Strong BAF = z > 6, at least 30 hets, and LRR in correct direction
        strong_baf_for_filterC = (
            baf_z > 6.0 and
            len(het_seg) >= 30 and
            ((event == "GAIN" and median_lrr > 0) or
             (event == "LOSS" and median_lrr < 0))
        )
        lrr_contrast_required = 0.08 if strong_baf_for_filterC else 0.12

        if len(background_data) > 50:
            background_median_lrr = np.nanmedian(background_data)
            lrr_contrast = median_lrr - background_median_lrr
            
            if DEBUG: print(f"  [{chrom}] Filter C (LRR contrast): contrast={lrr_contrast:.3f}, bg_median={background_median_lrr:.3f}, event_median={median_lrr:.3f}, required={lrr_contrast_required:.2f}, baf_z={baf_z:.1f}")

            if event == "GAIN" and lrr_contrast < lrr_contrast_required and not lrr_override_fired:
                if median_lrr >= 0.10:
                    confidence = "MEDIUM"
                    if DEBUG: print(f"  [{chrom}] Filter C: GAIN contrast weak but absolute median_lrr={median_lrr:.3f} confirms, downgraded to MEDIUM")
                elif delta > 0.15 and cf_baf > 0.3:
                    confidence = "MEDIUM"
                else:
                    if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER C: GAIN contrast {lrr_contrast:.3f} < {lrr_contrast_required:.2f}")
                    continue
            if event == "LOSS" and lrr_contrast > -lrr_contrast_required and not lrr_override_fired:
                if median_lrr <= -0.10:
                    confidence = "MEDIUM"
                    if DEBUG: print(f"  [{chrom}] Filter C: LOSS contrast weak but absolute median_lrr={median_lrr:.3f} confirms, downgraded to MEDIUM")
                elif delta > 0.15 and cf_baf > 0.3:
                    confidence = "MEDIUM"
                else:
                    if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER C: LOSS contrast {lrr_contrast:.3f} > {-lrr_contrast_required:.2f}")
                    continue
        else:
            if DEBUG: print(f"  [{chrom}] Filter C: skipped (only {len(background_data)} background points)")


        # Filter D: Chaos Filter
        if len(final_supporting_snps) > 5:
            ai_sd = final_supporting_snps["ai"].std()
            ai_mean = final_supporting_snps["ai"].mean()
            if DEBUG: print(f"  [{chrom}] Filter D (Chaos): ai_sd={ai_sd:.3f}, ai_mean={ai_mean:.3f}")
            if ai_mean < 0.08:
                if ai_sd > 0.08:
                    if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER D: low mean ({ai_mean:.3f}) + high sd ({ai_sd:.3f})")
                    continue
            if event in ("GAIN", "LOSS") and abs(median_lrr) > 0.15:
                chaos_max = 0.25
            elif event in ("GAIN", "LOSS") and abs(median_lrr) > 0.10:
                chaos_max = 0.20
            else:
                chaos_max = 0.15
            
            if ai_sd > chaos_max:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER D: ai_sd {ai_sd:.3f} > {chaos_max:.3f}")
                continue

        # Filter E: Peak-to-Median
        if len(final_supporting_snps) > 5:
            ai_values = final_supporting_snps["ai"]
            median_val = ai_values.median()
            peak_val = np.percentile(ai_values, 95)
            
            if DEBUG: print(f"  [{chrom}] Filter E (Peak-to-Median): median={median_val:.3f}, p95={peak_val:.3f}")

            if median_val < 0.04 and peak_val > 0.20 and span_mb < 10.0:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER E: spiky signal (median={median_val:.3f}, p95={peak_val:.3f})")
                continue
            if median_val < 0.03 and span_mb < 10.0:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER E: very low median ({median_val:.3f})")
                continue

        # Filter F: AI Contrast
        background_hets = d[~((d['position'] >= start_pos) & (d['position'] <= end_pos))]
        background_hets = background_hets[background_hets['is_het_like'] & background_hets['maf'].notna()]
        if len(background_hets) > 20:
            background_median_ai = background_hets['ai'].median()
            background_mad = np.median(np.abs(background_hets['ai'] - background_median_ai))
            background_sigma = 1.4826 * background_mad
            
            event_median_ai = het_seg['ai'].median()
            ai_contrast = event_median_ai - background_median_ai
            ai_ratio = event_median_ai / background_median_ai if background_median_ai > 0 else 99
            
            n_event_hets = len(het_seg)
            se_mean = background_sigma / np.sqrt(n_event_hets) if n_event_hets > 0 else background_sigma
            z_score = ai_contrast / se_mean if se_mean > 0 else 0
            
            passes_absolute = (ai_contrast > 0.05) and (ai_ratio > 1.5)
            passes_statistical = (z_score > 4.0) and (ai_contrast > 0.01) and (n_event_hets >= 30)

            lrr_confirms = (
                (event == "GAIN" and lrr_contrast is not None and lrr_contrast > 0.10) or
                (event == "LOSS" and lrr_contrast is not None and lrr_contrast < -0.10)
            )
            passes_lrr_rescue = lrr_confirms and (ai_ratio > 2.0) and (z_score > 3.0)
            
            if DEBUG: print(f"  [{chrom}] Filter F (AI contrast): contrast={ai_contrast:.3f}, ratio={ai_ratio:.1f}, z={z_score:.1f}, passes_abs={passes_absolute}, passes_stat={passes_statistical}")

            if not (passes_absolute or passes_statistical or passes_lrr_rescue):
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER F: AI not elevated vs background")
                continue
        else:
            if DEBUG: print(f"  [{chrom}] Filter F: skipped (only {len(background_hets)} background hets)")

        # Filter F genome-wide extension for large events
        chr_size_mb = chromosome_sizes.get(chrom, 250e6) / 1e6
        event_fraction = span_mb / chr_size_mb
        if span_mb > 15.0 and (len(background_hets) < 100 or event_fraction > 0.25):
            genome_bg = df[
                (df["chromosome"] != chrom) & 
                df["is_het_like"] & df["maf"].notna()
            ]
            if len(genome_bg) > 100:
                genome_bg_median_ai = genome_bg['ai'].median()
                genome_bg_mad = np.median(np.abs(genome_bg['ai'] - genome_bg_median_ai))
                genome_bg_sigma = 1.4826 * genome_bg_mad
                
                ai_contrast_genome = event_median_ai - genome_bg_median_ai
                se_genome = genome_bg_sigma / np.sqrt(n_event_hets)
                z_genome = ai_contrast_genome / se_genome if se_genome > 0 else 0
                
                if DEBUG: print(f"  [{chrom}] Filter F (genome-wide): "
                            f"contrast={ai_contrast_genome:.3f}, z={z_genome:.1f}, "
                            f"genome_bg={genome_bg_median_ai:.3f}")
                
                if ai_contrast_genome < 0.03 or z_genome < 5.0:
                    if DEBUG: print(f"  [{chrom}] ❌ KILLED: not elevated vs genome-wide background")
                    continue

        # Filter G: BAF-LRR agreement
        if np.isfinite(cls["lrr_fit_score"]) and cls["lrr_fit_score"] > 4.0:
            if event == "CN-LOH" and mean_ai > 0.15 and span_mb > 5:
                if DEBUG: print(f"  [{chrom}] Filter G: rescued CN-LOH (strong BAF), downgraded to MEDIUM")
                confidence = "MEDIUM"
            elif event in ("GAIN", "LOSS") and delta > 0.15 and lrr_contrast is not None:
                if DEBUG: print(f"  [{chrom}] Filter G: rescued GAIN or LOSS, downgraded to MEDIUM")
                lrr_direction_correct = (event == "GAIN" and lrr_contrast > 0.05) or \
                                        (event == "LOSS" and lrr_contrast < -0.05)
                if lrr_direction_correct:
                    confidence = "MEDIUM"
                else:
                    continue
            else:
                if DEBUG: print(f"  [{chrom}] ❌ KILLED at FILTER G: lrr_fit_score {cls['lrr_fit_score']:.2f} > 4.0")
                continue

        if DEBUG: print(f"  [{chrom}] ✅ PASSED ALL FILTERS — event={event}, span={span_mb:.1f} Mb, cf_baf={cf_baf_str}")

        # Record Result
        rec = {
            "chromosome": chrom, "start_pos": start_pos, "end_pos": end_pos, "span_mb": float(span_mb),
            "event": event, "cell_fraction_BAF": float(cf_baf) if np.isfinite(cf_baf) else np.nan,
            "cell_fraction_AI": float(cf_ai), 
            "cell_fraction_LRR": float(cf_lrr) if np.isfinite(cf_lrr) else np.nan,
            "cell_fraction_combined": float((cf_baf + cf_lrr) / 2) if (np.isfinite(cf_baf) and np.isfinite(cf_lrr)) else float(cf_baf) if np.isfinite(cf_baf) else np.nan,
            "detector": "BAF_AI",
            "n_hets": int(len(het_seg)), "mean_AI": float(mean_ai),
            "mean_LRR": float(mean_lrr) if np.isfinite(mean_lrr) else np.nan,
            "median_LRR": float(median_lrr) if np.isfinite(median_lrr) else np.nan,
            "median_LRR_corrected": float(median_lrr_corrected) if np.isfinite(median_lrr_corrected) else np.nan,
            "background_median_LRR": float(background_median_lrr_for_class),
            "mu1": float(mu1) if np.isfinite(mu1) else np.nan, "mu2": float(mu2) if np.isfinite(mu2) else np.nan,
            "delta_mu2_mu1": float(delta) if np.isfinite(delta) else np.nan,
            "enter_thresh": float(enter_thresh), "exit_thresh": float(exit_thresh),
            "window_hets": int(window_hets), "step_hets": int(step_hets),
            "lrr_sigma_chr": float(lrr_sigma),
            "lrr_expected_from_p_baf": float(cls.get("lrr_expected_from_p_baf", np.nan)),
            "lrr_fit_score": float(cls.get("lrr_fit_score", np.nan)),
            "flag_small_delta": bool(flag_small_delta), "flag_low_hets": bool(flag_low_hets),
            "flag_lrr_baf_mismatch": bool(flag_lrr_baf_mismatch), "flag_lrr_baf_mismatch_cls": bool(flag_lrr_baf_mismatch_cls),
            "flag_neutral_lrr_large_ai": bool(flag_neutral_lrr_large_ai), "confidence": confidence,
        }
        if "RSID" in seg_all.columns:
            rec["start_RSID"] = str(seg_all.iloc[0].get("RSID", "")); rec["end_RSID"] = str(seg_all.iloc[-1].get("RSID", ""))
        calls.append(rec)
            
    calls_df = pd.DataFrame(calls)
    if not calls_df.empty and "start_pos" in calls_df.columns:
        calls_df = calls_df.sort_values(["start_pos"]).reset_index(drop=True)
    
    return calls_df, wdf

def detect_arm_level_events(df, chrom, het_lo, het_hi, centromere_dict, 
                            min_z_score=8.0, min_hets=50):
    """
    Detect subtle arm-level events by comparing per-arm AI distributions.
    Bypasses the sliding window state machine entirely.
    """
    d = df[df["chromosome"] == chrom].sort_values("position").reset_index(drop=True)
    h = d[d["is_het_like"] & d["maf"].notna()].copy()
    
    if len(h) < 30 or chrom not in centromere_dict:
        return pd.DataFrame()
    
    cen_start = centromere_dict[chrom]['start']
    cen_end = centromere_dict[chrom]['end']

    # Calculate genome-wide background from all OTHER chromosomes
    other_chrom_hets = df[
        (df["chromosome"] != chrom) & 
        df["is_het_like"] & df["maf"].notna()
    ]
    genome_bg_ai = other_chrom_hets['ai'].values
    genome_bg_median = np.median(genome_bg_ai)
    genome_bg_sigma = 1.4826 * np.median(np.abs(genome_bg_ai - genome_bg_median))
    
    p_hets = h[h['position'] < cen_start]
    q_hets = h[h['position'] > cen_end]
    
    calls = []
    
    for arm_label, arm_hets, arm_start, arm_end in [
        ('p', p_hets, int(p_hets['position'].min()) if len(p_hets) > 0 else 0, cen_start),
        ('q', q_hets, cen_end, int(q_hets['position'].max()) if len(q_hets) > 0 else d['position'].max())
    ]:
        if len(arm_hets) < min_hets:
            continue
        
        other_hets = q_hets if arm_label == 'p' else p_hets
        
        arm_ai = arm_hets['ai'].values
        arm_median = np.median(arm_ai)
        
        # Use other arm as background if available, otherwise genome-wide
        if len(other_hets) >= 20:
            other_median = np.median(other_hets['ai'].values)
            bg_sigma = 1.4826 * np.median(np.abs(other_hets['ai'].values - other_median))
        else:
            other_median = genome_bg_median
            bg_sigma = genome_bg_sigma
        
        if bg_sigma <= 0:
            continue
        
        ai_contrast = arm_median - other_median
        se_mean = bg_sigma / np.sqrt(len(arm_ai))
        z_score = ai_contrast / se_mean if se_mean > 0 else 0
        
        if z_score < min_z_score or ai_contrast < 0.025:
            continue
        
        # LRR comparison between arms — use het LRR when available, fall back to ALL SNP LRR
        arm_lrr = arm_hets['lrr_filt'].dropna()
        other_lrr = other_hets['lrr_filt'].dropna()
        arm_median_lrr = np.nanmedian(arm_lrr) if len(arm_lrr) > 0 else np.nan
        other_median_lrr = np.nanmedian(other_lrr) if len(other_lrr) > 0 else np.nan
        
        # NEW: If het-based LRR is unavailable (100% events deplete hets), use ALL SNPs
        if np.isnan(arm_median_lrr):
            arm_all_snps = d[(d['position'] >= int(arm_start)) & (d['position'] <= int(arm_end))]
            arm_all_lrr = arm_all_snps['lrr_filt'].dropna()
            if len(arm_all_lrr) > 10:
                arm_median_lrr = np.nanmedian(arm_all_lrr)
        
        if np.isnan(other_median_lrr):
            # Use genome-wide background LRR as fallback
            other_chrom_lrr = df[(df["chromosome"] != chrom)]['lrr_filt'].dropna()
            if len(other_chrom_lrr) > 50:
                other_median_lrr = np.nanmedian(other_chrom_lrr)
            else:
                other_median_lrr = 0.0
        
        lrr_contrast = arm_median_lrr - other_median_lrr if np.isfinite(arm_median_lrr) and np.isfinite(other_median_lrr) else np.nan
        
        # Classify event using LRR contrast (arm vs other arm/background)
        # This prevents misclassification due to per-sample LRR drift
        if np.isfinite(lrr_contrast):
            # Use LRR contrast for classification (inherently background-corrected)
            if abs(lrr_contrast) <= 0.10:
                event = "CN-LOH"
            elif lrr_contrast >= 0.20:
                event = "GAIN"
            elif lrr_contrast <= -0.20:
                event = "LOSS"
            else:
                event = "CN-LOH"  # Default for AI-only detection
        elif np.isfinite(arm_median_lrr):
            # Fallback to absolute LRR if contrast unavailable
            if abs(arm_median_lrr) <= 0.10:
                event = "CN-LOH"
            elif arm_median_lrr >= 0.20:
                event = "GAIN"
            elif arm_median_lrr <= -0.20:
                event = "LOSS"
            else:
                event = "CN-LOH"
        else:
            # Last resort: use ALL SNP LRR from the arm region (not just hets)
            # This catches 100% loss/gain events where het depletion makes het-LRR unavailable
            arm_all_snps = d[(d['position'] >= int(arm_start)) & (d['position'] <= int(arm_end))]
            arm_all_lrr_vals = arm_all_snps['lrr_filt'].dropna()
            if len(arm_all_lrr_vals) > 10:
                all_snp_median_lrr = np.nanmedian(arm_all_lrr_vals)
                if all_snp_median_lrr <= -0.20:
                    event = "LOSS"
                elif all_snp_median_lrr >= 0.20:
                    event = "GAIN"
                else:
                    event = "CN-LOH"
            else:
                event = "CN-LOH"

        print(f"  [{chrom}] arm_median_lrr={arm_median_lrr:.3f}, event={event}")
        
        # Cell fraction from AI
        mean_ai = np.mean(arm_ai)
        if event == "CN-LOH":
            cf_baf = 2 * mean_ai
        elif event == "GAIN":
            cf_baf = 4 * mean_ai / (1 - 2 * mean_ai) if mean_ai < 0.5 else 1.0
        elif event == "LOSS":
            cf_baf = 4 * mean_ai / (1 + 2 * mean_ai)
        else:
            cf_baf = 2 * mean_ai
        cf_baf = float(np.clip(cf_baf, 0, 1))
        
        cf_baf = float(np.clip(cf_baf, 0, 1))
        
        # --- BOUNDARY REFINEMENT ---
        # Walk inward from both ends to find where signal actually starts/stops
        arm_hets_sorted = arm_hets.sort_values('position').reset_index(drop=True)
        refined_start = int(arm_hets_sorted.iloc[0]['position'])
        refined_end = int(arm_hets_sorted.iloc[-1]['position'])
        
        chunk_size = max(10, len(arm_hets_sorted) // 20)
        
        # Trim left: find first chunk with AI above background
        for k in range(0, len(arm_hets_sorted) - chunk_size, chunk_size // 2):
            chunk = arm_hets_sorted.iloc[k:k + chunk_size]
            chunk_median_ai = chunk['ai'].median()
            if chunk_median_ai > other_median + ai_contrast * 0.3:
                refined_start = int(chunk.iloc[0]['position'])
                break
        
        # Trim right: find last chunk with AI above background
        for k in range(len(arm_hets_sorted) - 1, chunk_size - 1, -(chunk_size // 2)):
            chunk = arm_hets_sorted.iloc[k - chunk_size + 1:k + 1]
            chunk_median_ai = chunk['ai'].median()
            if chunk_median_ai > other_median + ai_contrast * 0.3:
                refined_end = int(chunk.iloc[-1]['position'])
                break
        
        arm_start = refined_start
        arm_end = refined_end
        # --- END BOUNDARY REFINEMENT ---
        
        span_mb = (arm_end - arm_start) / 1e6
        
        confidence = "HIGH" if z_score > 6.0 else "MEDIUM"
        
        calls.append({
            'chromosome': chrom,
            'start_pos': int(arm_start),
            'end_pos': int(arm_end),
            'span_mb': span_mb,
            'event': event,
            'confidence': confidence,
            'cell_fraction_BAF': cf_baf,
            'cell_fraction_LRR': np.nan,
            'cell_fraction_combined': cf_baf,
            'n_hets': len(arm_hets),
            'n_total_snps': len(d[(d['position'] >= arm_start) & (d['position'] <= arm_end)]),
            'mean_AI': float(mean_ai),
            'median_LRR': float(arm_median_lrr) if np.isfinite(arm_median_lrr) else np.nan,
            'mean_LRR': float(np.nanmean(arm_lrr)) if len(arm_lrr) > 0 else np.nan,
            'detector': 'arm_level',
            'z_score': float(z_score),
            'ai_contrast': float(ai_contrast),
            'note': f'arm-level {arm_label}-arm (z={z_score:.1f})'
        })
    
    return pd.DataFrame(calls) if calls else pd.DataFrame()

def _check_extension_signal(df_by_chrom, chrom, ext_start, ext_end, ai_thresh=0.04, min_hets=5):
    """
    Check whether het SNPs in a boundary extension zone show real AI signal.
    
    Returns True if the extension has evidence of a real event (keep broader boundary),
    False if it looks like noise (use narrower boundary).
    """
    if chrom not in df_by_chrom:
        return False
    
    d = df_by_chrom[chrom]
    ext_hets = d[(d['position'] >= ext_start) & (d['position'] <= ext_end)
                 & d['is_het_like'] & d['ai'].notna()]
    
    if len(ext_hets) < min_hets:
        # Too few hets to tell — conservatively use narrower boundary
        return False
    
    median_ai = ext_hets['ai'].median()
    return median_ai >= ai_thresh


# ==========================================
# HELPER FUNCTIONS - TO CHECK HOW WELL DETECTED MCA MATCHES EXPECTED MCA
# ==========================================

def deduplicate_calls(calls, df_by_chrom=None):
    """Merge adjacent same-event calls, then deduplicate overlapping calls."""
    if len(calls) == 0:
        return pd.DataFrame()
    
    all_deduped = []
    
    for chrom in calls['chromosome'].unique():
        chrom_calls = calls[calls['chromosome'] == chrom].sort_values(
            ['start_pos', 'span_mb'], ascending=[True, False]
        ).reset_index(drop=True)
        
        # Step 1: Merge adjacent same-event calls
        merged = [chrom_calls.iloc[0].to_dict()]

        for i in range(1, len(chrom_calls)):
            current = merged[-1]
            next_call = chrom_calls.iloc[i].to_dict()
            
            gap_mb = (next_call['start_pos'] - current['end_pos']) / 1e6
            same_event = current.get('event') == next_call.get('event')

            merge = False

            if same_event:
                if gap_mb < 2:
                    merge = True
                elif gap_mb < 5:
                    lrr_diff = abs(current.get('median_LRR', 0) - next_call.get('median_LRR', 0))
                    if lrr_diff < 0.3:
                        merge = True

            if merge:
                calls_overlap = current['end_pos'] >= next_call['start_pos']

                if calls_overlap:
                    narrow_start = max(current['start_pos'], next_call['start_pos'])
                    narrow_end = min(current['end_pos'], next_call['end_pos'])
                    wide_start = min(current['start_pos'], next_call['start_pos'])
                    wide_end = max(current['end_pos'], next_call['end_pos'])
                    
                    cur_span = (current['end_pos'] - current['start_pos']) / 1e6
                    nxt_span = (next_call['end_pos'] - next_call['start_pos']) / 1e6
                    larger = current if cur_span >= nxt_span else next_call
                    larger_span = max(cur_span, nxt_span)
                    smaller_span = max(min(cur_span, nxt_span), 0.1)
                    
                    larger_is_strong = (
                        'BAF_AI' in larger.get('detector', '') and
                        larger.get('confidence', '') in ('HIGH', 'MEDIUM') and
                        larger_span > 15.0
                    )
                    size_dominates = (larger_span / smaller_span) > 3.0 and larger_span > 10.0
                    both_are_large = min(cur_span, nxt_span) > 10.0 and larger_span > 15.0
                    
                    if larger_is_strong or size_dominates or both_are_large:
                        current['start_pos'] = wide_start
                        current['end_pos'] = wide_end
                    else:
                        if df_by_chrom is not None and wide_start < narrow_start:
                            if _check_extension_signal(df_by_chrom, chrom, wide_start, narrow_start):
                                current['start_pos'] = wide_start
                            else:
                                current['start_pos'] = narrow_start
                        else:
                            current['start_pos'] = wide_start
                        
                        if df_by_chrom is not None and wide_end > narrow_end:
                            if _check_extension_signal(df_by_chrom, chrom, narrow_end, wide_end):
                                current['end_pos'] = wide_end
                            else:
                                current['end_pos'] = narrow_end
                        else:
                            current['end_pos'] = wide_end

                else:
                    current['start_pos'] = min(current['start_pos'], next_call['start_pos'])
                    current['end_pos'] = max(current['end_pos'], next_call['end_pos'])
                
                current['span_mb'] = (current['end_pos'] - current['start_pos']) / 1e6
                existing = set(current.get('detector', '').split('+'))
                new = set(next_call.get('detector', '').split('+'))
                current['detector'] = '+'.join(sorted(existing | new))
            else:
                merged.append(next_call)

        # # DEBUG: Show Step 1 output
        # for m in merged:
        #     print(f"  [DEDUP Step1] {chrom} {m.get('detector','?'):30s} {m.get('event','?'):8s} {m['start_pos']/1e6:8.1f}-{m['end_pos']/1e6:8.1f} Mb  span={m.get('span_mb',0):.1f}")

        # Step 2: Deduplicate overlapping calls
        merged = sorted(merged, key=lambda x: x['span_mb'], reverse=True)
        keep = []
        for call in merged:
            overlaps = False
            for i, kept in enumerate(keep):
                overlap_start = max(call['start_pos'], kept['start_pos'])
                overlap_end = min(call['end_pos'], kept['end_pos'])
                if overlap_start < overlap_end:
                    overlap_mb = (overlap_end - overlap_start) / 1e6
                    if overlap_mb > min(call['span_mb'], kept['span_mb']) * 0.5:
                        overlaps = True
                        
                        # --- EVENT TYPE ARBITRATION ---
                        if kept['event'] != call['event']:
                            kept_lrr = kept.get('median_LRR', np.nan)
                            call_lrr = call.get('median_LRR', np.nan)
                            
                            best_lrr = np.nan
                            best_lrr_source = None
                            for c in [kept, call]:
                                det = c.get('detector', '')
                                c_lrr = c.get('median_LRR', np.nan)
                                if np.isfinite(c_lrr):
                                    if 'homozygous' in det or 'lrr_baf' in det:
                                        best_lrr = c_lrr
                                        best_lrr_source = c
                                        break
                                    elif np.isnan(best_lrr):
                                        best_lrr = c_lrr
                                        best_lrr_source = c
                            
                            if np.isfinite(best_lrr):
                                if best_lrr <= -0.20:
                                    kept['event'] = 'LOSS'
                                elif best_lrr >= 0.20:
                                    kept['event'] = 'GAIN'
                                elif abs(best_lrr) <= 0.06:
                                    kept['event'] = 'CN-LOH'
                                elif abs(best_lrr) > 0.06:
                                    if 'lrr_baf' in call.get('detector', '') and 'lrr_baf' not in kept.get('detector', ''):
                                        kept['event'] = call['event']
                                    elif 'lrr_baf' in kept.get('detector', '') and 'lrr_baf' not in call.get('detector', ''):
                                        pass

                            # Sync call event after arbitration
                            call['event'] = kept['event']
                        
                        # --- BOUNDARY HANDLING ---
                        if kept['event'] == call['event']:
                            narrow_start = max(kept['start_pos'], call['start_pos'])
                            narrow_end = min(kept['end_pos'], call['end_pos'])
                            wide_start = min(kept['start_pos'], call['start_pos'])
                            wide_end = max(kept['end_pos'], call['end_pos'])
                            
                            larger = kept if kept['span_mb'] >= call['span_mb'] else call
                            smaller = call if kept['span_mb'] >= call['span_mb'] else kept
                            
                            larger_is_strong_baf_ai = (
                                'BAF_AI' in larger.get('detector', '') and
                                larger.get('confidence', '') in ('HIGH', 'MEDIUM') and
                                larger.get('span_mb', 0) > 15.0 and
                                not np.isnan(larger.get('delta_mu2_mu1', np.nan)) and
                                larger.get('delta_mu2_mu1', 0) > 0.15
                            )
                            size_ratio = larger.get('span_mb', 0) / max(smaller.get('span_mb', 1), 0.1)
                            larger_dominates = size_ratio > 3.0 and larger.get('span_mb', 0) > 10.0
                            both_are_large = min(larger.get('span_mb', 0), smaller.get('span_mb', 0)) > 10.0 and larger.get('span_mb', 0) > 15.0
                            
                            if larger_is_strong_baf_ai or larger_dominates or both_are_large:
                                kept['start_pos'] = wide_start
                                kept['end_pos'] = wide_end
                            else:
                                if df_by_chrom is not None and wide_start < narrow_start:
                                    if _check_extension_signal(df_by_chrom, chrom, wide_start, narrow_start):
                                        kept['start_pos'] = wide_start
                                    else:
                                        kept['start_pos'] = narrow_start
                                else:
                                    kept['start_pos'] = wide_start
                                
                                if df_by_chrom is not None and wide_end > narrow_end:
                                    if _check_extension_signal(df_by_chrom, chrom, narrow_end, wide_end):
                                        kept['end_pos'] = wide_end
                                    else:
                                        kept['end_pos'] = narrow_end
                                else:
                                    kept['end_pos'] = wide_end
                        else:
                            kept['start_pos'] = max(kept['start_pos'], call['start_pos'])
                            kept['end_pos'] = min(kept['end_pos'], call['end_pos'])

                        kept['span_mb'] = (kept['end_pos'] - kept['start_pos']) / 1e6
                        
                        existing_det = set(kept.get('detector', '').split('+'))
                        new_det = set(call.get('detector', '').split('+'))
                        kept['detector'] = '+'.join(sorted(existing_det | new_det))                        
                        break
            if not overlaps:
                keep.append(call)

        # DEBUG: Show Step 2 output
        # for k in keep:
        #     print(f"  [DEDUP Step2] {chrom} {k.get('detector','?'):30s} {k.get('event','?'):8s} {k['start_pos']/1e6:8.1f}-{k['end_pos']/1e6:8.1f} Mb  span={k.get('span_mb',0):.1f}")

        all_deduped.extend(keep)
    
    result = pd.DataFrame(all_deduped)
    if not result.empty:
        result = result.sort_values('start_pos').reset_index(drop=True)
    return result

def match_call_to_truth(call, truth, df, centromere_dict = centromere_dict, OVERLAP_THRESHOLD=0.5, debug = True, match_type = 'standard'):
    """
    Match detected call to truth with size-dependent length tolerance.
    """
    # Normalize chromosome names (remove 'chr' prefix)
    call_chrom = str(call['chromosome']).replace('chr', '')
    truth_chrom = str(truth['Chromosome']).replace('chr', '')
    
    if call_chrom != truth_chrom:
        return False, 0.0
    
    # Normalize event types (remove hyphens/underscores, uppercase)
    call_event = str(call['event']).replace('-', '').replace('_', '').upper()
    truth_event = str(truth['Type']).replace('-', '').replace('_', '').upper()
    
    if call_event != truth_event:
        return False, 0.0
    
    # Calculate overlap
    overlap_start = max(call['start_pos'], truth['Start_bp'])
    overlap_end = min(call['end_pos'], truth['End_bp'])
    
    if overlap_start >= overlap_end:
        return False, 0.0
    
    overlap_bp = overlap_end - overlap_start
    truth_length = truth['End_bp'] - truth['Start_bp']
    call_length = call['end_pos'] - call['start_pos']
    
    # Reciprocal overlap
    overlap_frac_truth = overlap_bp / truth_length
    overlap_frac_call = overlap_bp / call_length
    min_overlap = min(overlap_frac_truth, overlap_frac_call)
    
    if min_overlap < OVERLAP_THRESHOLD:
        return False, 0.0
    
    # === PANEL-AWARE WHOLE-ARM VALIDATION ===
    # Check for whole-arm events (multiple possible geometry names)
    geometry_str = str(truth.get('Geometry', '')).lower()
    is_whole_arm = any(keyword in geometry_str for keyword in ['whole', 'wholearm', 'whole_arm'])
    
    # Also check length category for "Whole_Arm"
    length_cat = str(truth.get('Length_Category', '')).lower()
    is_whole_arm = is_whole_arm or any(keyword in length_cat for keyword in ['whole', 'wholearm', 'whole_arm'])
    
    if debug:
        print(f"  Geometry: '{truth.get('Geometry', 'N/A')}'")
        print(f"  Length_Category: '{truth.get('Length_Category', 'N/A')}'")
        print(f"  Is whole-arm: {is_whole_arm}")
    
    if is_whole_arm and df is not None and centromere_dict is not None:
        chrom = call['chromosome']
        chrom_snps = df[df['chromosome'] == chrom].copy()
        
        if len(chrom_snps) > 0 and chrom in centromere_dict:
            first_snp = chrom_snps['position'].min()
            last_snp = chrom_snps['position'].max()
            centromere_start = centromere_dict[chrom]['start']
            centromere_end = centromere_dict[chrom]['end']
            
            # Determine if truth is p-arm or q-arm
            truth_center = (truth['Start_bp'] + truth['End_bp']) / 2
            
            if truth_center < centromere_start:  # p-arm
                # Find first SNP cluster on p-arm
                p_arm_snps = chrom_snps[chrom_snps['position'] < centromere_start]
                if len(p_arm_snps) > 0:
                    expected_start = p_arm_snps['position'].min()
                    expected_end = p_arm_snps['position'].max()
                else:
                    expected_start = first_snp
                    expected_end = centromere_start
            else:  # q-arm
                # Find first SNP cluster on q-arm (skip centromere gap)
                q_arm_snps = chrom_snps[chrom_snps['position'] > centromere_end]
                if len(q_arm_snps) > 0:
                    expected_start = q_arm_snps['position'].min()  # ← Real q-arm coverage start
                    expected_end = q_arm_snps['position'].max()
                else:
                    expected_start = centromere_end
                    expected_end = last_snp
            
            expected_length = expected_end - expected_start
            
            if debug:
                print(f"  Expected coverage: {expected_start/1e6:.1f} - {expected_end/1e6:.1f} Mb ({expected_length/1e6:.1f} Mb)")
                print(f"  Detected: {call['start_pos']/1e6:.1f} - {call['end_pos']/1e6:.1f} Mb ({call_length/1e6:.1f} Mb)")
            
            # Check if call covers actual SNP coverage (not theoretical boundaries)
            call_covers_start = call['start_pos'] <= expected_start + 5e6
            call_covers_end = call['end_pos'] >= expected_end - 5e6
            
            if debug:
                print(f"  Covers start: {call_covers_start}")
                print(f"  Covers end: {call_covers_end}")
            
            if call_covers_start and call_covers_end:
                truth_mb = expected_length / 1e6
                call_mb = call_length / 1e6
                size_diff_mb = abs(call_mb - truth_mb)
                allowed_diff_mb = max(truth_mb * 0.30, 5.0)
                
                if debug:
                    print(f"  Panel-aware: expected={truth_mb:.1f} Mb, detected={call_mb:.1f} Mb, diff={size_diff_mb:.1f} Mb, allowed={allowed_diff_mb:.1f} Mb")
                
                if size_diff_mb <= allowed_diff_mb:
                    if debug:
                        print(f"  ✅ MATCHED via panel-aware!")
                    return True, min_overlap
    
    # =======================================
    
    # Size-dependent length tolerance
    truth_mb = truth_length / 1e6
    call_mb = call_length / 1e6
    size_diff_mb = abs(call_mb - truth_mb) 
    
    if truth_mb <= 10:
        max_tolerance = 0.50
    elif truth_mb < 30:
        max_tolerance = 0.30
    else:
        max_tolerance = 0.20
    
    absolute_tolerance_mb = 5.0
    allowed_diff_mb = max(truth_mb * max_tolerance, absolute_tolerance_mb)
    if debug:
        print(f"  Standard check: truth={truth_mb:.1f} Mb, call={call_mb:.1f} Mb, diff={size_diff_mb:.1f} Mb, allowed={allowed_diff_mb:.1f} Mb")
    
    if size_diff_mb > allowed_diff_mb:
        if debug:
            print(f"  ❌ Failed standard size check")
        return False, 0.0
    
    if debug:
        print(f"  ✅ MATCHED via standard logic!")

    # Return a match result dict:
    truth_mb = truth_length / 1e6
    call_mb = call_length / 1e6
    size_ratio = call_mb / truth_mb  # <1 = undersized, >1 = oversized
    
    if size_ratio < 0.7:
        size_accuracy = 'undersized'
    elif size_ratio > 1.3:
        size_accuracy = 'oversized'
    else:
        size_accuracy = 'accurate (+/- 30%)'
    
    result = {
        'matched': True,
        'overlap': min_overlap,
        'truth_mb': truth_mb,
        'call_mb': call_mb,
        'size_ratio': size_ratio,
        'size_accuracy': size_accuracy,
        'match_type': match_type  # 'standard', 'partial', 'panel_aware'
    }
    return result

def partial_match_call_to_truth(call, truth, OVERLAP_MIN=0.20):
    """
    Check if a call partially matches truth — same chromosome with some overlap,
    but may differ in event type or size. Returns a dict with match details.
    
    This distinguishes "found the event but misclassified/mis-sized it" from
    "called something on a completely wrong chromosome" (genuine FP).
    """
    call_chrom = str(call['chromosome']).replace('chr', '')
    truth_chrom = str(truth['Chromosome']).replace('chr', '')
    
    if call_chrom != truth_chrom:
        return {'partial_match': False, 'reason': 'wrong_chromosome'}
    
    # Calculate overlap
    overlap_start = max(call['start_pos'], truth['Start_bp'])
    overlap_end = min(call['end_pos'], truth['End_bp'])
    
    if overlap_start >= overlap_end:
        return {'partial_match': False, 'reason': 'no_overlap'}
    
    overlap_bp = overlap_end - overlap_start
    truth_length = truth['End_bp'] - truth['Start_bp']
    call_length = call['end_pos'] - call['start_pos']
    
    overlap_frac_truth = overlap_bp / truth_length if truth_length > 0 else 0
    overlap_frac_call = overlap_bp / call_length if call_length > 0 else 0
    min_overlap = min(overlap_frac_truth, overlap_frac_call)
    
    if min_overlap < OVERLAP_MIN:
        return {'partial_match': False, 'reason': 'insufficient_overlap', 'overlap': min_overlap}
    
    # Check what went wrong
    call_event = str(call['event']).replace('-', '').replace('_', '').upper()
    truth_event = str(truth['Type']).replace('-', '').replace('_', '').upper()
    
    type_match = (call_event == truth_event)
    
    truth_mb = truth_length / 1e6
    call_mb = call_length / 1e6
    size_ratio = call_mb / truth_mb if truth_mb > 0 else 0
    
    # Determine mismatch reason
    if not type_match:
        reason = 'TYPE_MISMATCH'
    elif size_ratio < 0.5:
        reason = 'UNDERSIZED'
    elif size_ratio > 2.0:
        reason = 'OVERSIZED'
    else:
        reason = 'SIZE_MISMATCH'
    
    return {
        'partial_match': True,
        'reason': reason,
        'overlap': min_overlap,
        'type_match': type_match,
        'call_type': call_event,
        'truth_type': truth_event,
        'size_ratio': size_ratio,
        'call_mb': call_mb,
        'truth_mb': truth_mb
    }

def parse_sample_filename(filename):
    """
    Extract sample info from filename
    Format: SAMPLE_TYPE_X%_chrN_LENGTHMb_GEOMETRY_SNPs.txt
    """
    parts = filename.replace('_SNPs.txt', '').replace('_PON_normalised_read_depths_and_LRR.txt', '')
    parts = parts.split('_')
    
    # Find the simulation ID parts
    # Format like: CNTRL_182_s4_GAIN_0.1pct_chr1_2Mb_Interstitial
    
    sample = '_'.join(parts[:3])  # e.g., CNTRL_182_s4
    sim_id = '_'.join(parts[3:])   # e.g., GAIN_0.1pct_chr1_2Mb_Interstitial
    
    return sample, sim_id

# ==========================================
# MAIN PROCESSING
# ==========================================

def process_simulated_samples():
    """Run mCA caller on all simulated samples"""

    known_event_chroms= {'CNTRL_160_s8': ['chr7'],
                    'CNTRL_169_s7': ['chr18'],
                    'CNTRL_174_s3': ['chr4', 'chr7', 'chr13', 'chrX'],
                    'CNTRL_177_s4': ['chrX'],
                    'CNTRL_181_s7': ['chr16'],
                    'CNTRL_182_s4': ['chr2', 'chr3'],
                    'CNTRL_183_s6': ['chr1'],
                    'CNTRL_186_s4': ['chr17'],
                    'CNTRL_188_s7': ['chr7', 'chr17', 'chrX'],
                    'CNTRL_193_s2': ['chr9', 'chrX'],
                    'CNTRL_199_s7': ['chr5'],
                    'CNTRL_163_s6': ['chr5', 'chr16', 'chrX'],
                    'CNTRL_164_s6': ['chrX'],
                    'CNTRL_167_s2': ['chr3', 'chr8'],
                    'CNTRL_171_s2': ['chr2'],
                    'CNTRL_172_s7': ['chr7'],
                    'CNTRL_175_s3': ['chr21'],
                    'CNTRL_180_s3': ['chr14'],
                    'CNTRL_187_s2': ['chr10', 'chr17', 'chrX'],
                    'CNTRL_190_s4': ['chrX'],
                    'CNTRL_191_s7': ['chr4', 'chr14'],
                    'CNTRL_194_s8': ['chr6'],
                    'CNTRL_001_s10': ['chr2', 'chr13', 'chrX'],
                    'CNTRL_002_s8': ['chr1', 'chr8'],
                    'CNTRL_003_s9': ['chr1', 'chr10'],
                    'CNTRL_004_s10': ['chr19', 'chrX'],
                    'CNTRL_005_s9': ['chr6', 'chr9', 'chrX'],
                    'CNTRL_162_s5': ['chr19', 'chrX'],
                    'CNTRL_185_s4': ['chr15'],
                    'CNTRL_189_s4': ['chr6', 'chr15'],
                    'CNTRL_195_s3': ['chr3'],
                    'CNTRL_196_s8': ['chr22'],
                    'CNTRL_198_s2': ['chr3', 'chr6']}

    centromere_dict = load_centromeres('chromosome_ideogram_hg19.txt')
    print("Loaded centromeres for:", list(centromere_dict.keys()))
    
    print("="*80)
    print("mCA CALLER ON SIMULATED SAMPLES")
    print("="*80)
    
    # Load truth manifest
    print(f"\n📋 Loading truth manifest: {MANIFEST_FILE}")
    if not os.path.exists(MANIFEST_FILE):
        print(f"❌ Manifest not found: {MANIFEST_FILE}")
        return
    
    manifest_df = pd.read_csv(MANIFEST_FILE)
    print(f"✅ Loaded {len(manifest_df)} simulated events")

    # Get all simulated files
    print(f"\n🔍 Scanning {SIMULATED_DIR} for simulated samples...")

    # Filter to specific cell fractions (comment out to run all)
    # target_fractions = [0.10, 0.15, 0.20, 0.25, 0.50, 0.75, 1.00]
    # manifest_df = manifest_df[manifest_df['Fraction_Percent'].isin(target_fractions)]
    # print(f"🔬 Filtered to {len(manifest_df)} events at CF = {target_fractions}")
    
    snp_files = []
    for lrr_filename in manifest_df['File_Name'].unique():
        snp_filename = lrr_filename.replace('_PON_normalised_read_depths_and_LRR.txt', '_SNPs.txt')
        snp_path = os.path.join(SIMULATED_DIR, snp_filename)
        if os.path.exists(snp_path):
            snp_files.append(snp_filename)

    print(f"✅ Found {len(snp_files)} simulated samples matching filtered manifest")
    
    # Output directories
    PLOTS_FP_DIR = os.path.join(OUTPUT_DIR, "False_positive_mCA_plots")
    os.makedirs(PLOTS_FP_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)
    os.makedirs(PLOTS_UNDETECTED_DIR, exist_ok=True)

    # Process each sample
    all_detected_calls = []
    truth_comparison = []
    processed_count = 0

    # === PROGRESS TRACKING ===
    start_time = time.time()
    progress_interval = 20
    # =========================
    
    for i, snp_file in enumerate(snp_files, 1):
        # Parse filename
        sample_base = snp_file.replace('_SNPs.txt', '')
        lrr_file = sample_base + '_PON_normalised_read_depths_and_LRR.txt'
        
        snp_path = os.path.join(SIMULATED_DIR, snp_file)
        lrr_path = os.path.join(SIMULATED_DIR, lrr_file)
        
        if not os.path.exists(lrr_path):
            print(f"⚠️  Missing LRR file for {sample_base}, skipping...")
            continue
        
        if i % progress_interval == 0:
            print(f"\n[{i}/{len(snp_files)}] {sample_base}")    
        
        # Find corresponding truth entry
        truth_entry = manifest_df[manifest_df['File_Name'] == lrr_file]
        
        if truth_entry.empty:
            print(f"   ⚠️  No truth entry found in manifest")
            continue
        
        truth_row = truth_entry.iloc[0]
        
        try:
            # Load SNP data
            df_snp = pd.read_csv(snp_path, sep='\t')
            df_snp = df_snp.rename(columns={'REF': 'ref', 'ALT': 'alt', 'VAF': 'BAF'})

            # Load the LRR data
            df_lrr = pd.read_csv(lrr_path, sep='\t', skiprows=5, header=None, names=['chromosome', 'start', 'stop', 'band', 'mean_region_depth', 'normalised_region_depth', 'log2ratio', 'coefficient_of_variation_PON', 'p-value', 'mean_normalised_depth_PON', 'sample_depth_normalised_by_PON'])
            df_lrr = df_lrr[['chromosome', 'start', 'stop', 'log2ratio', 'p-value', 'coefficient_of_variation_PON']].copy()
            df_lrr['start'] = pd.to_numeric(df_lrr['start'], errors='coerce')
            df_lrr['stop'] = pd.to_numeric(df_lrr['stop'], errors='coerce')
            df_lrr = df_lrr.sort_values(['chromosome', 'start'])

            df_lrr['chromosome'] = df_lrr['chromosome'].astype(str)
            df_snp['chromosome'] = df_snp['chromosome'].astype(str)
            df_snp['position'] = pd.to_numeric(df_snp['position'], errors='coerce')

            # FAST MERGE
            merged_dfs = []
            
            for chrom in df_snp['chromosome'].unique():
                snps_c = df_snp[df_snp['chromosome'] == chrom].sort_values('position')
                lrr_c = df_lrr[df_lrr['chromosome'] == chrom]
                
                if lrr_c.empty or snps_c.empty:
                    continue

                snps_c['position'] = snps_c['position'].astype('int64')
                lrr_c = lrr_c.dropna(subset=['start', 'stop'])
                lrr_c['start'] = lrr_c['start'].astype('int64')
                lrr_c['stop']  = lrr_c['stop'].astype('int64')
                    
                m = pd.merge_asof(snps_c, lrr_c, left_on='position', right_on='start', 
                                direction='backward', suffixes=('', '_lrr'))
                m = m[m['position'] <= m['stop']]
                m = m.rename(columns={'coefficient_of_variation_PON': 'probe_CV', 'p-value': 'LRR_p_value', 'start': 'probe_start', 'stop': 'probe_stop', 'log2ratio': 'LRR'})
                m['BAF'] = pd.to_numeric(m['BAF'], errors='coerce')
                m['BAF_deviation'] = abs(0.5 - m['BAF'])
                
                merged_dfs.append(m)

            if not merged_dfs:
                return pd.DataFrame()
                
            df = pd.concat(merged_dfs, ignore_index=True)

            # FILTER NOISY PROBES
            df["probe_CV"] = pd.to_numeric(df["probe_CV"], errors="coerce")
            CV_THRESH = 0.3
            df = df[df["probe_CV"] <= CV_THRESH].reset_index(drop=True)

            # CLEAN DATAFRAME
            df = clean_dataframe(df, het_lo, het_hi)
            df_for_plotting = df.copy()
            
            # Group by chromosome for processing ALL chromosomes
            df_by_chrom = dict(tuple(df.groupby('chromosome')))
            
            # Run caller on ALL chromosomes
            all_calls_all_chroms = []
            
            _genome_baseline_lrr = None  # will be computed on first chromosome

            for chromosome in all_chromosomes:
                if chromosome not in df_by_chrom:
                    continue
                
                calls_100_hom = detect_large_homozygous(df_by_chrom[chromosome], chromosome, centromere_dict=centromere_dict)
                calls_100_hom_small = detect_small_homozygous(df_by_chrom[chromosome], chromosome, centromere_dict=centromere_dict)
                calls_100_hom_small_gains = detect_small_gains(df_by_chrom[chromosome], chromosome)

                # Run LRR-first detector with BAF confirmation
                # Compute genome-wide LRR baseline once per sample
                if _genome_baseline_lrr is None:
                    _genome_lrr_medians = []
                    for _gc, _gdf in df_by_chrom.items():
                        _glrr = _gdf['lrr_filt'].dropna()
                        if len(_glrr) > 50:
                            _genome_lrr_medians.append(np.nanmedian(_glrr))
                    _genome_baseline_lrr = np.nanmedian(_genome_lrr_medians) if len(_genome_lrr_medians) >= 5 else 0.0
                    # print(f"  Genome-wide LRR baseline: {_genome_baseline_lrr:.4f} (from {len(_genome_lrr_medians)} chromosomes)")
                
                debug_chrom=None

                calls_lrr_baf = detect_LRR_BAF(df_by_chrom[chromosome], chromosome, chromosome_sizes, centromere_dict=centromere_dict, debug_chrom=debug_chrom, genome_baseline_lrr=_genome_baseline_lrr)

                calls_ai, w = detect_BAF_AI(df_by_chrom[chromosome], chromosome, window_hets, step_hets, enter_thresh, exit_thresh, 
                                                    gain_lrr, loss_lrr, neutral_lrr, min_span_mb, ai_point_thresh, gap_mb, het_lo, het_hi, 
                                                    min_consecutive_hets, min_spanning_hets, min_supporting_fraction, support_ai_thresh,
                                                    centromere_dict=centromere_dict, chromosome_sizes=chromosome_sizes)
                
                # --- MULTI-SCALE BAF_AI PASSES ---
                calls_ai_multiscale = []
                if multiscale_configs:
                    for ms_window, ms_step, ms_min_span, ms_min_consec, ms_min_spanning in multiscale_configs:
                        chrom_data = df_by_chrom[chromosome]
                        n_chrom_hets = chrom_data['is_het_like'].sum() if 'is_het_like' in chrom_data.columns else 0
                        if n_chrom_hets < ms_window * 2:
                            continue
                        calls_ms, _ = detect_BAF_AI(
                            df_by_chrom[chromosome], chromosome,
                            window_hets=ms_window, step_hets=ms_step,
                            enter_thresh=enter_thresh, exit_thresh=exit_thresh,
                            gain_lrr=gain_lrr, loss_lrr=loss_lrr, neutral_lrr=neutral_lrr,
                            min_span_mb=ms_min_span, ai_point_thresh=ai_point_thresh,
                            gap_mb=gap_mb, het_lo=het_lo, het_hi=het_hi,
                            min_consecutive_hets=ms_min_consec,
                            min_spanning_hets=ms_min_spanning,
                            min_supporting_fraction=min_supporting_fraction,
                            support_ai_thresh=support_ai_thresh,
                            centromere_dict=centromere_dict,
                            chromosome_sizes=chromosome_sizes)
                        if not calls_ms.empty:
                            calls_ms['detector'] = calls_ms['detector'].astype(str) + f'_w{ms_window}'
                            calls_ai_multiscale.append(calls_ms)
                            print(f"  BAF_AI w{ms_window} calls: {calls_ms[['start_pos','end_pos','event','span_mb']].to_string()}")

                calls_arm_level = detect_arm_level_events(df_by_chrom[chromosome], chromosome, het_lo, het_hi, centromere_dict=centromere_dict, min_z_score=4.0, min_hets=50)

                # Combine calls for this chromosome
                chrom_calls = []
                if not calls_100_hom.empty:
                    calls_100_hom_filtered = calls_100_hom[calls_100_hom['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                    if not calls_100_hom_filtered.empty:
                        chrom_calls.append(calls_100_hom_filtered)

                if not calls_100_hom_small.empty:
                    calls_100_hom_small_filtered = calls_100_hom_small[calls_100_hom_small['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                    if not calls_100_hom_small_filtered.empty:
                        chrom_calls.append(calls_100_hom_small_filtered)

                if not calls_100_hom_small_gains.empty:
                    calls_100_hom_small_gains_filtered = calls_100_hom_small_gains[calls_100_hom_small_gains['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                    if not calls_100_hom_small_gains_filtered.empty:
                        chrom_calls.append(calls_100_hom_small_gains_filtered)

                if not calls_lrr_baf.empty:
                    calls_lrr_baf_filtered = calls_lrr_baf[calls_lrr_baf['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                    if not calls_lrr_baf_filtered.empty:
                        chrom_calls.append(calls_lrr_baf_filtered)

                for ms_calls in calls_ai_multiscale:
                    if not ms_calls.empty:
                        ms_filtered = ms_calls[ms_calls['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                        if not ms_filtered.empty:
                            chrom_calls.append(ms_filtered)

                if not calls_ai.empty:
                    calls_ai_filtered = calls_ai[calls_ai['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                    if not calls_ai_filtered.empty:
                        chrom_calls.append(calls_ai_filtered)

                if not calls_arm_level.empty:
                    calls_arm_level_filtered = calls_arm_level[calls_arm_level['confidence'].isin(['MEDIUM', 'HIGH'])].copy()
                    if not calls_arm_level_filtered.empty and len(chrom_calls) > 0:
                        existing = pd.concat(chrom_calls, ignore_index=True)
                        
                        keep = []
                        for _, arm_call in calls_arm_level_filtered.iterrows():
                            already_covered = False
                            for _, ex_call in existing.iterrows():
                                overlap = calculate_overlap(
                                    arm_call['start_pos'], arm_call['end_pos'],
                                    ex_call['start_pos'], ex_call['end_pos']
                                )
                                if overlap > 0.3:
                                    already_covered = True
                                    break
                            if not already_covered:
                                keep.append(arm_call)
                        
                        if keep:
                            chrom_calls.append(pd.DataFrame(keep))
                    
                    elif not calls_arm_level_filtered.empty:
                        chrom_calls.append(calls_arm_level_filtered)
                
                if len(chrom_calls) > 0:
                    chrom_detected = pd.concat(chrom_calls, ignore_index=True)
                    all_calls_all_chroms.append(chrom_detected)
            
            # Combine all calls from all chromosomes
            all_calls = []
            if len(all_calls_all_chroms) > 0:
                all_calls.append(pd.concat(all_calls_all_chroms, ignore_index=True))
            
            # ==========================================
            # CLASSIFY CALLS AND SAVE RESULTS
            # ==========================================
            
            # Get known germline events for this control
            sample_base_control = '_'.join(sample_base.split('_')[:3])  # e.g. 'CNTRL_190_s4'
            known_germline_chroms = set(known_event_chroms.get(sample_base_control, []))
            
            if len(all_calls) > 0:
                detected_calls = pd.concat(all_calls, ignore_index=True)
                # NEW: Snap telomeric boundaries before dedup
                detected_calls = snap_telomeric_boundaries(detected_calls, df_by_chrom, chromosome_sizes, centromere_dict=centromere_dict, max_gap_mb=3.0)
                final_calls = deduplicate_calls(detected_calls)
                final_calls['sample'] = sample_base
                final_calls_df = final_calls

                print(f"   📊 Detected {len(final_calls_df)} call(s)")

                # --- Try to match to truth ---
                best_match = None
                best_overlap = 0.0
                best_partial_match = None
                best_partial_info = None
                
                for idx, call in final_calls_df.iterrows():
                    result = match_call_to_truth(call, truth_row, df_for_plotting, centromere_dict=centromere_dict, OVERLAP_THRESHOLD=0.5)
                    if isinstance(result, dict):
                        is_match = result['matched']
                        overlap = result['overlap']
                    else:
                        is_match, overlap = result

                    if is_match and overlap > best_overlap:
                        best_match = call
                        best_overlap = overlap
                    
                    # If no full match, check for partial match (same chrom, some overlap, wrong type/size)
                    if not is_match:
                        pm = partial_match_call_to_truth(call, truth_row, OVERLAP_MIN=0.20)
                        if pm['partial_match']:
                            if best_partial_info is None or pm['overlap'] > best_partial_info['overlap']:
                                best_partial_match = call
                                best_partial_info = pm
                
                if best_match is None and best_partial_match is not None:
                    print(f"\n  ⚠️ PARTIAL MATCH on {truth_row['Chromosome']}: {best_partial_info['reason']}")
                    print(f"    Truth: {truth_row['Type']} {truth_row['Start_bp']:,} - {truth_row['End_bp']:,} ({best_partial_info['truth_mb']:.1f} Mb)")
                    print(f"    Detected: {best_partial_match['event']} {best_partial_match['start_pos']:,} - {best_partial_match['end_pos']:,} ({best_partial_info['call_mb']:.1f} Mb)")
                    print(f"    Overlap: {best_partial_info['overlap']:.2f}, Size ratio: {best_partial_info['size_ratio']:.2f}")
                
                # --- Store all detected calls ---
                for idx, call in final_calls_df.iterrows():
                    all_detected_calls.append({
                        'Sample': sample_base,
                        'Truth_Type': truth_row['Type'],
                        'Truth_Fraction': truth_row['Fraction_Percent'],
                        **call.to_dict()
                    })

                # --- Classify every call as TP, partial match, germline, or FP ---
                false_positive_count = 0
                genuine_fp_count = 0
                partial_fp_count = 0
                fp_calls = []
                genuine_fp_calls = []
                
                for idx, call in final_calls_df.iterrows():
                    # Skip the matched call (TP)
                    if best_match is not None:
                        if (call['start_pos'] == best_match['start_pos'] and 
                            call['chromosome'] == best_match['chromosome']):
                            continue
                    # Skip known germline events
                    if call['chromosome'] in known_germline_chroms:
                        continue
                    
                    # Check if this FP is a partial match to truth
                    pm = partial_match_call_to_truth(call, truth_row, OVERLAP_MIN=0.20)
                    
                    false_positive_count += 1
                    fp_calls.append(call)
                    
                    if pm['partial_match']:
                        partial_fp_count += 1
                    else:
                        genuine_fp_count += 1
                        genuine_fp_calls.append(call)

                # --- PLOT: Detected (TP) ---
                if best_match is not None:
                    print(f"   ✅ DETECTED: {best_match['event']} (overlap={best_overlap:.2f})")
                    try:
                        call_df = pd.DataFrame([best_match])
                        fig = plot_BAF_AI_windows(
                            df_for_plotting, best_match['chromosome'], call_df, None, enter_thresh, chromosome_sizes
                        )
                        plot_filename = (f"{sample_base}_{best_match['chromosome']}_"
                                       f"{best_match['event']}_{best_match['span_mb']:.1f}Mb.pdf")
                        plot_path = os.path.join(PLOTS_DIR, plot_filename)
                        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                        plt.close(fig)
                    except Exception as e:
                        print(f"      ⚠️  Could not save detected plot: {e}")
                
                # --- PLOT: Undetected (truth overlay) ---
                else:
                    print(f"   ❌ NOT DETECTED")
                    try:
                        truth_chrom = truth_row['Chromosome']
                        empty_calls = pd.DataFrame()
                        fig = plot_BAF_AI_windows(
                            df_for_plotting, truth_chrom, empty_calls, None, enter_thresh, chromosome_sizes
                        )
                        add_truth_region_overlay(
                            fig, truth_row['Start_bp'], truth_row['End_bp'], truth_row['Type'],
                            truth_fraction=truth_row.get('Fraction_Percent', None)
                        )
                        plot_filename = (f"{sample_base}_{truth_chrom}_"
                                       f"NOT_DETECTED_{truth_row['Type']}_{truth_row['Length_Category']}.pdf")
                        plot_path = os.path.join(PLOTS_UNDETECTED_DIR, plot_filename)
                        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                        plt.close(fig)
                    except Exception as e:
                        print(f"      ⚠️  Could not save undetected plot: {e}")
                
                # --- PLOT: False positives (non-TP, non-germline) ---
                if false_positive_count > 0:
                    fp_summary_strs = [f"{c['chromosome']}:{c['event']}:{c['span_mb']:.1f}Mb" for c in fp_calls]
                    if genuine_fp_count > 0 and partial_fp_count > 0:
                        print(f"      ⚠️  {false_positive_count} FP call(s): {genuine_fp_count} genuine + {partial_fp_count} partial match(es): {'; '.join(fp_summary_strs)}")
                    elif partial_fp_count > 0:
                        print(f"      🔶 {partial_fp_count} partial match FP(s) (boundary/type mismatch): {'; '.join(fp_summary_strs)}")
                    else:
                        print(f"      ⚠️  {genuine_fp_count} genuine false positive(s): {'; '.join(fp_summary_strs)}")
                    
                    for fp_call in fp_calls:
                        try:
                            # Check if this is a partial match
                            pm = partial_match_call_to_truth(fp_call, truth_row, OVERLAP_MIN=0.20)
                            fp_prefix = "PARTIAL" if pm['partial_match'] else "FP"
                            
                            fp_df = pd.DataFrame([fp_call])
                            fig_fp = plot_BAF_AI_windows(
                                df_for_plotting, fp_call['chromosome'], fp_df, None, enter_thresh, chromosome_sizes
                            )
                            
                            # Add truth overlay on partial matches
                            if pm['partial_match']:
                                add_truth_region_overlay(
                                    fig_fp, truth_row['Start_bp'], truth_row['End_bp'], truth_row['Type'],
                                    truth_fraction=truth_row.get('Fraction_Percent', None)
                                )
                            
                            plot_filename_fp = (f"{sample_base}_{fp_call['chromosome']}_"
                                              f"{fp_prefix}_{fp_call['event']}_{fp_call['span_mb']:.1f}Mb.pdf")
                            plot_path_fp = os.path.join(PLOTS_FP_DIR, plot_filename_fp)
                            fig_fp.savefig(plot_path_fp, dpi=300, bbox_inches='tight')
                            plt.close(fig_fp)
                        except Exception as e:
                            print(f"      ⚠️  Could not save FP plot: {e}")

                # --- RECORD TRUTH COMPARISON ---
                fp_summary_str = '; '.join([f"{c['chromosome']}:{c['event']}:{c['span_mb']:.1f}Mb" for c in fp_calls])
                
                if best_match is not None:
                    truth_comparison.append({
                        'Sample': truth_row['Sample'],
                        'Truth_Type': truth_row['Type'],
                        'Truth_Fraction': truth_row['Fraction_Percent'],
                        'Truth_Chromosome': truth_row['Chromosome'],
                        'Truth_Start': truth_row['Start_bp'],
                        'Truth_End': truth_row['End_bp'],
                        'Truth_Length_Mb': (truth_row['End_bp'] - truth_row['Start_bp']) / 1e6,
                        'Truth_Geometry': truth_row['Geometry'],
                        'Detected': True,
                        'Status': 'TRUE_POSITIVE',
                        'Detected_Type': best_match['event'],
                        'Detected_Chromosome': best_match['chromosome'],
                        'Detected_Start': best_match['start_pos'],
                        'Detected_End': best_match['end_pos'],
                        'Detected_Length_Mb': best_match['span_mb'],
                        'Detected_Cell_Fraction': best_match.get('cell_fraction_BAF', np.nan),
                        'Overlap': best_overlap,
                        'False_Positives': false_positive_count,
                        'Genuine_FPs': genuine_fp_count,
                        'Partial_Match_FPs': partial_fp_count,
                        'Total_Calls': len(final_calls_df),
                        'FP_Summary': fp_summary_str,
                        'Partial_Match_Reason': None
                    })
                elif best_partial_info is not None:
                    # Detected the event but with wrong type or boundary — PARTIAL MATCH
                    truth_comparison.append({
                        'Sample': truth_row['Sample'],
                        'Truth_Type': truth_row['Type'],
                        'Truth_Fraction': truth_row['Fraction_Percent'],
                        'Truth_Chromosome': truth_row['Chromosome'],
                        'Truth_Start': truth_row['Start_bp'],
                        'Truth_End': truth_row['End_bp'],
                        'Truth_Length_Mb': (truth_row['End_bp'] - truth_row['Start_bp']) / 1e6,
                        'Truth_Geometry': truth_row['Geometry'],
                        'Detected': False,
                        'Status': f"PARTIAL_MATCH_{best_partial_info['reason']}",
                        'Detected_Type': best_partial_match['event'],
                        'Detected_Chromosome': best_partial_match['chromosome'],
                        'Detected_Start': best_partial_match['start_pos'],
                        'Detected_End': best_partial_match['end_pos'],
                        'Detected_Length_Mb': best_partial_match['span_mb'],
                        'Detected_Cell_Fraction': best_partial_match.get('cell_fraction_BAF', np.nan),
                        'Overlap': best_partial_info['overlap'],
                        'False_Positives': false_positive_count,
                        'Genuine_FPs': genuine_fp_count,
                        'Partial_Match_FPs': partial_fp_count,
                        'Total_Calls': len(final_calls_df),
                        'FP_Summary': fp_summary_str,
                        'Partial_Match_Reason': best_partial_info['reason']
                    })
                else:
                    best_fp = None
                    if fp_calls:
                        fp_df_all = pd.DataFrame(fp_calls)
                        fp_same = fp_df_all[fp_df_all['chromosome'] == truth_row['Chromosome']]
                        if len(fp_same) > 0:
                            best_fp = fp_same.loc[fp_same['span_mb'].idxmax()]
                        else:
                            best_fp = fp_df_all.loc[fp_df_all['span_mb'].idxmax()]
                    
                    truth_comparison.append({
                        'Sample': truth_row['Sample'],
                        'Truth_Type': truth_row['Type'],
                        'Truth_Fraction': truth_row['Fraction_Percent'],
                        'Truth_Chromosome': truth_row['Chromosome'],
                        'Truth_Start': truth_row['Start_bp'],
                        'Truth_End': truth_row['End_bp'],
                        'Truth_Length_Mb': (truth_row['End_bp'] - truth_row['Start_bp']) / 1e6,
                        'Truth_Geometry': truth_row['Geometry'],
                        'Detected': False,
                        'Status': 'FALSE_NEGATIVE_WITH_FP' if false_positive_count > 0 else 'FALSE_NEGATIVE',
                        'Detected_Type': best_fp['event'] if best_fp is not None else None,
                        'Detected_Chromosome': best_fp['chromosome'] if best_fp is not None else None,
                        'Detected_Start': best_fp['start_pos'] if best_fp is not None else None,
                        'Detected_End': best_fp['end_pos'] if best_fp is not None else None,
                        'Detected_Length_Mb': best_fp['span_mb'] if best_fp is not None else None,
                        'Detected_Cell_Fraction': best_fp.get('cell_fraction_BAF', np.nan) if best_fp is not None else np.nan,
                        'Overlap': 0.0,
                        'False_Positives': false_positive_count,
                        'Genuine_FPs': genuine_fp_count,
                        'Partial_Match_FPs': partial_fp_count,
                        'Total_Calls': len(final_calls_df),
                        'FP_Summary': fp_summary_str,
                        'Partial_Match_Reason': None
                    })
            
            else:
                # === NO CALLS AT ALL ===
                print(f"   ❌ NOT DETECTED")
                
                try:
                    truth_chrom = truth_row['Chromosome']
                    empty_calls = pd.DataFrame()
                    fig = plot_BAF_AI_windows(
                        df_for_plotting, truth_chrom, empty_calls, None, enter_thresh, chromosome_sizes
                    )
                    add_truth_region_overlay(
                        fig, truth_row['Start_bp'], truth_row['End_bp'], truth_row['Type'],
                        truth_fraction=truth_row.get('Fraction_Percent', None)
                    )
                    plot_filename = (f"{sample_base}_{truth_chrom}_"
                                   f"NOT_DETECTED_{truth_row['Type']}_{truth_row['Length_Category']}.pdf")
                    plot_path = os.path.join(PLOTS_UNDETECTED_DIR, plot_filename)
                    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                    plt.close(fig)
                except Exception as e:
                    print(f"      ⚠️  Could not save undetected plot: {e}")

                truth_comparison.append({
                    'Sample': truth_row['Sample'],
                    'Truth_Type': truth_row['Type'],
                    'Truth_Fraction': truth_row['Fraction_Percent'],
                    'Truth_Chromosome': truth_row['Chromosome'],
                    'Truth_Start': truth_row['Start_bp'],
                    'Truth_End': truth_row['End_bp'],
                    'Truth_Length_Mb': (truth_row['End_bp'] - truth_row['Start_bp']) / 1e6,
                    'Truth_Geometry': truth_row['Geometry'],
                    'Detected': False,
                    'Status': 'FALSE_NEGATIVE',
                    'Detected_Type': None,
                    'Detected_Chromosome': None,
                    'Detected_Start': None,
                    'Detected_End': None,
                    'Detected_Length_Mb': None,
                    'Detected_Cell_Fraction': None,
                    'Overlap': None,
                    'False_Positives': 0,
                    'Genuine_FPs': 0,
                    'Partial_Match_FPs': 0,
                    'Total_Calls': 0,
                    'FP_Summary': '',
                    'Partial_Match_Reason': None
                })
        
        except Exception as e:
            print(f"   ❌ ERROR: {e}")
            truth_comparison.append({
                'Sample': truth_row['Sample'],
                'Truth_Type': truth_row['Type'],
                'Truth_Fraction': truth_row['Fraction_Percent'],
                'Truth_Chromosome': truth_row['Chromosome'],
                'Truth_Start': truth_row['Start_bp'],
                'Truth_End': truth_row['End_bp'],
                'Truth_Length_Mb': (truth_row['End_bp'] - truth_row['Start_bp']) / 1e6,
                'Truth_Geometry': truth_row['Geometry'],
                'Detected': False,
                'Status': 'ERROR',
                'Detected_Type': None,
                'Detected_Chromosome': None,
                'Detected_Start': None,
                'Detected_End': None,
                'Detected_Length_Mb': None,
                'Detected_Cell_Fraction': None,
                'Overlap': None,
                'False_Positives': None,
                'Genuine_FPs': None,
                'Partial_Match_FPs': None,
                'Total_Calls': None,
                'FP_Summary': '',
                'Partial_Match_Reason': None,
                'Error': str(e)
            })
            
        processed_count += 1

        # === PERIODIC PROGRESS SUMMARY ===
        if processed_count % progress_interval == 0 or i == len(snp_files):
            elapsed = time.time() - start_time
            detected_so_far = sum(1 for tc in truth_comparison if tc['Detected'])
            total_so_far = len(truth_comparison)
            
            if total_so_far > 0:
                success_rate = detected_so_far / total_so_far * 100
            else:
                success_rate = 0
            
            samples_per_sec = i / elapsed if elapsed > 0 else 0
            remaining_samples = len(snp_files) - i
            est_remaining_sec = remaining_samples / samples_per_sec if samples_per_sec > 0 else 0
            
            print("\n" + "="*60)
            print(f"📊 PROGRESS: {processed_count}/{len(snp_files)} samples processed ({processed_count/len(snp_files)*100:.1f}%)")
            print(f"✅ Detected: {detected_so_far}/{total_so_far} events ({success_rate:.1f}%)")
            print(f"⏱️  Elapsed: {elapsed/60:.1f} min | Est. remaining: {est_remaining_sec/60:.1f} min")
            print("="*60 + "\n")
        # ==================================
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    if all_detected_calls:
        df_all_calls = pd.DataFrame(all_detected_calls)
        out_path1 = os.path.join(OUTPUT_DIR, 'all_detected_mCAs.csv')
        df_all_calls.to_csv(out_path1, index=False)
        print(f"\n✅ Saved all detected calls ({len(df_all_calls)} calls): {out_path1}")
    else:
        print("\n⚠️  No calls detected in any sample")
    
    if truth_comparison:
        df_truth = pd.DataFrame(truth_comparison)
        out_path2 = os.path.join(OUTPUT_DIR, 'truth_comparison.csv')
        df_truth.to_csv(out_path2, index=False)
        print(f"✅ Saved truth comparison ({len(df_truth)} rows): {out_path2}")
        
        print("\n" + "="*80)
        print("SUMMARY STATISTICS")
        print("="*80)
        
        detected_count = df_truth['Detected'].sum()
        total_count = len(df_truth)
        partial_count = df_truth['Status'].str.startswith('PARTIAL_MATCH').sum()
        fn_count = df_truth['Status'].isin(['FALSE_NEGATIVE', 'FALSE_NEGATIVE_WITH_FP']).sum()
        
        print(f"\nOverall: {detected_count}/{total_count} events detected ({detected_count/total_count*100:.1f}%)")
        print(f"  True Positives:  {detected_count}")
        print(f"  Partial Matches: {partial_count} (detected but wrong type/size)")
        print(f"  False Negatives: {fn_count}")
        if partial_count > 0:
            print(f"  → Including partial matches as TP: {(detected_count+partial_count)}/{total_count} ({(detected_count+partial_count)/total_count*100:.1f}%)")
        
        # Partial match breakdown
        if partial_count > 0:
            print("\nPartial Match Breakdown:")
            for reason in df_truth[df_truth['Status'].str.startswith('PARTIAL_MATCH')]['Partial_Match_Reason'].value_counts().items():
                print(f"  {reason[0]:>20}: {reason[1]}")
        
        print("\nBy Cell Fraction:")
        for frac in sorted(df_truth['Truth_Fraction'].unique()):
            frac_df = df_truth[df_truth['Truth_Fraction'] == frac]
            detected = frac_df['Detected'].sum()
            partial = frac_df['Status'].str.startswith('PARTIAL_MATCH').sum()
            total = len(frac_df)
            extra = f" (+{partial} partial)" if partial > 0 else ""
            print(f"  {frac*100:>5.1f}%: {detected:>3}/{total:<3} detected ({detected/total*100:>5.1f}%){extra}")
        
        print("\nBy Event Type:")
        for event_type in sorted(df_truth['Truth_Type'].unique()):
            type_df = df_truth[df_truth['Truth_Type'] == event_type]
            detected = type_df['Detected'].sum()
            partial = type_df['Status'].str.startswith('PARTIAL_MATCH').sum()
            total = len(type_df)
            extra = f" (+{partial} partial)" if partial > 0 else ""
            print(f"  {event_type:>6}: {detected:>3}/{total:<3} detected ({detected/total*100:>5.1f}%){extra}")
        
        print("\nBy Geometry:")
        for geom in sorted(df_truth['Truth_Geometry'].unique()):
            geom_df = df_truth[df_truth['Truth_Geometry'] == geom]
            detected = geom_df['Detected'].sum()
            partial = geom_df['Status'].str.startswith('PARTIAL_MATCH').sum()
            total = len(geom_df)
            extra = f" (+{partial} partial)" if partial > 0 else ""
            print(f"  {geom:>15}: {detected:>3}/{total:<3} detected ({detected/total*100:>5.1f}%){extra}")
        
        # FP summary — distinguish genuine from partial match FPs
        total_fps = df_truth['False_Positives'].sum()
        samples_with_fps = (df_truth['False_Positives'] > 0).sum()
        
        genuine_fps_col = df_truth['Genuine_FPs'].dropna()
        partial_fps_col = df_truth['Partial_Match_FPs'].dropna()
        total_genuine_fps = int(genuine_fps_col.sum())
        total_partial_fps = int(partial_fps_col.sum())
        
        print(f"\nFalse Positives: {total_fps} total across {samples_with_fps} samples")
        print(f"  Genuine FPs (wrong chromosome / no truth overlap): {total_genuine_fps}")
        print(f"  Partial match FPs (correct region, wrong type/size): {total_partial_fps}")
    else:
        print("\n⚠️  No results to save")
    
    print("\n" + "="*80)
    print("COMPLETE")
    print("="*80)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run mCA caller on simulated samples")
    parser.add_argument("simulated_dir", help="Path to simulated samples directory")
    parser.add_argument("output_dir", help="Path to output directory")
    args = parser.parse_args()

    SIMULATED_DIR = args.simulated_dir
    MANIFEST_FILE = os.path.join(SIMULATED_DIR, "simulation_manifest.csv")
    OUTPUT_DIR = args.output_dir
    PLOTS_DIR = os.path.join(OUTPUT_DIR, "Detected_mCA_plots")
    PLOTS_UNDETECTED_DIR = os.path.join(OUTPUT_DIR, "Undetected_mCA_plots")
    PLOTS_FP_DIR = os.path.join(OUTPUT_DIR, "False_positive_mCA_plots")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)
    os.makedirs(PLOTS_UNDETECTED_DIR, exist_ok=True)
    os.makedirs(PLOTS_FP_DIR, exist_ok=True)

    process_simulated_samples()
