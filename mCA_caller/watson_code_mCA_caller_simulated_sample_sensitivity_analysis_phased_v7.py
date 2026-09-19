#!/usr/bin/env python3
"""
Phased mCA Detection: Complete Sensitivity Analysis

Reads the phased known-region results CSV and produces:
  1.  sensitivity_analysis_phased_all_sizes.pdf         (6-panel: overall, event, size, CF accuracy, geometry, limits)
  2.  sensitivity_analysis_phased_5Mb_plus.pdf           (same, ≥5Mb only)
  3.  chromosome_analysis_phased_all_sizes.pdf           (bar + heatmap by chromosome)
  4.  chromosome_analysis_phased_5Mb_plus.pdf            (same, ≥5Mb only)
  5.  sensitivity_heatmaps_phased_all_sizes.pdf          (event type × CF and size × CF heatmaps)
  6.  sensitivity_heatmaps_phased_5Mb_plus.pdf           (same, ≥5Mb only)
  7.  sensitivity_comparison_unphased_vs_phased.pdf      (3-panel overlay if truth_comparison.csv provided)
  8.  heatmap_comparison_unphased_vs_phased.pdf          (side-by-side size×CF heatmaps)
  9.  pass_usage_diagnostic.pdf                          (three-pass filter usage by event type and size)
  10. fpr_analysis_summary.pdf                           (4-panel: n_hets, event type, size, p-value distribution)
  11. sensitivity_vs_fpr_by_event.pdf                    (sensitivity + FPR overlay by event type)
  12. roc_curves_by_event_and_cf.pdf                     (ROC curves at selected CFs)

Requires:
  - sensitivity_results_all_ref_cfs.csv  (longitudinal pipeline output)
  - Within_family_fpr_results.csv        (within-family FPR null results)
  - truth_comparison.csv                 (optional, for unphased comparison figures)
  - phased_per_sample_results.csv        (optional, for ground-truth phased figures and
                                          three-way unphased vs ground-truth vs longitudinal
                                          comparison; useful to verify ref=100% longitudinal
                                          approximates perfect phasing)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Helvetica'
import os
import sys
import argparse

# === COLOR DEFINITIONS ===
c0 = (0.76, 0.76, 0.76)
c1 = (1.00, 0.18, 0.33); c2 = (1.00, 0.23, 0.19)
c3 = (1.00, 0.58, 0.00); c4 = (1.00, 0.80, 0.00)
c5 = (0.30, 0.85, 0.39); c6 = (0.35, 0.78, 0.98)
c7 = (0.20, 0.67, 0.86); c8 = (0.00, 0.48, 1.00)
c9 = (0.35, 0.34, 0.84); c10 = (0.00, 0.31, 0.57)

orange1='#feedde'; orange2='#fdbe85'; orange3='#fd8d3c'; orange4='#e6550d'; orange5='#a63603'
blue1='#eff3ff'; blue2='#bdd7e7'; blue3='#6baed6'; blue4='#3182bd'; blue5='#08519c'
green1='#edf8e9'; green2='#bae4b3'; green3='#74c476'; green4='#31a354'; green5='#006d2c'
grey1='#f7f7f7'; grey2='#cccccc'; grey3='#969696'; grey4='#636363'; grey5='#252525'
purple1='#f2f0f7'; purple2='#cbc9e2'; purple3='#9e9ac8'; purple4='#756bb1'; purple5='#54278f'
red1='#fee5d9'; red2='#fcae91'; red3='#fb6a4a'; red4='#de2d26'; red5='#a50f15'

# ═══════════════════════════════════════════════
# CONFIGURATION — defaults (overridden by CLI args)
# ═══════════════════════════════════════════════

DEFAULT_PHASED_CSV = 'CNV_panel_simulated_Samples_Feb2026/mCA_caller_simulation_results_TEST_set_v13_caller_all_cell_fractions/phased_per_sample_results.csv'
DEFAULT_LONGITUDINAL_CSV = 'CNV_panel_simulated_Samples_Feb2026/mCA_caller_simulation_results_TEST_set_v13_caller_all_cell_fractions/longitudinal_phased_results.csv'
DEFAULT_UNPHASED_CSV = 'CNV_panel_simulated_Samples_Feb2026/mCA_caller_simulation_results_TEST_set_v13_caller_all_cell_fractions/truth_comparison.csv'
DEFAULT_OUTPUT_DIR = 'CNV_panel_simulated_Samples_Feb2026/mCA_caller_simulation_results_TEST_set_v13_caller_all_cell_fractions/Sensitivity_analysis_phased'

# Chromosome ordering
CHROM_ORDER = [f'chr{i}' for i in range(1, 23)] + ['chrX']
CHROM_LABELS = [str(i) for i in range(1, 23)] + ['X']

# Size categories
ALL_SIZES = ['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
SIZES_5PLUS = ['5Mb', '10Mb', '20Mb', 'WholeArm']

# Shared styling
EVENT_COLORS = {'CN-LOH': '#2E86AB', 'GAIN': '#F77F00', 'LOSS': '#06A77D'}
SIZE_COLORS = {
    '2Mb': '#E63946', '3Mb': '#F77F00', '5Mb': '#FCBF49',
    '10Mb': '#06A77D', '20Mb': '#2E86AB', 'WholeArm': '#7B2D8E'
}
GEOM_COLORS = {'Interstitial': '#2E86AB', 'Telomeric': '#F77F00', 'Whole_Arm': '#06A77D'}

# ═══════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════

def sensitivity_by_group(df, group_col, cf_col='cf_true', sig_col='significant_raw'):
    """Compute sensitivity for each group × CF combination."""
    result = {}
    for grp in sorted(df[group_col].unique()):
        grp_df = df[df[group_col] == grp]
        sens = {}
        for cf in sorted(grp_df[cf_col].unique()):
            subset = grp_df[grp_df[cf_col] == cf]
            n = len(subset)
            n_sig = int(subset[sig_col].sum())
            sens[cf] = {'pct': n_sig / n * 100 if n > 0 else 0, 'det': n_sig, 'n': n}
        result[grp] = sens
    return result

def find_50pct_threshold(sens_dict):
    """Find CF (as %) at which sensitivity first reaches 50%."""
    for cf in sorted(sens_dict.keys()):
        if sens_dict[cf]['pct'] >= 50:
            return cf * 100
    return None

def load_unphased_sensitivity(csv_path):
    """Load unphased results and compute sensitivity by event and size."""
    tc = pd.read_csv(csv_path)
    event_map = {'CNLOH': 'CN-LOH', 'GAIN': 'GAIN', 'LOSS': 'LOSS'}
    tc['event'] = tc['Truth_Type'].map(event_map)
    tc['size_cat'] = pd.cut(
        tc['Truth_Length_Mb'],
        bins=[0, 2.5, 3.5, 5.5, 10.5, 20.5, 200],
        labels=['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    )

    by_event = {}
    for evt in ['CN-LOH', 'GAIN', 'LOSS']:
        evt_df = tc[tc['event'] == evt]
        sens = {}
        for cf in sorted(evt_df['Truth_Fraction'].unique()):
            subset = evt_df[evt_df['Truth_Fraction'] == cf]
            n = len(subset)
            n_tp = (subset['Status'] == 'TRUE_POSITIVE').sum()
            sens[cf] = n_tp / n * 100 if n > 0 else 0
        by_event[evt] = sens

    by_size = {}
    for size in ALL_SIZES:
        size_df = tc[tc['size_cat'] == size]
        sens = {}
        for cf in sorted(size_df['Truth_Fraction'].unique()):
            subset = size_df[size_df['Truth_Fraction'] == cf]
            n = len(subset)
            n_tp = (subset['Status'] == 'TRUE_POSITIVE').sum()
            sens[cf] = n_tp / n * 100 if n > 0 else 0
        by_size[size] = sens

    return by_event, by_size

def load_and_normalise(csv_path, longitudinal=False, ref_cf=None):
    """
    Load a phased results CSV and normalise columns so all plotting
    functions can use the same column names:
        cf_true, significant_raw, event, size, chromosome, cf_estimate, geometry

    Parameters
    ----------
    csv_path : str - path to CSV
    longitudinal : bool - if True, treat as longitudinal pipeline output
    ref_cf : float or None - for longitudinal data, filter to this reference CF
    """
    df = pd.read_csv(csv_path)

    if longitudinal:
        # Filter to tested rows only
        if 'status' in df.columns:
            df = df[df['status'] == 'tested'].copy()

        # Filter to requested reference CF
        if ref_cf is not None and 'ref_cf' in df.columns:
            df = df[np.isclose(df['ref_cf'], ref_cf)].copy()

        # Rename columns to match per-sample format
        rename_map = {}
        if 'target_cf' in df.columns and 'cf_true' not in df.columns:
            rename_map['target_cf'] = 'cf_true'
        if 'significant' in df.columns and 'significant_raw' not in df.columns:
            rename_map['significant'] = 'significant_raw'
        if 'n_hets' in df.columns and 'n_hets_inside' not in df.columns:
            rename_map['n_hets'] = 'n_hets_inside'
        if rename_map:
            df = df.rename(columns=rename_map)

    # Ensure required columns exist
    required = ['cf_true', 'significant_raw', 'event', 'size', 'chromosome']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns after normalisation: {missing}")

    return df

# ═══════════════════════════════════════════════
# FIGURE: 6-panel sensitivity analysis
# ═══════════════════════════════════════════════

def plot_sensitivity_6panel(df, title_suffix, output_name, fpr_df=None, ref_cf=None):
    """
    Layout (mirrors unphased figure):
      Row 1 (4 panels): Size | Overall | Event Type | Geometry
      Row 2 (2 panels): CF Estimation Accuracy | Chromosome × CF heatmap
    """
    fig = plt.figure(figsize=(26, 12))
    fig.suptitle(f'Phased Detection: Sensitivity Analysis ({title_suffix})',
                 fontsize=16, fontweight='bold', y=0.98)

    gs_top = fig.add_gridspec(1, 4, left=0.05, right=0.98, top=0.88, bottom=0.54,
                               wspace=0.35)
    gs_bot = fig.add_gridspec(1, 2, left=0.05, right=0.98, top=0.46, bottom=0.12,
                               wspace=0.35, width_ratios=[1.12, 3.7])

    ax_size  = fig.add_subplot(gs_top[0])
    ax_over  = fig.add_subplot(gs_top[1])
    ax_evt   = fig.add_subplot(gs_top[2])
    ax_geom  = fig.add_subplot(gs_top[3])
    ax_cf    = fig.add_subplot(gs_bot[0])
    ax_heat  = fig.add_subplot(gs_bot[1])

    cf_values = sorted(df['cf_true'].unique())

    def fmt_log_axis(ax):
        ax.set_xscale('log')
        ax.set_xticks([0.1, 0.5, 1, 5, 10, 50, 100])
        ax.get_xaxis().set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f'{x:g}'))
        ax.set_xlim(0.08, 130)

    # ── Panel a: Sensitivity by Event Size ──
    size_sens = sensitivity_by_group(df, 'size')
    size_order = [s for s in ALL_SIZES if s in size_sens]
    for size in size_order:
        cfs = sorted(size_sens[size].keys())
        ax_size.plot([c * 100 for c in cfs], [size_sens[size][c]['pct'] for c in cfs],
                     'o-', color=SIZE_COLORS.get(size, 'grey'), linewidth=2.5, markersize=8, label=size)
    ax_size.axhline(50, color='grey', linestyle='--', alpha=0.5)
    ax_size.legend(fontsize=9, loc='lower right', title='Event Size', title_fontsize=9)
    ax_size.set_xlabel('Cell Fraction (%)', fontsize=11)
    ax_size.set_ylabel('Sensitivity (%)', fontsize=11)
    ax_size.set_title('Sensitivity by Event Size', fontsize=12, fontweight='bold')
    fmt_log_axis(ax_size)
    ax_size.set_ylim(-2, 105)
    ax_size.grid(alpha=0.2, which="both")

    # ── Panel b: Overall Sensitivity ──
    event_sens = sensitivity_by_group(df, 'event')
    stats = []
    for cf in cf_values:
        subset = df[df['cf_true'] == cf]
        n = len(subset)
        det = int(subset['significant_raw'].sum())
        stats.append((cf * 100, det / n * 100 if n > 0 else 0, det, n))
    ax_over.plot([s[0] for s in stats], [s[1] for s in stats],
                 'o-', color='#2166AC', linewidth=2.5, markersize=9, zorder=3)
    for cp, s, d, n in stats:
        if s > 0 or cp >= 1:
            ax_over.annotate(f'{d}/{n}', (cp, s), textcoords='offset points',
                             xytext=(0, 12), fontsize=7, ha='center', fontweight='bold')
    ax_over.axhline(50, color='grey', linestyle='--', alpha=0.5)
    ax_over.set_xlabel('Cell Fraction (%)', fontsize=11)
    ax_over.set_ylabel('Sensitivity (%)', fontsize=11)
    ax_over.set_title('Overall Sensitivity vs Cell Fraction', fontsize=12, fontweight='bold')
    fmt_log_axis(ax_over)
    ax_over.set_ylim(-2, 105)
    ax_over.grid(alpha=0.2, which="both")
    ax_over.text(0.02, 0.95, '50% sensitivity', transform=ax_over.transAxes,
                 fontsize=9, color='grey', va='top')

    # ── Panel c: Sensitivity by Event Type ──
    for evt in ['CN-LOH', 'GAIN', 'LOSS']:
        if evt not in event_sens:
            continue
        cfs = sorted(event_sens[evt].keys())
        ax_evt.plot([c * 100 for c in cfs], [event_sens[evt][c]['pct'] for c in cfs],
                    'o-', color=EVENT_COLORS[evt], linewidth=2.5, markersize=8, label=evt)
    ax_evt.axhline(50, color='grey', linestyle='--', alpha=0.5)
    ax_evt.legend(fontsize=10, loc='lower right')
    ax_evt.set_xlabel('Cell Fraction (%)', fontsize=11)
    ax_evt.set_ylabel('Sensitivity (%)', fontsize=11)
    ax_evt.set_title('Sensitivity by Event Type', fontsize=12, fontweight='bold')
    fmt_log_axis(ax_evt)
    ax_evt.set_ylim(-2, 105)
    ax_evt.grid(alpha=0.2, which="both")

    # ── Panel d: Sensitivity by Geometry ──
    geom_sens = sensitivity_by_group(df, 'geometry')
    for geom in ['Interstitial', 'Telomeric', 'Whole_Arm']:
        if geom not in geom_sens:
            continue
        cfs = sorted(geom_sens[geom].keys())
        ax_geom.plot([c * 100 for c in cfs], [geom_sens[geom][c]['pct'] for c in cfs],
                     'o-', color=GEOM_COLORS.get(geom, 'grey'), linewidth=2.5, markersize=8, label=geom)
    ax_geom.axhline(50, color='grey', linestyle='--', alpha=0.5)
    ax_geom.legend(fontsize=10, loc='lower right')
    ax_geom.set_xlabel('Cell Fraction (%)', fontsize=11)
    ax_geom.set_ylabel('Sensitivity (%)', fontsize=11)
    ax_geom.set_title('Sensitivity by Geometry', fontsize=12, fontweight='bold')
    fmt_log_axis(ax_geom)
    ax_geom.set_ylim(-2, 105)
    ax_geom.grid(alpha=0.2, which="both")

    # ── Panel e: CF Estimation Accuracy — median ± IQR, no scatter ──
    sig_df = df[df['significant_raw'] == True].copy()
    if len(sig_df) > 0:
        ax_cf.plot([0.05, 120], [0.05, 120], 'k--', linewidth=1, alpha=0.5, label='Perfect')
        for evt in ['CN-LOH', 'GAIN', 'LOSS']:
            evt_df = sig_df[sig_df['event'] == evt]
            if len(evt_df) > 0:
                grp   = evt_df.groupby('cf_true')['cf_estimate']
                med   = grp.median() * 100
                q25   = grp.quantile(0.25) * 100
                q75   = grp.quantile(0.75) * 100
                x_vals = med.index * 100
                ax_cf.plot(x_vals, med.values, 'o-',
                           color=EVENT_COLORS[evt], markersize=8, linewidth=2.5,
                           label=evt, alpha=0.95, zorder=5)
                ax_cf.fill_between(x_vals, q25.values, q75.values,
                                   color=EVENT_COLORS[evt], alpha=0.18, zorder=4)
    ax_cf.legend(fontsize=8, loc='upper left')
    ax_cf.set_xlabel('True Cell Fraction (%)', fontsize=11)
    ax_cf.set_ylabel('Estimated Cell Fraction (%)', fontsize=11)
    ax_cf.set_title('CF Estimation Accuracy (median ± IQR)', fontsize=12, fontweight='bold')
    ax_cf.set_xscale('log'); ax_cf.set_yscale('log')
    ax_cf.set_xticks([0.1, 0.5, 1, 5, 10, 50, 100])
    ax_cf.set_yticks([0.1, 0.5, 1, 5, 10, 50, 100])
    ax_cf.get_xaxis().set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:g}'))
    ax_cf.get_yaxis().set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:g}'))
    ax_cf.set_xlim(0.05, 130); ax_cf.set_ylim(0.05, 130)
    ax_cf.grid(alpha=0.2, which='both')

    # ── Panel f: Sensitivity by chromosome × CF heatmap ──
    df_h = df.copy()
    df_h['chrom_clean'] = df_h['chromosome'].astype(str)
    if not df_h['chrom_clean'].str.startswith('chr').all():
        df_h['chrom_clean'] = 'chr' + df_h['chrom_clean'].str.replace('chr', '', regex=False)
    df_h = df_h[df_h['chrom_clean'].isin(CHROM_ORDER)]
    cf_heatmap = [cf for cf in sorted(df_h['cf_true'].unique()) if cf >= 0.005]
    cf_heatmap_labels = [f'{cf*100:g}%' for cf in cf_heatmap]

    matrix = np.zeros((len(cf_heatmap), len(CHROM_ORDER)))
    for j, chrom in enumerate(CHROM_ORDER):
        for i, cf in enumerate(cf_heatmap):
            subset = df_h[(df_h['chrom_clean'] == chrom) & (df_h['cf_true'] == cf)]
            n = len(subset)
            if n > 0:
                matrix[i, j] = subset['significant_raw'].sum() / n * 100

    # highest CF at top
    matrix_plot = matrix[::-1]
    cf_labels_plot = cf_heatmap_labels[::-1]

    im = ax_heat.imshow(matrix_plot, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
    for i in range(len(cf_heatmap)):
        for j in range(len(CHROM_ORDER)):
            val = matrix_plot[i, j]
            color = 'white' if val < 25 or val > 85 else 'black'
            ax_heat.text(j, i, f'{val:.0f}', ha='center', va='center',
                         fontsize=7, fontweight='bold', color=color)
    ax_heat.set_xticks(range(len(CHROM_ORDER)))
    ax_heat.set_xticklabels(CHROM_LABELS, fontsize=9)
    ax_heat.set_xlabel('Chromosome', fontsize=11)
    ax_heat.set_yticks(range(len(cf_heatmap)))
    ax_heat.set_yticklabels(cf_labels_plot, fontsize=10)
    ax_heat.set_ylabel('Cell Fraction', fontsize=11)
    cf_threshold = 10
    ax_heat.set_title(f'Sensitivity by Chromosome and Cell Fraction (≥{cf_threshold}%)',
                      fontsize=12, fontweight='bold')
    fig.colorbar(im, ax=ax_heat, fraction=0.015, pad=0.02).set_label('Sensitivity (%)', fontsize=10)

    for ax in [ax_size, ax_over, ax_evt, ax_geom, ax_cf, ax_heat]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)

    plt.savefig(os.path.join(OUTPUT_DIR, output_name + '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  ✅ {output_name}.pdf")


# ═══════════════════════════════════════════════
# FIGURE: Chromosome analysis (bar + heatmap)
# ═══════════════════════════════════════════════

def plot_chromosome_analysis(df, title_suffix, output_name):
    df = df.copy()
    df['chrom_clean'] = df['chromosome'].astype(str)
    if not df['chrom_clean'].str.startswith('chr').all():
        df['chrom_clean'] = 'chr' + df['chrom_clean'].str.replace('chr', '')
    df = df[df['chrom_clean'].isin(CHROM_ORDER)]

    cf_values = sorted(df['cf_true'].unique())
    cf_heatmap = [cf for cf in cf_values if cf >= 0.005]
    cf_heatmap_labels = [f'{cf*100:g}%' for cf in cf_heatmap]

    fig, (ax_bar, ax_heat) = plt.subplots(2, 1, figsize=(16, 10),
                                           gridspec_kw={'height_ratios': [1, 1.3]})
    fig.suptitle(f'Phased Detection: Sensitivity by Chromosome ({title_suffix})',
                 fontsize=16, fontweight='bold', y=0.98)

    # Bar chart
    bar_data = []
    for chrom in CHROM_ORDER:
        chrom_df = df[df['chrom_clean'] == chrom]
        n = len(chrom_df)
        n_sig = int(chrom_df['significant_raw'].sum()) if n > 0 else 0
        bar_data.append((n_sig / n * 100 if n > 0 else 0, n_sig, n))

    x_pos = range(len(CHROM_ORDER))
    sens_vals = [b[0] for b in bar_data]
    mean_sens = np.mean([b[0] for b in bar_data if b[2] > 0])

    ax_bar.bar(x_pos, sens_vals, color='#4C9ED9', edgecolor='white', linewidth=0.5)
    ax_bar.axhline(mean_sens, color='#D62728', linestyle='--', linewidth=1.5,
                   label=f'Mean: {mean_sens:.1f}%')
    for i, (s, n_sig, n) in enumerate(bar_data):
        if n > 0:
            ax_bar.text(i, s + 1.5, f'{n_sig}/{n}', ha='center', fontsize=7,
                       fontweight='bold', color='#333333')
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(CHROM_LABELS, fontsize=10)
    ax_bar.set_ylabel('Sensitivity (%)', fontsize=12)
    ax_bar.set_title('Sensitivity by Chromosome (All Cell Fractions)', fontsize=13, fontweight='bold')
    ax_bar.set_ylim(0, 105)
    ax_bar.legend(fontsize=11, loc='upper right')
    ax_bar.spines['top'].set_visible(False); ax_bar.spines['right'].set_visible(False)

    # Heatmap
    matrix = np.zeros((len(CHROM_ORDER), len(cf_heatmap)))
    for i, chrom in enumerate(CHROM_ORDER):
        for j, cf in enumerate(cf_heatmap):
            subset = df[(df['chrom_clean'] == chrom) & (df['cf_true'] == cf)]
            n = len(subset)
            if n > 0:
                matrix[i, j] = subset['significant_raw'].sum() / n * 100

    im = ax_heat.imshow(matrix.T, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
    for i in range(len(CHROM_ORDER)):
        for j in range(len(cf_heatmap)):
            val = matrix[i, j]
            color = 'white' if val < 25 or val > 85 else 'black'
            ax_heat.text(i, j, f'{val:.0f}', ha='center', va='center',
                        fontsize=7, fontweight='bold', color=color)

    ax_heat.set_xticks(range(len(CHROM_ORDER)))
    ax_heat.set_xticklabels(CHROM_LABELS, fontsize=10)
    ax_heat.set_xlabel('Chromosome', fontsize=12)
    ax_heat.set_yticks(range(len(cf_heatmap)))
    ax_heat.set_yticklabels(cf_heatmap_labels, fontsize=10)
    ax_heat.set_ylabel('Cell Fraction', fontsize=12)
    ax_heat.invert_yaxis()
    ax_heat.set_title('Sensitivity by Chromosome and Cell Fraction (≥0.5%)',
                      fontsize=13, fontweight='bold')
    
    for ax in [ax_bar, ax_heat]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)

    fig.colorbar(im, ax=ax_heat, fraction=0.02, pad=0.02).set_label('Sensitivity (%)', fontsize=11)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(OUTPUT_DIR, output_name + '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  ✅ {output_name}.pdf")


# ═══════════════════════════════════════════════
# FIGURE: Event type × CF and Size × CF heatmaps
# ═══════════════════════════════════════════════

def plot_sensitivity_heatmaps(df, title_suffix, output_name):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 4.5))

    cf_order = sorted(df['cf_true'].unique())
    cf_labels = [f'{cf*100:g}%' for cf in cf_order]

    # Event type × CF
    events_present = [e for e in ['CN-LOH', 'GAIN', 'LOSS'] if e in df['event'].unique()]
    matrix1 = np.zeros((len(events_present), len(cf_order)))
    for i, evt in enumerate(events_present):
        for j, cf in enumerate(cf_order):
            subset = df[(df['event'] == evt) & (df['cf_true'] == cf)]
            n = len(subset)
            if n > 0:
                matrix1[i, j] = subset['significant_raw'].sum() / n * 100

    # im1 = ax1.imshow(matrix1, cmap='RdYlGn', vmin=0, vmax=100, aspect='auto')
    im1 = ax1.imshow(matrix1, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
    for i in range(len(events_present)):
        for j in range(len(cf_order)):
            val = matrix1[i, j]
            color = 'white' if val < 25 or val > 85 else 'black'
            ax1.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=8, fontweight='bold', color=color)
    ax1.set_xticks(range(len(cf_labels)))
    ax1.set_xticklabels(cf_labels, fontsize=9, rotation=45, ha='right')
    ax1.set_yticks(range(len(events_present)))
    ax1.set_yticklabels(events_present, fontsize=11)
    ax1.set_xlabel('Cell Fraction', fontsize=11)
    ax1.set_title('Sensitivity by Event Type and Cell Fraction', fontsize=12, fontweight='bold')

    # Size × CF
    sizes_present = [s for s in ALL_SIZES if s in df['size'].unique()]
    matrix2 = np.zeros((len(sizes_present), len(cf_order)))
    for i, size in enumerate(sizes_present):
        for j, cf in enumerate(cf_order):
            subset = df[(df['size'] == size) & (df['cf_true'] == cf)]
            n = len(subset)
            if n > 0:
                matrix2[i, j] = subset['significant_raw'].sum() / n * 100

    im2 = ax2.imshow(matrix2, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
    for i in range(len(sizes_present)):
        for j in range(len(cf_order)):
            val = matrix2[i, j]
            color = 'white' if val < 25 or val > 85 else 'black'
            ax2.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=8, fontweight='bold', color=color)

    ax2.set_xticks(range(len(cf_labels)))
    ax2.set_xticklabels(cf_labels, fontsize=9, rotation=45, ha='right')
    ax2.set_yticks(range(len(sizes_present)))
    ax2.set_yticklabels(sizes_present, fontsize=11)
    ax2.set_xlabel('Cell Fraction', fontsize=11)
    ax2.set_title('Sensitivity by Event Size and Cell Fraction', fontsize=12, fontweight='bold')

    for ax in [ax1, ax2]:
        for axis in ['bottom','left', 'top', 'right']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)

    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    fig.colorbar(im2, cax=cbar_ax).set_label('Sensitivity (%)', fontsize=11)

    plt.savefig(os.path.join(OUTPUT_DIR, output_name + '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  ✅ {output_name}.pdf")


# ═══════════════════════════════════════════════
# FIGURE: Unphased vs Phased comparison
# ═══════════════════════════════════════════════

def plot_comparison_curves(phased_df, unphased_by_event, phased_label='Phased (ground truth)'):
    """3-panel: unphased vs phased sensitivity by event type."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)
    event_labels = {'CN-LOH': 'CN-LOH', 'GAIN': 'Gain', 'LOSS': 'Loss'}

    phased_sens = sensitivity_by_group(phased_df, 'event')

    for ax, evt in zip(axes, ['CN-LOH', 'GAIN', 'LOSS']):
        color = EVENT_COLORS[evt]

        # Unphased
        u_cfs = sorted(unphased_by_event[evt].keys())
        u_sens = [unphased_by_event[evt][cf] for cf in u_cfs]
        u_pct = [cf * 100 for cf in u_cfs]

        # Phased
        p_cfs = sorted(phased_sens[evt].keys())
        p_sens = [phased_sens[evt][cf]['pct'] for cf in p_cfs]
        p_pct = [cf * 100 for cf in p_cfs]

        # Fill between (common CFs)
        common = sorted(set(u_cfs) & set(p_cfs))
        common_pct = [cf * 100 for cf in common]
        ax.fill_between(common_pct,
                       [unphased_by_event[evt][cf] for cf in common],
                       [phased_sens[evt][cf]['pct'] for cf in common],
                       alpha=0.15, color=color)

        ax.plot(u_pct, u_sens, 'o-', color=color, alpha=0.4, linewidth=2, markersize=6, label='Unphased')
        ax.plot(p_pct, p_sens, 's-', color=color, linewidth=2.5, markersize=7, label=phased_label)
        ax.axhline(50, color='grey', linestyle=':', alpha=0.5)

        ax.set_title(event_labels[evt], fontsize=16, fontweight='bold', pad=10)
        ax.set_xlabel('Cell Fraction (%)', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Sensitivity (%)', fontsize=12)
        # ax.set_xlim(-1, 105); ax.set_ylim(-1, 105)
        # ax.set_xscale('symlog', linthresh=1)
        ax.set_xscale('log')
        ax.set_xticks([0.1, 0.5, 1, 5, 10, 25, 50, 100])
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.legend(fontsize=10, loc='upper left', framealpha=0.9)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)

    fig.suptitle('Detection Sensitivity: Unphased vs Phased Analysis',
                 fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_comparison_unphased_vs_phased.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ sensitivity_comparison_unphased_vs_phased.pdf")

def plot_comparison_heatmaps(phased_df, unphased_by_size, phased_title='Phased Detection (ground truth)'):
    """Side-by-side heatmaps: unphased vs phased, size × CF."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 4.5))

    cf_order = sorted(phased_df['cf_true'].unique())
    cf_labels = [f'{cf*100:g}%' for cf in cf_order]
    sizes_present = [s for s in ALL_SIZES if s in phased_df['size'].unique()]

    phased_size_sens = sensitivity_by_group(phased_df, 'size')

    for ax, (title, get_val) in zip(axes, [
        ('Unphased Detection', lambda s, c: unphased_by_size.get(s, {}).get(c, 0)),
        (phased_title,
         lambda s, c: phased_size_sens.get(s, {}).get(c, {}).get('pct', 0)
         if s in phased_size_sens and c in phased_size_sens.get(s, {}) else 0)
    ]):
        matrix = np.array([[get_val(size, cf) for cf in cf_order] for size in sizes_present])
        im = ax.imshow(matrix, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
        for i in range(len(sizes_present)):
            for j in range(len(cf_order)):
                val = matrix[i, j]
                color = 'white' if val < 30 or val > 85 else 'black'
                ax.text(j, i, f'{val:.0f}', ha='center', va='center',
                        fontsize=8, fontweight='bold', color=color)
        ax.set_xticks(range(len(cf_labels)))
        ax.set_xticklabels(cf_labels, fontsize=9, rotation=45, ha='right')
        ax.set_yticks(range(len(sizes_present)))
        ax.set_yticklabels(sizes_present, fontsize=10)
        ax.set_xlabel('Cell Fraction', fontsize=11)
        ax.set_title(title, fontsize=13, fontweight='bold', pad=10)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)

    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax).set_label('Sensitivity (%)', fontsize=11)

    plt.savefig(os.path.join(OUTPUT_DIR, 'heatmap_comparison_unphased_vs_phased.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print('  ✅ heatmap_comparison_unphased_vs_phased_vs_unphased.pdf')


# ═══════════════════════════════════════════════
# FIGURE: Three-way comparison (unphased vs ground-truth vs longitudinal)
# ═══════════════════════════════════════════════

def plot_threeway_curves(gt_df, long_df, unphased_by_event, ref_pct):
    """3-panel: unphased vs ground-truth vs longitudinal sensitivity by event type."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 6), sharey=True)

    gt_sens = sensitivity_by_group(gt_df, 'event')
    long_sens = sensitivity_by_group(long_df, 'event')

    for ax, evt in zip(axes, ['CN-LOH', 'GAIN', 'LOSS']):
        color = EVENT_COLORS[evt]

        # Unphased
        u_cfs = sorted(unphased_by_event[evt].keys())
        ax.plot([c * 100 for c in u_cfs],
                [unphased_by_event[evt][c] for c in u_cfs],
                'o-', color=color, alpha=0.3, linewidth=2, markersize=6, label='Unphased')

        # Ground-truth phased
        if evt in gt_sens:
            g_cfs = sorted(gt_sens[evt].keys())
            ax.plot([c * 100 for c in g_cfs],
                    [gt_sens[evt][c]['pct'] for c in g_cfs],
                    's-', color=color, linewidth=2.5, markersize=7, label='Phased (ground truth)')

        # Longitudinal
        if evt in long_sens:
            l_cfs = sorted(long_sens[evt].keys())
            ax.plot([c * 100 for c in l_cfs],
                    [long_sens[evt][c]['pct'] for c in l_cfs],
                    'D--', color=color, alpha=0.7, linewidth=2, markersize=6,
                    label=f'Longitudinal (index={ref_pct}%)')

        ax.axhline(50, color='grey', linestyle=':', alpha=0.5)
        ax.set_title(evt, fontsize=15, fontweight='bold')
        ax.set_xlabel('Cell Fraction (%)', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Sensitivity (%)', fontsize=12)
        ax.set_xscale('log')
        ax.set_xticks([0.1, 0.5, 1, 5, 10, 25, 50, 100])
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.set_ylim(-2, 105)
        ax.legend(fontsize=9, loc='upper left', framealpha=0.9)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.grid(axis='y', alpha=0.2)

    fig.suptitle('Detection Sensitivity: Unphased vs Ground-Truth Phased vs Longitudinal',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_threeway_comparison_by_event.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ sensitivity_threeway_comparison_by_event.pdf")

def plot_threeway_overall(gt_df, long_df, unphased_by_event, ref_pct):
    """Single panel: overall sensitivity for all three methods."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Unphased overall: average across events
    all_u_cfs = set()
    for evt in unphased_by_event:
        all_u_cfs.update(unphased_by_event[evt].keys())
    u_cfs = sorted(all_u_cfs)
    u_sens = []
    for cf in u_cfs:
        vals = [unphased_by_event[evt].get(cf, 0) for evt in unphased_by_event]
        u_sens.append(np.mean(vals))
    ax.plot([c * 100 for c in u_cfs], u_sens,
            'o-', color=grey4, linewidth=2, markersize=6, alpha=0.5, label='Unphased')

    # Ground-truth
    gt_cfs = sorted(gt_df['cf_true'].unique())
    gt_sens = []
    for cf in gt_cfs:
        sub = gt_df[gt_df['cf_true'] == cf]
        gt_sens.append(sub['significant_raw'].sum() / len(sub) * 100 if len(sub) > 0 else 0)
    ax.plot([c * 100 for c in gt_cfs], gt_sens,
            's-', color=blue5, linewidth=2.5, markersize=7, label='Phased (ground truth)')

    # Longitudinal
    l_cfs = sorted(long_df['cf_true'].unique())
    l_sens = []
    for cf in l_cfs:
        sub = long_df[long_df['cf_true'] == cf]
        l_sens.append(sub['significant_raw'].sum() / len(sub) * 100 if len(sub) > 0 else 0)
    ax.plot([c * 100 for c in l_cfs], l_sens,
            'D--', color=orange4, linewidth=2, markersize=6, label=f'Longitudinal (index={ref_pct}%)')

    ax.axhline(50, color='grey', linestyle=':', alpha=0.5)
    ax.set_xlabel('Cell Fraction (%)', fontsize=12)
    ax.set_ylabel('Sensitivity (%)', fontsize=12)
    ax.set_title('Overall Detection Sensitivity: Three-Way Comparison',
                 fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_xticks([0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 75])
    ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    ax.set_ylim(-2, 105)
    ax.legend(fontsize=11, loc='lower right', framealpha=0.9)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    for sp in ['bottom', 'left']:
        ax.spines[sp].set_linewidth(1.5)
        ax.spines[sp].set_color(grey3)
    ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
    ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
    ax.grid(axis='y', alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_threeway_comparison_overall.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ sensitivity_threeway_comparison_overall.pdf")

def plot_threeway_heatmaps(gt_df, long_df, unphased_by_size, ref_pct):
    """Side-by-side-by-side heatmaps: unphased vs ground-truth vs longitudinal."""
    fig, axes = plt.subplots(1, 3, figsize=(22, 4.5))

    gt_size_sens = sensitivity_by_group(gt_df, 'size')
    long_size_sens = sensitivity_by_group(long_df, 'size')

    # Use union of CFs present in both datasets
    cf_order_gt = sorted(gt_df['cf_true'].unique())
    cf_order_long = sorted(long_df['cf_true'].unique())
    cf_order = sorted(set(cf_order_gt) | set(cf_order_long))
    cf_labels = [f'{cf*100:g}%' for cf in cf_order]
    sizes_present = [s for s in ALL_SIZES if s in gt_df['size'].unique() or s in long_df['size'].unique()]

    panels = [
        ('Unphased', lambda s, c: unphased_by_size.get(s, {}).get(c, 0)),
        ('Phased (ground truth)',
         lambda s, c: gt_size_sens.get(s, {}).get(c, {}).get('pct', 0)
         if s in gt_size_sens and c in gt_size_sens.get(s, {}) else 0),
        (f'Longitudinal (index={ref_pct}%)',
         lambda s, c: long_size_sens.get(s, {}).get(c, {}).get('pct', 0)
         if s in long_size_sens and c in long_size_sens.get(s, {}) else 0),
    ]

    for ax, (title, get_val) in zip(axes, panels):
        matrix = np.array([[get_val(size, cf) for cf in cf_order] for size in sizes_present])
        im = ax.imshow(matrix, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
        for i in range(len(sizes_present)):
            for j in range(len(cf_order)):
                val = matrix[i, j]
                color = 'white' if val < 30 or val > 85 else 'black'
                ax.text(j, i, f'{val:.0f}', ha='center', va='center',
                        fontsize=7, fontweight='bold', color=color)
        ax.set_xticks(range(len(cf_labels)))
        ax.set_xticklabels(cf_labels, fontsize=8, rotation=45, ha='right')
        ax.set_yticks(range(len(sizes_present)))
        ax.set_yticklabels(sizes_present if ax == axes[0] else [], fontsize=10)
        ax.set_xlabel('Cell Fraction', fontsize=10)
        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)

    fig.suptitle('Detection Sensitivity by Event Size: Three-Way Comparison',
                 fontsize=15, fontweight='bold', y=1.03)
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.012, 0.7])
    fig.colorbar(im, cax=cbar_ax).set_label('Sensitivity (%)', fontsize=10)

    plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_threeway_comparison_heatmaps.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ sensitivity_threeway_comparison_heatmaps.pdf")

# ═══════════════════════════════════════════════
# FIGURE: Multi-reference CF comparison
# ═══════════════════════════════════════════════

REF_CF_COLORS = {1.0: blue5, 0.75: blue3, 0.5: orange3, 0.25: red4}
REF_CF_MARKERS = {1.0: 's', 0.75: 'D', 0.5: '^', 0.25: 'o'}

def plot_longitudinal_ref_comparison_overall(ref_dfs):
    """Overall sensitivity curve for each reference CF, overlaid."""
    fig, ax = plt.subplots(figsize=(8, 5.5))

    for ref_cf, df in sorted(ref_dfs.items(), reverse=True):
        cf_values = sorted(df['cf_true'].unique())
        sens = []
        for cf in cf_values:
            subset = df[df['cf_true'] == cf]
            n = len(subset)
            det = int(subset['significant_raw'].sum())
            sens.append(det / n * 100 if n > 0 else 0)

        label = f'index={ref_cf*100:.0f}% (N={len(df[df["cf_true"]==cf_values[0]])})'
        ax.plot([c * 100 for c in cf_values], sens,
                marker=REF_CF_MARKERS.get(ref_cf, 'o'), linestyle='-',
                color=REF_CF_COLORS.get(ref_cf, 'grey'),
                linewidth=2, markersize=6, label=label)

    ax.axhline(50, color='grey', linestyle='--', alpha=0.5)
    ax.set_xlabel('Target Cell Fraction (%)', fontsize=12)
    ax.set_ylabel('Sensitivity (%)', fontsize=12)
    ax.set_title('Longitudinal Detection: Effect of Index Sample Cell Fraction',
                 fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_xticks([0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 75])
    ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    ax.set_ylim(-2, 105)
    ax.legend(fontsize=10, loc='lower right')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(1.5)
        ax.spines[axis].set_color(grey3)
    ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
    ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
    ax.grid(axis='y', alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'longitudinal_ref_comparison_overall.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ longitudinal_ref_comparison_overall.pdf")

def plot_longitudinal_ref_comparison_by_event(ref_dfs):
    """3-panel: sensitivity by event type, one line per reference CF."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)

    for ax, evt in zip(axes, ['CN-LOH', 'GAIN', 'LOSS']):
        for ref_cf, df in sorted(ref_dfs.items(), reverse=True):
            evt_df = df[df['event'] == evt]
            if len(evt_df) == 0:
                continue
            cf_values = sorted(evt_df['cf_true'].unique())
            sens = []
            for cf in cf_values:
                subset = evt_df[evt_df['cf_true'] == cf]
                n = len(subset)
                det = int(subset['significant_raw'].sum())
                sens.append(det / n * 100 if n > 0 else 0)

            ax.plot([c * 100 for c in cf_values], sens,
                    marker=REF_CF_MARKERS.get(ref_cf, 'o'), linestyle='-',
                    color=REF_CF_COLORS.get(ref_cf, 'grey'),
                    linewidth=2, markersize=5, label=f'index={ref_cf*100:.0f}%')

        ax.axhline(50, color='grey', linestyle=':', alpha=0.5)
        ax.set_title(evt, fontsize=14, fontweight='bold')
        ax.set_xlabel('Target Cell Fraction (%)', fontsize=11)
        if ax == axes[0]:
            ax.set_ylabel('Sensitivity (%)', fontsize=12)
        ax.set_xscale('log')
        ax.set_xticks([0.1, 0.5, 1, 5, 10, 25, 50, 75])
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.set_ylim(-2, 105)
        ax.legend(fontsize=9, loc='lower right')
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.grid(axis='y', alpha=0.2)

    fig.suptitle('Longitudinal Detection by Event Type: Effect of Index Sample CF',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'longitudinal_ref_comparison_by_event.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ longitudinal_ref_comparison_by_event.pdf")

def plot_longitudinal_ref_comparison_by_size(ref_dfs):
    """Multi-panel: sensitivity by event size, one line per reference CF."""
    sizes = ['3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    fig, axes = plt.subplots(1, len(sizes), figsize=(4 * len(sizes), 5.5), sharey=True)

    for ax, size in zip(axes, sizes):
        for ref_cf, df in sorted(ref_dfs.items(), reverse=True):
            size_df = df[df['size'] == size]
            if len(size_df) == 0:
                continue
            cf_values = sorted(size_df['cf_true'].unique())
            sens = []
            for cf in cf_values:
                subset = size_df[size_df['cf_true'] == cf]
                n = len(subset)
                det = int(subset['significant_raw'].sum())
                sens.append(det / n * 100 if n > 0 else 0)

            ax.plot([c * 100 for c in cf_values], sens,
                    marker=REF_CF_MARKERS.get(ref_cf, 'o'), linestyle='-',
                    color=REF_CF_COLORS.get(ref_cf, 'grey'),
                    linewidth=2, markersize=5, label=f'index={ref_cf*100:.0f}%')

        ax.axhline(50, color='grey', linestyle=':', alpha=0.5)
        ax.set_title(size, fontsize=13, fontweight='bold')
        ax.set_xlabel('Target CF (%)', fontsize=10)
        if ax == axes[0]:
            ax.set_ylabel('Sensitivity (%)', fontsize=12)
        ax.set_xscale('log')
        ax.set_xticks([0.1, 1, 5, 10, 50])
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.set_ylim(-2, 105)
        ax.legend(fontsize=8, loc='lower right')
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.grid(axis='y', alpha=0.2)

    fig.suptitle('Longitudinal Detection by Event Size: Effect of Index Sample CF',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'longitudinal_ref_comparison_by_size.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ longitudinal_ref_comparison_by_size.pdf")

def plot_longitudinal_ref_heatmap_grid(ref_dfs):
    """Grid of size×CF heatmaps, one per reference CF."""
    ref_cfs_sorted = sorted(ref_dfs.keys(), reverse=True)
    n_refs = len(ref_cfs_sorted)
    fig, axes = plt.subplots(1, n_refs, figsize=(5.5 * n_refs, 4.5))
    if n_refs == 1:
        axes = [axes]

    for ax, ref_cf in zip(axes, ref_cfs_sorted):
        df = ref_dfs[ref_cf]
        cf_order = sorted(df['cf_true'].unique())
        cf_labels = [f'{cf*100:g}%' for cf in cf_order]
        sizes_present = [s for s in ALL_SIZES if s in df['size'].unique()]

        matrix = np.zeros((len(sizes_present), len(cf_order)))
        for i, size in enumerate(sizes_present):
            for j, cf in enumerate(cf_order):
                subset = df[(df['size'] == size) & (df['cf_true'] == cf)]
                n = len(subset)
                if n > 0:
                    matrix[i, j] = subset['significant_raw'].sum() / n * 100

        im = ax.imshow(matrix, cmap='RdYlBu', vmin=0, vmax=100, aspect='auto')
        for i in range(len(sizes_present)):
            for j in range(len(cf_order)):
                val = matrix[i, j]
                color = 'white' if val < 30 or val > 85 else 'black'
                ax.text(j, i, f'{val:.0f}', ha='center', va='center',
                        fontsize=7, fontweight='bold', color=color)

        ax.set_xticks(range(len(cf_labels)))
        ax.set_xticklabels(cf_labels, fontsize=8, rotation=45, ha='right')
        ax.set_yticks(range(len(sizes_present)))
        ax.set_yticklabels(sizes_present if ax == axes[0] else [], fontsize=10)
        ax.set_xlabel('Target CF', fontsize=10)
        ax.set_title(f'index={ref_cf*100:.0f}% CF', fontsize=12, fontweight='bold')

        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)

    fig.suptitle('Longitudinal Detection by Size: Effect of Index Sample CF',
                 fontsize=15, fontweight='bold', y=1.02)
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.012, 0.7])
    fig.colorbar(im, cax=cbar_ax).set_label('Sensitivity (%)', fontsize=10)

    plt.savefig(os.path.join(OUTPUT_DIR, 'longitudinal_ref_comparison_heatmaps.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅ longitudinal_ref_comparison_heatmaps.pdf")

# ═══════════════════════════════════════════════
# FIGURE: False Positive Rate (within-family design)
# ═══════════════════════════════════════════════

def load_fpr_data(fpr_csv, ref_cf=None):
    """Load FPR results, filter to tested rows and optional ref CF."""
    df = pd.read_csv(fpr_csv)
    df = df[df['status'] == 'tested'].copy()
    if ref_cf is not None and 'ref_cf' in df.columns:
        df = df[np.isclose(df['ref_cf'], ref_cf)]
    return df

def plot_fpr_summary(fpr_df, ref_pct):
    """
    4-panel single-row FPR summary:
      a. P-value distribution (uniform under null)
      b. FPR by index sample CF
      c. FPR by event type
      d. FPR by event size (≥5 Mb only)
    """
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fig.suptitle('False Positive Rate (Within-Family Null, phased caller)',
                 fontsize=14, fontweight='bold', y=1.02)
    alpha = 0.05

    # ── Panel a: P-value distribution ──
    ax = axes[0]
    pvals = fpr_df['p_onesided'].dropna()
    if len(pvals) > 0:
        ax.hist(pvals, bins=20, range=(0, 1), color=blue3, alpha=0.7,
                edgecolor='white', density=True, label=f'N = {len(pvals)}')
        ax.axhline(1.0, color='grey', linestyle='--', alpha=0.7,
                   label='Uniform (expected)')
        ax.axvline(alpha, color=red4, linestyle='-', linewidth=1.5, alpha=0.8,
                   label=f'\u03b1 = {alpha}')
    ax.set_xlabel('p-value (one-sided)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('P-value Distribution (null regions)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.2)

    # ── Panel b: FPR by index sample CF ──
    ax = axes[1]
    bar_vals = []
    if 'ref_cf' in fpr_df.columns:
        ref_cfs = sorted(fpr_df['ref_cf'].unique())
        bar_labels = []
        for rcf in ref_cfs:
            rcf_data = fpr_df[np.isclose(fpr_df['ref_cf'], rcf)]
            n = len(rcf_data)
            n_fp = int(rcf_data['significant'].sum())
            bar_vals.append(n_fp / n * 100 if n > 0 else 0)
            bar_labels.append(f'{rcf*100:.0f}%')
        x_pos = range(len(ref_cfs))
        ax.bar(x_pos, bar_vals, color=blue4, edgecolor='white')
        ax.axhline(alpha * 100, color='grey', linestyle='--', alpha=0.7,
                   label=f'\u03b1 = {alpha*100:.0f}%')
        for i, (rcf, val) in enumerate(zip(ref_cfs, bar_vals)):
            rcf_data = fpr_df[np.isclose(fpr_df['ref_cf'], rcf)]
            n_fp = int(rcf_data['significant'].sum())
            ax.text(i, val + 0.3, f'{n_fp}/{len(rcf_data)}\n({val:.1f}%)',
                    ha='center', fontsize=8, fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(bar_labels, fontsize=10)
        ax.set_xlabel('Index Sample CF', fontsize=11)
        ax.set_ylim(0, max(max(bar_vals) * 2, alpha * 100 * 3) if bar_vals else 15)
    else:
        ax.text(0.5, 0.5, 'ref_cf not available', ha='center', va='center',
                transform=ax.transAxes, fontsize=12, color='grey')
        ax.set_ylim(0, 15)
    ax.set_ylabel('FPR (%)', fontsize=11)
    ax.set_title('FPR by Index Sample CF', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.2)

    # ── Panel c: FPR by event type ──
    ax = axes[2]
    events = ['CN-LOH', 'GAIN', 'LOSS']
    bar_vals = []
    for evt in events:
        ed = fpr_df[fpr_df['event'] == evt]
        bar_vals.append(ed['significant'].sum() / len(ed) * 100 if len(ed) > 0 else 0)
    bar_colors = [EVENT_COLORS.get(e, grey4) for e in events]
    ax.bar(range(len(events)), bar_vals, color=bar_colors, edgecolor='white')
    ax.axhline(alpha * 100, color='grey', linestyle='--', alpha=0.7,
               label=f'\u03b1 = {alpha*100:.0f}%')
    for i, (evt, val) in enumerate(zip(events, bar_vals)):
        ed = fpr_df[fpr_df['event'] == evt]
        n_fp = int(ed['significant'].sum())
        ax.text(i, val + 0.3, f'{n_fp}/{len(ed)}\n({val:.1f}%)',
                ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(range(len(events)))
    ax.set_xticklabels(events, fontsize=11)
    ax.set_ylabel('FPR (%)', fontsize=11)
    ax.set_title('FPR by Event Type', fontsize=12, fontweight='bold')
    ax.set_ylim(0, max(max(bar_vals) * 2, alpha * 100 * 3) if bar_vals else 15)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.2)

    # ── Panel d: FPR by event size (≥5 Mb) ──
    ax = axes[3]
    sizes_present = [s for s in ALL_SIZES if s in fpr_df['size'].unique() and s not in ('2Mb', '3Mb')]
    bar_vals = []
    for sz in sizes_present:
        sd = fpr_df[fpr_df['size'] == sz]
        bar_vals.append(sd['significant'].sum() / len(sd) * 100 if len(sd) > 0 else 0)
    sz_colors = [SIZE_COLORS.get(s, grey4) for s in sizes_present]
    ax.bar(range(len(sizes_present)), bar_vals, color=sz_colors, edgecolor='white')
    ax.axhline(alpha * 100, color='grey', linestyle='--', alpha=0.7,
               label=f'\u03b1 = {alpha*100:.0f}%')
    for i, (sz, val) in enumerate(zip(sizes_present, bar_vals)):
        sd = fpr_df[fpr_df['size'] == sz]
        n_fp = int(sd['significant'].sum())
        ax.text(i, val + 0.3, f'{n_fp}/{len(sd)}',
                ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(range(len(sizes_present)))
    ax.set_xticklabels(sizes_present, fontsize=10)
    ax.set_ylabel('FPR (%)', fontsize=11)
    ax.set_title('FPR by Event Size', fontsize=12, fontweight='bold')
    ax.set_ylim(0, max(max(bar_vals) * 2, alpha * 100 * 3) if bar_vals else 15)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.2)

    for a in axes:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        for sp in ['bottom', 'left']:
            a.spines[sp].set_linewidth(1.5)
            a.spines[sp].set_color(grey3)
        a.xaxis.set_tick_params(width=1, color=grey3, length=6)
        a.yaxis.set_tick_params(width=1, color=grey3, length=6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fpr_analysis_summary.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  \u2705 fpr_analysis_summary.pdf")

def plot_sensitivity_vs_fpr(long_df, fpr_df, ref_pct):
    """
    3-panel: sensitivity curves with FPR shown as horizontal band, by event type.

    Within-family FPR has no target_cf dimension (null is independent of
    spike-in CF), so FPR is a single rate per event type shown as a
    horizontal line with 95% CI band.
    """
    fig, axes = plt.subplots(1, 3, figsize=(17, 6), sharey=True)
    alpha = 0.05

    for ax, evt in zip(axes, ['CN-LOH', 'GAIN', 'LOSS']):
        color = EVENT_COLORS[evt]

        # Sensitivity curve
        evt_tp = long_df[long_df['event'] == evt]
        if len(evt_tp) > 0:
            cf_vals = sorted(evt_tp['cf_true'].unique())
            sens = [evt_tp[evt_tp['cf_true'] == c]['significant_raw'].mean() * 100
                    for c in cf_vals]
            ax.plot([c * 100 for c in cf_vals], sens,
                    's-', color=color, linewidth=2.5, markersize=7,
                    label='Sensitivity (TPR)')

        # FPR as horizontal line — filter to matching index CF
        ref_cf_val = float(ref_pct) / 100.0
        fpr_df_filt = fpr_df[np.isclose(fpr_df['ref_cf'], ref_cf_val)] if 'ref_cf' in fpr_df.columns else fpr_df
        evt_fp = fpr_df_filt[fpr_df_filt['event'] == evt]
        if len(evt_fp) > 0:
            fpr_val = evt_fp['significant'].mean() * 100
            n_fp = int(evt_fp['significant'].sum())
            n_total = len(evt_fp)
            ax.axhline(fpr_val, color=red4, linestyle='-', linewidth=2, alpha=0.7,
                       label=f'FPR: {fpr_val:.1f}% ({n_fp}/{n_total})')

            # 95% CI band (Wilson interval)
            try:
                from statsmodels.stats.proportion import proportion_confint
                ci_lo, ci_hi = proportion_confint(n_fp, n_total, alpha=0.05,
                                                  method='wilson')
                ax.axhspan(ci_lo * 100, ci_hi * 100, color=red4, alpha=0.08)
            except Exception:
                pass

        ax.axhline(alpha * 100, color='grey', linestyle=':', alpha=0.5,
                   label=f'\u03b1 = {alpha*100:.0f}%')
        ax.set_title(evt, fontsize=15, fontweight='bold')
        ax.set_xlabel('Cell Fraction (%)', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Rate (%)', fontsize=12)
        ax.set_xscale('log')
        ax.set_xticks([0.1, 0.5, 1, 5, 10, 50, 100])
        ax.get_xaxis().set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f'{x:g}'))
        ax.set_xlim(0.08, 130)
        ax.set_ylim(-2, 105)
        ax.legend(fontsize=9, loc='center right')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.grid(axis='y', alpha=0.2)

    fig.suptitle(f'Sensitivity vs False Positive Rate (index={ref_pct}% CF)',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'sensitivity_vs_fpr_by_event_index{ref_pct}pct.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  \u2705 sensitivity_vs_fpr_by_event_index{ref_pct}pct.pdf")

def plot_roc_curves(long_df, fpr_df, ref_pct):
    """
    ROC curves: vary the p-value threshold and plot TPR vs FPR.

    Since FPR comes from within-family null tests (single pool of p-values)
    and TPR comes from longitudinal results (per CF), we get one ROC curve
    per CF level, per event type.

    3 panels (CN-LOH, GAIN, LOSS), each with ROC curves at selected CFs.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # P-value thresholds to sweep
    thresholds = np.concatenate([
        np.arange(0.0001, 0.001, 0.0001),
        np.arange(0.001, 0.01, 0.001),
        np.arange(0.01, 0.1, 0.01),
        np.arange(0.1, 1.01, 0.05),
    ])

    # Selected CFs to show (as fractions)
    show_cfs = [0.005, 0.01, 0.025, 0.05, 0.10]
    cf_colors = {0.005: '#d73027', 0.01: '#fc8d59', 0.025: '#fee090',
                 0.05: '#91bfdb', 0.10: '#4575b4'}
    cf_styles = {0.005: '-', 0.01: '-', 0.025: '-', 0.05: '-', 0.10: '-'}

    for ax, evt in zip(axes, ['CN-LOH', 'GAIN', 'LOSS']):

        # FPR p-values (null regions for this event type)
        null_pvals = fpr_df[fpr_df['event'] == evt]['p_onesided'].dropna().values

        if len(null_pvals) == 0:
            ax.text(0.5, 0.5, 'No null data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12)
            continue

        # Get available CFs for this event
        evt_tp = long_df[long_df['event'] == evt]
        available_cfs = sorted(evt_tp['cf_true'].unique())

        for cf in show_cfs:
            # Find closest available CF
            closest = min(available_cfs, key=lambda x: abs(x - cf))
            if abs(closest - cf) > 0.005:
                continue

            tp_pvals = evt_tp[evt_tp['cf_true'] == closest]['p_onesided'].dropna().values
            if len(tp_pvals) == 0:
                continue

            # Compute TPR and FPR at each threshold
            tpr_vals = []
            fpr_vals = []
            for t in thresholds:
                tpr_vals.append(np.mean(tp_pvals < t) * 100)
                fpr_vals.append(np.mean(null_pvals < t) * 100)

            ax.plot(fpr_vals, tpr_vals,
                    linestyle=cf_styles.get(cf, '-'),
                    color=cf_colors.get(cf, 'grey'),
                    linewidth=2, label=f'CF = {cf*100:g}%')

            # Mark the operating point at p=0.05
            tpr_at_05 = np.mean(tp_pvals < 0.05) * 100
            fpr_at_05 = np.mean(null_pvals < 0.05) * 100
            ax.plot(fpr_at_05, tpr_at_05, 'o', color=cf_colors.get(cf, 'grey'),
                    markersize=8, markeredgecolor='black', markeredgewidth=1,
                    zorder=5)

        # Diagonal (chance)
        ax.plot([0, 100], [0, 100], 'k--', linewidth=1, alpha=0.3, label='Chance')

        # Mark α = 5% line
        ax.axvline(np.mean(null_pvals < 0.05) * 100, color='grey',
                   linestyle=':', alpha=0.5)

        ax.set_title(evt, fontsize=15, fontweight='bold')
        ax.set_xlabel('False Positive Rate (%)', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('True Positive Rate / Sensitivity (%)', fontsize=12)
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_aspect('equal')
        ax.legend(fontsize=9, loc='lower right', title='Operating points at p=0.05',
                  title_fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.grid(alpha=0.15)

    fig.suptitle(f'ROC Curves by Event Type and Cell Fraction (index={ref_pct}% CF)',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'roc_curves_by_event_and_cf_index{ref_pct}pct.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  \u2705 roc_curves_by_event_and_cf_index{ref_pct}pct.pdf")


# ═══════════════════════════════════════════════
# FIGURE: Three-pass filter usage diagnostic
# ═══════════════════════════════════════════════

def plot_pass_usage(df, title_suffix, output_name):
    """
    2-panel bar chart showing how often each pass (1/2/3) of the three-pass
    het filter is triggered, broken down by event type and event size.
    Useful for documenting that pass 3 (near-hom extension) is only triggered
    for high-CF or small events where pass 1/2 find insufficient het SNPs.
    """
    if 'pass_used' not in df.columns:
        print(f"  ⚠ pass_used column not found, skipping {output_name}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(f'Three-Pass Het Filter Usage ({title_suffix})',
                 fontsize=14, fontweight='bold', y=1.01)

    pass_colors = {1: blue4, 2: orange3, 3: red4}
    pass_labels = {1: 'Pass 1 (0.2–0.8)', 2: 'Pass 2 (0.05–0.95)', 3: 'Pass 3 (0.001–0.999)'}

    for ax, (group_col, group_order) in zip(axes, [
        ('event', ['CN-LOH', 'GAIN', 'LOSS']),
        ('size',  [s for s in ALL_SIZES if s in df['size'].unique()])
    ]):
        groups = [g for g in group_order if g in df[group_col].unique()]
        x = np.arange(len(groups))
        width = 0.25

        for i, pass_n in enumerate([1, 2, 3]):
            counts = []
            for g in groups:
                grp_df = df[df[group_col] == g]
                n = len(grp_df)
                pct = (grp_df['pass_used'] == pass_n).sum() / n * 100 if n > 0 else 0
                counts.append(pct)
            bars = ax.bar(x + i * width, counts, width,
                          label=pass_labels[pass_n],
                          color=pass_colors[pass_n], alpha=0.85, edgecolor='white')

        ax.set_xticks(x + width)
        ax.set_xticklabels(groups, fontsize=11)
        ax.set_ylabel('% of calls using this pass', fontsize=11)
        ax.set_ylim(0, 108)
        ax.legend(fontsize=9, loc='upper right')
        ax.set_title(f'By {group_col.title()}', fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.grid(axis='y', alpha=0.2)

        # Annotate total N per group
        for i, g in enumerate(groups):
            n = len(df[df[group_col] == g])
            ax.text(x[i] + width, -5, f'N={n}', ha='center', fontsize=8, color=grey4)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, output_name + '.pdf'),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  ✅ {output_name}.pdf")


# ═══════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Phased mCA Detection: Sensitivity Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # Both + three-way comparison (recommended)
  python %(prog)s --phased phased_per_sample_results.csv --longitudinal longitudinal_phased_results.csv --unphased truth_comparison.csv

  # Ground-truth only
  python %(prog)s --phased phased_per_sample_results.csv --unphased truth_comparison.csv

  # Longitudinal only
  python %(prog)s --longitudinal longitudinal_phased_results.csv --unphased truth_comparison.csv

  # Longitudinal with specific reference CF for single-ref figures
  python %(prog)s --longitudinal results.csv --ref-cf 0.75

  # Custom output directory
  python %(prog)s --phased results.csv -o my_figures/

  # Include FPR analysis
  python %(prog)s --longitudinal longitudinal_phased_results.csv --fpr longitudinal_fpr_null_region_results.csv --unphased truth_comparison.csv
""")
    parser.add_argument('--phased', default=None,
                        help='Per-sample ground-truth phased results CSV (default: skip)')
    parser.add_argument('--longitudinal', default=None,
                        help='Longitudinal phased results CSV (default: skip)')
    parser.add_argument('--unphased', default=None,
                        help='Unphased truth_comparison CSV for comparison figures (default: skip)')
    parser.add_argument('-o', '--output-dir', default=None,
                        help='Output directory for figures')
    parser.add_argument('--ref-cf', type=float, default=1.0,
                        help='Reference CF for single-ref longitudinal figures (default: 1.0)')
    parser.add_argument('--ref-cfs', type=float, nargs='+', default=[1.0, 0.75, 0.5, 0.25],
                        help='Reference CFs for multi-ref comparison (default: 1.0 0.75 0.5 0.25)')
    parser.add_argument('--fpr', default=None,
                        help='Null-region FPR results CSV (from notebook cell)')
    parser.add_argument('--use-defaults', action='store_true',
                        help='Use built-in default paths (for running without arguments)')

    args = parser.parse_args()

    # ── Resolve paths ──
    if args.use_defaults:
        PHASED_CSV = args.phased or DEFAULT_PHASED_CSV
        LONGITUDINAL_CSV = args.longitudinal or DEFAULT_LONGITUDINAL_CSV
        UNPHASED_CSV = args.unphased or DEFAULT_UNPHASED_CSV
        OUTPUT_DIR = args.output_dir or DEFAULT_OUTPUT_DIR
    else:
        PHASED_CSV = args.phased
        LONGITUDINAL_CSV = args.longitudinal
        UNPHASED_CSV = args.unphased
        OUTPUT_DIR = args.output_dir or DEFAULT_OUTPUT_DIR

    FPR_CSV = args.fpr
    LONGITUDINAL_REF_CF = args.ref_cf
    LONGITUDINAL_REF_CFS = args.ref_cfs

    if not PHASED_CSV and not LONGITUDINAL_CSV and not FPR_CSV:
        print("Error: provide at least one of --phased, --longitudinal, or --fpr (or use --use-defaults)")
        parser.print_help()
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")

    # ── Load FPR data early so it can be passed to sensitivity figures ──
    fpr_df_global = None
    if FPR_CSV and os.path.exists(FPR_CSV):
        try:
            fpr_df_global = load_fpr_data(FPR_CSV)
            print(f"Pre-loaded FPR data: {len(fpr_df_global)} null tests")
        except Exception as e:
            print(f"Could not pre-load FPR data: {e}")

    # # ══════════════════════════════════════════════
    # # Per-sample ground-truth phased results
    # # ══════════════════════════════════════════════
    # if PHASED_CSV and os.path.exists(PHASED_CSV):
    #     df = load_and_normalise(PHASED_CSV, longitudinal=False)
    #     print(f"\nLoaded {len(df)} per-sample phased results from {PHASED_CSV}")
    #     print(f"  Events:      {df['event'].value_counts().to_dict()}")
    #     print(f"  Sizes:       {df['size'].value_counts().to_dict()}")
    #     print(f"  CF range:    {df['cf_true'].min():.4f} – {df['cf_true'].max():.4f}")
    #     print(f"  Chromosomes: {df['chromosome'].nunique()}")

    #     # All sizes
    #     print("\nAll sizes:")
    #     plot_sensitivity_6panel(df, 'All Sizes', 'sensitivity_analysis_phased_all_sizes')
    #     plot_chromosome_analysis(df, 'All Sizes', 'chromosome_analysis_phased_all_sizes')
    #     plot_sensitivity_heatmaps(df, 'All Sizes', 'sensitivity_heatmaps_phased_all_sizes')

    #     # 5Mb+
    #     df_5plus = df[df['size'].isin(SIZES_5PLUS)]
    #     print(f"\n≥5Mb subset: {len(df_5plus)} results")
    #     plot_sensitivity_6panel(df_5plus, '≥5 Mb', 'sensitivity_analysis_phased_5Mb_plus')
    #     plot_chromosome_analysis(df_5plus, '≥5 Mb', 'chromosome_analysis_phased_5Mb_plus')
    #     plot_sensitivity_heatmaps(df_5plus, '≥5 Mb', 'sensitivity_heatmaps_phased_5Mb_plus')

    #     # Comparison with unphased
    #     if UNPHASED_CSV and os.path.exists(UNPHASED_CSV):
    #         print("\nComparison figures (unphased vs phased ground-truth):")
    #         unphased_by_event, unphased_by_size = load_unphased_sensitivity(UNPHASED_CSV)
    #         plot_comparison_curves(df, unphased_by_event)
    #         plot_comparison_heatmaps(df, unphased_by_size)
    # elif PHASED_CSV:
    #     print(f"Skipping per-sample phased analysis (not found: {PHASED_CSV})")

    # ══════════════════════════════════════════════
    # Longitudinal phased results (single reference CF)
    # ══════════════════════════════════════════════
    if LONGITUDINAL_CSV and os.path.exists(LONGITUDINAL_CSV):
        ref_pct = f'{LONGITUDINAL_REF_CF*100:.0f}'
        df_long = load_and_normalise(LONGITUDINAL_CSV, longitudinal=True, ref_cf=LONGITUDINAL_REF_CF)
        print(f"\nLoaded {len(df_long)} longitudinal results (ref={ref_pct}% CF)")
        print(f"  Events:      {df_long['event'].value_counts().to_dict()}")
        print(f"  Sizes:       {df_long['size'].value_counts().to_dict()}")
        print(f"  CF range:    {df_long['cf_true'].min():.4f} – {df_long['cf_true'].max():.4f}")
        print(f"  Chromosomes: {df_long['chromosome'].nunique()}")

        suffix = f'ref{ref_pct}pct'

        # All sizes
        print(f"\nLongitudinal (ref={ref_pct}%), all sizes:")
        plot_sensitivity_6panel(df_long, f'Longitudinal index={ref_pct}%, All Sizes',
                                f'sensitivity_analysis_longitudinal_{suffix}_all_sizes',
                                fpr_df=fpr_df_global, ref_cf=LONGITUDINAL_REF_CF)
        plot_chromosome_analysis(df_long, f'Longitudinal index={ref_pct}%, All Sizes',
                                 f'chromosome_analysis_longitudinal_{suffix}_all_sizes')
        plot_sensitivity_heatmaps(df_long, f'Longitudinal index={ref_pct}%, All Sizes',
                                  f'sensitivity_heatmaps_longitudinal_{suffix}_all_sizes')
        plot_pass_usage(df_long, f'Longitudinal index={ref_pct}%, All Sizes',
                        f'pass_usage_diagnostic_{suffix}_all_sizes')

        # 5Mb+
        df_long_5plus = df_long[df_long['size'].isin(SIZES_5PLUS)]
        print(f"\nLongitudinal ≥5Mb subset: {len(df_long_5plus)} results")
        plot_sensitivity_6panel(df_long_5plus, f'Longitudinal index={ref_pct}%, ≥5 Mb',
                                f'sensitivity_analysis_longitudinal_{suffix}_5Mb_plus',
                                fpr_df=fpr_df_global, ref_cf=LONGITUDINAL_REF_CF)
        plot_chromosome_analysis(df_long_5plus, f'Longitudinal index={ref_pct}%, ≥5 Mb',
                                 f'chromosome_analysis_longitudinal_{suffix}_5Mb_plus')
        plot_sensitivity_heatmaps(df_long_5plus, f'Longitudinal index={ref_pct}%, ≥5 Mb',
                                  f'sensitivity_heatmaps_longitudinal_{suffix}_5Mb_plus')

        # Comparison with unphased
        if UNPHASED_CSV and os.path.exists(UNPHASED_CSV):
            print(f"\nComparison figures (unphased vs longitudinal ref={ref_pct}%):")
            unphased_by_event, unphased_by_size = load_unphased_sensitivity(UNPHASED_CSV)
            plot_comparison_curves(df_long, unphased_by_event,
                                  phased_label=f'Longitudinal (index={ref_pct}%)')
            # Rename to avoid overwriting the per-sample comparison
            os.rename(os.path.join(OUTPUT_DIR, 'sensitivity_comparison_unphased_vs_phased.pdf'),
                      os.path.join(OUTPUT_DIR, f'sensitivity_comparison_unphased_vs_longitudinal_{suffix}.pdf'))
            print(f"  → renamed to sensitivity_comparison_unphased_vs_longitudinal_{suffix}.pdf")

            plot_comparison_heatmaps(df_long, unphased_by_size,
                                    phased_title=f'Longitudinal (index={ref_pct}%)')
            os.rename(os.path.join(OUTPUT_DIR, 'heatmap_comparison_unphased_vs_phased.pdf'),
                      os.path.join(OUTPUT_DIR, f'heatmap_comparison_unphased_vs_longitudinal_{suffix}.pdf'))
            print(f"  → renamed to heatmap_comparison_unphased_vs_longitudinal_{suffix}.pdf")

        # # ── Three-way comparison (requires all three datasets) ──
        # if PHASED_CSV and os.path.exists(PHASED_CSV) and UNPHASED_CSV and os.path.exists(UNPHASED_CSV):
        #     print(f"\nThree-way comparison (unphased vs ground-truth vs longitudinal ref={ref_pct}%):")
        #     # Reload ground-truth if not already loaded
        #     try:
        #         gt_df = load_and_normalise(PHASED_CSV, longitudinal=False)
        #         unphased_by_event, unphased_by_size = load_unphased_sensitivity(UNPHASED_CSV)
        #         plot_threeway_overall(gt_df, df_long, unphased_by_event, ref_pct)
        #         plot_threeway_curves(gt_df, df_long, unphased_by_event, ref_pct)
        #         plot_threeway_heatmaps(gt_df, df_long, unphased_by_size, ref_pct)
        #     except Exception as e:
        #         print(f"  Three-way comparison failed: {e}")

        # ── Multi-reference CF comparison ──
        print(f"\n{'='*60}")
        print("Multi-reference CF comparison figures")
        print(f"{'='*60}")
        ref_dfs = {}
        for rcf in LONGITUDINAL_REF_CFS:
            try:
                rdf = load_and_normalise(LONGITUDINAL_CSV, longitudinal=True, ref_cf=rcf)
                if len(rdf) > 0:
                    ref_dfs[rcf] = rdf
                    print(f"  ref={rcf*100:.0f}%: {len(rdf)} results, "
                          f"{rdf['significant_raw'].sum():.0f} significant")
            except Exception as e:
                print(f"  ref={rcf*100:.0f}%: skipped ({e})")

        if len(ref_dfs) > 1:
            plot_longitudinal_ref_comparison_overall(ref_dfs)
            plot_longitudinal_ref_comparison_by_event(ref_dfs)
            plot_longitudinal_ref_comparison_by_size(ref_dfs)
            plot_longitudinal_ref_heatmap_grid(ref_dfs)
        else:
            print("  Need ≥2 reference CFs for comparison figures, skipping.")

    elif LONGITUDINAL_CSV:
        print(f"\nSkipping longitudinal analysis (not found: {LONGITUDINAL_CSV})")

    # ══════════════════════════════════════════════
    # FPR Analysis (null-region results)
    # ══════════════════════════════════════════════
    if FPR_CSV and os.path.exists(FPR_CSV):
        print(f"\n{'='*60}")
        print("False Positive Rate Analysis (within-family null)")
        print(f"{'='*60}")
        fpr_df = load_fpr_data(FPR_CSV)
        n_fp = int(fpr_df['significant'].sum())
        n_total = len(fpr_df)
        print(f"  Tested: {n_total}, False positives: {n_fp} ({n_fp/n_total*100:.2f}%)")

        print("\nFPR figures:")
        ref_pct = f'{LONGITUDINAL_REF_CF*100:.0f}'
        plot_fpr_summary(fpr_df, ref_pct)

        # Sensitivity vs FPR overlay (needs longitudinal data)
        if LONGITUDINAL_CSV and os.path.exists(LONGITUDINAL_CSV):
            try:
                long_for_fpr = load_and_normalise(
                    LONGITUDINAL_CSV, longitudinal=True, ref_cf=LONGITUDINAL_REF_CF
                )
                plot_sensitivity_vs_fpr(long_for_fpr, fpr_df, ref_pct)
                plot_roc_curves(long_for_fpr, fpr_df, ref_pct)
            except Exception as e:
                print(f"  Sensitivity vs FPR overlay failed: {e}")

    elif FPR_CSV:
        print(f"\nSkipping FPR analysis (not found: {FPR_CSV})")

    print(f"\nAll figures saved to {OUTPUT_DIR}/")
