#!/usr/bin/env python3
"""
Analyze mCA Caller Performance on Simulated Data

Calculates sensitivity, specificity, precision metrics and generates
publication-quality figures.

Includes detection breakdown analysis showing what happened to 'missed' events:
- Type misclassification vs boundary mismatch vs truly undetected
- Misclassification matrix (truth type vs called type)
- Effective sensitivity (strict TP vs TP + partial matches)
"""

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


import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from matplotlib.patches import Patch
import os

# Editable text in Illustrator / Inkscape
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

# Configuration - set in main()
RESULTS_DIR = None
TRUTH_FILE = None
OUTPUT_DIR = None

# Set plot style
# sns.set_style("whitegrid")
# sns.set_context("paper", font_scale=1.2)

# ============================================================
# OUTCOME CLASSIFICATION (for detection breakdown)
# ============================================================

def classify_outcome(row):
    """Classify each truth comparison row into a detection outcome category."""
    status = str(row['Status'])
    if status == 'TRUE_POSITIVE':
        return 'TP'
    elif status.startswith('PARTIAL_MATCH_TYPE'):
        return 'Type mismatch'
    elif status.startswith('PARTIAL_MATCH_OVER'):
        return 'Boundary mismatch'
    elif status.startswith('PARTIAL_MATCH_UNDER'):
        return 'Boundary mismatch'
    elif status.startswith('PARTIAL_MATCH_SIZE'):
        return 'Boundary mismatch'
    elif status.startswith('PARTIAL_MATCH'):
        return 'Boundary mismatch'
    else:
        return 'Not detected'

# Colours for detection breakdown
OUTCOME_COLORS = {
    'TP': '#2ca02c',
    'Type mismatch': '#ff7f0e',
    'Boundary mismatch': '#d4a017',
    'Not detected': '#d62728',
}
OUTCOME_ORDER = ['TP', 'Type mismatch', 'Boundary mismatch', 'Not detected']

# ============================================================
# ORIGINAL METRIC HELPERS
# ============================================================

def calculate_metrics(df):
    """Calculate sensitivity, precision, and other metrics"""
    total = len(df)
    detected = df['Detected'].sum()
    missed = total - detected
    fp_count = df['False_Positives'].sum() if 'False_Positives' in df.columns else 0
    
    sensitivity = detected / total if total > 0 else 0
    precision = detected / (detected + fp_count) if (detected + fp_count) > 0 else 0
    fdr = fp_count / (detected + fp_count) if (detected + fp_count) > 0 else 0
    fp_per_sample = fp_count / total if total > 0 else 0
    
    return {
        'Total': total,
        'True_Positives': detected,
        'False_Negatives': missed,
        'False_Positives': fp_count,
        'Sensitivity': sensitivity,
        'Precision': precision,
        'FP_per_sample': fp_per_sample,
        'FDR': fdr
    }

def metrics_by_group(df, group_col, group_label):
    """Calculate metrics for each group"""
    results = []
    for group in sorted(df[group_col].unique()):
        group_df = df[df[group_col] == group]
        metrics = calculate_metrics(group_df)
        metrics[group_label] = group
        results.append(metrics)
    return pd.DataFrame(results)

# ============================================================
# ORIGINAL FIGURES
# ============================================================

def plot_chromosome_analysis(df):
    """Sensitivity breakdown by chromosome"""
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), gridspec_kw={'height_ratios': [1, 1]})
    
    chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    
    chrom_data = []
    for chrom in chrom_order:
        chrom_df = df[df['Truth_Chromosome'] == chrom]
        if len(chrom_df) == 0:
            continue
        detected = chrom_df['Detected'].sum()
        total = len(chrom_df)
        sensitivity = detected / total * 100
        fp = chrom_df['False_Positives'].sum()
        chrom_data.append({
            'Chromosome': chrom.replace('chr', ''),
            'Sensitivity': sensitivity,
            'Total': total,
            'TP': detected,
            'FP': fp
        })
    
    chrom_df_plot = pd.DataFrame(chrom_data)
    
    bars = axes[0].bar(chrom_df_plot['Chromosome'], chrom_df_plot['Sensitivity'], 
                       color='#2E86AB', edgecolor='white', linewidth=0.5)
    axes[0].axhline(y=chrom_df_plot['Sensitivity'].mean(), color='#E63946', 
                    linestyle='--', linewidth=1.5, label=f"Mean: {chrom_df_plot['Sensitivity'].mean():.1f}%")
    
    for i, row in chrom_df_plot.iterrows():
            axes[0].text(i, row['Sensitivity'] + 3, 
                        f"{int(row['TP'])}/{int(row['Total'])}", 
                        ha='center', fontsize=6, rotation=0)
    
    axes[0].set_ylabel('Sensitivity (%)', fontweight='bold')
    axes[0].set_title('Sensitivity by Chromosome (All Cell Fractions)', fontweight='bold', fontsize=14)
    axes[0].set_ylim(0, 105)
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)
    
    detectable = df[df['Truth_Fraction'] >= 0.10]
    
    pivot = detectable.pivot_table(
        index='Truth_Fraction',
        columns='Truth_Chromosome',
        values='Detected',
        aggfunc='mean'
    ) * 100

    pivot.index = [f"{x*100:.0f}%" for x in pivot.index]
    pivot = pivot.reindex(columns=[c for c in chrom_order if c in pivot.columns])
    pivot.columns = [c.replace('chr', '') for c in pivot.columns]
    pivot = pivot.iloc[::-1]  # highest CF at top

    sns.heatmap(pivot, annot=True, fmt='.0f', cmap='RdYlBu', vmin=0, vmax=100,
                cbar_kws={'label': 'Sensitivity (%)'}, ax=axes[1], linewidths=0.5)
    axes[1].set_title('Sensitivity by Chromosome and Cell Fraction (≥10%)',
                     fontweight='bold', fontsize=14)
    axes[1].set_xlabel('Chromosome', fontweight='bold')
    axes[1].set_ylabel('Cell Fraction', fontweight='bold')

    for ax in axes.flatten():
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'mCA_simulations_chromosome_analysis.pdf'),
                dpi=300, bbox_inches='tight')
    print(f"✅ Saved chromosome analysis")
    plt.close()

def plot_sensitivity_curves(df):
    """Plot sensitivity vs cell fraction by event type and size"""
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 3, hspace=0.3, wspace=0.3)
    
    # 1. Overall sensitivity by fraction
    ax1 = plt.subplot(gs[0, 0])
    frac_metrics = metrics_by_group(df, 'Truth_Fraction', 'Fraction')
    frac_metrics['Fraction_Pct'] = frac_metrics['Fraction'] * 100
    
    ax1.plot(frac_metrics['Fraction_Pct'], frac_metrics['Sensitivity'] * 100, 
             'o-', linewidth=2, markersize=8, color='#2E86AB')
    ax1.axhline(y=50, color='gray', linestyle='--', alpha=0.5, label='50% sensitivity')
    ax1.set_xlabel('Cell Fraction (%)', fontweight='bold')
    ax1.set_ylabel('Sensitivity (%)', fontweight='bold')
    ax1.set_title('Overall Sensitivity vs Cell Fraction', fontweight='bold', fontsize=14)
    ax1.set_ylim(-5, 105)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    for _, row in frac_metrics.iterrows():
        ax1.text(row['Fraction_Pct'], row['Sensitivity']*100 + 3, 
                f"{int(row['True_Positives'])}/{int(row['Total'])}", 
                ha='center', fontsize=9, fontweight='bold')
    
    # 2. Sensitivity by event type
    ax2 = plt.subplot(gs[0, 1])
    for event_type in sorted(df['Truth_Type'].unique()):
        type_df = df[df['Truth_Type'] == event_type]
        frac_metrics = metrics_by_group(type_df, 'Truth_Fraction', 'Fraction')
        frac_metrics['Fraction_Pct'] = frac_metrics['Fraction'] * 100
        ax2.plot(frac_metrics['Fraction_Pct'], frac_metrics['Sensitivity'] * 100,
                'o-', linewidth=2, markersize=6, label=event_type)
    ax2.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Cell Fraction (%)', fontweight='bold')
    ax2.set_ylabel('Sensitivity (%)', fontweight='bold')
    ax2.set_title('Sensitivity by Event Type', fontweight='bold', fontsize=14)
    ax2.set_ylim(-5, 105)
    ax2.legend(frameon=True, fancybox=True)
    ax2.grid(True, alpha=0.3)

    # 3. Sensitivity by size
    ax3 = plt.subplot(gs[0, 2])
    size_order = ['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    colors = ['#E63946', '#F77F00', '#FCBF49', '#06A77D', '#2E86AB', '#7B2D8E']
    for size, color in zip(size_order, colors):
        size_df = df[df['Truth_Length_Category'] == size]
        if len(size_df) > 0:
            frac_metrics = metrics_by_group(size_df, 'Truth_Fraction', 'Fraction')
            frac_metrics['Fraction_Pct'] = frac_metrics['Fraction'] * 100
            ax3.plot(frac_metrics['Fraction_Pct'], frac_metrics['Sensitivity'] * 100,
                    'o-', linewidth=2, markersize=6, label=size, color=color)
    ax3.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Cell Fraction (%)', fontweight='bold')
    ax3.set_ylabel('Sensitivity (%)', fontweight='bold')
    ax3.set_title('Sensitivity by Event Size', fontweight='bold', fontsize=14)
    ax3.set_ylim(-5, 105)
    ax3.legend(frameon=True, fancybox=True)
    ax3.grid(True, alpha=0.3)
    
    # 4. Precision (left axis) + FP/sample (right axis) vs cell fraction
    ax4 = plt.subplot(gs[1, 0])
    frac_metrics = metrics_by_group(df, 'Truth_Fraction', 'Fraction')
    frac_metrics['Fraction_Pct'] = frac_metrics['Fraction'] * 100

    prec_color = '#06A77D'
    fp_color   = '#E63946'

    # Left axis: precision
    ax4.plot(frac_metrics['Fraction_Pct'], frac_metrics['Precision'] * 100,
             's-', linewidth=2, markersize=8, color=prec_color, label='Precision (PPV)', zorder=3)
    ax4.axhline(y=50, color='gray', linestyle='--', alpha=0.4)
    ax4.set_xlabel('Cell Fraction (%)', fontweight='bold')
    ax4.set_ylabel('Precision (PPV) (%)', fontweight='bold', color=prec_color)
    ax4.tick_params(axis='y', labelcolor=prec_color)
    ax4.set_ylim(-5, 105)
    ax4.set_title('Precision & False Positive Rate', fontweight='bold', fontsize=14)
    ax4.grid(True, alpha=0.3)

    # Right axis: FP per sample
    ax4b = ax4.twinx()
    ax4b.plot(frac_metrics['Fraction_Pct'], frac_metrics['FP_per_sample'],
              '^--', linewidth=2, markersize=8, color=fp_color, label='FP / sample', zorder=3)
    # Set a generous fixed ceiling so the FP line never fills the axis
    fp_max = frac_metrics['FP_per_sample'].max()
    ax4b.set_ylim(0, max(fp_max * 6, 1.0))
    ax4b.set_ylabel('False Positives per Sample', fontweight='bold', color=fp_color)
    ax4b.tick_params(axis='y', labelcolor=fp_color)
    # Annotate peak
    peak_row = frac_metrics.loc[frac_metrics['FP_per_sample'].idxmax()]
    ax4b.annotate(f"peak {peak_row['FP_per_sample']:.2f}",
                  xy=(peak_row['Fraction_Pct'], peak_row['FP_per_sample']),
                  xytext=(peak_row['Fraction_Pct'] + 6, peak_row['FP_per_sample'] * 1.5),
                  fontsize=8, fontweight='bold', color=fp_color,
                  arrowprops=dict(arrowstyle='->', color=fp_color, lw=1.2))

    # Combined legend
    lines1, labels1 = ax4.get_legend_handles_labels()
    lines2, labels2 = ax4b.get_legend_handles_labels()
    ax4.legend(lines1 + lines2, labels1 + labels2, frameon=True, fontsize=9, loc='center right')

    # 5. Sensitivity by geometry
    ax5 = plt.subplot(gs[1, 1])
    for geom in sorted(df['Truth_Geometry'].unique()):
        geom_df = df[df['Truth_Geometry'] == geom]
        frac_metrics = metrics_by_group(geom_df, 'Truth_Fraction', 'Fraction')
        frac_metrics['Fraction_Pct'] = frac_metrics['Fraction'] * 100
        ax5.plot(frac_metrics['Fraction_Pct'], frac_metrics['Sensitivity'] * 100,
                'o-', linewidth=2, markersize=6, label=geom)
    ax5.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    ax5.set_xlabel('Cell Fraction (%)', fontweight='bold')
    ax5.set_ylabel('Sensitivity (%)', fontweight='bold')
    ax5.set_title('Sensitivity by Geometry', fontweight='bold', fontsize=14)
    ax5.set_ylim(-5, 105)
    ax5.legend(frameon=True, fancybox=True, fontsize=9)
    ax5.grid(True, alpha=0.3)

    # Panel 6 intentionally left blank (gs[1,2] unused)

    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)
    # Style the twin axis right spine to match
    ax4b.spines['top'].set_visible(False)
    ax4b.spines['right'].set_linewidth(1.5)
    ax4b.spines['right'].set_color(grey3)
    ax4b.yaxis.set_tick_params(width=1, color=grey3, length=6)
    
    plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_analysis.pdf'), 
                dpi=300, bbox_inches='tight')
    print(f"✅ Saved sensitivity curves")
    plt.close()

def plot_heatmaps(df):
    """Create heatmaps showing detection by fraction and size/type"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    pivot1 = df.pivot_table(
        index='Truth_Type',
        columns='Truth_Fraction',
        values='Detected',
        aggfunc='mean'
    ) * 100
    pivot1.columns = [f"{x*100:.1f}%" for x in pivot1.columns]
    sns.heatmap(pivot1, annot=True, fmt='.0f', cmap='RdYlBu', vmin=0, vmax=100,
                cbar_kws={'label': 'Sensitivity (%)'}, ax=axes[0], linewidths=0.5)
    axes[0].set_title('Sensitivity by Event Type and Cell Fraction', 
                     fontweight='bold', fontsize=14)
    axes[0].set_xlabel('Cell Fraction', fontweight='bold')
    axes[0].set_ylabel('Event Type', fontweight='bold')
    
    pivot2 = df.pivot_table(
        index='Truth_Length_Category',
        columns='Truth_Fraction',
        values='Detected',
        aggfunc='mean'
    ) * 100
    pivot2.columns = [f"{x*100:.1f}%" for x in pivot2.columns]
    size_order = ['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    pivot2 = pivot2.reindex([s for s in size_order if s in pivot2.index])
    sns.heatmap(pivot2, annot=True, fmt='.0f', cmap='RdYlBu', vmin=0, vmax=100,
                cbar_kws={'label': 'Sensitivity (%)'}, ax=axes[1], linewidths=0.5)
    axes[1].set_title('Sensitivity by Event Size and Cell Fraction',
                     fontweight='bold', fontsize=14)
    axes[1].set_xlabel('Cell Fraction', fontweight='bold')
    axes[1].set_ylabel('Event Size', fontweight='bold')

    for ax in axes.flatten():
        for axis in ['bottom','left', 'top', 'right']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'mCA_simulations_sensitivity_heatmaps.pdf'),
                dpi=300, bbox_inches='tight')
    print(f"✅ Saved sensitivity heatmaps")
    plt.close()

def plot_summary_table(df):
    """Create summary statistics table with size breakdown within cell fraction"""
    size_order = ['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    
    table_data = []
    table_data.append(['Cell Fraction', 'Size', 'Total', 'TP', 'FN', 'Sensitivity', 'FP/sample'])
    
    fraction_row_indices = [0]
    row_idx = 1
    
    for frac in sorted(df['Truth_Fraction'].unique()):
        frac_df = df[df['Truth_Fraction'] == frac]
        frac_label = f"{frac*100:.1f}%"
        metrics = calculate_metrics(frac_df)
        table_data.append([
            frac_label, 'ALL',
            str(metrics['Total']),
            str(metrics['True_Positives']),
            str(metrics['False_Negatives']),
            f"{metrics['Sensitivity']*100:.1f}%",
            f"{metrics['False_Positives']/metrics['Total']:.2f}"
        ])
        fraction_row_indices.append(row_idx)
        row_idx += 1
        
        for size in size_order:
            size_df = frac_df[frac_df['Truth_Length_Category'] == size]
            if len(size_df) == 0:
                continue
            m = calculate_metrics(size_df)
            table_data.append([
                '', size,
                str(m['Total']),
                str(m['True_Positives']),
                str(m['False_Negatives']),
                f"{m['Sensitivity']*100:.1f}%",
                ''
            ])
            row_idx += 1
    
    fig, ax = plt.subplots(figsize=(14, max(8, len(table_data) * 0.35)))
    ax.axis('off')
    
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.12, 0.10, 0.10, 0.10, 0.10, 0.14, 0.12])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    for j in range(7):
        table[(0, j)].set_facecolor('#2E86AB')
        table[(0, j)].set_text_props(weight='bold', color='white')
    
    for ri in fraction_row_indices[1:]:
        for j in range(7):
            table[(ri, j)].set_facecolor('#E8F4F8')
            table[(ri, j)].set_text_props(weight='bold')
    
    ax.set_title('mCA Detection Summary by Cell Fraction and Event Size', 
                 fontweight='bold', fontsize=14, pad=20)
    
    plt.savefig(os.path.join(OUTPUT_DIR, 'mCA_simulations_summary_table.pdf'),
                dpi=300, bbox_inches='tight')
    print(f"✅ Saved summary table")
    plt.close()

# ============================================================
# NEW: DETECTION BREAKDOWN FIGURES
# ============================================================

def plot_detection_breakdown(df):
    """
    Stacked bar chart: for each CF × event type, show fraction that was
    TP / type mismatch / boundary mismatch / truly missed.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    fig.suptitle('Detection Breakdown: What happened to every simulated event?', 
                 fontsize=14, fontweight='bold', y=1.02)

    event_types = sorted(df['Truth_Type'].unique())
    cf_values = sorted(df['CF_pct'].unique())
    cf_values = [c for c in cf_values if c >= 10]

    for ax_idx, etype in enumerate(event_types):
        ax = axes[ax_idx]
        subset = df[df['Truth_Type'] == etype]
        
        bottoms = np.zeros(len(cf_values))
        
        for outcome in OUTCOME_ORDER:
            counts = []
            raw_counts = []
            for cf in cf_values:
                cf_data = subset[subset['CF_pct'] == cf]
                total = len(cf_data)
                n = len(cf_data[cf_data['Outcome_simple'] == outcome])
                pct = (n / total * 100) if total > 0 else 0
                counts.append(pct)
                raw_counts.append(n)
            
            bars = ax.bar([str(int(c)) + '%' for c in cf_values], counts, bottom=bottoms,
                           color=OUTCOME_COLORS[outcome], label=outcome, 
                           edgecolor='white', linewidth=0.5)
            
            # Add count labels for non-trivial segments
            for i, (cnt, bot, raw) in enumerate(zip(counts, bottoms, raw_counts)):
                if cnt > 5:
                    text_color = 'white' if outcome != 'Boundary mismatch' else 'black'
                    ax.text(i, bot + cnt/2, str(raw), ha='center', va='center', 
                           fontsize=8, fontweight='bold', color=text_color)
            
            bottoms += counts
        
        ax.set_title(etype, fontsize=13, fontweight='bold')
        ax.set_xlabel('Cell Fraction')
        ax.set_ylim(0, 105)
        ax.set_ylabel('% of simulated events' if ax_idx == 0 else '')
        ax.tick_params(axis='x', rotation=45)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)

    handles = [Patch(facecolor=OUTCOME_COLORS[o], label=o) for o in OUTCOME_ORDER]
    fig.legend(handles=handles, loc='upper center', ncol=4, 
               bbox_to_anchor=(0.5, 0.98), fontsize=11)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'detection_breakdown_by_type.pdf'),
                dpi=300, bbox_inches='tight')
    print("✅ Saved detection breakdown by type")
    plt.close()


def plot_misclassification_matrix(df):
    """
    Heatmap: truth event type vs called event type for all type-mismatched events.
    """
    type_mismatches = df[df['Outcome'] == 'Type mismatch'].copy()

    if len(type_mismatches) == 0:
        print("ℹ️  No type mismatches found — skipping misclassification matrix")
        return

    type_mismatches['Truth_Type_clean'] = type_mismatches['Truth_Type'].str.replace('-', '').str.upper()
    type_mismatches['Detected_Type_clean'] = type_mismatches['Detected_Type'].astype(str).str.replace('-', '').str.upper()
    
    labels = ['CNLOH', 'GAIN', 'LOSS']
    matrix = pd.DataFrame(0, index=labels, columns=labels)
    for _, row in type_mismatches.iterrows():
        t = row['Truth_Type_clean']
        d = row['Detected_Type_clean']
        if t in labels and d in labels:
            matrix.loc[t, d] += 1
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), 
                              gridspec_kw={'width_ratios': [1, 1.3]})
    fig.suptitle('Event Type Misclassification Analysis', fontsize=14, fontweight='bold')
    
    # Left: confusion matrix
    ax = axes[0]
    im = ax.imshow(matrix.values, cmap='RdYlBu', aspect='auto')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel('Called as', fontsize=13)
    ax.set_ylabel('Truth', fontsize=13)
    ax.set_title('Misclassification Matrix', fontsize=12, fontweight='bold')
    
    for i in range(len(labels)):
        for j in range(len(labels)):
            val = matrix.values[i, j]
            if val > 0:
                ax.text(j, i, str(int(val)), ha='center', va='center', 
                        fontsize=16, fontweight='bold',
                        color='white' if val > matrix.values.max()/2 else 'black')
    plt.colorbar(im, ax=ax, label='Count', shrink=0.8)
    
    # Right: type mismatches by cell fraction
    ax2 = axes[1]
    cf_values = sorted(type_mismatches['CF_pct'].unique())
    
    truth_types_present = [t for t in ['CN-LOH', 'GAIN', 'LOSS'] if t in type_mismatches['Truth_Type'].values]
    type_colors = {'CN-LOH': '#2E86AB', 'GAIN': '#F77F00', 'LOSS': '#06A77D'}
    
    bar_width = 0.25
    x = np.arange(len(cf_values))
    
    for i, tt in enumerate(truth_types_present):
        counts = []
        for cf in cf_values:
            n = len(type_mismatches[(type_mismatches['CF_pct'] == cf) & 
                                    (type_mismatches['Truth_Type'] == tt)])
            counts.append(n)
        ax2.bar(x + i * bar_width, counts, bar_width, 
                label=f'{tt} misclassified', color=type_colors.get(tt, 'grey'), alpha=0.8)
    
    ax2.set_xticks(x + bar_width * (len(truth_types_present) - 1) / 2)
    ax2.set_xticklabels([f'{int(c)}%' for c in cf_values], rotation=45)
    ax2.set_xlabel('Cell Fraction')
    ax2.set_ylabel('Count')
    ax2.set_title('Type Mismatches by Cell Fraction', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', alpha=0.3)

    for ax in axes.flatten():
        for axis in ['bottom','left', 'top', 'right']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'misclassification_matrix.pdf'),
                dpi=300, bbox_inches='tight')
    print("✅ Saved misclassification matrix")
    plt.close()

    # Print text summary
    print("\n  Type mismatch breakdown:")
    for _, row in type_mismatches.iterrows():
        cf = row['CF_pct']
        truth = row['Truth_Type']
        detected = row['Detected_Type']
        size = row.get('Truth_Length_Mb', '?')
        geom = row.get('Truth_Geometry', '?')
        chrom = row.get('Truth_Chromosome', '?')
        size_str = f"{size:.0f}" if isinstance(size, (int, float)) and not np.isnan(size) else '?'
        print(f"    {truth} → {detected}  ({chrom}, {size_str} Mb, {geom}, {cf:.0f}% CF)")


def plot_effective_sensitivity(df):
    """
    Side-by-side curves: strict TP sensitivity vs TP + partial match sensitivity.
    Shows the gap between 'detected correctly' and 'detected at all'.
    """
    event_types = sorted(df['Truth_Type'].unique())
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.suptitle('Sensitivity: strict TP vs. including partial matches', 
                 fontsize=14, fontweight='bold')

    for ax_idx, etype in enumerate(event_types):
        ax = axes[ax_idx]
        subset = df[df['Truth_Type'] == etype]
        
        cfs = sorted(subset['CF_pct'].unique())
        cfs = [c for c in cfs if c >= 5]
        
        strict_sens = []
        effective_sens = []
        
        for cf in cfs:
            cf_data = subset[subset['CF_pct'] == cf]
            total = len(cf_data)
            tp = len(cf_data[cf_data['Outcome_simple'] == 'TP'])
            partial = len(cf_data[cf_data['Outcome_simple'].isin(['Type mismatch', 'Boundary mismatch'])])
            
            strict_sens.append(tp / total * 100 if total > 0 else 0)
            effective_sens.append((tp + partial) / total * 100 if total > 0 else 0)
        
        ax.plot(cfs, strict_sens, 'o-', color='#2ca02c', linewidth=2, markersize=7, 
                label='Strict TP', zorder=3)
        ax.plot(cfs, effective_sens, 's--', color='#ff7f0e', linewidth=2, markersize=7, 
                label='TP + partial match', zorder=3)
        ax.fill_between(cfs, strict_sens, effective_sens, alpha=0.15, color='#ff7f0e')
        
        # Add gap annotations at key CFs
        for cf, s, e in zip(cfs, strict_sens, effective_sens):
            gap = e - s
            if gap >= 3:
                ax.annotate(f'+{gap:.0f}%', xy=(cf, (s + e) / 2), fontsize=7,
                           ha='left', va='center', color='#cc6600', fontweight='bold',
                           xytext=(5, 0), textcoords='offset points')
        
        ax.axhline(50, color='grey', linestyle=':', alpha=0.5, label='50% threshold')
        ax.set_title(etype, fontsize=13, fontweight='bold')
        ax.set_xlabel('Cell Fraction (%)')
        ax.set_ylabel('Sensitivity (%)' if ax_idx == 0 else '')
        ax.set_ylim(-2, 105)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color = grey3, length = 6)
        ax.yaxis.set_tick_params(width=1, color = grey3, length = 6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'effective_sensitivity_comparison.pdf'),
                dpi=300, bbox_inches='tight')
    print("✅ Saved effective sensitivity comparison")
    plt.close()


def plot_fpr_analysis(df):
    """
    FPR summary figure for the unphased caller, analogous to the phased caller's
    fpr_analysis_summary figure.  Uses within-sample false positives (calls on
    chromosomes where no mCA was simulated).

    Panels:
      top-left  : FPR (%) by cell fraction
      top-right : FPR (%) by event type
      bot-left  : FPR (%) by event size
      bot-right : FPR (%) by chromosome (bar chart; replaces p-value histogram
                  which does not apply to the deterministic unphased caller)
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('False Positive Rate – Unphased Caller (within-sample)',
                 fontsize=14, fontweight='bold')

    alpha_line_kw  = dict(color='gray', linestyle='--', linewidth=1.2, alpha=0.7)
    bar_label_kw   = dict(ha='center', va='bottom', fontsize=9, fontweight='bold')

    def _fpr_bar(ax, groups, labels, colors, xlabel, title):
        """Helper: draw a labelled FPR bar chart."""
        x = np.arange(len(groups))
        bars = ax.bar(x, [g['fpr'] for g in groups], color=colors,
                      edgecolor='white', linewidth=0.8)
        ax.axhline(5, **alpha_line_kw, label='α = 5%')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha='right')
        ax.set_ylabel('FPR (%)', fontweight='bold')
        ax.set_title(title, fontweight='bold', fontsize=12)
        ax.set_ylim(0, max(14, max(g['fpr'] for g in groups) * 1.4))
        ax.legend(fontsize=9)
        ax.grid(axis='y', alpha=0.3)
        for bar, g in zip(bars, groups):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.2,
                    f"{g['n_fp']}/{g['n_total']}\n({g['fpr']:.1f}%)",
                    **bar_label_kw)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # ── Panel 1: FPR by cell fraction ──────────────────────────────────────
    ax = axes[0, 0]
    cf_groups = []
    for cf in sorted(df['CF_pct'].unique()):
        sub = df[df['CF_pct'] == cf]
        n_fp    = sub['False_Positives'].sum()
        n_total = len(sub)
        cf_groups.append({'cf': cf, 'n_fp': n_fp, 'n_total': n_total,
                          'fpr': n_fp / n_total * 100 if n_total else 0})
    cf_labels = [f"{int(g['cf'])}%" for g in cf_groups]
    cf_colors = plt.cm.viridis_r(np.linspace(0.15, 0.85, len(cf_groups)))
    _fpr_bar(ax, cf_groups, cf_labels, cf_colors, 'Cell Fraction', 'FPR by Cell Fraction')
    ax.tick_params(axis='x', rotation=45)
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha('right')

    # ── Panel 2: FPR by event type ─────────────────────────────────────────
    ax = axes[0, 1]
    type_colors_map = {'CN-LOH': '#F77F00', 'GAIN': '#2E86AB', 'LOSS': '#06A77D'}
    type_groups = []
    type_labels = []
    type_colors = []
    for etype in sorted(df['Truth_Type'].unique()):
        sub = df[df['Truth_Type'] == etype]
        n_fp    = sub['False_Positives'].sum()
        n_total = len(sub)
        type_groups.append({'n_fp': n_fp, 'n_total': n_total,
                            'fpr': n_fp / n_total * 100 if n_total else 0})
        type_labels.append(etype)
        type_colors.append(type_colors_map.get(etype, '#888888'))
    _fpr_bar(ax, type_groups, type_labels, type_colors, 'Event Type', 'FPR by Event Type')

    # ── Panel 3: FPR by event size ─────────────────────────────────────────
    ax = axes[1, 0]
    size_order  = ['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    size_colors_list = ['#E63946', '#F77F00', '#FCBF49', '#06A77D', '#2E86AB', '#7B2D8E']
    size_groups = []
    size_labels = []
    size_cols   = []
    for size, col in zip(size_order, size_colors_list):
        sub = df[df['Truth_Length_Category'] == size]
        if len(sub) == 0:
            continue
        n_fp    = sub['False_Positives'].sum()
        n_total = len(sub)
        size_groups.append({'n_fp': n_fp, 'n_total': n_total,
                            'fpr': n_fp / n_total * 100 if n_total else 0})
        size_labels.append(size)
        size_cols.append(col)
    _fpr_bar(ax, size_groups, size_labels, size_cols, 'Event Size', 'FPR by Event Size')

    # ── Panel 4: FPR by chromosome ─────────────────────────────────────────
    ax = axes[1, 1]
    chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    chrom_data  = []
    for chrom in chrom_order:
        sub = df[df['Truth_Chromosome'] == chrom]
        if len(sub) == 0:
            continue
        n_fp    = sub['False_Positives'].sum()
        n_total = len(sub)
        chrom_data.append({
            'chrom': chrom.replace('chr', ''),
            'n_fp': n_fp, 'n_total': n_total,
            'fpr': n_fp / n_total * 100 if n_total else 0
        })
    chrom_labels = [d['chrom'] for d in chrom_data]
    chrom_fprs   = [d['fpr']   for d in chrom_data]
    mean_fpr = np.mean(chrom_fprs) if chrom_fprs else 0
    chrom_colors = ['#E63946' if f > 5 else '#2E86AB' for f in chrom_fprs]
    bars = ax.bar(chrom_labels, chrom_fprs, color=chrom_colors,
                  edgecolor='white', linewidth=0.5)
    ax.axhline(5,       **alpha_line_kw, label='α = 5%')
    ax.axhline(mean_fpr, color='#E63946', linestyle=':', linewidth=1.5,
               label=f'Mean: {mean_fpr:.1f}%')
    ax.set_ylabel('FPR (%)', fontweight='bold')
    ax.set_title('FPR by Chromosome', fontweight='bold', fontsize=12)
    ax.set_ylim(0, max(14, max(chrom_fprs) * 1.4) if chrom_fprs else 14)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.tick_params(axis='x', labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for ax in axes.flatten():
        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fpr_analysis_summary.pdf'),
                dpi=300, bbox_inches='tight')
    print("✅ Saved FPR analysis summary")
    plt.close()


def plot_roc_operating_points(df):
    """
    Sensitivity vs FPR operating-point plot for the unphased caller.

    Because the unphased caller is a deterministic pipeline (no p-value threshold
    to sweep), we cannot draw a continuous ROC curve.  Instead, each cell fraction
    defines one operating point.  Plotting these across CFs gives the trajectory
    through ROC space as clone size increases — equivalent to the phased ROC
    curves but showing discrete operating points rather than a swept curve.

    One panel per event type; points coloured by cell fraction.
    """
    event_types = sorted(df['Truth_Type'].unique())
    cf_values   = sorted(df['CF_pct'].unique())

    # Colour map: low CF = red/orange, high CF = blue
    cmap   = plt.cm.plasma
    norm   = plt.Normalize(vmin=min(cf_values), vmax=max(cf_values))

    fig, axes = plt.subplots(1, len(event_types),
                             figsize=(6 * len(event_types), 5), sharey=True)
    if len(event_types) == 1:
        axes = [axes]
    fig.suptitle('Sensitivity vs FPR Operating Points – Unphased Caller\n'
                 '(each point = one cell fraction; arrow shows increasing CF)',
                 fontsize=13, fontweight='bold')

    for ax, etype in zip(axes, event_types):
        sub = df[df['Truth_Type'] == etype]
        tprs = []
        fprs = []
        cfs_plot = []

        for cf in cf_values:
            cf_sub  = sub[sub['CF_pct'] == cf]
            n_total = len(cf_sub)
            if n_total == 0:
                continue
            n_tp = cf_sub['Detected'].sum()
            n_fp = cf_sub['False_Positives'].sum()
            # FPR denominator: total simulated events at this CF (same as sensitivity denom)
            tpr = n_tp / n_total * 100
            fpr = n_fp / n_total * 100
            tprs.append(tpr)
            fprs.append(fpr)
            cfs_plot.append(cf)

        # Diagonal chance line
        ax.plot([0, 100], [0, 100], color='gray', linestyle='--',
                linewidth=1, alpha=0.5, label='Chance', zorder=1)
        ax.axhline(50, color='gray', linestyle=':', linewidth=1, alpha=0.4)

        # Draw connecting line first (behind points)
        ax.plot(fprs, tprs, color='#cccccc', linewidth=1.5, zorder=2)

        # Annotate arrows along trajectory to show direction
        for i in range(len(fprs) - 1):
            dx = fprs[i+1] - fprs[i]
            dy = tprs[i+1] - tprs[i]
            if abs(dx) + abs(dy) > 1:
                ax.annotate('', xy=(fprs[i+1], tprs[i+1]),
                            xytext=(fprs[i], tprs[i]),
                            arrowprops=dict(arrowstyle='->', color='#aaaaaa',
                                            lw=1.2), zorder=3)

        # Scatter points coloured by CF
        sc = ax.scatter(fprs, tprs, c=cfs_plot, cmap=cmap, norm=norm,
                        s=80, zorder=4, edgecolors='white', linewidths=0.8)

        # Label each point with its CF
        for fpr_val, tpr_val, cf in zip(fprs, tprs, cfs_plot):
            offset_x = 1.5
            offset_y = 2
            ax.annotate(f'{cf:.0f}%',
                        xy=(fpr_val, tpr_val),
                        xytext=(fpr_val + offset_x, tpr_val + offset_y),
                        fontsize=7, color=cmap(norm(cf)),
                        fontweight='bold')

        ax.set_xlim(-2, max(max(fprs) * 1.5, 15))
        ax.set_ylim(-5, 105)
        ax.set_xlabel('False Positive Rate (%)', fontweight='bold')
        ax.set_title(etype, fontweight='bold', fontsize=13)
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(1.5)
            ax.spines[axis].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=6)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=6)

    axes[0].set_ylabel('Sensitivity (TPR) (%)', fontweight='bold')

    # Shared colourbar for CF
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.7, pad=0.02)
    cbar.set_label('Cell Fraction (%)', fontweight='bold')

    plt.savefig(os.path.join(OUTPUT_DIR, 'roc_operating_points.pdf'),
                dpi=300, bbox_inches='tight')
    print("✅ Saved ROC operating points")
    plt.close()


def plot_combined_sensitivity_figure(df_all, df_5mb):
    """
    Publication-ready 6-panel figure (a–f):

      Row 1:  a  b  c  d
      Row 2:  e (narrow)  |  f (wide chromosome heatmap)

    Panel a  – Sensitivity by event size          (ALL sizes)
    Panel b  – Overall sensitivity vs CF          (≥5 Mb)
    Panel c  – Sensitivity by event type          (≥5 Mb)
    Panel d  – Sensitivity by geometry            (≥5 Mb)
    Panel e  – Precision + FP/sample twin-axis    (≥5 Mb)
    Panel f  – Chromosome × CF sensitivity heatmap (≥5 Mb, ≥10% CF)

    All non-title text: Helvetica 6 pt
    Titles:             Helvetica 7 pt
    """
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42

    # ── font constants ──────────────────────────────────────────────────────
    SZ     = 10    # body: axis labels, tick labels, legends, annotations
    SZ_TTL = 12   # panel titles
    SZ_LBL = 13   # bold panel letters (a, b, c …)

    fig = plt.figure(figsize=(22, 10))
    gs  = gridspec.GridSpec(
        2, 4,
        figure=fig,
        hspace=0.38, wspace=0.35,
        height_ratios=[1, 1.05],
    )

    # ── helpers ─────────────────────────────────────────────────────────────
    def _style(ax):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for sp in ['bottom', 'left']:
            ax.spines[sp].set_linewidth(1.5)
            ax.spines[sp].set_color(grey3)
        ax.xaxis.set_tick_params(width=1, color=grey3, length=5)
        ax.yaxis.set_tick_params(width=1, color=grey3, length=5)
        ax.tick_params(labelsize=SZ)
        ax.xaxis.label.set_fontsize(SZ)
        ax.yaxis.label.set_fontsize(SZ)
        ax.title.set_fontsize(SZ_TTL)

    panel_labels = iter('abcdef')

    def _label(ax):
        ax.text(-0.12, 1.05, next(panel_labels),
                transform=ax.transAxes, fontsize=SZ_LBL,
                fontweight='bold', va='top', ha='left')

    # ── Panel a: Sensitivity by event size (ALL sizes) ──────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    size_order  = ['2Mb', '3Mb', '5Mb', '10Mb', '20Mb', 'WholeArm']
    size_colors = ['#E63946', '#F77F00', '#FCBF49', '#06A77D', '#2E86AB', '#7B2D8E']
    for size, color in zip(size_order, size_colors):
        sub = df_all[df_all['Truth_Length_Category'] == size]
        if len(sub) == 0:
            continue
        fm = metrics_by_group(sub, 'Truth_Fraction', 'Fraction')
        fm['Fraction_Pct'] = fm['Fraction'] * 100
        ax_a.plot(fm['Fraction_Pct'], fm['Sensitivity'] * 100,
                  'o-', linewidth=2, markersize=6, label=size, color=color)
    ax_a.axhline(50, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax_a.set_xlabel('Cell Fraction (%)', fontsize=SZ, fontweight='bold')
    ax_a.set_ylabel('Sensitivity (%)', fontsize=SZ, fontweight='bold')
    ax_a.set_title('Sensitivity by Event Size', fontsize=SZ_TTL, fontweight='bold')
    ax_a.set_ylim(-5, 105)
    ax_a.legend(frameon=True, fontsize=SZ, loc='upper left')
    ax_a.grid(True, alpha=0.3)
    _style(ax_a); _label(ax_a)

    # ── Panel b: Overall sensitivity vs CF (≥5 Mb) ──────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    fm_b = metrics_by_group(df_5mb, 'Truth_Fraction', 'Fraction')
    fm_b['Fraction_Pct'] = fm_b['Fraction'] * 100
    ax_b.plot(fm_b['Fraction_Pct'], fm_b['Sensitivity'] * 100,
              'o-', linewidth=2, markersize=8, color='#2E86AB')
    ax_b.axhline(50, color='gray', linestyle='--', alpha=0.5, linewidth=1,
                 label='50% sensitivity')
    for _, row in fm_b.iterrows():
        if row['Sensitivity'] > 0.01:
            ax_b.text(row['Fraction_Pct'], row['Sensitivity'] * 100 + 3,
                      f"{int(row['True_Positives'])}/{int(row['Total'])}",
                      ha='center', fontsize=SZ, fontweight='bold')
    ax_b.set_xlabel('Cell Fraction (%)', fontsize=SZ, fontweight='bold')
    ax_b.set_ylabel('Sensitivity (%)', fontsize=SZ, fontweight='bold')
    ax_b.set_title('Overall Sensitivity vs Cell Fraction', fontsize=SZ_TTL, fontweight='bold')
    ax_b.set_ylim(-5, 105)
    ax_b.legend(fontsize=SZ)
    ax_b.grid(True, alpha=0.3)
    _style(ax_b); _label(ax_b)

    # ── Panel c: Sensitivity by event type (≥5 Mb) ──────────────────────────
    ax_c = fig.add_subplot(gs[0, 2])
    type_colors = {'CNLOH': '#2E86AB', 'GAIN': '#F77F00', 'LOSS': '#06A77D'}
    for etype in sorted(df_5mb['Truth_Type'].unique()):
        sub = df_5mb[df_5mb['Truth_Type'] == etype]
        fm  = metrics_by_group(sub, 'Truth_Fraction', 'Fraction')
        fm['Fraction_Pct'] = fm['Fraction'] * 100
        ax_c.plot(fm['Fraction_Pct'], fm['Sensitivity'] * 100,
                  'o-', linewidth=2, markersize=6,
                  label=etype, color=type_colors.get(etype, '#888888'))
    ax_c.axhline(50, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax_c.set_xlabel('Cell Fraction (%)', fontsize=SZ, fontweight='bold')
    ax_c.set_ylabel('Sensitivity (%)', fontsize=SZ, fontweight='bold')
    ax_c.set_title('Sensitivity by Event Type', fontsize=SZ_TTL, fontweight='bold')
    ax_c.set_ylim(-5, 105)
    ax_c.legend(frameon=True, fontsize=SZ)
    ax_c.grid(True, alpha=0.3)
    _style(ax_c); _label(ax_c)

    # ── Panel d: Sensitivity by geometry (≥5 Mb) ────────────────────────────
    ax_d = fig.add_subplot(gs[0, 3])
    geom_colors = {'Interstitial': '#2E86AB', 'Telomeric': '#F77F00', 'Whole_Arm': '#06A77D'}
    for geom in sorted(df_5mb['Truth_Geometry'].unique()):
        sub = df_5mb[df_5mb['Truth_Geometry'] == geom]
        fm  = metrics_by_group(sub, 'Truth_Fraction', 'Fraction')
        fm['Fraction_Pct'] = fm['Fraction'] * 100
        ax_d.plot(fm['Fraction_Pct'], fm['Sensitivity'] * 100,
                  'o-', linewidth=2, markersize=6,
                  label=geom, color=geom_colors.get(geom, '#888888'))
    ax_d.axhline(50, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax_d.set_xlabel('Cell Fraction (%)', fontsize=SZ, fontweight='bold')
    ax_d.set_ylabel('Sensitivity (%)', fontsize=SZ, fontweight='bold')
    ax_d.set_title('Sensitivity by Geometry', fontsize=SZ_TTL, fontweight='bold')
    ax_d.set_ylim(-5, 105)
    ax_d.legend(frameon=True, fontsize=SZ)
    ax_d.grid(True, alpha=0.3)
    _style(ax_d); _label(ax_d)

    # ── Panel e: Precision + FP/sample twin-axis (≥5 Mb) ────────────────────
    ax_e = fig.add_subplot(gs[1, 0])
    fm_e = metrics_by_group(df_5mb, 'Truth_Fraction', 'Fraction')
    fm_e['Fraction_Pct'] = fm_e['Fraction'] * 100
    prec_color = '#06A77D'
    fp_color   = '#E63946'

    ax_e.plot(fm_e['Fraction_Pct'], fm_e['Precision'] * 100,
              's-', linewidth=2, markersize=8, color=prec_color,
              label='Precision (PPV)', zorder=3)
    ax_e.axhline(50, color='gray', linestyle='--', alpha=0.4, linewidth=1)
    ax_e.set_xlabel('Cell Fraction (%)', fontsize=SZ, fontweight='bold')
    ax_e.set_ylabel('Precision (PPV) (%)', fontsize=SZ, fontweight='bold', color=prec_color)
    ax_e.tick_params(axis='y', labelcolor=prec_color, labelsize=SZ)
    ax_e.set_ylim(-5, 105)
    ax_e.set_title('Precision & False Positive Rate', fontsize=SZ_TTL, fontweight='bold')
    ax_e.grid(True, alpha=0.3)

    ax_e2 = ax_e.twinx()
    ax_e2.plot(fm_e['Fraction_Pct'], fm_e['FP_per_sample'],
               '^--', linewidth=2, markersize=8, color=fp_color,
               label='FP / sample', zorder=3)
    fp_max = fm_e['FP_per_sample'].max()
    ax_e2.set_ylim(0, max(fp_max * 6, 1.0))
    ax_e2.set_ylabel('False Positives per Sample', fontsize=SZ, fontweight='bold',
                     color=fp_color)
    ax_e2.tick_params(axis='y', labelcolor=fp_color, labelsize=SZ)
    peak_row = fm_e.loc[fm_e['FP_per_sample'].idxmax()]
    ax_e2.annotate(f"peak {peak_row['FP_per_sample']:.2f}",
                   xy=(peak_row['Fraction_Pct'], peak_row['FP_per_sample']),
                   xytext=(peak_row['Fraction_Pct'] + 6,
                           peak_row['FP_per_sample'] * 1.8),
                   fontsize=SZ, fontweight='bold', color=fp_color,
                   arrowprops=dict(arrowstyle='->', color=fp_color, lw=1.0))

    lines1, labs1 = ax_e.get_legend_handles_labels()
    lines2, labs2 = ax_e2.get_legend_handles_labels()
    ax_e.legend(lines1 + lines2, labs1 + labs2,
                frameon=True, fontsize=SZ, loc='center right')

    _style(ax_e)
    ax_e2.spines['top'].set_visible(False)
    ax_e2.spines['right'].set_linewidth(1.5)
    ax_e2.spines['right'].set_color(grey3)
    ax_e2.yaxis.set_tick_params(width=1, color=grey3, length=5)
    _label(ax_e)

    # ── Panel f: Chromosome × CF heatmap (≥5 Mb, ≥10% CF) ──────────────────
    ax_f = fig.add_subplot(gs[1, 1:])

    chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    detectable  = df_5mb[df_5mb['Truth_Fraction'] >= 0.10]
    pivot = detectable.pivot_table(
        index='Truth_Fraction',
        columns='Truth_Chromosome',
        values='Detected',
        aggfunc='mean'
    ) * 100
    pivot.index = [f"{x*100:.0f}%" for x in pivot.index]
    pivot = pivot.reindex(columns=[c for c in chrom_order if c in pivot.columns])
    pivot.columns = [c.replace('chr', '') for c in pivot.columns]
    pivot = pivot.iloc[::-1]  # highest CF at top

    sns.heatmap(pivot, annot=True, fmt='.0f', cmap='RdYlBu',
                vmin=0, vmax=100,
                cbar_kws={'label': 'Sensitivity (%)', 'shrink': 0.8},
                ax=ax_f, linewidths=0.4,
                annot_kws={'size': SZ})
    ax_f.set_title('Sensitivity by Chromosome and Cell Fraction (≥10%)',
                   fontsize=SZ_TTL, fontweight='bold')
    ax_f.set_xlabel('Chromosome', fontsize=SZ, fontweight='bold')
    ax_f.set_ylabel('Cell Fraction', fontsize=SZ, fontweight='bold')
    ax_f.tick_params(axis='both', labelsize=SZ)
    ax_f.set_xticklabels(ax_f.get_xticklabels(), rotation=0)
    ax_f.set_yticklabels(ax_f.get_yticklabels(), rotation=0)
    ax_f.xaxis.set_tick_params(width=1, color=grey3, length=4)
    ax_f.yaxis.set_tick_params(width=1, color=grey3, length=4)

    # Colorbar font
    cbar = ax_f.collections[0].colorbar
    cbar.ax.tick_params(labelsize=SZ)
    cbar.set_label('Sensitivity (%)', fontsize=SZ)

    for sp in ['bottom', 'left', 'top', 'right']:
        ax_f.spines[sp].set_linewidth(1.5)
        ax_f.spines[sp].set_color(grey3)

    ax_f.text(-0.04, 1.05, next(panel_labels),
              transform=ax_f.transAxes, fontsize=SZ_LBL,
              fontweight='bold', va='top', ha='left')

    plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_analysis_combined.pdf'),
                dpi=300, bbox_inches='tight')
    print("✅ Saved combined sensitivity figure")
    plt.close()

    # Reset rcParams to defaults so other figures are unaffected
    mpl.rcParams.update(mpl.rcParamsDefault)
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype']  = 42


def print_cnloh_deep_dive(df):
    """Print detailed CN-LOH sensitivity breakdown to console."""
    cnloh = df[df['Truth_Type'] == 'CN-LOH'].copy()

    print("\n" + "=" * 60)
    print("CN-LOH SENSITIVITY DEEP DIVE")
    print("=" * 60)

    for cf in sorted(cnloh['CF_pct'].unique()):
        if cf < 10:
            continue
        cf_data = cnloh[cnloh['CF_pct'] == cf]
        total = len(cf_data)
        tp = len(cf_data[cf_data['Outcome_simple'] == 'TP'])
        type_mm = len(cf_data[cf_data['Outcome_simple'] == 'Type mismatch'])
        boundary = len(cf_data[cf_data['Outcome_simple'] == 'Boundary mismatch'])
        fn = len(cf_data[cf_data['Outcome_simple'] == 'Not detected'])

        print(f"\n  {cf:.0f}% CF: {tp}/{total} TP ({tp/total*100:.0f}%) | "
              f"{type_mm} type mismatch, {boundary} boundary, {fn} truly missed")

        tm = cf_data[cf_data['Outcome_simple'] == 'Type mismatch']
        if len(tm) > 0:
            for _, row in tm.iterrows():
                size = row.get('Truth_Length_Mb', '?')
                geom = row.get('Truth_Geometry', '?')
                size_str = f"{size:.0f}" if isinstance(size, (int, float)) and not np.isnan(size) else '?'
                print(f"    → called as {row['Detected_Type']} (truth: {size_str} Mb {geom})")

    # Size dependence
    print("\n  CN-LOH detection by size at key CFs:")
    for size_cat in sorted(cnloh['Truth_Length_Mb'].unique()):
        size_data = cnloh[cnloh['Truth_Length_Mb'] == size_cat]
        parts = []
        for cf in [25, 50, 75, 100]:
            cf_data = size_data[size_data['CF_pct'] == cf]
            if len(cf_data) == 0:
                continue
            tp = len(cf_data[cf_data['Outcome_simple'] == 'TP'])
            total = len(cf_data)
            partial = len(cf_data[cf_data['Outcome_simple'].isin(['Type mismatch', 'Boundary mismatch'])])
            parts.append(f"{cf}%: {tp}/{total} TP +{partial}p")
        if parts:
            print(f"    {size_cat:.0f} Mb: {' | '.join(parts)}")

    print("=" * 60)


# ============================================================
# MAIN
# ============================================================

def main():
    """Main analysis function"""
    
    print("=" * 80)
    print("ANALYZING mCA CALLER PERFORMANCE")
    print("=" * 80)
    
    # Load results
    print(f"\n📂 Loading: {TRUTH_FILE}")
    if not os.path.exists(TRUTH_FILE):
        print(f"❌ File not found: {TRUTH_FILE}")
        return
    
    df = pd.read_csv(TRUTH_FILE)
    print(f"✅ Loaded {len(df)} simulated events")
    
    # Add length in Mb if not present
    if 'Truth_Length_Mb' not in df.columns:
        df['Truth_Length_Mb'] = (df['Truth_End'] - df['Truth_Start']) / 1e6
    
    # Categorize length if not present
    if 'Truth_Length_Category' not in df.columns:
        def categorize_length(row):
            if row['Truth_Geometry'] == 'Whole_Arm':
                return 'WholeArm'
            length_mb = row['Truth_Length_Mb']
            if length_mb < 2.5:
                return '2Mb'
            elif length_mb < 3.5:
                return '3Mb'
            elif length_mb < 7.5:
                return '5Mb'
            elif length_mb < 15:
                return '10Mb'
            else:
                return '20Mb'
        df['Truth_Length_Category'] = df.apply(categorize_length, axis=1)
    
    # --- Add outcome classification columns for breakdown analysis ---
    df['CF_pct'] = df['Truth_Fraction'] * 100 if df['Truth_Fraction'].max() <= 1 else df['Truth_Fraction']
    if df['CF_pct'].max() <= 1.5:
        df['CF_pct'] = df['Truth_Fraction'] * 100
    df['Outcome'] = df.apply(classify_outcome, axis=1)
    df['Outcome_simple'] = df['Outcome'].replace({
        'Oversized': 'Boundary mismatch',
        'Undersized': 'Boundary mismatch',
        'Size mismatch': 'Boundary mismatch',
        'Other partial': 'Boundary mismatch',
    })

    # ==========================================
    # Generate figures — ALL events
    # ==========================================
    print("\n📊 Generating figures (all sizes)...")
    
    plot_sensitivity_curves(df)
    plot_heatmaps(df)
    plot_chromosome_analysis(df) 
    plot_summary_table(df)
    plot_detection_breakdown(df)
    plot_misclassification_matrix(df)
    plot_effective_sensitivity(df)
    plot_fpr_analysis(df)
    plot_roc_operating_points(df)
    
    # Rename all outputs to _all_sizes
    for f in os.listdir(OUTPUT_DIR):
        if f.endswith(('.pdf', '.png')) and '_5Mb_plus' not in f and '_all_sizes' not in f:
            base, ext = os.path.splitext(f)
            os.rename(os.path.join(OUTPUT_DIR, f), 
                      os.path.join(OUTPUT_DIR, f"{base}_all_sizes{ext}"))
    
    # ==========================================
    # Generate figures — ≥5 Mb only
    # ==========================================
    df_5mb = df[df['Truth_Length_Category'].isin(['5Mb', '10Mb', '20Mb', 'WholeArm'])]
    print(f"\n📊 Generating figures (≥5 Mb only): {len(df_5mb)} events...")
    
    plot_sensitivity_curves(df_5mb)
    plot_heatmaps(df_5mb)
    plot_chromosome_analysis(df_5mb)
    plot_summary_table(df_5mb)
    plot_detection_breakdown(df_5mb)
    plot_misclassification_matrix(df_5mb)
    plot_effective_sensitivity(df_5mb)
    plot_fpr_analysis(df_5mb)
    plot_roc_operating_points(df_5mb)
    plot_combined_sensitivity_figure(df, df_5mb)
    
    # Rename these to _5Mb_plus
    for f in os.listdir(OUTPUT_DIR):
        if f.endswith(('.pdf', '.png')) and '_all_sizes' not in f and '_5Mb_plus' not in f:
            base, ext = os.path.splitext(f)
            os.rename(os.path.join(OUTPUT_DIR, f), 
                      os.path.join(OUTPUT_DIR, f"{base}_5Mb_plus{ext}"))
    
    # ==========================================
    # CN-LOH deep dive (console output)
    # ==========================================
    print_cnloh_deep_dive(df)
    
    # ==========================================
    # Export summary statistics
    # ==========================================
    print("\n📝 Exporting summary statistics...")
    
    overall = calculate_metrics(df)
    frac_metrics = metrics_by_group(df, 'Truth_Fraction', 'Fraction')
    type_metrics = metrics_by_group(df, 'Truth_Type', 'Event Type')
    size_metrics = metrics_by_group(df, 'Truth_Length_Category', 'Size')
    
    # Partial match summary
    partial_summary = df.groupby(['Truth_Type', 'CF_pct', 'Outcome_simple']).size().reset_index(name='Count')
    
    with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'mCA_simulations_summary_metrics.xlsx')) as writer:
        pd.DataFrame([overall]).to_excel(writer, sheet_name='Overall', index=False)
        frac_metrics.to_excel(writer, sheet_name='By_Fraction', index=False)
        type_metrics.to_excel(writer, sheet_name='By_Type', index=False)
        size_metrics.to_excel(writer, sheet_name='By_Size', index=False)
        partial_summary.to_excel(writer, sheet_name='Detection_Breakdown', index=False)
    
    print(f"✅ Saved mCA_simulations_summary_metrics.xlsx (with Detection_Breakdown sheet)")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to: {OUTPUT_DIR}/")
    print("  Original figures:")
    print("    - sensitivity_analysis.pdf/png")
    print("    - mCA_simulations_sensitivity_heatmaps.pdf")
    print("    - mCA_simulations_chromosome_analysis.pdf")
    print("    - mCA_simulations_summary_table.pdf")
    print("  New breakdown figures:")
    print("    - detection_breakdown_by_type.pdf")
    print("    - misclassification_matrix.pdf")
    print("    - effective_sensitivity_comparison.pdf")
    print("    - fpr_analysis_summary.pdf")
    print("    - roc_operating_points.pdf")
    print("  Both _all_sizes and _5Mb_plus variants generated.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run mCA caller sensitivity analysis on simulated samples")
    parser.add_argument("results_dir", help="Path to mCA caller results for analysis")
    args = parser.parse_args()

    RESULTS_DIR = args.results_dir
    TRUTH_FILE = os.path.join(RESULTS_DIR, "truth_comparison.csv")
    OUTPUT_DIR = os.path.join(RESULTS_DIR, "Sensitivity_analysis_unphased")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    main()