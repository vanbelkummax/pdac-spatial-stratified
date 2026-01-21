#!/usr/bin/env python3
"""
PDAC Spatial Transcriptomics Analysis - Properly Stratified
============================================================

Patient Structure:
- YP12 (R): Pre + Post
- YP15 (R): Pre + Post
- YP03 (NR): Pre + Post
- YP04 (NR): Post only

Valid Comparisons:
- Post-treatment: R (n=2) vs NR (n=2) - best comparison
- Pre-treatment: R (n=2) vs NR (n=1) - limited
- Paired change: R (n=2) vs NR (n=1) - limited

NO P-VALUES with n=2 per group - descriptive only.
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_DIR = Path('/home/user/spatial-hackathon-2026/outputs/adata/polymathic')
OUTPUT_DIR = Path('/home/user/spatial-hackathon-2026-figures')
FIG_DIR = OUTPUT_DIR / 'figures'
TABLE_DIR = OUTPUT_DIR / 'tables'

# Patient-level metadata (NOT sample-level)
PATIENTS = {
    'YP12': {'response': 'R', 'pre': 'YP12A', 'post': 'YP12C'},
    'YP15': {'response': 'R', 'pre': 'YP15A', 'post': 'YP15C'},
    'YP03': {'response': 'NR', 'pre': 'YP03A', 'post': 'YP03C'},
    'YP04': {'response': 'NR', 'pre': None, 'post': 'YP04C'},  # NO PRE!
}

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


def load_sample(sample_id):
    """Load a single sample."""
    path = DATA_DIR / f'{sample_id}_polymathic.h5ad'
    if path.exists():
        return sc.read_h5ad(path)
    return None


def get_patient_metrics(patient_id, timepoint='post'):
    """Get metrics for a patient at specific timepoint."""
    info = PATIENTS[patient_id]
    sample_id = info[timepoint]
    if sample_id is None:
        return None

    adata = load_sample(sample_id)
    if adata is None:
        return None

    metrics = {
        'patient': patient_id,
        'response': info['response'],
        'timepoint': timepoint,
        'sample': sample_id,
        'n_spots': adata.n_obs,
    }

    # Cell type proportions
    if 'cell_type' in adata.obs.columns:
        ct_counts = adata.obs['cell_type'].value_counts()
        total = ct_counts.sum()
        for ct, count in ct_counts.items():
            metrics[f'pct_{ct}'] = 100 * count / total

    return metrics, adata


def create_figure1_overview():
    """Figure 1: Study overview and patient structure."""
    print("Creating Figure 1: Study Overview...")

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # A) Patient structure diagram
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)

    # Draw patient boxes
    positions = {'YP12': (2, 7), 'YP15': (2, 5), 'YP03': (7, 7), 'YP04': (7, 5)}

    for patient, (x, y) in positions.items():
        info = PATIENTS[patient]
        color = RESPONSE_COLORS[info['response']]

        # Patient box
        rect = plt.Rectangle((x-1, y-0.4), 2, 0.8, fill=True,
                             facecolor=color, alpha=0.3, edgecolor=color, linewidth=2)
        ax.add_patch(rect)
        ax.text(x, y, patient, ha='center', va='center', fontsize=12, fontweight='bold')

        # Pre/Post indicators
        if info['pre']:
            ax.text(x-0.5, y-0.6, 'Pre ✓', ha='center', fontsize=8, color='green')
        else:
            ax.text(x-0.5, y-0.6, 'Pre ✗', ha='center', fontsize=8, color='red')

        ax.text(x+0.5, y-0.6, 'Post ✓', ha='center', fontsize=8, color='green')

    # Labels
    ax.text(2, 8.5, 'RESPONDERS (n=2)', ha='center', fontsize=11, fontweight='bold', color=RESPONSE_COLORS['R'])
    ax.text(7, 8.5, 'NON-RESPONDERS (n=2)', ha='center', fontsize=11, fontweight='bold', color=RESPONSE_COLORS['NR'])

    ax.set_title('A) Patient Structure', fontsize=12)
    ax.axis('off')

    # B) Sample counts
    ax = axes[1]
    data = {
        'Pre-treatment': [2, 1],  # R, NR (YP04 has no pre)
        'Post-treatment': [2, 2],
    }
    x = np.arange(2)
    width = 0.35

    ax.bar(x - width/2, [2, 2], width, label='Responders', color=RESPONSE_COLORS['R'])
    ax.bar(x + width/2, [1, 2], width, label='Non-Responders', color=RESPONSE_COLORS['NR'])

    ax.set_ylabel('Number of Patients')
    ax.set_xticks(x)
    ax.set_xticklabels(['Pre-treatment', 'Post-treatment'])
    ax.legend()
    ax.set_title('B) Patients per Timepoint', fontsize=12)
    ax.set_ylim(0, 3)

    # Add note about YP04
    ax.text(0, 1.2, 'YP04 missing', ha='center', fontsize=8, style='italic', color='red')

    # C) Analysis approach
    ax = axes[2]
    ax.text(0.5, 0.9, 'Analysis Approach', ha='center', fontsize=12, fontweight='bold', transform=ax.transAxes)

    text = """
• Post-treatment comparison: R (n=2) vs NR (n=2)
  → Best powered comparison

• Pre-treatment comparison: R (n=2) vs NR (n=1)
  → Limited by missing YP04 Pre

• Paired analysis (Pre→Post change):
  R: 2 pairs, NR: 1 pair

⚠️ NO P-VALUES REPORTED
   n=2 per group insufficient for
   statistical inference

Descriptive statistics only:
  - Means ± SD
  - Effect sizes (Cohen's d)
  - Individual patient values shown
"""
    ax.text(0.1, 0.8, text, ha='left', va='top', fontsize=9, transform=ax.transAxes,
           family='monospace')
    ax.axis('off')
    ax.set_title('C) Statistical Considerations', fontsize=12)

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig01_study_overview.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved fig01_study_overview.png")


def create_figure2_celltype_post():
    """Figure 2: Cell type proportions - POST-TREATMENT only (n=2 vs n=2)."""
    print("Creating Figure 2: Cell Types (Post-treatment)...")

    # Collect post-treatment data
    data = []
    for patient_id in PATIENTS.keys():
        result = get_patient_metrics(patient_id, 'post')
        if result:
            metrics, _ = result
            data.append(metrics)

    df = pd.DataFrame(data)

    # Get cell type columns
    ct_cols = [c for c in df.columns if c.startswith('pct_')]

    # Prepare for plotting
    plot_data = []
    for _, row in df.iterrows():
        for col in ct_cols:
            if col in row and not pd.isna(row[col]):
                plot_data.append({
                    'patient': row['patient'],
                    'response': row['response'],
                    'cell_type': col.replace('pct_', ''),
                    'percentage': row[col]
                })

    plot_df = pd.DataFrame(plot_data)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # A) Stacked bar by patient
    ax = axes[0]
    patients_order = ['YP12', 'YP15', 'YP03', 'YP04']
    ct_pivot = plot_df.pivot(index='patient', columns='cell_type', values='percentage')
    ct_pivot = ct_pivot.reindex(patients_order)

    ct_pivot.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')

    # Color patient labels by response
    for i, patient in enumerate(patients_order):
        color = RESPONSE_COLORS[PATIENTS[patient]['response']]
        ax.get_xticklabels()[i].set_color(color)
        ax.get_xticklabels()[i].set_fontweight('bold')

    ax.set_ylabel('Percentage (%)')
    ax.set_xlabel('Patient')
    ax.set_title('A) Cell Type Composition (Post-treatment)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.set_xticklabels(patients_order, rotation=0)

    # B) R vs NR comparison for key cell types
    ax = axes[1]

    # Calculate means by response
    summary = plot_df.groupby(['response', 'cell_type'])['percentage'].agg(['mean', 'std']).reset_index()

    # Select top variable cell types
    ct_variance = plot_df.groupby('cell_type')['percentage'].var().sort_values(ascending=False)
    top_cts = ct_variance.head(8).index.tolist()

    summary_top = summary[summary['cell_type'].isin(top_cts)]

    # Grouped bar plot
    x = np.arange(len(top_cts))
    width = 0.35

    r_data = summary_top[summary_top['response'] == 'R'].set_index('cell_type').reindex(top_cts)
    nr_data = summary_top[summary_top['response'] == 'NR'].set_index('cell_type').reindex(top_cts)

    ax.bar(x - width/2, r_data['mean'], width, yerr=r_data['std'],
           label='R (n=2)', color=RESPONSE_COLORS['R'], capsize=3)
    ax.bar(x + width/2, nr_data['mean'], width, yerr=nr_data['std'],
           label='NR (n=2)', color=RESPONSE_COLORS['NR'], capsize=3)

    ax.set_ylabel('Percentage (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(top_cts, rotation=45, ha='right')
    ax.legend()
    ax.set_title('B) R vs NR (Post-treatment, mean ± SD)')

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig02_celltype_post.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved fig02_celltype_post.png")

    # Save data to Excel
    excel_data = {
        'patient': [], 'response': [], 'timepoint': [], 'cell_type': [], 'percentage': []
    }
    for _, row in plot_df.iterrows():
        excel_data['patient'].append(row['patient'])
        excel_data['response'].append(row['response'])
        excel_data['timepoint'].append('Post')
        excel_data['cell_type'].append(row['cell_type'])
        excel_data['percentage'].append(row['percentage'])

    pd.DataFrame(excel_data).to_excel(TABLE_DIR / 'T1_celltype_post.xlsx', index=False)
    print("  Saved T1_celltype_post.xlsx")


def create_figure3_celltype_pre():
    """Figure 3: Cell type proportions - PRE-TREATMENT (n=2 vs n=1)."""
    print("Creating Figure 3: Cell Types (Pre-treatment)...")

    # Collect pre-treatment data (YP04 excluded - no pre)
    data = []
    for patient_id in ['YP12', 'YP15', 'YP03']:  # Exclude YP04
        result = get_patient_metrics(patient_id, 'pre')
        if result:
            metrics, _ = result
            data.append(metrics)

    df = pd.DataFrame(data)

    # Get cell type columns
    ct_cols = [c for c in df.columns if c.startswith('pct_')]

    # Prepare for plotting
    plot_data = []
    for _, row in df.iterrows():
        for col in ct_cols:
            if col in row and not pd.isna(row[col]):
                plot_data.append({
                    'patient': row['patient'],
                    'response': row['response'],
                    'cell_type': col.replace('pct_', ''),
                    'percentage': row[col]
                })

    plot_df = pd.DataFrame(plot_data)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # A) Bar by patient
    ax = axes[0]
    patients_order = ['YP12', 'YP15', 'YP03']
    ct_pivot = plot_df.pivot(index='patient', columns='cell_type', values='percentage')
    ct_pivot = ct_pivot.reindex(patients_order)

    ct_pivot.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')

    for i, patient in enumerate(patients_order):
        color = RESPONSE_COLORS[PATIENTS[patient]['response']]
        ax.get_xticklabels()[i].set_color(color)
        ax.get_xticklabels()[i].set_fontweight('bold')

    ax.set_ylabel('Percentage (%)')
    ax.set_xlabel('Patient')
    ax.set_title('A) Cell Type Composition (Pre-treatment)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.set_xticklabels(patients_order, rotation=0)

    # B) Individual values with limited comparison
    ax = axes[1]

    # Show individual patient values for top cell types
    ct_variance = plot_df.groupby('cell_type')['percentage'].var().sort_values(ascending=False)
    top_cts = ct_variance.head(6).index.tolist()

    for i, ct in enumerate(top_cts):
        ct_data = plot_df[plot_df['cell_type'] == ct]

        for _, row in ct_data.iterrows():
            color = RESPONSE_COLORS[row['response']]
            marker = 'o' if row['response'] == 'R' else 's'
            ax.scatter(i, row['percentage'], c=color, marker=marker, s=100,
                      edgecolor='black', linewidth=1, zorder=5)
            ax.annotate(row['patient'], (i, row['percentage']),
                       xytext=(5, 0), textcoords='offset points', fontsize=8)

    ax.set_xticks(range(len(top_cts)))
    ax.set_xticklabels(top_cts, rotation=45, ha='right')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('B) Individual Patient Values\n(R: n=2, NR: n=1 - YP04 excluded)')

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=RESPONSE_COLORS['R'],
               markersize=10, label='R (n=2)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=RESPONSE_COLORS['NR'],
               markersize=10, label='NR (n=1)')
    ]
    ax.legend(handles=legend_elements)

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig03_celltype_pre.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved fig03_celltype_pre.png")

    # Save data
    pd.DataFrame(plot_data).to_excel(TABLE_DIR / 'T2_celltype_pre.xlsx', index=False)
    print("  Saved T2_celltype_pre.xlsx")


def create_figure4_paired_change():
    """Figure 4: Paired Pre→Post changes (n=2 R pairs, n=1 NR pair)."""
    print("Creating Figure 4: Paired Changes...")

    # Get paired data
    paired_data = []
    for patient_id in ['YP12', 'YP15', 'YP03']:  # Only patients with both
        pre_result = get_patient_metrics(patient_id, 'pre')
        post_result = get_patient_metrics(patient_id, 'post')

        if pre_result and post_result:
            pre_metrics, _ = pre_result
            post_metrics, _ = post_result

            ct_cols = [c for c in pre_metrics.keys() if c.startswith('pct_')]
            for col in ct_cols:
                if col in pre_metrics and col in post_metrics:
                    paired_data.append({
                        'patient': patient_id,
                        'response': PATIENTS[patient_id]['response'],
                        'cell_type': col.replace('pct_', ''),
                        'pre': pre_metrics[col],
                        'post': post_metrics[col],
                        'change': post_metrics[col] - pre_metrics[col]
                    })

    df = pd.DataFrame(paired_data)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # A) Slope plot for select cell types
    ax = axes[0]

    top_cts = df.groupby('cell_type')['change'].apply(lambda x: abs(x).mean()).nlargest(6).index.tolist()

    for ct in top_cts:
        ct_data = df[df['cell_type'] == ct]
        for _, row in ct_data.iterrows():
            color = RESPONSE_COLORS[row['response']]
            ax.plot([0, 1], [row['pre'], row['post']], 'o-', color=color,
                   markersize=8, linewidth=2, alpha=0.7)
            ax.annotate(f"{row['patient']}\n{ct}", (1.05, row['post']), fontsize=7)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Pre', 'Post'])
    ax.set_ylabel('Percentage (%)')
    ax.set_title('A) Pre→Post Changes (Top 6 Variable Cell Types)')
    ax.set_xlim(-0.2, 1.5)

    # B) Change by response
    ax = axes[1]

    # Calculate mean change per cell type per response
    change_summary = df.groupby(['response', 'cell_type'])['change'].mean().reset_index()
    change_summary = change_summary.pivot(index='cell_type', columns='response', values='change')

    # Sort by absolute difference
    change_summary['abs_diff'] = abs(change_summary['R'] - change_summary['NR'])
    change_summary = change_summary.sort_values('abs_diff', ascending=True)
    change_summary = change_summary.drop('abs_diff', axis=1)

    y = range(len(change_summary))

    ax.barh([i - 0.2 for i in y], change_summary['R'], height=0.35,
           label='R (n=2 pairs)', color=RESPONSE_COLORS['R'])
    ax.barh([i + 0.2 for i in y], change_summary['NR'], height=0.35,
           label='NR (n=1 pair)', color=RESPONSE_COLORS['NR'])

    ax.set_yticks(y)
    ax.set_yticklabels(change_summary.index)
    ax.set_xlabel('Change in % (Post - Pre)')
    ax.axvline(0, color='black', linewidth=0.5)
    ax.legend()
    ax.set_title('B) Mean Change by Response\n(Descriptive only)')

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig04_paired_change.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved fig04_paired_change.png")

    # Save data
    df.to_excel(TABLE_DIR / 'T3_paired_changes.xlsx', index=False)
    print("  Saved T3_paired_changes.xlsx")


def create_figure5_spatial_post():
    """Figure 5: Spatial metrics - POST-TREATMENT (n=2 vs n=2)."""
    print("Creating Figure 5: Spatial Metrics (Post-treatment)...")

    metrics_data = []

    for patient_id in PATIENTS.keys():
        info = PATIENTS[patient_id]
        sample_id = info['post']
        adata = load_sample(sample_id)

        if adata is None:
            continue

        row = {
            'patient': patient_id,
            'response': info['response'],
            'n_spots': adata.n_obs,
        }

        # Spatial entropy if available
        if 'spatial_entropy' in adata.obs.columns:
            row['spatial_entropy'] = adata.obs['spatial_entropy'].mean()

        # Cell type diversity
        if 'cell_type' in adata.obs.columns:
            ct_counts = adata.obs['cell_type'].value_counts()
            ct_props = ct_counts / ct_counts.sum()
            row['diversity_shannon'] = -np.sum(ct_props * np.log(ct_props + 1e-10))
            row['n_cell_types'] = len(ct_counts)

        # Tissue metrics from spatial coords
        if 'spatial' in adata.obsm:
            from scipy.spatial import ConvexHull
            from scipy.spatial.distance import cdist

            coords = adata.obsm['spatial']
            try:
                hull = ConvexHull(coords)
                row['tissue_area'] = hull.volume
            except:
                pass

            # NN distance
            dists = cdist(coords, coords)
            np.fill_diagonal(dists, np.inf)
            row['mean_nn_dist'] = np.mean(dists.min(axis=1))

        metrics_data.append(row)

    df = pd.DataFrame(metrics_data)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot each metric
    metrics_to_plot = ['n_spots', 'diversity_shannon', 'spatial_entropy', 'mean_nn_dist']
    titles = ['A) Number of Spots', 'B) Shannon Diversity',
              'C) Spatial Entropy', 'D) Mean NN Distance']

    for ax, metric, title in zip(axes.flat, metrics_to_plot, titles):
        if metric not in df.columns:
            ax.text(0.5, 0.5, f'{metric}\nnot available', ha='center', va='center')
            ax.set_title(title)
            continue

        for _, row in df.iterrows():
            color = RESPONSE_COLORS[row['response']]
            marker = 'o' if row['response'] == 'R' else 's'
            ax.scatter(row['response'], row[metric], c=color, marker=marker,
                      s=150, edgecolor='black', linewidth=1, zorder=5)
            ax.annotate(row['patient'], (row['response'], row[metric]),
                       xytext=(10, 0), textcoords='offset points', fontsize=9)

        # Add means
        for resp in ['R', 'NR']:
            vals = df[df['response'] == resp][metric].dropna()
            if len(vals) > 0:
                mean_val = vals.mean()
                ax.hlines(mean_val, -0.3, 0.3 if resp == 'R' else 0.7,
                         colors=RESPONSE_COLORS[resp], linestyles='--', linewidth=2)

        ax.set_title(title)
        ax.set_ylabel(metric.replace('_', ' ').title())

    plt.suptitle('Spatial Metrics (Post-treatment, n=2 vs n=2)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig05_spatial_post.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved fig05_spatial_post.png")

    df.to_excel(TABLE_DIR / 'T4_spatial_metrics.xlsx', index=False)
    print("  Saved T4_spatial_metrics.xlsx")


def create_summary_table():
    """Create summary statistics table."""
    print("Creating Summary Tables...")

    # Collect all metrics
    all_data = []

    for patient_id, info in PATIENTS.items():
        for timepoint in ['pre', 'post']:
            sample_id = info[timepoint]
            if sample_id is None:
                continue

            result = get_patient_metrics(patient_id, timepoint)
            if result:
                metrics, _ = result
                all_data.append(metrics)

    df = pd.DataFrame(all_data)

    # Summary by response and timepoint
    summary_rows = []

    for timepoint in ['pre', 'post']:
        tp_data = df[df['timepoint'] == timepoint]

        for response in ['R', 'NR']:
            resp_data = tp_data[tp_data['response'] == response]
            n = len(resp_data)

            if n == 0:
                continue

            row = {
                'timepoint': timepoint.capitalize(),
                'response': response,
                'n_patients': n,
                'patients': ', '.join(resp_data['patient'].tolist()),
            }

            # Add mean metrics
            for col in resp_data.columns:
                if col.startswith('pct_') or col in ['n_spots']:
                    vals = resp_data[col].dropna()
                    if len(vals) > 0:
                        row[f'{col}_mean'] = vals.mean()
                        if len(vals) > 1:
                            row[f'{col}_sd'] = vals.std()

            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_excel(TABLE_DIR / 'T0_summary.xlsx', index=False)
    print("  Saved T0_summary.xlsx")

    # Also save raw data
    df.to_excel(TABLE_DIR / 'T0_all_metrics.xlsx', index=False)
    print("  Saved T0_all_metrics.xlsx")


def main():
    print("="*60)
    print("PDAC SPATIAL ANALYSIS - PROPERLY STRATIFIED")
    print("="*60)
    print()
    print("Patient Structure:")
    print("  R:  YP12 (Pre+Post), YP15 (Pre+Post)")
    print("  NR: YP03 (Pre+Post), YP04 (Post only)")
    print()
    print("NO P-VALUES - descriptive statistics only (n=2 per group)")
    print("="*60)

    # Create all figures
    create_summary_table()
    create_figure1_overview()
    create_figure2_celltype_post()
    create_figure3_celltype_pre()
    create_figure4_paired_change()
    create_figure5_spatial_post()

    print()
    print("="*60)
    print("COMPLETE")
    print("="*60)
    print(f"Figures: {FIG_DIR}")
    print(f"Tables: {TABLE_DIR}")


if __name__ == '__main__':
    main()
