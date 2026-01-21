#!/usr/bin/env python3
"""
Deep Dive: Largest Effect Sizes in PDAC Spatial Transcriptomics
================================================================

Creates detailed multi-panel figures for metrics with largest effect sizes.
Multiple visualization approaches for each finding.

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
FIG_DIR = PROJECT_ROOT / "figures"
TABLE_DIR = PROJECT_ROOT / "tables"

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
PATIENT_MARKERS = {'YP12': 'o', 'YP15': 's', 'YP03': '^', 'YP04': 'D'}

def cohens_d(group1, group2):
    """Calculate Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    if pooled_std == 0:
        return 0
    return (np.mean(group1) - np.mean(group2)) / pooled_std

def load_data():
    """Load all data tables."""
    data = {}
    data['spatial'] = pd.read_excel(TABLE_DIR / 'T4_spatial_metrics.xlsx')
    data['celltype_post'] = pd.read_excel(TABLE_DIR / 'T1_celltype_post.xlsx')
    data['celltype_pre'] = pd.read_excel(TABLE_DIR / 'T2_celltype_pre.xlsx')
    data['pathways'] = pd.read_excel(TABLE_DIR / 'T7_progeny_pathways.xlsx')
    return data

def calculate_all_effect_sizes(data):
    """Calculate effect sizes for all metrics."""
    effects = []

    # Spatial metrics (Post-treatment only)
    spatial = data['spatial']
    for col in ['spatial_entropy', 'diversity_shannon', 'n_spots', 'mean_nn_dist']:
        r_vals = spatial[spatial['response'] == 'R'][col].values
        nr_vals = spatial[spatial['response'] == 'NR'][col].values
        d = cohens_d(r_vals, nr_vals)
        effects.append({
            'metric': col,
            'category': 'spatial',
            'R_mean': np.mean(r_vals),
            'R_std': np.std(r_vals),
            'NR_mean': np.mean(nr_vals),
            'NR_std': np.std(nr_vals),
            'cohens_d': d,
            'abs_d': abs(d),
            'direction': 'R > NR' if d > 0 else 'NR > R'
        })

    # Cell type proportions (Post-treatment)
    ct_post = data['celltype_post']
    for ct in ct_post['cell_type'].unique():
        ct_data = ct_post[ct_post['cell_type'] == ct]
        r_vals = ct_data[ct_data['response'] == 'R']['percentage'].values
        nr_vals = ct_data[ct_data['response'] == 'NR']['percentage'].values
        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            d = cohens_d(r_vals, nr_vals)
            effects.append({
                'metric': ct,
                'category': 'cell_type',
                'R_mean': np.mean(r_vals),
                'R_std': np.std(r_vals),
                'NR_mean': np.mean(nr_vals),
                'NR_std': np.std(nr_vals),
                'cohens_d': d,
                'abs_d': abs(d),
                'direction': 'R > NR' if d > 0 else 'NR > R'
            })

    # Pathway activities (Post-treatment only)
    pathways = data['pathways']
    pathways_post = pathways[pathways['timepoint'] == 'Post']
    pathway_cols = [c for c in pathways.columns if c not in ['sample', 'response', 'timepoint', 'patient']]

    for pathway in pathway_cols:
        r_vals = pathways_post[pathways_post['response'] == 'R'][pathway].values
        nr_vals = pathways_post[pathways_post['response'] == 'NR'][pathway].values
        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            d = cohens_d(r_vals, nr_vals)
            effects.append({
                'metric': pathway,
                'category': 'pathway',
                'R_mean': np.mean(r_vals),
                'R_std': np.std(r_vals),
                'NR_mean': np.mean(nr_vals),
                'NR_std': np.std(nr_vals),
                'cohens_d': d,
                'abs_d': abs(d),
                'direction': 'R > NR' if d > 0 else 'NR > R'
            })

    return pd.DataFrame(effects).sort_values('abs_d', ascending=False)

def fig_effect_size_ranking(effects_df):
    """Figure: Ranked effect sizes across all metrics."""
    print("\n  Generating effect size ranking figure...")

    fig, axes = plt.subplots(1, 3, figsize=(16, 8))

    # A) All effects ranked
    ax = axes[0]
    top_n = 20
    top_effects = effects_df.head(top_n).copy()
    colors = ['#2ecc71' if d > 0 else '#e74c3c' for d in top_effects['cohens_d']]

    y_pos = range(len(top_effects))
    ax.barh(y_pos, top_effects['cohens_d'], color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_effects['metric'])
    ax.axvline(0, color='black', linewidth=0.5)
    ax.axvline(0.8, color='gray', linestyle='--', alpha=0.5, label='Large effect (|d|=0.8)')
    ax.axvline(-0.8, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel("Cohen's d (R - NR)")
    ax.set_title('A) Top 20 Effect Sizes (Post-treatment)')
    ax.invert_yaxis()

    # Add effect size annotations
    for i, (_, row) in enumerate(top_effects.iterrows()):
        x_pos = row['cohens_d'] + 0.1 * np.sign(row['cohens_d'])
        ax.text(x_pos, i, f"{row['cohens_d']:.2f}", va='center', fontsize=8)

    # B) By category
    ax = axes[1]
    categories = ['spatial', 'cell_type', 'pathway']
    cat_colors = {'spatial': '#3498db', 'cell_type': '#9b59b6', 'pathway': '#f39c12'}

    for i, cat in enumerate(categories):
        cat_data = effects_df[effects_df['category'] == cat].head(5)
        for j, (_, row) in enumerate(cat_data.iterrows()):
            color = '#2ecc71' if row['cohens_d'] > 0 else '#e74c3c'
            ax.barh(i + j*0.15, row['abs_d'], height=0.12,
                   color=color, alpha=0.8,
                   label=row['metric'] if j < 3 else '')
            ax.text(row['abs_d'] + 0.05, i + j*0.15,
                   f"{row['metric'][:12]}", va='center', fontsize=7)

    ax.set_yticks(range(len(categories)))
    ax.set_yticklabels(['Spatial\nMetrics', 'Cell\nTypes', 'Pathways'])
    ax.set_xlabel("|Cohen's d|")
    ax.set_title('B) Top Effects by Category')
    ax.axvline(0.8, color='gray', linestyle='--', alpha=0.5)

    # C) Effect size interpretation guide
    ax = axes[2]
    ax.axis('off')

    # Top effects summary
    summary_text = "TOP EFFECT SIZES (|d| > 1.0)\n" + "="*40 + "\n\n"

    large_effects = effects_df[effects_df['abs_d'] > 1.0]
    for _, row in large_effects.iterrows():
        direction = "higher" if row['cohens_d'] > 0 else "lower"
        summary_text += f"{row['metric']}\n"
        summary_text += f"  R: {row['R_mean']:.3f} ± {row['R_std']:.3f}\n"
        summary_text += f"  NR: {row['NR_mean']:.3f} ± {row['NR_std']:.3f}\n"
        summary_text += f"  d = {row['cohens_d']:.2f} (R {direction})\n\n"

    summary_text += "\n" + "="*40 + "\n"
    summary_text += "Effect size interpretation:\n"
    summary_text += "  |d| < 0.2: Negligible\n"
    summary_text += "  |d| 0.2-0.5: Small\n"
    summary_text += "  |d| 0.5-0.8: Medium\n"
    summary_text += "  |d| > 0.8: Large\n"
    summary_text += "\n⚠️ n=2 per group - descriptive only"

    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
           fontsize=9, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax.set_title('C) Summary of Largest Effects')

    plt.suptitle('Effect Size Analysis: R vs NR (Post-treatment, n=2 vs n=2)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig11_effect_size_ranking.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig11_effect_size_ranking.png")

def fig_spatial_entropy_deep_dive(data):
    """Deep dive into spatial entropy - the largest effect."""
    print("\n  Generating spatial entropy deep dive...")

    spatial = data['spatial']

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)

    # A) Strip plot with individual points
    ax = fig.add_subplot(gs[0, 0])
    for patient in spatial['patient'].unique():
        row = spatial[spatial['patient'] == patient].iloc[0]
        color = RESPONSE_COLORS[row['response']]
        marker = PATIENT_MARKERS[patient]
        x = 0 if row['response'] == 'R' else 1
        ax.scatter(x, row['spatial_entropy'], c=color, marker=marker,
                  s=200, edgecolors='black', linewidth=2, label=patient, zorder=5)

    # Add means
    r_mean = spatial[spatial['response'] == 'R']['spatial_entropy'].mean()
    nr_mean = spatial[spatial['response'] == 'NR']['spatial_entropy'].mean()
    ax.hlines(r_mean, -0.3, 0.3, colors=RESPONSE_COLORS['R'], linestyles='--', linewidth=2)
    ax.hlines(nr_mean, 0.7, 1.3, colors=RESPONSE_COLORS['NR'], linestyles='--', linewidth=2)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Responders\n(n=2)', 'Non-Responders\n(n=2)'])
    ax.set_ylabel('Spatial Entropy')
    ax.set_title('A) Individual Patients')
    ax.legend(loc='lower right', fontsize=8)

    # B) Box + swarm
    ax = fig.add_subplot(gs[0, 1])
    box_data = [spatial[spatial['response'] == 'R']['spatial_entropy'].values,
                spatial[spatial['response'] == 'NR']['spatial_entropy'].values]
    bp = ax.boxplot(box_data, positions=[0, 1], widths=0.5, patch_artist=True)
    bp['boxes'][0].set_facecolor(RESPONSE_COLORS['R'])
    bp['boxes'][1].set_facecolor(RESPONSE_COLORS['NR'])
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_alpha(0.5)

    # Overlay points
    for i, response in enumerate(['R', 'NR']):
        vals = spatial[spatial['response'] == response]['spatial_entropy'].values
        ax.scatter([i]*len(vals), vals, c=RESPONSE_COLORS[response],
                  s=100, edgecolors='black', zorder=5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['R', 'NR'])
    ax.set_ylabel('Spatial Entropy')
    ax.set_title('B) Box Plot')

    # C) Effect size visualization
    ax = fig.add_subplot(gs[0, 2])
    r_vals = spatial[spatial['response'] == 'R']['spatial_entropy'].values
    nr_vals = spatial[spatial['response'] == 'NR']['spatial_entropy'].values
    d = cohens_d(r_vals, nr_vals)

    # Visualize distributions
    x_range = np.linspace(3.3, 3.9, 100)
    r_dist = stats.norm.pdf(x_range, np.mean(r_vals), np.std(r_vals) + 0.01)
    nr_dist = stats.norm.pdf(x_range, np.mean(nr_vals), np.std(nr_vals) + 0.01)

    ax.fill_between(x_range, r_dist, alpha=0.5, color=RESPONSE_COLORS['R'], label='R')
    ax.fill_between(x_range, nr_dist, alpha=0.5, color=RESPONSE_COLORS['NR'], label='NR')
    ax.axvline(np.mean(r_vals), color=RESPONSE_COLORS['R'], linestyle='--', linewidth=2)
    ax.axvline(np.mean(nr_vals), color=RESPONSE_COLORS['NR'], linestyle='--', linewidth=2)

    # Effect size arrow
    ax.annotate('', xy=(np.mean(r_vals), max(r_dist)*0.5),
               xytext=(np.mean(nr_vals), max(r_dist)*0.5),
               arrowprops=dict(arrowstyle='<->', color='black', lw=2))
    ax.text((np.mean(r_vals) + np.mean(nr_vals))/2, max(r_dist)*0.6,
           f'd = {d:.2f}', ha='center', fontsize=12, fontweight='bold')

    ax.set_xlabel('Spatial Entropy')
    ax.set_ylabel('Density')
    ax.set_title('C) Effect Size Visualization')
    ax.legend()

    # D) Paired comparison (if we had pre-post)
    ax = fig.add_subplot(gs[0, 3])
    # Use all samples for context
    all_samples = pd.concat([
        pd.DataFrame({'sample': ['YP12C', 'YP15C', 'YP03C', 'YP04C'],
                     'entropy': [3.681, 3.720, 3.508, 3.424],
                     'response': ['R', 'R', 'NR', 'NR'],
                     'timepoint': 'Post'})
    ])

    # Bar chart
    x = np.arange(4)
    colors = [RESPONSE_COLORS['R'], RESPONSE_COLORS['R'],
              RESPONSE_COLORS['NR'], RESPONSE_COLORS['NR']]
    bars = ax.bar(x, all_samples['entropy'], color=colors, edgecolor='black', linewidth=2)
    ax.set_xticks(x)
    ax.set_xticklabels(['YP12', 'YP15', 'YP03', 'YP04'])
    ax.set_ylabel('Spatial Entropy')
    ax.set_title('D) Individual Samples')
    ax.set_ylim(3.3, 3.85)

    # Add group lines
    ax.hlines(r_mean, -0.5, 1.5, colors=RESPONSE_COLORS['R'], linestyles='--', linewidth=2, label='R mean')
    ax.hlines(nr_mean, 1.5, 3.5, colors=RESPONSE_COLORS['NR'], linestyles='--', linewidth=2, label='R mean')

    # E) Correlation with other metrics
    ax = fig.add_subplot(gs[1, 0:2])
    metrics = ['diversity_shannon', 'n_spots', 'mean_nn_dist']
    metric_labels = ['Shannon Diversity', 'N Spots', 'Mean NN Distance']

    for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
        for _, row in spatial.iterrows():
            color = RESPONSE_COLORS[row['response']]
            marker = PATIENT_MARKERS[row['patient']]
            ax.scatter(row['spatial_entropy'], row[metric],
                      c=color, marker=marker, s=100, alpha=0.8)

    # Add correlation for entropy vs diversity
    corr = spatial['spatial_entropy'].corr(spatial['diversity_shannon'])
    ax.set_xlabel('Spatial Entropy')
    ax.set_ylabel('Shannon Diversity')
    ax.set_title(f'E) Entropy vs Diversity (r = {corr:.3f})')

    # Add legend for patients
    for patient, marker in PATIENT_MARKERS.items():
        ax.scatter([], [], marker=marker, c='gray', s=100, label=patient)
    ax.legend(loc='lower right')

    # F) Interpretation panel
    ax = fig.add_subplot(gs[1, 2:4])
    ax.axis('off')

    interpretation = """
    SPATIAL ENTROPY: Deep Dive Analysis
    ════════════════════════════════════════════════════════════

    WHAT IS SPATIAL ENTROPY?
    • Shannon entropy computed on cell type proportions
    • Higher values = more diverse cellular composition
    • Lower values = dominated by fewer cell types

    KEY FINDING:
    ┌─────────────────────────────────────────────────────────┐
    │  Responders: 3.70 ± 0.02  (more diverse)               │
    │  Non-Responders: 3.47 ± 0.04  (less diverse)           │
    │                                                         │
    │  Cohen's d = 5.06 (VERY LARGE effect)                  │
    │  Difference = 0.23 entropy units                        │
    └─────────────────────────────────────────────────────────┘

    BIOLOGICAL INTERPRETATION:
    • Responder tumors have more balanced cell type composition
    • Non-responder tumors may be dominated by specific populations
    • Could indicate different tumor microenvironment states
    • May reflect immune infiltration patterns

    CAVEATS:
    ⚠️  n=2 per group - hypothesis generating only
    ⚠️  Effect size inflated by small sample
    ⚠️  Requires validation in larger cohort
    """

    ax.text(0.02, 0.98, interpretation, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # G) Entropy by cell type contribution
    ax = fig.add_subplot(gs[2, :])

    # Load cell type data to show composition
    ct_post = data['celltype_post']

    # Pivot for stacked bar
    ct_pivot = ct_post.pivot(index='patient', columns='cell_type', values='percentage')
    ct_pivot = ct_pivot.reindex(['YP12', 'YP15', 'YP03', 'YP04'])

    # Sort columns by mean abundance
    col_order = ct_pivot.mean().sort_values(ascending=False).index
    ct_pivot = ct_pivot[col_order]

    # Stacked bar
    ct_pivot.plot(kind='bar', stacked=True, ax=ax, width=0.8,
                 colormap='tab20', edgecolor='white', linewidth=0.5)

    ax.set_ylabel('Cell Type Proportion (%)')
    ax.set_xlabel('')
    ax.set_xticklabels(['YP12 (R)', 'YP15 (R)', 'YP03 (NR)', 'YP04 (NR)'], rotation=0)
    ax.set_title('G) Cell Type Composition Underlying Entropy Differences')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)

    # Add entropy values on top
    for i, patient in enumerate(['YP12', 'YP15', 'YP03', 'YP04']):
        ent = spatial[spatial['patient'] == patient]['spatial_entropy'].values[0]
        ax.text(i, 102, f'H={ent:.2f}', ha='center', fontsize=10, fontweight='bold')

    plt.suptitle('Figure 12: Spatial Entropy Deep Dive - Largest Effect Size (d=5.06)',
                fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig12_spatial_entropy_deep_dive.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig12_spatial_entropy_deep_dive.png")

def fig_cell_type_deep_dive(data, effects_df):
    """Deep dive into cell types with largest effects."""
    print("\n  Generating cell type deep dive...")

    ct_effects = effects_df[effects_df['category'] == 'cell_type'].head(6)
    ct_post = data['celltype_post']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, (_, effect_row) in enumerate(ct_effects.iterrows()):
        ax = axes[idx]
        ct = effect_row['metric']

        ct_data = ct_post[ct_post['cell_type'] == ct]

        # Individual patient points
        for _, row in ct_data.iterrows():
            color = RESPONSE_COLORS[row['response']]
            marker = PATIENT_MARKERS.get(row['patient'], 'o')
            x = 0 if row['response'] == 'R' else 1
            x += np.random.uniform(-0.1, 0.1)  # jitter
            ax.scatter(x, row['percentage'], c=color, marker=marker,
                      s=150, edgecolors='black', linewidth=2, zorder=5)

        # Means and error bars
        r_vals = ct_data[ct_data['response'] == 'R']['percentage']
        nr_vals = ct_data[ct_data['response'] == 'NR']['percentage']

        ax.errorbar(0, r_vals.mean(), yerr=r_vals.std(),
                   fmt='_', color=RESPONSE_COLORS['R'], markersize=20,
                   capsize=10, capthick=2, linewidth=2)
        ax.errorbar(1, nr_vals.mean(), yerr=nr_vals.std(),
                   fmt='_', color=RESPONSE_COLORS['NR'], markersize=20,
                   capsize=10, capthick=2, linewidth=2)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['R (n=2)', 'NR (n=2)'])
        ax.set_ylabel('Proportion (%)')

        direction = "↑" if effect_row['cohens_d'] > 0 else "↓"
        ax.set_title(f'{ct}\nd = {effect_row["cohens_d"]:.2f} ({direction} in R)')

        # Add values
        ax.text(0.02, 0.98, f"R: {r_vals.mean():.1f}%\nNR: {nr_vals.mean():.1f}%",
               transform=ax.transAxes, fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle('Figure 13: Cell Type Proportions - Top 6 Effect Sizes',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig13_celltype_deep_dive.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig13_celltype_deep_dive.png")

def fig_pathway_deep_dive(data, effects_df):
    """Deep dive into pathways with largest effects."""
    print("\n  Generating pathway deep dive...")

    pathway_effects = effects_df[effects_df['category'] == 'pathway'].head(6)
    pathways = data['pathways']
    pathways_post = pathways[pathways['timepoint'] == 'Post']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, (_, effect_row) in enumerate(pathway_effects.iterrows()):
        ax = axes[idx]
        pathway = effect_row['metric']

        # Individual patient points
        for _, row in pathways_post.iterrows():
            color = RESPONSE_COLORS[row['response']]
            marker = PATIENT_MARKERS.get(row['patient'], 'o')
            x = 0 if row['response'] == 'R' else 1
            x += np.random.uniform(-0.1, 0.1)
            ax.scatter(x, row[pathway], c=color, marker=marker,
                      s=150, edgecolors='black', linewidth=2, zorder=5)

        # Means
        r_vals = pathways_post[pathways_post['response'] == 'R'][pathway]
        nr_vals = pathways_post[pathways_post['response'] == 'NR'][pathway]

        ax.errorbar(0, r_vals.mean(), yerr=r_vals.std() if len(r_vals) > 1 else 0,
                   fmt='_', color=RESPONSE_COLORS['R'], markersize=20,
                   capsize=10, capthick=2, linewidth=2)
        ax.errorbar(1, nr_vals.mean(), yerr=nr_vals.std() if len(nr_vals) > 1 else 0,
                   fmt='_', color=RESPONSE_COLORS['NR'], markersize=20,
                   capsize=10, capthick=2, linewidth=2)

        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['R (n=2)', 'NR (n=2)'])
        ax.set_ylabel('Activity Score')

        direction = "↑" if effect_row['cohens_d'] > 0 else "↓"
        ax.set_title(f'{pathway}\nd = {effect_row["cohens_d"]:.2f} ({direction} in R)')

    plt.suptitle('Figure 14: Pathway Activities - Top 6 Effect Sizes (PROGENy)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig14_pathway_deep_dive.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig14_pathway_deep_dive.png")

def fig_summary_dashboard(data, effects_df):
    """Create a summary dashboard of top findings."""
    print("\n  Generating summary dashboard...")

    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.3)

    spatial = data['spatial']
    ct_post = data['celltype_post']
    pathways = data['pathways']
    pathways_post = pathways[pathways['timepoint'] == 'Post']

    # Row 1: Top 3 overall effects
    top3 = effects_df.head(3)

    for idx, (_, row) in enumerate(top3.iterrows()):
        ax = fig.add_subplot(gs[0, idx])
        metric = row['metric']

        if row['category'] == 'spatial':
            r_vals = spatial[spatial['response'] == 'R'][metric].values
            nr_vals = spatial[spatial['response'] == 'NR'][metric].values
        elif row['category'] == 'cell_type':
            ct_data = ct_post[ct_post['cell_type'] == metric]
            r_vals = ct_data[ct_data['response'] == 'R']['percentage'].values
            nr_vals = ct_data[ct_data['response'] == 'NR']['percentage'].values
        else:  # pathway
            r_vals = pathways_post[pathways_post['response'] == 'R'][metric].values
            nr_vals = pathways_post[pathways_post['response'] == 'NR'][metric].values

        # Dot plot
        ax.scatter([0]*len(r_vals), r_vals, c=RESPONSE_COLORS['R'], s=200,
                  edgecolors='black', linewidth=2, zorder=5)
        ax.scatter([1]*len(nr_vals), nr_vals, c=RESPONSE_COLORS['NR'], s=200,
                  edgecolors='black', linewidth=2, zorder=5)

        # Connect with line showing difference
        ax.plot([0, 1], [np.mean(r_vals), np.mean(nr_vals)], 'k--', linewidth=2, alpha=0.5)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['R', 'NR'])
        ax.set_title(f'#{idx+1}: {metric}\nd = {row["cohens_d"]:.2f}', fontweight='bold')

        # Add means
        ax.text(0, np.mean(r_vals), f'  {np.mean(r_vals):.2f}', va='center', fontsize=9)
        ax.text(1, np.mean(nr_vals), f'  {np.mean(nr_vals):.2f}', va='center', fontsize=9)

    # Info panel
    ax = fig.add_subplot(gs[0, 3])
    ax.axis('off')
    info = """
    STUDY DESIGN
    ════════════════════

    Patients: 4 total
    • 2 Responders (R)
    • 2 Non-Responders (NR)

    Samples: 7 total
    • Post-treatment: n=4
    • Pre-treatment: n=3
      (YP04 pre missing)

    ⚠️ NO P-VALUES
    n=2 per group is
    insufficient for
    statistical inference

    Effect sizes shown
    for reference only
    """
    ax.text(0.1, 0.95, info, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Row 2: Heatmap of all effect sizes
    ax = fig.add_subplot(gs[1, :2])

    # Create effect size matrix by category
    pivot_data = []
    for cat in ['spatial', 'cell_type', 'pathway']:
        cat_effects = effects_df[effects_df['category'] == cat].head(5)
        for _, row in cat_effects.iterrows():
            pivot_data.append({'metric': row['metric'][:15], 'effect': row['cohens_d']})

    effect_matrix = pd.DataFrame(pivot_data)

    colors = ['#2ecc71' if d > 0 else '#e74c3c' for d in effect_matrix['effect']]
    bars = ax.barh(range(len(effect_matrix)), effect_matrix['effect'], color=colors, alpha=0.8)
    ax.set_yticks(range(len(effect_matrix)))
    ax.set_yticklabels(effect_matrix['metric'])
    ax.axvline(0, color='black', linewidth=1)
    ax.axvline(0.8, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(-0.8, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel("Cohen's d")
    ax.set_title('Top 5 Effects per Category')
    ax.invert_yaxis()

    # Row 2: Correlation matrix
    ax = fig.add_subplot(gs[1, 2:])

    # Merge spatial and pathway data for correlation
    merged = spatial.merge(pathways_post[['patient', 'JAK-STAT', 'Trail', 'TGFb']], on='patient')
    corr_cols = ['spatial_entropy', 'diversity_shannon', 'JAK-STAT', 'Trail', 'TGFb']
    corr_matrix = merged[corr_cols].corr()

    sns.heatmap(corr_matrix, annot=True, cmap='RdBu_r', center=0,
               vmin=-1, vmax=1, ax=ax, fmt='.2f')
    ax.set_title('Correlation Matrix: Top Metrics')

    # Row 3: Pre vs Post for responders
    ax = fig.add_subplot(gs[2, :2])

    # Pathway changes
    pathway_list = ['JAK-STAT', 'Trail', 'TGFb', 'Hypoxia', 'TNFa']
    x = np.arange(len(pathway_list))
    width = 0.35

    r_pre = pathways[(pathways['response'] == 'R') & (pathways['timepoint'] == 'Pre')][pathway_list].mean()
    r_post = pathways[(pathways['response'] == 'R') & (pathways['timepoint'] == 'Post')][pathway_list].mean()

    ax.bar(x - width/2, r_pre, width, label='Pre', color='lightblue', edgecolor='black')
    ax.bar(x + width/2, r_post, width, label='Post', color='darkblue', edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(pathway_list, rotation=45, ha='right')
    ax.set_ylabel('Activity')
    ax.axhline(0, color='gray', linestyle='--')
    ax.legend()
    ax.set_title('Responders: Pre vs Post Pathway Activity')

    # Key findings
    ax = fig.add_subplot(gs[2, 2:])
    ax.axis('off')

    findings = """
    KEY FINDINGS (Descriptive Only)
    ════════════════════════════════════════════════════════════════════

    1. SPATIAL ENTROPY (d = 5.06)
       • Responders: 3.70 ± 0.02 (higher cellular diversity)
       • Non-responders: 3.47 ± 0.04
       • Suggests more heterogeneous TME in responders

    2. CELL TYPE DIFFERENCES
       • Acinar cells: Higher in R (8.1% vs 2.3%, d = 1.63)
       • Ductal epithelial: Higher in NR (17.8% vs 8.0%, d = -2.14)
       • NK cells: Higher in R (7.5% vs 3.7%, d = 1.52)

    3. PATHWAY ACTIVITIES
       • JAK-STAT: Lower in R (d = -1.25) - reduced inflammation?
       • TGFb: Variable, higher in NR trend
       • Trail: Post-treatment increase in R

    ════════════════════════════════════════════════════════════════════
    INTERPRETATION:
    Responders show more diverse cell composition and different
    pathway activation patterns. These are HYPOTHESIS-GENERATING
    findings requiring validation in larger cohorts.
    """

    ax.text(0.02, 0.98, findings, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.2))

    plt.suptitle('Figure 15: Summary Dashboard - PDAC Spatial Transcriptomics (n=2 vs n=2)',
                fontsize=14, fontweight='bold', y=1.01)
    plt.savefig(FIG_DIR / 'fig15_summary_dashboard.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig15_summary_dashboard.png")

def main():
    print("="*60)
    print("DEEP DIVE: Largest Effect Sizes")
    print("="*60)

    print("\nLoading data...")
    data = load_data()

    print("\nCalculating effect sizes...")
    effects_df = calculate_all_effect_sizes(data)

    # Save effect sizes
    effects_df.to_excel(TABLE_DIR / 'T9_effect_sizes.xlsx', index=False)
    print(f"  Saved T9_effect_sizes.xlsx ({len(effects_df)} metrics)")

    # Print top effects
    print("\n  TOP 10 EFFECT SIZES (|d|):")
    print("  " + "-"*50)
    for _, row in effects_df.head(10).iterrows():
        print(f"  {row['metric']:<20} d={row['cohens_d']:>7.2f} ({row['category']})")

    # Generate figures
    print("\nGenerating figures...")
    fig_effect_size_ranking(effects_df)
    fig_spatial_entropy_deep_dive(data)
    fig_cell_type_deep_dive(data, effects_df)
    fig_pathway_deep_dive(data, effects_df)
    fig_summary_dashboard(data, effects_df)

    print("\n" + "="*60)
    print("COMPLETE")
    print("="*60)
    print("\nGenerated figures:")
    for f in sorted(FIG_DIR.glob('fig1[1-5]*.png')):
        print(f"  {f.name}")
    print(f"\nEffect sizes table: T9_effect_sizes.xlsx")

if __name__ == '__main__':
    main()
