#!/usr/bin/env python3
"""
Advanced Stratified Analysis - PDAC Spatial Transcriptomics
============================================================

Comprehensive analyses properly stratified by patient:
- Spatial Entropy (cell type diversity)
- Topology/Betti Curves (persistent homology)
- PROGENy Pathway Activity
- GSEA Hallmark Enrichment
- H&E + Spatial Transcriptome Visualization

Author: Max Van Belkum
Date: 2026-01-20
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import seaborn as sns
import scanpy as sc
import squidpy as sq
from scipy.stats import entropy
from scipy.spatial.distance import pdist, squareform
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Optional imports
try:
    import decoupler as dc
    HAS_DECOUPLER = True
except ImportError:
    HAS_DECOUPLER = False
    print("Warning: decoupler not available")

try:
    from gtda.homology import VietorisRipsPersistence
    from gtda.diagrams import BettiCurve
    HAS_GIOTTO = True
except ImportError:
    HAS_GIOTTO = False
    print("Warning: giotto-tda not available - topology analysis will be skipped")

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
OLD_REPO = Path("/home/user/spatial-hackathon-2026")
DATA_ROOT = Path("/home/user/data/hackathon/PDAC.Visium")
ADATA_DIR = OLD_REPO / "outputs" / "adata" / "polymathic"
FIG_DIR = PROJECT_ROOT / "figures"
TABLE_DIR = PROJECT_ROOT / "tables"

FIG_DIR.mkdir(exist_ok=True)
TABLE_DIR.mkdir(exist_ok=True)

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Color scheme
RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
TIMEPOINT_COLORS = {'Pre': '#3498db', 'Post': '#9b59b6'}

# Patient structure - THE GROUND TRUTH
PATIENTS = {
    'YP12': {'response': 'R', 'pre': 'YP12A', 'post': 'YP12C'},
    'YP15': {'response': 'R', 'pre': 'YP15A', 'post': 'YP15C'},
    'YP03': {'response': 'NR', 'pre': 'YP03A', 'post': 'YP03C'},
    'YP04': {'response': 'NR', 'pre': None, 'post': 'YP04C'},  # MISSING PRE!
}

SAMPLES = {
    "YP03A": {"response": "NR", "timepoint": "Pre", "patient": "YP03"},
    "YP03C": {"response": "NR", "timepoint": "Post", "patient": "YP03"},
    "YP04C": {"response": "NR", "timepoint": "Post", "patient": "YP04"},
    "YP12A": {"response": "R", "timepoint": "Pre", "patient": "YP12"},
    "YP12C": {"response": "R", "timepoint": "Post", "patient": "YP12"},
    "YP15A": {"response": "R", "timepoint": "Pre", "patient": "YP15"},
    "YP15C": {"response": "R", "timepoint": "Post", "patient": "YP15"},
}


def load_samples():
    """Load all polymathic h5ad files."""
    adatas = {}
    for sample in SAMPLES.keys():
        path = ADATA_DIR / f"{sample}_polymathic.h5ad"
        if path.exists():
            adata = sc.read_h5ad(path)
            adata.obs['sample'] = sample
            adata.obs['response'] = SAMPLES[sample]['response']
            adata.obs['timepoint'] = SAMPLES[sample]['timepoint']
            adata.obs['patient'] = SAMPLES[sample]['patient']
            adatas[sample] = adata
            print(f"  Loaded {sample}: {adata.n_obs} spots")
    return adatas


def get_hires_image(sample):
    """Get H&E hires image path for a sample."""
    img_path = DATA_ROOT / "spaceranger.out" / sample / "spatial" / "tissue_hires_image.png"
    if img_path.exists():
        return plt.imread(img_path)
    return None


def compute_spatial_entropy(adata, cell_type_col='cell_type'):
    """
    Compute spatial entropy based on cell type distribution.

    Higher entropy = more diverse cell type composition
    """
    if cell_type_col not in adata.obs.columns:
        # Try alternatives
        for col in ['cell_type_broad', 'leiden']:
            if col in adata.obs.columns:
                cell_type_col = col
                break
        else:
            return None

    ct_counts = adata.obs[cell_type_col].value_counts(normalize=True)
    # Shannon entropy (higher = more diverse)
    return entropy(ct_counts.values, base=2)


def compute_betti_features(coords, max_points=2000, seed=42):
    """
    Compute persistent homology features for tissue architecture.

    Returns Betti curves and summary statistics.
    """
    if not HAS_GIOTTO:
        return None

    np.random.seed(seed)

    # Subsample if needed (MaxMin for topology preservation)
    if len(coords) > max_points:
        # Simple random subsample for now
        idx = np.random.choice(len(coords), max_points, replace=False)
        coords = coords[idx]

    # Standardize coordinates
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)

    try:
        # Compute persistence
        persistence = VietorisRipsPersistence(
            homology_dimensions=[0, 1],
            max_edge_length=2.0,
            n_jobs=-1
        )
        diagrams = persistence.fit_transform(coords.reshape(1, -1, 2))

        # Compute Betti curves
        betti = BettiCurve(n_bins=100)
        curves = betti.fit_transform(diagrams)

        betti0_curve = curves[0, :, 0]  # Connected components
        betti1_curve = curves[0, :, 1]  # Holes/loops

        return {
            'betti0_curve': betti0_curve,
            'betti1_curve': betti1_curve,
            'betti0_auc': np.trapz(betti0_curve),
            'betti1_auc': np.trapz(betti1_curve),
            'betti0_max': betti0_curve.max(),
            'betti1_max': betti1_curve.max(),
            'topology_complexity': np.trapz(betti0_curve) + np.trapz(betti1_curve)
        }
    except Exception as e:
        print(f"    Topology error: {e}")
        return None


def run_progeny_analysis(adatas):
    """Run PROGENy pathway activity analysis with proper dtype conversion."""
    if not HAS_DECOUPLER:
        print("  decoupler not available, skipping PROGENy")
        return None

    print("\n  Loading PROGENy model...")
    progeny = dc.get_progeny(organism='human', top=500)
    # FIX: Convert pandas nullable dtypes to numpy
    progeny['weight'] = progeny['weight'].astype('float32')
    progeny['p_value'] = progeny['p_value'].astype('float32')

    pathway_results = {}

    for sample, adata in adatas.items():
        print(f"  Processing {sample}...")
        try:
            dc.run_mlm(
                mat=adata,
                net=progeny,
                source='source',
                target='target',
                weight='weight',
                verbose=False,
                use_raw=False
            )
            adata.obsm['pathway_activities'] = adata.obsm.pop('mlm_estimate')

            # Store mean activities
            mean_activities = adata.obsm['pathway_activities'].mean(axis=0)
            pathway_results[sample] = mean_activities.to_dict()
            print(f"    {len(mean_activities)} pathways computed")

        except Exception as e:
            print(f"    PROGENy failed: {e}")

    return pathway_results


def run_gsea_hallmark(adatas):
    """Run GSEA on Hallmark gene sets."""
    print("\n  Running GSEA Hallmark analysis...")

    try:
        # Load Hallmark gene sets
        msigdb = dc.get_resource('MSigDB')
        hallmarks = msigdb[msigdb['collection'] == 'hallmark']
        # FIX: Remove duplicate edges
        hallmarks = hallmarks.drop_duplicates(subset=['geneset', 'genesymbol'])
        print(f"    Loaded {len(hallmarks['geneset'].unique())} Hallmark gene sets")
    except Exception as e:
        print(f"    Failed to load MSigDB: {e}")
        return None

    gsea_results = {}

    # For each sample, compute enrichment scores
    for sample, adata in adatas.items():
        print(f"  Processing {sample}...")
        try:
            dc.run_ora(
                mat=adata,
                net=hallmarks,
                source='geneset',
                target='genesymbol',
                verbose=False,
                use_raw=False
            )

            if 'ora_estimate' in adata.obsm:
                mean_scores = adata.obsm['ora_estimate'].mean(axis=0)
                gsea_results[sample] = mean_scores.to_dict()
                print(f"    {len(mean_scores)} gene sets scored")
        except Exception as e:
            print(f"    GSEA failed: {e}")

    return gsea_results


def generate_figure_06_spatial_entropy(adatas, metrics_df):
    """
    Figure 6: Spatial Entropy Analysis

    Shows cell type diversity across samples, stratified by response and timepoint.
    """
    print("\n  Generating Figure 6: Spatial Entropy...")

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # A) Entropy by sample (bar chart with patient grouping)
    ax = fig.add_subplot(gs[0, 0])

    samples_ordered = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']
    samples_present = [s for s in samples_ordered if s in metrics_df['sample'].values]

    df_plot = metrics_df[metrics_df['sample'].isin(samples_present)].set_index('sample').loc[samples_present]

    colors = [RESPONSE_COLORS[SAMPLES[s]['response']] for s in samples_present]
    hatches = ['/' if SAMPLES[s]['timepoint'] == 'Pre' else '' for s in samples_present]

    bars = ax.bar(range(len(samples_present)), df_plot['spatial_entropy'], color=colors, edgecolor='black')
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)

    ax.set_xticks(range(len(samples_present)))
    ax.set_xticklabels(samples_present, rotation=45, ha='right')
    ax.set_ylabel('Spatial Entropy (bits)')
    ax.set_title('A) Spatial Entropy by Sample')

    # Add patient brackets
    ax.axvline(1.5, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(3.5, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(5.5, color='gray', linestyle='--', alpha=0.5)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=RESPONSE_COLORS['R'], label='Responder'),
        Patch(facecolor=RESPONSE_COLORS['NR'], label='Non-Responder'),
        Patch(facecolor='white', edgecolor='black', hatch='/', label='Pre-treatment'),
        Patch(facecolor='white', edgecolor='black', label='Post-treatment'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    # B) R vs NR comparison (Post-treatment only - n=2 vs n=2)
    ax = fig.add_subplot(gs[0, 1])

    post_r = metrics_df[(metrics_df['response'] == 'R') & (metrics_df['timepoint'] == 'Post')]['spatial_entropy']
    post_nr = metrics_df[(metrics_df['response'] == 'NR') & (metrics_df['timepoint'] == 'Post')]['spatial_entropy']

    # Individual points
    ax.scatter([0]*len(post_r), post_r, c=RESPONSE_COLORS['R'], s=100, zorder=5, label='R')
    ax.scatter([1]*len(post_nr), post_nr, c=RESPONSE_COLORS['NR'], s=100, zorder=5, label='NR')

    # Means
    ax.bar([0, 1], [post_r.mean(), post_nr.mean()],
           color=[RESPONSE_COLORS['R'], RESPONSE_COLORS['NR']], alpha=0.3, width=0.6)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Responders\n(n=2)', 'Non-Responders\n(n=2)'])
    ax.set_ylabel('Spatial Entropy (bits)')
    ax.set_title('B) Post-Treatment: R vs NR')

    # Effect size
    if len(post_r) > 1 and len(post_nr) > 1:
        pooled_std = np.sqrt((post_r.var() + post_nr.var()) / 2)
        if pooled_std > 0:
            d = (post_r.mean() - post_nr.mean()) / pooled_std
            ax.text(0.5, 0.95, f"Cohen's d = {d:.2f}", transform=ax.transAxes,
                   ha='center', fontsize=10, style='italic')

    ax.text(0.5, 0.02, 'NO P-VALUE (n=2)', transform=ax.transAxes,
           ha='center', fontsize=9, color='red', weight='bold')

    # C) Pre vs Post paired changes
    ax = fig.add_subplot(gs[0, 2])

    paired_data = []
    for patient, info in PATIENTS.items():
        if info['pre'] and info['post']:
            pre_val = metrics_df[metrics_df['sample'] == info['pre']]['spatial_entropy'].values
            post_val = metrics_df[metrics_df['sample'] == info['post']]['spatial_entropy'].values
            if len(pre_val) > 0 and len(post_val) > 0:
                paired_data.append({
                    'patient': patient,
                    'response': info['response'],
                    'pre': pre_val[0],
                    'post': post_val[0],
                    'change': post_val[0] - pre_val[0]
                })

    if paired_data:
        pdf = pd.DataFrame(paired_data)

        for _, row in pdf.iterrows():
            color = RESPONSE_COLORS[row['response']]
            ax.plot([0, 1], [row['pre'], row['post']], 'o-', color=color,
                   markersize=10, linewidth=2, label=f"{row['patient']} ({row['response']})")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Pre-treatment', 'Post-treatment'])
        ax.set_ylabel('Spatial Entropy (bits)')
        ax.set_title('C) Paired Pre→Post Changes')
        ax.legend(loc='best', fontsize=8)

    # D-F) Spatial plots showing cell type distributions for select samples
    for idx, sample in enumerate(['YP12C', 'YP15C', 'YP03C']):
        ax = fig.add_subplot(gs[1, idx])

        if sample in adatas:
            adata = adatas[sample]

            # Get cell type column
            ct_col = None
            for col in ['cell_type', 'cell_type_broad', 'leiden']:
                if col in adata.obs.columns:
                    ct_col = col
                    break

            if ct_col and 'spatial' in adata.obsm:
                coords = adata.obsm['spatial']
                categories = adata.obs[ct_col].astype('category')

                scatter = ax.scatter(coords[:, 0], coords[:, 1],
                                    c=categories.cat.codes,
                                    cmap='tab20', s=5, alpha=0.7)

                response = SAMPLES[sample]['response']
                entropy_val = metrics_df[metrics_df['sample'] == sample]['spatial_entropy'].values[0]

                ax.set_title(f"D{idx+1}) {sample} ({response})\nEntropy = {entropy_val:.3f}",
                            color=RESPONSE_COLORS[response])
                ax.set_xlabel('Spatial X')
                ax.set_ylabel('Spatial Y')
                ax.set_aspect('equal')

    plt.suptitle('Figure 6: Spatial Entropy Analysis\n(Cell Type Diversity)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig06_spatial_entropy.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig06_spatial_entropy.png")


def generate_figure_07_topology(adatas, metrics_df):
    """
    Figure 7: Topological Data Analysis (Persistent Homology)

    Betti curves and summary statistics.
    """
    if not HAS_GIOTTO:
        print("\n  Skipping Figure 7: giotto-tda not available")
        return

    print("\n  Generating Figure 7: Topology Analysis...")

    # Compute Betti features for all samples
    topology_data = {}
    for sample, adata in adatas.items():
        print(f"    Computing topology for {sample}...")
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            result = compute_betti_features(coords)
            if result:
                topology_data[sample] = result

    if not topology_data:
        print("    No topology data computed")
        return

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # A) Betti-0 curves (connected components)
    ax = fig.add_subplot(gs[0, 0])

    for sample, data in topology_data.items():
        color = RESPONSE_COLORS[SAMPLES[sample]['response']]
        linestyle = '--' if SAMPLES[sample]['timepoint'] == 'Pre' else '-'
        ax.plot(data['betti0_curve'], color=color, linestyle=linestyle,
               label=f"{sample}", alpha=0.8)

    ax.set_xlabel('Filtration Parameter')
    ax.set_ylabel('Betti-0 (Connected Components)')
    ax.set_title('A) Betti-0 Curves')
    ax.legend(fontsize=7, ncol=2)

    # B) Betti-1 curves (holes/loops)
    ax = fig.add_subplot(gs[0, 1])

    for sample, data in topology_data.items():
        color = RESPONSE_COLORS[SAMPLES[sample]['response']]
        linestyle = '--' if SAMPLES[sample]['timepoint'] == 'Pre' else '-'
        ax.plot(data['betti1_curve'], color=color, linestyle=linestyle,
               label=f"{sample}", alpha=0.8)

    ax.set_xlabel('Filtration Parameter')
    ax.set_ylabel('Betti-1 (Loops/Holes)')
    ax.set_title('B) Betti-1 Curves')
    ax.legend(fontsize=7, ncol=2)

    # C) Betti-1 AUC comparison (Post-treatment)
    ax = fig.add_subplot(gs[0, 2])

    post_samples = [s for s in topology_data.keys() if SAMPLES[s]['timepoint'] == 'Post']

    r_auc = [topology_data[s]['betti1_auc'] for s in post_samples if SAMPLES[s]['response'] == 'R']
    nr_auc = [topology_data[s]['betti1_auc'] for s in post_samples if SAMPLES[s]['response'] == 'NR']

    ax.scatter([0]*len(r_auc), r_auc, c=RESPONSE_COLORS['R'], s=100, zorder=5)
    ax.scatter([1]*len(nr_auc), nr_auc, c=RESPONSE_COLORS['NR'], s=100, zorder=5)
    ax.bar([0, 1], [np.mean(r_auc) if r_auc else 0, np.mean(nr_auc) if nr_auc else 0],
          color=[RESPONSE_COLORS['R'], RESPONSE_COLORS['NR']], alpha=0.3, width=0.6)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['R (n=2)', 'NR (n=2)'])
    ax.set_ylabel('Betti-1 AUC')
    ax.set_title('C) Post-Treatment Betti-1 AUC')
    ax.text(0.5, 0.02, 'NO P-VALUE (n=2)', transform=ax.transAxes,
           ha='center', fontsize=9, color='red', weight='bold')

    # D) Topology complexity by sample
    ax = fig.add_subplot(gs[1, 0])

    samples_ordered = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']
    samples_present = [s for s in samples_ordered if s in topology_data]

    complexity = [topology_data[s]['topology_complexity'] for s in samples_present]
    colors = [RESPONSE_COLORS[SAMPLES[s]['response']] for s in samples_present]
    hatches = ['/' if SAMPLES[s]['timepoint'] == 'Pre' else '' for s in samples_present]

    bars = ax.bar(range(len(samples_present)), complexity, color=colors, edgecolor='black')
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)

    ax.set_xticks(range(len(samples_present)))
    ax.set_xticklabels(samples_present, rotation=45, ha='right')
    ax.set_ylabel('Topology Complexity')
    ax.set_title('D) Topology Complexity by Sample')

    # E) Pre→Post paired changes in topology
    ax = fig.add_subplot(gs[1, 1])

    paired_topo = []
    for patient, info in PATIENTS.items():
        if info['pre'] and info['post']:
            if info['pre'] in topology_data and info['post'] in topology_data:
                paired_topo.append({
                    'patient': patient,
                    'response': info['response'],
                    'pre': topology_data[info['pre']]['betti1_auc'],
                    'post': topology_data[info['post']]['betti1_auc'],
                })

    if paired_topo:
        for row in paired_topo:
            color = RESPONSE_COLORS[row['response']]
            ax.plot([0, 1], [row['pre'], row['post']], 'o-', color=color,
                   markersize=10, linewidth=2, label=f"{row['patient']} ({row['response']})")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Pre-treatment', 'Post-treatment'])
        ax.set_ylabel('Betti-1 AUC')
        ax.set_title('E) Paired Pre→Post Betti-1 Changes')
        ax.legend(loc='best', fontsize=8)

    # F) Interpretation text box
    ax = fig.add_subplot(gs[1, 2])
    ax.axis('off')

    textstr = '''TOPOLOGICAL DATA ANALYSIS
========================

Method: Persistent Homology
(Vietoris-Rips filtration)

Betti Numbers:
• Betti-0: Connected components
  Higher = more fragmented tissue

• Betti-1: Holes/loops in structure
  Higher = more complex architecture

Interpretation:
• Tumor tissue often shows altered
  topological features vs normal
• Treatment may reorganize tissue
  architecture (change in loops)

NOTE: With n=2 per group, these
are DESCRIPTIVE observations only.
No statistical inference possible.'''

    ax.text(0.1, 0.9, textstr, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('Figure 7: Topological Data Analysis (Persistent Homology)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig07_topology_betti.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig07_topology_betti.png")

    # Return topology data for Excel export
    return topology_data


def generate_figure_08_progeny(pathway_results):
    """
    Figure 8: PROGENy Pathway Activity Analysis
    """
    if not pathway_results:
        print("\n  Skipping Figure 8: No pathway results")
        return

    print("\n  Generating Figure 8: PROGENy Pathway Activity...")

    # Convert to DataFrame
    pathway_df = pd.DataFrame(pathway_results).T
    pathway_df['sample'] = pathway_df.index
    pathway_df['response'] = pathway_df['sample'].map(lambda x: SAMPLES[x]['response'])
    pathway_df['timepoint'] = pathway_df['sample'].map(lambda x: SAMPLES[x]['timepoint'])
    pathway_df['patient'] = pathway_df['sample'].map(lambda x: SAMPLES[x]['patient'])

    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # Get pathway columns
    pathway_cols = [c for c in pathway_df.columns if c not in ['sample', 'response', 'timepoint', 'patient']]

    # A) Heatmap of pathway activities
    ax = fig.add_subplot(gs[0, 0])

    samples_ordered = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']
    samples_present = [s for s in samples_ordered if s in pathway_df.index]

    heatmap_data = pathway_df.loc[samples_present, pathway_cols].T

    sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, ax=ax,
               cbar_kws={'label': 'Activity'})
    ax.set_title('A) Pathway Activities by Sample')

    # Add response color annotation
    for i, sample in enumerate(samples_present):
        color = RESPONSE_COLORS[SAMPLES[sample]['response']]
        ax.add_patch(Rectangle((i, -0.5), 1, 0.5, color=color, clip_on=False))

    # B) R vs NR comparison (Post-treatment)
    ax = fig.add_subplot(gs[0, 1])

    post_r = pathway_df[(pathway_df['response'] == 'R') & (pathway_df['timepoint'] == 'Post')]
    post_nr = pathway_df[(pathway_df['response'] == 'NR') & (pathway_df['timepoint'] == 'Post')]

    r_means = post_r[pathway_cols].mean()
    nr_means = post_nr[pathway_cols].mean()

    x = np.arange(len(pathway_cols))
    width = 0.35

    ax.barh(x - width/2, r_means.values, width, label='R (n=2)', color=RESPONSE_COLORS['R'])
    ax.barh(x + width/2, nr_means.values, width, label='NR (n=2)', color=RESPONSE_COLORS['NR'])
    ax.set_yticks(x)
    ax.set_yticklabels(pathway_cols, fontsize=8)
    ax.set_xlabel('Mean Activity')
    ax.set_title('B) Post-Treatment: R vs NR')
    ax.axvline(0, color='black', linewidth=0.5)
    ax.legend()

    # C) Pre→Post changes (Responders)
    ax = fig.add_subplot(gs[1, 0])

    pre_r = pathway_df[(pathway_df['response'] == 'R') & (pathway_df['timepoint'] == 'Pre')]
    post_r = pathway_df[(pathway_df['response'] == 'R') & (pathway_df['timepoint'] == 'Post')]

    if len(pre_r) > 0 and len(post_r) > 0:
        r_change = post_r[pathway_cols].mean() - pre_r[pathway_cols].mean()
        r_change = r_change.sort_values()

        colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in r_change.values]
        ax.barh(range(len(r_change)), r_change.values, color=colors)
        ax.set_yticks(range(len(r_change)))
        ax.set_yticklabels(r_change.index, fontsize=8)
        ax.set_xlabel('Activity Change (Post - Pre)')
        ax.set_title('C) Responders: Pre→Post Pathway Changes')
        ax.axvline(0, color='black', linewidth=0.5)

    # D) Effect sizes
    ax = fig.add_subplot(gs[1, 1])

    effects = []
    for pathway in pathway_cols:
        r_vals = post_r[pathway].values
        nr_vals = post_nr[pathway].values

        if len(r_vals) > 1 and len(nr_vals) > 1:
            pooled_std = np.sqrt((np.var(r_vals) + np.var(nr_vals)) / 2)
            if pooled_std > 0:
                d = (np.mean(r_vals) - np.mean(nr_vals)) / pooled_std
            else:
                d = 0
            effects.append({'pathway': pathway, 'effect_size': d})

    if effects:
        eff_df = pd.DataFrame(effects).sort_values('effect_size')
        colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in eff_df['effect_size'].values]

        ax.barh(range(len(eff_df)), eff_df['effect_size'].values, color=colors)
        ax.set_yticks(range(len(eff_df)))
        ax.set_yticklabels(eff_df['pathway'].values, fontsize=8)
        ax.set_xlabel("Cohen's d (R - NR)")
        ax.set_title("D) Effect Sizes (Post-Treatment)")
        ax.axvline(0, color='black', linewidth=0.5)

        ax.text(0.5, 0.02, 'NO P-VALUES (n=2)', transform=ax.transAxes,
               ha='center', fontsize=9, color='red', weight='bold')

    plt.suptitle('Figure 8: PROGENy Pathway Activity Analysis', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig08_progeny_pathways.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig08_progeny_pathways.png")

    return pathway_df


def generate_figure_09_gsea(gsea_results):
    """
    Figure 9: GSEA Hallmark Enrichment
    """
    if not gsea_results:
        print("\n  Skipping Figure 9: No GSEA results")
        return

    print("\n  Generating Figure 9: GSEA Hallmark...")

    # Convert to DataFrame
    gsea_df = pd.DataFrame(gsea_results).T
    gsea_df['sample'] = gsea_df.index
    gsea_df['response'] = gsea_df['sample'].map(lambda x: SAMPLES[x]['response'])
    gsea_df['timepoint'] = gsea_df['sample'].map(lambda x: SAMPLES[x]['timepoint'])

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

    # Get gene set columns
    geneset_cols = [c for c in gsea_df.columns if c not in ['sample', 'response', 'timepoint']]

    # A) Heatmap
    ax = fig.add_subplot(gs[0, 0])

    samples_ordered = ['YP12A', 'YP12C', 'YP15A', 'YP15C', 'YP03A', 'YP03C', 'YP04C']
    samples_present = [s for s in samples_ordered if s in gsea_df.index]

    # Select top variable gene sets
    variances = gsea_df[geneset_cols].var()
    top_genesets = variances.nlargest(20).index.tolist()

    heatmap_data = gsea_df.loc[samples_present, top_genesets].T

    sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, ax=ax,
               cbar_kws={'label': 'Enrichment Score'})
    ax.set_title('A) Top 20 Variable Hallmark Gene Sets')
    ax.set_xticklabels(samples_present, rotation=45, ha='right')

    # B) R vs NR comparison (Post-treatment)
    ax = fig.add_subplot(gs[0, 1])

    post_r = gsea_df[(gsea_df['response'] == 'R') & (gsea_df['timepoint'] == 'Post')]
    post_nr = gsea_df[(gsea_df['response'] == 'NR') & (gsea_df['timepoint'] == 'Post')]

    if len(post_r) > 0 and len(post_nr) > 0:
        r_means = post_r[top_genesets].mean()
        nr_means = post_nr[top_genesets].mean()
        diff = r_means - nr_means
        diff = diff.sort_values()

        colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in diff.values]

        # Clean up names for display
        display_names = [n.replace('HALLMARK_', '').replace('_', ' ')[:25] for n in diff.index]

        ax.barh(range(len(diff)), diff.values, color=colors)
        ax.set_yticks(range(len(diff)))
        ax.set_yticklabels(display_names, fontsize=7)
        ax.set_xlabel('Enrichment Difference (R - NR)')
        ax.set_title('B) Post-Treatment R vs NR Difference')
        ax.axvline(0, color='black', linewidth=0.5)

        ax.text(0.5, 0.02, 'NO P-VALUES (n=2)', transform=ax.transAxes,
               ha='center', fontsize=9, color='red', weight='bold')

    plt.suptitle('Figure 9: GSEA Hallmark Gene Set Enrichment', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'fig09_gsea_hallmark.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig09_gsea_hallmark.png")

    return gsea_df


def generate_figure_10_he_spatial(adatas):
    """
    Figure 10: H&E + Spatial Transcriptome Visualization

    Per-patient panel showing:
    - H&E image
    - Spatial cell types
    - Gene expression (marker)
    """
    print("\n  Generating Figure 10: H&E + Spatial Transcriptome...")

    # Create figure with 4 patients x 4 columns
    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.3, wspace=0.2)

    patients_order = ['YP12', 'YP15', 'YP03', 'YP04']

    for row_idx, patient in enumerate(patients_order):
        info = PATIENTS[patient]

        # Get post-treatment sample (always available)
        post_sample = info['post']

        if post_sample not in adatas:
            continue

        adata = adatas[post_sample]
        response = info['response']

        # A) H&E image
        ax = fig.add_subplot(gs[row_idx, 0])
        he_img = get_hires_image(post_sample)
        if he_img is not None:
            ax.imshow(he_img)
            ax.set_title(f'{patient} H&E\n({response})', color=RESPONSE_COLORS[response], fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'H&E not\navailable', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{patient} ({response})', color=RESPONSE_COLORS[response])
        ax.axis('off')

        # B) Spatial cell types
        ax = fig.add_subplot(gs[row_idx, 1])

        ct_col = None
        for col in ['cell_type', 'cell_type_broad', 'leiden']:
            if col in adata.obs.columns:
                ct_col = col
                break

        if ct_col and 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            categories = adata.obs[ct_col].astype('category')

            scatter = ax.scatter(coords[:, 0], coords[:, 1],
                               c=categories.cat.codes,
                               cmap='tab20', s=3, alpha=0.8)
            ax.set_title(f'{post_sample} Cell Types')
            ax.set_aspect('equal')
        ax.axis('off')

        # C) Gene expression - KRT19 (ductal/tumor marker)
        ax = fig.add_subplot(gs[row_idx, 2])

        if 'spatial' in adata.obsm and 'KRT19' in adata.var_names:
            coords = adata.obsm['spatial']
            gene_expr = adata[:, 'KRT19'].X
            if hasattr(gene_expr, 'toarray'):
                gene_expr = gene_expr.toarray().flatten()
            else:
                gene_expr = np.array(gene_expr).flatten()

            scatter = ax.scatter(coords[:, 0], coords[:, 1],
                               c=gene_expr, cmap='Reds', s=3, alpha=0.8,
                               vmin=0)
            ax.set_title(f'KRT19 (Ductal)')
            ax.set_aspect('equal')
            plt.colorbar(scatter, ax=ax, shrink=0.5)
        ax.axis('off')

        # D) Gene expression - PTPRC (immune marker)
        ax = fig.add_subplot(gs[row_idx, 3])

        if 'spatial' in adata.obsm and 'PTPRC' in adata.var_names:
            coords = adata.obsm['spatial']
            gene_expr = adata[:, 'PTPRC'].X
            if hasattr(gene_expr, 'toarray'):
                gene_expr = gene_expr.toarray().flatten()
            else:
                gene_expr = np.array(gene_expr).flatten()

            scatter = ax.scatter(coords[:, 0], coords[:, 1],
                               c=gene_expr, cmap='Blues', s=3, alpha=0.8,
                               vmin=0)
            ax.set_title(f'PTPRC (Immune)')
            ax.set_aspect('equal')
            plt.colorbar(scatter, ax=ax, shrink=0.5)
        ax.axis('off')

    # Column labels
    fig.text(0.15, 0.98, 'H&E Histology', ha='center', fontsize=12, fontweight='bold')
    fig.text(0.38, 0.98, 'Cell Types', ha='center', fontsize=12, fontweight='bold')
    fig.text(0.61, 0.98, 'KRT19 Expression', ha='center', fontsize=12, fontweight='bold')
    fig.text(0.84, 0.98, 'PTPRC Expression', ha='center', fontsize=12, fontweight='bold')

    plt.suptitle('Figure 10: H&E and Spatial Transcriptome by Patient (Post-Treatment)', fontsize=14, y=1.01)
    plt.savefig(FIG_DIR / 'fig10_he_spatial.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("    Saved fig10_he_spatial.png")


def generate_excel_tables(metrics_df, topology_data, pathway_df, gsea_df):
    """Generate comprehensive Excel tables for all analyses."""
    print("\n  Generating Excel Tables...")

    # T5: Spatial Entropy
    if metrics_df is not None and 'spatial_entropy' in metrics_df.columns:
        entropy_df = metrics_df[['sample', 'patient', 'response', 'timepoint', 'spatial_entropy']].copy()
        entropy_df.to_excel(TABLE_DIR / 'T5_spatial_entropy.xlsx', index=False)
        print("    Saved T5_spatial_entropy.xlsx")

    # T6: Topology Metrics
    if topology_data:
        topo_rows = []
        for sample, data in topology_data.items():
            topo_rows.append({
                'sample': sample,
                'patient': SAMPLES[sample]['patient'],
                'response': SAMPLES[sample]['response'],
                'timepoint': SAMPLES[sample]['timepoint'],
                'betti0_auc': data['betti0_auc'],
                'betti1_auc': data['betti1_auc'],
                'betti0_max': data['betti0_max'],
                'betti1_max': data['betti1_max'],
                'topology_complexity': data['topology_complexity']
            })
        topo_df = pd.DataFrame(topo_rows)
        topo_df.to_excel(TABLE_DIR / 'T6_topology_metrics.xlsx', index=False)
        print("    Saved T6_topology_metrics.xlsx")

    # T7: Pathway Activities
    if pathway_df is not None:
        pathway_df.to_excel(TABLE_DIR / 'T7_progeny_pathways.xlsx', index=False)
        print("    Saved T7_progeny_pathways.xlsx")

    # T8: GSEA Hallmark
    if gsea_df is not None:
        gsea_df.to_excel(TABLE_DIR / 'T8_gsea_hallmark.xlsx', index=False)
        print("    Saved T8_gsea_hallmark.xlsx")


def main():
    print("="*70)
    print("ADVANCED STRATIFIED ANALYSES - PDAC SPATIAL TRANSCRIPTOMICS")
    print("="*70)
    print("\nPatient Structure:")
    print("  YP12 (R): Pre=YP12A, Post=YP12C")
    print("  YP15 (R): Pre=YP15A, Post=YP15C")
    print("  YP03 (NR): Pre=YP03A, Post=YP03C")
    print("  YP04 (NR): Pre=MISSING, Post=YP04C")
    print("\n⚠️  NO P-VALUES - n=2 per group insufficient for inference")

    print("\n" + "="*70)
    print("LOADING DATA")
    print("="*70)
    adatas = load_samples()

    if not adatas:
        print("ERROR: No samples loaded!")
        return

    # Compute basic metrics for all samples
    print("\n" + "="*70)
    print("COMPUTING METRICS")
    print("="*70)

    metrics_rows = []
    for sample, adata in adatas.items():
        print(f"  Processing {sample}...")

        row = {
            'sample': sample,
            'patient': SAMPLES[sample]['patient'],
            'response': SAMPLES[sample]['response'],
            'timepoint': SAMPLES[sample]['timepoint'],
            'n_spots': adata.n_obs,
        }

        # Spatial entropy
        entropy_val = compute_spatial_entropy(adata)
        if entropy_val is not None:
            row['spatial_entropy'] = entropy_val
            print(f"    Spatial entropy: {entropy_val:.3f}")

        metrics_rows.append(row)

    metrics_df = pd.DataFrame(metrics_rows)

    # Run advanced analyses
    print("\n" + "="*70)
    print("PATHWAY ANALYSIS (PROGENy)")
    print("="*70)
    pathway_results = run_progeny_analysis(adatas)

    print("\n" + "="*70)
    print("GSEA HALLMARK")
    print("="*70)
    gsea_results = run_gsea_hallmark(adatas)

    # Generate figures
    print("\n" + "="*70)
    print("GENERATING FIGURES")
    print("="*70)

    generate_figure_06_spatial_entropy(adatas, metrics_df)
    topology_data = generate_figure_07_topology(adatas, metrics_df)
    pathway_df = generate_figure_08_progeny(pathway_results)
    gsea_df = generate_figure_09_gsea(gsea_results)
    generate_figure_10_he_spatial(adatas)

    # Generate Excel tables
    generate_excel_tables(metrics_df, topology_data, pathway_df, gsea_df)

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)

    # List outputs
    print("\nGenerated Figures:")
    for f in sorted(FIG_DIR.glob('fig*.png')):
        print(f"  {f.name}")

    print("\nGenerated Tables:")
    for f in sorted(TABLE_DIR.glob('T*.xlsx')):
        print(f"  {f.name}")


if __name__ == '__main__':
    main()
