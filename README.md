# PDAC Spatial Transcriptomics - Properly Stratified Analysis

## Patient Structure

| Patient | Response | Pre-treatment | Post-treatment |
|---------|----------|---------------|----------------|
| YP12 | Responder | YP12A | YP12C |
| YP15 | Responder | YP15A | YP15C |
| YP03 | Non-Responder | YP03A | YP03C |
| YP04 | Non-Responder | **MISSING** | YP04C |

**Total: 4 patients (2 R, 2 NR)**

## Statistical Approach

**NO P-VALUES REPORTED** - With n=2 per group, statistical inference is not meaningful.

All analyses show:
- Individual patient values
- Means +/- SD where applicable
- Effect sizes (Cohen's d) for reference only

## Valid Comparisons

1. **Post-treatment**: R (n=2) vs NR (n=2) - Best comparison
2. **Pre-treatment**: R (n=2) vs NR (n=1) - Limited (YP04 missing)
3. **Paired change**: R (2 pairs) vs NR (1 pair) - Limited

---

## Figures

### Basic Analyses (fig01-05)

| Figure | Description |
|--------|-------------|
| fig01 | Study overview and patient structure |
| fig02 | Cell type proportions - Post-treatment (n=2 vs n=2) |
| fig03 | Cell type proportions - Pre-treatment (n=2 vs n=1) |
| fig04 | Paired Pre->Post changes |
| fig05 | Spatial metrics - Post-treatment |

### Advanced Analyses (fig06-10)

| Figure | Description | Methods |
|--------|-------------|---------|
| fig06 | Spatial Entropy (cell type diversity) | Shannon entropy on cell type distributions |
| fig07 | Topology Analysis (Betti curves) | Persistent homology via Vietoris-Rips filtration |
| fig08 | PROGENy Pathway Activity | 14 cancer-relevant pathways (decoupler) |
| fig09 | GSEA Hallmark Gene Sets | 50 Hallmark gene sets (MSigDB) |
| fig10 | H&E + Spatial Transcriptome | Per-patient visualization with KRT19/PTPRC markers |

### Deep Dive: Largest Effects (fig11-15)

| Figure | Description | Key Finding |
|--------|-------------|-------------|
| fig11 | Effect Size Ranking | Top 20 effects across all metrics |
| fig12 | Spatial Entropy Deep Dive | d=5.06, multiple visualizations |
| fig13 | Cell Type Deep Dive | Top 6 cell type differences |
| fig14 | Pathway Deep Dive | Top 6 pathway differences (PROGENy) |
| fig15 | Summary Dashboard | Integrated view of all findings |

---

## Key Findings (Descriptive Only)

### Top Effect Sizes (|d| > 1.0)

| Rank | Metric | Category | Cohen's d | Direction |
|------|--------|----------|-----------|-----------|
| 1 | Endothelial | Cell type | -5.35 | Higher in NR |
| 2 | Diversity/Entropy | Spatial | +5.06 | Higher in R |
| 3 | JAK-STAT | Pathway | -3.34 | Higher in NR |
| 4 | Unknown cells | Cell type | +3.17 | Higher in R |
| 5 | Low_Confidence | Cell type | -2.68 | Higher in NR |
| 6 | VEGF | Pathway | -2.65 | Higher in NR |
| 7 | Ductal_Epithelial | Cell type | -2.53 | Higher in NR |
| 8 | NK_cells | Cell type | +2.17 | Higher in R |
| 9 | TGFb | Pathway | -1.75 | Higher in NR |
| 10 | Acinar | Cell type | +1.58 | Higher in R |

### Spatial Entropy (Fig 6, 12)
- **Responders show higher spatial entropy** (more diverse cell type composition)
- R: 3.70 ± 0.02 vs NR: 3.47 ± 0.04 (Cohen's d = 5.06)
- Underlying composition: R have more balanced cell type distribution

### Cell Type Differences (Fig 13)
- **Ductal Epithelial**: Higher in NR (17.8% vs 8.0%, d = -2.53)
- **NK cells**: Higher in R (7.5% vs 3.7%, d = +2.17)
- **Acinar cells**: Higher in R (8.1% vs 2.3%, d = +1.58)
- Suggests immune infiltration (NK) in responders, epithelial dominance in NR

### Pathway Activity (Fig 8, 14)
- **JAK-STAT**: Higher in NR (d = -3.34) - active inflammation
- **VEGF**: Higher in NR (d = -2.65) - angiogenesis
- **TGFb**: Higher in NR (d = -1.75) - EMT/fibrosis
- **Hypoxia**: Higher in NR (d = -1.31)
- Non-responders show more aggressive pathway signatures

### Correlation Structure (Fig 15)
- Spatial entropy and JAK-STAT: r = -0.94 (inverse relationship)
- Higher diversity associated with lower inflammatory signaling

---

## Data Tables

### Basic Tables (T0-T4)

| Table | Contents |
|-------|----------|
| T0_summary.xlsx | Summary statistics by group |
| T0_all_metrics.xlsx | All patient-level metrics |
| T1_celltype_post.xlsx | Cell type % (Post) |
| T2_celltype_pre.xlsx | Cell type % (Pre) |
| T3_paired_changes.xlsx | Pre->Post changes |
| T4_spatial_metrics.xlsx | Spatial metrics |

### Advanced Tables (T5-T9)

| Table | Contents |
|-------|----------|
| T5_spatial_entropy.xlsx | Spatial entropy per sample |
| T6_topology_metrics.xlsx | Betti numbers and complexity |
| T7_progeny_pathways.xlsx | PROGENy pathway activities |
| T8_gsea_hallmark.xlsx | Hallmark gene set enrichment |
| T9_effect_sizes.xlsx | All effect sizes ranked (32 metrics) |

---

## Methods

### Spatial Entropy
Shannon entropy computed on cell type proportions. Higher values indicate more diverse cellular composition.

### Topology (Persistent Homology)
Vietoris-Rips filtration on spatial coordinates using giotto-tda.
- Betti-0: Connected components (higher = more fragmented)
- Betti-1: Loops/holes (higher = more complex architecture)

### PROGENy Pathway Activity
decoupler's MLM method on 14 cancer-relevant pathways:
Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, PI3K, TGFb, TNFa, Trail, VEGF, WNT, p53

### GSEA Hallmark
Over-representation analysis on MSigDB Hallmark collection (50 gene sets).

---

## Reproducibility

```bash
conda activate enact
python scripts/analysis_stratified.py      # Basic analyses (fig01-05)
python scripts/advanced_stratified.py      # Advanced analyses (fig06-10)
python scripts/deep_dive_effects.py        # Effect size deep dive (fig11-15)
```

---

## Acknowledgments

Data: PDAC Visium spatial transcriptomics cohort
Analysis: Spatial Biology Hackathon 2026
