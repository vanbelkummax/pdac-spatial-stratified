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

---

## Key Findings (Descriptive Only)

### Spatial Entropy (Fig 6)
- **Responders show higher spatial entropy** (more diverse cell type composition)
- R mean: 3.70 vs NR mean: 3.47 (Cohen's d = 5.06)
- Pre->Post: Responders increase entropy; NR (YP03) decreases

### Topology (Fig 7)
- Betti-1 curves (loops/holes in tissue architecture) show patient-specific patterns
- Treatment appears to reorganize tissue architecture differently between R vs NR

### Pathway Activity (Fig 8)
- **JAK-STAT**: Lower in R vs NR (effect d = -5.0)
- **PI3K, TRAIL**: Higher in R vs NR (effect d > 1.0)
- Responders show Trail pathway increase post-treatment

### GSEA Hallmark (Fig 9)
- **EMT, TNFa signaling, Hypoxia**: Higher in NR vs R
- Suggests more aggressive tumor phenotype in non-responders

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

### Advanced Tables (T5-T8)

| Table | Contents |
|-------|----------|
| T5_spatial_entropy.xlsx | Spatial entropy per sample |
| T6_topology_metrics.xlsx | Betti numbers and complexity |
| T7_progeny_pathways.xlsx | PROGENy pathway activities |
| T8_gsea_hallmark.xlsx | Hallmark gene set enrichment |

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
```

---

## Acknowledgments

Data: PDAC Visium spatial transcriptomics cohort
Analysis: Spatial Biology Hackathon 2026
