# PDAC Spatial Transcriptomics - Properly Stratified Analysis

## Patient Structure

| Patient | Response | Pre-treatment | Post-treatment |
|---------|----------|---------------|----------------|
| YP12 | Responder | ✓ YP12A | ✓ YP12C |
| YP15 | Responder | ✓ YP15A | ✓ YP15C |
| YP03 | Non-Responder | ✓ YP03A | ✓ YP03C |
| YP04 | Non-Responder | ✗ Missing | ✓ YP04C |

**Total: 4 patients (2 R, 2 NR)**

## Statistical Approach

⚠️ **NO P-VALUES REPORTED** - With n=2 per group, statistical inference is not meaningful.

All analyses show:
- Individual patient values
- Means ± SD where applicable
- Effect sizes (Cohen's d) for reference only

## Valid Comparisons

1. **Post-treatment**: R (n=2) vs NR (n=2) - Best comparison
2. **Pre-treatment**: R (n=2) vs NR (n=1) - Limited (YP04 missing)
3. **Paired change**: R (2 pairs) vs NR (1 pair) - Limited

## Figures

| Figure | Description |
|--------|-------------|
| fig01 | Study overview and patient structure |
| fig02 | Cell type proportions - Post-treatment (n=2 vs n=2) |
| fig03 | Cell type proportions - Pre-treatment (n=2 vs n=1) |
| fig04 | Paired Pre→Post changes |
| fig05 | Spatial metrics - Post-treatment |

## Data Tables

| Table | Contents |
|-------|----------|
| T0_summary.xlsx | Summary statistics by group |
| T0_all_metrics.xlsx | All patient-level metrics |
| T1_celltype_post.xlsx | Cell type % (Post) |
| T2_celltype_pre.xlsx | Cell type % (Pre) |
| T3_paired_changes.xlsx | Pre→Post changes |
| T4_spatial_metrics.xlsx | Spatial metrics |

## Reproducibility

```bash
conda activate enact
python scripts/analysis_stratified.py
```
