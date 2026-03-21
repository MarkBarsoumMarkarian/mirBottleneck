# mirBottleneck

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19121116.svg)](https://doi.org/10.5281/zenodo.19121116)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-submission-brightgreen)](https://github.com/Bioconductor/Contributions/issues/4203)
[![R CMD Check](https://github.com/MarkBarsoumMarkarian/mirBottleneck/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MarkBarsoumMarkarian/mirBottleneck/actions)

**miRNA Bottleneck Scoring for Pancreatic Cancer Survival Prediction**

Mark Barsoum Markarian · American University of Beirut · mkm25@mail.aub.edu

---

## What it does

Most miRNA biomarker tools ask which miRNAs correlate with a clinical outcome. mirBottleneck asks a different question: which miRNAs actually impose order on the transcriptome?

The package scores each miRNA along two independent axes:

| Score | What it captures |
|---|---|
| **VSS** (Variance Suppression Score) | How much a miRNA reduces expression variance across its validated targets |
| **Coherence Score** | How much a miRNA coordinates its targets into a coherent co-expression program |

miRNAs are classified into four archetypes based on where they fall in VSS-coherence space:

- **Dual** -- high VSS + high coherence. True transcriptome stabilizers.
- **Silencer** -- high VSS, low coherence. Broad noise suppressors (hsa-miR-21, hsa-miR-155).
- **Conductor** -- low VSS, high coherence. Program coordinators (hsa-miR-217, hsa-miR-802).
- **Weak** -- low on both axes.

A patient-level composite bottleneck index is built from the miRNAs most associated with survival and evaluated against overall survival using Cox proportional hazards models.

![Kaplan-Meier: High vs Low Bottleneck Score](inst/figures/km_plot.png)

Applied to 178 TCGA-PAAD patients: **HR = 6.55 (95% CI 3.42-12.55), p = 1.43 x 10^-8, C-index = 0.699.**

---

## One command, one report

```r
results <- run_mirBottleneck_project(
  mirna_log_rds      = "mirna_log.rds",
  rna_rds            = "rna_sym.rds",
  clinical_rds       = "clinical.rds",
  mirna_norm_map_rds = "norm_map.rds",
  out_dir            = "results/",
  report             = TRUE
)
# Saves results/mirBottleneck_report.html
# Self-contained, fully offline, opens in any browser
```

The HTML report includes the phase diagram, survival curves, sortable miRNA table, per-patient scores, and Cox model results.

---

## Installation

```r
# From Bioconductor (once accepted)
BiocManager::install("mirBottleneck")

# Development version from GitHub
remotes::install_github("MarkBarsoumMarkarian/mirBottleneck")
```

---

## Quick start with toy data

The package ships a small synthetic dataset (10 patients, 50 genes, 30 miRNAs) so you can run the full pipeline offline in under 60 seconds:

```r
library(mirBottleneck)

paths <- mirBottleneck:::.toy_paths()

results <- run_mirBottleneck_project(
  mirna_log_rds      = paths$mirna_log,
  rna_rds            = paths$rna_sym,
  clinical_rds       = paths$clinical_df,
  mirna_norm_map_rds = paths$mirna_norm_map,
  mirna_targets_rds  = paths$mirna_targets,  # precomputed, no internet needed
  query_network      = FALSE,
  out_dir            = tempdir(),
  cox_p_threshold    = 0.5,
  report             = TRUE
)
```

---

## Step-by-step usage

```r
library(mirBottleneck)

# 1. Build miRNA-target network (requires internet, queries miRTarBase)
targets <- build_network(
  mirna_ids   = rownames(mirna_log),
  rna_symbols = rownames(rna_sym)
)

# 2. Score miRNAs
vss    <- score_vss(mirna_log, rna_sym, targets, mirna_norm_map)
coher  <- score_coherence(mirna_log, rna_sym, targets, mirna_norm_map,
                          scored_mirnas = vss$mirna)

# 3. Classify into archetypes
classified <- classify_bottleneck(vss, coher)

# 4. Build patient-level composite score
composite <- composite_score(mirna_log, classified, clinical_df, mirna_norm_map)

# 5. Survival analysis
surv <- survival_model(composite$patient_scores, clinical_df)

# 6. Generate HTML report
generate_report(
  results     = c(composite, surv, list(combined_bottleneck = classified,
                  vss_scores = vss, coherence_scores = coher)),
  out_dir     = "results/",
  cohort_name = "My Cohort"
)
```

---

## Package structure

```
mirBottleneck/
├── R/
│   ├── build_network.R        # miRNA-target network from miRTarBase
│   ├── score_vss.R            # Variance Suppression Score
│   ├── score_coherence.R      # Coherence Induction Score
│   ├── classify.R             # miRNA classification + bottleneck index
│   ├── composite_score.R      # Patient-level composite score
│   ├── survival_model.R       # Cox models + log-rank test
│   ├── generate_report.R      # HTML report generation
│   ├── fetch_data.R           # Download full dataset from Zenodo
│   ├── toy_data.R             # Toy dataset paths
│   └── utils.R                # Barcode harmonization, normalization helpers
│
├── inst/
│   ├── extdata/               # Toy datasets + reference results
│   │   ├── toy_mirna_log.rds
│   │   ├── toy_rna_sym.rds
│   │   ├── toy_clinical_df.rds
│   │   ├── toy_mirna_norm_map.rds
│   │   ├── toy_mirna_targets.rds
│   │   └── README.txt
│   ├── report/
│   │   └── mirBottleneck_report.Rmd   # HTML report template
│   └── figures/
│       ├── km_plot.png
│       └── vss_coherence_classification.png
│
├── vignettes/
│   └── mirBottleneck-workflow.Rmd
│
├── tests/testthat/
│   ├── test-classify-bottleneck.R
│   └── test-mirna-expr-helpers.R
│
├── DESCRIPTION
├── NAMESPACE
└── LICENSE
```

---

## Full TCGA-PAAD dataset

The complete harmonized dataset (~51 MB, 178 patients, 60,660 genes, 1,881 miRNAs) is archived on Zenodo. It is not bundled in the package but can be downloaded on demand:

**DOI:** https://doi.org/10.5281/zenodo.19121116

```r
# Downloads and caches to tools::R_user_dir("mirBottleneck", "cache")
paths <- fetch_mirBottleneck_data()

rna  <- readRDS(paths[["rna_log"]])    # 60,660 x 178
mirna <- readRDS(paths[["mirna_log"]]) # 1,881 x 178
```

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("multiMiR", "biomaRt"))

# CRAN
install.packages(c("dplyr", "survival", "rmarkdown", "plotly", "DT"))
```

---

## Citation

If you use this package or the harmonized dataset, please cite:

> Markarian MB. mirBottleneck: a dual-score framework for identifying transcriptome-stabilizing miRNAs and predicting survival in pancreatic adenocarcinoma. Manuscript in preparation. 2026.
> https://github.com/MarkBarsoumMarkarian/mirBottleneck

> Markarian MB. Harmonized multi-omics TCGA-PAAD dataset. Zenodo. 2026.
> https://doi.org/10.5281/zenodo.19121116

---

## License

MIT (c) Mark Barsoum Markarian
