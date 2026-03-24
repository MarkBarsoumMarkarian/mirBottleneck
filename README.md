# mirBottleneck

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19121116.svg)](https://doi.org/10.5281/zenodo.19121116)

> **miRNA Bottleneck Scoring for Pancreatic Cancer Survival Prediction**  
> Mark Barsoum Markarian · American University of Beirut · mkm25@mail.aub.edu

---

## Overview

`mirBottleneck` is an R package that identifies miRNAs acting as **transcriptome stabilizers** in pancreatic adenocarcinoma (TCGA-PAAD). It scores each miRNA along two complementary axes:

| Score | What it measures |
|---|---|
| **VSS** — Variance Suppression Score | How much a miRNA reduces expression variance across its validated targets |
| **Coherence Score** | How much a miRNA coordinates its targets into coherent co-expression programs |

miRNAs are then classified into four functional archetypes — **Silencer**, **Conductor**, **Dual**, and **Weak** — and a patient-level **composite bottleneck index** is built from the top survival-associated bottleneck miRNAs. This index is evaluated against overall survival using Cox proportional hazards models.

![Kaplan-Meier: High vs Low Bottleneck Score](inst/figures/km_plot.png)

> **New (v0.99.1):** Hallmark ssGSEA entropy analysis added — high bottleneck score patients show significantly elevated transcriptome disorder (Spearman rho = 0.261, p = 0.0007). See `inst/scripts/08_hallmark_entropy.R`.

---

## Repository Structure

```
mirBottleneck/
├── R/                          # Package source code
│   ├── build_network.R         # miRNA–target network from miRTarBase
│   ├── score_vss.R             # Variance Suppression Score
│   ├── score_coherence.R       # Coherence Induction Score
│   ├── classify.R              # miRNA classification + bottleneck index
│   ├── composite_score.R       # Patient-level composite score
│   ├── survival_model.R        # Cox models + log-rank test
│   └── utils.R                 # Barcode harmonization, normalization helpers
│
├── data/                       # Harmonized TCGA-PAAD analysis objects
│   ├── clinical.rds            # Patient clinical metadata (184 × 9)
│   ├── master_clinical.rds     # Extended clinical table
│   ├── rna_matrix.rds          # Raw RNA-seq counts (60,660 genes × 178 patients)
│   ├── rna_log.rds             # Log2-transformed RNA-seq matrix
│   ├── rna_final.rds           # Final filtered RNA matrix (gene symbols)
│   ├── mirna_matrix.rds        # Raw miRNA counts (1,881 × 178 patients)
│   ├── mirna_log.rds           # Log2-transformed miRNA matrix
│   ├── mirna_final.rds         # Final filtered miRNA matrix
│   ├── mirna_norm_map.rds      # Precursor ID → normalized ID mapping
│   ├── gene_symbol_map.rds     # ENSEMBL ID → HGNC symbol mapping
│   ├── maf_filtered.rds        # Filtered somatic mutations (MAF format)
│   ├── mirna_targets.rds       # miRNA–target interactions (all)
│   ├── mirna_targets_strict.rds# miRNA–target interactions (strict filter)
│   ├── mirna_target_interactions.rds # Full interaction table
│   ├── interactions_filtered.rds     # Filtered interaction network
│   ├── interactions_strict.rds       # Strict interaction network
│   ├── shared_patients.rds     # Patients with complete multi-omics data (n=167)
│   ├── bottleneck_scores.rds   # Per-miRNA VSS + coherence scores
│   ├── coherence_scores.rds    # Coherence induction scores
│   ├── combined_bottleneck.rds # Classified miRNA table with bottleneck index
│   ├── cox_individual.rds      # Individual miRNA Cox results
│   ├── cox_model_final.rds     # Final composite Cox model
│   ├── cox_composite_final.rds # Full composite model output
│   ├── survival_df.rds         # Patient survival + bottleneck score (167 × 11)
│   └── survival_composite.rds  # Composite survival analysis results
│
├── figures/
│   └── km_plot.png             # Kaplan-Meier: High vs Low bottleneck score
│
├── DESCRIPTION                 # R package metadata
├── NAMESPACE                   # Exported functions
└── LICENSE                     # MIT License
```

---

## R Package Functions

### `build_network(mirna_ids, rna_symbols, ...)`
Queries miRTarBase via `multiMiR` to retrieve validated miRNA–target interactions. Filters to targets present in the RNA-seq dataset and batches queries to avoid server timeouts.

### `score_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)`
Fits a linear model of each target gene on miRNA expression across patients. VSS = mean R² across all validated targets per miRNA.

### `score_coherence(mirna_log, rna_sym, mirna_targets, ...)`
Splits patients by median miRNA expression, then measures whether target genes become more correlated in the high-expression group. Coherence = Δ mean pairwise correlation (high − low), with empirical p-values via permutation.

### `classify_bottleneck(vss_scores, coherence_scores)`
Normalizes both scores to [0,1], computes a composite bottleneck index, and classifies each miRNA:

| Class | VSS | Coherence |
|---|---|---|
| **Dual** (strongest) | High | High |
| **Silencer** | High | Low |
| **Conductor** | Low | High |
| **Weak** | Low | Low |

### `composite_score(mirna_log, classified_df, clinical_df, ...)`
Builds a direction-aware weighted patient-level score using survival-associated bottleneck miRNAs. miRNAs with HR > 1 push the score up; protective miRNAs push it down. Weights are proportional to |Cox coefficient|.

### `survival_model(patient_scores, clinical_df)`
Fits three Cox models (score alone, clinical alone, combined) and a log-rank test comparing High vs Low score groups.

### Utilities (`utils.R`)
- `harmonize_barcode()` — standardize TCGA barcodes to 12-character patient level
- `normalize_mirna_id()` — convert mature miRNA IDs to lowercase precursor format
- `normalize_01()` — min-max normalization

---

## Data

All `.rds` objects in `data/` are derived from **TCGA-PAAD** (pancreatic adenocarcinoma) raw data downloaded via `TCGAbiolinks`. Raw data is not included due to size and access restrictions. The harmonized objects represent the analysis-ready output of the preprocessing pipeline.

**Cohort:** 178 patients (RNA-seq + miRNA) · 167 with complete multi-omics data  
**Source:** The Cancer Genome Atlas (TCGA) — https://www.cancer.gov/tcga

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("multiMiR", "biomaRt", "TCGAbiolinks", "survminer"))

# CRAN
install.packages(c("dplyr", "survival"))
```

---

## Quick Start

```r
# Load harmonized data
mirna_log  <- readRDS("data/mirna_log.rds")
rna_final  <- readRDS("data/rna_final.rds")
mirna_norm <- readRDS("data/mirna_norm_map.rds")
clinical   <- readRDS("data/clinical.rds")

# Source functions
sapply(list.files("R", full.names = TRUE), source)

# 1. Build miRNA-target network
targets <- build_network(rownames(mirna_log), rownames(rna_final))

# 2. Score miRNAs
vss    <- score_vss(mirna_log, rna_final, targets, mirna_norm)
coher  <- score_coherence(mirna_log, rna_final, targets, mirna_norm)

# 3. Classify
classified <- classify_bottleneck(vss, coher)

# 4. Patient-level composite score
composite  <- composite_score(mirna_log, classified, clinical, mirna_norm)

# 5. Survival analysis
surv <- survival_model(composite$patient_scores, clinical)
```

---


---

## Hallmark Entropy Analysis

Script: `inst/scripts/08_hallmark_entropy.R`

Tests the hypothesis that bottleneck miRNAs maintain transcriptional order by correlating the patient-level composite bottleneck score with pathway-level Shannon entropy (ssGSEA over 50 MSigDB Hallmark gene sets).

| Result | Value |
|---|---|
| Spearman rho (bottleneck score ~ entropy) | 0.261 |
| p-value | 0.0007 |
| Top disordered pathway (high-bottleneck) | HALLMARK_PROTEIN_SECRETION |
| Second | HALLMARK_DNA_REPAIR |

High bottleneck score patients have significantly higher transcriptome entropy, supporting the interpretation that the survival signal reflects loss of regulatory order rather than activation of a single oncogenic pathway.

Results are in `inst/figures/hallmark_entropy_full.csv` and `inst/figures/hallmark_pathway_disorder.csv`.


---

## Hallmark Entropy Analysis

Script: [`inst/scripts/08_hallmark_entropy.R`](inst/scripts/08_hallmark_entropy.R)

Tests the hypothesis that bottleneck miRNAs maintain transcriptional order by correlating the composite bottleneck score with pathway-level Shannon entropy (ssGSEA over 50 MSigDB Hallmark gene sets, n = 167 TCGA-PAAD patients).

| Result | Value |
|---|---|
| Spearman rho (bottleneck score ~ transcriptome entropy) | **0.261** |
| p-value | **0.0007** |
| Top disordered pathway in high-bottleneck patients | `HALLMARK_PROTEIN_SECRETION` |
| Second | `HALLMARK_DNA_REPAIR` |

High bottleneck score patients have significantly higher transcriptome entropy, supporting the interpretation that the survival signal reflects loss of regulatory order across oncogenic programmes rather than activation of a single pathway.

Result tables: [`inst/figures/hallmark_entropy_full.csv`](inst/figures/hallmark_entropy_full.csv) · [`inst/figures/hallmark_pathway_disorder.csv`](inst/figures/hallmark_pathway_disorder.csv)

## Citation

If you use this code or data, please cite:

> Markarian, M.B. (2026). *mirBottleneck: miRNA Bottleneck Scoring for Pancreatic Cancer Survival Prediction*. GitHub. https://github.com/MarkBarsoumMarkarian/mirBottleneck

---


---

## Full Dataset (Zenodo)

The complete harmonized TCGA-PAAD dataset (~51 MB) is archived on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.19121116

Download directly from R (cached locally, never required for package checks):

```r
paths <- fetch_mirBottleneck_data()
rna   <- readRDS(paths[["rna_log"]])   # 60,660 genes x 178 patients
```

## License

MIT © Mark Barsoum Markarian
