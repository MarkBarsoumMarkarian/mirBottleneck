mirBottleneck — inst/extdata contents
======================================

TOY DATASETS (shipped with package — used for examples, tests, vignettes)
--------------------------------------------------------------------------
toy_mirna_log.rds       log2(count+1) miRNA matrix, 30 miRNAs x 10 patients (synthetic)
toy_rna_sym.rds         log2(count+1) RNA-seq matrix, 50 genes x 10 patients (synthetic)
toy_clinical_df.rds     Clinical data frame, 10 patients x 9 columns (synthetic)
toy_mirna_norm_map.rds  miRNA precursor ID normalization map (synthetic)
toy_mirna_targets.rds   Precomputed miRNA-target network, 30 miRNAs (synthetic)

All toy files are fully synthetic (set.seed(42)) and contain no real patient data.
They are designed to run the full mirBottleneck pipeline offline in < 60 seconds.

REAL TCGA-PAAD ANALYSIS OUTPUTS (reference results only)
---------------------------------------------------------
clinical.rds            Patient clinical metadata (184 x 9)
master_clinical.rds     Extended clinical table
mirna_norm_map.rds      Full miRNA normalization map (TCGA-PAAD)
gene_symbol_map.rds     ENSEMBL ID to HGNC symbol mapping
shared_patients.rds     167 patients with complete multi-omics data
bottleneck_scores.rds   Per-miRNA VSS + coherence scores
coherence_scores.rds    Coherence induction scores
combined_bottleneck.rds Classified miRNA table with bottleneck index
cox_individual.rds      Individual miRNA Cox results
cox_model_final.rds     Final composite Cox model
cox_composite_final.rds Full composite model output
survival_df.rds         Patient survival + bottleneck score (167 x 11)
survival_composite.rds  Composite survival analysis results

FULL DATASET (not in package — too large)
------------------------------------------
The full TCGA-PAAD harmonized dataset (~51 MB) is available on Zenodo:
https://doi.org/10.5281/zenodo.19121116

Download with:
  mirBottleneck::fetch_mirBottleneck_data()  # caches to tools::R_user_dir()
