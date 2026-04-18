# mirBottleneck — Complete Beginner Explanation (What Was Done, Why, and How)

## 1) One-sentence summary
This project built and packaged a full R-based workflow that finds microRNAs (miRNAs) that seem to control how "stable" or "organized" gene expression is in pancreatic cancer, then tested whether those miRNA patterns are linked to patient survival.

---

## 2) Why this project exists
Pancreatic ductal adenocarcinoma (PDAC) is very deadly and needs better biomarkers.
Most biomarker workflows test one gene/miRNA at a time.
This project asks a different question:

- Instead of only asking "Is this miRNA high or low?"
- Ask "Does this miRNA make its target genes behave in a more controlled and coherent way across patients?"

So the package was designed to measure **regulatory influence at network scale**, not only single-marker association.

---

## 3) Core idea in plain language
Each miRNA gets **two scores**:

1. **VSS (Variance Suppression Score)**  
   "When this miRNA changes across patients, how much of the variation in its target genes does it explain?"

2. **Coherence Score**  
   "When this miRNA is high, do its target genes become more synchronized with each other?"

Then these two scores are combined to classify miRNAs into 4 behavior types:

- **Dual**: high VSS + high coherence
- **Silencer**: high VSS + low coherence
- **Conductor**: low VSS + high coherence
- **Weak**: low VSS + low coherence

After that, a **patient-level composite risk score** is built from survival-associated miRNAs and tested with survival models.

---

## 4) Data used and cohort construction
The project uses TCGA-PAAD data (RNA-seq, miRNA-seq, and clinical survival info).

Main processing decisions described in the manuscript and encoded in the pipeline:

- Keep primary tumor samples.
- Harmonize TCGA barcodes to patient-level IDs (first 12 chars).
- Keep only patients present across required data types.
- Use log2(count + 1) transformed expression matrices for analysis.
- Build final shared cohort for scoring/survival analysis.

Reported cohort values in manuscript:

- 178 patients with matched RNA + miRNA profiles
- 167 with complete survival data
- 164 in final multivariable model after missing-covariate exclusions

Clinical covariates used in Cox models:

- age
- stage_clean
- OS_days
- OS_status

---

## 5) Exactly what the package code does (end-to-end)
The orchestrator is `run_mirBottleneck_project()`.
It performs the full analysis in this sequence:

1. Read RDS inputs (`mirna_log`, `rna_sym`, `clinical_df`, `mirna_norm_map`).
2. Validate input structures and required columns.
3. Restrict to shared patients across datasets.
4. Build or load miRNA-target network.
5. Compute VSS for each miRNA.
6. Keep top miRNAs by VSS (parameterized) for coherence scoring.
7. Compute coherence scores + permutation p-values.
8. Merge/classify miRNAs into archetypes and compute bottleneck index.
9. Fit individual miRNA Cox models (adjusted for age + stage).
10. Select contributing miRNAs by nominal p-value threshold.
11. Build direction-aware weighted patient composite score.
12. Fit survival models (score-only, clinical-only, combined).
13. Save all intermediate/final outputs (RDS + CSV).
14. Optionally generate a complete offline HTML report.

---

## 6) Network construction details
Implemented in `build_network()`.

What it does:

- Queries validated miRNA-target interactions from multiMiR / miRTarBase.
- Uses batched queries to reduce timeout risk.
- Performs rescue queries for key known miRNAs.
- Normalizes miRNA IDs.
- Keeps only targets that exist in RNA expression data.
- Groups into a per-miRNA target list.
- Requires at least a minimum target count (default 5).

Output columns:

- `mirna`
- `targets` (list column)
- `n_targets`

---

## 7) VSS scoring details
Implemented in `score_vss()`.

For each miRNA:

- Get its expression vector across patients.
- For each validated target gene, fit linear model: `target_expression ~ miRNA_expression`.
- Take model R² for each target.
- VSS = mean R² across targets.

Interpretation:

- Higher VSS means miRNA expression tracks and explains more of target variation.

Practical filter in code:

- Skip miRNAs with fewer than 5 present targets.

---

## 8) Coherence scoring details
Implemented in `score_coherence()`.

For each miRNA:

- Split patients into high vs low miRNA expression (median split).
- Compute mean pairwise target-gene correlation in each group.
- Coherence score = (high-group mean correlation) - (low-group mean correlation).
- Estimate empirical p-value via permutations.

Practical filter in code:

- Skip/NA if fewer than 10 present targets.

Interpretation:

- Positive score: targets are more coordinated when miRNA is high.
- Negative score: targets are less coordinated when miRNA is high.

---

## 9) Classification + bottleneck index
Implemented in `classify_bottleneck()`.

Steps:

- Join VSS and coherence results by miRNA.
- Min-max normalize both metrics to [0,1].
- Compute bottleneck index = average(normalized VSS, normalized coherence).
- Use medians of normalized metrics to assign class (dual/silencer/conductor/weak).

This gives both:

- a continuous score (`bottleneck_index`)
- an interpretable categorical class (`class`)

---

## 10) Patient composite score (survival-oriented)
Implemented in `composite_score()`.

Key logic:

- Fit one Cox model per miRNA: survival ~ miRNA + age + stage.
- Keep miRNAs with nominal `p < cox_p_threshold` (default 0.15).
- Use Cox coefficient magnitude as weight.
- Use coefficient sign as direction (risk-increasing or risk-decreasing).
- For each patient, compute weighted mean of selected miRNA expressions.

So each patient gets one number representing the net weighted influence of selected bottleneck miRNAs.

---

## 11) Survival modeling
Implemented in `survival_model()`.

Three models are fit:

1. composite score only
2. clinical covariates only
3. composite + clinical

Also:

- Median split into high/low score groups
- Kaplan–Meier and log-rank test
- C-index and AIC reported

Purpose:

- test whether the molecular composite adds prognostic value beyond basic clinical variables.

---

## 12) Report generation
Implemented in `generate_report()` using `inst/report/mirBottleneck_report.Rmd`.

It produces a self-contained HTML report that includes:

- cohort summary
- network statistics
- VSS and coherence distributions
- phase diagram (class space)
- ranked miRNA tables
- survival metrics and KM curve
- per-patient composite score plot
- session info

Design goal: one file that a collaborator can open offline in a browser.

---

## 13) Extra analysis included in repository
Script: `inst/scripts/08_hallmark_entropy.R`

What it adds:

- Runs ssGSEA on Hallmark pathways.
- Computes per-patient Shannon entropy of pathway activity.
- Correlates entropy with bottleneck score.
- Tests entropy in Cox model.
- Computes pathway-level "disorder" association with bottleneck score.
- Saves result tables and figures.

This explores whether bottleneck score tracks broader transcriptome disorder.

---

## 14) Main manuscript results (plain-language recap)
From `mirBottleneck_manuscript_v19(1).docx`:

- Framework applied to TCGA-PAAD.
- 98 miRNAs retained after filtering.
- Known PDAC oncomiRs (like miR-21 and miR-155) rank highly in silencer behavior.
- Patient risk score significantly associated with survival in multivariable Cox model.
- Reported discrimination is moderate after optimism correction.
- Authors explicitly note this is discovery-stage and needs external validation.

Important interpretation note emphasized in manuscript:

- in-sample performance can be optimistic,
- so independent-cohort validation is the key next step.

---

## 15) What files are central in this repository
- `R/build_network.R`: validated target network creation
- `R/score_vss.R`: VSS computation
- `R/score_coherence.R`: coherence computation + permutation p-values
- `R/classify.R`: archetype assignment + bottleneck index
- `R/composite_score.R`: patient score construction from Cox-selected miRNAs
- `R/survival_model.R`: survival comparisons and model summaries
- `R/run_mirBottleneck_project.R`: full pipeline wrapper
- `R/generate_report.R`: HTML report renderer
- `inst/report/mirBottleneck_report.Rmd`: report template
- `vignettes/mirBottleneck-workflow.Rmd`: runnable tutorial with toy data
- `inst/extdata/*`: packaged toy and precomputed analysis objects
- `inst/scripts/08_hallmark_entropy.R`: entropy extension analysis

---

## 16) What “success” means in this project
The project is not just "fit one survival model".
It delivers a reproducible framework with:

- data harmonization assumptions,
- network-aware scoring,
- interpretable archetypes,
- patient-level risk aggregation,
- and automated reporting.

So success is:

1. reproducible end-to-end run,
2. biologically plausible miRNA ordering/classification,
3. statistically meaningful survival association,
4. transparent documentation of limitations and need for external validation.

---

## 17) Limitations (important for beginners)
The manuscript and code acknowledge several limitations:

- miRTarBase is literature-biased (well-studied miRNAs may be advantaged).
- Bulk RNA-seq mixes cell populations.
- Liberal selection threshold (p < 0.15) is discovery-oriented.
- Composite weights are learned in-sample.
- Therefore external validation is mandatory for clinical claims.

This is responsible practice: strong discovery signal, cautious interpretation.

---

## 18) If you are new: minimal mental model
Think of the pipeline like this:

1. Build each miRNA’s target neighborhood.
2. Measure whether each miRNA seems to control noise (VSS) and coordination (coherence).
3. Put miRNAs on a 2D map and assign behavior class.
4. Combine the most survival-relevant miRNAs into a patient score.
5. Test whether that score separates better vs worse survival.

That is the full scientific and computational story in one flow.

---

## 19) Deliverable created now
This file was created as a beginner-oriented "explain everything" companion document based on:

- the manuscript: `mirBottleneck_manuscript_v19(1).docx`
- and the repository source code / workflow files.

