#' Run mirBottleneck on a real project (RDS inputs)
#'
#' Reads user-provided `.rds` objects (no demos) and writes results to `out_dir`.
#' Supports either a user-supplied miRNA-target network (offline deterministic)
#' or online network building via `multiMiR` (miRTarBase) through `build_network()`.
#'
#' Expected objects:
#' - mirna_log: numeric matrix, miRNAs x patients (rownames miRNA IDs, colnames patient IDs)
#' - rna_sym: numeric matrix, genes x patients (rownames gene symbols, colnames patient IDs)
#' - clinical_df: data.frame with patient, OS_days, OS_status, age, stage_clean
#' - mirna_norm_map: data.frame with original, norm
#'
#' @param mirna_log_rds Path to `.rds` containing `mirna_log`
#' @param rna_rds Path to `.rds` containing `rna_sym`
#' @param clinical_rds Path to `.rds` containing `clinical_df`
#' @param mirna_norm_map_rds Path to `.rds` containing `mirna_norm_map`
#' @param out_dir Output directory (created if needed)
#' @param mirna_targets_rds Optional path to `.rds` containing `mirna_targets`
#'   as returned by `build_network()` (columns: mirna, targets(list), n_targets)
#' @param query_network If TRUE and `mirna_targets_rds` is NULL, build network online.
#' @param batch_size Passed to `build_network()`
#' @param min_targets Passed to `build_network()`
#' @param rescue_mirnas Passed to `build_network()`
#' @param vss_keep_top Integer, keep top-N miRNAs by VSS for coherence scoring.
#' @param n_perm Passed to `compute_cis()`
#' @param max_targets Passed to `compute_cis()`.
#' @param BPPARAM_vss Optional BiocParallel parameter for `compute_vss()` (serial fallback).
#' @param BPPARAM_cis Optional BiocParallel parameter for `compute_cis()` (serial fallback).
#' @param cox_p_threshold p-value threshold for composite_score().
#' @param report Logical. If TRUE, generate an HTML report.
#' @param cohort_name Character label used in the HTML report.
#'
#' @return Invisible list of key objects.
#' @export
#' @examples
#' paths <- mirBottleneck:::.toy_paths()
#' run_mirBottleneck_project(
#'   mirna_log_rds      = paths$mirna_log,
#'   rna_rds            = paths$rna_sym,
#'   clinical_rds       = paths$clinical_df,
#'   mirna_norm_map_rds = paths$mirna_norm_map,
#'   mirna_targets_rds  = paths$mirna_targets,
#'   query_network      = FALSE,
#'   out_dir            = tempdir(),
#'   cox_p_threshold    = 0.5,
#'   report             = FALSE
#' )
run_mirBottleneck_project <- function(
  mirna_log_rds,
  rna_rds,
  clinical_rds,
  mirna_norm_map_rds,
  out_dir,
  mirna_targets_rds = NULL,
  query_network = TRUE,
  batch_size = 50,
  min_targets = 5,
  rescue_mirnas = c("hsa-miR-21-5p",  "hsa-miR-21-3p",
                    "hsa-miR-155-5p", "hsa-miR-155-3p",
                    "hsa-miR-196a-5p","hsa-miR-196b-5p",
                    "hsa-miR-210-3p"),
  vss_keep_top = 300,
  n_perm = 200,
  max_targets = 200,
  BPPARAM_vss = NULL,
  BPPARAM_cis = NULL,
  cox_p_threshold = 0.15,
  report = FALSE,
  cohort_name = "TCGA-PAAD"
) {
  if (missing(out_dir) || is.null(out_dir) || !nzchar(out_dir)) {
    stop("out_dir is required.")
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Load inputs ----
  mirna_log      <- readRDS(mirna_log_rds)
  rna_sym        <- readRDS(rna_rds)
  clinical_df    <- readRDS(clinical_rds)
  mirna_norm_map <- readRDS(mirna_norm_map_rds)

  # Basic validation
  if (!is.matrix(mirna_log)) stop("mirna_log must be a matrix.")
  if (!is.matrix(rna_sym)) stop("rna_sym must be a matrix.")
  req_cols <- c("patient", "OS_days", "OS_status", "age", "stage_clean")
  if (!all(req_cols %in% colnames(clinical_df))) {
    stop("clinical_df must contain columns: ", paste(req_cols, collapse = ", "))
  }
  if (!all(c("original", "norm") %in% colnames(mirna_norm_map))) {
    stop("mirna_norm_map must contain columns: original, norm")
  }

  # Restrict to shared patients
  pats <- intersect(colnames(mirna_log), colnames(rna_sym))
  pats <- intersect(pats, clinical_df$patient)
  if (length(pats) < 10) stop("Too few shared patients after intersection.")
  mirna_log   <- mirna_log[, pats, drop = FALSE]
  rna_sym     <- rna_sym[, pats, drop = FALSE]
  clinical_df <- clinical_df[clinical_df$patient %in% pats, , drop = FALSE]

  # ---- Network (C3) ----
  if (!is.null(mirna_targets_rds)) {
    mirna_targets <- readRDS(mirna_targets_rds)
  } else {
    if (!isTRUE(query_network)) {
      stop("mirna_targets_rds is NULL and query_network=FALSE. Provide a network or enable query_network.")
    }
    mirna_ids   <- rownames(mirna_log)
    rna_symbols <- rownames(rna_sym)
    mirna_targets <- build_network(
      mirna_ids = mirna_ids,
      rna_symbols = rna_symbols,
      batch_size = batch_size,
      min_targets = min_targets,
      rescue_mirnas = rescue_mirnas
    )
  }

  # Persist network used
  saveRDS(mirna_targets, file = file.path(out_dir, "mirna_targets.rds"))

  # ---- Scoring pipeline ----
  vss_scores <- compute_vss(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map,
    BPPARAM = BPPARAM_vss
  )
  saveRDS(vss_scores, file = file.path(out_dir, "vss_scores.rds"))
  utils::write.csv(vss_scores, file = file.path(out_dir, "vss_scores.csv"), row.names = FALSE)

  scored_mirnas <- head(vss_scores$mirna, n = min(vss_keep_top, nrow(vss_scores)))

  coherence_scores <- compute_cis(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map,
    scored_mirnas = scored_mirnas,
    n_perm = n_perm,
    max_targets = max_targets,
    BPPARAM = BPPARAM_cis
  )
  saveRDS(coherence_scores, file = file.path(out_dir, "coherence_scores.rds"))
  utils::write.csv(coherence_scores, file = file.path(out_dir, "coherence_scores.csv"), row.names = FALSE)

  combined <- classify_archetypes(vss_scores = vss_scores, coherence_scores = coherence_scores)
  saveRDS(combined, file = file.path(out_dir, "mirna_bottleneck.rds"))
  utils::write.csv(combined, file = file.path(out_dir, "mirna_bottleneck.csv"), row.names = FALSE)

  comp <- composite_score(
    mirna_log = mirna_log,
    classified_df = combined,
    clinical_df = clinical_df,
    mirna_norm_map = mirna_norm_map,
    cox_p_threshold = cox_p_threshold
  )
  saveRDS(comp, file = file.path(out_dir, "composite_score.rds"))
  utils::write.csv(comp$contributing_mirnas,
                   file = file.path(out_dir, "contributing_mirnas.csv"),
                   row.names = FALSE)
  utils::write.csv(comp$cox_individual,
                   file = file.path(out_dir, "cox_individual.csv"),
                   row.names = FALSE)
  utils::write.csv(data.frame(patient = names(comp$patient_scores),
                              composite_score = as.numeric(comp$patient_scores)),
                   file = file.path(out_dir, "patient_scores.csv"),
                   row.names = FALSE)

  surv <- survival_model(patient_scores = comp$patient_scores, clinical_df = clinical_df)
  saveRDS(surv, file = file.path(out_dir, "survival_results.rds"))
  utils::write.csv(surv$surv_df, file = file.path(out_dir, "survival_df.csv"), row.names = FALSE)

  results <- list(
    combined_bottleneck = combined,
    cox_individual      = comp$cox_individual,
    cox_final           = surv$cox_combined,
    survival_df         = surv$surv_df,
    mirna_targets       = mirna_targets,
    contributing_mirnas = comp$contributing_mirnas,
    clinical_df         = clinical_df,
    vss_scores          = vss_scores,
    coherence_scores    = coherence_scores,
    patient_scores      = comp$patient_scores,
    composite           = comp,
    survival            = surv
  )

  if (isTRUE(report)) {
    generate_report(
      results     = results,
      out_dir     = out_dir,
      cohort_name = cohort_name,
      open        = FALSE
    )
  }

  invisible(results)
}
