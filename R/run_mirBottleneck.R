#' Run the full mirBottleneck pipeline for one cancer type
#'
#' Orchestrates the complete mirBottleneck analysis: network building, VSS,
#' CIS, classification, composite scoring, and survival modeling. Designed
#' to be called in a loop across cancer types. Results are checkpointed to
#' disk so Colab sessions can be resumed without recomputing completed cancers.
#'
#' @param cancer_type Character. TCGA project code, e.g. "TCGA-LUAD".
#' @param save_dir Character. Root directory for all outputs. A subdirectory
#'   named after the cancer type will be created automatically.
#' @param cohort_data Optional. Pre-fetched output of fetch_cohort_data().
#'   If NULL (default), fetch_cohort_data() is called automatically.
#' @param force_rerun Logical. If FALSE (default), skip cancers that already
#'   have a completed result RDS in save_dir.
#' @param cox_p_threshold Numeric. p-value cutoff for composite score miRNA
#'   selection (default 0.15).
#' @param n_perm Integer. Permutations for CIS empirical p-values (default 500).
#'   Use 200 for fast runs, 1000 for publication.
#' @param min_targets Integer. Minimum validated targets for a miRNA to be
#'   included (default 5).
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @return Invisibly returns a list with all pipeline outputs. Also writes
#'   an RDS file to \code{save_dir/<cancer_type>/result.rds} containing:
#'   \itemize{
#'     \item \code{cancer_type}
#'     \item \code{n_patients}
#'     \item \code{bottleneck_scores} - VSS + CIS per miRNA
#'     \item \code{classified}        - archetypes + bottleneck index
#'     \item \code{composite}         - patient scores + contributing miRNAs
#'     \item \code{survival}          - Cox models + log-rank + summary metrics
#'     \item \code{entropy}           - per-patient Shannon entropy
#'     \item \code{entropy_cor}       - Spearman rho + p for entropy vs score
#'     \item \code{clinical}          - clinical data frame used
#'     \item \code{timestamp}
#'   }
#'
#' @export
#' @examples
#' \donttest{
#' result <- run_mirBottleneck("TCGA-LUAD", save_dir = "/tmp/mirBottleneck_pancancer")
#' result$survival$summary
#' }
run_mirBottleneck <- function(cancer_type,
                              save_dir         = "mirBottleneck_results",
                              cohort_data      = NULL,
                              force_rerun      = FALSE,
                              cox_p_threshold  = 0.15,
                              n_perm           = 500,
                              min_targets      = 5,
                              verbose          = TRUE) {

  .vcat <- function(...) if (verbose) message("[", cancer_type, "] ", ...)

  out_dir  <- file.path(save_dir, cancer_type)
  rds_path <- file.path(out_dir, "result.rds")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------ #
  # Checkpoint: skip if already done                                    #
  # ------------------------------------------------------------------ #
  if (!force_rerun && file.exists(rds_path)) {
    .vcat("Checkpoint found. Loading cached result.")
    return(invisible(readRDS(rds_path)))
  }

  tryCatch({

    # ---------------------------------------------------------------- #
    # Step 1: Data                                                      #
    # ---------------------------------------------------------------- #
    if (is.null(cohort_data)) {
      .vcat("Fetching cohort data...")
      cohort_data <- fetch_cohort_data(cancer_type,
                                       save_dir = file.path(save_dir, "_gdc_cache",
                                                            cancer_type),
                                       verbose  = verbose)
    }

    mirna_log      <- cohort_data$mirna_log
    rna_log        <- cohort_data$rna_log
    clinical       <- cohort_data$clinical
    mirna_norm_map <- cohort_data$mirna_norm_map

    .vcat("Patients: ", cohort_data$n_patients,
          " | miRNAs: ", nrow(mirna_log),
          " | Genes: ", nrow(rna_log))

    # ---------------------------------------------------------------- #
    # Step 2: Build miRNA-target network                                #
    # ---------------------------------------------------------------- #
    .vcat("Building miRNA-target network...")
    network_path <- file.path(out_dir, "network.rds")

    if (file.exists(network_path)) {
      .vcat("  Loading cached network.")
      mirna_targets <- readRDS(network_path)
    } else {
      mirna_targets <- build_network(
        mirna_ids   = rownames(mirna_log),
        rna_symbols = rownames(rna_log),
        min_targets = min_targets
      )
      saveRDS(mirna_targets, network_path)
    }

    .vcat("miRNAs with sufficient targets: ", nrow(mirna_targets))

    if (nrow(mirna_targets) == 0) {
      stop("No miRNAs passed the min_targets filter. Check miRNA ID format.")
    }

    # ---------------------------------------------------------------- #
    # Step 3: VSS                                                       #
    # ---------------------------------------------------------------- #
    .vcat("Computing VSS...")
    vss_path <- file.path(out_dir, "vss.rds")

    if (file.exists(vss_path)) {
      vss <- readRDS(vss_path)
    } else {
      vss <- score_vss(mirna_log, rna_log, mirna_targets, mirna_norm_map)
      saveRDS(vss, vss_path)
    }

    # ---------------------------------------------------------------- #
    # Step 4: CIS                                                       #
    # ---------------------------------------------------------------- #
    .vcat("Computing CIS (n_perm = ", n_perm, ")...")
    cis_path <- file.path(out_dir, "cis.rds")

    if (file.exists(cis_path)) {
      cis <- readRDS(cis_path)
    } else {
      scored_mirnas <- vss$mirna
      cis <- score_coherence(mirna_log, rna_log, mirna_targets,
                             mirna_norm_map, scored_mirnas,
                             n_perm = n_perm)
      saveRDS(cis, cis_path)
    }

    # ---------------------------------------------------------------- #
    # Step 5: Classify archetypes                                       #
    # ---------------------------------------------------------------- #
    .vcat("Classifying archetypes...")
    classified <- classify_bottleneck(vss, cis)

    archetype_counts <- table(classified$class)
    .vcat("Archetypes: ",
          paste(names(archetype_counts), archetype_counts, sep = "=",
                collapse = " | "))

    # ---------------------------------------------------------------- #
    # Step 6: Composite score                                           #
    # ---------------------------------------------------------------- #
    .vcat("Building composite bottleneck score...")
    composite <- composite_score(
      mirna_log       = mirna_log,
      classified_df   = classified,
      clinical_df     = clinical,
      mirna_norm_map  = mirna_norm_map,
      cox_p_threshold = cox_p_threshold
    )

    .vcat("Contributing miRNAs: ", nrow(composite$contributing_mirnas))

    if (nrow(composite$contributing_mirnas) == 0) {
      warning(cancer_type, ": no miRNAs passed cox_p_threshold = ",
              cox_p_threshold, ". Survival results will be NA.")
    }

    # ---------------------------------------------------------------- #
    # Step 7: Survival analysis                                         #
    # ---------------------------------------------------------------- #
    .vcat("Running survival analysis...")
    survival <- tryCatch(
      survival_model(composite$patient_scores, clinical),
      error = function(e) {
        warning(cancer_type, " survival_model failed: ", e$message)
        list(summary = c(c_index_score = NA, logrank_p = NA,
                         hr = NA, aic_score = NA))
      }
    )

    # ---------------------------------------------------------------- #
    # Step 8: Entropy                                                   #
    # ---------------------------------------------------------------- #
    .vcat("Computing transcriptome entropy...")
    entropy_result <- compute_entropy(
      rna_log        = rna_log,
      patient_scores = composite$patient_scores,
      clinical       = clinical
    )

    # ---------------------------------------------------------------- #
    # Assemble and save result                                          #
    # ---------------------------------------------------------------- #
    result <- list(
      cancer_type       = cancer_type,
      n_patients        = cohort_data$n_patients,
      bottleneck_scores = merge(vss, cis, by = "mirna", all = FALSE),
      classified        = classified,
      composite         = composite,
      survival          = survival,
      entropy           = entropy_result$entropy_per_patient,
      entropy_cor       = entropy_result$spearman,
      clinical          = clinical,
      timestamp         = Sys.time()
    )

    saveRDS(result, rds_path)
    .vcat("Result saved to: ", rds_path)

    invisible(result)

  }, error = function(e) {
    message("[", cancer_type, "] PIPELINE FAILED: ", e$message)
    fail_path <- file.path(out_dir, "FAILED.txt")
    writeLines(c(as.character(Sys.time()), e$message), fail_path)
    invisible(NULL)
  })
}


#' Run mirBottleneck across all 33 TCGA cancer types
#'
#' Loops over all TCGA cancer types with miRNA-seq data, calls
#' run_mirBottleneck() for each, and aggregates results into a
#' pan-cancer summary table. Designed for Colab: checkpointing
#' means interrupted sessions resume from where they left off.
#'
#' @param save_dir Character. Root directory for all outputs.
#' @param cancer_types Character vector. Defaults to all 33 TCGA types
#'   with miRNA-seq and survival data.
#' @param n_perm Integer. CIS permutations (default 500).
#' @param cox_p_threshold Numeric. Composite score p cutoff (default 0.15).
#' @param force_rerun Logical. Rerun even if checkpoint exists (default FALSE).
#' @param verbose Logical (default TRUE).
#'
#' @return Data frame with one row per cancer type and columns:
#'   cancer_type, n_patients, c_index, hr, logrank_p, entropy_rho,
#'   entropy_p, n_dual, n_silencer, n_conductor, n_weak, status
#'
#' @export
run_pancancer <- function(save_dir        = "mirBottleneck_pancancer",
                          cancer_types    = .tcga_mirna_types(),
                          n_perm          = 500,
                          cox_p_threshold = 0.15,
                          force_rerun     = FALSE,
                          verbose         = TRUE) {

  message("=== mirBottleneck Pan-Cancer Run ===")
  message("Cancer types: ", length(cancer_types))
  message("Save dir: ", save_dir)
  message("Checkpoint-resumable: TRUE")
  message("")

  results_list <- vector("list", length(cancer_types))
  names(results_list) <- cancer_types

  for (i in seq_along(cancer_types)) {
    ct <- cancer_types[i]
    message("\n[", i, "/", length(cancer_types), "] Starting: ", ct)

    res <- run_mirBottleneck(
      cancer_type     = ct,
      save_dir        = save_dir,
      force_rerun     = force_rerun,
      cox_p_threshold = cox_p_threshold,
      n_perm          = n_perm,
      verbose         = verbose
    )

    results_list[[ct]] <- res
  }

  # Aggregate into summary table
  summary_df <- .aggregate_pancancer(results_list, cancer_types, save_dir)

  summary_path <- file.path(save_dir, "pancancer_summary.rds")
  saveRDS(summary_df, summary_path)
  message("\nPan-cancer summary saved to: ", summary_path)

  summary_df
}


# ------------------------------------------------------------------ #
# Internal: aggregate results                                         #
# ------------------------------------------------------------------ #
#' @keywords internal
.aggregate_pancancer <- function(results_list, cancer_types, save_dir) {

  rows <- lapply(cancer_types, function(ct) {

    # Try loaded result first, then disk
    res <- results_list[[ct]]
    if (is.null(res)) {
      rds_path <- file.path(save_dir, ct, "result.rds")
      if (file.exists(rds_path)) res <- readRDS(rds_path)
    }

    if (is.null(res)) {
      return(data.frame(cancer_type = ct, status = "FAILED",
                        stringsAsFactors = FALSE))
    }

    surv_sum <- res$survival$summary
    ec       <- res$entropy_cor
    cls      <- table(res$classified$class)

    data.frame(
      cancer_type  = ct,
      n_patients   = res$n_patients,
      c_index      = .safe_get(surv_sum, "c_index_score"),
      hr           = .safe_get(surv_sum, "hr"),
      logrank_p    = .safe_get(surv_sum, "logrank_p"),
      aic_score    = .safe_get(surv_sum, "aic_score"),
      entropy_rho  = if (!is.null(ec)) ec$estimate else NA_real_,
      entropy_p    = if (!is.null(ec)) ec$p.value  else NA_real_,
      n_dual       = .safe_table(cls, "dual"),
      n_silencer   = .safe_table(cls, "silencer"),
      n_conductor  = .safe_table(cls, "conductor"),
      n_weak       = .safe_table(cls, "weak"),
      status       = "OK",
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

#' @keywords internal
.safe_get <- function(x, key) {
  if (is.null(x) || !key %in% names(x)) return(NA_real_)
  as.numeric(x[[key]])
}

#' @keywords internal
.safe_table <- function(tbl, key) {
  if (key %in% names(tbl)) as.integer(tbl[[key]]) else 0L
}

#' @keywords internal
.tcga_mirna_types <- function() {
  c(
    "TCGA-ACC",  "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
    "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM",  "TCGA-HNSC",
    "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
    "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
    "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
    "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
    "TCGA-UCEC", "TCGA-UCS",  "TCGA-UVM"
  )
}
