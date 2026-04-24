#' Bootstrap optimism correction for the composite bottleneck C-index
#'
#' Estimates optimism-corrected model performance using the Harrell bootstrap
#' optimism correction method. For each bootstrap sample: fits the composite
#' score model, evaluates apparent performance on the bootstrap sample and
#' test performance on the original data, then subtracts the mean optimism
#' from the apparent (training) C-index.
#'
#' @param mirna_log Numeric matrix (miRNAs x patients), log2-transformed.
#' @param classified_df Data frame from classify_bottleneck().
#' @param clinical_df Data frame with patient, OS_days, OS_status, age,
#'   stage_clean.
#' @param mirna_norm_map Data frame with original and norm columns.
#' @param cox_p_threshold Numeric. p-value cutoff for miRNA selection
#'   (default 0.15). Must match what was used in composite_score().
#' @param n_boot Integer. Number of bootstrap replicates (default 200).
#'   Use 1000 for publication.
#' @param seed Integer. Random seed for reproducibility (default 42).
#' @param verbose Logical (default TRUE).
#'
#' @return List with:
#'   \itemize{
#'     \item \code{c_index_apparent}   - naive C-index on full data
#'     \item \code{c_index_corrected}  - optimism-corrected C-index
#'     \item \code{optimism}           - mean optimism across bootstrap samples
#'     \item \code{optimism_sd}        - SD of optimism
#'     \item \code{boot_results}       - data frame of per-bootstrap metrics
#'     \item \code{n_boot}             - number of successful bootstrap runs
#'   }
#'
#' @export
#' @importFrom survival coxph Surv
#' @importFrom stats median
#'
#' @examples
#' \donttest{
#' cohort <- fetch_cohort_data("TCGA-PAAD", save_dir = "/tmp/tcga")
#' classified <- classify_bottleneck(
#'   compute_vss(cohort$mirna_log, cohort$rna_log,
#'               build_network(rownames(cohort$mirna_log),
#'                             rownames(cohort$rna_log)),
#'               cohort$mirna_norm_map),
#'   compute_cis(cohort$mirna_log, cohort$rna_log,
#'               build_network(rownames(cohort$mirna_log),
#'                             rownames(cohort$rna_log)),
#'               cohort$mirna_norm_map,
#'               n_perm = 200)
#' )
#' boot_result <- bootstrap_validate(
#'   mirna_log      = cohort$mirna_log,
#'   classified_df  = classified,
#'   clinical_df    = cohort$clinical,
#'   mirna_norm_map = cohort$mirna_norm_map,
#'   n_boot         = 200
#' )
#' boot_result$c_index_corrected
#' }
bootstrap_validate <- function(mirna_log,
                                classified_df,
                                clinical_df,
                                mirna_norm_map,
                                cox_p_threshold = 0.15,
                                n_boot          = 200,
                                seed            = 42,
                                verbose         = TRUE) {

  set.seed(seed)

  patients <- intersect(colnames(mirna_log), clinical_df$patient)
  patients <- patients[
    !is.na(clinical_df$OS_days[match(patients, clinical_df$patient)]) &
    clinical_df$OS_days[match(patients, clinical_df$patient)] > 0
  ]
  n <- length(patients)

  if (n < 30) warning("n=", n, " is small for bootstrap validation.")

  .vcat <- function(...) if (verbose) message(...)

  .vcat("Bootstrap optimism correction: n=", n,
        " patients, ", n_boot, " replicates")

  # --- Apparent performance on full data ---
  full_composite <- composite_score(
    mirna_log       = mirna_log[, patients, drop = FALSE],
    classified_df   = classified_df,
    clinical_df     = clinical_df[clinical_df$patient %in% patients, ],
    mirna_norm_map  = mirna_norm_map,
    cox_p_threshold = cox_p_threshold
  )

  if (nrow(full_composite$contributing_mirnas) == 0) {
    stop("No miRNAs selected with cox_p_threshold=", cox_p_threshold,
         ". Try a more lenient threshold.")
  }

  full_surv    <- survival_model(full_composite$patient_scores,
                                  clinical_df[clinical_df$patient %in% patients, ])
  c_apparent   <- full_surv$summary["c_index_score"]
  contrib_mirs <- full_composite$contributing_mirnas

  .vcat("Apparent C-index: ", round(c_apparent, 4))

  # --- Bootstrap loop ---
  boot_results <- vector("list", n_boot)
  n_success    <- 0L

  pb_step <- max(1L, n_boot %/% 10L)

  for (b in seq_len(n_boot)) {

    if (verbose && b %% pb_step == 0)
      message("  Bootstrap ", b, "/", n_boot, "...")

    boot_idx  <- sample(n, n, replace = TRUE)
    boot_pts  <- patients[boot_idx]
    orig_pts  <- patients  # for test evaluation

    boot_mirna <- mirna_log[, boot_idx, drop = FALSE]
    colnames(boot_mirna) <- boot_pts

    boot_clin <- clinical_df[match(boot_pts, clinical_df$patient), ]
    boot_clin$patient <- boot_pts

    # Deduplicate for composite_score (it needs unique patient IDs)
    # We aggregate by taking mean for duplicate patients
    boot_mirna_dedup <- .dedup_boot_matrix(boot_mirna)
    boot_clin_dedup  <- boot_clin[!duplicated(boot_clin$patient), ]

    result_b <- tryCatch({

      # Fit model on bootstrap sample
      comp_b <- composite_score(
        mirna_log       = boot_mirna_dedup,
        classified_df   = classified_df,
        clinical_df     = boot_clin_dedup,
        mirna_norm_map  = mirna_norm_map,
        cox_p_threshold = cox_p_threshold
      )

      if (nrow(comp_b$contributing_mirnas) == 0) return(NULL)

      # Apparent performance on bootstrap sample
      surv_b_app <- survival_model(comp_b$patient_scores, boot_clin_dedup)
      c_boot_app <- surv_b_app$summary["c_index_score"]

      # Test performance: apply bootstrap model weights to ORIGINAL data
      # Use locked weights from bootstrap model
      c_boot_test <- .apply_locked_score(
        contributing_mirnas = comp_b$contributing_mirnas,
        mirna_log           = mirna_log[, orig_pts, drop = FALSE],
        clinical_df         = clinical_df[clinical_df$patient %in% orig_pts, ],
        mirna_norm_map      = mirna_norm_map
      )

      data.frame(
        boot         = b,
        c_apparent   = c_boot_app,
        c_test       = c_boot_test,
        optimism     = c_boot_app - c_boot_test,
        n_mirnas     = nrow(comp_b$contributing_mirnas),
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      if (verbose) message("  Bootstrap ", b, " failed: ", e$message)
      NULL
    })

    if (!is.null(result_b)) {
      n_success <- n_success + 1L
      boot_results[[b]] <- result_b
    }
  }

  boot_df <- dplyr::bind_rows(Filter(Negate(is.null), boot_results))

  if (nrow(boot_df) == 0)
    stop("All bootstrap replicates failed. Check inputs.")

  mean_optimism <- mean(boot_df$optimism, na.rm = TRUE)
  sd_optimism   <- sd(boot_df$optimism,   na.rm = TRUE)
  c_corrected   <- c_apparent - mean_optimism

  .vcat("\nBootstrap optimism correction results:")
  .vcat("  Successful replicates : ", n_success, "/", n_boot)
  .vcat("  Mean optimism         : ", round(mean_optimism, 4))
  .vcat("  SD optimism           : ", round(sd_optimism,   4))
  .vcat("  Apparent C-index      : ", round(c_apparent,    4))
  .vcat("  Corrected C-index     : ", round(c_corrected,   4))

  list(
    c_index_apparent  = as.numeric(c_apparent),
    c_index_corrected = as.numeric(c_corrected),
    optimism          = mean_optimism,
    optimism_sd       = sd_optimism,
    boot_results      = boot_df,
    n_boot            = n_success
  )
}


# ------------------------------------------------------------------ #
# Internal helpers                                                    #
# ------------------------------------------------------------------ #

#' @keywords internal
.dedup_boot_matrix <- function(mat) {
  # Average expression for duplicate patient IDs (from bootstrap resampling)
  pts <- colnames(mat)
  if (!anyDuplicated(pts)) return(mat)

  unique_pts <- unique(pts)
  out <- vapply(unique_pts, function(p) {
    idx <- which(pts == p)
    if (length(idx) == 1L) return(mat[, idx])
    rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(mat)))

  out <- t(out)
  rownames(out) <- rownames(mat)
  out
}

#' @keywords internal
.apply_locked_score <- function(contributing_mirnas, mirna_log,
                                 clinical_df, mirna_norm_map) {

  # Normalize IDs
  ext_norm <- normalize_mirna(rownames(mirna_log))
  rownames(mirna_log) <- ext_norm

  train_mirnas <- normalize_mirna(contributing_mirnas$mirna)
  found        <- intersect(train_mirnas, ext_norm)

  if (length(found) == 0) return(NA_real_)

  cont_sub       <- contributing_mirnas
  cont_sub$mirna_norm <- train_mirnas
  cont_sub       <- cont_sub[cont_sub$mirna_norm %in% found, ]

  shared_pts     <- intersect(colnames(mirna_log), clinical_df$patient)
  mirna_sub      <- mirna_log[cont_sub$mirna_norm, shared_pts, drop = FALSE]

  patient_scores <- apply(mirna_sub, 2, function(expr_vec) {
    weights  <- abs(cont_sub$coef)
    directed <- expr_vec * sign(cont_sub$coef)
    stats::weighted.mean(directed, w = weights, na.rm = TRUE)
  })

  clin_sub <- clinical_df[clinical_df$patient %in% shared_pts, ]

  cox_test <- tryCatch(
    survival::coxph(
      survival::Surv(OS_days, OS_status) ~ composite,
      data = data.frame(
        composite  = patient_scores[clin_sub$patient],
        OS_days    = clin_sub$OS_days,
        OS_status  = clin_sub$OS_status
      )
    ),
    error = function(e) NULL
  )

  if (is.null(cox_test)) return(NA_real_)
  as.numeric(summary(cox_test)$concordance[1])
}
