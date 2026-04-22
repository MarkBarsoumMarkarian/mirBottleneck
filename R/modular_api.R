#' Compute Variance Suppression Scores (modular API)
#'
#' Modular entry point for VSS computation. Accepts either classic matrix/data.frame
#' inputs (`mirna_log`, `rna_sym`) or a `SummarizedExperiment` carrying miRNA and
#' mRNA assays.
#'
#' @param mirna_log Numeric matrix/data.frame of miRNA expression (miRNAs x samples).
#' @param rna_sym Numeric matrix/data.frame of mRNA expression (genes x samples).
#' @param mirna_targets Data frame with columns `mirna`, `targets` (list), and `n_targets`.
#' @param mirna_norm_map Data frame with columns `original`, `norm`.
#' @param se Optional `SummarizedExperiment` with miRNA and mRNA assays.
#' @param mirna_assay Assay name for miRNA expression when `se` is used.
#' @param mrna_assay Assay name for mRNA expression when `se` is used.
#' @param BPPARAM Optional BiocParallel parameter object. Currently reserved;
#'   computation falls back to serial execution.
#' @return Data frame with columns: `mirna`, `vss`, `n_targets`, ordered by decreasing `vss`.
#' @export
compute_vss <- function(mirna_log = NULL,
                        rna_sym = NULL,
                        mirna_targets,
                        mirna_norm_map,
                        se = NULL,
                        mirna_assay = "mirna",
                        mrna_assay = "mrna",
                        BPPARAM = NULL) {
  resolved <- .resolve_expression_inputs(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    se = se,
    mirna_assay = mirna_assay,
    mrna_assay = mrna_assay,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map
  )
  mirna_log <- resolved$mirna_log
  rna_sym <- resolved$rna_sym

  if (!is.null(BPPARAM)) {
    message("BPPARAM provided; compute_vss currently executes in serial mode.")
  }

  mirna_expr <- .build_mirna_expr(mirna_log, mirna_targets$mirna, mirna_norm_map)

  n <- nrow(mirna_targets)
  scores <- numeric(n)
  n_tested <- integer(n)

  message("Computing VSS for ", n, " miRNAs...")

  for (i in seq_len(n)) {
    mid <- mirna_targets$mirna[i]
    targets <- mirna_targets$targets[[i]]

    if (!mid %in% names(mirna_expr)) next
    x <- mirna_expr[[mid]]

    targets_present <- targets[targets %in% rownames(rna_sym)]
    if (length(targets_present) < 5) next

    target_mat <- rna_sym[targets_present, , drop = FALSE]

    r2_vals <- apply(target_mat, 1, function(y) {
      tryCatch(summary(lm(y ~ x))$r.squared, error = function(e) NA_real_)
    })

    scores[i] <- mean(r2_vals, na.rm = TRUE)
    n_tested[i] <- length(targets_present)

    if (i %% 100 == 0) message("  ", i, " / ", n)
  }

  result <- data.frame(
    mirna = mirna_targets$mirna,
    vss = scores,
    n_targets = n_tested,
    stringsAsFactors = FALSE
  )
  result <- result[result$vss > 0, , drop = FALSE]
  result[order(-result$vss), , drop = FALSE]
}

#' Compute Coherence Induction Scores (modular API)
#'
#' Modular entry point for coherence scoring. Accepts either classic matrix/data.frame
#' inputs (`mirna_log`, `rna_sym`) or a `SummarizedExperiment` carrying miRNA and
#' mRNA assays.
#'
#' @param mirna_log Numeric matrix/data.frame of miRNA expression (miRNAs x samples).
#' @param rna_sym Numeric matrix/data.frame of mRNA expression (genes x samples).
#' @param mirna_targets Data frame with columns `mirna`, `targets` (list), and `n_targets`.
#' @param mirna_norm_map Data frame with columns `original`, `norm`.
#' @param scored_mirnas Character vector of miRNA IDs to score.
#' @param n_perm Number of permutations for empirical p-values.
#' @param max_targets Maximum targets per miRNA used in coherence scoring.
#' @param se Optional `SummarizedExperiment` with miRNA and mRNA assays.
#' @param mirna_assay Assay name for miRNA expression when `se` is used.
#' @param mrna_assay Assay name for mRNA expression when `se` is used.
#' @param BPPARAM Optional BiocParallel parameter object. Currently reserved;
#'   computation falls back to serial execution.
#' @return Data frame with columns: `mirna`, `coherence_score`, `coherence_p`.
#' @export
compute_cis <- function(mirna_log = NULL,
                        rna_sym = NULL,
                        mirna_targets,
                        mirna_norm_map,
                        scored_mirnas,
                        n_perm = 200,
                        max_targets = 200,
                        se = NULL,
                        mirna_assay = "mirna",
                        mrna_assay = "mrna",
                        BPPARAM = NULL) {
  resolved <- .resolve_expression_inputs(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    se = se,
    mirna_assay = mirna_assay,
    mrna_assay = mrna_assay,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map,
    scored_mirnas = scored_mirnas
  )
  mirna_log <- resolved$mirna_log
  rna_sym <- resolved$rna_sym

  if (!is.null(BPPARAM)) {
    message("BPPARAM provided; compute_cis currently executes in serial mode.")
  }

  mirna_expr <- .build_mirna_expr(mirna_log, scored_mirnas, mirna_norm_map)

  n <- length(scored_mirnas)
  coherence_score <- numeric(n)
  coherence_p <- numeric(n)

  message("Computing coherence scores for ", n, " miRNAs...")

  for (i in seq_len(n)) {
    mid <- scored_mirnas[i]

    tryCatch({
      if (!mid %in% names(mirna_expr)) next
      x <- mirna_expr[[mid]]

      targ_row <- mirna_targets[mirna_targets$mirna == mid, , drop = FALSE]
      if (nrow(targ_row) == 0) next
      targets <- targ_row$targets[[1]]
      targets_present <- targets[targets %in% rownames(rna_sym)]

      if (length(targets_present) < 10) {
        coherence_score[i] <- NA_real_
        next
      }

      if (length(targets_present) > max_targets) {
        tv <- apply(rna_sym[targets_present, , drop = FALSE], 1, var)
        targets_present <- names(sort(tv, decreasing = TRUE))[seq_len(max_targets)]
      }

      target_mat <- rna_sym[targets_present, , drop = FALSE]
      if (any(dim(target_mat) == 0)) next

      high_idx <- which(x >= median(x))
      low_idx <- which(x < median(x))

      ch <- mean_pairwise_cor(target_mat[, high_idx, drop = FALSE])
      cl <- mean_pairwise_cor(target_mat[, low_idx, drop = FALSE])
      if (is.na(ch) || is.na(cl)) next

      obs_delta <- ch - cl

      perm_deltas <- numeric(n_perm)
      for (p in seq_len(n_perm)) {
        px <- sample(x)
        ph <- mean_pairwise_cor(target_mat[, which(px >= median(px)), drop = FALSE])
        pl <- mean_pairwise_cor(target_mat[, which(px < median(px)), drop = FALSE])
        if (!is.na(ph) && !is.na(pl)) perm_deltas[p] <- ph - pl
      }

      coherence_score[i] <- obs_delta
      coherence_p[i] <- mean(perm_deltas >= obs_delta, na.rm = TRUE)
    }, error = function(e) {
      message("  Skipping ", mid, ": ", e$message)
    })

    if (i %% 100 == 0) message("  ", i, " / ", n)
  }

  result <- data.frame(
    mirna = scored_mirnas,
    coherence_score = coherence_score,
    coherence_p = coherence_p,
    stringsAsFactors = FALSE
  )

  result <- result[!(result$coherence_score == 0 &
                     result$coherence_p == 0 &
                     !is.na(result$coherence_score)), , drop = FALSE]
  result <- result[!is.na(result$coherence_score), , drop = FALSE]
  result[order(-result$coherence_score), , drop = FALSE]
}

#' Classify miRNA archetypes (modular API)
#'
#' Joins VSS and coherence outputs, computes normalized bottleneck index,
#' and assigns each miRNA to one of four archetypes.
#'
#' @param vss_scores Data frame from `compute_vss()` or `score_vss()`.
#' @param coherence_scores Data frame from `compute_cis()` or `score_coherence()`.
#' @return Data frame containing scores, normalized components, bottleneck index,
#'   and assigned class.
#' @export
classify_archetypes <- function(vss_scores, coherence_scores) {
  n_vss <- nrow(vss_scores)
  n_coh <- nrow(coherence_scores)
  combined <- dplyr::inner_join(vss_scores, coherence_scores, by = "mirna")
  if (nrow(combined) < min(n_vss, n_coh)) {
    message("Dropped ", min(n_vss, n_coh) - nrow(combined),
            " miRNAs not shared between VSS and coherence inputs.")
  }
  combined <- combined |>
    dplyr::mutate(
      vss_norm = normalize_01(vss),
      coherence_norm = normalize_01(coherence_score),
      bottleneck_index = (vss_norm + coherence_norm) / 2
    )

  vss_med <- median(combined$vss_norm, na.rm = TRUE)
  coh_med <- median(combined$coherence_norm, na.rm = TRUE)

  combined <- combined |>
    dplyr::mutate(
      class = dplyr::case_when(
        vss_norm >= vss_med & coherence_norm >= coh_med ~ "dual",
        vss_norm >= vss_med & coherence_norm < coh_med ~ "silencer",
        vss_norm < vss_med & coherence_norm >= coh_med ~ "conductor",
        TRUE ~ "weak"
      )
    ) |>
    dplyr::arrange(dplyr::desc(bottleneck_index))

  message("Class distribution:")
  print(table(combined$class))
  combined
}
