#' Compute Coherence Induction Scores
#'
#' Splits patients by median miRNA expression and measures whether target genes
#' become more correlated in the high-expression group. The coherence score is
#' the difference in mean pairwise target correlation between high and low groups.
#' An empirical p-value is computed via permutation.
#'
#' @param mirna_log Numeric matrix of log2-transformed miRNA expression
#'   (miRNAs x patients).
#' @param rna_sym Numeric matrix of log2-transformed RNA-seq expression
#'   (gene symbols x patients).
#' @param mirna_targets Data frame from \code{build_network()}.
#' @param mirna_norm_map Data frame mapping precursor IDs to normalized IDs.
#' @param scored_mirnas Character vector of miRNA IDs to score (typically those
#'   passing VSS threshold).
#' @param n_perm Integer. Number of permutations for p-value (default 200).
#' @param max_targets Integer. Maximum targets per miRNA for speed (default 200).
#' @return Data frame with columns: mirna, coherence_score, coherence_p
#' @export
#' @examples
#' set.seed(1)
#' mirna_log <- matrix(rnorm(4 * 8), nrow = 4)
#' rownames(mirna_log) <- c("hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-miR-1-3p", "hsa-miR-2-3p")
#' colnames(mirna_log) <- paste0("P", 1:8)
#'
#' rna_sym <- matrix(rnorm(6 * 8), nrow = 6)
#' rownames(rna_sym) <- c("KRAS", "TP53", "SMAD4", "CDKN2A", "EGFR", "BRCA1")
#' colnames(rna_sym) <- paste0("P", 1:8)
#'
#' mirna_targets <- data.frame(
#'   mirna = c("hsa-miR-21-5p", "hsa-miR-155-5p"),
#'   targets = I(list(c("KRAS", "TP53"), c("EGFR", "BRCA1"))),
#'   n_targets = c(2L, 2L)
#' )
#'
#' mirna_norm_map <- data.frame(
#'   original = rownames(mirna_log),
#'   norm = rownames(mirna_log)
#' )
#'
#' scored_mirnas <- c("hsa-miR-21-5p", "hsa-miR-155-5p")
#'
#' score_coherence(
#'   mirna_log, rna_sym, mirna_targets, mirna_norm_map,
#'   scored_mirnas = scored_mirnas,
#'   n_perm = 5,
#'   max_targets = 10
#' )
score_coherence <- function(mirna_log, rna_sym, mirna_targets, mirna_norm_map,
                            scored_mirnas, n_perm = 200, max_targets = 200) {

  mirna_expr <- .build_mirna_expr(mirna_log, scored_mirnas, mirna_norm_map)

  n               <- length(scored_mirnas)
  coherence_score <- numeric(n)
  coherence_p     <- numeric(n)

  message("Computing coherence scores for ", n, " miRNAs...")

  for (i in seq_len(n)) {
    mid <- scored_mirnas[i]

    tryCatch({
      if (!mid %in% names(mirna_expr)) next
      x <- mirna_expr[[mid]]

      targ_row <- mirna_targets[mirna_targets$mirna == mid, ]
      if (nrow(targ_row) == 0) next
      targets         <- targ_row$targets[[1]]
      targets_present <- targets[targets %in% rownames(rna_sym)]

      if (length(targets_present) < 10) {
        coherence_score[i] <- NA
        next
      }

      if (length(targets_present) > max_targets) {
        tv              <- apply(rna_sym[targets_present, ], 1, var)
        targets_present <- names(sort(tv, decreasing = TRUE))[seq_len(max_targets)]
      }

      target_mat <- rna_sym[targets_present, , drop = FALSE]
      if (any(dim(target_mat) == 0)) next

      high_idx <- which(x >= median(x))
      low_idx  <- which(x <  median(x))

      ch <- mean_pairwise_cor(target_mat[, high_idx, drop = FALSE])
      cl <- mean_pairwise_cor(target_mat[, low_idx,  drop = FALSE])
      if (is.na(ch) || is.na(cl)) next

      obs_delta <- ch - cl

      perm_deltas <- numeric(n_perm)
      for (p in seq_len(n_perm)) {
        px <- sample(x)
        ph <- mean_pairwise_cor(target_mat[, which(px >= median(px)), drop = FALSE])
        pl <- mean_pairwise_cor(target_mat[, which(px <  median(px)), drop = FALSE])
        if (!is.na(ph) && !is.na(pl)) perm_deltas[p] <- ph - pl
      }

      coherence_score[i] <- obs_delta
      coherence_p[i]     <- mean(perm_deltas >= obs_delta, na.rm = TRUE)

    }, error = function(e) {
      message("  Skipping ", mid, ": ", e$message)
    })

    if (i %% 100 == 0) message("  ", i, " / ", n)
  }

  result <- data.frame(mirna           = scored_mirnas,
                       coherence_score = coherence_score,
                       coherence_p     = coherence_p,
                       stringsAsFactors = FALSE)

  # Remove degenerate entries
  result <- result[!(result$coherence_score == 0 &
                     result$coherence_p     == 0 &
                     !is.na(result$coherence_score)), ]
  result <- result[!is.na(result$coherence_score), ]
  result[order(-result$coherence_score), ]
}
