#' Compute Variance Suppression Scores
#'
#' For each miRNA, fits a linear model of each target gene's expression on
#' the miRNA's expression across patients. The VSS is the mean R-squared
#' across all validated targets, quantifying how much transcriptome variance
#' the miRNA explains.
#'
#' @param mirna_log Numeric matrix of log2-transformed miRNA expression
#'   (miRNAs x patients). Rownames must be in precursor format (hsa-mir-21).
#' @param rna_sym Numeric matrix of log2-transformed RNA-seq expression
#'   (gene symbols x patients).
#' @param mirna_targets Data frame from \code{build_network()} with columns
#'   mirna, targets (list), n_targets.
#' @param mirna_norm_map Data frame with columns original (precursor ID) and
#'   norm (normalized precursor ID), used to map isoforms.
#' @return Data frame with columns: mirna, vss, n_targets, sorted descending by vss
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
#' score_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)
score_vss <- function(mirna_log, rna_sym, mirna_targets, mirna_norm_map) {

  # Build miRNA expression list (average across isoforms)
  mirna_expr <- .build_mirna_expr(mirna_log, mirna_targets$mirna, mirna_norm_map)

  n       <- nrow(mirna_targets)
  scores  <- numeric(n)
  n_tested <- integer(n)

  message("Computing VSS for ", n, " miRNAs...")

  for (i in seq_len(n)) {
    mid     <- mirna_targets$mirna[i]
    targets <- mirna_targets$targets[[i]]

    if (!mid %in% names(mirna_expr)) next
    x <- mirna_expr[[mid]]

    targets_present <- targets[targets %in% rownames(rna_sym)]
    if (length(targets_present) < 5) next

    target_mat <- rna_sym[targets_present, , drop = FALSE]

    r2_vals <- apply(target_mat, 1, function(y) {
      tryCatch(summary(lm(y ~ x))$r.squared, error = function(e) NA)
    })

    scores[i]    <- mean(r2_vals, na.rm = TRUE)
    n_tested[i]  <- length(targets_present)

    if (i %% 100 == 0) message("  ", i, " / ", n)
  }

  result <- data.frame(mirna     = mirna_targets$mirna,
                       vss       = scores,
                       n_targets = n_tested,
                       stringsAsFactors = FALSE)
  result <- result[result$vss > 0, ]
  result[order(-result$vss), ]
}
