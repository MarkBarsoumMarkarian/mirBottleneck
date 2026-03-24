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
#' @examples
#' set.seed(1)
#' mirnas <- c("hsa-mir-1", "hsa-mir-2")
#' patients <- paste0("P", 1:10)
#' genes <- paste0("G", 1:6)
#'
#' mirna_log <- matrix(rnorm(length(mirnas) * length(patients)),
#'   nrow = length(mirnas),
#'   dimnames = list(mirnas, patients)
#' )
#' rna_sym <- matrix(rnorm(length(genes) * length(patients)),
#'   nrow = length(genes),
#'   dimnames = list(genes, patients)
#' )
#' mirna_targets <- data.frame(
#'   mirna = mirnas,
#'   targets = I(list(genes[1:3], genes[4:6])),
#'   n_targets = c(3, 3)
#' )
#' mirna_norm_map <- data.frame(original = mirnas, norm = mirnas)
#' score_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)
#' @export
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
