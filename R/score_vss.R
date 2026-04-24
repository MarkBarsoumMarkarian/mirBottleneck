#' Compute Variance Suppression Scores
#'
#' Backward-compatible wrapper around `compute_vss()`.
#'
#' @inheritParams compute_vss
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
  compute_vss(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map
  )
}
