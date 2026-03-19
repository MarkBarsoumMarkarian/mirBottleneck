#' Harmonize TCGA barcodes to patient level
#'
#' Standardizes any TCGA barcode to the first 12 characters (patient-level),
#' removing sample type, vial, portion, and analyte suffixes.
#'
#' @param x Character vector of TCGA barcodes
#' @return Character vector of 12-character patient-level barcodes
#' @export
#' @examples
#' harmonize_barcode("TCGA-3A-A9I7-01A-21R-A38N-13")
#' # Returns "TCGA-3A-A9I7"
harmonize_barcode <- function(x) {
  substr(trimws(x), 1, 12)
}

#' Normalize miRNA ID to lowercase precursor format
#'
#' Converts mature miRNA IDs (hsa-miR-21-5p) to lowercase precursor format
#' (hsa-mir-21) by lowercasing and stripping arm suffixes.
#'
#' @param x Character vector of miRNA IDs
#' @return Normalized character vector
#' @keywords internal
normalize_mirna <- function(x) {
  x <- tolower(x)
  x <- gsub("-5p$|-3p$", "", x)
  x
}

#' Normalize a numeric vector to 0-1 range
#'
#' @param x Numeric vector
#' @return Numeric vector scaled to [0, 1]; returns all-zero for zero-variance
#'   input, and all-NA for non-finite input.
#' @keywords internal
normalize_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(rep(NA_real_, length(x)))
  if (rng[2] - rng[1] == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

#' Compute mean pairwise correlation of a matrix
#'
#' @param mat Numeric matrix (genes x samples)
#' @return Mean of upper triangle of correlation matrix, or NA
#' @keywords internal
mean_pairwise_cor <- function(mat) {
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) return(NA)
  if (any(dim(mat) == 0)) return(NA)
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[row_vars > 0, , drop = FALSE]
  if (nrow(mat) < 2) return(NA)
  cor_mat <- cor(t(mat), use = "pairwise.complete.obs")
  upper <- cor_mat[upper.tri(cor_mat)]
  if (length(upper) == 0) return(NA)
  mean(upper, na.rm = TRUE)
}
