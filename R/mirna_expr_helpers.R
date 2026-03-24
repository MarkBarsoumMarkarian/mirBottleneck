# Internal helpers for miRNA expression averaging

#' Build a named list of per-patient miRNA expression vectors
#'
#' For each unique normalized miRNA ID, averages expression across all
#' isoform rows in \code{mirna_log}, then returns a named list where each
#' element is a numeric vector of length = number of patients.
#'
#' @param mirna_log Numeric matrix (miRNAs x patients).
#' @param norm_ids Character vector of normalized IDs, parallel to rownames(mirna_log).
#' @param mirna_norm_map Data frame with columns \code{original} and \code{norm}.
#' @return Named list of numeric vectors (one per unique normalized miRNA).
#' @examples
#' set.seed(1)
#' mirna_log <- matrix(rnorm(12), nrow = 3,
#'   dimnames = list(c("hsa-mir-1","hsa-mir-2","hsa-mir-3"), paste0("P", 1:4))
#' )
#' mirna_norm_map <- data.frame(
#'   original = rownames(mirna_log),
#'   norm = rownames(mirna_log),
#'   stringsAsFactors = FALSE
#' )
#' x <- .build_mirna_expr(mirna_log, norm_ids = rownames(mirna_log), mirna_norm_map = mirna_norm_map)
#' names(x)
#' @keywords internal
.build_mirna_expr <- function(mirna_log, norm_ids, mirna_norm_map) {

  if (missing(norm_ids) || is.null(norm_ids)) {
    idx      <- match(rownames(mirna_log), mirna_norm_map$original)
    norm_ids <- ifelse(is.na(idx), rownames(mirna_log), mirna_norm_map$norm[idx])
  }

  unique_ids <- unique(norm_ids)

  result <- lapply(unique_ids, function(mid) {
    rows <- which(norm_ids == mid)
    if (length(rows) == 1L) {
      mirna_log[rows, ]
    } else {
      colMeans(mirna_log[rows, , drop = FALSE])
    }
  })

  names(result) <- unique_ids
  result
}

#' Build an averaged miRNA expression matrix
#'
#' Collapses isoform rows by averaging and optionally subsets to specific
#' patients, returning a matrix with one row per unique normalized miRNA ID.
#'
#' @param mirna_log Numeric matrix (miRNAs x patients).
#' @param norm_ids Character vector of normalized IDs, parallel to rownames(mirna_log).
#' @param mirna_norm_map Data frame with columns \code{original} and \code{norm}.
#' @param patients Optional character vector of patient IDs to subset columns.
#' @return Numeric matrix (unique miRNAs x patients).
#' @examples
#' set.seed(1)
#' mirna_log <- matrix(rnorm(12), nrow = 3,
#'   dimnames = list(c("hsa-mir-1","hsa-mir-2","hsa-mir-3"), paste0("P", 1:4))
#' )
#' mirna_norm_map <- data.frame(
#'   original = rownames(mirna_log),
#'   norm = rownames(mirna_log),
#'   stringsAsFactors = FALSE
#' )
#' mat <- .build_mirna_expr_mat(
#'   mirna_log,
#'   norm_ids = rownames(mirna_log),
#'   mirna_norm_map = mirna_norm_map,
#'   patients = colnames(mirna_log)
#' )
#' dim(mat)
#' @keywords internal
.build_mirna_expr_mat <- function(mirna_log, norm_ids, mirna_norm_map,
                                  patients = NULL) {
  expr_list <- .build_mirna_expr(mirna_log, norm_ids, mirna_norm_map)
  mat       <- do.call(rbind, expr_list)

  if (!is.null(patients)) {
    patients <- intersect(patients, colnames(mat))
    mat      <- mat[, patients, drop = FALSE]
  }

  mat
}
