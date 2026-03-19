# Internal helpers for miRNA expression averaging

#' Build a named list of per-patient miRNA expression vectors
#'
#' For each unique normalized miRNA ID, averages expression across all
#' isoform rows in \code{mirna_log}, then returns a named list where each
#' element is a numeric vector of length = number of patients.
#'
#' @param mirna_log Numeric matrix (miRNAs x patients). Rownames are precursor
#'   IDs in original (possibly isoform-suffixed) format.
#' @param norm_ids Character vector of normalized IDs, same length as
#'   \code{nrow(mirna_log)}, parallel to rownames(mirna_log).
#' @param mirna_norm_map Data frame with columns \code{original} and
#'   \code{norm} mapping raw precursor IDs to normalized IDs.
#' @return Named list; each element is a numeric vector (length = ncol(mirna_log))
#'   representing averaged expression for that normalized miRNA across patients.
#' @keywords internal
.build_mirna_expr <- function(mirna_log, norm_ids, mirna_norm_map) {

  # Map rownames to normalized IDs via mirna_norm_map if norm_ids not supplied
  if (missing(norm_ids) || is.null(norm_ids)) {
    idx      <- match(rownames(mirna_log), mirna_norm_map$original)
    norm_ids <- ifelse(is.na(idx), rownames(mirna_log),
                       mirna_norm_map$norm[idx])
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
#' Collapses isoform rows by averaging, returning a matrix with one row per
#' unique normalized miRNA ID.
#'
#' @inheritParams .build_mirna_expr
#' @return Numeric matrix (unique miRNAs x patients).
#' @keywords internal
.build_mirna_expr_mat <- function(mirna_log, norm_ids, mirna_norm_map) {
  expr_list <- .build_mirna_expr(mirna_log, norm_ids, mirna_norm_map)
  do.call(rbind, expr_list)
}
