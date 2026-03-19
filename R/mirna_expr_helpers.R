#' Build per-miRNA expression vectors (averaging isoforms)
#'
#' Returns a named list where each element is a numeric vector across patients.
#'
#' @keywords internal
.build_mirna_expr <- function(mirna_log, norm_ids, mirna_norm_map) {
  if (!is.matrix(mirna_log)) stop("mirna_log must be a matrix.")
  if (!all(c("original", "norm") %in% colnames(mirna_norm_map))) {
    stop("mirna_norm_map must contain columns: original, norm")
  }

  norm_ids <- unique(norm_ids)
  out <- setNames(vector("list", length(norm_ids)), norm_ids)

  for (id in norm_ids) {
    original_ids <- mirna_norm_map$original[mirna_norm_map$norm == id]
    original_ids <- intersect(original_ids, rownames(mirna_log))

    if (length(original_ids) == 0) {
      out[[id]] <- NULL
      next
    }

    m <- mirna_log[original_ids, , drop = FALSE]
    x <- if (nrow(m) == 1) as.numeric(m[1, ]) else colMeans(m, na.rm = TRUE)
    names(x) <- colnames(mirna_log)

    out[[id]] <- x
  }

  out[!vapply(out, is.null, logical(1))]
}

#' Build miRNA expression matrix (rows=normalized miRNA IDs, cols=patients)
#'
#' @keywords internal
.build_mirna_expr_mat <- function(mirna_log, norm_ids, mirna_norm_map, patients) {
  if (!is.matrix(mirna_log)) stop("mirna_log must be a matrix.")
  if (!all(c("original", "norm") %in% colnames(mirna_norm_map))) {
    stop("mirna_norm_map must contain columns: original, norm")
  }

  norm_ids <- unique(norm_ids)
  patients <- intersect(patients, colnames(mirna_log))

  expr_mat <- matrix(
    NA_real_,
    nrow = length(norm_ids),
    ncol = length(patients),
    dimnames = list(norm_ids, patients)
  )

  for (i in seq_along(norm_ids)) {
    id <- norm_ids[i]
    original_ids <- mirna_norm_map$original[mirna_norm_map$norm == id]
    original_ids <- intersect(original_ids, rownames(mirna_log))
    if (length(original_ids) == 0) next

    m <- mirna_log[original_ids, patients, drop = FALSE]
    expr_mat[i, ] <- if (nrow(m) == 1) as.numeric(m[1, ]) else colMeans(m, na.rm = TRUE)
  }

  expr_mat
}