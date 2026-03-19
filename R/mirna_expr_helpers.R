#' @keywords internal
.build_mirna_expr <- function(mirna_log, norm_ids, mirna_norm_map) {
  if (!all(c("original", "norm") %in% colnames(mirna_norm_map))) {
    stop("Required columns missing in mirna_norm_map")
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

    # If only one isoform row, keep it as a vector; else average isoforms per patient
    out[[id]] <- if (nrow(m) == 1) as.numeric(m[1, ]) else colMeans(m, na.rm = TRUE)
    names(out[[id]]) <- colnames(mirna_log)
  }

  # drop NULL entries
  out[!vapply(out, is.null, logical(1))]
}

#' @keywords internal
.build_mirna_expr_mat <- function(mirna_log, norm_ids, mirna_norm_map, patients) {
  if (!all(c("original", "norm") %in% colnames(mirna_norm_map))) {
    stop("Required columns missing in mirna_norm_map")
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
