#' @keywords internal
.build_mirna_expr <- function(mirna_log, norm_ids, mirna_norm_map) {
    if (!all(c('original', 'norm') %in% colnames(mirna_norm_map))) {
        stop('Required columns missing in mirna_norm_map')
    }

    # Match normalized IDs with original miRNA IDs
    avg_expr <- sapply(unique(norm_ids), function(id) {
        original_ids <- mirna_norm_map$original[mirna_norm_map$norm == id]
        if (length(original_ids) == 0) return(NA)
        rowMeans(mirna_log[rownames(mirna_log) %in% original_ids, , drop = FALSE], na.rm = TRUE)
    })

    names(avg_expr) <- unique(norm_ids)
    return(avg_expr)
}

#' @keywords internal
.build_mirna_expr_mat <- function(mirna_log, norm_ids, mirna_norm_map, patients) {
    if (!all(c('original', 'norm') %in% colnames(mirna_norm_map))) {
        stop('Required columns missing in mirna_norm_map')
    }

    expr_mat <- matrix(NA, nrow = length(unique(norm_ids)), ncol = length(patients), 
                       dimnames = list(unique(norm_ids), patients))

    for (i in seq_along(unique(norm_ids))) {
        original_ids <- mirna_norm_map$original[mirna_norm_map$norm == unique(norm_ids)[i]]
        expr_mat[i, ] <- rowMeans(mirna_log[rownames(mirna_log) %in% original_ids, patients, drop = FALSE], na.rm = TRUE)
    }
    return(expr_mat)
}