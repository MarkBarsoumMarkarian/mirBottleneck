# Internal input validation and extraction helpers

#' Validate and extract matrix inputs
#'
#' @param mirna_log Numeric matrix/data.frame (miRNAs x patients)
#' @param rna_sym Numeric matrix/data.frame (genes x patients)
#' @param mirna_targets Data frame with columns mirna, targets, n_targets
#' @param mirna_norm_map Data frame with columns original, norm
#' @param scored_mirnas Optional character vector of miRNA IDs
#' @keywords internal
.validate_inputs_matrix <- function(mirna_log,
                                    rna_sym,
                                    mirna_targets,
                                    mirna_norm_map,
                                    scored_mirnas = NULL) {
  if (!is.matrix(mirna_log) && !is.data.frame(mirna_log)) {
    stop("mirna_log must be a matrix or data.frame.", call. = FALSE)
  }
  if (!is.matrix(rna_sym) && !is.data.frame(rna_sym)) {
    stop("rna_sym must be a matrix or data.frame.", call. = FALSE)
  }

  mirna_log <- as.matrix(mirna_log)
  rna_sym <- as.matrix(rna_sym)

  if (!is.numeric(mirna_log)) stop("mirna_log must be numeric.", call. = FALSE)
  if (!is.numeric(rna_sym)) stop("rna_sym must be numeric.", call. = FALSE)

  if (is.null(rownames(mirna_log)) || anyNA(rownames(mirna_log))) {
    stop("mirna_log must have non-missing rownames (miRNA IDs).", call. = FALSE)
  }
  if (is.null(rownames(rna_sym)) || anyNA(rownames(rna_sym))) {
    stop("rna_sym must have non-missing rownames (gene symbols).", call. = FALSE)
  }

  if (is.null(colnames(mirna_log)) || is.null(colnames(rna_sym))) {
    stop("mirna_log and rna_sym must have colnames (sample IDs).", call. = FALSE)
  }

  mir_samples <- colnames(mirna_log)
  mrna_samples <- colnames(rna_sym)
  if (!setequal(mir_samples, mrna_samples)) {
    stop(
      "Sample IDs in mirna_log and rna_sym do not match. ",
      "Provide aligned matrices or a SummarizedExperiment with matching assay columns.",
      call. = FALSE
    )
  }
  if (!identical(mir_samples, mrna_samples)) {
    message("Reordering rna_sym columns to match mirna_log sample order.")
    rna_sym <- rna_sym[, mir_samples, drop = FALSE]
  }

  if (!is.data.frame(mirna_targets)) {
    stop("mirna_targets must be a data.frame.", call. = FALSE)
  }
  if (!all(c("mirna", "targets", "n_targets") %in% colnames(mirna_targets))) {
    stop("mirna_targets must contain columns: mirna, targets, n_targets.", call. = FALSE)
  }

  if (!is.data.frame(mirna_norm_map) ||
      !all(c("original", "norm") %in% colnames(mirna_norm_map))) {
    stop("mirna_norm_map must contain columns: original, norm.", call. = FALSE)
  }

  if (!is.null(scored_mirnas) && !is.character(scored_mirnas)) {
    stop("scored_mirnas must be a character vector.", call. = FALSE)
  }

  list(
    mirna_log = mirna_log,
    rna_sym = rna_sym
  )
}

#' Extract miRNA and mRNA assays from SummarizedExperiment
#'
#' @param se SummarizedExperiment instance
#' @param mirna_assay Character. miRNA assay name.
#' @param mrna_assay Character. mRNA assay name.
#' @keywords internal
.extract_from_se <- function(se, mirna_assay = "mirna", mrna_assay = "mrna") {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object.", call. = FALSE)
  }

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required for se input.", call. = FALSE)
  }

  assay_names <- SummarizedExperiment::assayNames(se)
  if (!mirna_assay %in% assay_names) {
    stop(
      "Missing miRNA assay '", mirna_assay, "'. Available assays: ",
      paste(assay_names, collapse = ", "),
      call. = FALSE
    )
  }
  if (!mrna_assay %in% assay_names) {
    stop(
      "Missing mRNA assay '", mrna_assay, "'. Available assays: ",
      paste(assay_names, collapse = ", "),
      call. = FALSE
    )
  }

  mirna_log <- SummarizedExperiment::assay(se, mirna_assay)
  rna_sym <- SummarizedExperiment::assay(se, mrna_assay)

  if (is.null(colnames(mirna_log)) || is.null(colnames(rna_sym))) {
    stop("Both assays must have sample IDs as colnames.", call. = FALSE)
  }
  if (!identical(colnames(mirna_log), colnames(rna_sym))) {
    stop(
      "Sample alignment mismatch between assays '", mirna_assay,
      "' and '", mrna_assay, "'. Ensure identical sample IDs and order.",
      call. = FALSE
    )
  }

  list(
    mirna_log = mirna_log,
    rna_sym = rna_sym
  )
}

#' Resolve expression inputs from matrices or SummarizedExperiment
#'
#' @keywords internal
.resolve_expression_inputs <- function(mirna_log = NULL,
                                       rna_sym = NULL,
                                       se = NULL,
                                       mirna_assay = "mirna",
                                       mrna_assay = "mrna",
                                       mirna_targets,
                                       mirna_norm_map,
                                       scored_mirnas = NULL) {
  if (!is.null(se)) {
    extracted <- .extract_from_se(
      se = se,
      mirna_assay = mirna_assay,
      mrna_assay = mrna_assay
    )
    mirna_log <- extracted$mirna_log
    rna_sym <- extracted$rna_sym
  }

  validated <- .validate_inputs_matrix(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map,
    scored_mirnas = scored_mirnas
  )

  validated
}
