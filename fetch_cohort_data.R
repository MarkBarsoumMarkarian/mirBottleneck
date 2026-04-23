#' Fetch and preprocess TCGA miRNA + RNA + clinical data for any cancer type
#'
#' Downloads miRNA-seq, RNA-seq, and clinical data from TCGA via TCGAbiolinks,
#' preprocesses them into the matrix format expected by mirBottleneck functions,
#' and returns a list ready to pipe into run_mirBottleneck().
#'
#' @param cancer_type Character. TCGA project code, e.g. "TCGA-PAAD", "TCGA-BRCA".
#' @param save_dir Character. Directory to cache downloaded GDC files.
#'   Defaults to a subdirectory named after the cancer type inside tempdir().
#' @param min_os_days Numeric. Minimum OS_days to retain a patient (default 1).
#' @param log2_pseudocount Numeric. Pseudocount before log2 transform (default 1).
#' @param min_mirna_cpm Numeric. Minimum CPM for miRNA filtering (default 1).
#' @param min_gene_cpm Numeric. Minimum CPM for RNA filtering (default 1).
#' @param min_mirna_prev Numeric. Minimum sample prevalence for miRNA (default 0.2).
#' @param min_gene_prev Numeric. Minimum sample prevalence for RNA (default 0.2).
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{mirna_log}   - log2 miRNA matrix (miRNAs x patients)
#'     \item \code{rna_log}     - log2 RNA matrix (genes x patients)
#'     \item \code{clinical}    - data frame: patient, OS_days, OS_status, age, stage_clean
#'     \item \code{mirna_norm_map} - data frame: original, norm
#'     \item \code{cancer_type} - the input cancer_type string
#'     \item \code{n_patients}  - number of patients with complete data
#'   }
#'
#' @export
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare getTCGATable
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom dplyr select mutate filter rename left_join coalesce
#' @importFrom stringr str_extract str_to_upper
#'
#' @examples
#' \donttest{
#' cohort <- fetch_cohort_data("TCGA-LUAD", save_dir = "/tmp/tcga_cache")
#' str(cohort)
#' }
fetch_cohort_data <- function(cancer_type,
                              save_dir        = file.path(tempdir(), cancer_type),
                              min_os_days     = 1,
                              log2_pseudocount = 1,
                              min_mirna_cpm   = 1,
                              min_gene_cpm    = 1,
                              min_mirna_prev  = 0.2,
                              min_gene_prev   = 0.2,
                              verbose         = TRUE) {

  .vcat <- function(...) if (verbose) message(...)

  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------ #
  # 1. Download miRNA-seq                                               #
  # ------------------------------------------------------------------ #
  .vcat("[", cancer_type, "] Querying miRNA-seq...")

  q_mirna <- TCGAbiolinks::GDCquery(
    project           = cancer_type,
    data.category     = "Transcriptome Profiling",
    data.type         = "miRNA Expression Quantification",
    workflow.type     = "BCGSC miRNA Profiling",
    legacy            = FALSE
  )
  TCGAbiolinks::GDCdownload(q_mirna, directory = save_dir, method = "api")
  se_mirna <- TCGAbiolinks::GDCprepare(q_mirna, directory = save_dir)

  mirna_raw <- SummarizedExperiment::assay(se_mirna, "reads_per_million_miRNA_mapped")
  mirna_barcodes <- harmonize_barcode(colnames(mirna_raw))
  colnames(mirna_raw) <- mirna_barcodes

  # Keep tumor samples only (barcode digit 14-15 == "01")
  tumor_idx <- grepl("-01$|-01[A-Z]", colnames(se_mirna))
  mirna_raw <- mirna_raw[, tumor_idx, drop = FALSE]
  colnames(mirna_raw) <- harmonize_barcode(colnames(mirna_raw))

  # Deduplicate patients (keep first occurrence)
  mirna_raw <- mirna_raw[, !duplicated(colnames(mirna_raw)), drop = FALSE]

  .vcat("[", cancer_type, "] miRNA matrix: ", nrow(mirna_raw), " x ", ncol(mirna_raw))

  # ------------------------------------------------------------------ #
  # 2. Download RNA-seq                                                 #
  # ------------------------------------------------------------------ #
  .vcat("[", cancer_type, "] Querying RNA-seq...")

  q_rna <- TCGAbiolinks::GDCquery(
    project           = cancer_type,
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts",
    legacy            = FALSE
  )
  TCGAbiolinks::GDCdownload(q_rna, directory = save_dir, method = "api")
  se_rna <- TCGAbiolinks::GDCprepare(q_rna, directory = save_dir)

  rna_raw <- SummarizedExperiment::assay(se_rna, "unstranded")
  tumor_idx_rna <- grepl("-01$|-01[A-Z]", colnames(se_rna))
  rna_raw <- rna_raw[, tumor_idx_rna, drop = FALSE]
  colnames(rna_raw) <- harmonize_barcode(colnames(rna_raw))
  rna_raw <- rna_raw[, !duplicated(colnames(rna_raw)), drop = FALSE]

  # Map ENSEMBL to gene symbols
  gene_info <- as.data.frame(SummarizedExperiment::rowData(se_rna))
  symbol_col <- intersect(c("gene_name", "external_gene_name", "gene_id"), colnames(gene_info))[1]
  gene_symbols <- gene_info[[symbol_col]]
  names(gene_symbols) <- rownames(rna_raw)

  # Remove duplicated symbols and NAs
  valid <- !is.na(gene_symbols) & gene_symbols != "" & !duplicated(gene_symbols)
  rna_raw <- rna_raw[valid, , drop = FALSE]
  rownames(rna_raw) <- gene_symbols[valid]

  .vcat("[", cancer_type, "] RNA matrix: ", nrow(rna_raw), " x ", ncol(rna_raw))

  # ------------------------------------------------------------------ #
  # 3. Clinical data                                                    #
  # ------------------------------------------------------------------ #
  .vcat("[", cancer_type, "] Fetching clinical data...")

  clin_raw <- TCGAbiolinks::GDCquery_clinic(cancer_type, type = "clinical")

  clinical <- .parse_clinical(clin_raw, cancer_type)

  # ------------------------------------------------------------------ #
  # 4. Find shared patients                                             #
  # ------------------------------------------------------------------ #
  shared <- Reduce(intersect, list(
    colnames(mirna_raw),
    colnames(rna_raw),
    clinical$patient
  ))

  clinical <- clinical[clinical$patient %in% shared, ]
  clinical <- clinical[!is.na(clinical$OS_days) & clinical$OS_days >= min_os_days, ]
  shared   <- intersect(shared, clinical$patient)

  mirna_raw <- mirna_raw[, shared, drop = FALSE]
  rna_raw   <- rna_raw[,   shared, drop = FALSE]

  .vcat("[", cancer_type, "] Shared patients with complete data: ", length(shared))

  if (length(shared) < 20) {
    warning(cancer_type, ": only ", length(shared),
            " shared patients. Results may be unreliable.")
  }

  # ------------------------------------------------------------------ #
  # 5. Filter + log2 transform                                          #
  # ------------------------------------------------------------------ #
  mirna_log <- .filter_and_log(mirna_raw, min_cpm = min_mirna_cpm,
                                min_prev = min_mirna_prev,
                                pseudocount = log2_pseudocount,
                                label = "miRNA")

  rna_log <- .filter_and_log(rna_raw, min_cpm = min_gene_cpm,
                              min_prev = min_gene_prev,
                              pseudocount = log2_pseudocount,
                              label = "RNA")

  # ------------------------------------------------------------------ #
  # 6. Build mirna_norm_map                                             #
  # ------------------------------------------------------------------ #
  mirna_norm_map <- data.frame(
    original = rownames(mirna_log),
    norm     = normalize_mirna(rownames(mirna_log)),
    stringsAsFactors = FALSE
  )

  .vcat("[", cancer_type, "] Done. miRNA: ", nrow(mirna_log),
        " | RNA: ", nrow(rna_log), " | patients: ", ncol(mirna_log))

  list(
    mirna_log      = mirna_log,
    rna_log        = rna_log,
    clinical       = clinical,
    mirna_norm_map = mirna_norm_map,
    cancer_type    = cancer_type,
    n_patients     = length(shared)
  )
}


# ------------------------------------------------------------------ #
# Internal helpers                                                    #
# ------------------------------------------------------------------ #

#' @keywords internal
.parse_clinical <- function(clin_raw, cancer_type) {

  clin <- as.data.frame(clin_raw)

  # Patient barcode
  id_col <- intersect(c("submitter_id", "bcr_patient_barcode", "case_submitter_id"),
                      colnames(clin))[1]
  clin$patient <- harmonize_barcode(clin[[id_col]])

  # OS days: prefer days_to_death, fall back to days_to_last_follow_up
  os_days <- suppressWarnings(as.numeric(clin$days_to_death))
  fu_days <- suppressWarnings(as.numeric(clin$days_to_last_follow_up))
  clin$OS_days <- ifelse(!is.na(os_days), os_days, fu_days)

  # OS status: 1 = dead, 0 = alive/censored
  clin$OS_status <- as.integer(
    grepl("dead|deceased", tolower(as.character(clin$vital_status)))
  )

  # Age
  age_col <- intersect(c("age_at_index", "age_at_diagnosis", "bcr_patient_age"),
                       colnames(clin))[1]
  clin$age <- suppressWarnings(as.numeric(clin[[age_col]]))

  # Stage: normalize to a clean factor
  stage_cols <- intersect(c("ajcc_pathologic_stage", "tumor_stage",
                             "clinical_stage", "pathologic_stage"), colnames(clin))
  if (length(stage_cols) > 0) {
    raw_stage <- as.character(clin[[stage_cols[1]]])
    clin$stage_clean <- .normalize_stage(raw_stage)
  } else {
    clin$stage_clean <- NA_character_
  }

  clin[, c("patient", "OS_days", "OS_status", "age", "stage_clean")]
}

#' @keywords internal
.normalize_stage <- function(x) {
  x <- toupper(trimws(x))
  x[grepl("^STAGE\\s*I[^IV]|^STAGE\\s*I$|^I$|^1$", x)] <- "I"
  x[grepl("^STAGE\\s*II[^I]|^STAGE\\s*II$|^II$|^2$", x)] <- "II"
  x[grepl("^STAGE\\s*III[^I]|^STAGE\\s*III$|^III$|^3$", x)] <- "III"
  x[grepl("^STAGE\\s*IV|^IV$|^4$", x)] <- "IV"
  x[!x %in% c("I", "II", "III", "IV")] <- NA_character_
  x
}

#' @keywords internal
.filter_and_log <- function(mat, min_cpm, min_prev, pseudocount, label) {
  lib_sizes <- colSums(mat, na.rm = TRUE)
  cpm_mat   <- t(t(mat) / lib_sizes * 1e6)

  # Filter by CPM threshold in at least min_prev fraction of samples
  keep <- rowSums(cpm_mat >= min_cpm, na.rm = TRUE) >= (min_prev * ncol(mat))
  mat  <- mat[keep, , drop = FALSE]
  message("  ", label, ": kept ", sum(keep), " / ", length(keep),
          " features after CPM filter")

  log2(mat + pseudocount)
}

#' @keywords internal
.build_mirna_expr_mat <- function(mirna_log, norm_ids, mirna_norm_map, patients) {
  # Map norm IDs back to row names in mirna_log
  lookup <- setNames(mirna_norm_map$original, mirna_norm_map$norm)
  raw_ids <- lookup[norm_ids]

  present <- raw_ids[raw_ids %in% rownames(mirna_log)]
  mat <- mirna_log[present, patients, drop = FALSE]
  rownames(mat) <- names(present)
  mat
}
