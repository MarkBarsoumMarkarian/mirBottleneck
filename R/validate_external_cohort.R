#' Validate mirBottleneck score in an external cohort
#'
#' Applies the contributing miRNAs and their Cox weights derived from TCGA-PAAD
#' to a new cohort (e.g. ICGC-PACA-AU or a GEO dataset), computes the composite
#' bottleneck score, and evaluates survival performance. No refitting is done:
#' this is a true locked-model external validation.
#'
#' @param ext_mirna_log Numeric matrix. log2 miRNA expression for the external
#'   cohort (miRNAs x patients). Row names must be in the same ID format as
#'   the training cohort (precursor or mature, will be normalized).
#' @param ext_clinical Data frame with columns: patient, OS_days, OS_status.
#'   age and stage_clean are optional but recommended.
#' @param contributing_mirnas Data frame from the PAAD composite_score() output
#'   (\code{composite$contributing_mirnas}). Must contain: mirna, coef, hr.
#' @param cohort_name Character. Label for this cohort in output (default "External").
#' @param mirna_norm_map Optional data frame mapping original to normalized miRNA IDs
#'   in the external cohort. If NULL, normalize_mirna() is applied directly.
#'
#' @return List with:
#'   \itemize{
#'     \item \code{patient_scores} - named numeric vector of composite scores
#'     \item \code{survival}       - output of survival_model()
#'     \item \code{n_mirnas_found} - how many training miRNAs were in ext cohort
#'     \item \code{n_mirnas_total} - how many training miRNAs were expected
#'     \item \code{cohort_name}
#'   }
#' @export
#' @examples
#' \donttest{
#' # Load PAAD result
#' paad <- readRDS("mirBottleneck_results/TCGA-PAAD/result.rds")
#' # Load external cohort matrices
#' icgc <- fetch_icgc_paca()
#' val  <- validate_external_cohort(
#'   ext_mirna_log      = icgc$mirna_log,
#'   ext_clinical       = icgc$clinical,
#'   contributing_mirnas = paad$composite$contributing_mirnas,
#'   cohort_name        = "ICGC-PACA-AU"
#' )
#' val$survival$summary
#' }
validate_external_cohort <- function(ext_mirna_log,
                                     ext_clinical,
                                     contributing_mirnas,
                                     cohort_name   = "External",
                                     mirna_norm_map = NULL) {

  message("[", cohort_name, "] External validation: locked-model scoring")

  # Normalize external miRNA IDs
  ext_norm <- normalize_mirna(rownames(ext_mirna_log))
  rownames(ext_mirna_log) <- ext_norm

  # Normalize training miRNA IDs
  train_mirnas <- normalize_mirna(contributing_mirnas$mirna)

  # Find overlap
  found <- intersect(train_mirnas, ext_norm)
  missing <- setdiff(train_mirnas, ext_norm)

  message("[", cohort_name, "] Training miRNAs found: ",
          length(found), " / ", length(train_mirnas))

  if (length(missing) > 0) {
    message("[", cohort_name, "] Missing: ", paste(missing, collapse = ", "))
  }

  if (length(found) == 0) {
    stop("No training miRNAs found in external cohort. Check ID format.")
  }

  # Subset contributing_mirnas to those found
  cont_sub <- contributing_mirnas
  cont_sub$mirna_norm <- train_mirnas
  cont_sub <- cont_sub[cont_sub$mirna_norm %in% found, ]

  # Align patients
  shared_pts <- intersect(colnames(ext_mirna_log), ext_clinical$patient)
  ext_mirna_sub <- ext_mirna_log[cont_sub$mirna_norm, shared_pts, drop = FALSE]

  # Compute composite score (same locked weights from training)
  patient_scores <- apply(ext_mirna_sub, 2, function(expr_vec) {
    weights  <- abs(cont_sub$coef)
    directed <- expr_vec * sign(cont_sub$coef)
    stats::weighted.mean(directed, w = weights, na.rm = TRUE)
  })

  message("[", cohort_name, "] Patients scored: ", length(patient_scores))

  # Survival analysis
  clin_sub <- ext_clinical[ext_clinical$patient %in% shared_pts, ]

  # Add dummy age/stage if missing (survival_model needs them for full model)
  if (!"age" %in% colnames(clin_sub)) clin_sub$age <- NA_real_
  if (!"stage_clean" %in% colnames(clin_sub)) clin_sub$stage_clean <- NA_character_

  surv_result <- tryCatch(
    survival_model(patient_scores, clin_sub),
    error = function(e) {
      warning("[", cohort_name, "] survival_model failed: ", e$message)
      list(summary = c(c_index_score = NA, logrank_p = NA,
                       c_index_combined = NA))
    }
  )

  message("[", cohort_name, "] C-index: ",
          round(surv_result$summary["c_index_score"], 4),
          " | Log-rank p: ",
          signif(surv_result$summary["logrank_p"], 3))

  list(
    patient_scores  = patient_scores,
    survival        = surv_result,
    n_mirnas_found  = length(found),
    n_mirnas_total  = length(train_mirnas),
    cohort_name     = cohort_name
  )
}


#' Fetch ICGC-PACA-AU miRNA + clinical data
#'
#' Downloads ICGC PACA-AU donor expression and clinical data from the ICGC
#' Data Portal. Requires prior ICGC registration and DACO approval for
#' controlled-access files. Open-access miRNA data is available without login.
#'
#' @param data_dir Character. Directory to store downloaded files (default tempdir()).
#' @param use_cached Logical. If TRUE and files exist in data_dir, skip download.
#'
#' @return List with:
#'   \itemize{
#'     \item \code{mirna_log}  - log2 miRNA matrix (miRNAs x patients)
#'     \item \code{clinical}   - data frame: patient, OS_days, OS_status, age
#'     \item \code{n_patients}
#'   }
#' @export
fetch_icgc_paca <- function(data_dir   = file.path(tempdir(), "ICGC_PACA"),
                             use_cached = TRUE) {

  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  # ICGC open-access miRNA expression endpoint
  mirna_url  <- paste0(
    "https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/",
    "PACA-AU/miRNA-Seq/open/exp_mirna.PACA-AU.tsv.gz"
  )
  donor_url  <- paste0(
    "https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/",
    "PACA-AU/donor.PACA-AU.tsv.gz"
  )

  mirna_file <- file.path(data_dir, "exp_mirna.PACA-AU.tsv.gz")
  donor_file <- file.path(data_dir, "donor.PACA-AU.tsv.gz")

  if (!use_cached || !file.exists(mirna_file)) {
    message("Downloading ICGC PACA-AU miRNA expression...")
    tryCatch(
      utils::download.file(mirna_url, mirna_file, quiet = FALSE, method = "curl"),
      error = function(e) stop("ICGC download failed: ", e$message,
                               "\nIf access is restricted, download manually from ",
                               "https://dcc.icgc.org/projects/PACA-AU")
    )
  }

  if (!use_cached || !file.exists(donor_file)) {
    message("Downloading ICGC PACA-AU donor clinical data...")
    utils::download.file(donor_url, donor_file, quiet = FALSE, method = "curl")
  }

  # Parse miRNA expression
  message("Parsing miRNA expression...")
  mirna_raw <- read.delim(gzfile(mirna_file), check.names = FALSE,
                           stringsAsFactors = FALSE)

  # ICGC format: rows = donors, cols = miRNA IDs (or transposed)
  # Detect orientation
  if ("icgc_donor_id" %in% colnames(mirna_raw) ||
      "donor_id" %in% colnames(mirna_raw)) {
    # Rows are donors, pivot to miRNA x patient
    id_col <- grep("donor_id", colnames(mirna_raw), value = TRUE)[1]
    mir_col <- grep("mirna_id|mir_id|mirbase_id", colnames(mirna_raw),
                    ignore.case = TRUE, value = TRUE)[1]
    val_col <- grep("normalized|rpm|read_count|expression",
                    colnames(mirna_raw), ignore.case = TRUE, value = TRUE)[1]

    mirna_wide <- tidyr::pivot_wider(
      mirna_raw[, c(id_col, mir_col, val_col)],
      names_from  = id_col,
      values_from = val_col,
      id_cols     = mir_col
    )
    rownames_mat <- mirna_wide[[mir_col]]
    mirna_mat <- as.matrix(mirna_wide[, -1])
    rownames(mirna_mat) <- rownames_mat
  } else {
    # Assume miRNAs as rows already
    mirna_mat <- as.matrix(mirna_raw[, -1])
    rownames(mirna_mat) <- mirna_raw[[1]]
  }

  storage.mode(mirna_mat) <- "numeric"

  # Log2 transform
  mirna_log <- log2(mirna_mat + 1)

  # Parse clinical
  message("Parsing clinical data...")
  donor_raw <- read.delim(gzfile(donor_file), check.names = FALSE,
                           stringsAsFactors = FALSE)

  patient_col <- grep("donor_id|icgc_donor", colnames(donor_raw), value = TRUE)[1]
  donor_raw$patient <- donor_raw[[patient_col]]

  # Survival: ICGC uses donor_survival_time (days) and donor_vital_status
  os_days <- suppressWarnings(as.numeric(donor_raw$donor_survival_time))
  fu_days <- suppressWarnings(as.numeric(donor_raw$donor_interval_of_last_followup))
  donor_raw$OS_days <- ifelse(!is.na(os_days), os_days, fu_days)
  donor_raw$OS_status <- as.integer(
    grepl("deceased|dead", tolower(donor_raw$donor_vital_status))
  )

  age_col <- grep("donor_age|age_at", colnames(donor_raw),
                  ignore.case = TRUE, value = TRUE)[1]
  donor_raw$age <- suppressWarnings(as.numeric(donor_raw[[age_col]]))

  clinical <- donor_raw[, c("patient", "OS_days", "OS_status", "age")]
  clinical <- clinical[!is.na(clinical$OS_days) & clinical$OS_days > 0, ]

  # Shared patients
  shared <- intersect(colnames(mirna_log), clinical$patient)
  mirna_log <- mirna_log[, shared, drop = FALSE]
  clinical  <- clinical[clinical$patient %in% shared, ]

  message("ICGC-PACA-AU: ", length(shared), " patients with complete data")

  list(
    mirna_log  = mirna_log,
    clinical   = clinical,
    n_patients = length(shared)
  )
}


#' Fetch a GEO dataset as mirBottleneck-compatible matrices
#'
#' Downloads a GEO series (GSE) containing miRNA expression and survival data,
#' parses it into the standard mirBottleneck matrix format.
#'
#' @param geo_id Character. GEO series accession, e.g. "GSE71729".
#' @param mirna_platform_pattern Character. Regex to match the miRNA GPL platform
#'   name (default "miRNA|mirna").
#' @param os_days_col Character. Column name in pData containing OS days.
#' @param os_status_col Character. Column name in pData containing OS status
#'   (1=event, 0=censored).
#' @param data_dir Character. Cache directory.
#'
#' @return List with mirna_log, clinical, n_patients.
#' @export
#' @importFrom GEOquery getGEO
fetch_geo_cohort <- function(geo_id,
                              mirna_platform_pattern = "miRNA|mirna",
                              os_days_col   = NULL,
                              os_status_col = NULL,
                              data_dir      = file.path(tempdir(), geo_id)) {

  if (!requireNamespace("GEOquery", quietly = TRUE))
    stop("Install GEOquery: BiocManager::install('GEOquery')")

  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  message("Downloading GEO series: ", geo_id, "...")
  gse <- GEOquery::getGEO(geo_id, destdir = data_dir, GSEMatrix = TRUE,
                           AnnotGPL = FALSE)

  # If multiple platforms, find the miRNA one
  if (length(gse) > 1) {
    platform_titles <- sapply(gse, function(g)
      GEOquery::Meta(g)$platform_title[1])
    mirna_idx <- grep(mirna_platform_pattern, platform_titles,
                      ignore.case = TRUE)[1]
    if (is.na(mirna_idx)) {
      warning("No miRNA platform found matching '", mirna_platform_pattern,
              "'. Using first platform.")
      mirna_idx <- 1L
    }
    gse <- gse[[mirna_idx]]
  } else {
    gse <- gse[[1]]
  }

  # Expression matrix
  expr_mat <- Biobase::exprs(gse)
  pdata    <- Biobase::pData(gse)
  fdata    <- Biobase::fData(gse)

  # Map probe IDs to miRNA names
  mir_col <- grep("miRNA|mirna|mir_id|accession",
                  colnames(fdata), ignore.case = TRUE, value = TRUE)[1]

  if (!is.na(mir_col)) {
    rownames(expr_mat) <- make.unique(as.character(fdata[[mir_col]]))
  }

  # Normalize IDs
  rownames(expr_mat) <- normalize_mirna(rownames(expr_mat))

  # Ensure log2 scale (GEO data is often already log2)
  if (max(expr_mat, na.rm = TRUE) > 50) {
    message("Values appear linear, applying log2(x+1) transform.")
    expr_mat <- log2(expr_mat + 1)
  }

  # Clinical data from pData
  pdata$patient <- rownames(pdata)

  # Auto-detect OS columns if not provided
  if (is.null(os_days_col)) {
    os_days_col <- grep("os.*day|surv.*day|days.*surv|overall.*surv.*day|time",
                        colnames(pdata), ignore.case = TRUE, value = TRUE)[1]
    message("Auto-detected OS days column: ", os_days_col)
  }
  if (is.null(os_status_col)) {
    os_status_col <- grep("os.*stat|vital|dead|event|censor",
                          colnames(pdata), ignore.case = TRUE, value = TRUE)[1]
    message("Auto-detected OS status column: ", os_status_col)
  }

  if (is.na(os_days_col) || is.na(os_status_col)) {
    stop("Could not detect OS columns. Specify os_days_col and os_status_col manually.\n",
         "Available columns: ", paste(colnames(pdata), collapse = ", "))
  }

  os_days <- suppressWarnings(as.numeric(pdata[[os_days_col]]))
  os_status_raw <- pdata[[os_status_col]]

  # Parse OS status: handles 0/1, alive/dead, yes/no
  os_status <- suppressWarnings(as.integer(os_status_raw))
  if (anyNA(os_status)) {
    os_status <- as.integer(grepl("dead|deceased|1|yes|event",
                                  tolower(os_status_raw)))
  }

  clinical <- data.frame(
    patient   = pdata$patient,
    OS_days   = os_days,
    OS_status = os_status,
    stringsAsFactors = FALSE
  )
  clinical <- clinical[!is.na(clinical$OS_days) & clinical$OS_days > 0, ]

  shared <- intersect(colnames(expr_mat), clinical$patient)
  expr_mat <- expr_mat[, shared, drop = FALSE]
  clinical <- clinical[clinical$patient %in% shared, ]

  message(geo_id, ": ", length(shared), " patients with complete miRNA + OS data")

  list(
    mirna_log  = expr_mat,
    clinical   = clinical,
    n_patients = length(shared)
  )
}


#' Meta-analysis forest plot across multiple validation cohorts
#'
#' Takes validation results from validate_external_cohort() and the TCGA-PAAD
#' training result, fits a random-effects meta-analysis on log(HR), and
#' produces a publication-ready forest plot.
#'
#' @param results_list Named list. Each element is output of
#'   validate_external_cohort() or run_mirBottleneck(). Names become cohort labels.
#' @param output_file Character. Path to save the PDF forest plot (default NULL,
#'   just returns the plot).
#' @param title Character. Plot title.
#'
#' @return Invisibly returns the meta-analysis result (meta::metagen object).
#' @export
plot_meta_forest <- function(results_list,
                              output_file = NULL,
                              title       = "mirBottleneck Score: Multi-Cohort Survival Meta-Analysis") {

  if (!requireNamespace("meta", quietly = TRUE))
    stop("Install meta: install.packages('meta')")

  # Extract HR and SE from each cohort's Cox model
  cohort_rows <- lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]

    # Handle both run_mirBottleneck and validate_external_cohort output formats
    cox_model <- if (!is.null(res$survival$cox_score)) {
      res$survival$cox_score
    } else NULL

    if (is.null(cox_model)) {
      return(data.frame(cohort = nm, log_hr = NA, se = NA, n = NA,
                        stringsAsFactors = FALSE))
    }

    s     <- summary(cox_model)
    coef  <- s$coefficients[1, "coef"]
    se    <- s$coefficients[1, "se(coef)"]
    n_pat <- if (!is.null(res$n_patients)) res$n_patients else nrow(s$conf.int)

    data.frame(
      cohort = nm,
      log_hr = coef,
      se     = se,
      n      = n_pat,
      stringsAsFactors = FALSE
    )
  })

  meta_df <- dplyr::bind_rows(cohort_rows)
  meta_df <- meta_df[!is.na(meta_df$log_hr), ]

  if (nrow(meta_df) < 2) {
    warning("Need at least 2 cohorts for meta-analysis. Returning single cohort result.")
  }

  message("Running random-effects meta-analysis on ", nrow(meta_df), " cohorts...")
  ma <- meta::metagen(
    TE       = meta_df$log_hr,
    seTE     = meta_df$se,
    studlab  = meta_df$cohort,
    data     = meta_df,
    sm       = "HR",
    common   = FALSE,
    random   = TRUE,
    title    = title
  )

  message("Pooled HR: ", round(exp(ma$TE.random), 3),
          " [", round(exp(ma$lower.random), 3),
          "-", round(exp(ma$upper.random), 3), "]",
          " | I² = ", round(ma$I2 * 100, 1), "%",
          " | p = ", signif(ma$pval.random, 3))

  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = max(4, nrow(meta_df) + 3))
    meta::forest.meta(ma,
                      sortvar     = meta_df$log_hr,
                      print.tau2  = TRUE,
                      print.I2    = TRUE,
                      col.diamond = "#2166AC",
                      col.square  = "#4DAF4A",
                      xlab        = "Hazard Ratio (log scale)")
    dev.off()
    message("Forest plot saved to: ", output_file)
  } else {
    meta::forest.meta(ma,
                      sortvar     = meta_df$log_hr,
                      print.tau2  = TRUE,
                      print.I2    = TRUE,
                      col.diamond = "#2166AC",
                      col.square  = "#4DAF4A",
                      xlab        = "Hazard Ratio (log scale)")
  }

  invisible(ma)
}
