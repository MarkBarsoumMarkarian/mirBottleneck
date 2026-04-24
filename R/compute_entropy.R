#' Compute transcriptome Shannon entropy per patient
#'
#' Converts gene expression to a probability distribution per patient and
#' computes Shannon entropy. High entropy = more disordered/noisy transcriptome.
#' Tests the hypothesis that high bottleneck load drives transcriptome entropy.
#'
#' @param rna_log Numeric matrix (genes x patients), log2-transformed.
#' @param patient_scores Named numeric vector from composite_score().
#' @param clinical Data frame with patient, OS_days, OS_status.
#' @param base Numeric. Logarithm base for entropy (default 2, bits).
#'
#' @return List with:
#'   \itemize{
#'     \item \code{entropy_per_patient} - named numeric vector of entropy values
#'     \item \code{spearman}  - result of cor.test() between score and entropy
#'     \item \code{entropy_df} - data frame: patient, entropy, bottleneck_score,
#'                               OS_days, OS_status
#'   }
#' @export
#' @examples
#' set.seed(1)
#' rna_log <- matrix(abs(rnorm(200 * 30)), nrow = 200)
#' rownames(rna_log) <- paste0("GENE", 1:200)
#' colnames(rna_log) <- paste0("P", 1:30)
#' patient_scores <- rnorm(30); names(patient_scores) <- paste0("P", 1:30)
#' clinical <- data.frame(patient = paste0("P", 1:30),
#'                        OS_days = rexp(30)*365, OS_status = rbinom(30,1,0.5))
#' compute_entropy(rna_log, patient_scores, clinical)
compute_entropy <- function(rna_log, patient_scores, clinical, base = 2) {

  # Convert log2 back to linear space, compute probability distribution
  # per patient, then Shannon entropy
  linear_mat <- 2^rna_log

  entropy_vec <- apply(linear_mat, 2, function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0) return(NA_real_)
    p <- x / sum(x)
    -sum(p * log(p, base = base), na.rm = TRUE)
  })

  names(entropy_vec) <- colnames(rna_log)

  # Align with patient scores
  shared <- intersect(names(patient_scores), names(entropy_vec))
  scores_shared  <- patient_scores[shared]
  entropy_shared <- entropy_vec[shared]

  spearman <- tryCatch(
    cor.test(scores_shared, entropy_shared, method = "spearman",
             exact = FALSE),
    error = function(e) NULL
  )

  if (!is.null(spearman)) {
    message("Entropy ~ Bottleneck Score: Spearman rho = ",
            round(spearman$estimate, 4),
            ", p = ", signif(spearman$p.value, 3))
  }

  # Build data frame
  entropy_df <- data.frame(
    patient          = shared,
    entropy          = entropy_shared,
    bottleneck_score = scores_shared,
    stringsAsFactors = FALSE
  )

  if (!is.null(clinical)) {
    clin_sub <- clinical[clinical$patient %in% shared,
                         c("patient", "OS_days", "OS_status")]
    entropy_df <- merge(entropy_df, clin_sub, by = "patient", all.x = TRUE)
  }

  list(
    entropy_per_patient = entropy_vec,
    spearman            = spearman,
    entropy_df          = entropy_df
  )
}


#' Compute TF activity via DoRothEA + VIPER
#'
#' Estimates transcription factor activity scores per patient using the
#' DoRothEA regulon database and the VIPER algorithm. Tests whether
#' TF activity dysregulation correlates with bottleneck score, providing
#' a mechanistic link between miRNA bottleneck load and transcriptome entropy.
#'
#' @param rna_log Numeric matrix (genes x patients), log2-transformed.
#' @param patient_scores Named numeric vector from composite_score().
#' @param confidence_levels Character vector. DoRothEA confidence levels to
#'   include (default c("A","B","C")).
#' @param organism Character. "human" (default) or "mouse".
#' @param top_n_tfs Integer. Report Spearman correlations for top N TFs by
#'   absolute correlation with bottleneck score (default 20).
#'
#' @return List with:
#'   \itemize{
#'     \item \code{tf_activity}  - matrix of TF activity scores (TFs x patients)
#'     \item \code{tf_cor}       - data frame of TF-score Spearman correlations
#'     \item \code{tf_entropy_cor} - Spearman: mean |TF activity| vs entropy
#'   }
#' @export
compute_tf_activity <- function(rna_log,
                                patient_scores,
                                confidence_levels = c("A", "B", "C"),
                                organism          = "human",
                                top_n_tfs         = 20) {

  if (!requireNamespace("decoupleR", quietly = TRUE))
    stop("Install decoupleR: BiocManager::install('decoupleR')")
  if (!requireNamespace("dorothea", quietly = TRUE))
    stop("Install dorothea: BiocManager::install('dorothea')")

  # Load DoRothEA regulons
  message("Loading DoRothEA regulons (confidence: ",
          paste(confidence_levels, collapse = ","), ")...")

  regulons <- dorothea::dorothea_hs  # human; swap to dorothea_mm for mouse
  regulons <- regulons[regulons$confidence %in% confidence_levels, ]

  # Run VIPER via decoupleR
  message("Running VIPER on ", ncol(rna_log), " patients...")
  tf_res <- decoupleR::run_viper(
    mat      = rna_log,
    network  = regulons,
    .source  = "tf",
    .target  = "target",
    .mor     = "mor",
    minsize  = 5,
    verbose  = FALSE
  )

  # Pivot to matrix: TFs x patients
  tf_mat <- tf_res |>
    dplyr::filter(statistic == "viper") |>
    tidyr::pivot_wider(id_cols = "source",
                       names_from = "condition",
                       values_from = "score") |>
    tibble::column_to_rownames("source") |>
    as.matrix()

  # Correlate each TF activity with bottleneck score
  shared <- intersect(colnames(tf_mat), names(patient_scores))
  tf_sub <- tf_mat[, shared, drop = FALSE]
  sc_sub <- patient_scores[shared]

  tf_cor <- apply(tf_sub, 1, function(act) {
    ct <- tryCatch(
      cor.test(act, sc_sub, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(ct)) return(c(rho = NA, p = NA))
    c(rho = ct$estimate, p = ct$p.value)
  })

  tf_cor_df <- as.data.frame(t(tf_cor))
  tf_cor_df$tf  <- rownames(tf_cor_df)
  tf_cor_df$p_adj <- p.adjust(tf_cor_df$p, method = "BH")
  tf_cor_df <- tf_cor_df[order(abs(tf_cor_df$rho), decreasing = TRUE), ]

  message("Top ", min(top_n_tfs, nrow(tf_cor_df)),
          " TFs correlated with bottleneck score:")
  print(head(tf_cor_df[, c("tf", "rho", "p", "p_adj")], top_n_tfs))

  list(
    tf_activity    = tf_mat,
    tf_cor         = tf_cor_df
  )
}


#' Immune deconvolution via TIMER2 REST API
#'
#' Calls the TIMER2.0 web API to estimate immune cell infiltration fractions
#' for TCGA samples, then tests their correlation with bottleneck score.
#' No local compute required; TIMER2 does all the heavy lifting server-side.
#'
#' @param cancer_type Character. TCGA cancer type abbreviation, e.g. "PAAD".
#'   Do NOT include the "TCGA-" prefix here.
#' @param patient_barcodes Character vector of TCGA patient barcodes.
#' @param patient_scores Named numeric vector from composite_score().
#' @param method Character. Deconvolution method (default "TIMER").
#'   Options: "TIMER", "CIBERSORT", "CIBERSORT-ABS", "EPIC", "MCP-counter".
#'
#' @return List with:
#'   \itemize{
#'     \item \code{immune_df}   - data frame of immune fractions per patient
#'     \item \code{immune_cor}  - Spearman correlation of each cell type vs score
#'   }
#' @export
compute_immune_deconvolution <- function(cancer_type,
                                         patient_barcodes,
                                         patient_scores,
                                         method = "TIMER") {

  base_url <- "http://timer.cistrome.org"

  # TIMER2 expects the cancer type without TCGA- prefix
  ct_short <- gsub("^TCGA-", "", toupper(cancer_type))

  message("Querying TIMER2 for ", ct_short, " (method: ", method, ")...")

  url <- paste0(base_url, "/infiltration_estimation_for_tcga.zip?",
                "cancers[]=", ct_short,
                "&methods[]=", tolower(method))

  tmp <- tempfile(fileext = ".zip")
  tryCatch(
    utils::download.file(url, tmp, quiet = TRUE, method = "curl"),
    error = function(e) stop("TIMER2 download failed: ", e$message)
  )

  # Unzip and read
  ex_dir <- tempfile()
  utils::unzip(tmp, exdir = ex_dir)
  csv_files <- list.files(ex_dir, pattern = "\\.csv$", full.names = TRUE,
                          recursive = TRUE)

  if (length(csv_files) == 0)
    stop("TIMER2 response contained no CSV files.")

  immune_raw <- read.csv(csv_files[1], check.names = FALSE)

  # Harmonize barcodes
  barcode_col <- grep("sample|barcode|tumor", colnames(immune_raw),
                      ignore.case = TRUE, value = TRUE)[1]
  immune_raw$patient <- harmonize_barcode(immune_raw[[barcode_col]])
  immune_raw <- immune_raw[immune_raw$patient %in% patient_barcodes, ]

  # Cell type columns (exclude metadata cols)
  meta_cols   <- c(barcode_col, "patient", "cancer_type", "purity")
  cell_cols   <- setdiff(colnames(immune_raw), meta_cols)
  cell_cols   <- cell_cols[grepl("B_cell|CD4|CD8|Macrophage|Neutrophil|NK|DC|Treg",
                                 cell_cols, ignore.case = TRUE)]

  immune_df <- immune_raw[, c("patient", cell_cols), drop = FALSE]

  # Correlate each cell type with bottleneck score
  shared <- intersect(immune_df$patient, names(patient_scores))
  sc_sub <- patient_scores[shared]

  immune_cor <- lapply(cell_cols, function(ct_col) {
    vals <- setNames(immune_df[[ct_col]], immune_df$patient)[shared]
    ct_r <- tryCatch(
      cor.test(vals, sc_sub, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(ct_r)) return(NULL)
    data.frame(cell_type = ct_col,
               rho = ct_r$estimate,
               p   = ct_r$p.value,
               stringsAsFactors = FALSE)
  })

  immune_cor_df <- dplyr::bind_rows(immune_cor)
  immune_cor_df$p_adj <- p.adjust(immune_cor_df$p, method = "BH")
  immune_cor_df <- immune_cor_df[order(immune_cor_df$p), ]

  message("Immune deconvolution complete. Top correlations with bottleneck score:")
  print(head(immune_cor_df, 10))

  list(
    immune_df  = immune_df,
    immune_cor = immune_cor_df
  )
}


#' Full entropy mechanism analysis
#'
#' Convenience wrapper that runs entropy computation, TF activity, and
#' immune deconvolution together and returns a unified mechanism report.
#'
#' @param rna_log Matrix (genes x patients).
#' @param patient_scores Named numeric vector.
#' @param clinical Data frame with patient, OS_days, OS_status.
#' @param cancer_type Character. TCGA project code (for TIMER2).
#' @param run_tf Logical. Run TF activity analysis (requires decoupleR +
#'   dorothea, default TRUE).
#' @param run_immune Logical. Run TIMER2 immune deconvolution (requires
#'   internet, default TRUE).
#'
#' @return List combining entropy, TF, and immune outputs.
#' @export
entropy_mechanism <- function(rna_log,
                              patient_scores,
                              clinical,
                              cancer_type,
                              run_tf     = TRUE,
                              run_immune = TRUE) {

  message("=== Entropy Mechanism Analysis ===")

  # Entropy
  entropy_res <- compute_entropy(rna_log, patient_scores, clinical)

  # TF activity
  tf_res <- NULL
  if (run_tf) {
    tf_res <- tryCatch(
      compute_tf_activity(rna_log, patient_scores),
      error = function(e) {
        warning("TF activity failed: ", e$message); NULL
      }
    )
  }

  # Immune deconvolution
  immune_res <- NULL
  if (run_immune) {
    immune_res <- tryCatch(
      compute_immune_deconvolution(
        cancer_type      = cancer_type,
        patient_barcodes = clinical$patient,
        patient_scores   = patient_scores
      ),
      error = function(e) {
        warning("Immune deconvolution failed: ", e$message); NULL
      }
    )
  }

  # Mediation: does entropy mediate bottleneck -> survival?
  mediation_res <- NULL
  if (!is.null(entropy_res$entropy_df) &&
      "OS_days" %in% colnames(entropy_res$entropy_df)) {
    mediation_res <- tryCatch(
      .run_mediation(entropy_res$entropy_df),
      error = function(e) {
        warning("Mediation analysis failed: ", e$message); NULL
      }
    )
  }

  list(
    entropy    = entropy_res,
    tf         = tf_res,
    immune     = immune_res,
    mediation  = mediation_res
  )
}


#' @keywords internal
.run_mediation <- function(entropy_df) {
  # Simple mediation: bottleneck_score -> entropy -> OS (Cox)
  # Using Baron-Kenny approach with log-rank on entropy terciles
  df <- entropy_df[complete.cases(entropy_df), ]

  if (nrow(df) < 30) return(NULL)

  # Path a: bottleneck -> entropy
  path_a <- lm(entropy ~ bottleneck_score, data = df)

  # Path b + c': bottleneck + entropy -> OS
  path_bc <- tryCatch(
    survival::coxph(survival::Surv(OS_days, OS_status) ~
                      bottleneck_score + entropy,
                    data = df),
    error = function(e) NULL
  )

  # Path c: bottleneck -> OS (without mediator)
  path_c <- tryCatch(
    survival::coxph(survival::Surv(OS_days, OS_status) ~
                      bottleneck_score,
                    data = df),
    error = function(e) NULL
  )

  list(
    path_a_summary  = summary(path_a),
    path_bc_summary = if (!is.null(path_bc)) summary(path_bc) else NULL,
    path_c_summary  = if (!is.null(path_c))  summary(path_c)  else NULL,
    note = paste(
      "Baron-Kenny mediation: path_a = bottleneck->entropy,",
      "path_c = bottleneck->OS (total),",
      "path_bc = bottleneck+entropy->OS (direct).",
      "Mediation supported if entropy coef significant in path_bc",
      "and bottleneck coef attenuates vs path_c."
    )
  )
}
