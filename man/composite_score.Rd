#' Compute patient-level composite bottleneck score
#'
#' Builds a direction-aware weighted composite score per patient using the
#' miRNAs most associated with survival. miRNAs with HR > 1 push the score
#' up; protective miRNAs push it down. Weights are proportional to the
#' absolute Cox coefficient.
#'
#' @param mirna_log Numeric matrix of log2-transformed miRNA expression
#'   (miRNAs x patients).
#' @param classified_df Data frame from \code{classify_mirnas()}.
#' @param clinical_df Data frame with columns: patient, OS_days, OS_status,
#'   age, stage_clean.
#' @param mirna_norm_map Data frame mapping precursor IDs to normalized IDs.
#' @param cox_p_threshold Numeric. Nominal p-value cutoff for miRNA inclusion
#'   in composite score (default 0.15).
#' @return List with:
#'   \itemize{
#'     \item \code{patient_scores}: Named numeric vector of composite scores
#'     \item \code{contributing_mirnas}: Data frame of miRNAs used in score
#'     \item \code{cox_individual}: Full individual Cox results
#'   }
#' @examples
#' set.seed(1)
#' mirnas <- c("hsa-mir-1", "hsa-mir-2", "hsa-mir-3")
#' patients <- paste0("P", 1:8)
#'
#' mirna_log <- matrix(
#'   rnorm(length(mirnas) * length(patients)),
#'   nrow = length(mirnas),
#'   dimnames = list(mirnas, patients)
#' )
#'
#' classified_df <- data.frame(
#'   mirna = mirnas,
#'   class = c("weak", "dual", "silencer"),
#'   stringsAsFactors = FALSE
#' )
#'
#' clinical_df <- data.frame(
#'   patient = patients,
#'   OS_days = c(100, 200, 150, 300, 400, 120, 250, 500),
#'   OS_status = c(1, 1, 0, 1, 0, 1, 0, 0),
#'   age = c(60, 55, 70, 50, 65, 62, 58, 49),
#'   stage_clean = c("II", "I", "III", "II", "IV", "II", "I", "III"),
#'   stringsAsFactors = FALSE
#' )
#'
#' mirna_norm_map <- data.frame(original = mirnas, norm = mirnas, stringsAsFactors = FALSE)
#'
#' out <- composite_score(
#'   mirna_log = mirna_log,
#'   classified_df = classified_df,
#'   clinical_df = clinical_df,
#'   mirna_norm_map = mirna_norm_map,
#'   cox_p_threshold = 1
#' )
#' names(out)
#' @export
#' @importFrom survival coxph Surv
#' @importFrom dplyr filter mutate arrange
composite_score <- function(mirna_log, classified_df, clinical_df,
                            mirna_norm_map, cox_p_threshold = 0.15) {

  shared_patients <- clinical_df$patient
  mirna_expr_mat  <- .build_mirna_expr_mat(mirna_log,
                                           norm_ids     = classified_df$mirna,
                                           mirna_norm_map = mirna_norm_map,
                                           patients     = shared_patients)

  base_df <- clinical_df[!is.na(clinical_df$OS_days) & clinical_df$OS_days > 0, ]

  # Individual Cox models adjusted for age and stage
  message("Running ", nrow(classified_df), " individual Cox models...")

  cox_results <- lapply(classified_df$mirna, function(mid) {
    expr <- mirna_expr_mat[mid, base_df$patient]
    if (all(is.na(expr))) return(NULL)
    if (var(expr, na.rm = TRUE) == 0) return(NULL)

    df <- base_df
    df$mirna_expr <- as.numeric(expr)

    tryCatch({
      fit <- coxph(Surv(OS_days, OS_status) ~ mirna_expr + age + stage_clean,
                   data = df)
      s   <- summary(fit)
      data.frame(
        mirna     = mid,
        coef      = s$coefficients["mirna_expr", "coef"],
        hr        = s$coefficients["mirna_expr", "exp(coef)"],
        se        = s$coefficients["mirna_expr", "se(coef)"],
        z         = s$coefficients["mirna_expr", "z"],
        p         = s$coefficients["mirna_expr", "Pr(>|z|)"],
        c_index   = s$concordance[1],
        class     = classified_df$class[classified_df$mirna == mid],
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  })

  cox_df <- dplyr::bind_rows(cox_results)
  cox_df$p_adj <- p.adjust(cox_df$p, method = "BH")
  cox_df <- cox_df[order(cox_df$p), ]

  # Select contributing miRNAs
  contributing <- cox_df[cox_df$p < cox_p_threshold, ]
  message(nrow(contributing), " miRNAs pass p < ", cox_p_threshold)

  # Build composite score
  sel_mat <- mirna_expr_mat[contributing$mirna, shared_patients, drop = FALSE]

  patient_scores <- apply(sel_mat, 2, function(expr_vec) {
    weights  <- abs(contributing$coef)
    directed <- expr_vec * sign(contributing$coef)
    weighted.mean(directed, w = weights, na.rm = TRUE)
  })

  list(
    patient_scores      = patient_scores,
    contributing_mirnas = contributing,
    cox_individual      = cox_df
  )
}
