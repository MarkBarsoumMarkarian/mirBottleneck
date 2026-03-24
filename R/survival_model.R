#' Fit survival models for the composite bottleneck score
#'
#' Fits three Cox proportional hazards models (score alone, clinical alone,
#' combined) and a log-rank test comparing high vs low score groups.
#'
#' @param patient_scores Named numeric vector from \code{composite_score()}.
#' @param clinical_df Data frame with columns: patient, OS_days, OS_status,
#'   age, stage_clean.
#' @return List with:
#'   \itemize{
#'     \item \code{surv_df}: Patient-level data frame with scores and survival
#'     \item \code{cox_score}: Cox model, score alone
#'     \item \code{cox_clinical}: Cox model, clinical alone
#'     \item \code{cox_combined}: Cox model, score + clinical
#'     \item \code{logrank_p}: Log-rank p-value
#'     \item \code{summary}: Named vector of key metrics
#'   }
#' @export
#' @importFrom survival coxph Surv survfit survdiff
#' @importFrom stats pchisq
#' @examples
#' set.seed(1)
#' # patient_scores must be a *named numeric vector*
#' patient_scores <- rnorm(20)
#' names(patient_scores) <- paste0("P", 1:20)
#'
#' # clinical_df must contain: patient, OS_days, OS_status, age, stage_clean
#' clinical_df <- data.frame(
#'   patient = paste0("P", 1:20),
#'   OS_days = as.numeric(rexp(20, rate = 0.1) * 365),
#'   OS_status = sample(0:1, 20, replace = TRUE),
#'   age = sample(40:80, 20, replace = TRUE),
#'   stage_clean = sample(c("I", "II", "III", "IV"), 20, replace = TRUE)
#' )
#'
#' res <- survival_model(patient_scores, clinical_df)
#' names(res)
survival_model <- function(patient_scores, clinical_df) {

  surv_df <- clinical_df
  surv_df$composite <- patient_scores[match(surv_df$patient,
                                             names(patient_scores))]
  surv_df <- surv_df[!is.na(surv_df$composite) &
                     !is.na(surv_df$OS_days) &
                     surv_df$OS_days > 0, ]

  surv_df$score_group <- ifelse(
    surv_df$composite >= median(surv_df$composite, na.rm = TRUE),
    "High", "Low"
  )

  message("Patients in survival analysis: ", nrow(surv_df))

  cox_score    <- coxph(Surv(OS_days, OS_status) ~ composite,
                        data = surv_df)
  cox_clinical <- coxph(Surv(OS_days, OS_status) ~ age + stage_clean,
                        data = surv_df)
  cox_combined <- coxph(Surv(OS_days, OS_status) ~ composite + age + stage_clean,
                        data = surv_df)

  lr     <- survdiff(Surv(OS_days, OS_status) ~ score_group, data = surv_df)
  lr_p   <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  km_fit <- survfit(Surv(OS_days, OS_status) ~ score_group, data = surv_df)

  c1 <- summary(cox_score)$concordance[1]
  c2 <- summary(cox_clinical)$concordance[1]
  c3 <- summary(cox_combined)$concordance[1]

  message("\nSurvival results:")
  message("  C-index (score)   : ", round(c1, 4))
  message("  C-index (clinical): ", round(c2, 4))
  message("  C-index (combined): ", round(c3, 4))
  message("  AIC (score)       : ", round(AIC(cox_score),    2))
  message("  AIC (clinical)    : ", round(AIC(cox_clinical), 2))
  message("  AIC (combined)    : ", round(AIC(cox_combined), 2))
  message("  Log-rank p        : ", signif(lr_p, 3))

  list(
    surv_df      = surv_df,
    cox_score    = cox_score,
    cox_clinical = cox_clinical,
    cox_combined = cox_combined,
    km_fit       = km_fit,
    logrank_p    = lr_p,
    summary      = c(
      c_index_score    = c1,
      c_index_clinical = c2,
      c_index_combined = c3,
      aic_score        = AIC(cox_score),
      aic_clinical     = AIC(cox_clinical),
      aic_combined     = AIC(cox_combined),
      logrank_p        = lr_p
    )
  )
}
