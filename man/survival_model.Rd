\name{survival_model}
\alias{survival_model}
\title{Fit survival models for the composite bottleneck score}
\usage{
survival_model(patient_scores, clinical_df)
}
\arguments{
\item{patient_scores}{Named numeric vector from \code{composite_score()}.}
\item{clinical_df}{Data frame with columns: patient, OS_days, OS_status, age, stage_clean.}
}
\value{
List with: surv_df (patient-level data frame with scores and survival),
cox_score (Cox model, score alone), cox_clinical (Cox model, clinical alone),
cox_combined (Cox model, score + clinical), logrank_p (log-rank p-value),
summary (named vector of key metrics).
}
\description{
Fits three Cox proportional hazards models (score alone, clinical alone,
combined) and a log-rank test comparing high vs low score groups.
}
\examples{
patients <- paste0("P", 1:8)

patient_scores <- setNames(c(-1.2, -0.3, 0.2, 0.8, 1.1, -0.6, 0.4, 0.9), patients)

clinical_df <- data.frame(
  patient = patients,
  OS_days = c(100, 200, 150, 300, 400, 120, 250, 500),
  OS_status = c(1, 1, 0, 1, 0, 1, 0, 0),
  age = c(60, 55, 70, 50, 65, 62, 58, 49),
  stage_clean = c("II", "I", "III", "II", "IV", "II", "I", "III")
)

fit <- survival_model(patient_scores, clinical_df)
names(fit)
}
