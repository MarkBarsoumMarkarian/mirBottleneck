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
