\name{composite_score}
\alias{composite_score}
\title{Compute patient-level composite bottleneck score}
\usage{
composite_score(
  mirna_log,
  classified_df,
  clinical_df,
  mirna_norm_map,
  cox_p_threshold = 0.15
)
}
\arguments{
\item{mirna_log}{Numeric matrix of log2-transformed miRNA expression (miRNAs x patients).}
\item{classified_df}{Data frame from \code{classify_bottleneck()}.}
\item{clinical_df}{Data frame with columns: patient, OS_days, OS_status, age, stage_clean.}
\item{mirna_norm_map}{Data frame mapping precursor IDs to normalized IDs.}
\item{cox_p_threshold}{Numeric. Nominal p-value cutoff for miRNA inclusion (default 0.15).}
}
\value{
List with: patient_scores (named numeric vector), contributing_mirnas (data frame),
cox_individual (full individual Cox results).
}
\description{
Builds a direction-aware weighted composite score per patient using the
miRNAs most associated with survival. Weights are proportional to |Cox coefficient|.
}
