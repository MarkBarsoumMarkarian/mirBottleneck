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
\examples{
set.seed(1)

mirnas <- c("hsa-mir-1", "hsa-mir-2", "hsa-mir-3")
patients <- paste0("P", 1:6)

mirna_log <- matrix(rnorm(length(mirnas) * length(patients)),
                    nrow = length(mirnas),
                    dimnames = list(mirnas, patients))

classified_df <- data.frame(
  mirna = mirnas,
  vss = c(0.2, 0.7, 0.4),
  coherence_score = c(0.1, 0.5, -0.2),
  bottleneck_index = c(0.3, 0.9, 0.2),
  class = c("weak", "dual", "silencer")
)

clinical_df <- data.frame(
  patient = patients,
  OS_days = c(200, 450, 120, 800, 500, 300),
  OS_status = c(1, 0, 1, 0, 1, 1),
  age = c(60, 55, 70, 50, 65, 62),
  stage_clean = c("II", "I", "III", "II", "IV", "II")
)

mirna_norm_map <- data.frame(
  original = mirnas,
  norm = mirnas
)

out <- composite_score(
  mirna_log = mirna_log,
  classified_df = classified_df,
  clinical_df = clinical_df,
  mirna_norm_map = mirna_norm_map,
  cox_p_threshold = 1
)

str(out)
}
