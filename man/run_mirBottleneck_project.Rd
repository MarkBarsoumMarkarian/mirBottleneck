\name{run_mirBottleneck_project}
\alias{run_mirBottleneck_project}
\title{Run the full mirBottleneck analysis pipeline}
\usage{
run_mirBottleneck_project(
  mirna_log_rds,
  rna_rds,
  clinical_rds,
  mirna_norm_map_rds,
  out_dir,
  mirna_targets_rds = NULL,
  query_network = TRUE,
  batch_size = 50,
  min_targets = 5,
  rescue_mirnas = c("hsa-miR-21-5p", "hsa-miR-21-3p", "hsa-miR-155-5p",
    "hsa-miR-155-3p", "hsa-miR-196a-5p", "hsa-miR-196b-5p", "hsa-miR-210-3p"),
  vss_keep_top = 300,
  n_perm = 200,
  max_targets = 200,
  cox_p_threshold = 0.15,
  report = FALSE,
  cohort_name = "TCGA-PAAD"
)
}
\arguments{
\item{mirna_log_rds}{Path to RDS file containing log2-transformed miRNA matrix (miRNAs x patients).}
\item{rna_rds}{Path to RDS file containing log2-transformed RNA-seq matrix (genes x patients).}
\item{clinical_rds}{Path to RDS file containing clinical data frame.}
\item{mirna_norm_map_rds}{Path to RDS file containing miRNA normalization map.}
\item{out_dir}{Output directory for results (created if absent).}
\item{mirna_targets_rds}{Optional path to pre-computed target interactions RDS.}
\item{query_network}{Logical. Whether to query miRTarBase (default TRUE).}
\item{batch_size}{Integer. Batch size for network queries (default 50).}
\item{min_targets}{Integer. Minimum targets per miRNA (default 5).}
\item{rescue_mirnas}{Character vector of mature miRNA IDs to rescue explicitly.}
\item{vss_keep_top}{Integer. Number of top VSS miRNAs to pass to coherence scoring (default 300).}
\item{n_perm}{Integer. Permutations for coherence p-values (default 200).}
\item{max_targets}{Integer. Max targets per miRNA for coherence scoring (default 200).}
\item{cox_p_threshold}{Numeric. p-value cutoff for composite score inclusion (default 0.15).}
\item{report}{Logical. Generate a self-contained HTML report after analysis (default FALSE).}
\item{cohort_name}{Character. Cohort label shown in the HTML report header (default "TCGA-PAAD").}
}
\value{
Invisibly returns a named list of all analysis outputs. Results are also
saved as RDS files to \code{out_dir}.
}
\description{
End-to-end wrapper that runs network construction, VSS scoring, coherence scoring,
miRNA classification, patient-level composite score, and survival modelling in sequence.
}
