% Generated manually
\name{compute_vss}
\alias{compute_vss}
\title{Compute Variance Suppression Scores (modular API)}
\usage{
compute_vss(
  mirna_log = NULL,
  rna_sym = NULL,
  mirna_targets,
  mirna_norm_map,
  se_mirna = NULL,
  se_mrna = NULL,
  mirna_assay = "mirna",
  mrna_assay = "mrna",
  BPPARAM = NULL
)
}
\arguments{
\item{mirna_log}{Numeric matrix/data.frame of miRNA expression (miRNAs x samples).}
\item{rna_sym}{Numeric matrix/data.frame of mRNA expression (genes x samples).}
\item{mirna_targets}{Data frame with columns \code{mirna}, \code{targets} (list), and \code{n_targets}.}
\item{mirna_norm_map}{Data frame with columns \code{original}, \code{norm}.}
\item{se_mirna}{Optional \code{SummarizedExperiment} containing the miRNA assay.}
\item{se_mrna}{Optional \code{SummarizedExperiment} containing the mRNA assay.}
\item{mirna_assay}{Assay name for miRNA expression when \code{se_mirna} is used.}
\item{mrna_assay}{Assay name for mRNA expression when \code{se_mrna} is used.}
\item{BPPARAM}{Optional BiocParallel parameter object. Currently reserved; computation falls back to serial execution.}
}
\value{Data frame with columns: \code{mirna}, \code{vss}, \code{n_targets}, ordered by decreasing \code{vss}.}
\description{Modular entry point for VSS computation supporting matrix/data.frame or \code{SummarizedExperiment} input.}
