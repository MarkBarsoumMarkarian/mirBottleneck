\name{score_vss}
\alias{score_vss}
\title{Compute Variance Suppression Scores}
\usage{
score_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)
}
\arguments{
\item{mirna_log}{Numeric matrix of log2-transformed miRNA expression (miRNAs x patients).
Rownames must be in precursor format (hsa-mir-21).}
\item{rna_sym}{Numeric matrix of log2-transformed RNA-seq expression (gene symbols x patients).}
\item{mirna_targets}{Data frame from \code{build_network()} with columns mirna, targets (list), n_targets.}
\item{mirna_norm_map}{Data frame with columns original (precursor ID) and norm (normalized precursor ID).}
}
\value{
Data frame with columns: mirna, vss, n_targets, sorted descending by vss.
}
\description{
For each miRNA, fits a linear model of each target gene's expression on
the miRNA's expression across patients. The VSS is the mean R-squared
across all validated targets, quantifying how much transcriptome variance
the miRNA explains.
}
\examples{
set.seed(1)

mirnas <- c("hsa-mir-1", "hsa-mir-2")
patients <- paste0("P", 1:8)
genes <- c("G1", "G2", "G3", "G4")

mirna_log <- matrix(rnorm(length(mirnas) * length(patients)),
                    nrow = length(mirnas),
                    dimnames = list(mirnas, patients))

rna_sym <- matrix(rnorm(length(genes) * length(patients)),
                  nrow = length(genes),
                  dimnames = list(genes, patients))

mirna_targets <- data.frame(
  mirna = mirnas,
  targets = I(list(c("G1", "G2"), c("G3", "G4"))),
  n_targets = c(2, 2)
)

mirna_norm_map <- data.frame(original = mirnas, norm = mirnas)

vss <- score_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)
vss
}
