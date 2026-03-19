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
