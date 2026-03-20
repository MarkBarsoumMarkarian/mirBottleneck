\name{score_coherence}
\alias{score_coherence}
\title{Compute Coherence Induction Scores}
\usage{
score_coherence(
  mirna_log,
  rna_sym,
  mirna_targets,
  mirna_norm_map,
  scored_mirnas,
  n_perm = 200,
  max_targets = 200
)
}
\arguments{
\item{mirna_log}{Numeric matrix of log2-transformed miRNA expression (miRNAs x patients).}
\item{rna_sym}{Numeric matrix of log2-transformed RNA-seq expression (gene symbols x patients).}
\item{mirna_targets}{Data frame from \code{build_network()}.}
\item{mirna_norm_map}{Data frame mapping precursor IDs to normalized IDs.}
\item{scored_mirnas}{Character vector of miRNA IDs to score.}
\item{n_perm}{Integer. Number of permutations for p-value (default 200).}
\item{max_targets}{Integer. Maximum targets per miRNA for speed (default 200).}
}
\value{
Data frame with columns: mirna, coherence_score, coherence_p.
}
\description{
Splits patients by median miRNA expression and measures whether target genes
become more correlated in the high-expression group. The coherence score is
the difference in mean pairwise target correlation between high and low groups.
An empirical p-value is computed via permutation.
}
