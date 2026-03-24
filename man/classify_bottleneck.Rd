\name{classify_bottleneck}
\alias{classify_bottleneck}
\title{Classify miRNAs and compute combined bottleneck index}
\usage{
classify_bottleneck(vss_scores, coherence_scores)
}
\arguments{
\item{vss_scores}{Data frame from \code{score_vss()} with columns: mirna, vss, n_targets.}
\item{coherence_scores}{Data frame from \code{score_coherence()} with columns: mirna, coherence_score, coherence_p.}
}
\value{
Data frame with columns: mirna, vss, coherence_score, bottleneck_index, class.
Classes are: dual, silencer, conductor, weak.
}
\description{
Normalizes VSS and coherence scores to [0,1], computes a composite bottleneck
index, and classifies each miRNA into one of four archetypes based on median splits.
}
\examples{
vss_scores <- data.frame(
  mirna = c("hsa-mir-1", "hsa-mir-2", "hsa-mir-3", "hsa-mir-4"),
  vss = c(0.1, 0.8, 0.4, 0.6),
  n_targets = c(10, 12, 8, 15)
)

coh_scores <- data.frame(
  mirna = c("hsa-mir-1", "hsa-mir-2", "hsa-mir-3", "hsa-mir-4"),
  coherence_score = c(-0.2, 0.9, 0.1, 0.5),
  coherence_p = c(0.7, 0.01, 0.5, 0.04)
)

classified <- classify_bottleneck(vss_scores, coh_scores)
classified
}
