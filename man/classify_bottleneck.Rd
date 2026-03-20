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
