% Generated manually
\name{plot_archetype_landscape}
\alias{plot_archetype_landscape}
\title{Plot miRNA archetype landscape}
\usage{
plot_archetype_landscape(classified_df)
}
\arguments{
\item{classified_df}{Data frame from \code{classify_archetypes()} or \code{classify_bottleneck()}. Must contain columns \code{vss}, \code{coherence_score}, and \code{class}.}
}
\value{A \code{ggplot2} object.}
\description{Creates a VSS-vs-coherence scatter plot colored by archetype class and returns a customizable \code{ggplot2} object.}
