\name{harmonize_barcode}
\alias{harmonize_barcode}
\title{Harmonize TCGA barcodes to patient level}
\usage{
harmonize_barcode(x)
}
\arguments{
\item{x}{Character vector of TCGA barcodes.}
}
\value{
Character vector of 12-character patient-level barcodes.
}
\description{
Standardizes any TCGA barcode to the first 12 characters (patient-level),
removing sample type, vial, portion, and analyte suffixes.
}
\examples{
harmonize_barcode("TCGA-3A-A9I7-01A-21R-A38N-13")
}
