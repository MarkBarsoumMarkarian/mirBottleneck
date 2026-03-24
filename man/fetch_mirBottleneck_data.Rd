\name{fetch_mirBottleneck_data}
\alias{fetch_mirBottleneck_data}
\title{Download the full mirBottleneck TCGA-PAAD dataset from Zenodo}
\usage{
fetch_mirBottleneck_data(
  dest = tools::R_user_dir("mirBottleneck", "cache"),
  overwrite = FALSE,
  checksum = TRUE
)
}
\arguments{
\item{dest}{Character. Directory to download files into. Defaults to a
persistent user cache directory.}
\item{overwrite}{Logical. Re-download even if files already exist (default FALSE).}
\item{checksum}{Logical. Verify checksums after download (default TRUE).}
}
\value{
Invisibly returns a named character vector of paths to downloaded files.
}
\description{
Downloads the complete harmonized TCGA-PAAD multi-omics dataset (~51 MB)
from Zenodo (DOI: 10.5281/zenodo.19121116) and caches it locally.
This function is optional and is never called during package checks or tests.
}
\examples{
\dontrun{
  paths <- fetch_mirBottleneck_data()
  rna   <- readRDS(paths[["rna_log"]])
  dim(rna)
}
}

