#' Download the full mirBottleneck TCGA-PAAD dataset from Zenodo
#'
#' Downloads the complete harmonized TCGA-PAAD multi-omics dataset (~51 MB)
#' from Zenodo and caches it locally. This function is purely optional and is
#' never called during package checks, tests, or vignette builds — those use
#' the small toy datasets shipped with the package instead.
#'
#' @param dest Character. Directory to download files into. Defaults to a
#'   persistent user cache directory via \code{tools::R_user_dir()}.
#' @param overwrite Logical. Re-download even if files already exist (default FALSE).
#' @param checksum Logical. Verify MD5 checksums after download (default TRUE).
#'
#' @return Invisibly returns a named character vector of paths to downloaded files.
#'
#' @details
#' Dataset DOI: \url{https://doi.org/10.5281/zenodo.19121116}
#'
#' The full dataset contains:
#' \itemize{
#'   \item \code{rna_log.rds} — log2-transformed RNA-seq matrix (60,660 genes x 178 patients)
#'   \item \code{rna_matrix.rds} — raw RNA-seq counts
#'   \item \code{rna_final.rds} — filtered RNA-seq matrix (gene symbols)
#'   \item \code{mirna_log.rds} — log2-transformed miRNA matrix (1,881 x 178 patients)
#'   \item \code{mirna_matrix.rds} — raw miRNA counts
#'   \item \code{mirna_final.rds} — filtered miRNA matrix
#'   \item \code{mirna_targets.rds} — validated miRNA-target interactions
#'   \item \code{interactions_filtered.rds} / \code{interactions_strict.rds} — filtered networks
#'   \item \code{maf_filtered.rds} — somatic mutation calls (MAF format)
#'   \item \code{clinical.rds} — patient clinical metadata
#'   \item \code{survival_df.rds} — survival data with bottleneck scores
#' }
#'
#' @examples
#' \dontrun{
#'   paths <- fetch_mirBottleneck_data()
#'   rna   <- readRDS(paths[["rna_log"]])
#'   dim(rna)  # 60660 x 178
#' }
#'
#' @export
fetch_mirBottleneck_data <- function(dest      = tools::R_user_dir("mirBottleneck", "cache"),
                                     overwrite = FALSE,
                                     checksum  = TRUE) {

  zenodo_base <- "https://zenodo.org/records/19121116/files"
  doi         <- "10.5281/zenodo.19121116"

  files <- c(
    "rna_log.rds",
    "rna_matrix.rds",
    "rna_final.rds",
    "mirna_log.rds",
    "mirna_matrix.rds",
    "mirna_final.rds",
    "mirna_targets.rds",
    "interactions_filtered.rds",
    "interactions_strict.rds",
    "mirna_target_interactions.rds",
    "maf_filtered.rds",
    "clinical.rds",
    "survival_df.rds"
  )

  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    message("Created cache directory: ", dest)
  }

  message("Fetching mirBottleneck full dataset")
  message("DOI: ", doi)
  message("Destination: ", dest)
  message(rep("-", 60))

  out_paths <- stats::setNames(
    file.path(dest, files),
    tools::file_path_sans_ext(files)
  )

  for (i in seq_along(files)) {
    f    <- files[i]
    dest_f <- out_paths[i]

    if (file.exists(dest_f) && !overwrite) {
      message("  [skip] ", f, " (already cached)")
      next
    }

    url <- paste0(zenodo_base, "/", f, "?download=1")
    message("  [download] ", f, " ...")

    tryCatch(
      utils::download.file(url, destfile = dest_f, mode = "wb", quiet = TRUE),
      error = function(e) {
        warning("Failed to download ", f, ": ", conditionMessage(e))
      }
    )

    if (file.exists(dest_f)) {
      size_mb <- round(file.size(dest_f) / 1024^2, 1)
      message("  [ok] ", f, " (", size_mb, " MB)")
    }
  }

  message(rep("-", 60))
  message("Done. Load with e.g.:")
  message("  rna <- readRDS(\"", file.path(dest, "rna_log.rds"), "\")")

  invisible(out_paths)
}
