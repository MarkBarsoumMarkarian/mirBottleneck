#' Paths to packaged toy datasets
#'
#' Returns a named list of absolute paths to the small synthetic datasets
#' shipped with mirBottleneck. These are used in examples, tests, and
#' vignettes so that no network access or large downloads are required.
#'
#' The toy datasets contain 10 synthetic patients, 50 genes, and 30 miRNAs,
#' and a precomputed interaction network — enough to demonstrate the full
#' pipeline deterministically in seconds.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{mirna_log}{Log2-transformed miRNA expression matrix (30 x 10).}
#'     \item{rna_sym}{Log2-transformed RNA-seq matrix with gene symbols (50 x 10).}
#'     \item{clinical_df}{Clinical data frame (10 x 9).}
#'     \item{mirna_norm_map}{miRNA precursor ID normalization map.}
#'     \item{mirna_targets}{Precomputed miRNA-target interaction network.}
#'   }
#' @keywords internal
#' @examples
#' paths <- mirBottleneck:::.toy_paths()
#' mirna_log <- readRDS(paths$mirna_log)
#' dim(mirna_log) # 30 x 10
.toy_paths <- function() {
  p <- function(f) system.file("extdata", f, package = "mirBottleneck",
                               mustWork = TRUE)
  list(
    mirna_log      = p("toy_mirna_log.rds"),
    rna_sym        = p("toy_rna_sym.rds"),
    clinical_df    = p("toy_clinical_df.rds"),
    mirna_norm_map = p("toy_mirna_norm_map.rds"),
    mirna_targets  = p("toy_mirna_targets.rds")
  )
}
