#' Build miRNA-target interaction network
#'
#' Retrieves validated miRNA-target interactions from miRTarBase via multiMiR,
#' filters to targets present in the RNA-seq dataset, and returns a target list
#' per miRNA. Queries are batched to avoid server timeouts.
#'
#' @param mirna_ids Character vector of miRNA IDs (precursor format, e.g. hsa-mir-21)
#' @param rna_symbols Character vector of gene symbols present in RNA-seq data
#' @param batch_size Integer. Number of miRNAs per multiMiR query (default 50)
#' @param min_targets Integer. Minimum validated targets to retain a miRNA (default 5)
#' @param rescue_mirnas Character vector of additional mature miRNA IDs to query
#'   explicitly (e.g. known oncomiRs missed by batch query). Default includes
#'   hsa-miR-21-5p, hsa-miR-155-5p.
#' @return A data frame with columns: mirna, targets (list column), n_targets
#' @export
#' @importFrom multiMiR get_multimir
#' @importFrom dplyr filter mutate select group_by summarise bind_rows distinct
#' @examples
#' \donttest{
#' # build_network() may query a remote miRNA-target resource; keep out of routine checks
#' if (requireNamespace("multiMiR", quietly = TRUE)) {
#'   mirna_ids <- c("hsa-miR-21-5p", "hsa-miR-155-5p")
#'   rna_symbols <- c("KRAS", "TP53", "SMAD4", "CDKN2A", "EGFR")
#'   net <- build_network(
#'     mirna_ids = mirna_ids,
#'     rna_symbols = rna_symbols,
#'     batch_size = 50,
#'     min_targets = 1
#'   )
#'   head(net)
#' }
#' }
build_network <- function(mirna_ids,
                          rna_symbols,
                          batch_size   = 50,
                          min_targets  = 5,
                          rescue_mirnas = c("hsa-miR-21-5p",  "hsa-miR-21-3p",
                                            "hsa-miR-155-5p", "hsa-miR-155-3p",
                                            "hsa-miR-196a-5p","hsa-miR-196b-5p",
                                            "hsa-miR-210-3p")) {

  message("Querying miRTarBase in batches of ", batch_size, "...")
  n_batches <- ceiling(length(mirna_ids) / batch_size)
  results   <- vector("list", n_batches)

  for (i in seq_len(n_batches)) {
    idx   <- ((i - 1) * batch_size + 1):min(i * batch_size, length(mirna_ids))
    batch <- mirna_ids[idx]
    message("  Batch ", i, " of ", n_batches)
    tryCatch({
      res <- multiMiR::get_multimir(org = "hsa", mirna = batch,
                                    table = "validated", summary = FALSE)
      if (!is.null(res) && nrow(res@data) > 0) results[[i]] <- res@data
    }, error = function(e) message("  Batch ", i, " failed: ", e$message))
    Sys.sleep(2)
  }

  # Rescue query for known important miRNAs
  if (length(rescue_mirnas) > 0) {
    message("Running rescue query for ", length(rescue_mirnas), " mature IDs...")
    tryCatch({
      res <- multiMiR::get_multimir(org = "hsa", mirna = rescue_mirnas,
                                    table = "validated", summary = FALSE)
      if (!is.null(res) && nrow(res@data) > 0)
        results[[n_batches + 1]] <- res@data
    }, error = function(e) message("Rescue query failed: ", e$message))
  }

  interactions <- dplyr::bind_rows(results)
  message("Total interactions retrieved: ", nrow(interactions))

  # Restrict to miRTarBase only
  interactions <- dplyr::filter(interactions, database == "mirtarbase")

  # Normalize miRNA IDs and filter to our targets
  interactions <- interactions |>
    dplyr::mutate(mirna_norm = normalize_mirna(mature_mirna_id)) |>
    dplyr::filter(target_symbol %in% rna_symbols)

  # Build per-miRNA target list
  mirna_targets <- interactions |>
    dplyr::select(mirna = mirna_norm, target = target_symbol) |>
    dplyr::distinct() |>
    dplyr::group_by(mirna) |>
    dplyr::summarise(targets   = list(target),
                     n_targets = dplyr::n(),
                     .groups   = "drop") |>
    dplyr::filter(n_targets >= min_targets)

  message("miRNAs with >= ", min_targets, " validated targets: ", nrow(mirna_targets))
  mirna_targets
}
