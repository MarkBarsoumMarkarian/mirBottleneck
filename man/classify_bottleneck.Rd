# classify.R - Bottleneck classification and combined index

#' Classify miRNAs and compute combined bottleneck index
#'
#' @param vss_scores dataframe from score_vss()
#' @param coherence_scores dataframe from score_coherence()
#' @return dataframe with mirna, vss, coherence_score, bottleneck_index, class
#' @examples
#' vss_scores <- data.frame(
#'   mirna = c("hsa-mir-1", "hsa-mir-2"),
#'   vss = c(0.2, 0.8),
#'   n_targets = c(10, 12)
#' )
#' coherence_scores <- data.frame(
#'   mirna = c("hsa-mir-1", "hsa-mir-2"),
#'   coherence_score = c(0.1, 0.6),
#'   coherence_p = c(0.5, 0.02)
#' )
#' classify_bottleneck(vss_scores, coherence_scores)
#' @export
classify_bottleneck <- function(vss_scores, coherence_scores) {
  combined <- dplyr::inner_join(vss_scores, coherence_scores, by = "mirna")
  combined <- combined |>
    dplyr::mutate(
      vss_norm         = normalize_01(vss),
      coherence_norm   = normalize_01(coherence_score),
      bottleneck_index = (vss_norm + coherence_norm) / 2
    )
  vss_med <- median(combined$vss_norm,       na.rm = TRUE)
  coh_med <- median(combined$coherence_norm, na.rm = TRUE)
  combined <- combined |>
    dplyr::mutate(
      class = dplyr::case_when(
        vss_norm >= vss_med & coherence_norm >= coh_med ~ "dual",
        vss_norm >= vss_med & coherence_norm <  coh_med ~ "silencer",
        vss_norm <  vss_med & coherence_norm >= coh_med ~ "conductor",
        TRUE ~ "weak"
      )
    ) |>
    dplyr::arrange(dplyr::desc(bottleneck_index))
  message("Class distribution:")
  print(table(combined$class))
  combined
}
