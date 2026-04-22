# classify.R - Bottleneck classification and combined index

#' Classify miRNAs and compute combined bottleneck index
#'
#' Backward-compatible wrapper around `classify_archetypes()`.
#'
#' @param vss_scores dataframe from score_vss()
#' @param coherence_scores dataframe from score_coherence()
#' @return dataframe with mirna, vss, coherence_score, bottleneck_index, class
#' @export
#' @examples
#' vss_scores <- data.frame(
#'   mirna = c("hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-miR-1-3p"),
#'   vss = c(0.9, 0.4, 0.1)
#' )
#' coherence_scores <- data.frame(
#'   mirna = c("hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-miR-1-3p"),
#'   coherence_score = c(0.8, 0.3, 0.2)
#' )
#' classify_bottleneck(vss_scores, coherence_scores)
classify_bottleneck <- function(vss_scores, coherence_scores) {
  classify_archetypes(vss_scores = vss_scores, coherence_scores = coherence_scores)
}
