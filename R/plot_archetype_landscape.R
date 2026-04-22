#' Plot miRNA archetype landscape
#'
#' Creates a VSS-vs-coherence scatter plot colored by archetype class and returns
#' a standard `ggplot2` object for downstream customization.
#'
#' @param classified_df Data frame from `classify_archetypes()` or `classify_bottleneck()`.
#'   Must contain columns `vss`, `coherence_score`, and `class`.
#' @return A `ggplot2` object.
#' @export
plot_archetype_landscape <- function(classified_df) {
  req_cols <- c("vss", "coherence_score", "class")
  if (!all(req_cols %in% colnames(classified_df))) {
    stop(
      "classified_df must contain columns: ",
      paste(req_cols, collapse = ", "),
      call. = FALSE
    )
  }

  med_vss <- median(classified_df$vss, na.rm = TRUE)
  med_coh <- median(classified_df$coherence_score, na.rm = TRUE)

  ggplot2::ggplot(
    classified_df,
    ggplot2::aes(x = vss, y = coherence_score, color = class)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::geom_vline(xintercept = med_vss, linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = med_coh, linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      x = "Variance Suppression Score (VSS)",
      y = "Coherence Induction Score",
      color = "Archetype"
    ) +
    ggplot2::theme_minimal()
}
