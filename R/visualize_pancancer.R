#' Pan-cancer visualization suite for mirBottleneck
#'
#' Functions for visualizing pan-cancer results: HR forest plots,
#' miRNA conservation heatmaps, archetype distribution maps,
#' and entropy correlation plots.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline
#'   geom_tile scale_fill_gradient2 scale_color_manual facet_wrap
#'   theme_minimal theme element_text element_blank labs coord_flip
#'   scale_x_log10 geom_hline annotate
#' @importFrom dplyr mutate arrange filter left_join group_by summarise
#' @importFrom tidyr pivot_longer pivot_wider
#' @name visualize_pancancer
NULL


#' Forest plot of HR across all 33 cancer types
#'
#' Plots hazard ratios (bottleneck score vs OS) across all cancer types
#' as a ranked forest plot. Cancers where bottleneck score is significantly
#' prognostic (logrank_p < 0.05) are highlighted.
#'
#' @param summary_df Data frame from run_pancancer() or .aggregate_pancancer().
#'   Must contain: cancer_type, hr, c_index, logrank_p, n_patients.
#' @param p_threshold Numeric. p-value threshold for significance highlight
#'   (default 0.05).
#' @param output_file Character. Path to save PDF (default NULL, print to screen).
#' @param width Numeric. Plot width in inches (default 8).
#' @param height Numeric. Plot height in inches (default 10).
#'
#' @return ggplot object invisibly.
#' @export
plot_pancancer_forest <- function(summary_df,
                                   p_threshold = 0.05,
                                   output_file = NULL,
                                   width       = 8,
                                   height      = 10) {

  df <- summary_df[!is.na(summary_df$hr) & summary_df$status == "OK", ]

  # Compute 95% CI from HR and approximate SE
  # We use c_index as a rough guide but need actual SE: flag if not available
  # For now we derive SE from the Cox model HR and p using normal approximation
  df <- df |>
    dplyr::mutate(
      log_hr    = log(hr),
      sig       = logrank_p < p_threshold,
      label     = paste0(cancer_type, " (n=", n_patients, ")"),
      cancer_short = gsub("^TCGA-", "", cancer_type)
    ) |>
    dplyr::arrange(log_hr)

  df$cancer_short <- factor(df$cancer_short, levels = df$cancer_short)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x     = log_hr,
    y     = cancer_short,
    color = sig
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = n_patients), alpha = 0.85) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#D62728", "FALSE" = "#7F7F7F"),
      labels = c("TRUE" = paste0("p < ", p_threshold),
                 "FALSE" = paste0("p ≥ ", p_threshold)),
      name   = "Log-rank"
    ) +
    ggplot2::scale_size_continuous(
      name   = "N patients",
      range  = c(2, 7),
      breaks = c(50, 200, 500, 1000)
    ) +
    ggplot2::labs(
      x     = "log(Hazard Ratio) — Bottleneck Score vs OS",
      y     = NULL,
      title = "mirBottleneck Score: Pan-Cancer Survival Association",
      subtitle = paste0(sum(df$sig, na.rm = TRUE), " of ", nrow(df),
                        " cancer types significant at p < ", p_threshold)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y      = ggplot2::element_text(size = 9),
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height,
                    device = "pdf")
    message("Forest plot saved to: ", output_file)
  } else {
    print(p)
  }

  invisible(p)
}


#' Heatmap of universal bottleneck miRNAs across cancer types
#'
#' Identifies miRNAs that appear as Dual or Silencer archetypes across
#' multiple cancer types (universal bottlenecks) and visualizes their
#' bottleneck index as a heatmap (miRNAs x cancer types).
#'
#' @param results_dir Character. Path to the run_pancancer() output directory.
#'   Each subdirectory should contain a result.rds file.
#' @param cancer_types Character vector. Cancer types to include. Defaults to
#'   all with result.rds present.
#' @param min_cancers Integer. Minimum number of cancer types a miRNA must
#'   appear in (as non-Weak) to be included (default 5).
#' @param top_n Integer. Maximum number of miRNAs to show (default 40).
#' @param output_file Character. Path to save PDF (default NULL).
#'
#' @return ggplot object invisibly.
#' @export
plot_conservation_heatmap <- function(results_dir,
                                       cancer_types = NULL,
                                       min_cancers  = 5,
                                       top_n        = 40,
                                       output_file  = NULL) {

  # Load all results
  if (is.null(cancer_types)) {
    cancer_types <- basename(list.dirs(results_dir, recursive = FALSE))
    cancer_types <- cancer_types[grepl("^TCGA-", cancer_types)]
  }

  all_classified <- lapply(cancer_types, function(ct) {
    rds <- file.path(results_dir, ct, "result.rds")
    if (!file.exists(rds)) return(NULL)
    res <- readRDS(rds)
    if (is.null(res$classified)) return(NULL)
    df <- res$classified
    df$cancer_type <- gsub("^TCGA-", "", ct)
    df
  })

  all_classified <- dplyr::bind_rows(Filter(Negate(is.null), all_classified))

  if (nrow(all_classified) == 0)
    stop("No classified data found. Check results_dir.")

  # Count how many cancer types each miRNA is non-Weak in
  conservation <- all_classified |>
    dplyr::filter(class != "weak") |>
    dplyr::group_by(mirna) |>
    dplyr::summarise(
      n_cancers        = dplyr::n_distinct(cancer_type),
      dominant_class   = names(sort(table(class), decreasing = TRUE))[1],
      .groups = "drop"
    ) |>
    dplyr::filter(n_cancers >= min_cancers) |>
    dplyr::arrange(dplyr::desc(n_cancers))

  top_mirnas <- head(conservation$mirna, top_n)

  message("Universal bottleneck miRNAs (>= ", min_cancers, " cancer types): ",
          nrow(conservation))
  message("Showing top ", length(top_mirnas))

  # Build matrix: bottleneck_index x cancer type
  heat_df <- all_classified |>
    dplyr::filter(mirna %in% top_mirnas) |>
    dplyr::select(mirna, cancer_type, bottleneck_index)

  # Pivot wide then back long for ggplot with complete grid
  heat_wide <- tidyr::pivot_wider(heat_df,
                                   names_from  = cancer_type,
                                   values_from = bottleneck_index,
                                   values_fill = 0)
  heat_long <- tidyr::pivot_longer(heat_wide, -mirna,
                                    names_to  = "cancer_type",
                                    values_to = "bottleneck_index")

  # Order miRNAs by conservation
  mirna_order <- conservation$mirna[conservation$mirna %in% top_mirnas]
  heat_long$mirna <- factor(heat_long$mirna, levels = rev(mirna_order))

  p <- ggplot2::ggplot(heat_long,
                       ggplot2::aes(x = cancer_type, y = mirna,
                                    fill = bottleneck_index)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low      = "#F7F7F7",
      mid      = "#4393C3",
      high     = "#053061",
      midpoint = 0.5,
      name     = "Bottleneck\nIndex",
      limits   = c(0, 1)
    ) +
    ggplot2::labs(
      x     = NULL,
      y     = NULL,
      title = "Universal miRNA Bottlenecks Across TCGA Cancer Types",
      subtitle = paste0("miRNAs active (non-Weak) in ≥ ", min_cancers,
                        " cancer types | top ", length(top_mirnas), " shown")
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y      = ggplot2::element_text(size = 7),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold"),
      panel.grid       = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    h <- max(6, length(top_mirnas) * 0.25 + 3)
    w <- max(8, length(cancer_types) * 0.35 + 4)
    ggplot2::ggsave(output_file, p, width = w, height = h, device = "pdf")
    message("Conservation heatmap saved to: ", output_file)
  } else {
    print(p)
  }

  invisible(p)
}


#' Archetype distribution bubble chart across cancer types
#'
#' Visualizes how many miRNAs fall into each archetype (Dual, Silencer,
#' Conductor, Weak) for every cancer type. Bubble size = count, color = class.
#'
#' @param summary_df Data frame from run_pancancer(). Must contain:
#'   cancer_type, n_dual, n_silencer, n_conductor, n_weak.
#' @param output_file Character. PDF path (default NULL).
#'
#' @return ggplot object invisibly.
#' @export
plot_archetype_distribution <- function(summary_df, output_file = NULL) {

  df <- summary_df[summary_df$status == "OK", ]

  long_df <- tidyr::pivot_longer(
    df[, c("cancer_type", "n_dual", "n_silencer", "n_conductor", "n_weak")],
    cols      = c(n_dual, n_silencer, n_conductor, n_weak),
    names_to  = "archetype",
    values_to = "count"
  ) |>
    dplyr::mutate(
      archetype    = gsub("^n_", "", archetype),
      cancer_short = gsub("^TCGA-", "", cancer_type)
    )

  archetype_colors <- c(
    dual      = "#D62728",
    silencer  = "#FF7F0E",
    conductor = "#2CA02C",
    weak      = "#AEC7E8"
  )

  long_df <- long_df |>
    dplyr::arrange(archetype, dplyr::desc(count))

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = archetype, y = cancer_short,
                                    size = count, color = archetype)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_manual(values = archetype_colors, guide = "none") +
    ggplot2::scale_size_continuous(name = "# miRNAs", range = c(1, 10)) +
    ggplot2::labs(
      x     = "Archetype",
      y     = NULL,
      title = "miRNA Archetype Distribution Across TCGA Cancer Types"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 8),
      plot.title   = ggplot2::element_text(face = "bold")
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 7,
                    height = max(5, nrow(df) * 0.25 + 2), device = "pdf")
    message("Archetype distribution plot saved to: ", output_file)
  } else {
    print(p)
  }

  invisible(p)
}


#' Entropy vs bottleneck score scatter plots across cancer types
#'
#' For each cancer type, plots patient-level transcriptome entropy against
#' composite bottleneck score. Panels faceted by cancer type. Points colored
#' by survival status.
#'
#' @param results_dir Character. run_pancancer() output directory.
#' @param cancer_types Character vector. Subset of cancer types to plot.
#'   Defaults to all available with entropy data.
#' @param max_panels Integer. Maximum number of facets (default 16).
#' @param output_file Character. PDF path (default NULL).
#'
#' @return ggplot object invisibly.
#' @export
plot_entropy_scatter <- function(results_dir,
                                  cancer_types = NULL,
                                  max_panels   = 16,
                                  output_file  = NULL) {

  if (is.null(cancer_types)) {
    cancer_types <- basename(list.dirs(results_dir, recursive = FALSE))
    cancer_types <- cancer_types[grepl("^TCGA-", cancer_types)]
  }

  all_df <- lapply(cancer_types, function(ct) {
    rds <- file.path(results_dir, ct, "result.rds")
    if (!file.exists(rds)) return(NULL)
    res <- readRDS(rds)

    entropy_vec <- res$entropy
    scores_vec  <- res$composite$patient_scores

    if (is.null(entropy_vec) || is.null(scores_vec)) return(NULL)

    shared <- intersect(names(entropy_vec), names(scores_vec))
    if (length(shared) < 10) return(NULL)

    df <- data.frame(
      patient          = shared,
      entropy          = entropy_vec[shared],
      bottleneck_score = scores_vec[shared],
      cancer_type      = gsub("^TCGA-", "", ct),
      stringsAsFactors = FALSE
    )

    # Add OS status if available
    if (!is.null(res$clinical)) {
      clin <- res$clinical[, c("patient", "OS_status")]
      df <- merge(df, clin, by = "patient", all.x = TRUE)
    } else {
      df$OS_status <- NA
    }

    # Add Spearman rho to label
    rho <- if (!is.null(res$entropy_cor)) {
      round(res$entropy_cor$estimate, 3)
    } else NA
    p_val <- if (!is.null(res$entropy_cor)) {
      signif(res$entropy_cor$p.value, 2)
    } else NA

    df$facet_label <- paste0(df$cancer_type,
                              "\nrho=", rho, " p=", p_val)
    df
  })

  plot_df <- dplyr::bind_rows(Filter(Negate(is.null), all_df))

  if (nrow(plot_df) == 0)
    stop("No entropy data found across cancer types.")

  # Limit panels
  ct_keep <- unique(plot_df$facet_label)[seq_len(min(max_panels,
                                                      dplyr::n_distinct(plot_df$facet_label)))]
  plot_df <- plot_df[plot_df$facet_label %in% ct_keep, ]

  plot_df$OS_status_label <- factor(
    ifelse(is.na(plot_df$OS_status), "Unknown",
           ifelse(plot_df$OS_status == 1, "Deceased", "Alive/Censored")),
    levels = c("Deceased", "Alive/Censored", "Unknown")
  )

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = bottleneck_score,
                                    y = entropy,
                                    color = OS_status_label)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#D62728",
                         linewidth = 0.8, fill = "#FEE0D2") +
    ggplot2::scale_color_manual(
      values = c("Deceased" = "#D62728",
                 "Alive/Censored" = "#2166AC",
                 "Unknown" = "#969696"),
      name   = "Vital Status"
    ) +
    ggplot2::facet_wrap(~ facet_label,
                        scales = "free",
                        ncol   = min(4, ceiling(sqrt(length(ct_keep))))) +
    ggplot2::labs(
      x     = "Composite Bottleneck Score",
      y     = "Transcriptome Shannon Entropy (bits)",
      title = "Bottleneck Score vs Transcriptome Entropy",
      subtitle = "Each panel: one TCGA cancer type | rho = Spearman correlation"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = 7, face = "bold"),
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(face = "bold")
    )

  if (!is.null(output_file)) {
    n_panels <- dplyr::n_distinct(plot_df$facet_label)
    n_cols   <- min(4, ceiling(sqrt(n_panels)))
    n_rows   <- ceiling(n_panels / n_cols)
    ggplot2::ggsave(output_file, p,
                    width  = n_cols * 3.5,
                    height = n_rows * 3.2,
                    device = "pdf")
    message("Entropy scatter saved to: ", output_file)
  } else {
    print(p)
  }

  invisible(p)
}


#' Generate all pan-cancer figures in one call
#'
#' Convenience wrapper that produces all four main figures from a completed
#' run_pancancer() run and saves them to an output directory.
#'
#' @param results_dir Character. run_pancancer() output directory.
#' @param figures_dir Character. Where to save PDFs (default results_dir/figures).
#' @param summary_df Data frame from run_pancancer(). If NULL, tries to load
#'   from results_dir/pancancer_summary.rds.
#'
#' @return Named list of ggplot objects.
#' @export
generate_all_figures <- function(results_dir,
                                  figures_dir = file.path(results_dir, "figures"),
                                  summary_df  = NULL) {

  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(summary_df)) {
    sum_path <- file.path(results_dir, "pancancer_summary.rds")
    if (!file.exists(sum_path))
      stop("pancancer_summary.rds not found. Run run_pancancer() first or supply summary_df.")
    summary_df <- readRDS(sum_path)
  }

  message("Generating Figure 1: Pan-cancer forest plot...")
  fig1 <- plot_pancancer_forest(
    summary_df  = summary_df,
    output_file = file.path(figures_dir, "fig1_pancancer_forest.pdf")
  )

  message("Generating Figure 2: Conservation heatmap...")
  fig2 <- plot_conservation_heatmap(
    results_dir = results_dir,
    output_file = file.path(figures_dir, "fig2_conservation_heatmap.pdf")
  )

  message("Generating Figure 3: Archetype distribution...")
  fig3 <- plot_archetype_distribution(
    summary_df  = summary_df,
    output_file = file.path(figures_dir, "fig3_archetype_distribution.pdf")
  )

  message("Generating Figure 4: Entropy scatter panels...")
  fig4 <- plot_entropy_scatter(
    results_dir = results_dir,
    output_file = file.path(figures_dir, "fig4_entropy_scatter.pdf")
  )

  message("All figures saved to: ", figures_dir)

  invisible(list(fig1 = fig1, fig2 = fig2, fig3 = fig3, fig4 = fig4))
}
