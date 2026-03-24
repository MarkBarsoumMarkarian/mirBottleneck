#' Generate a self-contained HTML analysis report
#'
#' Renders a self-contained interactive HTML report from the output of
#' \code{run_mirBottleneck_project()}. The report includes cohort summary,
#' network statistics, VSS and coherence distributions, the bottleneck
#' phase diagram, a sortable miRNA table, Kaplan-Meier survival curves,
#' Cox model results, per-patient scores, and session info.
#'
#' The report is fully offline — no server, no Shiny, no external calls.
#' It opens in any browser and can be shared as a single ~3MB HTML file.
#'
#' @param results List returned by \code{run_mirBottleneck_project()}.
#' @param out_dir Character. Directory to write the HTML file into.
#'   Defaults to the current working directory.
#' @param cohort_name Character. Cohort label displayed in the report header
#'   (default "TCGA-PAAD").
#' @param filename Character. Output filename (default
#'   "mirBottleneck_report.html").
#' @param open Logical. Open the report in the default browser after
#'   rendering (default FALSE).
#'
#' @return Invisibly returns the path to the generated HTML file.
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("DT", quietly = TRUE) &&
#'     requireNamespace("rmarkdown", quietly = TRUE)) {
#'   paths   <- mirBottleneck:::.toy_paths()
#'   results <- run_mirBottleneck_project(
#'     mirna_log_rds      = paths$mirna_log,
#'     rna_rds            = paths$rna_sym,
#'     clinical_rds       = paths$clinical_df,
#'     mirna_norm_map_rds = paths$mirna_norm_map,
#'     mirna_targets_rds  = paths$mirna_targets,
#'     query_network      = FALSE,
#'     out_dir            = tempdir(),
#'     cox_p_threshold    = 0.5
#'   )
#'   generate_report(results, out_dir = tempdir(), open = FALSE)
#' }
#'
#' @export
generate_report <- function(results,
                            out_dir     = getwd(),
                            cohort_name = "TCGA-PAAD",
                            filename    = "mirBottleneck_report.html",
                            open        = FALSE) {

  # ---- check dependencies --------------------------------------------------
  for (pkg in c("rmarkdown", "plotly", "DT", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required to generate the report. ",
           "Install with: install.packages('", pkg, "')")
    }
  }

  # ---- locate template -----------------------------------------------------
  tmpl <- system.file("report", "mirBottleneck_report.Rmd",
                      package = "mirBottleneck", mustWork = TRUE)

  # ---- validate results list -----------------------------------------------
  required <- c("combined_bottleneck", "cox_individual", "cox_final",
                "survival_df", "mirna_targets", "clinical_df",
                "vss_scores", "coherence_scores",
                "patient_scores", "contributing_mirnas")

  missing_keys <- setdiff(required, names(results))
  if (length(missing_keys) > 0) {
    stop("results is missing required elements: ",
         paste(missing_keys, collapse = ", "),
         "\nMake sure you are passing the full list from run_mirBottleneck_project().")
  }

  # ---- output path ---------------------------------------------------------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_path <- file.path(out_dir, filename)

  # ---- render --------------------------------------------------------------
  message("Rendering mirBottleneck HTML report...")
  message("  Template : ", tmpl)
  message("  Output   : ", out_path)

  rmarkdown::render(
    input       = tmpl,
    output_file = out_path,
    params      = list(
      results     = results,
      cohort_name = cohort_name,
      pkg_version = as.character(utils::packageVersion("mirBottleneck"))
    ),
    envir       = new.env(parent = globalenv()),
    quiet       = TRUE
  )

  message("Done! Report saved to: ", out_path)

  if (open) utils::browseURL(out_path)

  invisible(out_path)
}
