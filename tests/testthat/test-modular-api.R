test_that("compute_vss returns expected structure", {
  set.seed(11)
  patients <- paste0("P", 1:12)

  mirna_log <- matrix(rnorm(20 * length(patients)), nrow = 20,
                      dimnames = list(paste0("mir", 1:20), patients))
  rna_sym <- matrix(rnorm(20 * length(patients)), nrow = 20,
                    dimnames = list(paste0("g", 1:20), patients))
  mirna_targets <- data.frame(
    mirna = c("mir1", "mir2"),
    targets = I(list(paste0("g", 1:10), paste0("g", 11:20))),
    n_targets = c(10L, 10L),
    stringsAsFactors = FALSE
  )
  mirna_norm_map <- data.frame(
    original = rownames(mirna_log),
    norm = rownames(mirna_log),
    stringsAsFactors = FALSE
  )

  res <- compute_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)

  expect_true(all(c("mirna", "vss", "n_targets") %in% colnames(res)))
  expect_true(nrow(res) >= 1)
})

test_that("compute_cis returns expected structure", {
  set.seed(12)
  patients <- paste0("P", 1:12)

  mirna_log <- matrix(rnorm(2 * length(patients)), nrow = 2,
                      dimnames = list(c("mir1", "mir2"), patients))
  rna_sym <- matrix(rnorm(20 * length(patients)), nrow = 20,
                    dimnames = list(paste0("g", 1:20), patients))
  mirna_targets <- data.frame(
    mirna = c("mir1", "mir2"),
    targets = I(list(paste0("g", 1:10), paste0("g", 11:20))),
    n_targets = c(10L, 10L),
    stringsAsFactors = FALSE
  )
  mirna_norm_map <- data.frame(
    original = rownames(mirna_log),
    norm = rownames(mirna_log),
    stringsAsFactors = FALSE
  )

  res <- compute_cis(
    mirna_log = mirna_log,
    rna_sym = rna_sym,
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map,
    scored_mirnas = mirna_targets$mirna,
    n_perm = 5,
    max_targets = 20
  )

  expect_true(all(c("mirna", "coherence_score", "coherence_p") %in% colnames(res)))
})

test_that("classify_archetypes assigns known quadrants", {
  vss_scores <- data.frame(
    mirna = c("a", "b", "c", "d"),
    vss = c(0.9, 0.8, 0.2, 0.1),
    n_targets = c(10L, 10L, 10L, 10L),
    stringsAsFactors = FALSE
  )
  coherence_scores <- data.frame(
    mirna = c("a", "b", "c", "d"),
    coherence_score = c(0.9, 0.1, 0.8, 0.1),
    coherence_p = c(0.01, 0.20, 0.03, 0.60),
    stringsAsFactors = FALSE
  )

  res <- classify_archetypes(vss_scores, coherence_scores)

  expect_equal(nrow(res), 4)
  expect_equal(res$class[res$mirna == "a"], "dual")
  expect_equal(res$class[res$mirna == "b"], "silencer")
  expect_equal(res$class[res$mirna == "c"], "conductor")
  expect_equal(res$class[res$mirna == "d"], "weak")
})

test_that("run_mirBottleneck_project smoke test works on toy data", {
  toy_path <- function(f) system.file("extdata", f, package = "mirBottleneck", mustWork = FALSE)
  paths <- list(
    mirna_log = toy_path("toy_mirna_log.rds"),
    rna_sym = toy_path("toy_rna_sym.rds"),
    clinical_df = toy_path("toy_clinical_df.rds"),
    mirna_norm_map = toy_path("toy_mirna_norm_map.rds"),
    mirna_targets = toy_path("toy_mirna_targets.rds")
  )
  skip_if_not(all(nzchar(unlist(paths))), "Toy datasets not available")

  out <- file.path(tempdir(), "mirbottleneck_smoke")
  res <- run_mirBottleneck_project(
    mirna_log_rds = paths$mirna_log,
    rna_rds = paths$rna_sym,
    clinical_rds = paths$clinical_df,
    mirna_norm_map_rds = paths$mirna_norm_map,
    mirna_targets_rds = paths$mirna_targets,
    query_network = FALSE,
    out_dir = out,
    n_perm = 5,
    cox_p_threshold = 0.5,
    report = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("vss_scores", "coherence_scores", "combined_bottleneck") %in% names(res)))
})

test_that("SummarizedExperiment path has shape parity with matrix path", {
  skip_if_not_installed("SummarizedExperiment")

  set.seed(13)
  patients <- paste0("P", 1:12)
  mirna_log <- matrix(rnorm(20 * length(patients)), nrow = 20,
                      dimnames = list(paste0("mir", 1:20), patients))
  rna_sym <- matrix(rnorm(20 * length(patients)), nrow = 20,
                    dimnames = list(paste0("g", 1:20), patients))

  mirna_targets <- data.frame(
    mirna = c("mir1", "mir2"),
    targets = I(list(paste0("g", 1:10), paste0("g", 11:20))),
    n_targets = c(10L, 10L),
    stringsAsFactors = FALSE
  )
  mirna_norm_map <- data.frame(
    original = rownames(mirna_log),
    norm = rownames(mirna_log),
    stringsAsFactors = FALSE
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(mirna = mirna_log, mrna = rna_sym)
  )

  mat_vss <- compute_vss(mirna_log, rna_sym, mirna_targets, mirna_norm_map)
  se_vss <- compute_vss(
    mirna_targets = mirna_targets,
    mirna_norm_map = mirna_norm_map,
    se = se
  )

  expect_equal(ncol(mat_vss), ncol(se_vss))
  expect_equal(nrow(mat_vss), nrow(se_vss))
})

test_that("validation errors for missing assays and mismatched sample IDs are informative", {
  skip_if_not_installed("SummarizedExperiment")

  patients <- paste0("P", 1:12)
  mirna_log <- matrix(rnorm(2 * length(patients)), nrow = 2,
                      dimnames = list(c("mir1", "mir2"), patients))
  rna_sym <- matrix(rnorm(20 * length(patients)), nrow = 20,
                    dimnames = list(paste0("g", 1:20), patients))
  mirna_targets <- data.frame(
    mirna = c("mir1", "mir2"),
    targets = I(list(paste0("g", 1:10), paste0("g", 11:20))),
    n_targets = c(10L, 10L),
    stringsAsFactors = FALSE
  )
  mirna_norm_map <- data.frame(
    original = rownames(mirna_log),
    norm = rownames(mirna_log),
    stringsAsFactors = FALSE
  )

  se_missing <- SummarizedExperiment::SummarizedExperiment(
    assays = list(only_mirna = mirna_log)
  )

  expect_error(
    compute_vss(
      mirna_targets = mirna_targets,
      mirna_norm_map = mirna_norm_map,
      se = se_missing
    ),
    "Missing miRNA assay|Missing mRNA assay"
  )

  rna_mismatch <- rna_sym
  colnames(rna_mismatch) <- paste0("Q", 1:12)
  expect_error(
    compute_vss(
      mirna_log = mirna_log,
      rna_sym = rna_mismatch,
      mirna_targets = mirna_targets,
      mirna_norm_map = mirna_norm_map
    ),
    "Sample IDs in mirna_log and rna_sym do not match"
  )
})


test_that("plot_archetype_landscape returns a ggplot object", {
  classified <- data.frame(
    mirna = c("a", "b"),
    vss = c(0.8, 0.2),
    coherence_score = c(0.7, 0.1),
    class = c("dual", "weak"),
    stringsAsFactors = FALSE
  )

  p <- plot_archetype_landscape(classified)
  expect_s3_class(p, "ggplot")
})
