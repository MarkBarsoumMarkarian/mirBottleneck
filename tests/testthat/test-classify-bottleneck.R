test_that("classify_bottleneck returns required columns and class values", {
  vss_scores <- data.frame(
    mirna     = c("a", "b", "c", "d"),
    vss       = c(0.9, 0.1, 0.8, 0.2),
    n_targets = c(10,  10,  10,  10),
    stringsAsFactors = FALSE
  )

  coherence_scores <- data.frame(
    mirna           = c("a", "b", "c", "d"),
    coherence_score = c(0.8, 0.7, 0.1, 0.2),
    coherence_p     = c(0.01, 0.02, 0.4, 0.5),
    stringsAsFactors = FALSE
  )

  result <- classify_bottleneck(vss_scores, coherence_scores)

  # Required columns present
  expect_true(all(c("mirna", "vss", "coherence_score",
                    "bottleneck_index", "class") %in% colnames(result)))

  # All 4 miRNAs returned
  expect_equal(nrow(result), 4)

  # Valid class values only (lowercase)
  valid_classes <- c("dual", "silencer", "conductor", "weak")
  expect_true(all(result$class %in% valid_classes))

  # bottleneck_index in [0, 1]
  expect_true(all(result$bottleneck_index >= 0))
  expect_true(all(result$bottleneck_index <= 1))

  # "a": high VSS + high coherence → dual
  expect_equal(result$class[result$mirna == "a"], "dual")

  # "c": high VSS + low coherence → silencer
  expect_equal(result$class[result$mirna == "c"], "silencer")
})

test_that("classify_bottleneck handles single-miRNA input", {
  vss       <- data.frame(mirna = "x", vss = 0.5, n_targets = 5)
  coherence <- data.frame(mirna = "x", coherence_score = 0.5, coherence_p = 0.1)
  result    <- classify_bottleneck(vss, coherence)
  expect_equal(nrow(result), 1)
  expect_true(result$class %in% c("dual", "silencer", "conductor", "weak"))
})
