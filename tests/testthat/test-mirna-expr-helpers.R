test_that(".build_mirna_expr returns named list of patient vectors", {
  mirna_log <- matrix(
    c(1, 2, 3,
      4, 5, 6,
      10, 20, 30),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("hsa-mir-21-1", "hsa-mir-21-2", "hsa-mir-155"),
      c("P1", "P2", "P3")
    )
  )

  mirna_norm_map <- data.frame(
    original = c("hsa-mir-21-1", "hsa-mir-21-2", "hsa-mir-155"),
    norm     = c("hsa-mir-21",   "hsa-mir-21",   "hsa-mir-155"),
    stringsAsFactors = FALSE
  )

  norm_ids <- mirna_norm_map$norm[match(rownames(mirna_log),
                                        mirna_norm_map$original)]

  result <- mirBottleneck:::.build_mirna_expr(
    mirna_log     = mirna_log,
    norm_ids      = norm_ids,
    mirna_norm_map = mirna_norm_map
  )

  expect_type(result, "list")
  expect_named(result, c("hsa-mir-21", "hsa-mir-155"), ignore.order = TRUE)

  # hsa-mir-21: mean of rows 1 & 2 per patient → c(2.5, 3.5, 4.5)
  expect_equal(result[["hsa-mir-21"]], c(P1 = 2.5, P2 = 3.5, P3 = 4.5))

  # hsa-mir-155: single row, no averaging
  expect_equal(result[["hsa-mir-155"]], c(P1 = 10, P2 = 20, P3 = 30))
})

test_that("harmonize_barcode truncates to 12 characters", {
  expect_equal(harmonize_barcode("TCGA-3A-A9I7-01A-21R-A38N-13"), "TCGA-3A-A9I7")
  expect_equal(harmonize_barcode("TCGA-AB-1234"), "TCGA-AB-1234")
  expect_equal(nchar(harmonize_barcode("TCGA-XX-XXXX-01A")), 12)
})

test_that("normalize_01 maps to [0, 1]", {
  x      <- c(2, 4, 6, 8, 10)
  result <- mirBottleneck:::normalize_01(x)
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
  expect_equal(length(result), length(x))
})

test_that("normalize_01 returns all zeros for zero-variance input", {
  x      <- c(5, 5, 5)
  result <- mirBottleneck:::normalize_01(x)
  expect_true(all(result == 0))
})
