# Gate G-S3c: BACE input cleaning. At lambda = 1 with fixed thresholds (n = 300, seed 105, L'Ecuyer RNG as in
# rubin_cell.R) the ordinal trait has an empty level and BACE as shipped stops at once ("Mixed model equations
# singular"). fit_bace_mi() drops the empty level, records it, and the fit completes; the shipped datasets keep
# df_miss's columns and levels.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({ source(file.path(root, "script", "campaign_gnn_off_lib.R")); source(file.path(root, "script", "rubin_bace.R")) })
old_rng <- RNGkind()[1]
RNGkind("L'Ecuyer-CMRG")
cell <- make_cell("types_mixed", 300L, 105L, lambda = 1, rho = 0, thresholds = "fixed", driver = TRUE)

testthat::test_that("the case really has an empty level (positive control for the fix)", {
  testthat::expect_true(any(table(cell$df_miss$ord) == 0))
})

testthat::test_that("fit_bace_mi drops the empty level, records it, and the shipped datasets stay aligned", {
  set.seed(1)
  fb <- fit_bace_mi(cell, M = 3L, nitt = 3000L, burnin = 1000L, thin = 10L, runs = 2L)
  testthat::expect_true(any(grepl("^ord: dropped", fb$diag$input_fix)))
  s <- mi_bace_shipped(fb$outb, cell$df_miss)
  testthat::expect_length(s, 3L)
  testthat::expect_identical(names(s[[1]]), names(cell$df_miss))
  testthat::expect_identical(levels(s[[1]]$ord), levels(cell$df_miss$ord))
  testthat::expect_false(anyNA(s[[1]]$c1))
})

RNGkind(old_rng)  # do not leak the RNG kind into later test files (make_cell() datasets depend on it)
