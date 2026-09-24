# Gate G-S3b: the chained BACE arm (mi_bace_chain) breaks the shared anchor that BACE's own final step
# has (Meng review B2). Negative control: in the shipped final runs, c1's design matrix is identical in
# every run, because c1 is fitted first and its predictors come from the one converged dataset. In the
# chained runs those predictors change from run to run.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "script", "campaign_gnn_off_lib.R"))
source(file.path(root, "script", "rubin_bace.R"))

cell <- make_cell("types_mixed", n = 40, seed = 1L, lambda = 0.7, rho = 0.5)
M <- 4L
set.seed(11L)
fb <- fit_bace_mi(cell, M = M, nitt = 3000L, burnin = 1000L, thin = 10L, runs = 2L)
ch <- mi_bace_chain(fb, cell$df_miss, M = M)
message(sprintf("[test-bace-chain] shipped fit %.1fs, chain %.1fs (M = %d)", fb$wall_s, ch$wall_s, M))

testthat::test_that("chained arm returns M datasets with observed cells untouched", {
  testthat::expect_length(ch$datasets, M)
  obs <- !is.na(cell$df_miss)
  for (d in ch$datasets) {
    for (v in names(cell$df_miss)) {
      o <- obs[, v]
      testthat::expect_identical(as.character(d[[v]][o]), as.character(cell$df_miss[[v]][o]))
    }
  }
})

testthat::test_that("shipped runs share c1's design matrix; chained runs do not", {
  X_ship <- lapply(fb$outb$final_results$all_models, function(m) as.matrix(m[["c1"]]$X))
  X_chain <- lapply(ch$models, function(m) as.matrix(m[["c1"]]$X))
  same <- function(L) all(vapply(L[-1], function(x) isTRUE(all.equal(x, L[[1]])), logical(1)))
  testthat::expect_true(same(X_ship))      # negative control: the shared anchor (B2)
  testthat::expect_false(same(X_chain))    # the chain varies c1's predictors across runs
})
