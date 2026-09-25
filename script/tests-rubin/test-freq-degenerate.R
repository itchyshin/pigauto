# Gate G-S2d: numerical failures in the frequentist draws are caught, not returned as garbage.
# (1) A parameter set whose joint covariance is not positive definite makes cond_draw() stop with a
#     "degenerate conditional distribution" error (before 2026-09-24 a silent ginv() fallback returned
#     draws of 1e6 to 1e67 on about 1 in 300 bootstrap refits).
# (2) mi_freq_A() redraws the bootstrap sample when a conditional draw fails, counts it in n_degenerate,
#     and still returns M datasets.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({ source(file.path(root, "script", "campaign_gnn_off_lib.R")); source(file.path(root, "script", "rubin_freq.R")) })
RNGkind("L'Ecuyer-CMRG")
cell <- make_cell("types_mixed", 60L, 3L, lambda = 0.7, rho = 0.5, thresholds = "fixed", driver = TRUE)
bt <- default_block_traits(cell)
pars <- fit_block(cell$df_miss, cell$tree, bt, cell$trait_types)
Yt <- transform_block(cell$df_miss, bt, pars$is_prp)

testthat::test_that("a non-positive-definite parameter set raises the degenerate error", {
  bad <- pars; bad$Sigma_p[1, 2] <- bad$Sigma_p[2, 1] <- 3 * sqrt(bad$Sigma_p[1, 1] * bad$Sigma_p[2, 2])
  testthat::expect_error(cond_draw(Yt, bad, cell$tree, 2L), "degenerate conditional distribution")
  testthat::expect_silent(cond_draw(Yt, pars, cell$tree, 2L))        # the fitted parameters still pass
})

testthat::test_that("mi_freq_A redraws after a failed conditional draw and counts it", {
  real <- cond_draw; calls <- 0L
  assign("cond_draw", function(...) { calls <<- calls + 1L; if (calls %in% c(2L, 5L)) stop("degenerate conditional distribution: test") else real(...) },
         envir = globalenv())
  on.exit(assign("cond_draw", real, envir = globalenv()))
  set.seed(4)
  res <- mi_freq_A(cell, 6L)
  testthat::expect_identical(res$n_degenerate, 2L)
  testthat::expect_length(Filter(Negate(is.null), res$datasets), 6L)
})
