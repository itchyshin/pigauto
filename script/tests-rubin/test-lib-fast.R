# The fast eigenbasis PGLS slope must reproduce est_pgls_slope() (nlme::gls + corPagel, REML, bounded
# fallback): estimate, variance and lambda, across interior and boundary lambda.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({ source(file.path(root, "script", "campaign_gnn_off_lib.R")); source(file.path(root, "script", "rubin_lib.R")) })

grid <- expand.grid(lambda = c(0.3, 0.7, 1), n = c(60L, 150L), seed = 1:6)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  cell <- make_cell("types_mixed", g$n, 500L + g$seed, lambda = g$lambda, rho = 0.5, thresholds = "fixed", driver = TRUE)
  a <- est_pgls_slope(cell$truth, cell$tree); b <- est_pgls_slope_fast(cell$truth, cell$tree)
  data.frame(g, est_a = a$estimate, est_b = b$estimate, var_a = a$variance, var_b = b$variance,
             lam_a = a$lambda_hat, lam_b = b$lambda_hat)
}))

testthat::test_that("fast PGLS matches nlme::gls REML with corPagel (36 datasets, lambda 0.3 to 1)", {
  message(sprintf("[test-lib-fast] max |d est| %.2e, max rel d var %.2e, max |d lambda| %.2e",
                  max(abs(res$est_a - res$est_b)), max(abs(res$var_b / res$var_a - 1)), max(abs(res$lam_a - res$lam_b))))
  testthat::expect_lt(max(abs(res$est_a - res$est_b)), 1e-4)
  testthat::expect_lt(max(abs(res$var_b / res$var_a - 1)), 1e-4)
  testthat::expect_lt(max(abs(res$lam_a - res$lam_b)), 1e-4)
})

testthat::test_that("a supplied eigendecomposition gives the same answer", {
  cell <- make_cell("types_mixed", 80L, 9L, lambda = 0.7, rho = 0.5)
  e <- pagel_eigen(cell$tree, rownames(cell$truth))
  testthat::expect_equal(est_pgls_slope_fast(cell$truth, cell$tree, eig = e)$estimate,
                         est_pgls_slope_fast(cell$truth, cell$tree)$estimate, tolerance = 1e-12)
})
