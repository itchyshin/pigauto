# script/tests-rubin/test-lib-estimands.R
# Gate G-S1c for est_pgls_slope() / est_phylo_cor(). 100 reps of
# make_cell("types_mixed", n = 100, seed = i, lambda = 0.7, rho = 0.5) on the COMPLETE truth (no
# missingness): the population PGLS slope of c2 ~ c1 and the population phylogenetic correlation
# are both rho = 0.5 (see script/campaign_gnn_off_lib.R header comment on sim_latents()).
#
# campaign_gnn_off_lib.R (read-only; sourced here only to call make_cell()) needs ape + rcoal etc,
# already a pigauto dependency.

source(normalizePath(file.path(testthat::test_path(), "..", "rubin_lib.R")))
source(normalizePath(file.path(testthat::test_path(), "..", "campaign_gnn_off_lib.R")))

testthat::test_that("est_pgls_slope / est_phylo_cor recover rho = 0.5 on average over 100 reps", {
  n_rep <- 100L
  rho_true <- 0.5
  t0 <- Sys.time()

  slopes <- numeric(n_rep)
  ses_slope <- numeric(n_rep)
  rs <- numeric(n_rep)
  ses_r <- numeric(n_rep)
  converged <- logical(n_rep)

  for (i in seq_len(n_rep)) {
    cell <- make_cell("types_mixed", n = 100, seed = i, lambda = 0.7, rho = rho_true)
    truth <- cell$truth; tree <- cell$tree

    slope_fit <- est_pgls_slope(truth, tree)
    converged[i] <- slope_fit$converged
    slopes[i] <- slope_fit$estimate
    ses_slope[i] <- sqrt(slope_fit$variance)

    cor_fit <- est_phylo_cor(truth, tree, lambda = slope_fit$lambda_hat)
    rs[i] <- cor_fit$r
    ses_r[i] <- sqrt(cor_fit$variance)
  }

  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  message(sprintf(
    "G-S1c: mean slope = %.4f, mean r = %.4f, %d/%d converged, elapsed = %.1fs",
    mean(slopes), mean(rs), sum(converged), n_rep, elapsed))

  testthat::expect_true(all(converged))
  testthat::expect_true(all(is.finite(ses_slope)))
  testthat::expect_true(all(is.finite(ses_r)))
  testthat::expect_equal(mean(slopes), rho_true, tolerance = 0.05, scale = 1)
  testthat::expect_equal(mean(rs), rho_true, tolerance = 0.05, scale = 1)
  testthat::expect_lt(elapsed, 180)
})

testthat::test_that("est_phylo_cor returns NA fields, not an error, when lambda is not finite (Meng B1)", {
  cell <- make_cell("types_mixed", n = 40, seed = 3L, lambda = 0.7, rho = 0.5)
  r <- est_phylo_cor(cell$truth, cell$tree, lambda = NA_real_)
  testthat::expect_true(is.na(r$r) && is.na(r$z))
  testthat::expect_equal(r$variance, 1 / 37)
})
