# script/tests-rubin/test-freq-draws.R
# Gate G-S2b for cond_draw()'s sampling step. n = 30, p = 2, known (not fitted) parameters; 5,000
# cond_draw() draws. Checks: (1) observed cells are returned exactly unchanged in every draw;
# (2) the empirical covariance of the drawn missing cells matches the analytic conditional
# covariance joint_cov()/cond_draw() computes (V_mm - V_mo solve(V_oo) V_om), max abs diff < 0.05.
# This is a pure sampling-fidelity check -- G-S2a already establishes that the conditional-mean/
# covariance machinery matches Rphylopars' own answer on real fitted data.

source(normalizePath(file.path(testthat::test_path(), "..", "rubin_freq.R")))

testthat::test_that("cond_draw reproduces the analytic conditional mean/covariance under known params", {
  set.seed(20260924)
  n <- 30L; p <- 2L
  tree <- ape::rcoal(n)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  sp <- tree$tip.label

  mu <- c(t1 = 0, t2 = 1)
  Sigma_p <- matrix(c(1, 0.4, 0.4, 0.8), 2, 2, dimnames = list(c("t1", "t2"), c("t1", "t2")))
  lambda <- 0.6
  pars <- list(mu = mu, Sigma_p = Sigma_p, Sigma_e = NULL, lambda = lambda,
               block_traits = c("t1", "t2"), is_prp = c(FALSE, FALSE), species = sp)

  jc <- joint_cov(pars, tree)
  Lfull <- chol(jc$V)
  y_full <- jc$mu_vec + as.vector(crossprod(Lfull, stats::rnorm(n * p)))
  Yfull <- matrix(y_full, n, p, dimnames = list(sp, c("t1", "t2")))

  Y <- Yfull
  Y[1:10, "t2"] <- NA
  Y[11:15, "t1"] <- NA
  testthat::expect_equal(sum(is.na(Y)), 15L)

  M <- 5000L
  cd <- cond_draw(Y, pars, tree, M)
  testthat::expect_length(cd$mis_idx, 15L)

  obs_ok <- vapply(cd$draws, function(m) isTRUE(all.equal(m[!is.na(Y)], Y[!is.na(Y)])), logical(1))
  testthat::expect_true(all(obs_ok))

  draws_mat <- vapply(cd$draws, function(m) as.vector(m)[cd$mis_idx], numeric(length(cd$mis_idx)))
  emp_mean <- rowMeans(draws_mat)
  emp_cov <- stats::cov(t(draws_mat))

  max_mean_diff <- max(abs(emp_mean - cd$cond_mean))
  max_cov_diff <- max(abs(emp_cov - cd$cond_cov))
  cat(sprintf("G-S2b max abs diff: mean %.4f, cov %.4f (n_mis = %d, M = %d)\n",
              max_mean_diff, max_cov_diff, length(cd$mis_idx), M))

  testthat::expect_lt(max_cov_diff, 0.05)
})
