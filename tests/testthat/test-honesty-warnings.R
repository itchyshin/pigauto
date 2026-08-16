# Silent-fallback honesty warnings (arc/silent-fallback-honesty).
#
# 1. P1-8: user covariates never reach the joint MVN / threshold-joint
#    baselines -- warn instead of silently ignoring them.
# 2. Small-validation conformal ceiling: at n_val < 19 the 95% split-
#    conformal target is arithmetically unreachable (ceiling is
#    n_val / (n_val + 1)); the warning must say so, not just "noisy".

test_that("[P1-8] joint MVN baseline warns that covariates are ignored", {
  skip_if_not_installed("Rphylopars")
  set.seed(42L)
  n <- 24L
  tree <- ape::rcoal(n)
  # >= 2 continuous traits => use_continuous_joint fires (single-obs,
  # no multi_proportion, default lambda_mode).
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   row.names = tree$tip.label)
  df$x1[c(2L, 9L)] <- NA
  df$x2[c(4L, 11L)] <- NA
  covs <- data.frame(env = stats::rnorm(n), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, covariates = covs)

  expect_warning(
    fit_baseline(pd, tree),
    regexp = "covariates.*ignored by the joint"
  )
})

test_that("[P1-8] joint MVN baseline is silent without covariates", {
  skip_if_not_installed("Rphylopars")
  set.seed(43L)
  n <- 24L
  tree <- ape::rcoal(n)
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   row.names = tree$tip.label)
  df$x1[c(2L, 9L)] <- NA
  df$x2[c(4L, 11L)] <- NA
  pd <- preprocess_traits(df, tree)

  expect_no_warning(fit_baseline(pd, tree))
})

test_that("[conformal-ceiling] calibrate_gates names the n_val < 19 ceiling", {
  set.seed(2033L)
  n <- 30L
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L,
                  mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  truth <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  val   <- matrix(FALSE, n, 1L)
  val[1:12, 1L] <- TRUE   # 12 val cells: above min_val_cells (10), below 19

  # 12 >= min_val_cells so the "Small validation set" branch does NOT
  # fire; the unreachable-ceiling branch must fire on its own.
  expect_warning(
    pigauto:::calibrate_gates(
      trait_map = tm, mu_cal = mu, delta_cal = delta,
      X_truth_r = truth, val_mask_mat = val,
      gate_grid = seq(0, 1, 0.1), gate_cap = 1,
      safety_floor = TRUE,
      mean_baseline_per_col = c(x = 0),
      simplex_step = 0.1,
      latent_names = "x", verbose = FALSE, seed = 2033L),
    regexp = "NOT achievable.*n_val / \\(n_val \\+ 1\\)"
  )
})

test_that("[conformal-ceiling] no ceiling warning at n_val >= 19", {
  set.seed(2034L)
  n <- 60L
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L,
                  mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  truth <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  val   <- matrix(FALSE, n, 1L)
  val[1:20, 1L] <- TRUE   # 20 val cells: ceiling 20/21 > 0.95

  expect_no_warning(
    pigauto:::calibrate_gates(
      trait_map = tm, mu_cal = mu, delta_cal = delta,
      X_truth_r = truth, val_mask_mat = val,
      gate_grid = seq(0, 1, 0.1), gate_cap = 1,
      safety_floor = TRUE,
      mean_baseline_per_col = c(x = 0),
      simplex_step = 0.1,
      latent_names = "x", verbose = FALSE, seed = 2034L)
  )
})

test_that("[conformal-ceiling] binary traits do not trigger the ceiling warning", {
  # Conformal scores are not computed for binary/categorical, so the
  # ceiling message would be misleading there.
  set.seed(2035L)
  n <- 30L
  tm <- list(list(name = "b", type = "binary", latent_cols = 1L,
                  levels = c("no", "yes"), mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "b"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "b"))
  truth <- matrix(rbinom(n, 1L, 0.5), n, 1L, dimnames = list(NULL, "b"))
  val   <- matrix(FALSE, n, 1L)
  val[1:12, 1L] <- TRUE   # 12 val cells, binary

  expect_no_warning(
    pigauto:::calibrate_gates(
      trait_map = tm, mu_cal = mu, delta_cal = delta,
      X_truth_r = truth, val_mask_mat = val,
      gate_grid = seq(0, 1, 0.1), gate_cap = 1,
      safety_floor = TRUE,
      mean_baseline_per_col = c(b = 0),
      simplex_step = 0.1,
      latent_names = "b", verbose = FALSE, seed = 2035L)
  )
})
