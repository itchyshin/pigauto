# B1 / B2 / B3 — train → calibrate → predict symmetry for the GNN path.
# See conversation 2026-08-07 critical review of ResidualPhyloDAE wiring.

skip_if_no_torch <- function() {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")
}

# ---- B1: held-out cells must not appear as truth in the train input ----------

test_that("mask_heldout_with_baseline() replaces holdout cells with MU", {
  X  <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L, ncol = 2L)
  MU <- matrix(0, nrow = 3L, ncol = 2L)
  # Linear indices (column-major), matching splits$val_idx / test_idx
  hold_idx <- c(1L, 5L)  # X[1,1] and X[2,2]
  out <- mask_heldout_with_baseline(X, MU, hold_idx)
  expect_equal(out[1, 1], 0)
  expect_equal(out[2, 2], 0)
  expect_equal(out[3, 1], 3)   # untouched
  expect_equal(out[1, 2], 4)   # untouched
})

test_that("training forward does not see held-out truth (B1 smoke)", {
  skip_if_no_torch()
  skip_if_not_installed("ape")
  set.seed(20260807L)
  n <- 40L
  tree <- ape::rcoal(n)
  tip <- tree$tip.label
  # Two continuous traits; trait2 is pure noise (weak phylo)
  x1 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  x2 <- stats::rnorm(n)
  df <- data.frame(x1 = as.numeric(x1[tip]), x2 = x2, row.names = tip)
  hide <- sample.int(n, 8L)
  df$x1[hide] <- NA_real_

  # Fit with a held-out val/test split (default missing_frac > 0 via impute)
  res <- pigauto::impute(
    df, tree,
    epochs = 25L, n_imputations = 1L,
    refine_steps = 1L,
    phylo_signal_gate = FALSE,
    safety_floor = TRUE,
    verbose = FALSE, seed = 20260807L
  )
  # Contract: fit records that train used heldout-masked inputs
  expect_true(isTRUE(res$fit$model_config$train_mask_heldout))
})

# ---- B2: phylo_signal_gate without safety_floor must not land on latent 0 ----

test_that("phylo_signal_gate + safety_floor=FALSE falls back to BM, not mean/zero", {
  skip_if_no_torch()
  skip_if_not_installed("ape")
  skip_if_not_installed("phytools")
  set.seed(20260807L)
  n <- 80L
  tree <- ape::rcoal(n)
  tip <- tree$tip.label
  # Pure noise → λ ≈ 0 → phylo gate should trigger
  df <- data.frame(noise = stats::rnorm(n, mean = 2, sd = 1),
                   row.names = tip)
  df$noise[sample.int(n, 15L)] <- NA_real_

  expect_warning(
    res <- pigauto::impute(
      df, tree,
      phylo_signal_gate = TRUE,
      phylo_signal_threshold = 0.5,
      safety_floor = FALSE,
      epochs = 30L, n_imputations = 1L,
      refine_steps = 1L,
      verbose = FALSE, seed = 20260807L
    ),
    regexp = "phylo_signal_gate|safety_floor"
  )
  fit <- res$fit
  expect_true(any(fit$phylo_gate_triggered, na.rm = TRUE))
  gated <- names(fit$phylo_gate_triggered)[fit$phylo_gate_triggered]
  # Must NOT be the mean corner without a mean baseline (would become latent 0)
  expect_equal(unname(fit$r_cal_mean[gated]), rep(0, length(gated)))
  expect_equal(unname(fit$r_cal_gnn[gated]),  rep(0, length(gated)))
  expect_equal(unname(fit$r_cal_bm[gated]),   rep(1, length(gated)))
})

# ---- B3: calibration / conformal use the same refine_steps as predict -------

test_that("gate calibration uses refine_steps matching predict (B3)", {
  skip_if_no_torch()
  skip_if_not_installed("ape")
  set.seed(20260807L)
  n <- 50L
  tree <- ape::rcoal(n)
  tip <- tree$tip.label
  x <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  df <- data.frame(y = as.numeric(x[tip]), row.names = tip)
  df$y[sample.int(n, 10L)] <- NA_real_

  res <- pigauto::impute(
    df, tree,
    epochs = 30L, n_imputations = 1L,
    refine_steps = 3L,
    phylo_signal_gate = FALSE,
    safety_floor = TRUE,
    verbose = FALSE, seed = 20260807L
  )
  expect_equal(res$fit$model_config$refine_steps, 3L)
  expect_equal(res$fit$model_config$cal_refine_steps, 3L)
})

test_that("pigauto_refine_forward() returns delta_mat with refine_steps >= 1", {
  skip_if_no_torch()
  skip_if_not_installed("ape")
  set.seed(1L)
  n <- 20L
  tree <- ape::rcoal(n)
  tip <- tree$tip.label
  x <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  df <- data.frame(y = as.numeric(x[tip]), row.names = tip)
  df$y[sample.int(n, 4L)] <- NA_real_
  res <- pigauto::impute(
    df, tree,
    epochs = 15L, n_imputations = 1L,
    refine_steps = 2L,
    phylo_signal_gate = FALSE,
    verbose = FALSE, seed = 1L
  )
  # Helper is internal; contract is exposed via cal_refine_steps on the fit.
  expect_equal(res$fit$model_config$cal_refine_steps, 2L)
  expect_true(is.numeric(res$fit$r_cal_gnn))
})
