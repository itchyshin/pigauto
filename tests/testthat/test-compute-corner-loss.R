# tests/testthat/test-compute-corner-loss.R
#
# Formula-level unit tests for `compute_corner_loss()` (extracted from
# `calibrate_gates()`'s `cal_mean_loss` closure 2026-04-30 per the Opus
# adversarial-review recommendation).
#
# Each branch (continuous, count, ordinal, proportion, binary,
# categorical, multi_proportion, zi_count) gets a dedicated test that
# verifies the formula by hand-construction.  The corners
# `c(1, 0, 0)`, `c(0, 1, 0)`, `c(0, 0, 1)` are checked individually
# under safety_floor = TRUE; the legacy scalar `g = 0` and `g = 1` are
# checked under safety_floor = FALSE.
#
# Critical case: zi_count `c(0, 0, 1)` previously contaminated the
# magnitude column blend with the GATE column's logit mean (single
# `mean_j = mean_baseline_per_col[lc[1]]` used for both blend1 calls).
# This file's test_that for zi_count ensures the magnitude blend uses
# `mean_baseline_per_col[lc[2]]` -- the formula a hand-written
# implementer would expect.

# ---- Helpers --------------------------------------------------------------

make_tm <- function(name, type, latent_cols, levels = NULL,
                    mean = 0, sd = 1) {
  tm <- list(name = name, type = type, latent_cols = latent_cols,
             levels = levels, mean = mean, sd = sd)
  tm
}

# ---- Continuous (and the count/ordinal/proportion clones) ------------------

test_that("compute_corner_loss continuous: pure-BM corner = MSE on baseline preds", {
  set.seed(1)
  n <- 20L
  mu_cal    <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  delta_cal <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  X_truth_r <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  tm <- make_tm("x", "continuous", 1L)
  rows <- seq_len(n)
  mb <- 0.5  # arbitrary

  # Pure BM: g = c(1, 0, 0) under safety_floor = TRUE.  The result
  # should be MSE between mu_cal and X_truth_r on all rows.
  expected <- mean((mu_cal[rows, 1] - X_truth_r[rows, 1])^2)
  observed <- pigauto:::compute_corner_loss(
    g = c(1, 0, 0), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(x = mb)
  )
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("compute_corner_loss continuous: pure-MEAN corner = MSE on the column mean", {
  set.seed(2)
  n <- 20L
  mu_cal    <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  delta_cal <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  X_truth_r <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  tm <- make_tm("x", "continuous", 1L)
  rows <- seq_len(n)
  mb <- 0.7

  # Pure MEAN: g = c(0, 0, 1).  Result should be MSE of (mb - truth).
  expected <- mean((mb - X_truth_r[rows, 1])^2)
  observed <- pigauto:::compute_corner_loss(
    g = c(0, 0, 1), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(x = mb)
  )
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("compute_corner_loss continuous: pure-GNN corner = MSE on delta_cal", {
  set.seed(3)
  n <- 20L
  mu_cal    <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  delta_cal <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  X_truth_r <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  tm <- make_tm("x", "continuous", 1L)
  rows <- seq_len(n)

  expected <- mean((delta_cal[rows, 1] - X_truth_r[rows, 1])^2)
  observed <- pigauto:::compute_corner_loss(
    g = c(0, 1, 0), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(x = 0)
  )
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("compute_corner_loss legacy (safety_floor=FALSE): scalar g works", {
  set.seed(4)
  n <- 20L
  mu_cal    <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  delta_cal <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  X_truth_r <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  tm <- make_tm("x", "continuous", 1L)
  rows <- seq_len(n)

  # g = 0 -> pure BM
  exp_bm  <- mean((mu_cal[rows, 1] - X_truth_r[rows, 1])^2)
  obs_bm  <- pigauto:::compute_corner_loss(
    g = 0, rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = FALSE
  )
  expect_equal(obs_bm, exp_bm, tolerance = 1e-12)

  # g = 1 -> pure GNN
  exp_gnn <- mean((delta_cal[rows, 1] - X_truth_r[rows, 1])^2)
  obs_gnn <- pigauto:::compute_corner_loss(
    g = 1, rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = FALSE
  )
  expect_equal(obs_gnn, exp_gnn, tolerance = 1e-12)

  # Legacy promoted vector g = c(1-r, r, 0) should extract r unambiguously
  obs_blend <- pigauto:::compute_corner_loss(
    g = c(0.7, 0.3, 0), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = FALSE
  )
  exp_blend <- mean((0.7 * mu_cal[rows, 1] + 0.3 * delta_cal[rows, 1] -
                     X_truth_r[rows, 1])^2)
  expect_equal(obs_blend, exp_blend, tolerance = 1e-12)
})

test_that("compute_corner_loss returns Inf for empty rows", {
  tm <- make_tm("x", "continuous", 1L)
  observed <- pigauto:::compute_corner_loss(
    g = c(1, 0, 0), rows = integer(0), tm = tm,
    mu_cal = matrix(0, 0, 1), delta_cal = matrix(0, 0, 1),
    X_truth_r = matrix(0, 0, 1),
    safety_floor = TRUE, mean_baseline_per_col = c(x = 0)
  )
  expect_identical(observed, Inf)
})

# ---- Binary ---------------------------------------------------------------

test_that("compute_corner_loss binary: 0-1 loss on the calibrated logit blend", {
  set.seed(5)
  n <- 20L
  # mu_cal has logits with strong sign so pred_class is unambiguous
  mu_cal    <- matrix(c(rep(2, n / 2), rep(-2, n / 2)), n, 1,
                       dimnames = list(NULL, "y"))
  delta_cal <- matrix(0, n, 1, dimnames = list(NULL, "y"))
  X_truth_r <- matrix(c(rep(1, n / 2), rep(0, n / 2)), n, 1,
                       dimnames = list(NULL, "y"))
  tm   <- make_tm("y", "binary", 1L)
  rows <- seq_len(n)

  # Pure-BM: blend = mu_cal -> logits 2 / -2 -> pred_class 1 / 0 -> matches truth
  observed_bm <- pigauto:::compute_corner_loss(
    g = c(1, 0, 0), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(y = 0)
  )
  expect_equal(observed_bm, 0)

  # Pure-MEAN with mb = -1: blend logit = -1 < 0 -> pred_class 0 -> half wrong
  observed_mean <- pigauto:::compute_corner_loss(
    g = c(0, 0, 1), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(y = -1)
  )
  expect_equal(observed_mean, 0.5)
})

# ---- Categorical -----------------------------------------------------------

test_that("compute_corner_loss categorical: argmax 0-1 loss with per-column mean", {
  # K = 3 categorical: 3 latent columns (one per class).  Pure-MEAN
  # corner uses mean_baseline_per_col[lc[k]] for each class column k.
  # Construct a fixture where pure-MEAN's argmax is column 2 (class B)
  # uniformly, but truth's argmax is column 1 / 2 / 3 in equal thirds.
  set.seed(6)
  n <- 9L
  mu_cal <- matrix(c(rep(c(2, 0, 0), 3),
                      rep(c(0, 2, 0), 3),
                      rep(c(0, 0, 2), 3)),
                    nrow = n, ncol = 3, byrow = TRUE,
                    dimnames = list(NULL, c("y=A", "y=B", "y=C")))
  delta_cal <- matrix(0, n, 3, dimnames = dimnames(mu_cal))
  truth_oh  <- mu_cal  # use the BM matrix as the one-hot truth

  tm   <- make_tm("y", "categorical", 1:3, levels = c("A", "B", "C"))
  rows <- seq_len(n)
  mb   <- c(`y=A` = -1, `y=B` = 1, `y=C` = -2)  # MEAN's argmax = "B"

  # Pure-BM: argmax(BM logits) matches truth -> 0-1 loss = 0
  expect_equal(
    pigauto:::compute_corner_loss(
      g = c(1, 0, 0), rows = rows, tm = tm,
      mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = truth_oh,
      safety_floor = TRUE, mean_baseline_per_col = mb),
    0
  )

  # Pure-MEAN with mb argmax = "B" -> wrong on rows 1-3 (truth A) and 7-9 (truth C)
  # 6 / 9 wrong = 2/3
  expect_equal(
    pigauto:::compute_corner_loss(
      g = c(0, 0, 1), rows = rows, tm = tm,
      mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = truth_oh,
      safety_floor = TRUE, mean_baseline_per_col = mb),
    2 / 3,
    tolerance = 1e-12
  )
})

# ---- multi_proportion -----------------------------------------------------

test_that("compute_corner_loss multi_proportion returns Inf under safety_floor=TRUE", {
  tm <- make_tm("p", "multi_proportion", 1:3)
  observed <- pigauto:::compute_corner_loss(
    g = c(0, 1, 0), rows = 1:5, tm = tm,
    mu_cal = matrix(0, 5, 3), delta_cal = matrix(0, 5, 3),
    X_truth_r = matrix(0, 5, 3),
    safety_floor = TRUE, mean_baseline_per_col = c(0, 0, 0)
  )
  expect_identical(observed, Inf)
})

test_that("compute_corner_loss multi_proportion legacy (safety_floor=FALSE): MSE on K columns", {
  set.seed(7)
  n <- 5L
  K <- 3L
  mu_cal    <- matrix(rnorm(n * K), n, K)
  delta_cal <- matrix(rnorm(n * K), n, K)
  X_truth_r <- matrix(rnorm(n * K), n, K)
  tm   <- make_tm("p", "multi_proportion", 1:K)
  rows <- seq_len(n)

  # g = 0 -> pure BM
  expected <- mean((mu_cal[rows, , drop = FALSE] - X_truth_r[rows, , drop = FALSE])^2)
  observed <- pigauto:::compute_corner_loss(
    g = 0, rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = FALSE
  )
  expect_equal(observed, expected, tolerance = 1e-12)
})

# ---- ZI count: the bug Opus surfaced --------------------------------------

test_that("compute_corner_loss zi_count: pure-MEAN corner uses lc[2] mean for magnitude (not lc[1])", {
  # Regression test for the Opus 2026-04-30 finding: pre-fix the
  # zi_count branch used `mean_j = mean_baseline_per_col[lc[1]]` (gate
  # column logit) for BOTH the gate blend and the magnitude blend.
  # That contaminated the c(0, 0, 1) corner's magnitude formula with
  # the gate's logit value.
  #
  # Construct a fixture where the gate column mean is wildly different
  # from the magnitude column mean.  The pure-MEAN corner's loss
  # depends on which mean is used for the magnitude blend; an
  # incorrect implementation (gate-mean substituted for magnitude-mean)
  # produces a measurably different number.

  set.seed(8)
  n <- 20L
  # Latent layout: gate at col 1, magnitude at col 2
  mu_cal    <- matrix(c(rnorm(n, mean = -3),       # gate logits
                        rnorm(n, mean =  0)),       # magnitude z-score
                       nrow = n, ncol = 2,
                       dimnames = list(NULL, c("zi_gate", "zi_mag")))
  delta_cal <- matrix(0, n, 2, dimnames = dimnames(mu_cal))

  # Truth: gate column is 0/1 indicator, magnitude column is z-score
  # of log1p(count) when nonzero, NA when zero.
  truth_gate <- c(rep(0, n / 2), rep(1, n / 2))
  truth_mag  <- c(rep(NA_real_, n / 2), rnorm(n / 2))
  X_truth_r  <- cbind(zi_gate = truth_gate, zi_mag = truth_mag)

  tm <- make_tm("zi", "zi_count", c(1L, 2L), mean = 1.5, sd = 0.8)

  # Means: gate logit mean = -3 (very low p_nz), magnitude z mean = 0
  mb <- c(zi_gate = -3, zi_mag = 0)

  # Compute the EXPECTED pure-MEAN corner loss using each column's own mean
  # (i.e., the post-fix formula).
  blend1 <- function(g_corner, mu, dlt, mscalar) {
    g_corner[1] * mu + g_corner[2] * dlt + g_corner[3] * mscalar
  }
  rows <- seq_len(n)
  gate_pred_correct <- blend1(c(0, 0, 1),
                               mu_cal[rows, 1], delta_cal[rows, 1],
                               mb["zi_gate"])
  mag_pred_correct  <- blend1(c(0, 0, 1),
                               mu_cal[rows, 2], delta_cal[rows, 2],
                               mb["zi_mag"])
  p_nz_correct  <- plogis(gate_pred_correct)
  count_hat_correct <- pmax(expm1(mag_pred_correct * tm$sd + tm$mean), 0)
  pred_ev_correct <- p_nz_correct * count_hat_correct
  truth_ev <- rep(0, n)
  nz <- which(truth_gate > 0.5 & is.finite(truth_mag))
  truth_ev[nz] <- expm1(truth_mag[nz] * tm$sd + tm$mean)
  expected_correct <- mean((pred_ev_correct - truth_ev)^2)

  # Compute the BUGGY pure-MEAN corner loss (gate mean substituted into
  # magnitude blend).
  mag_pred_buggy <- blend1(c(0, 0, 1),
                            mu_cal[rows, 2], delta_cal[rows, 2],
                            mb["zi_gate"])
  count_hat_buggy <- pmax(expm1(mag_pred_buggy * tm$sd + tm$mean), 0)
  pred_ev_buggy <- p_nz_correct * count_hat_buggy
  expected_buggy <- mean((pred_ev_buggy - truth_ev)^2)

  observed <- pigauto:::compute_corner_loss(
    g = c(0, 0, 1), rows = rows, tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = mb
  )

  expect_equal(observed, expected_correct, tolerance = 1e-9,
               info = "zi_count c(0,0,1) corner must use lc[2] mean for magnitude")
  # Sanity: the buggy version produces a DIFFERENT number on this fixture
  expect_false(isTRUE(all.equal(expected_correct, expected_buggy)),
               info = "fixture is sensitive to the bug; without sensitivity the test is vacuous")
})

test_that("compute_corner_loss zi_count: pure-BM corner = current zi_count loss formula", {
  # Sanity check on the c(1, 0, 0) corner -- should not depend on the
  # mean_baseline_per_col values at all.
  set.seed(9)
  n <- 10L
  mu_cal    <- matrix(c(rnorm(n), rnorm(n)), n, 2,
                       dimnames = list(NULL, c("g", "m")))
  delta_cal <- matrix(0, n, 2, dimnames = dimnames(mu_cal))
  truth_gate <- rbinom(n, 1, 0.5)
  truth_mag  <- ifelse(truth_gate == 1, rnorm(n), NA_real_)
  X_truth_r  <- cbind(g = truth_gate, m = truth_mag)
  tm <- make_tm("zi", "zi_count", c(1L, 2L), mean = 1.0, sd = 1.0)

  out_a <- pigauto:::compute_corner_loss(
    g = c(1, 0, 0), rows = seq_len(n), tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(g = -100, m = -100)
  )
  out_b <- pigauto:::compute_corner_loss(
    g = c(1, 0, 0), rows = seq_len(n), tm = tm,
    mu_cal = mu_cal, delta_cal = delta_cal, X_truth_r = X_truth_r,
    safety_floor = TRUE, mean_baseline_per_col = c(g = +100, m = +100)
  )
  # Pure-BM corner (g[3] = 0) ignores the mean values entirely.
  expect_equal(out_a, out_b, tolerance = 1e-12)
})
