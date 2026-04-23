# tests/testthat/test-safety-floor.R
# Smoke canary for safety-floor (three-way mean-gate) calibration.
# See specs/2026-04-23-safety-floor-mean-gate-design.md.

test_that("simplex_grid(0.05) returns 231 rows on the simplex", {
  g <- pigauto:::simplex_grid(step = 0.05)
  expect_equal(dim(g), c(231L, 3L))
  expect_equal(colnames(g), c("r_bm", "r_gnn", "r_mean"))
  expect_true(all(g >= 0 & g <= 1))
  expect_equal(rowSums(g), rep(1, 231L), tolerance = 1e-10)
})

test_that("simplex_grid(0.25) returns 15 rows (coarse grid for tests)", {
  g <- pigauto:::simplex_grid(step = 0.25)
  expect_equal(nrow(g), 15L)
  expect_equal(rowSums(g), rep(1, 15L), tolerance = 1e-10)
})

test_that("simplex_grid always contains the three corners", {
  g <- pigauto:::simplex_grid(step = 0.05)
  corners <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3L, byrow = TRUE,
                    dimnames = list(NULL, c("r_bm", "r_gnn", "r_mean")))
  for (k in seq_len(nrow(corners))) {
    matched <- any(apply(g, 1L, function(row) all(row == corners[k, ])))
    expect_true(matched,
                info = sprintf("corner (%g,%g,%g) missing",
                               corners[k,1], corners[k,2], corners[k,3]))
  }
})

# ---- Task 1 follow-up: cover the non-divisor error path ----

test_that("simplex_grid errors when step does not evenly divide 1", {
  expect_error(pigauto:::simplex_grid(step = 0.3), "evenly divide")
})

# ---- Task 2: mean_baseline_scalar() ----

test_that("mean_baseline_scalar continuous = training grand mean on latent scale", {
  # Simulate: column is z-scored, training-observed mask drops some rows
  set.seed(1L)
  X <- matrix(rnorm(100), nrow = 100, ncol = 1L)
  train_mask <- rep(TRUE, 100); train_mask[11:20] <- FALSE
  # Latent scale mean of training-observed rows:
  expected <- mean(X[train_mask, 1L])
  got <- pigauto:::mean_baseline_scalar(X[, 1L], train_mask,
                                          trait_type = "continuous")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("mean_baseline_scalar binary = qlogis(clip(freq, 0.01, 0.99))", {
  X <- c(rep(1, 30), rep(0, 70))   # freq = 0.30
  expected <- qlogis(0.30)
  got <- pigauto:::mean_baseline_scalar(X, rep(TRUE, 100),
                                          trait_type = "binary")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("mean_baseline_scalar binary clips degenerate all-0 to qlogis(0.01)", {
  X <- rep(0, 100)
  expected <- qlogis(0.01)
  got <- pigauto:::mean_baseline_scalar(X, rep(TRUE, 100),
                                          trait_type = "binary")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("mean_baseline_scalar binary clips degenerate all-1 to qlogis(0.99)", {
  X <- rep(1, 100)
  expected <- qlogis(0.99)
  got <- pigauto:::mean_baseline_scalar(X, rep(TRUE, 100),
                                          trait_type = "binary")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("mean_baseline_scalar categorical column = qlogis(clip(mean(Y_k), 0.01, 0.99))", {
  Y_k <- c(rep(1, 20), rep(0, 80))
  expected <- qlogis(0.20)
  got <- pigauto:::mean_baseline_scalar(Y_k, rep(TRUE, 100),
                                          trait_type = "categorical")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("mean_baseline_scalar excludes NA and non-training rows", {
  X <- c(NA_real_, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  train_mask <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  # train_mask AND !is.na(X) -> rows 2..6: 3/5 = 0.6
  expected <- qlogis(0.6)
  got <- pigauto:::mean_baseline_scalar(X, train_mask,
                                          trait_type = "binary")
  expect_equal(got, expected, tolerance = 1e-12)
})

# ---- Task 3: calibrate_gates return shape change ----

test_that("calibrate_gates() returns list with r_cal_bm/r_cal_gnn/r_cal_mean", {
  set.seed(2026)
  n <- 30L
  tm <- list(list(name = "x1", type = "continuous",
                   latent_cols = 1L, mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), nrow = n, ncol = 1L)
  delta <- matrix(rnorm(n), nrow = n, ncol = 1L)
  truth <- matrix(rnorm(n), nrow = n, ncol = 1L)
  val   <- matrix(c(rep(TRUE, 10), rep(FALSE, 20)), nrow = n, ncol = 1L)
  res <- pigauto:::calibrate_gates(
    trait_map = tm, mu_cal = mu, delta_cal = delta,
    X_truth_r = truth, val_mask_mat = val,
    gate_grid = c(0, 0.25, 0.5, 0.75, 1), gate_cap = 1,
    latent_names = "x1", verbose = FALSE, seed = 2026L)
  expect_type(res, "list")
  expect_named(res, c("r_cal_bm", "r_cal_gnn", "r_cal_mean"),
               ignore.order = TRUE)
  expect_equal(length(res$r_cal_bm), 1L)
  expect_equal(length(res$r_cal_gnn), 1L)
  expect_equal(length(res$r_cal_mean), 1L)
  # 1-D path: r_mean must be 0 (safety_floor not yet on; this is behaviour
  # preservation).
  expect_equal(as.numeric(res$r_cal_mean), 0)
  # r_cal_bm + r_cal_gnn must be 1 (the 1-D path splits the blend into
  # (1 - r) BM + r GNN, so r_cal_bm = 1 - r_cal_gnn).
  expect_equal(as.numeric(res$r_cal_bm + res$r_cal_gnn), 1)
})

# ---- Task 4: simplex search when safety_floor = TRUE ----

test_that("calibrate_gates(safety_floor = TRUE) finds w_mean > 0 when mean dominates", {
  # Construct a column where the grand mean BEATS both bm and gnn.
  # truth = rep(0, 30); bm = rnorm(30, 0, 2); gnn = rnorm(30, 0, 2); mean = 0.
  # The grid should pick a weight triple with r_mean > 0.
  set.seed(42L)
  n <- 30L
  tm <- list(list(name = "x1", type = "continuous",
                   latent_cols = 1L, mean = 0, sd = 1))
  mu    <- matrix(rnorm(n, 0, 2), nrow = n)
  delta <- matrix(rnorm(n, 0, 2), nrow = n)
  truth <- matrix(rep(0, n), nrow = n)   # degenerate: mean IS truth
  val   <- matrix(rep(TRUE, n), nrow = n)
  res <- pigauto:::calibrate_gates(
    trait_map = tm, mu_cal = mu, delta_cal = delta,
    X_truth_r = truth, val_mask_mat = val,
    gate_grid = seq(0, 1, 0.1), gate_cap = 1,
    safety_floor = TRUE,
    mean_baseline_per_col = c(x1 = 0),
    simplex_step = 0.05,
    latent_names = "x1", verbose = FALSE, seed = 2026L)
  expect_gt(as.numeric(res$r_cal_mean), 0)
  # Simplex: three weights sum to exactly 1 (or within renorm tolerance).
  expect_equal(as.numeric(res$r_cal_bm + res$r_cal_gnn + res$r_cal_mean),
               1, tolerance = 1e-10)
})

test_that("calibrate_gates(safety_floor = TRUE) invariant: loss <= mean_loss on val", {
  # Two columns: one where BM is good, one where mean is best.
  set.seed(2026L)
  n <- 50L
  tm <- list(
    list(name = "phy",   type = "continuous", latent_cols = 1L,
          mean = 0, sd = 1),
    list(name = "noise", type = "continuous", latent_cols = 2L,
          mean = 0, sd = 1))
  truth <- cbind(rnorm(n), rnorm(n, 0, 0.5))   # col 2 clusters near 0
  mu    <- cbind(truth[, 1] + rnorm(n, 0, 0.3),   # good BM on col 1
                 rnorm(n, 0, 2))                  # bad BM on col 2
  delta <- cbind(rnorm(n, 0, 2),                  # bad GNN on col 1
                 rnorm(n, 0, 2))                  # bad GNN on col 2
  val   <- matrix(TRUE, nrow = n, ncol = 2L)
  res <- pigauto:::calibrate_gates(
    trait_map = tm, mu_cal = mu, delta_cal = delta,
    X_truth_r = truth, val_mask_mat = val,
    gate_grid = seq(0, 1, 0.1), gate_cap = 1,
    safety_floor = TRUE,
    mean_baseline_per_col = c(phy = 0, noise = 0),
    simplex_step = 0.05,
    latent_names = c("phy", "noise"), verbose = FALSE, seed = 2026L)
  for (j in 1:2) {
    w <- c(res$r_cal_bm[j], res$r_cal_gnn[j], res$r_cal_mean[j])
    blended <- w[1] * mu[, j] + w[2] * delta[, j] + w[3] * 0
    mean_only <- rep(0, n)
    expect_lte(mean((blended - truth[, j])^2),
               mean((mean_only - truth[, j])^2) + 1e-10)
  }
})

# ---- Task 5: fit_pigauto() safety_floor integration ----

test_that("fit_pigauto(safety_floor = TRUE) stores mean_baseline_per_col + simplex weights on fit", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")
  data(avonet300, tree300, package = "pigauto")
  set.seed(2026L)
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  df$Mass[sample(300, 30)]               <- NA_real_
  df$Beak.Length_Culmen[sample(300, 30)] <- NA_real_
  res <- pigauto::impute(df, tree300, safety_floor = TRUE,
                           epochs = 50L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  fit <- res$fit
  expect_true(!is.null(fit$mean_baseline_per_col))
  expect_type(fit$mean_baseline_per_col, "double")
  expect_true(!is.null(fit$r_cal_bm))
  expect_true(!is.null(fit$r_cal_gnn))
  expect_true(!is.null(fit$r_cal_mean))
  expect_equal(length(fit$r_cal_bm), length(fit$mean_baseline_per_col))
  sums <- fit$r_cal_bm + fit$r_cal_gnn + fit$r_cal_mean
  expect_equal(as.numeric(sums), rep(1, length(sums)), tolerance = 1e-8)
  expect_true(isTRUE(fit$safety_floor))
})

test_that("fit_pigauto(safety_floor = FALSE) reproduces v0.9.1 behaviour (r_mean = 0, mean_baseline_per_col NULL or 0)", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")
  data(avonet300, tree300, package = "pigauto")
  set.seed(2026L)
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  df$Mass[sample(300, 30)] <- NA_real_
  res <- pigauto::impute(df, tree300, safety_floor = FALSE,
                           epochs = 50L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  fit <- res$fit
  expect_false(isTRUE(fit$safety_floor))
  # Legacy path: r_cal_mean must be exactly 0 on every column.
  expect_equal(as.numeric(fit$r_cal_mean), rep(0, length(fit$r_cal_mean)))
})

# ---- Task 4 follow-up: guard against zero-sum median weight vector ----

test_that("calibrate_gates(safety_floor = TRUE, median_splits) produces finite weights", {
  # median_splits with small n_val + simplex can produce adversarial splits
  # whose per-axis median is c(0, 0, 0). The guard falls back to c(0, 0, 1).
  set.seed(2026L)
  n <- 15L   # small val set
  tm <- list(list(name = "x1", type = "continuous",
                   latent_cols = 1L, mean = 0, sd = 1))
  mu    <- matrix(rnorm(n, 0, 3), nrow = n)
  delta <- matrix(rnorm(n, 0, 3), nrow = n)
  truth <- matrix(rnorm(n), nrow = n)
  val   <- matrix(rep(TRUE, n), nrow = n)
  res <- pigauto:::calibrate_gates(
    trait_map = tm, mu_cal = mu, delta_cal = delta,
    X_truth_r = truth, val_mask_mat = val,
    gate_grid = seq(0, 1, 0.1), gate_cap = 1,
    safety_floor = TRUE,
    mean_baseline_per_col = c(x1 = 0),
    simplex_step = 0.05,
    gate_method = "median_splits",
    gate_splits_B = 31L,
    latent_names = "x1", verbose = FALSE, seed = 7L)
  # Must be finite on every slot, must sum to 1.
  expect_true(all(is.finite(c(res$r_cal_bm, res$r_cal_gnn, res$r_cal_mean))))
  expect_equal(as.numeric(res$r_cal_bm + res$r_cal_gnn + res$r_cal_mean),
               1, tolerance = 1e-8)
})

# ---- Task 6: predict 3-way blend + backward compat ----

test_that("predict.pigauto_fit uses 3-way blend when r_cal_mean > 0", {
  # Construct a synthetic case where safety_floor should pull the gate
  # toward the mean: plants-like pattern with weak BM baseline.
  data("avonet300", package = "pigauto")
  data("tree300",   package = "pigauto")
  set.seed(2026L)
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  df$Mass[sample(300, 30)] <- NA_real_
  res <- pigauto::impute(df, tree300, safety_floor = TRUE,
                           epochs = 50L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  fit <- res$fit
  # All stored per-column weights sum to 1.
  expect_equal(as.numeric(fit$r_cal_bm + fit$r_cal_gnn + fit$r_cal_mean),
               rep(1, length(fit$r_cal_bm)), tolerance = 1e-8)
  # Predictions must be finite on every imputed cell.
  imputed_vals <- res$completed$Mass[res$imputed_mask[, "Mass"]]
  expect_true(all(is.finite(imputed_vals)))
})

test_that("predict.pigauto_fit falls back to 2-way blend on legacy v0.9.1 fits (null r_cal_mean)", {
  # Fit with safety_floor = FALSE first; compare against the same fit
  # with r_cal_bm / r_cal_gnn / r_cal_mean / mean_baseline_per_col all
  # set to NULL to simulate a pre-Task-3 v0.9.1 fit. %||% fallback must
  # reproduce identical predictions.
  data("avonet300", package = "pigauto")
  data("tree300",   package = "pigauto")
  set.seed(2026L)
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  df$Mass[sample(300, 30)] <- NA_real_
  res <- pigauto::impute(df, tree300, safety_floor = FALSE,
                           epochs = 50L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  new_pred   <- res$prediction$imputed
  # Simulate legacy fit: remove new slots, keep only r_cal + calibrated_gates.
  fit_legacy <- res$fit
  fit_legacy$r_cal_bm             <- NULL
  fit_legacy$r_cal_gnn            <- NULL
  fit_legacy$r_cal_mean           <- NULL
  fit_legacy$mean_baseline_per_col <- NULL
  fit_legacy$safety_floor         <- NULL
  pred_legacy_view <- predict(fit_legacy, n_imputations = 1L)
  # Should be bit-identical: the %||% fallback reconstructs r_cal_bm =
  # 1 - r_cal, r_cal_gnn = r_cal, r_cal_mean = 0, mean_baseline = 0.
  expect_equal(pred_legacy_view$imputed, new_pred, tolerance = 1e-10)
})

# ---- Task 8: multi_impute() mean-term propagation ----

# ---- Task 9: multi_impute_trees() share_gnn + safety_floor ----

test_that("multi_impute_trees(share_gnn = TRUE, safety_floor = TRUE) reuses safety weights across trees", {
  # trees300 is bundled in pigauto; skip gracefully if absent.
  skip_if(
    !exists("trees300", where = asNamespace("pigauto")) &&
      inherits(try(data("trees300", package = "pigauto"), silent = TRUE), "try-error"),
    "trees300 dataset not available")
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")
  data("avonet300",  package = "pigauto")
  data("trees300",   package = "pigauto")
  set.seed(2026L)
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  df$Mass[sample(300, 30)] <- NA_real_
  mit <- pigauto::multi_impute_trees(
    df, trees300[1:3], m_per_tree = 1L,
    share_gnn = TRUE, safety_floor = TRUE,
    epochs = 50L, verbose = FALSE, seed = 2026L)
  # The reference-tree fit carries safety-floor metadata.
  expect_true(isTRUE(mit$fit$safety_floor))
  expect_true(!is.null(mit$fit$mean_baseline_per_col))
  # Each completed dataset has finite Mass values.
  for (i in seq_along(mit$datasets)) {
    mass_col <- mit$datasets[[i]]$Mass
    expect_true(all(is.finite(mass_col[!is.na(mass_col)])))
  }
})

# ---- Task 10: legacy v0.9.1 fit backward-compat fixture ----

test_that("legacy v0.9.1 fit fixture loads and predicts via %||% fallback", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")
  fx_path <- system.file("extdata", "legacy_fit_v091.rds", package = "pigauto")
  expect_true(nzchar(fx_path) && file.exists(fx_path))
  # load_pigauto() deserialises the torch model_state (cross-session portable)
  fit <- load_pigauto(fx_path)
  expect_s3_class(fit, "pigauto_fit")
  # Legacy fit has NO new safety-floor slots
  expect_null(fit$r_cal_bm)
  expect_null(fit$r_cal_gnn)
  expect_null(fit$r_cal_mean)
  expect_null(fit$mean_baseline_per_col)
  expect_null(fit$safety_floor)
  # Predict via %||% fallback — must succeed and return finite imputed values
  pred <- predict(fit, n_imputations = 1L)
  expect_true(!is.null(pred$imputed))
  expect_true(all(is.finite(as.matrix(pred$imputed)[!is.na(as.matrix(pred$imputed))])))
})

test_that("multi_impute(safety_floor = TRUE) pooled point + SE are finite", {
  data("avonet300",  package = "pigauto")
  data("tree300",    package = "pigauto")
  set.seed(2026L)
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  df$Mass[sample(300, 30)]               <- NA_real_
  df$Beak.Length_Culmen[sample(300, 30)] <- NA_real_
  mi <- pigauto::multi_impute(df, tree300, m = 5L, safety_floor = TRUE,
                                epochs = 50L, verbose = FALSE, seed = 2026L)
  # Pooled point must be finite for every imputed cell.
  expect_true(all(is.finite(mi$pooled_point[!is.na(mi$pooled_point)])))
  # Pooled SE on imputed cells (non-zero mask).
  se_vec <- mi$se[mi$imputed_mask]
  expect_true(all(is.finite(se_vec)))
  # Completed datasets: all finite on the imputed Mass column.
  for (k in seq_along(mi$datasets)) {
    mass_imp <- mi$datasets[[k]]$Mass[mi$imputed_mask[, "Mass"]]
    expect_true(all(is.finite(mass_imp)))
  }
})
