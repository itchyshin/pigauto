# Regression test for the BACE-compat CI evaluator
# (script/gha/_bace_compat.R::.bace_eval_per_imputation).
#
# Bug history (2026-05-17): the evaluator scored each of the m stochastic
# imputation draws produced by .sample_conformal_draw() against the truth.
# For categorical traits, each draw is sampled proportional to the per-cell
# probability vector — so per-draw accuracy is bounded by E[p_max] across
# cells and is much lower than argmax-of-pooled-probability accuracy.
#
# Fix: pool_method = "auto" (new default) routes discrete-type accuracy
# through mi$pooled_point (which contains argmax of averaged probability
# from predict.pigauto_fit + pool_imputations) while keeping continuous
# traits on per-draw RMSE for Rubin-style noise propagation.

test_that("BACE-compat eval pools discrete-type accuracy via pooled_point", {
  skip_if_not_installed("ape")

  # Source the BACE-compat helpers from script/gha — they are not exported
  # functions of pigauto but live alongside the package source.
  bace_compat_path <- file.path("..", "..", "script", "gha", "_bace_compat.R")
  if (!file.exists(bace_compat_path)) {
    skip(sprintf("script/gha/_bace_compat.R not found at %s", bace_compat_path))
  }
  source(bace_compat_path, local = TRUE)

  # Synthetic truth + per-draw predictions where the model is consistent
  # (always predicts class A on the masked cells, which matches truth). The
  # m stochastic draws ALSO predict mostly A but sometimes class B or C
  # (mimicking categorical sampling noise). pooled_point holds the model's
  # argmax = always A.
  truth_long <- data.frame(
    species_tip = paste0("sp", 1:10),
    trait       = "cat1",
    true_value  = "A",
    stringsAsFactors = FALSE
  )
  trait_types_lookup <- c(cat1 = "categorical")

  # m=5 stochastic draws: 60% A, 40% B/C noise
  set.seed(2026L)
  make_draw <- function(p_A = 0.6) {
    classes <- sample(c("A", "B", "C"), 10L, replace = TRUE,
                       prob = c(p_A, (1 - p_A) / 2, (1 - p_A) / 2))
    data.frame(cat1 = factor(classes, levels = c("A", "B", "C")),
                row.names = paste0("sp", 1:10))
  }
  draws <- lapply(seq_len(5L), function(i) make_draw())

  # pooled_point: model's argmax = always A (correct for all 10)
  pooled <- data.frame(cat1 = factor(rep("A", 10L), levels = c("A", "B", "C")),
                        row.names = paste0("sp", 1:10))

  # --- OLD behaviour (per_draw): expected accuracy ~ 0.6
  ev_old <- .bace_eval_per_imputation(draws, truth_long, trait_types_lookup,
                                        t_fit_sec = 0,
                                        pool_method = "per_draw")
  acc_per_draw <- median(ev_old$accuracy, na.rm = TRUE)
  expect_lt(acc_per_draw, 0.85)   # stochastic noise drags it below 1.0
  expect_gt(acc_per_draw, 0.35)   # but still above random (1/3)

  # --- NEW behaviour (auto, the default): pooled accuracy = 1.0
  ev_new <- .bace_eval_per_imputation(draws, truth_long, trait_types_lookup,
                                        t_fit_sec = 0,
                                        pooled_point = pooled)
  expect_true(nrow(ev_new) >= 1L)
  # Discrete trait → one row with imputation_idx = 0 carrying the
  # pooled-point accuracy. Since pooled_point matches truth perfectly,
  # this should be 1.0.
  pooled_row <- ev_new[ev_new$trait == "cat1" & ev_new$imputation_idx == 0L, ,
                         drop = FALSE]
  expect_equal(nrow(pooled_row), 1L)
  expect_equal(pooled_row$accuracy, 1.0)

  # Sanity: per-draw rows for discrete trait should NOT be emitted under
  # auto mode (they were emitted under per_draw mode).
  per_draw_rows <- ev_new[ev_new$imputation_idx > 0L & ev_new$type == "categorical", ,
                            drop = FALSE]
  expect_equal(nrow(per_draw_rows), 0L)
})

test_that("BACE-compat eval computes coverage_95 + interval_width from conformal bounds", {
  skip_if_not_installed("ape")
  bace_compat_path <- file.path("..", "..", "script", "gha", "_bace_compat.R")
  if (!file.exists(bace_compat_path)) {
    skip(sprintf("script/gha/_bace_compat.R not found at %s", bace_compat_path))
  }
  source(bace_compat_path, local = TRUE)

  truth_long <- data.frame(
    species_tip = paste0("sp", 1:10),
    trait       = "cont1",
    true_value  = as.character(seq_len(10L)),
    stringsAsFactors = FALSE
  )
  trait_types_lookup <- c(cont1 = "continuous")

  set.seed(2026L)
  draws <- list(data.frame(cont1 = (1:10) + rnorm(10L, 0, 0.1),
                            row.names = paste0("sp", 1:10)))
  pooled <- data.frame(cont1 = (1:10) + rnorm(10L, 0, 0.1),
                         row.names = paste0("sp", 1:10))

  # Conformal bounds: tight intervals covering 9 of 10 truth values
  truth_vals <- 1:10
  pred_vals <- pooled$cont1
  half_width <- abs(pred_vals - truth_vals) + c(rep(0.1, 9), -0.5)   # last cell uncovered
  conformal_lower <- matrix(pred_vals - half_width, ncol = 1L,
                              dimnames = list(paste0("sp", 1:10), "cont1"))
  conformal_upper <- matrix(pred_vals + half_width, ncol = 1L,
                              dimnames = list(paste0("sp", 1:10), "cont1"))

  ev <- .bace_eval_per_imputation(draws, truth_long, trait_types_lookup,
                                    t_fit_sec = 0,
                                    pooled_point    = pooled,
                                    conformal_lower = conformal_lower,
                                    conformal_upper = conformal_upper)
  cont_row <- ev[ev$type == "continuous", , drop = FALSE]
  expect_equal(nrow(cont_row), 1L)
  # 9 of 10 truths inside their conformal interval
  expect_equal(cont_row$coverage_95, 0.9)
  expect_true(is.finite(cont_row$interval_width))
  expect_true(cont_row$interval_width > 0)
})

test_that("BACE-compat eval computes Brier score from per-class probabilities", {
  skip_if_not_installed("ape")
  bace_compat_path <- file.path("..", "..", "script", "gha", "_bace_compat.R")
  if (!file.exists(bace_compat_path)) {
    skip(sprintf("script/gha/_bace_compat.R not found at %s", bace_compat_path))
  }
  source(bace_compat_path, local = TRUE)

  truth_long <- data.frame(
    species_tip = paste0("sp", 1:6),
    trait       = "cat1",
    true_value  = c("A","A","A","B","B","C"),
    stringsAsFactors = FALSE
  )
  trait_types_lookup <- c(cat1 = "categorical")

  # Perfect probabilities: 1.0 on the true class, 0 elsewhere → Brier = 0
  probs_perfect <- matrix(0, nrow = 6L, ncol = 3L,
                            dimnames = list(paste0("sp", 1:6),
                                            c("A","B","C")))
  for (i in seq_along(truth_long$true_value)) {
    probs_perfect[i, truth_long$true_value[i]] <- 1
  }
  pooled <- data.frame(cat1 = factor(truth_long$true_value,
                                       levels = c("A","B","C")),
                         row.names = paste0("sp", 1:6))

  ev <- .bace_eval_per_imputation(list(pooled), truth_long, trait_types_lookup,
                                    t_fit_sec = 0,
                                    pooled_point  = pooled,
                                    probabilities = list(cat1 = probs_perfect),
                                    trait_levels  = list(cat1 = c("A","B","C")))
  cat_row <- ev[ev$type == "categorical", , drop = FALSE]
  expect_equal(cat_row$accuracy, 1.0)
  expect_equal(cat_row$brier, 0.0, tolerance = 1e-10)

  # Uniform probabilities (all 1/3): Brier should be > 0.
  probs_uniform <- matrix(1/3, nrow = 6L, ncol = 3L,
                            dimnames = dimnames(probs_perfect))
  ev2 <- .bace_eval_per_imputation(list(pooled), truth_long, trait_types_lookup,
                                     t_fit_sec = 0,
                                     pooled_point  = pooled,
                                     probabilities = list(cat1 = probs_uniform),
                                     trait_levels  = list(cat1 = c("A","B","C")))
  cat_row2 <- ev2[ev2$type == "categorical", , drop = FALSE]
  expect_gt(cat_row2$brier, 0.3)
})

test_that("BACE-compat eval pools continuous RMSE via pooled_point under auto mode", {
  skip_if_not_installed("ape")
  bace_compat_path <- file.path("..", "..", "script", "gha", "_bace_compat.R")
  if (!file.exists(bace_compat_path)) {
    skip(sprintf("script/gha/_bace_compat.R not found at %s", bace_compat_path))
  }
  source(bace_compat_path, local = TRUE)

  # Continuous truth + 3 draws with mild Gaussian noise around truth.
  # pooled_point is tighter to truth than any individual draw -- that's
  # the whole reason BACE reports RMSE-of-mean rather than mean-of-RMSE.
  truth_long <- data.frame(
    species_tip = paste0("sp", 1:10),
    trait       = "cont1",
    true_value  = as.character(seq_len(10L)),
    stringsAsFactors = FALSE
  )
  trait_types_lookup <- c(cont1 = "continuous")

  set.seed(2026L)
  make_draw <- function() {
    data.frame(cont1 = (1:10) + rnorm(10L, 0, 0.5),
                row.names = paste0("sp", 1:10))
  }
  draws <- lapply(seq_len(3L), function(i) make_draw())
  pooled <- data.frame(cont1 = (1:10) + rnorm(10L, 0, 0.1),
                         row.names = paste0("sp", 1:10))

  # ---- AUTO: continuous goes through pooled_point (BACE-equivalent)
  ev_auto <- .bace_eval_per_imputation(draws, truth_long, trait_types_lookup,
                                         t_fit_sec = 0,
                                         pooled_point = pooled)
  cont_rows <- ev_auto[ev_auto$type == "continuous", , drop = FALSE]
  expect_equal(nrow(cont_rows), 1L)
  expect_equal(cont_rows$imputation_idx, 0L)
  expect_true(!is.na(cont_rows$rmse))
  expect_true(is.na(cont_rows$accuracy))

  # The pooled RMSE should equal the direct RMSE-on-pooled-point.
  direct_rmse <- sqrt(mean((pooled$cont1 - (1:10))^2))
  expect_equal(cont_rows$rmse, direct_rmse, tolerance = 1e-10)

  # ---- PER_DRAW: 3 rows, mean-of-RMSE pessimism visible
  ev_per_draw <- .bace_eval_per_imputation(draws, truth_long, trait_types_lookup,
                                             t_fit_sec = 0,
                                             pool_method = "per_draw")
  pd_rows <- ev_per_draw[ev_per_draw$type == "continuous", , drop = FALSE]
  expect_equal(nrow(pd_rows), 3L)
  expect_true(all(pd_rows$rmse > direct_rmse))   # Jensen
})
