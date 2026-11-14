# C-Blk1 (P0): conformal MI for zi_count must score / draw on magnitude
# latent scale, not pred$se of E[X]. Gate stays Bernoulli.

test_that("compute_conformal_scores scores zi_count on magnitude latent, not gate", {
  n <- 8L
  trait_map <- list(list(
    name = "zi",
    type = "zi_count",
    latent_cols = 1:2,
    mean = 0,
    sd = 1
  ))
  # Col 1 = gate, col 2 = log1p-z magnitude. Pure-BM blend: pred = mu_cal.
  mu_cal <- matrix(c(rep(0, n), rep(1, n)), n, 2)
  delta_cal <- matrix(c(rep(10, n), rep(9, n)), n, 2)
  X_truth_r <- matrix(c(rep(1, n), rep(1, n)), n, 2)
  val_mask <- matrix(TRUE, n, 2)

  scores <- pigauto:::compute_conformal_scores(
    trait_map = trait_map,
    calibrated_gates = c(0, 0),
    mu_cal = mu_cal,
    delta_cal = delta_cal,
    X_truth_r = X_truth_r,
    val_mask_mat = val_mask,
    r_cal_bm = c(1, 1),
    r_cal_gnn = c(0, 0),
    r_cal_mean = c(0, 0),
    mean_baseline_per_col = c(0, 0)
  )

  # Magnitude residual is |1 - 1| = 0. Scoring the gate would give |1 - 0| = 1.
  expect_equal(unname(scores["zi"]), 0, tolerance = 1e-12)
})

test_that("compute_conformal_scores zi_count quantile matches mag residuals", {
  n <- 10L
  trait_map <- list(list(
    name = "zi",
    type = "zi_count",
    latent_cols = 1:2,
    mean = 0,
    sd = 1
  ))
  mu_cal <- matrix(c(rep(0, n), rep(2, n)), n, 2)
  delta_cal <- matrix(0, n, 2)
  X_truth_r <- matrix(c(rep(0, n), rep(1, n)), n, 2)
  val_mask <- matrix(TRUE, n, 2)

  scores <- pigauto:::compute_conformal_scores(
    trait_map = trait_map,
    calibrated_gates = c(0, 0),
    mu_cal = mu_cal,
    delta_cal = delta_cal,
    X_truth_r = X_truth_r,
    val_mask_mat = val_mask
  )

  # All mag residuals = 1; split quantile at 95% is 1.
  expect_equal(unname(scores["zi"]), 1, tolerance = 1e-12)
})

test_that(".sample_conformal_draw zi_count does not use E[X] SE as count|nz SD", {
  n <- 80L
  tm <- list(
    name = "parasites",
    type = "zi_count",
    latent_cols = 1:2,
    mean = 0,
    sd = 1
  )
  count_hat <- 10
  z_mag <- log1p(count_hat)
  pred <- list(
    imputed = data.frame(parasites = as.integer(rep(count_hat, n))),
    imputed_latent = cbind(rep(8, n), rep(z_mag, n)),
    se = matrix(100, n, 1, dimnames = list(NULL, "parasites")),
    se_latent = cbind(rep(NA_real_, n), rep(0.05, n)),
    probabilities = list(parasites = rep(1, n)),
    conformal_scores = c(parasites = 0.2)
  )
  imputed_mask <- matrix(TRUE, n, 1, dimnames = list(NULL, "parasites"))

  draw <- pigauto:::.sample_conformal_draw(
    pred,
    imputed_mask,
    list(tm),
    seed_i = 42L
  )
  s <- stats::sd(as.numeric(draw$parasites))

  # Wrong path: N(cond_mu, pred$se=100) → SD tens–hundreds after round.
  # Right path: mag z ~ N(log1p(10), 0.2/1.96) → count SD ~ 1.
  expect_lt(s, 10)
  expect_gt(s, 0.2)
  expect_true(all(draw$parasites >= 0L))
  expect_true(all(draw$parasites > 0L))
})

test_that(".sample_conformal_draw zi_count fallback uses mag latent SE not E[X] SE", {
  n <- 80L
  tm <- list(
    name = "parasites",
    type = "zi_count",
    latent_cols = 1:2,
    mean = 0,
    sd = 1
  )
  count_hat <- 10
  z_mag <- log1p(count_hat)
  pred <- list(
    imputed = data.frame(parasites = as.integer(rep(count_hat, n))),
    imputed_latent = cbind(rep(8, n), rep(z_mag, n)),
    se = matrix(100, n, 1, dimnames = list(NULL, "parasites")),
    se_latent = cbind(rep(NA_real_, n), rep(0.08, n)),
    probabilities = list(parasites = rep(1, n)),
    conformal_scores = c(parasites = NA_real_)
  )
  imputed_mask <- matrix(TRUE, n, 1, dimnames = list(NULL, "parasites"))

  draw <- pigauto:::.sample_conformal_draw(
    pred,
    imputed_mask,
    list(tm),
    seed_i = 7L
  )
  s <- stats::sd(as.numeric(draw$parasites))
  expect_lt(s, 10)
  expect_gt(s, 0.1)
})

test_that(".sample_conformal_draw zi_count gate stays Bernoulli, not Normal(E[X])", {
  n <- 60L
  tm <- list(
    name = "parasites",
    type = "zi_count",
    latent_cols = 1:2,
    mean = 0,
    sd = 1
  )
  pred <- list(
    imputed = data.frame(parasites = as.integer(rep(8L, n))),
    imputed_latent = cbind(rep(0, n), rep(log1p(8), n)),
    se = matrix(50, n, 1, dimnames = list(NULL, "parasites")),
    se_latent = cbind(rep(NA_real_, n), rep(0.05, n)),
    probabilities = list(parasites = c(rep(0, n / 2), rep(1, n / 2))),
    conformal_scores = c(parasites = 0.15)
  )
  imputed_mask <- matrix(TRUE, n, 1, dimnames = list(NULL, "parasites"))

  draw <- pigauto:::.sample_conformal_draw(
    pred,
    imputed_mask,
    list(tm),
    seed_i = 99L
  )
  expect_identical(draw$parasites[seq_len(n / 2)], rep(0L, n / 2))
  expect_true(all(draw$parasites[(n / 2 + 1L):n] > 0L))
})

test_that("multi_impute conformal with trait_types zi_count stores mag-scale scores", {
  skip_if_not_installed("torch")
  set.seed(20260808L)
  n <- 48L
  tree <- ape::rtree(n)
  counts <- stats::rpois(n, lambda = 8)
  counts[sample.int(n, round(n * 0.2))] <- 0L
  df <- data.frame(parasites = as.integer(counts), row.names = tree$tip.label)
  miss <- c(2L, 5L, 9L, 14L, 20L, 28L, 35L, 42L)
  df$parasites[miss] <- NA_integer_

  mi <- suppressWarnings(multi_impute(
    df,
    tree,
    m = 5L,
    draws_method = "conformal",
    trait_types = c(parasites = "zi_count"),
    epochs = 12L,
    eval_every = 6L,
    patience = 4L,
    verbose = FALSE,
    seed = 20260808L
  ))

  expect_equal(mi$draws_method, "conformal")
  expect_true("parasites" %in% names(mi$fit$conformal_scores))
  expect_true(is.finite(mi$fit$conformal_scores[["parasites"]]))

  draw_mat <- vapply(
    mi$datasets,
    function(d) as.integer(d$parasites[miss]),
    integer(length(miss))
  )
  per_cell_sd <- apply(draw_mat, 1L, stats::sd)
  expect_true(
    any(per_cell_sd > 0),
    info = "zi_count conformal MI should vary missing cells across draws"
  )
  expect_true(all(unlist(lapply(mi$datasets, function(d) d$parasites)) >= 0L))
})

test_that("multi_impute mc_dropout with trait_types zi_count still returns m datasets", {
  skip_if_not_installed("torch")
  set.seed(20260809L)
  n <- 28L
  tree <- ape::rtree(n)
  counts <- stats::rpois(n, lambda = 6)
  counts[sample.int(n, round(n * 0.3))] <- 0L
  df <- data.frame(parasites = as.integer(counts), row.names = tree$tip.label)
  df$parasites[c(3L, 8L, 15L, 22L)] <- NA_integer_

  mi <- suppressWarnings(multi_impute(
    df,
    tree,
    m = 3L,
    draws_method = "mc_dropout",
    trait_types = c(parasites = "zi_count"),
    epochs = 10L,
    eval_every = 5L,
    patience = 4L,
    verbose = FALSE,
    seed = 20260809L
  ))

  expect_equal(mi$draws_method, "mc_dropout")
  expect_equal(length(mi$datasets), 3L)
  for (d in mi$datasets) {
    expect_equal(nrow(d), n)
    expect_false(any(is.na(d$parasites)))
    expect_true(all(d$parasites >= 0L))
  }
})
