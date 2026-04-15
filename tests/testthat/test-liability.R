test_that("liability_info() returns schema for each trait type", {
  # continuous: 1 liability, no threshold
  info <- liability_info(list(type = "continuous", n_latent = 1L))
  expect_equal(info$n_liability, 1L)
  expect_equal(info$kind, "continuous")

  # binary: 1 liability, single threshold at 0
  info <- liability_info(list(type = "binary", n_latent = 1L, levels = c("A","B")))
  expect_equal(info$n_liability, 1L)
  expect_equal(info$kind, "threshold")

  # categorical K=4: K liabilities, argmax observation
  info <- liability_info(list(type = "categorical", n_latent = 4L,
                              levels = c("a","b","c","d")))
  expect_equal(info$n_liability, 4L)
  expect_equal(info$kind, "argmax")

  # ordinal K=5: 1 liability, K-1 thresholds
  info <- liability_info(list(type = "ordinal", n_latent = 1L,
                              levels = as.character(1:5)))
  expect_equal(info$n_liability, 1L)
  expect_equal(info$kind, "ordered_threshold")
  expect_length(info$thresholds, 4L)  # K-1

  # count: 1 liability = log1p-scale
  info <- liability_info(list(type = "count", n_latent = 1L))
  expect_equal(info$n_liability, 1L)
  expect_equal(info$kind, "continuous")

  # proportion: 1 liability = logit-scale
  info <- liability_info(list(type = "proportion", n_latent = 1L))
  expect_equal(info$n_liability, 1L)
  expect_equal(info$kind, "continuous")

  # multi_proportion K=12: K liabilities with sum-zero constraint
  info <- liability_info(list(type = "multi_proportion", n_latent = 12L))
  expect_equal(info$n_liability, 12L)
  expect_equal(info$kind, "sum_zero")

  # zi_count: 2 liabilities (gate threshold + magnitude continuous)
  info <- liability_info(list(type = "zi_count", n_latent = 2L))
  expect_equal(info$n_liability, 2L)
  expect_equal(info$kind, "mixed_2")
})


test_that("liability_info() errors on unknown type", {
  # Catches accidental omission when a new trait type is added upstream
  # but the liability contract is not updated.
  expect_error(
    liability_info(list(type = "not_a_real_type", n_latent = 1L)),
    "No liability contract defined for type 'not_a_real_type'",
    fixed = TRUE
  )
})


test_that("estep_liability_binary returns posterior mean/var on truncated Gaussian", {
  # If y = 1, liability ~ truncated N(0, 1) from below at 0 → mean = sqrt(2/π) ≈ 0.7979
  res <- estep_liability_binary(y = 1, mu_prior = 0, sd_prior = 1)
  expect_equal(res$mean, sqrt(2 / pi), tolerance = 1e-6)
  expect_lt(res$var, 1)  # truncation always reduces variance

  # y = 0 → mean = -sqrt(2/π)
  res0 <- estep_liability_binary(y = 0, mu_prior = 0, sd_prior = 1)
  expect_equal(res0$mean, -sqrt(2 / pi), tolerance = 1e-6)

  # Very positive prior + y = 1: posterior mean close to prior mean (not much truncation)
  res2 <- estep_liability_binary(y = 1, mu_prior = 3, sd_prior = 1)
  expect_lt(abs(res2$mean - 3), 0.05)  # truncation barely bites
  expect_lt(abs(res2$var - 1), 0.02)   # analytic value ≈ 0.987; 0-pt is 3 SDs away
})


test_that("estep_liability_ordinal returns interval-truncated posterior", {
  thresholds <- c(-1, 0, 1)  # K = 4 classes, 3 thresholds

  # Class k=2 means liability in (-1, 0]. Standard N(0,1): should be moderately negative.
  res <- estep_liability_ordinal(k = 2L, thresholds = thresholds,
                                 mu_prior = 0, sd_prior = 1)
  expect_lt(res$mean, 0)
  expect_gt(res$mean, -1)
  expect_lt(res$var, 1)

  # Edge class k=1 (liability < -1): posterior mean should be less than -1
  res1 <- estep_liability_ordinal(k = 1L, thresholds = thresholds,
                                  mu_prior = 0, sd_prior = 1)
  expect_lt(res1$mean, -1)

  # Edge class k=4 (liability > 1): posterior mean > 1
  res4 <- estep_liability_ordinal(k = 4L, thresholds = thresholds,
                                  mu_prior = 0, sd_prior = 1)
  expect_gt(res4$mean, 1)
})

test_that("estep_liability_categorical plug-in: observed class gets boosted mean", {
  mu_prior <- c(0, 0, 0, 0)   # K=4 liabilities, prior means
  sd_prior <- c(1, 1, 1, 1)
  res <- estep_liability_categorical(k = 2L, mu_prior = mu_prior,
                                     sd_prior = sd_prior)
  expect_length(res$mean, 4L)
  expect_length(res$var, 4L)
  expect_gt(res$mean[2], res$mean[1])   # observed class is highest
  expect_gt(res$mean[2], res$mean[3])
  expect_gt(res$mean[2], res$mean[4])
  # Sum-to-zero projection: the predictions operate on a CLR-like space
  expect_equal(sum(res$mean), 0, tolerance = 1e-6)
})

test_that("estep_liability dispatches correctly on trait type", {
  # Binary
  tm_bin <- list(type = "binary", n_latent = 1L, levels = c("A","B"))
  # Observed "B" (2nd level) -> y = 1 in encoded form
  res <- estep_liability(tm_bin, observed = 1, mu_prior = 0, sd_prior = 1)
  expect_equal(res$mean, sqrt(2 / pi), tolerance = 1e-6)

  # Continuous: liability = observed value directly
  tm_cont <- list(type = "continuous", n_latent = 1L)
  res_c <- estep_liability(tm_cont, observed = 1.5, mu_prior = 0, sd_prior = 1)
  expect_equal(res_c$mean, 1.5)
  expect_equal(res_c$var, 0)  # observed directly -> zero posterior variance

  # NA observed (truly missing): returns prior
  res_na <- estep_liability(tm_cont, observed = NA_real_,
                            mu_prior = 0, sd_prior = 1)
  expect_equal(res_na$mean, 0)
  expect_equal(res_na$var, 1)
})


test_that("whole trait_map can be processed via estep_liability", {
  skip_on_cran()
  set.seed(1)
  tree <- ape::rcoal(30)

  # Mixed-type synthetic data
  df <- data.frame(
    mass    = stats::rnorm(30),
    clutch  = rpois(30, 5),
    diet    = factor(sample(c("insect", "plant", "fish"), 30, replace = TRUE)),
    threat  = ordered(sample(c("LC", "NT", "VU"), 30, replace = TRUE),
                      levels = c("LC","NT","VU"))
  )
  rownames(df) <- tree$tip.label

  pd <- preprocess_traits(df, tree, trait_types = c(clutch = "count"))

  # For each trait-map entry, call estep_liability with a prior and a
  # sample observed row.
  X <- pd$X_scaled
  for (tm in pd$trait_map) {
    lc <- tm$latent_cols
    observed <- X[1, lc]
    mu_prior <- rep(0, length(lc))
    sd_prior <- rep(1, length(lc))
    res <- estep_liability(tm, observed, mu_prior, sd_prior)
    expect_length(res$mean, length(lc))
    expect_length(res$var,  length(lc))
    expect_true(all(is.finite(res$mean)))
    expect_true(all(res$var >= 0))
  }
})
