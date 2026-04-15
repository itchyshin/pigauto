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
