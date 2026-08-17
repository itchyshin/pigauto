skip_if_no_libtorch()

# zi_count conditional-on-nonzero conformal intervals (predict.pigauto_fit).
#
# compute_conformal_scores() (fit_helpers.R) has always scored zi_count on
# the magnitude (log1p-z) latent column, but predict.pigauto_fit() dropped
# zi_count from the type filter that builds conformal_lower / conformal_upper,
# so the score was computed but never surfaced. This restores
# conformal_lower / conformal_upper for zi_count as CONDITIONAL intervals for
# count | nonzero (back-transformed magnitude +/- score), NOT intervals for
# E[X] = P(nonzero) * E[count | nonzero]. The gate probability P(nonzero) is
# available separately in pred$probabilities.

make_zi_fixture <- function(n = 60L, seed = 21L) {
  set.seed(seed)
  tr <- ape::rtree(n)
  sp <- tr$tip.label
  # ~40% structural zeros, the rest small counts
  cnt <- ifelse(stats::runif(n) < 0.4, 0L, stats::rpois(n, 4) + 1L)
  df <- data.frame(row.names = sp, z = as.integer(cnt),
                   x = as.numeric(ape::rTraitCont(tr, model = "BM")))
  list(tree = tr, df = df)
}

test_that("zi_count fit produces finite conformal bounds where magnitude is defined", {
  fx <- make_zi_fixture()
  res <- impute(fx$df, fx$tree, trait_types = c(z = "zi_count"),
                epochs = 20L, verbose = FALSE, seed = 9L)
  pred <- res$prediction

  expect_true("z" %in% colnames(pred$conformal_lower))
  expect_true("z" %in% colnames(pred$conformal_upper))
  cs <- pred$conformal_scores["z"]
  skip_if(is.na(cs), "conformal score for z unavailable for this fixture/seed")

  expect_true(all(is.finite(pred$conformal_lower[, "z"])))
  expect_true(all(is.finite(pred$conformal_upper[, "z"])))
  expect_true(all(pred$conformal_lower[, "z"] <= pred$conformal_upper[, "z"]))
  # Conditional-on-nonzero counts are non-negative.
  expect_true(all(pred$conformal_lower[, "z"] >= 0))
})

test_that("zi_count conformal bounds back-transform monotonically with the score width", {
  fx <- make_zi_fixture()
  res <- impute(fx$df, fx$tree, trait_types = c(z = "zi_count"),
                epochs = 20L, verbose = FALSE, seed = 9L)
  fit <- res$fit
  skip_if(is.na(fit$conformal_scores["z"]),
          "conformal score for z unavailable for this fixture/seed")

  fit_narrow <- fit
  fit_narrow$conformal_scores["z"] <- 0.1
  fit_wide <- fit
  fit_wide$conformal_scores["z"] <- 2.0

  pred_narrow <- predict(fit_narrow, return_se = TRUE)
  pred_wide   <- predict(fit_wide, return_se = TRUE)

  width_narrow <- pred_narrow$conformal_upper[, "z"] - pred_narrow$conformal_lower[, "z"]
  width_wide   <- pred_wide$conformal_upper[, "z"] - pred_wide$conformal_lower[, "z"]

  # expm1() is monotonic increasing, so a wider latent-scale half-width must
  # not produce a narrower original-scale interval anywhere, and must
  # produce a strictly wider one somewhere.
  expect_true(all(width_wide >= width_narrow - 1e-8))
  expect_true(any(width_wide > width_narrow + 1e-8))

  # Point-estimate magnitude is unaffected by the score-width swap.
  expect_equal(pred_narrow$imputed$z, pred_wide$imputed$z)
})

test_that("zi_count conformal interval is conditional -- not an E[X] interval -- and the gate stays separately available", {
  fx <- make_zi_fixture()
  res <- impute(fx$df, fx$tree, trait_types = c(z = "zi_count"),
                epochs = 20L, verbose = FALSE, seed = 9L)
  pred <- res$prediction
  skip_if(is.na(pred$conformal_scores["z"]),
          "conformal score for z unavailable for this fixture/seed")

  # Gate probability P(nonzero) is available separately, distinct from the
  # conditional interval.
  expect_true("z" %in% names(pred$probabilities))
  p_nz <- pred$probabilities[["z"]]
  expect_true(all(p_nz >= 0 & p_nz <= 1))

  # The conditional interval brackets the conditional point estimate
  # (magnitude latent back-transformed, i.e. count | nonzero) -- it is
  # built from imputed_latent column 2 of the zi_count trait, not from the
  # decoded E[X] = P(nonzero) * E[count | nonzero] point estimate
  # `pred$imputed$z`.
  tm <- Filter(function(t) identical(t$name, "z"), pred$trait_map)[[1]]
  lc <- tm$latent_cols
  mag_latent <- pred$imputed_latent[, lc[2L]]
  expected_point <- pmax(expm1(mag_latent * tm$sd + tm$mean), 0)
  expect_true(all(pred$conformal_lower[, "z"] <= expected_point + 1e-6))
  expect_true(all(pred$conformal_upper[, "z"] >= expected_point - 1e-6))
  # And the conditional point estimate generally differs from decoded E[X]
  # (equal only where P(nonzero) happens to be 1 for every row).
  expect_true(any(abs(expected_point - pred$imputed$z) > 1e-6) || all(p_nz > 0.999))
})
