# discrete_lambda = c("estimate", "fixed_1") (default "estimate") estimates Pagel's lambda on binary /
# ordinal liability columns in fit_joint_threshold_baseline(), which also reaches categorical traits
# through the one-vs-rest fits. "fixed_1" must reproduce the previous behaviour exactly.

sim_discrete_lambda <- function(n = 120, lambda = 0.3, seed = 7) {
  set.seed(seed)
  tr <- ape::rcoal(n)
  C <- ape::vcv(tr, corr = TRUE)
  V <- lambda * C + (1 - lambda) * diag(n)
  L <- t(chol(V)) %*% matrix(stats::rnorm(n * 3), n)
  rownames(L) <- tr$tip.label
  d <- data.frame(c1 = L[, 1],
                  bin = factor(ifelse(L[, 2] > 0, "yes", "no")),
                  cat3 = factor(cut(L[, 3], stats::qnorm(c(0, 1/3, 2/3, 1)), labels = c("a", "b", "c"))),
                  row.names = tr$tip.label)
  for (v in names(d)) d[sample(n, round(0.25 * n)), v] <- NA
  list(tree = tr, d = d)
}

test_that("discrete_lambda = 'fixed_1' leaves the threshold baseline at lambda = 1 (byte-identical to the old path)", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda()
  pd <- preprocess_traits(s$d, s$tree)
  a <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL, discrete_lambda = "fixed_1")
  b <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL, discrete_lambda = "fixed_1")
  expect_identical(a$mu_liab, b$mu_liab)
  bin_col <- which(a$liab_types[a$fit_cols_idx] == "binary")
  expect_true(all(a$lambda_per_trait_fit[bin_col] == 1))
  # the default is "estimate", so it must differ from fixed_1 at low signal
  d <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL)
  expect_false(identical(a$mu_liab, d$mu_liab))
})

test_that("the default estimates lambda on the binary liability column", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda(lambda = 0.1)
  pd <- preprocess_traits(s$d, s$tree)
  fit <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL)
  bin_col <- which(fit$liab_types[fit$fit_cols_idx] == "binary")
  lam <- fit$lambda_per_trait_fit[bin_col]
  expect_true(all(is.finite(lam)))
  expect_true(all(lam >= 0 & lam <= 1))
  expect_true(any(lam < 1))
  expect_true(all(is.finite(fit$mu_liab[, fit$fit_cols_idx])))
})

test_that("impute() fills every missing cell and records discrete_lambda in model_config", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda()
  res <- suppressWarnings(impute(s$d, s$tree, verbose = FALSE))
  expect_false(anyNA(res$completed$bin))
  expect_false(anyNA(res$completed$cat3))
  expect_identical(res$fit$model_config$discrete_lambda, "estimate")
  res1 <- suppressWarnings(impute(s$d, s$tree, verbose = FALSE, discrete_lambda = "fixed_1"))
  expect_identical(res1$fit$model_config$discrete_lambda, "fixed_1")
  expect_error(impute(s$d, s$tree, verbose = FALSE, discrete_lambda = "bogus"))
})

test_that("a lambda_fixed rebuild from the fitted lambda_per_trait reproduces the baseline", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda(lambda = 0.2)
  pd <- preprocess_traits(s$d, s$tree)
  bl <- fit_baseline(pd, s$tree, splits = NULL, predict_method = "exact")
  expect_identical(bl$discrete_lambda, "estimate")
  bin_nm <- grep("^bin", names(bl$lambda_per_trait), value = TRUE)
  expect_true(length(bin_nm) >= 1L)
  expect_true(any(bl$lambda_per_trait[bin_nm] < 1))
  rb <- fit_baseline(pd, s$tree, splits = NULL, predict_method = "exact",
                     lambda_fixed = bl$lambda_per_trait)
  expect_equal(rb$mu, bl$mu, tolerance = 1e-6)
  expect_equal(rb$lambda_per_trait, bl$lambda_per_trait, tolerance = 1e-6)
})
