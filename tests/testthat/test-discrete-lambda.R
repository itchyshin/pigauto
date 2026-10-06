# Prototype (feat/discrete-lambda): options(pigauto.discrete_lambda = "estimate") estimates Pagel's
# lambda on binary / ordinal liability columns in fit_joint_threshold_baseline(), which also reaches
# categorical traits through the one-vs-rest fits. Off (the default) must leave everything unchanged.

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

test_that("discrete lambda is off by default and leaves the threshold baseline unchanged", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda()
  pd <- preprocess_traits(s$d, s$tree)
  old <- options(pigauto.discrete_lambda = NULL); on.exit(options(old), add = TRUE)
  a <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL)
  options(pigauto.discrete_lambda = "fixed_1")
  b <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL)
  expect_identical(a$mu_liab, b$mu_liab)
  expect_identical(a$se_liab, b$se_liab)
  bin_col <- which(a$liab_types == "binary")
  expect_true(all(a$lambda_per_trait_fit[bin_col] == 1))
})

test_that("discrete lambda = 'estimate' estimates lambda on the binary liability column", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda(lambda = 0.1)
  pd <- preprocess_traits(s$d, s$tree)
  old <- options(pigauto.discrete_lambda = "estimate"); on.exit(options(old), add = TRUE)
  fit <- fit_joint_threshold_baseline(pd, s$tree, splits = NULL)
  bin_col <- which(fit$liab_types[fit$fit_cols_idx] == "binary")
  lam <- fit$lambda_per_trait_fit[bin_col]
  expect_true(all(is.finite(lam)))
  expect_true(all(lam >= 0 & lam <= 1))
  expect_true(any(lam < 1))
  expect_true(all(is.finite(fit$mu_liab[, fit$fit_cols_idx])))
})

test_that("impute() runs end to end with discrete lambda on and fills every missing cell", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda()
  old <- options(pigauto.discrete_lambda = "estimate"); on.exit(options(old), add = TRUE)
  res <- suppressWarnings(impute(s$d, s$tree, verbose = FALSE))
  expect_false(anyNA(res$completed$bin))
  expect_false(anyNA(res$completed$cat3))
})
