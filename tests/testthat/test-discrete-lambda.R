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

# ---- "auto": per-trait choice by validation Brier ---------------------------

test_that("select_discrete_lambda: lower Brier wins, ties and NA keep fixed_1", {
  bf <- c(a = 0.20, b = 0.20, c = 0.30, d = NA, e = 0.25)
  be <- c(a = 0.15, b = 0.20, c = 0.35, d = 0.10, e = NA)
  out <- select_discrete_lambda(bf, be)
  expect_identical(out, c(a = "estimate", b = "fixed_1", c = "fixed_1",
                          d = "fixed_1", e = "fixed_1"))
})

test_that("discrete_val_brier matches a hand calculation", {
  p <- matrix(c(0.8, 0.2, 0.5, 0.5, 0.1, 0.9), ncol = 2, byrow = TRUE)
  y <- matrix(c(1, 0, 0, 1, 0, 1), ncol = 2, byrow = TRUE)
  # rows 1 and 2: (0.04+0.04) and (0.25+0.25) -> mean 0.29
  expect_equal(discrete_val_brier(p, y, c(1L, 2L)), 0.29)
  expect_true(is.na(discrete_val_brier(p, y, integer(0))))
})

test_that("auto records one choice per binary/categorical trait", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda(n = 150, lambda = 0.3)
  pd <- preprocess_traits(s$d, s$tree)
  set.seed(1)
  sp <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, val_frac = 0.5)
  old <- options(pigauto.discrete_lambda = "auto"); on.exit(options(old), add = TRUE)
  b <- suppressWarnings(fit_baseline(pd, s$tree, splits = sp))
  ch <- b$discrete_lambda_chosen
  expect_false(is.null(ch))
  expect_setequal(names(ch), c("bin", "cat3"))
  expect_length(ch, 2L)
  expect_true(all(ch %in% c("fixed_1", "estimate")))
})

test_that("auto is a no-op relative to fixed_1 when it picks fixed_1 everywhere, and unset has no field", {
  skip_if_not(joint_mvn_available())
  s <- sim_discrete_lambda(n = 150, lambda = 1)
  pd <- preprocess_traits(s$d, s$tree)
  set.seed(2)
  sp <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, val_frac = 0.5)
  old <- options(pigauto.discrete_lambda = NULL); on.exit(options(old), add = TRUE)
  b0 <- suppressWarnings(fit_baseline(pd, s$tree, splits = sp))
  expect_null(b0$discrete_lambda_chosen)
})

test_that("impute() under auto fills every discrete cell at strong and weak signal", {
  skip_if_not(joint_mvn_available())
  old <- options(pigauto.discrete_lambda = "auto"); on.exit(options(old), add = TRUE)
  for (lam in c(1, 0.1)) {
    s <- sim_discrete_lambda(n = 150, lambda = lam, seed = 11)
    res <- suppressWarnings(impute(s$d, s$tree, verbose = FALSE))
    expect_false(anyNA(res$completed$bin))
    expect_false(anyNA(res$completed$cat3))
  }
})
