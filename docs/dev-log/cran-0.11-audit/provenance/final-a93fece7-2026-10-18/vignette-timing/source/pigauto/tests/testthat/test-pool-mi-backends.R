# End-to-end pool_mi() backend tests: multi_impute(draws_method = "posterior")
# -> with_imputations() -> pool_mi(), so the provenance path is exercised
# (the backend tests in test-multi-impute.R pass bare lists of fits).

pmb_mi <- function(n = 30L, m = 3L) {
  set.seed(21L)
  tree <- ape::rtree(n)
  R <- stats::cov2cor(ape::vcv(tree))
  Z <- t(chol(R + diag(1e-10, n))) %*% matrix(stats::rnorm(2L * n), n) %*%
    chol(matrix(c(1, 0.6, 0.6, 1), 2))
  df <- data.frame(mass = Z[, 1] * 2 + 10, wing = Z[, 2] * 2 + 10,
                   row.names = tree$tip.label)
  df$mass[c(2, 5, 9, 14)] <- NA
  df$wing[c(5, 7, 20, 25, 28)] <- NA
  df <- df[sample.int(n), , drop = FALSE]
  ctl <- list(n_chains = 2L, burnin = 60L, n_iter = 120L, keep_draws = 40L,
              auto_extend = FALSE)
  suppressWarnings(multi_impute(df, tree, m = m, draws_method = "posterior",
                                posterior_control = ctl, seed = 1,
                                verbose = FALSE))
}

pmb_check <- function(pooled, terms) {
  expect_s3_class(pooled, "data.frame")
  expect_true(all(terms %in% pooled$term))
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(is.finite(pooled$std.error)))
  expect_true(all(pooled$std.error > 0))
  # Between-imputation variance must enter: pooling only one fit, or fits on
  # identical datasets, would give riv = 0.
  expect_true(all(pooled$riv[pooled$term %in% terms] > 0))
}

test_that("pool_mi() pools glmmTMB fits from with_imputations()", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  mi <- pmb_mi()
  fits <- with_imputations(mi, function(d) glmmTMB::glmmTMB(mass ~ wing, data = d))
  pmb_check(pool_mi(fits), c("(Intercept)", "wing"))
})

test_that("pool_mi() pools lme4::lmer fits from with_imputations()", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  mi <- pmb_mi()
  fits <- with_imputations(mi, function(d) {
    d$grp <- factor(rep(seq_len(6L), length.out = nrow(d)))
    lme4::lmer(mass ~ wing + (1 | grp), data = d)
  })
  pmb_check(pool_mi(fits), c("(Intercept)", "wing"))
})

test_that("pool_mi() refuses brmsfit lists and points to brm_multiple()", {
  fake <- structure(list(), class = "brmsfit")
  expect_error(pool_mi(list(fake, fake)), "brm_multiple")
})
