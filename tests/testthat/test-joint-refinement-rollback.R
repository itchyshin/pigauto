test_that("divergence returns the last accepted joint refinement", {
  set.seed(59)
  tree <- ape::rtree(16)
  i <- seq_len(16)
  traits <- data.frame(x = sin(i * .7), y = cos(i * .4), row.names = tree$tip.label)
  traits$x[c(2, 4, 6, 8, 10)] <- NA_real_
  traits$y[c(3, 5, 7, 9, 11)] <- NA_real_
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous", y = "continuous"),
                            log_transform = FALSE)
  L <- data$X_scaled
  fit <- function(iter) fit_mvn_bm_inhouse(L, tree = tree, predict_method = "per_column",
                                            max_iter = iter, tol = 1e-14)
  initial <- fit(0L)
  accepted <- fit(1L)
  stopped <- fit(8L)
  expect_false(accepted$diverged)
  expect_identical(accepted$n_iter, 1L)
  expect_true(stopped$diverged)
  expect_identical(stopped$n_iter, 1L)
  expect_gt(max(abs(accepted$anc_recon - initial$anc_recon)), .1)
  expect_gt(norm(accepted$pars$phylocov - initial$pars$phylocov, "F"), .01)

  # Independently calculate both steps for this two-trait case: the normal
  # conditional for the other trait, inverse-variance mean update, and dense
  # correlation-matrix covariance M-step. Conservative variance is retained.
  R <- stats::cov2cor(ape::vcv.phylo(tree))[rownames(L), rownames(L)]
  step <- function(mean, variance, Sigma) {
    next_mean <- mean
    for (j in seq_len(ncol(L))) {
      other <- setdiff(seq_len(ncol(L)), j)
      missing <- which(is.na(L[, j]))
      inverse <- solve(Sigma[other, other, drop = FALSE] + diag(1e-8, length(other)))
      cross_mean <- as.vector(mean[missing, other, drop = FALSE] %*%
                              inverse %*% Sigma[other, j, drop = FALSE])
      cross_var <- max(Sigma[j, j] - as.numeric(Sigma[j, other, drop = FALSE] %*%
                                                inverse %*% Sigma[other, j, drop = FALSE]), 1e-8)
      bm_var <- pmax(variance[missing, j], 1e-8)
      next_mean[missing, j] <- (mean[missing, j] / bm_var + cross_mean / cross_var) /
                                (1 / bm_var + 1 / cross_var)
    }
    next_sigma <- (crossprod(next_mean, solve(R, next_mean)) +
                   diag(colSums(variance), ncol(L))) / nrow(L)
    list(mean = next_mean, variance = variance, Sigma = next_sigma)
  }
  first <- step(initial$anc_recon, initial$anc_var, initial$pars$phylocov)
  second <- step(first$mean, first$variance, first$Sigma)
  delta1 <- norm(first$Sigma - initial$pars$phylocov, "F") /
    max(norm(initial$pars$phylocov, "F"), 1e-8)
  delta2 <- norm(second$Sigma - first$Sigma, "F") / max(norm(first$Sigma, "F"), 1e-8)
  expect_true(is.finite(delta1) && is.finite(delta2))
  expect_gt(delta2, delta1)
  expect_equal(accepted$anc_recon, first$mean, tolerance = 1e-7)
  expect_equal(unname(accepted$pars$phylocov), unname(first$Sigma), tolerance = 1e-7)
  expect_equal(stopped$anc_recon, accepted$anc_recon, tolerance = 1e-8)
  expect_equal(stopped$anc_var, accepted$anc_var, tolerance = 1e-8)
  expect_equal(unname(stopped$pars$phylocov), unname(accepted$pars$phylocov), tolerance = 1e-8)
  expect_equal(unname(stopped$anc_recon[!is.na(L)]), unname(L[!is.na(L)]))

  baseline <- function(iter) fit_baseline(data, tree, predict_method = "per_column",
                                        lambda_mode = "fixed_1", joint_refine_iter = iter)
  zero <- baseline(0L)
  one <- baseline(1L)
  many <- baseline(8L)
  expect_identical(unname(many$path), c("joint_mvn", "joint_mvn"))
  expect_gt(max(abs(one$mu - zero$mu)), .1)
  expect_equal(many$mu, one$mu, tolerance = 1e-8)
  expect_equal(many$se, one$se, tolerance = 1e-8)
  expect_equal(unname(many$mu[!is.na(L)]), unname(L[!is.na(L)]))
})

test_that("divergence retains the accepted pooled variance state", {
  set.seed(59)
  tree <- ape::rtree(16)
  i <- seq_len(16)
  traits <- data.frame(x = sin(i * .7), y = cos(i * .4), row.names = tree$tip.label)
  traits$x[c(2, 4, 6, 8, 10)] <- NA_real_
  traits$y[c(3, 5, 7, 9, 11)] <- NA_real_
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous", y = "continuous"),
                            log_transform = FALSE)
  fit <- function(iter) fit_mvn_bm_inhouse(data$X_scaled, tree = tree,
    predict_method = "per_column", max_iter = iter, tol = 1e-14,
    refine_variance = "pooled")
  initial <- fit(0L)
  accepted <- fit(1L)
  stopped <- fit(8L)
  expect_true(stopped$diverged)
  expect_identical(stopped$n_iter, 1L)
  expect_gt(max(abs(accepted$anc_var - initial$anc_var)), .1)
  expect_equal(stopped$anc_recon, accepted$anc_recon, tolerance = 1e-8)
  expect_equal(stopped$anc_var, accepted$anc_var, tolerance = 1e-8)
  expect_equal(unname(stopped$pars$phylocov), unname(accepted$pars$phylocov), tolerance = 1e-8)
  expect_equal(unname(stopped$anc_recon[!is.na(data$X_scaled)]),
               unname(data$X_scaled[!is.na(data$X_scaled)]))
})
