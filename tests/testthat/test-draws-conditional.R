# Joint conditional multivariate-BM draws (prototype, R/draws_conditional.R).
#
# The load-bearing test is the first one: the empirical mean AND FULL
# covariance of many joint draws from draw_conditional_bm() must match a
# dense Sigma-kron-R conditional computed directly, by Monte Carlo. If the
# draws were independent across tips (the multi_impute() conformal bug this
# prototype exists to fix), the empirical off-diagonal covariance among
# missing cells would be ~0 while the dense reference is not -- this test
# would catch that.

# Dense reference over TIP cells only (vec(L) ~ MVN(0, Sigma %x% R),
# column-major / trait-major cell order, matching test-exact-conditional.R's
# .dense_ref() convention). Returns the FULL conditional covariance at the
# missing cells, not just its diagonal.
.dcb_dense_ref <- function(L, Sigma, R) {
  V <- kronecker(Sigma, R)
  v <- as.numeric(L)
  oi <- which(!is.na(v)); mi <- which(is.na(v))
  Voo <- V[oi, oi, drop = FALSE]
  mu  <- as.numeric(V[mi, oi, drop = FALSE] %*% solve(Voo, v[oi]))
  cov <- V[mi, mi, drop = FALSE] -
    V[mi, oi, drop = FALSE] %*% solve(Voo, V[oi, mi, drop = FALSE])
  list(mu = mu, cov = cov, mi = mi)
}

.dcb_toy <- function(n = 60L, K = 2L, seed = 21L, miss = 0.3, rho = 0.7) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  Sigma_true <- matrix(c(1, rho, rho, 1), K, K)
  R <- stats::cov2cor(ape::vcv(tree))
  L <- t(chol(R)) %*% matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma_true)
  colnames(L) <- c("x", "y")
  rownames(L) <- tree$tip.label
  set.seed(seed + 1L)
  L[matrix(stats::runif(n * K) < miss, n, K)] <- NA

  traits <- as.data.frame(L)
  pd <- preprocess_traits(traits, tree, log_transform = FALSE)
  list(tree = tree, traits = traits, data = pd, L = L)
}

test_that("[draws] joint mean and FULL covariance match a dense Sigma-kron-R conditional (Monte Carlo)", {
  skip_if_not_installed("Matrix")
  d <- .dcb_toy()
  fr <- list(data = d$data, tree = d$tree)

  m <- 4000L
  res <- draw_conditional_bm(fr, m = m, seed = 99L)
  expect_s3_class(res, "pigauto_draws_conditional_bm")
  expect_length(res$datasets, m)
  expect_length(res$latent_draws, m)

  # Dense reference computed from the SAME Sigma the function estimated
  # internally, so this isolates correctness of the joint sampling step
  # (Sigma estimation itself is fit_mvn_bm_inhouse()'s job, tested in
  # test-joint-solver.R).
  R <- stats::cov2cor(ape::vcv(d$tree))
  ref <- .dcb_dense_ref(d$data$X_scaled, res$sigma_hat, R)

  # Stack draws: m x n_miss matrix, column order matches as.numeric(L)
  # (trait-major, species within trait), i.e. ref$mi.
  draw_mat <- t(vapply(res$latent_draws, function(Xk) as.numeric(Xk)[ref$mi],
                       numeric(length(ref$mi))))

  emp_mu  <- colMeans(draw_mat)
  emp_cov <- stats::cov(draw_mat)

  # Monte Carlo tolerance: SE(mean) ~ sqrt(diag(cov)/m); use a generous
  # multiple since this also has to tolerate the sparse Cholesky's own
  # floating-point error.
  se_mu <- sqrt(pmax(diag(ref$cov), 1e-8) / m)
  expect_true(all(abs(emp_mu - ref$mu) < pmax(8 * se_mu, 0.03)))

  # Covariance: compare the full matrix (this is the point of the test --
  # independent-across-tips draws would zero out the off-diagonal here).
  expect_equal(emp_cov, ref$cov, tolerance = 0.15)
})

test_that("[draws] observed cells are untouched in every draw", {
  skip_if_not_installed("Matrix")
  d <- .dcb_toy(seed = 5L)
  fr <- list(data = d$data, tree = d$tree)
  res <- draw_conditional_bm(fr, m = 5L, seed = 1L)

  obs <- !is.na(d$data$X_scaled)
  for (Xk in res$latent_draws) {
    expect_equal(Xk[obs], d$data$X_scaled[obs])
  }
  for (df in res$datasets) {
    expect_equal(as.matrix(df)[obs], as.matrix(d$traits)[obs])
  }
})

test_that("[draws] returns m completed datasets with all cells filled", {
  skip_if_not_installed("Matrix")
  d <- .dcb_toy(seed = 7L)
  fr <- list(data = d$data, tree = d$tree)
  m <- 6L
  res <- draw_conditional_bm(fr, m = m, seed = 2L)

  expect_equal(res$m, m)
  expect_length(res$datasets, m)
  for (df in res$datasets) {
    expect_false(anyNA(df))
    expect_equal(colnames(df), colnames(d$traits))
    expect_equal(rownames(df), rownames(d$traits))
  }
  # Different draws should differ from each other at missing cells.
  expect_false(isTRUE(all.equal(res$datasets[[1]], res$datasets[[2]])))
})

test_that("[draws] errors clearly on non-continuous trait types", {
  skip_if_not_installed("Matrix")
  set.seed(3L)
  tree <- ape::rcoal(20L)
  df <- data.frame(
    x = stats::rnorm(20L),
    b = factor(sample(c("a", "b"), 20L, TRUE)),
    row.names = tree$tip.label
  )
  df$x[1:4] <- NA
  pd <- preprocess_traits(df, tree)
  expect_error(draw_conditional_bm(list(data = pd, tree = tree), m = 3L),
              "continuous traits only")
})

test_that("[draws] errors clearly on multi-obs data", {
  skip_if_not_installed("Matrix")
  set.seed(9L)
  tree <- ape::rcoal(10L)
  df <- data.frame(
    species = rep(tree$tip.label, each = 2L),
    x = stats::rnorm(20L)
  )
  df$x[1:3] <- NA
  pd <- preprocess_traits(df, tree, species_col = "species")
  expect_error(draw_conditional_bm(list(data = pd, tree = tree), m = 3L),
              "single-observation-per-species")
})

test_that("[draws] errors clearly on a malformed fit_or_result", {
  skip_if_not_installed("Matrix")
  expect_error(draw_conditional_bm(list(foo = 1), m = 3L),
              "must expose \\$data")
})
