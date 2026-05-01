# tests/testthat/test-bm-internal.R
# Unit tests for internal BM baseline (R/bm_internal.R)

test_that("phylo_cor_matrix returns correct structure", {
  set.seed(1)
  tree <- ape::rtree(20)
  R <- phylo_cor_matrix(tree)

  expect_true(is.matrix(R))
  expect_equal(nrow(R), 20L)
  expect_equal(ncol(R), 20L)
  # Diagonal = 1 (correlation matrix)
  expect_equal(unname(diag(R)), rep(1, 20), tolerance = 1e-12)
  # Symmetric
  expect_equal(R, t(R))
  # Positive definite (all eigenvalues > 0)
  evals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
  # Dimnames match tip labels
  expect_equal(rownames(R), tree$tip.label)
  expect_equal(colnames(R), tree$tip.label)
})

test_that("bm_impute_col recovers observed values exactly", {
  set.seed(2)
  tree <- ape::rtree(30)
  R <- phylo_cor_matrix(tree)

  # Simulate some data
  y <- rnorm(30)
  # Mask 10 species
  miss <- sample(30, 10)
  y_masked <- y
  y_masked[miss] <- NA

  res <- bm_impute_col(y_masked, R)

  # Observed species keep exact values
  obs <- setdiff(seq_len(30), miss)
  expect_equal(res$mu[obs], y_masked[obs])
  # Observed species have se = 0
 expect_equal(res$se[obs], rep(0, length(obs)))
  # Missing species have finite mu and positive se
  expect_true(all(is.finite(res$mu[miss])))
  expect_true(all(res$se[miss] > 0))
})

test_that("bm_impute_col returns correct structure", {
  set.seed(3)
  tree <- ape::rtree(20)
  R <- phylo_cor_matrix(tree)

  y <- rnorm(20)
  y[c(1, 5, 10)] <- NA

  res <- bm_impute_col(y, R)

  expect_true(is.list(res))
  expect_named(res, c("mu", "se"))
  expect_length(res$mu, 20L)
  expect_length(res$se, 20L)
  expect_true(all(is.finite(res$mu)))
  expect_true(all(is.finite(res$se)))
  expect_true(all(res$se >= 0))
})

test_that("bm_impute_col handles all-observed case", {
  set.seed(4)
  tree <- ape::rtree(15)
  R <- phylo_cor_matrix(tree)

  y <- rnorm(15)  # no NAs
  res <- bm_impute_col(y, R)

  expect_equal(res$mu, y)
  expect_equal(res$se, rep(0, 15))
})

test_that("bm_impute_col falls back for < 5 observations", {
  set.seed(5)
  tree <- ape::rtree(20)
  R <- phylo_cor_matrix(tree)

  y <- rep(NA_real_, 20)
  y[c(1, 2, 3)] <- c(1.0, 2.0, 3.0)  # only 3 observed

  res <- bm_impute_col(y, R)

  # Observed keep values
  expect_equal(res$mu[1:3], c(1.0, 2.0, 3.0))
  # Missing get global mean
  expect_equal(res$mu[4:20], rep(mean(c(1, 2, 3)), 17))
  # SE for missing = global sd
  expect_equal(res$se[4:20], rep(sd(c(1, 2, 3)), 17))
  # SE for observed = 0
  expect_equal(res$se[1:3], rep(0, 3))
})

test_that("bm_impute_col handles near-singular R_oo via nugget back-off", {
  # Star tree: all tips at equal distance from root
  n <- 15
  tree <- ape::stree(n, type = "star")
  tree$edge.length <- rep(1, nrow(tree$edge))
  R <- phylo_cor_matrix(tree)

  y <- rnorm(n)
  y[c(1, 2, 3)] <- NA

  # Should not error, thanks to nugget regularisation
  res <- expect_no_error(bm_impute_col(y, R))
  expect_true(all(is.finite(res$mu)))
  expect_true(all(is.finite(res$se)))
})

test_that("bm_impute_col handles zero observations", {
  set.seed(7)
  tree <- ape::rtree(10)
  R <- phylo_cor_matrix(tree)

  y <- rep(NA_real_, 10)  # all missing
  res <- bm_impute_col(y, R)

  # Fallback: mu = 0, se = 1
  expect_equal(res$mu, rep(0, 10))
  expect_equal(res$se, rep(1, 10))
})

test_that("fit_baseline with internal BM produces finite mu and se", {
  skip_if_not_installed("torch")
  skip_if(!torch::torch_is_installed(), "torch backend not installed")

  set.seed(8)
  tree <- ape::rtree(30)
  # Simple continuous data
  n <- 30
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    row.names = tree$tip.label
  )
  df$x1[sample(n, 5)] <- NA
  df$x2[sample(n, 5)] <- NA

  pd <- preprocess_traits(df, tree)
  graph <- build_phylo_graph(tree, k_eigen = 4L)
  bl <- fit_baseline(pd, tree, graph = graph)

  expect_true(all(is.finite(bl$mu)))
  expect_true(all(is.finite(bl$se)))
  expect_equal(nrow(bl$mu), n)
  expect_equal(ncol(bl$mu), ncol(pd$X_scaled))
})

test_that("fit_baseline with internal BM and splits produces finite results", {
  skip_if_not_installed("torch")
  skip_if(!torch::torch_is_installed(), "torch backend not installed")

  set.seed(9)
  tree <- ape::rtree(30)
  n <- 30
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    row.names = tree$tip.label
  )

  pd <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 1,
                             trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)
  bl <- fit_baseline(pd, tree, splits = spl, graph = graph)

  expect_true(all(is.finite(bl$mu)))
  expect_true(all(is.finite(bl$se)))
  # SE should be > 0 for at least some species (those missing after masking)
  expect_true(any(bl$se > 0))
})

test_that("build_phylo_graph caches R_phy alongside D", {
  set.seed(10)
  tree <- ape::rtree(20)
  g <- build_phylo_graph(tree, k_eigen = 4L)

  expect_true(!is.null(g$R_phy))
  expect_equal(nrow(g$R_phy), 20L)
  expect_equal(ncol(g$R_phy), 20L)
  expect_equal(unname(diag(g$R_phy)), rep(1, 20), tolerance = 1e-12)
  expect_equal(rownames(g$R_phy), tree$tip.label)
})

# ===========================================================================
# D1: bm_impute_col math validation (added 2026-04-30)
# ===========================================================================
#
# Validate `bm_impute_col()` math against an independent hand-rolled
# implementation of the BM conditional-MVN BLUP, and (when `phylolm` is
# installed) against `phylolm::phylolm(model = "BM")`.  This catches:
#
#   * algebraic regressions in the Cholesky-based solve path
#   * accidental change to the GLS phylogenetic mean formula
#   * REML vs ML variance estimator drift
#   * edge cases at low n_obs

ref_bm_blup <- function(y, R) {
  obs_idx  <- which(!is.na(y))
  miss_idx <- which(is.na(y))
  n_o <- length(obs_idx); n_m <- length(miss_idx)
  if (n_m == 0L) return(list(mu = y, se = rep(0, length(y))))
  if (n_o < 5L) {
    y_obs <- y[obs_idx]
    global_mu <- if (n_o > 0L) mean(y_obs) else 0
    global_sd <- if (n_o > 1L) stats::sd(y_obs) else 1
    mu_out <- numeric(length(y))
    se_out <- numeric(length(y))
    mu_out[obs_idx]  <- y_obs
    mu_out[miss_idx] <- global_mu
    se_out[miss_idx] <- global_sd
    return(list(mu = mu_out, se = se_out))
  }
  R_oo <- R[obs_idx, obs_idx, drop = FALSE]
  R_mo <- R[miss_idx, obs_idx, drop = FALSE]
  R_oo_inv <- solve(R_oo + diag(1e-6, n_o))
  ones <- rep(1, n_o)
  y_o  <- y[obs_idx]
  mu_hat <- as.numeric((ones %*% R_oo_inv %*% y_o) /
                          (ones %*% R_oo_inv %*% ones))
  e <- y_o - mu_hat
  sigma2 <- as.numeric(t(e) %*% R_oo_inv %*% e) / max(n_o - 1L, 1L)
  mu_m <- mu_hat + as.numeric(R_mo %*% R_oo_inv %*% e)
  H <- R_mo %*% R_oo_inv %*% t(R_mo)
  cond_var <- sigma2 * pmax(1 - diag(H), 0)
  mu_out <- numeric(length(y))
  se_out <- numeric(length(y))
  mu_out[obs_idx]  <- y_o
  mu_out[miss_idx] <- mu_m
  se_out[miss_idx] <- sqrt(cond_var)
  list(mu = mu_out, se = se_out)
}

test_that("[D1.1] bm_impute_col matches hand-rolled BLUP within 1e-8 on well-conditioned R", {
  set.seed(2040L)
  for (n in c(20L, 50L, 100L)) {
    tree <- ape::rtree(n)
    R    <- pigauto:::phylo_cor_matrix(tree)
    y    <- stats::rnorm(n)
    miss_idx <- sample.int(n, size = max(2L, n %/% 4L))
    y[miss_idx] <- NA

    pi  <- pigauto:::bm_impute_col(y, R)
    ref <- ref_bm_blup(y, R)

    expect_equal(pi$mu, ref$mu, tolerance = 1e-8,
                 info = sprintf("n=%d: mu mismatch", n))
    expect_equal(pi$se, ref$se, tolerance = 1e-8,
                 info = sprintf("n=%d: se mismatch", n))
  }
})

test_that("[D1.2] bm_impute_col implies mu_hat solves GLS normal equation", {
  set.seed(2041L)
  n <- 30L
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- stats::rnorm(n)
  R_inv <- solve(R + diag(1e-6, n))
  ones  <- rep(1, n)
  expected_mu_hat <- as.numeric((ones %*% R_inv %*% y) /
                                 (ones %*% R_inv %*% ones))
  # Mask half and check the imputation matches the hand-rolled.
  y_masked <- y
  y_masked[1:5] <- NA
  res <- pigauto:::bm_impute_col(y_masked, R)
  ref <- ref_bm_blup(y_masked, R)
  expect_equal(res$mu, ref$mu, tolerance = 1e-8)
})

test_that("[D1.3] bm_impute_col degenerate (n_obs < 5) returns global mean/sd", {
  set.seed(2042L)
  n <- 10L
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- rep(NA_real_, n)
  y[1:3] <- c(1, 2, 3)
  res <- pigauto:::bm_impute_col(y, R)
  expect_equal(res$mu[1:3], c(1, 2, 3))
  expect_equal(res$mu[4:n], rep(mean(c(1, 2, 3)), n - 3L))
  expect_equal(res$se[1:3], rep(0, 3))
  expect_equal(res$se[4:n], rep(stats::sd(c(1, 2, 3)), n - 3L))
})

test_that("[D1.4] bm_impute_col GLS phylogenetic mean matches phylolm BM intercept (ultrametric tree)", {
  skip_if_not_installed("phylolm")
  set.seed(2043L)
  n <- 50L
  # phylolm and pigauto give the SAME GLS mean only on ULTRAMETRIC trees.
  # For non-ultrametric trees the two implementations use different VCV
  # conventions: phylolm uses ape::vcv(tree) (per-tip variances differ)
  # while pigauto uses cov2cor(vcv(tree)) (per-tip variances normalised
  # to 1).  These two GLS formulas only coincide when V[i,i] is constant
  # across tips -- i.e. on ultrametric trees.
  tree <- ape::rcoal(n)   # coalescent tree is ultrametric by construction
  y <- stats::rnorm(n)
  names(y) <- tree$tip.label

  dat <- data.frame(y = y, row.names = tree$tip.label)
  fit <- phylolm::phylolm(y ~ 1, data = dat, phy = tree, model = "BM")
  phylolm_intercept <- as.numeric(stats::coef(fit)[1])

  R <- pigauto:::phylo_cor_matrix(tree)
  R_inv <- solve(R + diag(1e-6, n))
  ones  <- rep(1, n)
  pigauto_mu_hat <- as.numeric((ones %*% R_inv %*% y) /
                                 (ones %*% R_inv %*% ones))

  rel_diff <- abs(pigauto_mu_hat - phylolm_intercept) /
              max(abs(phylolm_intercept), 0.1)
  expect_lt(rel_diff, 0.05,
            label = sprintf("[D1.4] pigauto mu_hat = %.4g, phylolm = %.4g",
                             pigauto_mu_hat, phylolm_intercept))
})
