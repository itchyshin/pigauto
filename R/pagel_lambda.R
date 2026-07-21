# Pagel's lambda BM baseline (v0.10).
#
# Spec: specs/2026-05-18-pagel-lambda-baseline-design.md.
#
# Pagel's lambda is a single-parameter generalization of Brownian motion
# that shrinks the off-diagonal of the phylogenetic correlation matrix
# toward zero:
#
#   R(lambda) = lambda * R + (1 - lambda) * I
#
# This is equivalent to scaling all INTERNAL branch lengths by lambda
# while leaving terminal (tip-incident) edges unchanged: under that
# transform, off-diagonal entries of vcv(tree) scale by lambda while the
# diagonal (root-to-tip total) is preserved. The equivalence is exact for
# ultrametric trees and a close approximation for non-ultrametric trees;
# we rely on the tree-transform form throughout because it keeps the
# downstream sparse Hadfield-Nakagawa machinery working unchanged.
#
# Reference: Pagel (1999) Nature 401: 877-884.

# ---- bayes_lambda_for_col -------------------------------------------------
#
# Bayesian model-averaged Pagel's lambda (v0.11). Returns the posterior
# distribution over lambda under a flat prior, computed deterministically
# from the profile log-likelihood evaluated on a fine grid.
#
# Why: ML (`ml_lambda_for_col`) and k-fold CV (`cv_lambda_for_col`) both
# collapse to a single point lambda_hat. On traits with weak phylogenetic
# signal, the profile likelihood is flat across a wide interior region
# and a single point estimate is mask-sensitive — different masks give
# wildly different lambda_hat values. Bayesian model averaging integrates
# over the whole likelihood-weighted lambda distribution, so the
# downstream prediction is stable.
#
# Mathematically: w_i propto exp(-nll(lambda_i)), normalised. Under a
# flat prior over [0, 1] this is exactly the marginal posterior on a
# coarse grid; the eigendecomp cache makes evaluation cheap enough that
# a 41-point grid is essentially free.
#
# This is the deterministic analogue of what BACE's MCMCglmm does via
# Markov-chain samples from the phylo / residual variance-component
# posterior. The integration is identical; the sampling is just a fast
# numerical Riemann sum rather than a chain.
#
# Spec: specs/2026-05-18-pagel-lambda-baseline-design.md (v0.11 update).
#
# Returns a list with:
#   $lambda_grid          numeric vector of grid points
#   $weights              numeric vector of posterior weights (sums to 1)
#   $lambda_post_mean     scalar -- posterior mean of lambda
#   $lambda_post_entropy  scalar -- Shannon entropy of the posterior
#                         (in nats; high entropy = wide posterior, low =
#                         concentrated)
#
# Degenerate cases (n_obs < 10 or all-NLL non-finite) return a point
# mass at lambda = 1 with entropy 0, matching the ml_lambda fallback.
# @noRd
bayes_lambda_for_col <- function(y, R, nugget = 1e-6,
                                   lambda_grid = seq(0.005, 0.995,
                                                       length.out = 41L)) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  point_mass_1 <- function() {
    list(lambda_grid = 1.0,
         weights = 1.0,
         lambda_post_mean = 1.0,
         lambda_post_entropy = 0.0)
  }
  if (n_o < 10L) return(point_mass_1())

  cache <- build_pagel_nll_cache(y, R, nugget = nugget)
  nll_vals <- vapply(lambda_grid, cache$nll, numeric(1L))
  ok <- is.finite(nll_vals)
  if (!any(ok)) return(point_mass_1())

  # Softmax with max-subtract for numerical stability. Non-finite NLLs
  # (Cholesky failures at extreme lambdas) get zero weight.
  log_w <- -nll_vals
  log_w[!ok] <- -Inf
  log_w <- log_w - max(log_w[ok])
  w <- exp(log_w)
  w_sum <- sum(w)
  if (!is.finite(w_sum) || w_sum <= 0) return(point_mass_1())
  w <- w / w_sum

  # Posterior summaries.
  post_mean <- sum(w * lambda_grid)
  # Shannon entropy in nats; mask near-zero weights to avoid 0*log(0).
  w_safe <- w[w > 1e-12]
  post_entropy <- -sum(w_safe * log(w_safe))

  list(lambda_grid = lambda_grid,
       weights = w,
       lambda_post_mean = post_mean,
       lambda_post_entropy = post_entropy)
}

# ---- cv_lambda_for_col ----------------------------------------------------
#
# Cross-validated Pagel's lambda selection (v0.11). Direct fix for the
# weak-signal ML pathology: ML picks lambda near 0 when the profile
# likelihood is flat there, even when the RMSE-optimal interior lambda
# would generalise better. CV optimises out-of-fold prediction error
# directly, so it recovers the right interior lambda.
#
# Spec: specs/2026-05-18-cv-lambda-selection-design.md.
#
# Falls back to lambda = 1 (matches ml_lambda_for_col convention) when
# n_obs < 50 (folds become too small to be useful).
#
# Returns a scalar lambda_hat in [0, 1].
# @noRd
cv_lambda_for_col <- function(y, R, nugget = 1e-6,
                                lambda_grid = seq(0.01, 0.99, by = 0.05),
                                k = 5L, seed = NULL) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  if (n_o < 50L) return(1.0)

  # Seeded fold assignment over observed cells.
  if (!is.null(seed)) set.seed(seed)
  fold_id <- sample(rep(seq_len(k), length.out = n_o))

  # Per-lambda cumulative squared error across folds.
  cv_sse <- numeric(length(lambda_grid))
  for (f in seq_len(k)) {
    held <- obs[fold_id == f]
    if (length(held) == 0L) next
    y_train <- y
    y_train[held] <- NA
    truth <- y[held]
    for (i in seq_along(lambda_grid)) {
      pred <- bm_impute_col(y_train, R, nugget = nugget,
                              lambda = lambda_grid[i])
      cv_sse[i] <- cv_sse[i] + sum((pred$mu[held] - truth)^2)
    }
  }
  if (!any(is.finite(cv_sse))) return(1.0)
  lambda_grid[which.min(cv_sse)]
}

# ---- build_pagel_nll_cache --------------------------------------------------
#
# One-time eigendecomposition cache for the Pagel's lambda profile NLL.
#
# The profile-REML negative log-likelihood for Pagel's lambda on a column y
# observed at cells `obs`, with phylo correlation R restricted to R_oo, is:
#
#   nll(lambda) = 0.5 * ( (n_o - 1) * log(sigma2_hat(lambda))
#                       + log|R_oo(lambda)| )
#
# where R_oo(lambda) = lambda * R_oo + (1 - lambda) * I.
#
# Eigendecompose R_oo = U Lambda U^T once (O(n_o^3)). Then for any lambda,
# d_i(lambda) = lambda * Lambda_i + (1 - lambda), and R_oo(lambda)^{-1} =
# U diag(1/d) U^T. Each NLL evaluation reduces to O(n_o) once
# c_y = U^T y_o and c_1 = U^T 1 are cached.
#
# Spec: specs/2026-05-18-pagel-lambda-eigendecomp-speedup-design.md.
#
# Returns a list with:
#   $nll(lambda)  -- closure that evaluates NLL at any lambda in [0, 1]
#   $n_o          -- number of observed cells (for callers that need it)
# @noRd
build_pagel_nll_cache <- function(y, R, nugget = 1e-6) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  if (n_o < 2L) {
    # Degenerate: NLL undefined. Return a closure that always reports +Inf
    # so optimisers correctly back off.
    return(list(
      nll = function(lambda) .Machine$double.xmax,
      n_o = n_o
    ))
  }
  R_oo <- R[obs, obs, drop = FALSE]
  y_o  <- y[obs]
  eig <- eigen(R_oo, symmetric = TRUE)
  U <- eig$vectors
  # Clip tiny negative eigenvalues from numerical noise; bound away from
  # zero so adding (1 - lambda) keeps d strictly positive even at lambda = 1.
  evals <- pmax(eig$values, 0)
  c_y <- as.numeric(crossprod(U, y_o))
  c_1 <- as.numeric(crossprod(U, rep(1, n_o)))

  nll_at <- function(lambda) {
    # d_i(lambda) = lambda * Lambda_i + (1 - lambda), plus a tiny nugget on
    # the diagonal of R_oo(lambda) for parity with the dense-Cholesky path.
    d <- lambda * evals + (1 - lambda) + nugget
    if (any(d <= 0)) return(.Machine$double.xmax)
    inv_d <- 1 / d
    # GLS phylogenetic mean. In the original basis:
    #   sum_a = 1^T R_oo(lambda)^{-1} 1
    #         = 1^T U diag(1/d) U^T 1
    #         = c_1^T diag(1/d) c_1 = sum(c_1^2 * inv_d)
    #   sum_b = 1^T R_oo(lambda)^{-1} y_o
    #         = c_1^T diag(1/d) c_y = sum(c_1 * inv_d * c_y)
    sum_a <- sum(c_1 * c_1 * inv_d)
    sum_b <- sum(c_1 * c_y * inv_d)
    mu_hat <- sum_b / sum_a
    # Residual rotated to the eigenbasis: c_e = U^T (y - mu*1) = c_y - mu*c_1.
    c_e <- c_y - mu_hat * c_1
    # REML variance.
    sigma2 <- sum(c_e * c_e * inv_d) / max(n_o - 1L, 1L)
    if (!is.finite(sigma2) || sigma2 <= 0) return(.Machine$double.xmax)
    log_det <- sum(log(d))
    0.5 * ((n_o - 1L) * log(sigma2) + log_det)
  }

  list(nll = nll_at, n_o = n_o)
}

# Scale internal edges of a phylo tree by lambda, with compensation
# on terminal edges that preserves root-to-tip distances for
# ultrametric trees.
#
# Concretely: each terminal edge t (parent height h_p, child = tip)
# is replaced by
#   t' = t + (1 - lambda) * h_p
# and each internal edge i is scaled by lambda. The result is that
# cov2cor(vcv(out)) equals lambda * cov2cor(vcv(tree)) + (1-lambda)*I
# exactly (when the tree is ultrametric).
#
# For non-ultrametric trees the formula is a close approximation but
# not exact; the diagonal of the new correlation matrix still equals
# one, while off-diagonals are scaled approximately by lambda. Users
# with strongly non-ultrametric trees should be aware.
#
# @param tree object of class "phylo".
# @param lambda numeric scalar in [0, 1].
# @return a phylo with edge lengths transformed.
# @noRd
transform_tree_pagel <- function(tree, lambda) {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0 || lambda > 1) {
    stop("'lambda' must be a numeric scalar in [0, 1]; got: ",
         paste(lambda, collapse = ", "), call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  # node.depth.edgelength: distance from root for every node, ordered
  # so node ids 1..n_tip are tips and n_tip+1..n_tip+n_int are internal.
  node_depths <- ape::node.depth.edgelength(tree)
  parents <- tree$edge[, 1L]
  children <- tree$edge[, 2L]
  is_terminal <- children <= n_tip

  out <- tree
  # Internal edges: scale by lambda.
  out$edge.length[!is_terminal] <-
    lambda * tree$edge.length[!is_terminal]
  # Terminal edges: t' = t + (1-lambda) * parent_depth, where
  # parent_depth = root-to-parent distance under the original tree.
  parent_depth_term <- node_depths[parents[is_terminal]]
  out$edge.length[is_terminal] <-
    tree$edge.length[is_terminal] + (1 - lambda) * parent_depth_term
  out
}
