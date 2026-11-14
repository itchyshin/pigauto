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
# a 41-point grid is essentially free. Since 2026-09-24, nll(lambda) is
# full REML (the GLS mean's determinant correction is included), so this
# is now a REML posterior rather than a profile-likelihood posterior --
# the weighting still integrates the same likelihood-weighted surface,
# just with the mean-uncertainty penalty folded in.
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
# The full-REML (restricted likelihood, GLS mean profiled out) negative
# log-likelihood for Pagel's lambda on a column y observed at cells `obs`,
# with phylo correlation R restricted to R_oo, is:
#
#   nll(lambda) = 0.5 * ( (n_o - 1) * log(sigma2_hat(lambda))
#                       + log|R_oo(lambda)|
#                       + log(1' R_oo(lambda)^-1 1) )
#
# where R_oo(lambda) = lambda * R_oo + (1 - lambda) * I. The third term is
# the REML determinant correction (Harville 1977) for profiling out the
# scalar GLS mean; omitting it (as earlier versions of this file did)
# biases lambda_hat low at weak signal. Measured (200 seeds, n = 300, 30%
# MCAR, BM DGP, observed-mean shift +1.5): bias -0.061 -> -0.022 at
# lambda = 0.3, -0.038 -> -0.027 at lambda = 0.7, unchanged (-0.005) at
# lambda = 1 (2026-09-24).
#
# Eigendecompose R_oo = U Lambda U^T once (O(n_o^3)). Then for any lambda,
# d_i(lambda) = lambda * Lambda_i + (1 - lambda), and R_oo(lambda)^{-1} =
# U diag(1/d) U^T. Each NLL evaluation reduces to O(n_o) once
# c_y = U^T y_o and c_1 = U^T 1 are cached.
#
# Spec: specs/2026-05-18-pagel-lambda-eigendecomp-speedup-design.md.
#
# Optional design matrix X (n x p, full length, same row order as y):
# when supplied, the GLS coefficients beta are profiled out inside the
# NLL instead of fixing the model to an intercept-only mean, so a
# covariate-aware caller (bm_impute_col_with_cov) can estimate lambda
# jointly with beta. X = NULL reproduces the intercept-only closure
# below exactly -- that branch is untouched by the X-aware addition.
#
# Returns a list with:
#   $nll(lambda)  -- closure that evaluates NLL at any lambda in [0, 1]
#   $mu_at(lambda) -- closure returning the GLS phylogenetic mean at that
#                     lambda (intercept-only branch only; reuses the same
#                     eigenbasis as $nll, so a caller that needs both the
#                     ML lambda and its GLS mean -- e.g. the joint
#                     baseline's Henderson-centering step -- does not pay
#                     for a second O(n_o^3) factorisation of R_oo. See
#                     `ml_lambda_and_mu_for_col()` in R/bm_internal.R.
#                     NULL on the X-aware branch (no scalar mean there).
#   $n_o          -- number of observed cells (for callers that need it)
# @noRd
build_pagel_nll_cache <- function(y, R, nugget = 1e-6, X = NULL) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  if (n_o < 2L) {
    # Degenerate: NLL undefined. Return a closure that always reports +Inf
    # so optimisers correctly back off.
    return(list(
      nll = function(lambda) .Machine$double.xmax,
      mu_at = function(lambda) NA_real_,
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

  if (is.null(X)) {
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
      # Full-REML determinant correction for profiling out the scalar GLS
      # mean: 0.5 * log(1' R_oo(lambda)^-1 1) = 0.5 * log(sum_a). Guard
      # against sum_a <= 0 (should not happen once d > 0, but mirrors the
      # other singular guards in this file).
      if (!is.finite(sum_a) || sum_a <= 0) return(.Machine$double.xmax)
      0.5 * ((n_o - 1L) * log(sigma2) + log_det + log(sum_a))
    }
    mu_at <- function(lambda) {
      d <- lambda * evals + (1 - lambda) + nugget
      if (any(d <= 0)) return(NA_real_)
      inv_d <- 1 / d
      sum_a <- sum(c_1 * c_1 * inv_d)
      sum_b <- sum(c_1 * c_y * inv_d)
      sum_b / sum_a
    }
    return(list(nll = nll_at, mu_at = mu_at, n_o = n_o))
  }

  # X-aware branch: same eigenbasis, but profile a p-column GLS regression
  # instead of a scalar GLS mean. c_X = U^T X_o rotates the design matrix
  # into the same basis as c_y / c_1 above, so beta(lambda), like mu_hat(lambda)
  # above, is O(n_o) per evaluation after this one-off O(n_o * p) rotation.
  X_o <- X[obs, , drop = FALSE]
  p_x <- ncol(X_o)
  c_X <- crossprod(U, X_o)  # n_o x p_x

  nll_at_x <- function(lambda) {
    d <- lambda * evals + (1 - lambda) + nugget
    if (any(d <= 0)) return(.Machine$double.xmax)
    inv_d <- 1 / d
    # beta(lambda) = (c_X' D^-1 c_X)^-1 c_X' D^-1 c_y, D = diag(d).
    XtDinvX <- crossprod(c_X, c_X * inv_d)
    XtDinvy <- crossprod(c_X, c_y * inv_d)
    beta <- tryCatch(solve(XtDinvX, XtDinvy), error = function(e) NULL)
    if (is.null(beta) || !all(is.finite(beta))) return(.Machine$double.xmax)
    # Residual rotated to the eigenbasis: c_e = c_y - c_X %*% beta.
    c_e <- c_y - as.numeric(c_X %*% beta)
    sigma2 <- sum(c_e * c_e * inv_d) / max(n_o - p_x, 1L)
    if (!is.finite(sigma2) || sigma2 <= 0) return(.Machine$double.xmax)
    log_det <- sum(log(d))
    # Full-REML NLL: adds the log|X' R(lambda)^-1 X| = log|XtDinvX|
    # determinant correction for profiling out beta (Harville 1977), the
    # matrix analogue of the intercept-only branch's log(sum_a) term above.
    log_det_x <- tryCatch(
      determinant(XtDinvX, logarithm = TRUE)$modulus,
      error = function(e) NA_real_
    )
    if (!is.finite(log_det_x)) return(.Machine$double.xmax)
    0.5 * ((n_o - p_x) * log(sigma2) + log_det + as.numeric(log_det_x))
  }
  list(nll = nll_at_x, mu_at = NULL, n_o = n_o)
}

# ---- .pagel_lambda_from_cache -----------------------------------------------
#
# Grid-then-refine ML search over lambda in [0, 1], given an
# ALREADY-BUILT NLL cache (build_pagel_nll_cache()). Factored out of
# `ml_lambda_and_mu_for_col()` (R/bm_internal.R) so a caller that has
# already built the cache for another reason (e.g.
# `.mvn_resolve_lambda()`'s lambda_block search over the same columns,
# R/joint_mvn_solver.R) can reuse it for the per-column estimate too,
# instead of paying for a second O(n_o^3) eigendecomposition of the
# same R_oo (Rose review, 2026-09-23, "SPEED").
#
# Returns list(lambda_hat, mu_hat). mu_hat is NA_real_ when
# `cache$mu_at` is NULL (the X-aware cache branch has no scalar mean)
# or the cache is degenerate (its $nll is always +Inf).
# @noRd
.pagel_lambda_from_cache <- function(cache, lambda_grid = c(0.005, 0.995)) {
  nll   <- cache$nll
  mu_at <- cache$mu_at
  grid <- seq(lambda_grid[1L], lambda_grid[2L], length.out = 11L)
  grid_nll <- vapply(grid, nll, numeric(1L))
  if (!any(is.finite(grid_nll))) {
    return(list(lambda_hat = 1.0, mu_hat = NA_real_))
  }
  best_idx <- which.min(grid_nll)
  lo <- grid[max(1L, best_idx - 1L)]
  hi <- grid[min(length(grid), best_idx + 1L)]
  if (lo == hi) {
    lambda_hat <- grid[best_idx]
  } else {
    opt <- tryCatch(
      stats::optimise(nll, interval = c(lo, hi), tol = 1e-4),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$minimum)) {
      lambda_hat <- grid[best_idx]
    } else if (opt$objective <= grid_nll[best_idx] + 1e-8) {
      lambda_hat <- max(0, min(1, opt$minimum))
    } else {
      lambda_hat <- grid[best_idx]
    }
  }
  list(lambda_hat = lambda_hat,
       mu_hat = if (is.null(mu_at)) NA_real_ else mu_at(lambda_hat))
}

# Scale internal edges of a phylo tree by lambda, with compensation
# on terminal edges that preserves root-to-tip distances.
#
# Concretely: each terminal edge t (parent height h_p, child = tip)
# is replaced by
#   t' = t + (1 - lambda) * h_p
# and each internal edge i is scaled by lambda.
#
# This is exact on the CORRELATION scale for any tree, ultrametric or
# not: for tips i, j with MRCA depth d_ij, every edge on the root-to-MRCA
# path is internal, so vcv(out)[i, j] = lambda * d_ij = lambda * A[i, j].
# Root-to-tip depth is preserved by construction (lambda * h_p + t' =
# h_p + t = A[i, i]), so the diagonal stays at A[i, i] and
#   cov2cor(vcv(out)) = lambda * cov2cor(vcv(tree)) + (1 - lambda) * I
# exactly, for any tree topology or branch lengths.
#
# On the COVARIANCE scale the identity only holds approximately for a
# non-ultrametric tree: vcv(out) = lambda * A + (1 - lambda) * diag(A),
# not lambda * A + (1 - lambda) * I, because diag(A) is not constant
# when tips are not equidistant from the root. The solver here always
# works on the correlation scale, so this does not affect it, but a
# caller reading vcv(out) directly on a non-ultrametric tree should be
# aware of the distinction.
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
