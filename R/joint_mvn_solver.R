# R/joint_mvn_solver.R
#
# In-house ML estimator of the K x K trait covariance matrix Sigma
# under a matrix-normal Brownian motion model:
#
#   vec(L) ~ MVN(0, Sigma %x% R)
#
# where L is the n x K liability matrix (rows = tips, columns =
# liability dimensions; NAs allowed) and R is the n x n phylogenetic
# correlation matrix (cov2cor(vcv(tree))). Replaces
# Rphylopars::phylopars() for pigauto's joint MVN baseline,
# threshold-joint baseline, and OVR categorical fits.
#
# Output contract matches what phylopars exposes via the fields
# pigauto actually consumes:
#
#   $anc_recon : n x K matrix of posterior cell means (tip rows only;
#                rownames match input rownames).
#   $anc_var   : n x K matrix of posterior cell variances. 0 at
#                observed cells, > 0 at imputed cells.
#   $pars$phylocov : K x K Sigma estimate.
#
# Algorithm (Fix B: Fisher-ML on observed cells + EM for L_hat).
#
#   1. Initialise Sigma_hat from observed-cell pairwise-complete sample
#      covariance. Captures cross-trait correlation from iteration 0.
#   2. Maximise the observed-data marginal log-likelihood with respect
#      to a Cholesky parameterisation of Sigma. The likelihood treats
#      each species's observed traits as iid MVN(0, Sigma_obs,obs);
#      phylogeny (R) enters separately at the cell-imputation step
#      below. This is the matrix-normal Sigma MLE that phylopars
#      implements without the full Kronecker pruning machinery.
#   3. Per-column BM init of L_hat / L_var (uses R via bm_impute_col).
#   4. EM cell-refinement: inverse-variance pool between per-column BM
#      posterior and cross-trait conditional MVN posterior under the
#      Fisher-ML Sigma.
#
# Why Fisher-ML beats my previous naive (1/n) L_hat^T R^-1 L_hat init:
# the previous init imputed each column independently before forming
# Sigma_hat, so the imputed values had zero cross-trait correlation by
# construction. Sigma_hat off-diagonals were biased toward zero, decoded
# OVR-categorical class probabilities clustered near 1/K, accuracy
# collapsed to near-random. Fisher-ML on observed cells avoids that
# imputation pollution entirely.

# Force a symmetric matrix to be positive-definite via eigen-clip.
.mvn_ensure_pd <- function(M, eps = 1e-8) {
  if (!isSymmetric(M, tol = 100 * eps)) M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  evals <- pmax(ev$values, eps)
  K <- length(evals)
  ev$vectors %*% diag(evals, nrow = K) %*% t(ev$vectors)
}

# ---- Fisher-ML Sigma (opt-in `sigma_method = "fisher_ml"`) ---------------
#
# Ported (adapted, not cherry-picked) from commit e7ca41c (2026-05-16,
# "autoresearch exp-4 (B-1): Fisher-ML Sigma via Cholesky-optim
# observed-data NLL"), never merged -- see docs/dev-log/
# 2026-08-16-continuous-gap-diagnosis.md for why it was revived: the
# default single-pass solver loses 0.14-1.27 z-RMSE to converged
# Rphylopars REML on AVONET300. `sigma_method = "single_pass"` (default)
# is byte-identical to pre-existing behaviour; the machinery below only
# executes under `sigma_method = "fisher_ml"`.
#
# The prototype's own Sigma-init strategy (species-iid pairwise-complete
# sample covariance) is exactly what the v0.9.2 "Bug fixes: in-house
# Sigma solver" NEWS entry later replaced with the R-credited Kronecker
# M-step (`.mvn_sigma_kron_M()` below), because species-iid double-counts
# correlated tips and produced near-singular Sigma on strongly
# phylo-conserved data. That fix owns the single_pass default; fisher_ml
# is a separate, self-contained alternative a caller opts into, not a
# resurrection of the pre-fix default. It leaves the *other* two v0.9.2
# fixes untouched: the per-column BM init (`.mvn_init_per_column()`)
# still goes through Henderson `cor_scale = TRUE` when available, and
# `max_iter` still defaults to 0L (EM cell-refinement stays opt-in for
# both sigma_method values, for the same near-sister-tip divergence
# reason documented in NEWS.md).

# Pack/unpack the K(K+1)/2 free parameters of a K x K SPD Sigma via its
# lower-triangular Cholesky factor LC, with LC[j, j] = exp(par[diag_idx])
# so diagonals stay positive without constrained optimisation.
.mvn_par_to_chol <- function(par, K) {
  LC <- matrix(0, K, K)
  idx <- 1L
  for (j in seq_len(K)) {
    LC[j, j] <- exp(par[idx]); idx <- idx + 1L
    if (j < K) for (i in (j + 1L):K) {
      LC[i, j] <- par[idx]; idx <- idx + 1L
    }
  }
  LC
}

.mvn_chol_to_par <- function(LC) {
  K <- ncol(LC)
  par <- numeric(K * (K + 1L) / 2L)
  idx <- 1L
  for (j in seq_len(K)) {
    par[idx] <- log(max(LC[j, j], 1e-12)); idx <- idx + 1L
    if (j < K) for (i in (j + 1L):K) {
      par[idx] <- LC[i, j]; idx <- idx + 1L
    }
  }
  par
}

# Observed-data negative log-likelihood for L | Sigma assuming each
# species is iid MVN(0, Sigma) (phylogenetic R coupling is dropped from
# this term; it re-enters at the per-column-BM init and EM cell-
# refinement steps). Drops constant terms. Lower is better.
.mvn_obs_nll <- function(par, L, K, eps = 1e-8) {
  LC <- .mvn_par_to_chol(par, K)
  Sigma <- LC %*% t(LC)
  total <- 0
  for (i in seq_len(nrow(L))) {
    obs <- !is.na(L[i, ])
    n_obs <- sum(obs)
    if (n_obs == 0L) next
    S_oo <- Sigma[obs, obs, drop = FALSE] + diag(eps, n_obs)
    chol_S <- tryCatch(chol(S_oo), error = function(e) NULL)
    if (is.null(chol_S)) return(.Machine$double.xmax)
    log_det <- 2 * sum(log(diag(chol_S)))
    li <- as.numeric(L[i, obs])
    inv_li <- chol2inv(chol_S) %*% li
    quad <- sum(li * inv_li)
    total <- total + 0.5 * (log_det + quad)
  }
  total
}

# Fisher-ML refinement of Sigma starting from Sigma_start. Uses optim()
# over the Cholesky-parameterised free parameters (K is typically 2-8 in
# pigauto's joint paths, so the K(K+1)/2 free parameters keep optim
# cheap). `fallback_fn()` is the zero-argument "single_pass" Sigma this
# call site would otherwise have produced; it is invoked -- with a
# warning -- when optim() errors or fails to converge, so fisher_ml can
# never leave the caller with a worse Sigma than the default.
.mvn_sigma_fisher_ml <- function(L, Sigma_start, fallback_fn, eps = 1e-8,
                                  optim_maxit = 50L) {
  K <- ncol(L)
  Sigma_start_pd <- .mvn_ensure_pd(Sigma_start, eps = eps)
  LC_init <- t(chol(Sigma_start_pd))
  par_init <- .mvn_chol_to_par(LC_init)
  fit <- tryCatch(
    stats::optim(par_init, .mvn_obs_nll, L = L, K = K, eps = eps,
                 method = "BFGS", control = list(maxit = optim_maxit)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0L) {
    warning("fit_mvn_bm_inhouse: sigma_method = \"fisher_ml\" optim() did ",
            "not converge; falling back to the single_pass Sigma estimate.",
            call. = FALSE)
    return(fallback_fn())
  }
  LC <- .mvn_par_to_chol(fit$par, K)
  LC %*% t(LC)
}

# Proper matrix-normal Kronecker M-step using Hadfield-Nakagawa sparse
# R^{-1}. Replaces the previous iid-species approximation when a tree
# is available. Given the current E-step's posterior mean L_hat and
# per-cell variance L_var:
#
#   Sigma_hat = (1/n) (L_hat^T R^{-1} L_hat + V_correction)
#
# where V_correction = diag(colSums of cell-level posterior variances).
# This is the closed-form M-step for matrix-normal under the standard
# EM decomposition with vec(L) ~ MVN(0, Sigma %x% R). Using Henderson
# Q_S in place of dense R^{-1} keeps it O(K * n) instead of O(n^3).
.mvn_sigma_kron_M <- function(L_hat, L_var, henderson) {
  K <- ncol(L_hat); n <- nrow(L_hat)
  # Joint MVN solver works on the correlation scale R = cov2cor(vcv(tree)),
  # so request R^{-1} (not A^{-1}) from Henderson. The cor_scale = TRUE
  # rescale uses sqrt(diag(A)) stored in `henderson`.
  Rinv_L <- henderson_R_inv_apply(L_hat, henderson, cor_scale = TRUE)
  M <- crossprod(L_hat, Rinv_L) + diag(colSums(L_var), nrow = K)
  M / n
}

# Observed-cell GLS phylogenetic mean under R(lambda), for one column:
#
#   mu_hat(lambda) = (1' R_oo(lambda)^{-1} 1)^{-1} 1' R_oo(lambda)^{-1} y_o
#
# This is exactly bm_impute_col()'s own `mu_hat` (bm_internal.R, the
# "1. GLS phylogenetic mean" step), duplicated here rather than exposed
# from R/bm_internal.R or R/pagel_lambda.R (both owned by a concurrent
# lane in this worktree) so this file stays self-contained. Used only to
# give henderson_bm_predict() the same free-mean model that
# bm_impute_col() and the lambda_k profile use, since henderson_bm_predict
# itself assumes a zero root state (mean-model consistency, S2
# correction / Rose review section 7 B(i)).
.mvn_gls_mean_at_lambda <- function(y, R, lambda, nugget = 1e-8) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  if (n_o < 1L) return(0)
  R_oo <- R[obs, obs, drop = FALSE]
  if (lambda < 1) {
    diag_R <- diag(R_oo)
    R_oo <- lambda * R_oo
    diag(R_oo) <- lambda * diag_R + (1 - lambda)
  }
  R_oo <- R_oo + diag(nugget, n_o)
  chol_R <- tryCatch(chol(R_oo), error = function(e) NULL)
  if (is.null(chol_R)) return(mean(y[obs]))
  chol_solve <- function(b) backsolve(chol_R, forwardsolve(t(chol_R), b))
  ones <- rep(1, n_o)
  a <- chol_solve(ones)
  b <- chol_solve(y[obs])
  sum(b) / sum(a)
}

# Per-column BM init: exact univariate-BM conditional MVN.
# Returns L_hat (NAs replaced) and L_var (0 at observed cells, > 0 at
# imputed cells). For columns with < 2 observed values, leaves the
# column at zero with unit variance (uninformative prior).
#
# Two code paths share the same numerical answer to ~1e-3:
#   - bm_impute_col(yj, R, lambda = lambda_vec[j])  legacy dense
#     O(n_obs^3) per column, lambda-aware, free GLS mean throughout.
#   - henderson_bm_predict(yj, H)  sparse O(n) per column via Hadfield-
#     Nakagawa (2010) eq 29, where H is already built on the Pagel-
#     transformed tree T(lambda_vec[j]). This path assumes a ZERO root
#     state, so when lambda_vec[j] != 1 the column is first centered at
#     its own .mvn_gls_mean_at_lambda() before the Henderson call and
#     the mean is added back to the result -- otherwise the sparse and
#     dense paths disagree by the observed-subset mean as lambda falls
#     (measured by Rose's review, section 7 B(i): max |delta mu| grows
#     from ~0.001 at lambda=1 to ~0.25 at lambda=0.3 on non-centered
#     data). At lambda_vec[j] == 1 no centering happens, matching the
#     pre-lane path exactly. Cuts init time at large n from minutes to
#     seconds (~18x at n=2000) and avoids forming dense (n x n) R^{-1}.
#
# `henderson_for_col(j)` returns the Henderson object to use for column
# j (already built at that column's lambda), or NULL to force the dense
# fallback (K == 1, no tree, Matrix unavailable, or use_henderson =
# FALSE). `lambda_vec` (length K) is consulted on both paths.
#
# `mu_hat_vec` (optional, length K, NA where not precomputed): the GLS
# mean at column j's OWN lambda_j, when the caller already has it from
# the same eigendecomposition used to estimate lambda_j (S6 perf fix,
# Rose review 2026-09-23: `.mvn_resolve_lambda()`'s "estimate" branch
# gets this for free from `ml_lambda_and_mu_for_col()`). When
# `mu_hat_vec[j]` is finite it is used as-is; otherwise this function
# falls back to its own `.mvn_gls_mean_at_lambda()` call, unchanged.
.mvn_init_per_column <- function(L, R, eps = 1e-8, henderson_for_col = NULL,
                                  lambda_vec = NULL, mu_hat_vec = NULL) {
  n <- nrow(L); K <- ncol(L)
  if (is.null(lambda_vec)) lambda_vec <- rep(1, K)
  if (is.null(mu_hat_vec)) mu_hat_vec <- rep(NA_real_, K)
  L_hat <- L; L_var <- matrix(0, n, K, dimnames = dimnames(L))
  for (j in seq_len(K)) {
    yj <- L[, j]
    obs <- !is.na(yj)
    if (sum(obs) < 2L) {
      L_hat[, j] <- ifelse(is.na(yj), 0, yj)
      L_var[!obs, j] <- 1
      next
    }
    henderson_j <- if (!is.null(henderson_for_col)) henderson_for_col(j) else NULL
    lam_j <- lambda_vec[j]
    if (!is.null(henderson_j) && lam_j != 1) {
      mu_gls <- if (is.finite(mu_hat_vec[j])) {
        mu_hat_vec[j]
      } else {
        .mvn_gls_mean_at_lambda(yj, R, lam_j, nugget = eps)
      }
      res <- henderson_bm_predict(yj - mu_gls, henderson_j, eps = eps,
                                    cor_scale = TRUE)
      res$mu <- res$mu + mu_gls
    } else if (!is.null(henderson_j)) {
      res <- henderson_bm_predict(yj, henderson_j, eps = eps, cor_scale = TRUE)
    } else {
      res <- bm_impute_col(yj, R, nugget = eps, lambda = lam_j)
    }
    L_hat[, j] <- res$mu
    L_var[, j] <- res$se^2
    L_var[obs, j] <- 0
  }
  list(L_hat = L_hat, L_var = L_var)
}

# ---- Pagel's lambda resolution (S2, feat/joint-lambda-default) -----------
#
# Design: docs/dev-log/2026-09-22-joint-lambda-alignment.md, sections 2-7.
# `lambda_cols` names the "continuous-family" columns (continuous, count,
# ordinal, proportion, zi magnitude); NULL means all columns are eligible.
# Discrete liability columns (binary, zi gate, ordinal-via-OVR synthetic
# columns) are excluded by the caller and stay at lambda = 1 -- section 7's
# post-review decision (B iii cut): the spec's decision 4 describes a
# lambda estimated over the WHOLE liability block, which this lane does
# not build, and `tests/testthat/test-lambda-per-type.R` already locks in
# lambda = 1 on the discrete path (the August arc/lambda-per-type fix).
# `lambda_block` (the block value below) is used only by the Sigma
# M-step, the opt-in exact conditional, and the opt-in EM refine -- never
# to override a non-lambda_cols column's own lambda.
.mvn_resolve_lambda_cols <- function(lambda_cols, K) {
  if (is.null(lambda_cols)) return(seq_len(K))
  if (is.logical(lambda_cols)) {
    if (length(lambda_cols) != K) {
      stop("fit_mvn_bm_inhouse: logical 'lambda_cols' must have length ",
           "ncol(L).", call. = FALSE)
    }
    return(which(lambda_cols))
  }
  lambda_cols <- as.integer(lambda_cols)
  if (anyNA(lambda_cols) || any(lambda_cols < 1L) || any(lambda_cols > K)) {
    stop("fit_mvn_bm_inhouse: 'lambda_cols' indices out of range.",
         call. = FALSE)
  }
  lambda_cols
}

# Resolves `lambda` ("fixed_1", "estimate", or a numeric scalar in
# [0, 1]) into a per-column vector, a single block value, and a tag
# recording which mode ran. "fixed_1" takes a literal shortcut (no
# caches, no optim()) so it is identical to the pre-lane solver.
#
# "estimate": lambda_bar (the block value) is the argmin over [0.01,
# 0.99] of the SUM of the profile-REML NLL caches of the lambda_cols
# columns (build_pagel_nll_cache(), reused from R/pagel_lambda.R). Each
# lambda_cols column with >= 10 observed cells then gets its own
# per-trait lambda_k via ml_lambda_for_col() (R/bm_internal.R; reused,
# not duplicated); a lambda_cols column with < 10 observed cells falls
# back to lambda_bar (too little data for its own estimate). Every
# column OUTSIDE lambda_cols is fixed at lambda = 1 (section 7's B iii
# cut) -- lambda_bar is never imposed on a column the caller did not
# name.
.mvn_resolve_lambda <- function(lambda, L, R, lambda_cols_idx, eps = 1e-8) {
  K <- ncol(L)
  col_names <- colnames(L)
  if (identical(lambda, "fixed_1")) {
    return(list(lambda_vec = stats::setNames(rep(1, K), col_names),
                lambda_block = 1,
                lambda_mode_used = "fixed_1"))
  }
  if (is.numeric(lambda) && length(lambda) == 1L) {
    if (!is.finite(lambda) || lambda < 0 || lambda > 1) {
      stop("fit_mvn_bm_inhouse: numeric 'lambda' must be a scalar in ",
           "[0, 1]; got: ", lambda, call. = FALSE)
    }
    return(list(lambda_vec = stats::setNames(rep(lambda, K), col_names),
                lambda_block = lambda,
                lambda_mode_used = "numeric"))
  }
  if (is.numeric(lambda) && length(lambda) == K) {
    # S4 (dispatcher, predict-time lambda_fixed rebuild): a full per-column
    # vector supplied by the caller. Every entry is used as-is -- lambda_cols
    # is irrelevant here because the caller has already decided every
    # column's value (typically by replaying a previous fit's own
    # $lambda_per_trait). `lambda_block` is a diagnostic mean, only
    # consulted by the opt-in Sigma M-step / exact conditional / EM refine.
    if (!all(is.finite(lambda)) || any(lambda < 0) || any(lambda > 1)) {
      stop("fit_mvn_bm_inhouse: numeric 'lambda' vector must have all ",
           "entries in [0, 1].", call. = FALSE)
    }
    lambda_vec <- stats::setNames(as.numeric(lambda), col_names)
    return(list(lambda_vec = lambda_vec,
                lambda_block = mean(lambda_vec),
                lambda_mode_used = "numeric_vector"))
  }
  if (!identical(lambda, "estimate")) {
    stop("fit_mvn_bm_inhouse: 'lambda' must be \"fixed_1\", \"estimate\", ",
         "or a numeric scalar in [0, 1]; got: ", paste(lambda, collapse = ", "),
         call. = FALSE)
  }

  # S6 perf fix (Rose review, 2026-09-23, "SPEED"): build each
  # lambda_cols_idx column's eigendecomposition cache ONCE and reuse it
  # for BOTH the lambda_block search below and this column's own
  # per-trait lambda_k / GLS mean, via `.pagel_lambda_from_cache()`
  # (R/pagel_lambda.R). Before this fix, `ml_lambda_for_col()` rebuilt
  # the same O(n_o^3) cache a second time per column, and
  # `.mvn_init_per_column()`'s Henderson-centering step rebuilt it a
  # third time (as a dense Cholesky) via `.mvn_gls_mean_at_lambda()`.
  caches <- lapply(lambda_cols_idx, function(j) {
    build_pagel_nll_cache(L[, j], R, nugget = eps)
  })
  lambda_block <- if (length(caches) == 0L) {
    1.0
  } else {
    block_obj <- function(lam) {
      sum(vapply(caches, function(cc) cc$nll(lam), numeric(1L)))
    }
    opt <- tryCatch(stats::optimize(block_obj, interval = c(0.01, 0.99)),
                     error = function(e) NULL)
    if (is.null(opt) || !is.finite(opt$objective) || !is.finite(opt$minimum)) {
      1.0
    } else {
      opt$minimum
    }
  }

  # Every column starts at lambda = 1 (the B iii cut: nothing outside
  # lambda_cols is ever touched), then lambda_cols columns are
  # overwritten with their own estimate or, for a low-n column, with
  # lambda_bar. `mu_hat_vec` carries each column's GLS mean AT its own
  # lambda_k, read off the same cache, for `.mvn_init_per_column()` to
  # reuse (NA where not computed, e.g. low-n columns using lambda_block).
  lambda_vec <- rep(1, K)
  mu_hat_vec <- rep(NA_real_, K)
  for (idx in seq_along(lambda_cols_idx)) {
    j <- lambda_cols_idx[idx]
    n_obs_j <- sum(!is.na(L[, j]))
    if (n_obs_j >= 10L) {
      lm_j <- .pagel_lambda_from_cache(caches[[idx]])
      lambda_vec[j] <- lm_j$lambda_hat
      mu_hat_vec[j] <- lm_j$mu_hat
    } else {
      lambda_vec[j] <- lambda_block
    }
  }
  names(lambda_vec) <- col_names

  list(lambda_vec = lambda_vec, lambda_block = lambda_block,
       lambda_mode_used = "estimate", mu_hat_vec = mu_hat_vec)
}

# E-step: refine imputations at originally-NA cells via inverse-variance
# pooling between the per-column BM posterior and the cross-trait
# conditional MVN posterior (uses Sigma off-diagonals via
# build_conditional_prior).
#
# `refine_variance` controls what happens to L_var when the mean moves:
#
#   "conservative" (default): update the MEAN only; leave L_var at the
#     per-column BM posterior variance. Rationale: `prec_bm + prec_cross`
#     is the precision of two INDEPENDENT estimates, but the BM posterior
#     and the cross-trait conditional posterior are both functions of the
#     same observed cells -- they are strongly dependent, so adding their
#     precisions double-counts the information and understates the
#     variance. Measured (2026-08-17 recovery sim,
#     docs/dev-log/2026-08-17-sigma-recovery-results.md): the pooled rule
#     drives 95% SE coverage from 0.925 down to 0.856 at one iteration and
#     0.618 at three, while the refined MEAN genuinely improves RMSE. So
#     the mean update is kept and the variance is held at a value we can
#     defend. Because the refined mean is more accurate at the same
#     variance, this is conservative (coverage moves UP, not down).
#
#   "pooled": the historical `1 / (prec_bm + prec_cross)` rule. Retained
#     only for reproducing pre-2026-08-17 behaviour and for the recovery
#     sim's comparison arm. Not recommended: its intervals are
#     overconfident by construction.
#
# The exact conditional variance under vec(L) ~ MVN(0, Sigma %x% R) is what
# a full joint solve returns; neither rule here computes it. "conservative"
# errs toward over-coverage, which is the safe direction for a package
# whose headline UQ claim is interval validity.
.mvn_estep_refine <- function(L_obs_mask, L_hat, L_var, Sigma, eps = 1e-8,
                               refine_variance = c("conservative", "pooled")) {
  refine_variance <- match.arg(refine_variance)
  cross <- build_conditional_prior(Sigma, L_hat, eps = eps)
  L_hat_new <- L_hat
  L_var_new <- L_var
  K <- ncol(L_hat)
  for (j in seq_len(K)) {
    idx <- which(!L_obs_mask[, j])
    if (!length(idx)) next
    v_bm    <- pmax(L_var[idx, j], eps)
    m_bm    <- L_hat[idx, j]
    v_cross <- pmax(cross$sd_prior[idx, j]^2, eps)
    m_cross <- cross$mu_prior[idx, j]
    prec_bm    <- 1 / v_bm
    prec_cross <- 1 / v_cross
    prec_tot   <- prec_bm + prec_cross
    L_hat_new[idx, j] <- (m_bm * prec_bm + m_cross * prec_cross) / prec_tot
    if (identical(refine_variance, "pooled")) {
      L_var_new[idx, j] <- 1 / prec_tot
    }
    # "conservative": L_var_new[idx, j] stays at the BM value.
  }
  list(L_hat = L_hat_new, L_var = L_var_new)
}

# Main entry point. Inputs:
#   L     : n x K liability matrix with NAs at unobserved cells.
#           rownames must match tip labels of `tree` (or of `R`).
#   tree  : phylo. Optional if R is supplied.
#   R     : n x n phylogenetic correlation matrix. Computed from
#           `tree` if missing.
#   max_iter : maximum EM iterations for the L_hat refinement loop
#              (default 5).
#   tol   : relative-Frobenius stopping criterion on Sigma.
#   eps   : ridge added to covariance solves.
#   sigma_method : "single_pass" (default, byte-identical to
#              pre-existing behaviour) or "fisher_ml" (opt-in; see the
#              "Fisher-ML Sigma" comment block above `.mvn_par_to_chol`
#              for provenance and the fallback-on-non-convergence
#              contract).
#   lambda : "fixed_1" (default, identical to the pre-lane solver -- no
#              new code path executes), "estimate" (per-trait Pagel's
#              lambda via profile REML on `lambda_cols`; every other
#              column stays at lambda = 1; see `.mvn_resolve_lambda()`),
#              a numeric scalar in [0, 1] applied to every column, or a
#              numeric vector of length `ncol(L)` giving each column's
#              own fixed value directly (S4 / dispatcher predict-time
#              rebuild: typically a previous fit's own
#              `$lambda_per_trait`, replayed rather than re-estimated).
#   lambda_cols : integer/logical index into columns of `L` naming the
#              continuous-family columns eligible for their own
#              lambda_k under `lambda = "estimate"`. NULL (default)
#              means all columns. Ignored for "fixed_1" and any numeric
#              `lambda` (every column's value is already decided).
#
# Returns a list with the phylopars-compatible fields described above,
# plus diagnostics ($n_iter, $converged) and $lambda_per_trait,
# $lambda_block, $lambda_mode_used (see docs/dev-log/
# 2026-09-22-joint-lambda-alignment.md sections 4, 6 and 7).
# $lambda_block is a diagnostic only: it feeds the Sigma M-step, the
# opt-in `predict_method = "exact"`, and the opt-in `max_iter > 0` EM
# refine, never a `lambda_cols` decision.
fit_mvn_bm_inhouse <- function(L, tree = NULL, R = NULL,
                                max_iter = 0L, tol = 1e-4, eps = 1e-8,
                                use_henderson = TRUE,
                                sigma_method = c("single_pass", "fisher_ml"),
                                refine_variance = c("conservative",
                                                    "pooled"),
                                predict_method = c("per_column", "exact"),
                                lambda = "fixed_1",
                                lambda_cols = NULL) {
  sigma_method <- match.arg(sigma_method)
  refine_variance <- match.arg(refine_variance)
  predict_method <- match.arg(predict_method)
  if (is.null(R)) {
    if (is.null(tree)) stop("fit_mvn_bm_inhouse: either tree or R must be supplied.")
    R <- phylo_cor_matrix(tree)
  }
  spp <- rownames(L)
  if (is.null(spp)) stop("fit_mvn_bm_inhouse: L must have rownames.")
  R <- R[spp, spp]
  n <- nrow(L); K <- ncol(L)

  L_obs_mask <- !is.na(L)

  # ---- Pagel's lambda resolution ------------------------------------------
  lambda_cols_idx <- .mvn_resolve_lambda_cols(lambda_cols, K)
  lam <- .mvn_resolve_lambda(lambda, L, R, lambda_cols_idx, eps = eps)
  lambda_vec <- lam$lambda_vec
  lambda_block <- lam$lambda_block
  lambda_mode_used <- lam$lambda_mode_used
  mu_hat_vec <- lam$mu_hat_vec

  # Build Hadfield-Nakagawa (2010) sparse S^{-1} on demand, keyed by
  # lambda rounded to 1e-3 so at most K distinct builds happen (usually
  # far fewer -- most columns share lambda = 1 or lambda_block). The
  # build itself is O(n^2) in time and memory (build_henderson_S_inv
  # opens with `diag(ape::vcv(tree))`, which forms the dense n x n
  # covariance just to read its diagonal); it is the downstream SOLVE
  # that is O(n) via the sparse Q. Both are cheap relative to the
  # O(n^3) dense R^{-1} this replaces (>= 10x speed-up at n=2000).
  # Falls back to dense when tree is unset or Matrix package is
  # unavailable.
  #
  # K=1 stays on the legacy dense path (bm_impute_col, which estimates
  # a free GLS mean directly). For K>=2, henderson_bm_predict's own
  # zero-root assumption is corrected below by centering each column at
  # its own lambda-dependent GLS mean before the Henderson call and
  # adding it back after (mean-model consistency, section 7 B(i)).
  #
  # lambda_val == 1 reuses `tree` untransformed rather than routing
  # through transform_tree_pagel(tree, 1) -- both are mathematically
  # identical (internal edges * 1, terminal edges + 0 * depth), but the
  # literal reuse is what makes `lambda = "fixed_1"` identical to the
  # pre-lane solver rather than merely numerically close.
  henderson_cache <- new.env(parent = emptyenv())
  get_henderson_at <- function(lambda_val) {
    if (!isTRUE(use_henderson) || K < 2L || is.null(tree) ||
        !requireNamespace("Matrix", quietly = TRUE)) {
      return(NULL)
    }
    key <- sprintf("%.3f", round(lambda_val, 3))
    if (exists(key, envir = henderson_cache, inherits = FALSE)) {
      return(get(key, envir = henderson_cache, inherits = FALSE))
    }
    tr <- if (lambda_val == 1) tree else transform_tree_pagel(tree, lambda_val)
    h <- tryCatch(build_henderson_S_inv(tr), error = function(e) NULL)
    assign(key, h, envir = henderson_cache)
    h
  }

  # Shared Henderson object at lambda_block: used by the Sigma M-step,
  # the opt-in exact conditional, and the opt-in EM refine -- all three
  # need one common R(lambda) per the hybrid design (section 4 of the
  # dev-log). At lambda = "fixed_1" this is the same build as before.
  henderson_bar <- get_henderson_at(lambda_block)

  # ---- L_hat init: per-column BM, each column at its own lambda_k --------
  init <- .mvn_init_per_column(L, R, eps = eps,
                                henderson_for_col = function(j) {
                                  get_henderson_at(lambda_vec[j])
                                },
                                lambda_vec = lambda_vec,
                                mu_hat_vec = mu_hat_vec)
  L_hat <- init$L_hat
  L_var <- init$L_var

  # ---- Sigma init: closed-form Kronecker MLE on per-column-BM L_hat ------
  # Treating species as iid for the Sigma MLE (pairwise sample covariance,
  # complete-rows covariance, or the optim-based marginal Fisher-ML
  # `.mvn_sigma_fisher_ml`) badly underestimates Sigma when phylogenetic
  # signal is strong, because correlated species are double-counted as
  # independent observations. On the K=5 BM-on-coalescent-tree synthetic,
  # pairwise sample cov gives diag = [0.49, 0.07, 0.07, 2.63, 0.33] vs the
  # true diag = [1.47, 0.84, 0.39, 1.60, 1.06] -- two diagonals are 10x
  # too small. That tiny-diagonal Sigma is near-singular; build_conditional
  # _prior() then produces an ill-conditioned inverse that explodes L_hat
  # at missing cells in the first E-step, and Sigma blows up to ~30x truth.
  #
  # The closed-form matrix-normal MLE Sigma_hat = (1/n) L_hat^T R^{-1} L_hat
  # on per-column-BM-imputed L_hat properly credits R for the species
  # correlation, recovering diag = [1.09, 0.51, 0.24, 1.01, 0.90] -- the
  # right order of magnitude. Cross-trait off-diagonals are biased toward
  # zero in this init (per-column BM doesn't see cross-trait correlation),
  # but EM iterations below recover them via the cross-prior E-step.
  # `single_pass_fallback()` is exactly the Sigma this init block would
  # produce for K >= 2 without fisher_ml -- reused both as fisher_ml's
  # own optim() starting value and as its non-convergence fallback.
  single_pass_fallback <- function() {
    if (!is.null(henderson_bar)) {
      .mvn_sigma_kron_M(L_hat, L_var, henderson_bar)
    } else {
      S0 <- stats::cov(L, use = "pairwise.complete.obs")
      S0[!is.finite(S0)] <- 0
      diag(S0) <- pmax(diag(S0), eps)
      .mvn_ensure_pd(S0, eps = eps)
    }
  }

  Sigma <- if (K == 1L) {
    Sig <- stats::var(L[, 1L], na.rm = TRUE)
    if (!is.finite(Sig) || Sig <= 0) Sig <- 1
    matrix(Sig, 1L, 1L)
  } else if (identical(sigma_method, "fisher_ml")) {
    # Prototype step 1: pairwise-complete sample covariance seeds the
    # Cholesky optim (captures cross-trait sign/structure from observed
    # data before any per-column-BM imputation pollution). Step 2: the
    # optim() itself, in .mvn_sigma_fisher_ml() above.
    Sigma0 <- stats::cov(L, use = "pairwise.complete.obs")
    Sigma0[!is.finite(Sigma0)] <- 0
    diag(Sigma0) <- pmax(diag(Sigma0), eps)
    .mvn_sigma_fisher_ml(L, Sigma0, single_pass_fallback, eps = eps)
  } else if (!is.null(henderson_bar)) {
    .mvn_sigma_kron_M(L_hat, L_var, henderson_bar)
  } else {
    # No tree: closed-form M-step needs dense R^{-1}. Fall back to
    # pairwise-complete (PD-stabilised) and accept the species-iid bias.
    S0 <- stats::cov(L, use = "pairwise.complete.obs")
    S0[!is.finite(S0)] <- 0
    diag(S0) <- pmax(diag(S0), eps)
    .mvn_ensure_pd(S0, eps = eps)
  }
  Sigma <- .mvn_ensure_pd(Sigma, eps = eps)

  # ---- exact matrix-normal conditional (opt-in) --------------------------
  # predict_method = "exact" replaces the PREDICTION with the exact
  # conditional mean/variance under vec(L) ~ MVN(0, Sigma %x% R), computed
  # in the precision form via Hadfield & Nakagawa (2010) sparse S^{-1}.
  # Sigma estimation is unchanged (the 2026-08-17 recovery sim showed the
  # existing estimate is sound; it was the prediction step that discarded
  # it). Returns NULL on any failure -- oversized problem, singular Sigma,
  # non-finite output -- and we fall through to the per-column path rather
  # than degrade silently.
  # Derivation + gates: docs/dev-log/2026-08-17-exact-conditional-design.md
  if (identical(predict_method, "exact") && K >= 2L && !is.null(henderson_bar)) {
    ec <- exact_conditional_mvn(L, Sigma, henderson_bar, cor_scale = TRUE,
                                 eps = eps)
    if (!is.null(ec)) {
      return(list(
        anc_recon = ec$mu,
        anc_var   = ec$var,
        diverged  = FALSE,
        pars      = list(phylocov = Sigma),
        n_iter    = 0L,
        converged = TRUE,
        predict_method = "exact",
        lambda_per_trait = lambda_vec,
        lambda_block = lambda_block,
        lambda_mode_used = lambda_mode_used
      ))
    }
    warning("fit_mvn_bm_inhouse: predict_method = \"exact\" was not usable ",
            "here (problem too large, or a numerically unusable Sigma); ",
            "falling back to the per-column path.", call. = FALSE)
  }

  if (K == 1L || max_iter <= 0L) {
    # K=1: per-column BM is exact ML; nothing to refine.
    # max_iter=0: caller has opted out of the cross-trait EM refinement
    # (introduced 2026-05-17 after the EM was found to diverge on data
    # with strong phylogenetic signal at near-sister tips, because the
    # cross-trait conditional posterior in `build_conditional_prior()`
    # ignores R-mediated cross-row correlation; the resulting L_hat at
    # missing cells differs sharply from observed sister-tip values,
    # and R^{-1} amplifies that discrepancy in the next M-step into a
    # multiplicative Sigma blow-up). With max_iter = 0, returns the
    # per-column BM imputation and the closed-form L_hat^T R^{-1} L_hat
    # / n Sigma estimate -- both consistent under matrix-normal BM, and
    # empirically much better than the EM-refined version on synthetic
    # K=5 BM data (0.93 vs 0.53 argmax accuracy at n=100, miss=0.30).
    return(list(
      anc_recon = L_hat,
      anc_var   = L_var,
      pars      = list(phylocov = Sigma),
      n_iter    = 0L,
      converged = TRUE,
      lambda_per_trait = lambda_vec,
      lambda_block = lambda_block,
      lambda_mode_used = lambda_mode_used
    ))
  }

  converged <- FALSE
  diverged <- FALSE
  iter <- 0L
  # Divergence guard (2026-08-17). The EM was disabled outright on
  # 2026-05-17 because Sigma could blow up multiplicatively on
  # strong-signal data with near-sister tips. Rather than forbid the loop,
  # track the Sigma step size: it must SHRINK. If an iteration's relative
  # change exceeds the previous one, the loop is running away -- stop and
  # return the last good iterate rather than the diverged one. `prev_delta`
  # starts at Inf so the first iteration can never trip the guard.
  prev_delta <- Inf
  L_hat_prev <- L_hat; L_var_prev <- L_var; Sigma_prev <- Sigma
  for (k in seq_len(max(max_iter, 1L))) {
    iter <- k
    refined <- .mvn_estep_refine(L_obs_mask, L_hat, L_var, Sigma, eps = eps,
                                  refine_variance = refine_variance)
    # single_pass: proper Kronecker M-step when Henderson sparse R^{-1}
    # is available; fall back to closed-form sample-cov (drops R)
    # otherwise. The Kronecker M-step is the closed-form ML for
    # matrix-normal under EM and uses the phylogenetic R correctly.
    # fisher_ml: prototype step 4, M-step re-run of the observed-data
    # NLL optim on the *refined* (fully-completed) L_hat, seeded from
    # the previous iteration's Sigma -- matches e7ca41c's EM loop, which
    # passes `refined$L_hat` (not the original NA-carrying L) into the
    # optim at every M-step.
    loop_fallback <- function() {
      if (!is.null(henderson_bar)) {
        .mvn_sigma_kron_M(refined$L_hat, refined$L_var, henderson_bar)
      } else {
        .mvn_sigma_ml(refined$L_hat, refined$L_var,
                       R_inv = solve(R + diag(eps, n)))
      }
    }
    Sigma_new <- if (identical(sigma_method, "fisher_ml")) {
      .mvn_sigma_fisher_ml(refined$L_hat, Sigma, loop_fallback, eps = eps)
    } else if (!is.null(henderson_bar)) {
      .mvn_sigma_kron_M(refined$L_hat, refined$L_var, henderson_bar)
    } else {
      .mvn_sigma_ml(refined$L_hat, refined$L_var,
                     R_inv = solve(R + diag(eps, n)))
    }
    delta <- norm(Sigma_new - Sigma, "F") / max(norm(Sigma, "F"), eps)
    if (!is.finite(delta) || (k > 1L && delta > prev_delta)) {
      # Running away (or non-finite): discard this iterate, keep the last
      # good one. This is the guard that makes max_iter > 0 safe to expose.
      diverged <- TRUE
      L_hat <- L_hat_prev; L_var <- L_var_prev; Sigma <- Sigma_prev
      iter <- k - 1L
      break
    }
    prev_delta <- delta
    L_hat_prev <- L_hat; L_var_prev <- L_var; Sigma_prev <- Sigma
    L_hat <- refined$L_hat
    L_var <- refined$L_var
    Sigma <- Sigma_new
    if (delta < tol) { converged <- TRUE; break }
  }

  list(
    anc_recon = L_hat,
    anc_var   = L_var,
    diverged  = diverged,
    pars      = list(phylocov = Sigma),
    n_iter    = iter,
    converged = converged,
    lambda_per_trait = lambda_vec,
    lambda_block = lambda_block,
    lambda_mode_used = lambda_mode_used
  )
}

# Back-compat shim for .mvn_sigma_ml (no longer called directly by
# fit_mvn_bm_inhouse, but retained so any external test that imports it
# keeps working).
.mvn_sigma_ml <- function(L_hat, L_var, R_inv) {
  n <- nrow(L_hat)
  K <- ncol(L_hat)
  S <- crossprod(L_hat, R_inv %*% L_hat) + diag(colSums(L_var), nrow = K)
  S / n
}

# Run Rphylopars::phylopars() on the liability matrix L and return the
# fields pigauto consumes, in the same shape as fit_mvn_bm_inhouse():
# $anc_recon (tip + internal rows, K cols), $anc_var (variance, not SE),
# $pars$phylocov (K x K BM rate matrix).
#
# `model = "lambda"` lets phylopars estimate its own (single, shared)
# Pagel's lambda; this is NOT read back into $lambda_per_trait /
# $lambda_block in this slice (fit_joint_solver() sets both to
# NA_real_ on the rphylopars path) -- see that dispatcher's roxygen.
.fit_mvn_bm_rphylopars <- function(L, tree, model = "BM") {
  spp <- rownames(L)
  df <- data.frame(species = spp, L, stringsAsFactors = FALSE)
  Rphylopars::phylopars(df, tree = tree, model = model,
                        phylo_correlated = TRUE, pheno_correlated = TRUE,
                        REML = TRUE)
}

#' Dispatch the joint MVN Sigma solver
#'
#' Diagnosed in \code{docs/dev-log/2026-08-16-continuous-gap-diagnosis.md}:
#' on AVONET300 the in-house single-pass solver
#' (\code{fit_mvn_bm_inhouse()}, \code{max_iter = 0}) loses 0.14-1.27
#' z-RMSE to \code{Rphylopars::phylopars()}'s converged REML fit on
#' comparable data. This dispatcher lets callers opt in to the phylopars
#' solver while leaving the default (and byte-identical output) on the
#' in-house path.
#'
#' \code{joint_solver = "rphylopars"} is wrapped in a robustness guard:
#' if the phylopars call errors, or if any tip row of \code{$anc_recon}
#' is non-finite, this falls back to the in-house solver with a warning.
#' Rphylopars' own \code{solve(): system is singular} warnings are left
#' to propagate (an approximate-but-finite solution is not itself a
#' failure); only a hard error or non-finite output triggers fallback.
#'
#' @param L n x K liability matrix (rownames = tip labels, NAs allowed).
#' @param tree phylo.
#' @param joint_solver character, \code{"inhouse"} (default) or
#'   \code{"rphylopars"}.
#' @param sigma_method character, \code{"single_pass"} (default) or
#'   \code{"fisher_ml"}. Only consulted when \code{joint_solver =
#'   "inhouse"}; see \code{fit_mvn_bm_inhouse()}'s \code{sigma_method}
#'   argument for the algorithm and fallback contract.
#' @param joint_refine_iter integer, default \code{0L}. Forwarded to
#'   \code{fit_mvn_bm_inhouse()}'s \code{max_iter} argument (cross-trait
#'   EM cell-refinement using the estimated Sigma). \code{0L} preserves
#'   the single-pass, byte-identical default. Only consulted when
#'   \code{joint_solver = "inhouse"} (including the fallback path when
#'   \code{joint_solver = "rphylopars"} fails).
#' @param lambda \code{"fixed_1"} (default), \code{"estimate"}, or a
#'   numeric scalar in [0, 1]. Forwarded to \code{fit_mvn_bm_inhouse()}
#'   when \code{joint_solver = "inhouse"} (including its fallback path).
#'   When \code{joint_solver = "rphylopars"}, only whether \code{lambda}
#'   is \code{"fixed_1"} matters: it selects \code{model = "BM"} vs
#'   \code{model = "lambda"} in \code{.fit_mvn_bm_rphylopars()}.
#'   Phylopars' own lambda estimate is not read back in this slice --
#'   \code{$lambda_per_trait} and \code{$lambda_block} are set to
#'   \code{NA_real_} on that path.
#' @param lambda_cols integer/logical index into columns of \code{L}, or
#'   NULL (default, all columns). Forwarded to \code{fit_mvn_bm_inhouse()}
#'   when \code{joint_solver = "inhouse"}; see its \code{lambda_cols}
#'   argument. Not consulted on the \code{"rphylopars"} path.
#' @return list with the phylopars-compatible fields described in the
#'   file header comment above (\code{$anc_recon}, \code{$anc_var},
#'   \code{$pars$phylocov}), plus \code{$lambda_per_trait} and
#'   \code{$lambda_block}.
#' @keywords internal
#' @noRd
fit_joint_solver <- function(L, tree, joint_solver = "inhouse",
                             predict_method = "per_column",
                              sigma_method = "single_pass",
                              joint_refine_iter = 0L,
                              lambda = "fixed_1",
                              lambda_cols = NULL) {
  if (identical(joint_solver, "inhouse")) {
    return(fit_mvn_bm_inhouse(L = L, tree = tree, sigma_method = sigma_method,
                              predict_method = predict_method,
                               max_iter = joint_refine_iter,
                               lambda = lambda, lambda_cols = lambda_cols))
  }

  model <- if (identical(lambda, "fixed_1")) "BM" else "lambda"
  fit <- tryCatch(.fit_mvn_bm_rphylopars(L, tree, model = model),
                   error = function(e) e)
  ok <- !inherits(fit, "error")
  msg <- "non-finite tip prediction in $anc_recon"
  if (ok) {
    spp <- rownames(L)
    tip_rows <- match(spp, rownames(fit$anc_recon))
    ok <- !anyNA(tip_rows) && all(is.finite(fit$anc_recon[tip_rows, , drop = FALSE]))
    # Plausibility guard (feat/joint-lambda-default, 2026-09-23). Under
    # model = "lambda", phylopars occasionally returns finite but explosive
    # tip predictions: in the 18-cell simulation, 1 of 200 seeds at lambda 0.7,
    # n 100 and 6-8 of ~198 at lambda 1, n 1000 reached z-RMSE of 10^3 to 10^5,
    # while the in-house solver never exceeded 1.3. A prediction more than 10x
    # beyond the observed range of its column is not a usable imputation.
    if (ok) {
      obs_max <- max(abs(L), na.rm = TRUE)
      pred_max <- max(abs(fit$anc_recon[tip_rows, , drop = FALSE]))
      if (is.finite(obs_max) && pred_max > 10 * max(obs_max, 1)) {
        ok <- FALSE
        msg <- sprintf("implausible tip prediction (max |pred| %.3g vs max |observed| %.3g)",
                       pred_max, obs_max)
      }
    }
  }
  if (!ok) {
    if (inherits(fit, "error")) msg <- conditionMessage(fit)
    warning("fit_joint_solver: joint_solver = \"rphylopars\" failed (",
            msg, "); falling back to the in-house solver.", call. = FALSE)
    return(fit_mvn_bm_inhouse(L = L, tree = tree, max_iter = joint_refine_iter,
                              predict_method = predict_method,
                              lambda = lambda, lambda_cols = lambda_cols))
  }
  fit$lambda_per_trait <- NA_real_
  fit$lambda_block <- NA_real_
  fit
}
