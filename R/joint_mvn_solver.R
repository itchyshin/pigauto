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

# Per-column BM init: exact univariate-BM conditional MVN.
# Returns L_hat (NAs replaced) and L_var (0 at observed cells, > 0 at
# imputed cells). For columns with < 2 observed values, leaves the
# column at zero with unit variance (uninformative prior).
#
# Two code paths share the same numerical answer to ~1e-3:
#   - bm_impute_col(yj, R)         legacy dense O(n_obs^3) per column.
#   - henderson_bm_predict(yj, H)  sparse O(n) per column via Hadfield-
#     Nakagawa (2010) eq 29. Cuts init time at large n from minutes to
#     seconds (~18x at n=2000) and avoids forming dense (n x n) R^{-1}.
.mvn_init_per_column <- function(L, R, eps = 1e-8, henderson = NULL) {
  n <- nrow(L); K <- ncol(L)
  L_hat <- L; L_var <- matrix(0, n, K, dimnames = dimnames(L))
  for (j in seq_len(K)) {
    yj <- L[, j]
    obs <- !is.na(yj)
    if (sum(obs) < 2L) {
      L_hat[, j] <- ifelse(is.na(yj), 0, yj)
      L_var[!obs, j] <- 1
      next
    }
    res <- if (!is.null(henderson)) {
      henderson_bm_predict(yj, henderson, eps = eps, cor_scale = TRUE)
    } else {
      bm_impute_col(yj, R, nugget = eps)
    }
    L_hat[, j] <- res$mu
    L_var[, j] <- res$se^2
    L_var[obs, j] <- 0
  }
  list(L_hat = L_hat, L_var = L_var)
}

# E-step: refine imputations at originally-NA cells via inverse-variance
# pooling between the per-column BM posterior and the cross-trait
# conditional MVN posterior (uses Sigma off-diagonals via
# build_conditional_prior).
.mvn_estep_refine <- function(L_obs_mask, L_hat, L_var, Sigma, eps = 1e-8) {
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
    L_var_new[idx, j] <- 1 / prec_tot
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
#
# Returns a list with the phylopars-compatible fields described above,
# plus diagnostics ($n_iter, $converged).
fit_mvn_bm_inhouse <- function(L, tree = NULL, R = NULL,
                                max_iter = 0L, tol = 1e-4, eps = 1e-8,
                                use_henderson = TRUE,
                                sigma_method = c("single_pass", "fisher_ml")) {
  sigma_method <- match.arg(sigma_method)
  if (is.null(R)) {
    if (is.null(tree)) stop("fit_mvn_bm_inhouse: either tree or R must be supplied.")
    R <- phylo_cor_matrix(tree)
  }
  spp <- rownames(L)
  if (is.null(spp)) stop("fit_mvn_bm_inhouse: L must have rownames.")
  R <- R[spp, spp]
  n <- nrow(L); K <- ncol(L)

  L_obs_mask <- !is.na(L)

  # Build Hadfield-Nakagawa (2010) sparse S^{-1} once; reuse across the
  # K per-column BM imputations. O(n) build time vs O(n^3) for dense
  # R^{-1}; >= 10x speed-up at n=2000 with identical mean predictions
  # to ~1e-3 vs the dense path. Falls back to dense when tree is unset
  # or Matrix package is unavailable.
  #
  # K=1 stays on the legacy dense path: bm_impute_col estimates a GLS
  # phylogenetic mean from data, while henderson_bm_predict assumes a
  # zero root state. For joint (K>=2) liability matrices the inputs
  # are centered by build_liability_matrix so the zero-root assumption
  # is correct; for single-trait back-compat the GLS-mean estimate
  # matters.
  henderson <- if (isTRUE(use_henderson) && K >= 2L && !is.null(tree) &&
                    requireNamespace("Matrix", quietly = TRUE)) {
    tryCatch(build_henderson_S_inv(tree), error = function(e) NULL)
  } else NULL

  # ---- L_hat init: per-column BM (uses phylo R correctly per column) -----
  init <- .mvn_init_per_column(L, R, eps = eps, henderson = henderson)
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
    if (!is.null(henderson)) {
      .mvn_sigma_kron_M(L_hat, L_var, henderson)
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
  } else if (!is.null(henderson)) {
    .mvn_sigma_kron_M(L_hat, L_var, henderson)
  } else {
    # No tree: closed-form M-step needs dense R^{-1}. Fall back to
    # pairwise-complete (PD-stabilised) and accept the species-iid bias.
    S0 <- stats::cov(L, use = "pairwise.complete.obs")
    S0[!is.finite(S0)] <- 0
    diag(S0) <- pmax(diag(S0), eps)
    .mvn_ensure_pd(S0, eps = eps)
  }
  Sigma <- .mvn_ensure_pd(Sigma, eps = eps)

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
      converged = TRUE
    ))
  }

  converged <- FALSE
  iter <- 0L
  for (k in seq_len(max(max_iter, 1L))) {
    iter <- k
    refined <- .mvn_estep_refine(L_obs_mask, L_hat, L_var, Sigma, eps = eps)
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
      if (!is.null(henderson)) {
        .mvn_sigma_kron_M(refined$L_hat, refined$L_var, henderson)
      } else {
        .mvn_sigma_ml(refined$L_hat, refined$L_var,
                       R_inv = solve(R + diag(eps, n)))
      }
    }
    Sigma_new <- if (identical(sigma_method, "fisher_ml")) {
      .mvn_sigma_fisher_ml(refined$L_hat, Sigma, loop_fallback, eps = eps)
    } else if (!is.null(henderson)) {
      .mvn_sigma_kron_M(refined$L_hat, refined$L_var, henderson)
    } else {
      .mvn_sigma_ml(refined$L_hat, refined$L_var,
                     R_inv = solve(R + diag(eps, n)))
    }
    delta <- norm(Sigma_new - Sigma, "F") / max(norm(Sigma, "F"), eps)
    L_hat <- refined$L_hat
    L_var <- refined$L_var
    Sigma <- Sigma_new
    if (delta < tol) { converged <- TRUE; break }
  }

  list(
    anc_recon = L_hat,
    anc_var   = L_var,
    pars      = list(phylocov = Sigma),
    n_iter    = iter,
    converged = converged
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
.fit_mvn_bm_rphylopars <- function(L, tree) {
  spp <- rownames(L)
  df <- data.frame(species = spp, L, stringsAsFactors = FALSE)
  Rphylopars::phylopars(df, tree = tree, model = "BM",
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
#' @return list with the phylopars-compatible fields described in the
#'   file header comment above (\code{$anc_recon}, \code{$anc_var},
#'   \code{$pars$phylocov}).
#' @keywords internal
#' @noRd
fit_joint_solver <- function(L, tree, joint_solver = "inhouse",
                              sigma_method = "single_pass") {
  if (identical(joint_solver, "inhouse")) {
    return(fit_mvn_bm_inhouse(L = L, tree = tree, sigma_method = sigma_method))
  }

  fit <- tryCatch(.fit_mvn_bm_rphylopars(L, tree), error = function(e) e)
  ok <- !inherits(fit, "error")
  if (ok) {
    spp <- rownames(L)
    tip_rows <- match(spp, rownames(fit$anc_recon))
    ok <- !anyNA(tip_rows) && all(is.finite(fit$anc_recon[tip_rows, , drop = FALSE]))
  }
  if (!ok) {
    msg <- if (inherits(fit, "error")) conditionMessage(fit) else
      "non-finite tip prediction in $anc_recon"
    warning("fit_joint_solver: joint_solver = \"rphylopars\" failed (",
            msg, "); falling back to the in-house solver.", call. = FALSE)
    return(fit_mvn_bm_inhouse(L = L, tree = tree))
  }
  fit
}
