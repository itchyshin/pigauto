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
#                rownames match input rownames). For observed cells
#                this is the input value (BM is exact on observed
#                data). For NA cells this is the EM estimate.
#   $anc_var   : n x K matrix of posterior cell variances. 0 at
#                observed cells, > 0 at imputed cells.
#   $pars$phylocov : K x K Sigma estimate.
#
# Algorithm (closed-form ML with EM for missing cells):
#
#   1. Per-column initialise L_hat, L_var via bm_impute_col() — exact
#      univariate-BM conditional MVN, using R only (no cross-trait info).
#   2. Matrix-normal MLE of Sigma from completed L_hat with variance
#      correction for the imputed cells:
#         Sigma_hat = (L_hat^T R^{-1} L_hat + diag(colSums(L_var))) / n
#   3. E-step: for each originally-NA cell, refine the imputation by
#      combining the per-column BM posterior with the cross-trait
#      conditional MVN posterior (build_conditional_prior). Both
#      contribute via inverse-variance pooling.
#   4. M-step: re-estimate Sigma from the refined L_hat.
#   5. Iterate 3-4 until rel-Frobenius change in Sigma < tol.
#
# No external dependencies. ~150 LOC including docs.

# Closed-form matrix-normal MLE of Sigma given completed L and R.
#
# Sigma_hat = (1/n) * (L^T R^{-1} L + V_correction)
# where V_correction = diag(colSums of cell-level posterior variances).
# When L is fully observed and V_correction = 0 this is the classical
# MLE; when cells are imputed with variance V_ij the correction term
# is the trace of the E[L^T R^{-1} L] expectation under the posterior.
.mvn_sigma_ml <- function(L_hat, L_var, R_inv) {
  n <- nrow(L_hat)
  # `diag(numeric_scalar)` is ambiguous in R (treated as dim, not value), so
  # be explicit about nrow = K to handle the K=1 single-trait case correctly.
  K <- ncol(L_hat)
  S <- crossprod(L_hat, R_inv %*% L_hat) + diag(colSums(L_var), nrow = K)
  S / n
}

# Per-column BM init: exact univariate-BM conditional MVN. Uses R only.
# Returns L_hat (NAs replaced) and L_var (0 at observed cells, > 0 at
# imputed cells). For columns with < 2 observed values, leaves the
# column at zero with unit variance (uninformative prior).
.mvn_init_per_column <- function(L, R, eps = 1e-8) {
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
    res <- bm_impute_col(yj, R, nugget = eps)
    L_hat[, j] <- res$mu
    L_var[, j] <- res$se^2
    # Observed cells are exact for BM: no posterior variance
    L_var[obs, j] <- 0
  }
  list(L_hat = L_hat, L_var = L_var)
}

# E-step: refine imputations at originally-NA cells via inverse-variance
# pooling between the per-column BM posterior and the cross-trait
# conditional MVN posterior (which uses Sigma off-diagonals).
.mvn_estep_refine <- function(L_obs_mask, L_hat, L_var, Sigma, eps = 1e-8) {
  # Cross-trait conditional prior at every cell (treats species iid,
  # conditions cell (i, j) on observed L_hat[i, -j]).
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
#   max_iter : maximum EM iterations (default 5).
#   tol   : relative-Frobenius stopping criterion on Sigma.
#   eps   : ridge added to covariance solves.
#
# Returns a list with the phylopars-compatible fields described above,
# plus diagnostics ($n_iter, $converged).
fit_mvn_bm_inhouse <- function(L, tree = NULL, R = NULL,
                                max_iter = 5L, tol = 1e-3, eps = 1e-8) {
  if (is.null(R)) {
    if (is.null(tree)) stop("fit_mvn_bm_inhouse: either tree or R must be supplied.")
    R <- phylo_cor_matrix(tree)
  }
  spp <- rownames(L)
  if (is.null(spp)) stop("fit_mvn_bm_inhouse: L must have rownames.")
  R <- R[spp, spp]
  n <- nrow(L); K <- ncol(L)

  L_obs_mask <- !is.na(L)
  init <- .mvn_init_per_column(L, R, eps = eps)
  L_hat <- init$L_hat
  L_var <- init$L_var

  R_reg <- R + diag(eps, n)
  R_chol <- tryCatch(chol(R_reg), error = function(e) NULL)
  if (is.null(R_chol)) {
    # Tiny additional ridge for nearly-singular trees.
    R_reg  <- R + diag(eps * 10, n)
    R_chol <- chol(R_reg)
  }
  R_inv <- chol2inv(R_chol)

  Sigma <- .mvn_sigma_ml(L_hat, L_var, R_inv)

  # K = 1 has no off-diagonals to exploit. Return per-column BM result
  # directly so the single-trait case is bit-identical to the legacy
  # per-column BM path that callers and tests rely on.
  if (K == 1L) {
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
  for (k in seq_len(max(max_iter - 1L, 0L))) {
    iter <- k
    refined <- .mvn_estep_refine(L_obs_mask, L_hat, L_var, Sigma, eps = eps)
    Sigma_new <- .mvn_sigma_ml(refined$L_hat, refined$L_var, R_inv)
    delta <- norm(Sigma_new - Sigma, "F") /
              max(norm(Sigma, "F"), eps)
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
