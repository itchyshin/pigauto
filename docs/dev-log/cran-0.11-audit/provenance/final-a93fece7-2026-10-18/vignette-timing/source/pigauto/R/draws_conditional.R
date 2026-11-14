# R/draws_conditional.R
#
# Prototype: JOINT proper draws from the multivariate-BM baseline posterior,
# for the mi-gls-attenuation investigation (arc/mi-gls-attenuation).
#
# Background (measured this session; see useful/mondrian-mi-se-justification.md
# "## Result" and script/mondrian_confirmation/13_mi_gls_attenuation_diag.R):
# with bivariate BM (rho = 0.7), x partly missing, y observed,
# multi_impute(draws_method = "conformal") halves the pooled PGLS slope of
# y ~ x (0.35 vs the true 0.70). Two compounding causes:
#
#   (a) CENTRE: the default predict_method = "per_column" baseline point
#       prediction is (near-)the x-only conditional mean, not the joint
#       (x, y) conditional mean -- residual SD 0.325 vs the joint-oracle
#       0.239.
#   (b) SPREAD: multi_impute()'s conformal-width draws sample each missing
#       cell independently (R/multi_impute.R, .sample_conformal_draw()) --
#       conformal draw SD 0.293 vs the proper joint-conditional SD 0.218,
#       AND draws are independent across tips, so downstream GLS sees no
#       phylogenetically-correlated imputation uncertainty at all.
#
# draw_conditional_bm() fixes both by sampling the missing latent cells
# JOINTLY (across tips AND traits at once) from the EXACT Gaussian
# conditional of vec(L) ~ MVN(0, Sigma %x% R) given all observed cells,
# where Sigma is pigauto's own joint-MVN trait-covariance estimate
# (fit_mvn_bm_inhouse(), R/joint_mvn_solver.R -- read-only reuse, not
# re-estimated here) and R = cov2cor(vcv(tree)) via the sparse
# Henderson-Nakagawa precision (R/henderson_s_inv.R).
#
# This is explicitly the r_cal = 0 (baseline-only) proper draw: no GNN
# correction is applied. pigauto's r_cal = 0 fallback is deliberately
# preserved as a safety net elsewhere in the package (see fit_pigauto.R /
# AGENTS.md); this prototype is that same fallback used as the reference
# "proper MI" baseline for the attenuation investigation, NOT a general
# replacement for multi_impute(). It is NOT wired into multi_impute().
#
# The mean computation mirrors exact_conditional_mvn() (R/exact_conditional.R)
# exactly (same index conventions, same sparse Kronecker precision), and is
# extended here to draw a full JOINT sample (not just report the marginal
# per-cell variance) via Matheron's conditional-simulation identity on the
# sparse Cholesky factor: for precision P_uu = L L^T (perm = FALSE, so no
# permutation bookkeeping is needed), x = L^{-T} w with w ~ N(0, I) satisfies
# Cov(x) = L^{-T} L^{-1} = P_uu^{-1} exactly, so mu_cond + x is a proper draw
# from N(mu_cond, P_uu^{-1}).

# ---- Internal: pull $data / $tree out of a pigauto_result-like object -----
.dcb_extract_data_tree <- function(fit_or_result) {
  data <- fit_or_result$data
  tree <- fit_or_result$tree
  if (is.null(data) || !inherits(data, "pigauto_data") || is.null(tree)) {
    stop("draw_conditional_bm(): `fit_or_result` must expose $data (a ",
         "pigauto_data object from preprocess_traits()) and $tree -- e.g. ",
         "the object returned by impute(). Got an object of class ",
         paste(class(fit_or_result), collapse = "/"), ".", call. = FALSE)
  }
  list(data = data, tree = tree)
}

# ---- Internal: decode one z-scored latent matrix back to original scale ---
# Continuous-only decode: raw = z * sd + mean, then exp() if log_transform.
# Mirrors the continuous branch of .sample_conformal_draw() (R/multi_impute.R).
.dcb_decode_latent <- function(Xk, trait_map) {
  nms <- vapply(trait_map, `[[`, character(1), "name")
  out <- as.data.frame(matrix(NA_real_, nrow(Xk), length(trait_map)))
  names(out) <- nms
  rownames(out) <- rownames(Xk)
  for (tm in trait_map) {
    z   <- Xk[, tm$latent_cols]
    raw <- z * tm$sd + tm$mean
    if (isTRUE(tm$log_transform)) raw <- exp(raw)
    out[[tm$name]] <- raw
  }
  out
}

#' Joint conditional multivariate-BM draws (prototype, internal)
#'
#' Draws `m` proper completions of the missing continuous-trait cells,
#' jointly across tips and traits, from the exact Gaussian conditional of
#' the multivariate-BM baseline `vec(L) ~ MVN(0, Sigma \%x\% R)`. See the
#' file header of R/draws_conditional.R for the investigation background.
#'
#' Continuous traits only (errors clearly otherwise); single observation
#' per species only (errors on multi-obs data).
#'
#' @param fit_or_result An object exposing `$data` (a `pigauto_data` from
#'   [preprocess_traits()]) and `$tree` -- typically the return value of
#'   [impute()].
#' @param m integer, number of joint draws (default `20L`).
#' @param seed optional integer; when supplied, the whole batch of `m`
#'   draws is generated from one reproducible `set.seed(seed)` call (draws
#'   differ by column of the underlying random matrix, not by re-seeding
#'   per draw).
#' @param sigma_method `"em"` (default) estimates the trait covariance by
#'   EM maximum likelihood with missing cells; `"inhouse"` reuses the
#'   plug-in estimate from `fit_mvn_bm_inhouse(max_iter = 0)`, which shrinks
#'   cross-trait covariance when several traits have missing cells.
#' @param eps numeric ridge added to `Sigma` and to edge lengths before
#'   inversion (default `1e-8`).
#' @param max_cells integer safety cap on the number of unknown
#'   (trait, extended-tree-node) cells; refuses rather than hangs on an
#'   oversized problem (default `60000L`).
#'
#' @return A list of class `"pigauto_draws_conditional_bm"`:
#'   \describe{
#'     \item{`datasets`}{list of `m` data.frames, original scale, same
#'       columns/rownames as `data$X_original`; observed cells preserved,
#'       missing cells filled with the joint draw.}
#'     \item{`latent_draws`}{list of `m` `n_species x K` z-scored latent
#'       matrices (pre-decode); useful for testing/diagnostics against the
#'       exact Gaussian conditional directly.}
#'     \item{`mu_cond`}{`n_species x K` conditional mean (z-scored latent;
#'       the r_cal = 0 baseline-only point prediction); observed cells
#'       echoed.}
#'     \item{`sigma_hat`}{the `K x K` trait covariance estimate used
#'       (pigauto's own `fit_mvn_bm_inhouse()` estimate).}
#'     \item{`imputed_mask`}{`n_species x K` logical, `TRUE` where a cell
#'       was originally missing.}
#'     \item{`m`}{number of draws.}
#'   }
#' @keywords internal
#' @noRd
draw_conditional_bm <- function(fit_or_result, m = 20L, seed = NULL,
                                 eps = 1e-8, max_cells = 60000L,
                                 sigma_method = c("em", "inhouse")) {
  sigma_method <- match.arg(sigma_method)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("draw_conditional_bm() requires the 'Matrix' package.", call. = FALSE)
  }
  m <- as.integer(m)
  if (!is.finite(m) || m < 1L) {
    stop("`m` must be a positive integer. Got m = ", m, ".", call. = FALSE)
  }

  dt   <- .dcb_extract_data_tree(fit_or_result)
  data <- dt$data
  tree <- dt$tree

  if (isTRUE(data$multi_obs)) {
    stop("draw_conditional_bm() supports single-observation-per-species ",
         "data only (multi_obs = FALSE). This prototype does not handle ",
         "the multi-obs aggregation path.", call. = FALSE)
  }

  trait_map  <- data$trait_map
  all_types  <- vapply(trait_map, `[[`, character(1), "type")
  bad_types  <- setdiff(unique(all_types), "continuous")
  if (length(bad_types) > 0L) {
    stop("draw_conditional_bm() supports continuous traits only. Found ",
         "non-continuous trait type(s): ", paste(bad_types, collapse = ", "),
         ". This is a prototype for the mi-gls-attenuation investigation, ",
         "not a general-purpose replacement for multi_impute().",
         call. = FALSE)
  }

  X <- data$X_scaled
  if (is.null(rownames(X))) rownames(X) <- data$species_names
  K <- ncol(X)
  n <- nrow(X)

  henderson <- build_henderson_S_inv(tree)
  if (henderson$n_tips != n) {
    stop("draw_conditional_bm(): tree tip count does not match the ",
         "preprocessed data (internal error).", call. = FALSE)
  }

  # ---- Sigma: reuse pigauto's own joint-MVN estimator, read-only ---------
  # fit_mvn_bm_inhouse() (R/joint_mvn_solver.R) is exactly what
  # fit_joint_mvn_baseline() calls internally. max_iter = 0L / predict_method
  # = "per_column" is its cheapest mode: it still returns the closed-form
  # matrix-normal Sigma_hat = (1/n) L_hat^T R^{-1} L_hat -- a GLS/REML-style
  # estimate against R = cov2cor(vcv(tree)) -- without running the (slower,
  # occasionally divergent) cross-trait EM refinement. The mean/variance it
  # also returns is NOT used here; the joint (not just per-cell) conditional
  # is built separately below because it needs the Cholesky factor itself,
  # not just its diagonal.
  if (identical(sigma_method, "inhouse")) {
    sig_fit <- fit_mvn_bm_inhouse(L = X, tree = tree, R = NULL, max_iter = 0L,
                                   predict_method = "per_column")
    Sigma <- sig_fit$pars$phylocov
  } else {
    # EM ML estimate (see .dcb_sigma_em): the plug-in "inhouse" estimate
    # shrinks cross-trait covariance when several traits have missing cells.
    R_tip <- stats::cov2cor(ape::vcv(tree))[rownames(X), rownames(X)]
    Sigma <- .dcb_sigma_em(X, R_tip)
  }
  Sigma <- (Sigma + t(Sigma)) / 2

  Sig_inv <- tryCatch(solve(Sigma + diag(eps, K)), error = function(e) NULL)
  if (is.null(Sig_inv)) {
    stop("draw_conditional_bm(): Sigma is numerically singular.",
         call. = FALSE)
  }
  Sig_inv <- (Sig_inv + t(Sig_inv)) / 2

  N     <- henderson$N
  Q     <- henderson$Q
  T_idx <- henderson$tip_idx
  I_idx <- henderson$int_idx
  cell_idx <- function(t, i) (t - 1L) * N + i

  obs_cells <- integer(0); obs_vals <- numeric(0)
  unk_cells <- integer(0)
  miss_row  <- integer(0); miss_col <- integer(0)

  for (t in seq_len(K)) {
    y  <- X[, t]
    o  <- which(!is.na(y))
    mi <- which(is.na(y))
    if (length(o)) {
      yo <- y[o] * henderson$tip_sqrt_d[o]
      obs_cells <- c(obs_cells, cell_idx(t, T_idx[o]))
      obs_vals  <- c(obs_vals, yo)
    }
    if (length(mi)) {
      unk_cells <- c(unk_cells, cell_idx(t, T_idx[mi]))
      miss_row  <- c(miss_row, mi)
      miss_col  <- c(miss_col, rep(t, length(mi)))
    }
  }
  n_miss_tip <- length(unk_cells)

  mu_full <- X
  if (n_miss_tip == 0L) {
    latent_list <- replicate(m, X, simplify = FALSE)
    datasets <- lapply(latent_list, function(Xk) {
      decoded <- .dcb_decode_latent(Xk, trait_map)
      build_completed(data$X_original, decoded, species_col = NULL,
                      input_row_order = NULL)$completed
    })
    return(structure(
      list(datasets = datasets, latent_draws = latent_list, mu_cond = mu_full,
           sigma_hat = Sigma, imputed_mask = matrix(FALSE, n, K,
                                                     dimnames = dimnames(X)),
           m = m),
      class = "pigauto_draws_conditional_bm"
    ))
  }
  if (!length(obs_cells)) {
    stop("draw_conditional_bm(): no observed cells to condition on.",
         call. = FALSE)
  }

  for (t in seq_len(K)) unk_cells <- c(unk_cells, cell_idx(t, I_idx))
  n_unk <- length(unk_cells)
  if (n_unk > max_cells) {
    stop("draw_conditional_bm(): problem too large (", n_unk, " unknown ",
         "cells > max_cells = ", max_cells, "). Reduce K or n, or raise ",
         "max_cells.", call. = FALSE)
  }

  P <- Matrix::kronecker(Matrix::Matrix(Sig_inv, sparse = FALSE), Q)
  P_uu <- P[unk_cells, unk_cells, drop = FALSE]
  P_uo <- P[unk_cells, obs_cells, drop = FALSE]

  # perm = FALSE: skip the AMD/METIS fill-reducing permutation so the sparse
  # Cholesky factor can be used directly (system = "Lt") for Matheron's
  # conditional-simulation trick without a separate un-permute step. A
  # deliberate simplicity-over-speed trade-off for this prototype; fine at
  # the regime sizes here (n <= 1000, K = 2).
  chol_uu <- tryCatch(
    Matrix::Cholesky(Matrix::forceSymmetric(P_uu), LDL = FALSE, perm = FALSE),
    error = function(e) NULL)
  if (is.null(chol_uu)) {
    stop("draw_conditional_bm(): sparse Cholesky of the conditional ",
         "precision failed (Sigma or the tree may be degenerate).",
         call. = FALSE)
  }

  # Conditional mean -- identical formula to exact_conditional_mvn()'s mu.
  mu_v   <- as.numeric(Matrix::solve(chol_uu, -as.matrix(P_uo %*% obs_vals)))
  mu_tip <- mu_v[seq_len(n_miss_tip)]

  # m joint draws at once: W is n_unk x m, one column per draw, one
  # factorisation reused for all of them.
  if (!is.null(seed)) set.seed(as.integer(seed))
  W      <- matrix(stats::rnorm(n_unk * m), n_unk, m)
  Xdraw  <- as.matrix(Matrix::solve(chol_uu, W, system = "Lt"))
  draw_tip <- Xdraw[seq_len(n_miss_tip), , drop = FALSE]

  sd_r      <- henderson$tip_sqrt_d[miss_row]
  mu_orig   <- mu_tip / sd_r
  draw_orig <- draw_tip / sd_r   # recycled row-wise: nrow(draw_tip) == length(sd_r)

  if (!all(is.finite(mu_orig)) || !all(is.finite(draw_orig))) {
    stop("draw_conditional_bm(): non-finite draw produced (numerical ",
         "failure).", call. = FALSE)
  }

  mu_full[cbind(miss_row, miss_col)] <- mu_orig

  latent_list <- vector("list", m)
  for (k in seq_len(m)) {
    Xk <- X
    Xk[cbind(miss_row, miss_col)] <- mu_orig + draw_orig[, k]
    latent_list[[k]] <- Xk
  }

  datasets <- lapply(latent_list, function(Xk) {
    decoded <- .dcb_decode_latent(Xk, trait_map)
    build_completed(data$X_original, decoded, species_col = NULL,
                    input_row_order = NULL)$completed
  })

  imputed_mask <- matrix(FALSE, n, K, dimnames = dimnames(X))
  imputed_mask[cbind(miss_row, miss_col)] <- TRUE

  structure(
    list(datasets = datasets, latent_draws = latent_list, mu_cond = mu_full,
         sigma_hat = Sigma, imputed_mask = imputed_mask, m = m),
    class = "pigauto_draws_conditional_bm"
  )
}

# ---------------------------------------------------------------------------
# .dcb_sigma_em(L, R): EM maximum-likelihood estimate of the K x K trait
# covariance Sigma for vec(L) ~ MVN(0, Sigma %x% R) with missing cells.
#
# Why (2026-09-23, arc/mi-gls-attenuation): the plug-in estimate from
# fit_mvn_bm_inhouse(max_iter = 0) fills each column with its own per-column
# conditional mean, which carries no cross-trait information, and omits the
# conditional-covariance term. With 30% of both traits missing it shrank a
# true phylogenetic correlation of 0.67 to 0.46 (six trees, n = 300), and its
# EM option returned the same value. This estimator is the textbook EM for a
# matrix-normal model with missing entries:
#   E-step: mu_m = C_mo C_oo^{-1} l_o and V_m = C_mm - C_mo C_oo^{-1} C_om,
#           with C = Sigma %x% R;
#   M-step: Sigma_jk = (1/n) [ Lhat_j' R^{-1} Lhat_k
#                              + sum_{a,b} R^{-1}_ab V[(a,j), (b,k)] ].
# Dense linear algebra: intended for the prototype's sizes (n <= ~1000,
# K small). On complete data it returns the closed form in one step.
# ---------------------------------------------------------------------------
.dcb_sigma_em <- function(L, R, max_iter = 200L, tol = 1e-7) {
  L <- as.matrix(L); n <- nrow(L); K <- ncol(L)
  Rinv <- solve(R)
  miss <- is.na(L)
  Lf <- L
  for (j in seq_len(K)) Lf[miss[, j], j] <- 0
  Sigma <- crossprod(Lf, Rinv %*% Lf) / n
  if (!any(miss)) return((Sigma + t(Sigma)) / 2)
  diag(Sigma) <- pmax(diag(Sigma), 1e-6)
  v <- as.vector(L)
  m_idx <- which(is.na(v)); o_idx <- which(!is.na(v))
  tip_of <- ((seq_along(v) - 1L) %% n) + 1L
  trait_of <- ((seq_along(v) - 1L) %/% n) + 1L
  for (it in seq_len(max_iter)) {
    C <- kronecker(Sigma, R)
    C_oo <- C[o_idx, o_idx]; C_mo <- C[m_idx, o_idx]
    A <- t(solve(C_oo, t(C_mo)))                 # C_mo C_oo^{-1}
    mu_m <- as.vector(A %*% v[o_idx])
    V_m <- C[m_idx, m_idx] - A %*% t(C_mo)
    Lhat <- L; Lhat[miss] <- mu_m
    S_new <- crossprod(Lhat, Rinv %*% Lhat)
    ta <- tip_of[m_idx]; tr <- trait_of[m_idx]
    W <- Rinv[ta, ta] * V_m                      # R^{-1}_ab V[(a,j),(b,k)]
    for (j in seq_len(K)) for (k in seq_len(K)) {
      S_new[j, k] <- S_new[j, k] + sum(W[tr == j, tr == k, drop = FALSE])
    }
    S_new <- (S_new + t(S_new)) / (2 * n)
    if (max(abs(S_new - Sigma)) < tol * max(1, max(abs(Sigma)))) {
      Sigma <- S_new; break
    }
    Sigma <- S_new
  }
  Sigma
}
