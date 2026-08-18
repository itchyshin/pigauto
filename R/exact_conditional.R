# R/exact_conditional.R
#
# Exact matrix-normal conditional for the joint BM baseline.
#
# Model:  vec(L) ~ MVN(0, Sigma %x% R),  L is n x K, R = cov2cor(vcv(tree)).
#
# The correct baseline prediction is the conditional mean/variance of the
# missing cells given the observed cells. Doing that in the COVARIANCE form
# needs Sigma_OO^{-1}, which is |O| x |O| dense (|O| ~ 0.7 n K) -- cubic and
# hopeless at real n. In the PRECISION form it is cheap, because for a
# Gaussian with joint precision P:
#
#   precision(L_unk | L_obs) = P[unk, unk]
#   mean(L_unk | L_obs)      = -P[unk, unk]^{-1} P[unk, obs] L_obs
#
# with NO inversion of the observed block. And the joint precision here is
#
#   P = Sigma^{-1} %x% S^{-1}
#
# where S^{-1} is the Hadfield-Nakagawa sparse precision over the EXTENDED
# tree (tips + non-root internals) from build_henderson_S_inv(): O(n)
# nonzeros. So P is sparse and a sparse Cholesky of P[unk, unk] is the whole
# computation.
#
# This is the K-trait Kronecker lift of henderson_bm_predict()
# (R/henderson_s_inv.R), which already does exactly this for K = 1 and is in
# production. Conventions (cor_scale via tip_sqrt_d, zero root state,
# internals treated as unobserved) are inherited from it deliberately -- do
# not re-derive them here.
#
# Design + pre-registered gates:
#   docs/dev-log/2026-08-17-exact-conditional-design.md

# Cell index in the lifted space: trait t (1..K), extended node i (1..N)
# maps to (t - 1) * N + i. This matches kronecker(Sigma_inv, Q), whose
# (t,u) block is Sigma_inv[t,u] * Q.
.ec_cell <- function(t, i, N) (t - 1L) * N + i

# Solve for the diagonal of P_uu^{-1} at the requested positions, in
# blocks of `block` unit right-hand sides. Exact (not stochastic); the
# blocking only bounds peak memory. `pos` are positions within the `unk`
# ordering.
.ec_inv_diag <- function(chol_uu, n_unk, pos, block = 256L) {
  out <- numeric(length(pos))
  if (!length(pos)) return(out)
  starts <- seq(1L, length(pos), by = block)
  for (s in starts) {
    idx <- s:min(s + block - 1L, length(pos))
    rhs <- Matrix::sparseMatrix(i = pos[idx], j = seq_along(idx),
                                x = rep(1, length(idx)),
                                dims = c(n_unk, length(idx)))
    V <- as.matrix(Matrix::solve(chol_uu, rhs))
    out[idx] <- V[cbind(pos[idx], seq_along(idx))]
  }
  out
}

# Exact conditional mean/variance for the missing cells of L.
#
# @param L        n x K matrix, NA at missing cells, rownames = tip labels
#                 in the tree's tip order.
# @param Sigma    K x K trait covariance (the existing estimate; this
#                 function changes PREDICTION, not estimation).
# @param henderson  output of build_henderson_S_inv(tree).
# @param cor_scale  work on R = cov2cor(A) (TRUE) or raw A (FALSE).
#                   TRUE matches the joint solver's convention.
# @param max_cells  refuse (return NULL) above this many unknown cells, so
#                   the caller can fall back rather than hang. The
#                   variance step is the binding cost.
# @return list(mu, var) each n x K (observed cells echo L, var 0 there),
#         or NULL if the problem is too large / inputs unusable.
exact_conditional_mvn <- function(L, Sigma, henderson, cor_scale = TRUE,
                                   eps = 1e-8, max_cells = 20000L,
                                   block = 256L) {
  if (!requireNamespace("Matrix", quietly = TRUE)) return(NULL)
  n <- nrow(L); K <- ncol(L)
  if (is.null(henderson) || henderson$n_tips != n) return(NULL)

  N <- henderson$N
  Q <- henderson$Q
  T_idx <- henderson$tip_idx
  I_idx <- henderson$int_idx

  # Sigma^{-1}, symmetrised and ridged.
  Sig <- (Sigma + t(Sigma)) / 2
  Sig_inv <- tryCatch(solve(Sig + diag(eps, K)), error = function(e) NULL)
  if (is.null(Sig_inv)) return(NULL)
  Sig_inv <- (Sig_inv + t(Sig_inv)) / 2

  # ---- index sets in the lifted (trait, node) space ----------------------
  obs_cells <- integer(0); obs_vals <- numeric(0)
  unk_cells <- integer(0)
  # positions of the missing TIP cells inside `unk_cells`, plus where they
  # belong in the output matrix
  miss_row <- integer(0); miss_col <- integer(0)

  for (t in seq_len(K)) {
    y <- L[, t]
    o <- which(!is.na(y)); m <- which(is.na(y))
    if (length(o)) {
      yo <- y[o]
      if (isTRUE(cor_scale)) yo <- yo * henderson$tip_sqrt_d[o]
      obs_cells <- c(obs_cells, .ec_cell(t, T_idx[o], N))
      obs_vals  <- c(obs_vals, yo)
    }
    if (length(m)) {
      unk_cells <- c(unk_cells, .ec_cell(t, T_idx[m], N))
      miss_row  <- c(miss_row, m)
      miss_col  <- c(miss_col, rep(t, length(m)))
    }
  }
  n_miss_tip <- length(unk_cells)          # these come FIRST in `unk`
  if (n_miss_tip == 0L) {
    return(list(mu = L, var = matrix(0, n, K, dimnames = dimnames(L))))
  }
  if (!length(obs_cells)) return(NULL)
  # every internal-node cell, all traits, is unobserved
  for (t in seq_len(K)) unk_cells <- c(unk_cells, .ec_cell(t, I_idx, N))
  n_unk <- length(unk_cells)
  if (n_unk > max_cells) return(NULL)

  # ---- joint precision and the conditional solve -------------------------
  P <- tryCatch(Matrix::kronecker(Matrix::Matrix(Sig_inv, sparse = FALSE), Q),
                error = function(e) NULL)
  if (is.null(P)) return(NULL)
  P_uu <- P[unk_cells, unk_cells, drop = FALSE]
  P_uo <- P[unk_cells, obs_cells, drop = FALSE]

  chol_uu <- tryCatch(
    Matrix::Cholesky(Matrix::forceSymmetric(P_uu), LDL = FALSE),
    error = function(e) NULL)
  if (is.null(chol_uu)) return(NULL)

  v <- tryCatch(
    as.numeric(Matrix::solve(chol_uu, -as.matrix(P_uo %*% obs_vals))),
    error = function(e) NULL)
  if (is.null(v) || !all(is.finite(v[seq_len(n_miss_tip)]))) return(NULL)

  var_unit <- tryCatch(
    .ec_inv_diag(chol_uu, n_unk, seq_len(n_miss_tip), block = block),
    error = function(e) NULL)
  if (is.null(var_unit)) return(NULL)

  # ---- write back, undoing the correlation scaling -----------------------
  mu <- L
  var_out <- matrix(0, n, K, dimnames = dimnames(L))
  mu_m <- v[seq_len(n_miss_tip)]
  var_m <- pmax(var_unit, 0)
  if (isTRUE(cor_scale)) {
    sd_r <- henderson$tip_sqrt_d[miss_row]
    mu_m <- mu_m / sd_r
    var_m <- var_m / (sd_r^2)
  }
  mu[cbind(miss_row, miss_col)] <- mu_m
  var_out[cbind(miss_row, miss_col)] <- var_m
  if (!all(is.finite(mu[cbind(miss_row, miss_col)]))) return(NULL)
  list(mu = mu, var = var_out)
}
