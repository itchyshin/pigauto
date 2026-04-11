# Internal helpers for univariate phylogenetic BM imputation.
#
# These implement the standard conditional multivariate-normal approach:
#   R = cov2cor(vcv(tree))       -- phylogenetic correlation matrix
#   mu_hat = GLS phylogenetic mean
#   sigma2 = REML residual variance
#   E[y_m | y_o] = mu_hat + R_mo R_oo^{-1} (y_o - mu_hat)
#   Var[y_m | y_o] = sigma2 * (1 - diag(R_mo R_oo^{-1} R_om))
#
# Each trait column is imputed independently (univariate, not joint).
# This is the standard MCMCglmm/glmmTMB approach: separate the
# phylogenetic correlation structure from the residual variance.
#
# @references
# Goolsby, E.W., Bruggeman, J. & Ane, C. (2017) Rphylopars: fast
#   multivariate phylogenetic comparative methods for missing data and
#   within-species variation. Methods in Ecology and Evolution, 8,
#   22-27.
# @keywords internal

# ---- phylo_cor_matrix -------------------------------------------------------

#' Compute the phylogenetic correlation matrix from a tree
#'
#' @param tree object of class \code{"phylo"}.
#' @return Symmetric positive-definite matrix (\eqn{n \times n}), with
#'   diagonal = 1 and dimnames from \code{tree$tip.label}.
#' @keywords internal
phylo_cor_matrix <- function(tree) {
  C <- ape::vcv(tree)
  stats::cov2cor(C)
}

# ---- bm_impute_col ----------------------------------------------------------

#' Univariate BM conditional imputation for one latent column
#'
#' Given a phylogenetic correlation matrix and a partially observed trait
#' vector, compute the GLS phylogenetic mean, REML variance, and
#' conditional mean/SE for missing species.
#'
#' @param y numeric vector (length \code{n_species}). \code{NA} = missing.
#' @param R phylogenetic correlation matrix (\code{n_species x n_species}).
#'   Must be ordered consistently with \code{y}.
#' @param nugget numeric scalar added to diagonal of \eqn{R_{oo}} for
#'   regularisation (default \code{1e-6}).
#' @return A list with:
#'   \describe{
#'     \item{mu}{numeric vector (length \code{n_species}). Observed species
#'       keep their values; missing species get the conditional mean.}
#'     \item{se}{numeric vector (length \code{n_species}). Observed species
#'       get 0; missing species get
#'       \eqn{\sqrt{\sigma^2 (1 - h_i)}}.}
#'   }
#' @keywords internal
bm_impute_col <- function(y, R, nugget = 1e-6) {
  n <- length(y)
  mu_out <- numeric(n)
  se_out <- numeric(n)

  obs_idx  <- which(!is.na(y))
  miss_idx <- which(is.na(y))
  n_o <- length(obs_idx)
  n_m <- length(miss_idx)

 # Edge case: all observed ------------------------------------------------
  if (n_m == 0L) {
    return(list(mu = y, se = se_out))
  }

  # Edge case: too few observations for GLS --------------------------------
  if (n_o < 5L) {
    y_obs <- y[obs_idx]
    global_mu <- if (n_o > 0L) mean(y_obs) else 0
    global_sd <- if (n_o > 1L) stats::sd(y_obs) else 1
    mu_out[obs_idx]  <- y_obs
    mu_out[miss_idx] <- global_mu
    se_out[miss_idx] <- global_sd
    return(list(mu = mu_out, se = se_out))
  }

  # Normal path: conditional MVN -------------------------------------------
  y_o <- y[obs_idx]
  R_oo <- R[obs_idx, obs_idx, drop = FALSE]
  R_mo <- R[miss_idx, obs_idx, drop = FALSE]

  # Cholesky with nugget back-off
  L <- NULL
  nug <- nugget
  for (attempt in seq_len(6L)) {
    R_oo_reg <- R_oo + diag(nug, n_o)
    L <- tryCatch(chol(R_oo_reg), error = function(e) NULL)
    if (!is.null(L)) break
    nug <- nug * 2
  }
  if (is.null(L)) {
    # Last resort: fallback to global mean/sd
    warning("BM baseline: Cholesky of R_oo failed after nugget back-off; ",
            "using global mean fallback.", call. = FALSE)
    mu_out[obs_idx]  <- y_o
    mu_out[miss_idx] <- mean(y_o)
    se_out[miss_idx] <- stats::sd(y_o)
    return(list(mu = mu_out, se = se_out))
  }

  # Solve via Cholesky: solve(R_oo_reg, b) = backsolve(L, forwardsolve(t(L), b))
  chol_solve <- function(b) {
    backsolve(L, forwardsolve(t(L), b))
  }

  # 1. GLS phylogenetic mean
  ones  <- rep(1, n_o)
  a     <- chol_solve(ones)        # R_oo^{-1} 1
  b     <- chol_solve(y_o)         # R_oo^{-1} y_o
  mu_hat <- sum(b) / sum(a)        # (1'R_oo^{-1} y_o) / (1'R_oo^{-1} 1)

  # 2. REML variance
  e       <- y_o - mu_hat
  e_solve <- chol_solve(e)         # R_oo^{-1} (y_o - mu_hat)
  sigma2  <- as.numeric(crossprod(e, e_solve)) / max(n_o - 1L, 1L)

  # 3. Conditional mean for missing species
  mu_m <- mu_hat + as.numeric(R_mo %*% e_solve)

  # 4. Conditional variance for missing species
  #    h_i = diag(R_mo R_oo^{-1} R_om) computed without forming the full product
  #    W = L^{-T} R_om  =>  h_i = colSums(W^2)
  W <- forwardsolve(t(L), t(R_mo))   # (n_o x n_m)
  h <- colSums(W^2)
  cond_var <- sigma2 * pmax(1 - h, 0)

  # Assemble output
  mu_out[obs_idx]  <- y_o
  mu_out[miss_idx] <- mu_m
  se_out[miss_idx] <- sqrt(cond_var)

  list(mu = mu_out, se = se_out)
}
