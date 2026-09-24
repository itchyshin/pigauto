# R/mi_posterior.R
#
# Posterior multiple imputation for continuous traits:
# multi_impute(draws_method = "posterior").
#
# Specification: docs/dev-log/mi-posterior/design.md (sections 1 to 4).
#
# MODEL (continuous traits on pigauto's latent scale, one row per species)
#
#   y_ik = mu_k + a_ik + e_ik
#   vec(a) ~ N(0, Sigma_P %x% Qc^{-1})    over all N non-root nodes
#   e_i    ~ N(0, Sigma_E)                independently across tips
#   mu     flat
#
# so that vec(Y) ~ N(1 mu', Sigma_P %x% R + Sigma_E %x% I_n) with
# R = cov2cor(vcv(tree)). Qc = D Q D, where Q is the Hadfield & Nakagawa
# (2010, J. Evol. Biol. 23: 494-508) sparse node precision from
# build_henderson_S_inv() (root excluded, so mu is the root state) and
# D = diag(tip_sqrt_d at tips, 1 at internal nodes). The tip-block Schur
# complement of Qc is R^{-1} exactly, so tip values of `a` are in the units
# of y and the internal nodes sit on an irrelevant rescaled scale. Qc is the
# only node precision used anywhere below.
#
# GIBBS SWEEP (design.md section 2.2)
#
#   1. (a, mu, y_mis) | y_obs, Sigma_P, Sigma_E as ONE block: first
#      (a, mu) | y_obs with y_mis integrated out (sparse precision
#      Qc %x% Sigma_P^{-1} plus per-tip likelihood blocks
#      P_i = E_O (Sigma_E[O, O])^{-1} E_O'), then y_mis | a, mu, y_obs per
#      tip. Sparse Cholesky with a FIXED pattern: the symbolic analysis is
#      done once, and every sweep only rewrites the numeric values and calls
#      Matrix::update().
#   2. xi = a diag(1 / alpha)                     (working effects)
#   3. alpha | xi, mu, y, Sigma_E                 (Gaussian, conjugate)
#   4. Sigma_W | xi ~ IW(nu_W + N, S_W + xi' Qc xi)
#      Sigma_P = diag(alpha) Sigma_W diag(alpha)
#   5. Sigma_E | e ~ IW(nu_E + n, S_E + e'e)
#
# PARAMETER EXPANSION. Steps 2 to 4 are the parameter-expanded
# data-augmentation scheme for the phylogenetic covariance as used for
# MCMCglmm's G structures (alpha.mu = 0, alpha.V = V_alpha I): Gelman (2006,
# Bayesian Analysis 1: 515-534) for the scaled-F prior it implies on each
# standard deviation, and Hadfield (2010, J. Stat. Softw. 33(2)) for the
# multivariate form. Priors: Sigma_W ~ IW(K + 1, I_K), alpha ~ N(0, 1000 I_K),
# Sigma_E ~ IW(K + 1, 0.01 diag(observed latent variances)). alpha and
# Sigma_W are not identified; only Sigma_P, Sigma_E, lambda and mu are
# reported and diagnosed.
#
# DIAGNOSTICS. Rank-normalised split R-hat (the maximum of the bulk and the
# folded-tail versions) and rank-normalised bulk ESS with Geyer's initial
# monotone sequence, following Vehtari, Gelman, Simpson, Carpenter & Buerkner
# (2021, Bayesian Analysis 16: 667-718). Implemented here to avoid a new
# dependency; they are computed on every post-burn-in sweep of the
# parameters (not only on the thinned, kept sweeps).
#
# IMPROPER MODE (param_uncertainty = "none" or "both", validation only).
# "both" returns the proper draws exactly as "full" does, plus the plug-in
# draws below from the SAME chain run (drawn after the chains, so the proper
# results are identical to "full" with the same seed). The full
# sampler runs first (same chains and lengths; its diagnostics are the ones
# reported). Sigma_P and Sigma_E are then fixed at their posterior means and
# only step 1 runs; with the covariances fixed the precision does not change,
# so it is factorised once and every kept draw is an exact independent draw
# of (a, mu, y_mis) | y_obs. mu is still drawn (flat prior), as in a REML
# plug-in. The posterior mean was chosen over #187's REML estimate because
# #187 estimates only the per-trait lambda (diagonals), not the two full
# K x K covariances, so the posterior mean is the only plug-in that matches
# the full model and isolates parameter uncertainty in the comparison.

# ---------------------------------------------------------------------------
# Problem set-up
# ---------------------------------------------------------------------------

# Qc = D Q D (design.md 2.0). Zero-length edges are floored at 1e-6 of the
# tree height before Q is built so that 1 / edge length stays finite. The
# tip depths are passed in from ape::node.depth.edgelength() (O(n)), so the
# dense n x n vcv(tree) is never formed here; the values are the same as
# diag(vcv(tree)).
.mip_build_Qc <- function(tree) {
  tip_depths <- NULL
  if (!is.null(tree$edge.length)) {
    h <- max(ape::node.depth.edgelength(tree))
    floor_len <- 1e-6 * if (is.finite(h) && h > 0) h else 1
    tree$edge.length <- pmax(tree$edge.length, floor_len)
    tip_depths <- ape::node.depth.edgelength(tree)[
      seq_len(length(tree$tip.label))]
  }
  hs <- build_henderson_S_inv(tree, tip_depths = tip_depths)
  d <- c(hs$tip_sqrt_d, rep(1, hs$n_internal_nonroot))
  Dm <- Matrix::Diagonal(x = d)
  Qc <- Matrix::forceSymmetric(Dm %*% hs$Q %*% Dm, uplo = "U")
  list(Qc = Qc, n = hs$n_tips, N = hs$N, tip_labels = hs$tip_labels,
       tip_sqrt_d = hs$tip_sqrt_d)
}

# Group tips by their observed-trait pattern.
.mip_patterns <- function(Y) {
  obs <- !is.na(Y)
  key <- apply(obs, 1L, function(r) paste(as.integer(r), collapse = ""))
  ukey <- unique(key)
  lapply(ukey, function(kk) {
    rows <- which(key == kk)
    o <- which(obs[rows[1L], ])
    list(rows = rows, obs = o, miss = setdiff(seq_len(ncol(Y)), o))
  })
}

# Build the problem object from a latent matrix Y (n x K, NA = missing,
# rows aligned to tree$tip.label) and a tree.
.mip_problem <- function(Y, tree) {
  Y <- as.matrix(Y)
  qc <- .mip_build_Qc(tree)
  if (!is.null(rownames(Y)) && !identical(rownames(Y), qc$tip_labels)) {
    stop("mi_posterior: rows of the latent matrix are not in tree tip order ",
         "(internal error).", call. = FALSE)
  }
  if (nrow(Y) != qc$n) {
    stop("mi_posterior: latent matrix and tree disagree on the number of ",
         "tips (internal error).", call. = FALSE)
  }
  n_obs_k <- colSums(!is.na(Y))
  if (any(n_obs_k < 2L)) {
    stop("draws_method = \"posterior\" needs at least 2 observed values per ",
         "trait; too few in: ",
         paste(colnames(Y)[n_obs_k < 2L], collapse = ", "), ".", call. = FALSE)
  }
  pats <- .mip_patterns(Y)
  miss <- which(is.na(Y), arr.ind = TRUE)
  miss <- miss[order(miss[, 2L], miss[, 1L]), , drop = FALSE]
  list(Y = Y, Y0 = replace(Y, is.na(Y), 0), K = ncol(Y), n = qc$n, N = qc$N,
       Qc = qc$Qc, tip_sqrt_d = qc$tip_sqrt_d, patterns = pats,
       miss = unname(miss),
       obs_var = apply(Y, 2L, stats::var, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# The (a, mu) | y_obs precision with a fixed sparsity pattern
# ---------------------------------------------------------------------------
#
# System ordering: node-major a (index (v - 1) K + k, v = 1..N), then mu
# (N K + k) when include_mu. Every stored entry of the upper triangle is a
# linear function of theta = c(vec(Sigma_P^{-1}), vec(P_1), ..., vec(P_G)),
# where P_g is the likelihood block of observed-pattern g. The map is the
# sparse matrix A (nnz x length(theta)), built once: P@x <- A %*% theta.
.mip_template <- function(prob, include_mu = TRUE) {
  K <- prob$K; N <- prob$N
  n_sys <- N * K + if (include_mu) K else 0L
  Qt <- Matrix::mat2triplet(Matrix::triu(prob$Qc))
  up_q <- Qt$i <= Qt$j
  qi <- as.integer(Qt$i[up_q]); qj <- as.integer(Qt$j[up_q]); qx <- Qt$x[up_q]
  kk <- rep(seq_len(K), times = K); ll <- rep(seq_len(K), each = K)
  th_sig <- (ll - 1L) * K + kk
  # Prior: Qc[v, w] * SigPinv[k, l] at ((v-1)K+k, (w-1)K+l), upper only.
  ri <- rep((qi - 1L) * K, each = K * K) + rep(kk, length(qi))
  ci <- rep((qj - 1L) * K, each = K * K) + rep(ll, length(qi))
  ti <- rep(th_sig, length(qi))
  cf <- rep(qx, each = K * K)
  keep <- ri <= ci
  ri <- ri[keep]; ci <- ci[keep]; ti <- ti[keep]; cf <- cf[keep]
  # Likelihood blocks.
  off <- K * K
  for (g in seq_along(prob$patterns)) {
    pg <- prob$patterns[[g]]
    o <- pg$obs
    if (!length(o)) next
    base <- off + (g - 1L) * K * K
    ok <- rep(o, times = length(o)); ol <- rep(o, each = length(o))
    thg <- base + (ol - 1L) * K + ok
    tips <- pg$rows
    # (a_i k, a_i l), k <= l
    su <- ok <= ol
    ri <- c(ri, rep((tips - 1L) * K, each = sum(su)) + rep(ok[su], length(tips)))
    ci <- c(ci, rep((tips - 1L) * K, each = sum(su)) + rep(ol[su], length(tips)))
    ti <- c(ti, rep(thg[su], length(tips)))
    cf <- c(cf, rep(1, sum(su) * length(tips)))
    if (include_mu) {
      # (a_i k, mu l), all k, l
      ri <- c(ri, rep((tips - 1L) * K, each = length(ok)) + rep(ok, length(tips)))
      ci <- c(ci, rep(N * K + ol, length(tips)))
      ti <- c(ti, rep(thg, length(tips)))
      cf <- c(cf, rep(1, length(ok) * length(tips)))
      # (mu k, mu l), k <= l, summed over the pattern's tips
      ri <- c(ri, N * K + ok[su]); ci <- c(ci, N * K + ol[su])
      ti <- c(ti, thg[su]); cf <- c(cf, rep(length(tips), sum(su)))
    }
  }
  P <- Matrix::sparseMatrix(i = ri, j = ci, x = 1, dims = c(n_sys, n_sys),
                            symmetric = TRUE, repr = "C")
  colj <- rep(seq_len(n_sys), diff(P@p))
  skey <- (colj - 1) * n_sys + (P@i + 1)
  tkey <- (ci - 1) * n_sys + ri
  pos <- match(tkey, skey)
  if (anyNA(pos)) stop("mi_posterior: template indexing failed (internal).",
                       call. = FALSE)
  A <- Matrix::sparseMatrix(i = pos, j = ti, x = cf,
                            dims = c(length(P@x),
                                     K * K * (1L + length(prob$patterns))))
  list(P = P, A = A, include_mu = include_mu, n_sys = n_sys, L = NULL)
}

# Likelihood blocks P_g (K x K, zero outside the observed set).
.mip_lik_blocks <- function(prob, Sigma_E) {
  K <- prob$K
  lapply(prob$patterns, function(pg) {
    Pg <- matrix(0, K, K)
    o <- pg$obs
    if (length(o)) Pg[o, o] <- .mip_solve_spd(Sigma_E[o, o, drop = FALSE])
    Pg
  })
}

.mip_solve_spd <- function(S) {
  S <- (S + t(S)) / 2
  chol2inv(chol(S))
}

# Refresh the numeric values of the precision (and its factor).
.mip_refactor <- function(tpl, Sigma_P, lik_blocks) {
  theta <- c(as.vector(.mip_solve_spd(Sigma_P)),
             unlist(lapply(lik_blocks, as.vector), use.names = FALSE))
  tpl$P@x <- as.numeric(tpl$A %*% theta)
  tpl$L <- if (is.null(tpl$L)) {
    Matrix::Cholesky(tpl$P, perm = TRUE, LDL = FALSE, super = FALSE)
  } else {
    Matrix::update(tpl$L, tpl$P)
  }
  tpl
}

# Linear term: sum over tips of P_i (y_i - mu_fixed) (mu_fixed = 0 when mu is
# part of the system).
.mip_linear <- function(prob, tpl, lik_blocks, mu_fixed = NULL) {
  K <- prob$K; N <- prob$N
  Yc <- prob$Y0
  if (!is.null(mu_fixed)) {
    Yc <- sweep(prob$Y, 2L, mu_fixed)
    Yc[is.na(Yc)] <- 0
  }
  B <- matrix(0, N, K)
  for (g in seq_along(prob$patterns)) {
    pg <- prob$patterns[[g]]
    if (!length(pg$obs)) next
    B[pg$rows, ] <- Yc[pg$rows, , drop = FALSE] %*% lik_blocks[[g]]
  }
  b <- as.vector(t(B))
  if (tpl$include_mu) b <- c(b, colSums(B))
  b
}

# Posterior mean of (a, mu) and n_draws joint draws around it.
.mip_amu <- function(prob, tpl, b, n_draws = 1L) {
  K <- prob$K; N <- prob$N
  m <- as.numeric(Matrix::solve(tpl$L, b, system = "A"))
  split <- function(v) {
    list(a = matrix(v[seq_len(N * K)], N, K, byrow = TRUE),
         mu = if (tpl$include_mu) v[N * K + seq_len(K)] else NULL)
  }
  if (n_draws == 0L) return(split(m))
  W <- matrix(stats::rnorm(tpl$n_sys * n_draws), tpl$n_sys, n_draws)
  Z <- Matrix::solve(tpl$L, Matrix::solve(tpl$L, W, system = "Lt"),
                     system = "Pt")
  Z <- as.matrix(Z) + m
  if (n_draws == 1L) return(split(Z[, 1L]))
  lapply(seq_len(n_draws), function(d) split(Z[, d]))
}

# y_mis | a, mu, y_obs, Sigma_E per tip (step 1b). With draw = FALSE returns
# the conditional mean (no noise).
.mip_fill_ymis <- function(prob, a_tip, mu, Sigma_E, draw = TRUE) {
  Y <- prob$Y
  for (pg in prob$patterns) {
    M <- pg$miss
    if (!length(M)) next
    rows <- pg$rows; O <- pg$obs
    base <- a_tip[rows, M, drop = FALSE] +
      matrix(mu[M], length(rows), length(M), byrow = TRUE)
    if (length(O)) {
      eO <- Y[rows, O, drop = FALSE] - a_tip[rows, O, drop = FALSE] -
        matrix(mu[O], length(rows), length(O), byrow = TRUE)
      Bm <- Sigma_E[M, O, drop = FALSE] %*% .mip_solve_spd(Sigma_E[O, O, drop = FALSE])
      base <- base + eO %*% t(Bm)
      Cm <- Sigma_E[M, M, drop = FALSE] - Bm %*% Sigma_E[O, M, drop = FALSE]
    } else {
      Cm <- Sigma_E[M, M, drop = FALSE]
    }
    if (draw) {
      Cm <- (Cm + t(Cm)) / 2
      Z <- matrix(stats::rnorm(length(rows) * length(M)), length(rows))
      base <- base + Z %*% chol(Cm)
    }
    Y[rows, M] <- base
  }
  Y
}

# Exact conditional mean of y_mis | y_obs at fixed parameters (mu drawn or
# fixed). Used by the G2b gate and the tests.
.mip_cond_mean <- function(prob, Sigma_P, Sigma_E, mu = NULL) {
  tpl <- .mip_template(prob, include_mu = is.null(mu))
  lb <- .mip_lik_blocks(prob, Sigma_E)
  tpl <- .mip_refactor(tpl, Sigma_P, lb)
  s <- .mip_amu(prob, tpl, .mip_linear(prob, tpl, lb, mu_fixed = mu), 0L)
  mu_use <- if (is.null(mu)) s$mu else mu
  Yc <- .mip_fill_ymis(prob, s$a[seq_len(prob$n), , drop = FALSE], mu_use,
                       Sigma_E, draw = FALSE)
  Yc[prob$miss]
}

# n_draws exact independent draws of y_mis | y_obs at fixed Sigma_P, Sigma_E
# (and mu when supplied). Returns list(ymis = n_mis x n_draws, mu = n_draws x K).
.mip_fixed_draws <- function(prob, Sigma_P, Sigma_E, n_draws, mu = NULL,
                             block = 500L) {
  tpl <- .mip_template(prob, include_mu = is.null(mu))
  lb <- .mip_lik_blocks(prob, Sigma_E)
  tpl <- .mip_refactor(tpl, Sigma_P, lb)
  b <- .mip_linear(prob, tpl, lb, mu_fixed = mu)
  n_mis <- nrow(prob$miss)
  out <- matrix(NA_real_, n_mis, n_draws)
  mus <- matrix(NA_real_, n_draws, prob$K)
  done <- 0L
  while (done < n_draws) {
    nb <- min(block, n_draws - done)
    dr <- .mip_amu(prob, tpl, b, nb)
    if (nb == 1L) dr <- list(dr)
    for (d in seq_len(nb)) {
      mu_d <- if (is.null(mu)) dr[[d]]$mu else mu
      Yc <- .mip_fill_ymis(prob, dr[[d]]$a[seq_len(prob$n), , drop = FALSE],
                           mu_d, Sigma_E)
      out[, done + d] <- Yc[prob$miss]
      mus[done + d, ] <- mu_d
    }
    done <- done + nb
  }
  list(ymis = out, mu = mus)
}

# ---------------------------------------------------------------------------
# Collapsed log-likelihood log p(y_obs | Sigma_P, Sigma_E), a and mu
# integrated out (flat mu prior), up to a constant free of the parameters:
#
#   -1/2 sum_i log|Sigma_E[O_i, O_i]| - N/2 log|Sigma_P| - 1/2 log|Lambda|
#   -1/2 (sum_i y_i' P_i y_i - b' Lambda^{-1} b)
#
# where Lambda and b are the precision and linear term of (a, mu) | y_obs.
# (log|Qc %x% Sigma_P^{-1}| = K log|Qc| - N log|Sigma_P|; the first term is
# constant.) Returns list(ll, tpl) because the factor is refreshed in place.
# ---------------------------------------------------------------------------
.mip_collapsed_ll <- function(prob, tpl, Sigma_P, Sigma_E) {
  lb <- .mip_lik_blocks(prob, Sigma_E)
  tpl <- .mip_refactor(tpl, Sigma_P, lb)
  b <- .mip_linear(prob, tpl, lb)
  m <- as.numeric(Matrix::solve(tpl$L, b, system = "A"))
  ld_Lam <- as.numeric(Matrix::determinant(tpl$L, logarithm = TRUE,
                                           sqrt = FALSE)$modulus)
  ld_E <- 0; q_y <- 0
  for (g in seq_along(prob$patterns)) {
    pg <- prob$patterns[[g]]
    o <- pg$obs
    if (!length(o)) next
    ld_E <- ld_E + length(pg$rows) *
      as.numeric(determinant(Sigma_E[o, o, drop = FALSE])$modulus)
    Yg <- prob$Y0[pg$rows, , drop = FALSE]
    q_y <- q_y + sum((Yg %*% lb[[g]]) * Yg)
  }
  ld_P <- as.numeric(determinant(Sigma_P)$modulus)
  ll <- -0.5 * ld_E - 0.5 * prob$N * ld_P - 0.5 * ld_Lam -
    0.5 * (q_y - sum(b * m))
  list(ll = ll, tpl = tpl)
}

.mip_log_iw <- function(S, nu, Psi) {
  K <- nrow(S)
  -0.5 * (nu + K + 1) * as.numeric(determinant(S)$modulus) -
    0.5 * sum(diag(Psi %*% .mip_solve_spd(S)))
}

# Collapsed Metropolis moves (added to the design.md sweep, before step 1):
# for each trait k, three multiplicative moves on the expanded state with
# (a, mu, y_mis) integrated out:
#   "alpha": alpha_k -> c alpha_k                (log c ~ N(0, s))
#   "E"    : row/column k of Sigma_E scaled by d (log d ~ N(0, s))
#   "ridge": both, with log c = -log d           (trades Sigma_P[k,k]
#            against Sigma_E[k,k], i.e. moves lambda_k)
# Each is a group move with a symmetric proposal on the log scale, so the
# acceptance ratio carries the Jacobians |c| (alpha_k) and d^(K+1) (the
# symmetric K x K Sigma_E). Target: p(alpha, Sigma_E | Sigma_W, y_obs), the
# (a, mu, y_mis)-marginal. Because step 1 then redraws (a, mu, y_mis) from
# its full conditional, the sweep is a valid partially collapsed Gibbs
# sampler (van Dyk & Park 2008, JASA 103: 790-796). Why: near lambda = 1 the
# Gibbs update Sigma_E | a is nearly deterministic (Sigma_E small, a tied to
# y), and near lambda = 0 the same holds for Sigma_P; the expansion covers
# only Sigma_P. Measured at n = 1000 without these moves: bulk ESS 20 to 41
# for lambda from 4 x 10,000 sweeps at truth lambda = 1 and 0.95.
.mip_mh_moves <- function(prob, tpl, alpha, Sigma_W, Sigma_E, hyper, step,
                          step_off, cur_ll = NULL) {
  K <- prob$K
  sP <- function(al) .mip_px_sigma(al, Sigma_W)
  if (is.null(cur_ll)) {
    cc <- .mip_collapsed_ll(prob, tpl, sP(alpha), Sigma_E)
    cur_ll <- cc$ll; tpl <- cc$tpl
  }
  acc <- matrix(0, K, 3L)
  sdA <- sqrt(hyper$V_alpha)
  for (k in seq_len(K)) {
    for (mv in 1:3) {
      z <- stats::rnorm(1L, 0, step[k, mv])
      lc <- switch(mv, z, 0, z)
      ld <- switch(mv, 0, z, -z)
      al_new <- alpha
      al_new[k] <- alpha[k] * exp(lc)
      SE_new <- Sigma_E
      if (ld != 0) {
        sc <- rep(1, K); sc[k] <- exp(ld)
        SE_new <- Sigma_E * tcrossprod(sc)
      }
      pr <- .mip_collapsed_ll(prob, tpl, sP(al_new), SE_new)
      tpl <- pr$tpl
      log_r <- pr$ll - cur_ll +
        stats::dnorm(al_new[k], 0, sdA, log = TRUE) -
        stats::dnorm(alpha[k], 0, sdA, log = TRUE) + lc +
        .mip_log_iw(SE_new, hyper$nu_E, hyper$S_E) -
        .mip_log_iw(Sigma_E, hyper$nu_E, hyper$S_E) + (K + 1) * ld
      if (is.finite(log_r) && log(stats::runif(1L)) < log_r) {
        alpha <- al_new; Sigma_E <- SE_new; cur_ll <- pr$ll
        acc[k, mv] <- 1
      }
    }
  }
  # Prior-only rebalance per trait k: row/column k of Sigma_W scaled by d
  # and alpha_k by 1 / d, so Sigma_P (hence the likelihood) is unchanged and
  # no factorisation is needed. Jacobian d^(K+1) (Sigma_W) times 1/d
  # (alpha_k). It keeps the expanded chain irreducible in Sigma_W's
  # diagonal when the Gibbs steps are switched off (test hook).
  for (k in seq_len(K)) {
    ld <- stats::rnorm(1L, 0, 0.5)
    sc <- rep(1, K); sc[k] <- exp(ld)
    SW_new <- Sigma_W * tcrossprod(sc)
    al_new <- alpha; al_new[k] <- alpha[k] / exp(ld)
    log_r <- .mip_log_iw(SW_new, hyper$nu_W, hyper$S_W) -
      .mip_log_iw(Sigma_W, hyper$nu_W, hyper$S_W) +
      stats::dnorm(al_new[k], 0, sdA, log = TRUE) -
      stats::dnorm(alpha[k], 0, sdA, log = TRUE) + K * ld
    if (is.finite(log_r) && log(stats::runif(1L)) < log_r) {
      Sigma_W <- SW_new; alpha <- al_new
    }
  }

  # Off-diagonal moves, one set per trait pair (k < l):
  #   1: correlation of Sigma_E[k,l] on the Fisher-z scale,
  #      atanh(r) -> atanh(r) + eps, diagonal fixed; Jacobian
  #      (1 - r_new^2) / (1 - r_old^2)
  #   2: the same for Sigma_W[k,l]
  #   3: covariance ridge: Sigma_P[k,l] += delta and Sigma_E[k,l] -= delta,
  #      delta = eps sqrt(v_k v_l), v = diag(Sigma_P + Sigma_E); in the
  #      expanded state Sigma_W[k,l] += delta / (alpha_k alpha_l). The total
  #      cross-trait covariance is well identified but its phylogenetic /
  #      residual split is not, which Gibbs steps 4-5 traverse slowly
  #      (measured: bulk ESS 58 to 140 for Sigma_P[1,2] at n = 1000).
  # eps ~ N(0, s). Each move leaves the diagonals and alpha unchanged and
  # its scale depends only on those, so the proposal is symmetric in the
  # moved coordinate (move 3 is a translation with unit Jacobian); a
  # non-positive-definite proposal has zero prior density and is rejected.
  n_pair <- K * (K - 1L) / 2L
  acc_off <- matrix(0, max(0L, n_pair), 3L)
  if (K > 1L) {
    pairs <- which(upper.tri(diag(K)), arr.ind = TRUE)
    not_pd <- function(S) {
      inherits(tryCatch(chol(S), error = function(e) e), "error")
    }
    for (p in seq_len(nrow(pairs))) {
      k <- pairs[p, 1L]; l <- pairs[p, 2L]
      for (which_m in 1:3) {
        z <- stats::rnorm(1L, 0, step_off[p, which_m])
        SE_new <- Sigma_E; SW_new <- Sigma_W
        log_jac <- 0
        if (which_m <= 2L) {
          S <- if (which_m == 1L) Sigma_E else Sigma_W
          sc <- sqrt(S[k, k] * S[l, l])
          r_old <- S[k, l] / sc
          if (abs(r_old) >= 1) next
          r_new <- tanh(atanh(r_old) + z)
          log_jac <- log1p(-r_new^2) - log1p(-r_old^2)
          if (which_m == 1L) {
            SE_new[k, l] <- SE_new[l, k] <- r_new * sc
          } else {
            SW_new[k, l] <- SW_new[l, k] <- r_new * sc
          }
        } else {
          SPc <- sP(alpha)
          v <- diag(SPc) + diag(Sigma_E)
          delta <- z * sqrt(v[k] * v[l])
          SE_new[k, l] <- SE_new[l, k] <- Sigma_E[k, l] - delta
          SW_new[k, l] <- SW_new[l, k] <-
            Sigma_W[k, l] + delta / (alpha[k] * alpha[l])
        }
        if (not_pd(SE_new) || not_pd(SW_new)) next
        pr <- .mip_collapsed_ll(prob, tpl, .mip_px_sigma(alpha, SW_new),
                                SE_new)
        tpl <- pr$tpl
        log_r <- pr$ll - cur_ll +
          .mip_log_iw(SE_new, hyper$nu_E, hyper$S_E) -
          .mip_log_iw(Sigma_E, hyper$nu_E, hyper$S_E) +
          .mip_log_iw(SW_new, hyper$nu_W, hyper$S_W) -
          .mip_log_iw(Sigma_W, hyper$nu_W, hyper$S_W) + log_jac
        if (is.finite(log_r) && log(stats::runif(1L)) < log_r) {
          Sigma_E <- SE_new; Sigma_W <- SW_new
          cur_ll <- pr$ll
          acc_off[p, which_m] <- 1
        }
      }
    }
  }
  list(alpha = alpha, Sigma_E = Sigma_E, Sigma_W = Sigma_W, tpl = tpl,
       acc = acc, acc_off = acc_off, ll = cur_ll)
}

# ---------------------------------------------------------------------------
# Inverse-Wishart and the PX steps
# ---------------------------------------------------------------------------

# PX mapping (design.md 2.1): Sigma_P = diag(alpha) Sigma_W diag(alpha) and
# a = xi diag(alpha) (row i of a is diag(alpha) xi_i).
.mip_px_sigma <- function(alpha, Sigma_W) {
  S <- Sigma_W * tcrossprod(alpha)
  (S + t(S)) / 2
}
.mip_px_effects <- function(xi, alpha) sweep(xi, 2L, alpha, "*")

# Sigma ~ IW(nu, S): density proportional to
# |Sigma|^{-(nu + K + 1) / 2} exp(-tr(S Sigma^{-1}) / 2).
.mip_riwish <- function(nu, S) {
  S <- (S + t(S)) / 2
  W <- stats::rWishart(1L, df = nu, Sigma = .mip_solve_spd(S))[, , 1L]
  out <- .mip_solve_spd(W)
  (out + t(out)) / 2
}

# Step 3: alpha | xi, mu, y (completed), Sigma_E. At tips
# y_i - mu = diag(xi_i) alpha + e_i, e_i ~ N(0, Sigma_E).
.mip_draw_alpha <- function(xi_tip, Rres, Sigma_E, V_alpha) {
  K <- ncol(xi_tip)
  SEi <- .mip_solve_spd(Sigma_E)
  Pa <- SEi * crossprod(xi_tip) + diag(1 / V_alpha, K)
  rhs <- colSums(xi_tip * (Rres %*% SEi))
  Ua <- chol((Pa + t(Pa)) / 2)
  m <- backsolve(Ua, forwardsolve(t(Ua), rhs))
  as.numeric(m + backsolve(Ua, stats::rnorm(K)))
}

# ---------------------------------------------------------------------------
# One chain
# ---------------------------------------------------------------------------
#
# start: list(Sigma_P, Sigma_E). hyper: list(nu_W, S_W, V_alpha, nu_E, S_E).
# Returns the post-burn-in parameter trace (every sweep) and the kept sweeps
# (every thin-th) with y_mis and the full matrices.
.mip_run_chain <- function(prob, start, hyper, burnin, n_iter, thin,
                           seed = NULL, mh = TRUE, gibbs = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  K <- prob$K; n <- prob$n; N <- prob$N
  tpl <- .mip_template(prob, include_mu = TRUE)
  Sigma_E <- start$Sigma_E
  alpha <- rep(1, K)
  Sigma_W <- start$Sigma_P
  Sigma_P <- Sigma_W
  step <- matrix(0.5, K, 3L)
  acc_win <- matrix(0, K, 3L)
  acc_kept <- matrix(0, K, 3L)
  n_pair <- K * (K - 1L) / 2L
  step_off <- matrix(0.2, max(n_pair, 0L), 3L)
  acc_win_off <- matrix(0, max(n_pair, 0L), 3L)
  acc_kept_off <- matrix(0, max(n_pair, 0L), 3L)
  adapt <- function(st, rate, it) {
    st <- st * exp(ifelse(rate > 0.44, 1, -1) * min(0.5, 5 / sqrt(it / 50)))
    pmin(pmax(st, 1e-3), 5)
  }
  n_kept <- n_iter %/% thin
  up <- which(upper.tri(diag(K), diag = TRUE))
  trace <- matrix(NA_real_, n_iter, 2L * length(up) + 2L * K)
  kept_ymis <- matrix(NA_real_, nrow(prob$miss), n_kept)
  kept_SP <- array(NA_real_, c(K, K, n_kept))
  kept_SE <- array(NA_real_, c(K, K, n_kept))
  kept_mu <- matrix(NA_real_, n_kept, K)
  qc <- prob$Qc
  for (it in seq_len(burnin + n_iter)) {
    # 0. Collapsed Metropolis moves on (alpha_k, Sigma_E) (see .mip_mh_moves).
    if (mh) {
      mv <- .mip_mh_moves(prob, tpl, alpha, Sigma_W, Sigma_E, hyper, step,
                          step_off)
      alpha <- mv$alpha; Sigma_E <- mv$Sigma_E; Sigma_W <- mv$Sigma_W
      tpl <- mv$tpl
      if (it <= burnin) {
        # Adapt step sizes in windows of 50 sweeps during burn-in only
        # (target acceptance 0.44 for one-dimensional moves); frozen after,
        # so the kept phase is a fixed, valid Metropolis-within-Gibbs kernel.
        acc_win <- acc_win + mv$acc
        acc_win_off <- acc_win_off + mv$acc_off
        if (it %% 50L == 0L) {
          step <- adapt(step, acc_win / 50, it)
          if (n_pair > 0L) step_off <- adapt(step_off, acc_win_off / 50, it)
          acc_win[] <- 0; acc_win_off[] <- 0
        }
      } else {
        acc_kept <- acc_kept + mv$acc
        acc_kept_off <- acc_kept_off + mv$acc_off
      }
      Sigma_P <- .mip_px_sigma(alpha, Sigma_W)
    }
    # 1. (a, mu, y_mis) | y_obs as one block.
    lb <- .mip_lik_blocks(prob, Sigma_E)
    tpl <- .mip_refactor(tpl, Sigma_P, lb)
    s <- .mip_amu(prob, tpl, .mip_linear(prob, tpl, lb), 1L)
    a <- s$a; mu <- s$mu
    Yc <- .mip_fill_ymis(prob, a[seq_len(n), , drop = FALSE], mu, Sigma_E)
    if (!gibbs) {
      # Test hook: collapsed moves + step 1 only (steps 2-5 skipped).
      Rres <- sweep(Yc, 2L, mu)
    } else {
    # 2. working effects.
    xi <- sweep(a, 2L, alpha, "/")
    # 3. alpha.
    Rres <- sweep(Yc, 2L, mu)
    alpha <- .mip_draw_alpha(xi[seq_len(n), , drop = FALSE], Rres, Sigma_E,
                             hyper$V_alpha)
    a <- .mip_px_effects(xi, alpha)
    # 4. Sigma_W | xi, then Sigma_P.
    SS <- as.matrix(crossprod(xi, qc %*% xi))
    Sigma_W <- .mip_riwish(hyper$nu_W + N, hyper$S_W + SS)
    Sigma_P <- .mip_px_sigma(alpha, Sigma_W)
    # 5. Sigma_E | e.
    e <- Rres - a[seq_len(n), , drop = FALSE]
    Sigma_E <- .mip_riwish(hyper$nu_E + n, hyper$S_E + crossprod(e))
    }
    if (it > burnin) {
      j <- it - burnin
      lam <- diag(Sigma_P) / (diag(Sigma_P) + diag(Sigma_E))
      trace[j, ] <- c(Sigma_P[up], Sigma_E[up], lam, mu)
      if (j %% thin == 0L) {
        kk <- j %/% thin
        kept_ymis[, kk] <- Yc[prob$miss]
        kept_SP[, , kk] <- Sigma_P
        kept_SE[, , kk] <- Sigma_E
        kept_mu[kk, ] <- mu
      }
    }
  }
  list(trace = trace, ymis = kept_ymis, Sigma_P = kept_SP, Sigma_E = kept_SE,
       mu = kept_mu, mh_step = list(scale = step, off = step_off),
       mh_accept = if (mh) list(scale = acc_kept / n_iter,
                                off = acc_kept_off / n_iter) else NULL)
}

.mip_param_names <- function(K) {
  up <- which(upper.tri(diag(K), diag = TRUE), arr.ind = TRUE)
  up <- up[order(up[, 2L], up[, 1L]), , drop = FALSE]
  ij <- sprintf("[%d,%d]", up[, 1L], up[, 2L])
  c(paste0("Sigma_P", ij), paste0("Sigma_E", ij),
    sprintf("lambda[%d]", seq_len(K)), sprintf("mu[%d]", seq_len(K)))
}

# ---------------------------------------------------------------------------
# Convergence diagnostics (Vehtari et al. 2021)
# ---------------------------------------------------------------------------

# sims: iterations x chains. Returns iterations/2 x (2 * chains).
.mip_split_chains <- function(sims) {
  n <- nrow(sims)
  h <- n %/% 2L
  if (h < 2L) return(sims)
  cbind(sims[seq_len(h), , drop = FALSE],
        sims[n - h + seq_len(h), , drop = FALSE])
}

.mip_rank_normalise <- function(sims) {
  r <- rank(as.vector(sims), ties.method = "average")
  S <- length(r)
  matrix(stats::qnorm((r - 3 / 8) / (S + 1 / 4)), nrow(sims), ncol(sims))
}

.mip_rhat_basic <- function(sims) {
  n <- nrow(sims)
  cm <- colMeans(sims)
  cv <- apply(sims, 2L, stats::var)
  W <- mean(cv)
  B <- n * stats::var(cm)
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

# Biased (Geyer 1992) autocovariance via FFT.
.mip_autocov <- function(x) {
  N <- length(x)
  M2 <- 2L * stats::nextn(N)
  yc <- c(x - mean(x), rep(0, M2 - N))
  f <- stats::fft(yc)
  Re(stats::fft(Conj(f) * f, inverse = TRUE))[seq_len(N)] / (M2 * N)
}

.mip_ess_raw <- function(sims) {
  chains <- ncol(sims); n <- nrow(sims)
  if (n < 4L) return(NA_real_)
  acov <- vapply(seq_len(chains), function(c) .mip_autocov(sims[, c]),
                 numeric(n))
  acov <- matrix(acov, n, chains)
  mean_var <- mean(acov[1L, ]) * n / (n - 1)
  if (!is.finite(mean_var) || mean_var <= 0) return(NA_real_)
  var_plus <- mean_var * (n - 1) / n
  if (chains > 1L) var_plus <- var_plus + stats::var(colMeans(sims))
  rho <- rep(0, n)
  t <- 0L
  rho_even <- 1
  rho[1L] <- rho_even
  rho_odd <- 1 - (mean_var - mean(acov[2L, ])) / var_plus
  rho[2L] <- rho_odd
  while (t < n - 5L && !is.nan(rho_even + rho_odd) &&
         (rho_even + rho_odd > 0)) {
    t <- t + 2L
    rho_even <- 1 - (mean_var - mean(acov[t + 1L, ])) / var_plus
    rho_odd <- 1 - (mean_var - mean(acov[t + 2L, ])) / var_plus
    if ((rho_even + rho_odd) >= 0) {
      rho[t + 1L] <- rho_even
      rho[t + 2L] <- rho_odd
    }
  }
  max_t <- t
  if (rho_even > 0) rho[max_t + 1L] <- rho_even
  # Geyer's initial monotone sequence.
  t <- 0L
  while (t <= max_t - 4L) {
    t <- t + 2L
    if (rho[t + 1L] + rho[t + 2L] > rho[t - 1L] + rho[t]) {
      rho[t + 1L] <- (rho[t - 1L] + rho[t]) / 2
      rho[t + 2L] <- rho[t + 1L]
    }
  }
  S <- chains * n
  tau <- -1 + 2 * sum(rho[seq_len(max_t)]) + rho[max_t + 1L]
  tau <- max(tau, 1 / log10(S))
  S / tau
}

# Rank-normalised split R-hat: max of bulk and folded (tail) versions.
.mip_rhat <- function(sims) {
  sims <- as.matrix(sims)
  if (stats::var(as.vector(sims)) == 0) return(NA_real_)
  sp <- .mip_split_chains(sims)
  bulk <- .mip_rhat_basic(.mip_rank_normalise(sp))
  fold <- abs(sp - stats::median(sp))
  tail <- .mip_rhat_basic(.mip_rank_normalise(fold))
  max(bulk, tail, na.rm = TRUE)
}

.mip_ess_bulk <- function(sims) {
  sims <- as.matrix(sims)
  if (stats::var(as.vector(sims)) == 0) return(NA_real_)
  .mip_ess_raw(.mip_rank_normalise(.mip_split_chains(sims)))
}

# chains: list of per-chain trace matrices (iterations x parameters).
.mip_diagnostics <- function(traces, K, rhat_max = 1.05, ess_min = 400) {
  nm <- .mip_param_names(K)
  diag_idx <- grep("^(Sigma_P|Sigma_E|lambda)", nm)
  out <- data.frame(parameter = nm[diag_idx], rhat = NA_real_,
                    ess_bulk = NA_real_, stringsAsFactors = FALSE)
  for (r in seq_along(diag_idx)) {
    sims <- vapply(traces, function(tr) tr[, diag_idx[r]],
                   numeric(nrow(traces[[1L]])))
    sims <- matrix(sims, ncol = length(traces))
    out$rhat[r] <- .mip_rhat(sims)
    out$ess_bulk[r] <- .mip_ess_bulk(sims)
  }
  attr(out, "converged") <- all(is.finite(out$rhat)) &&
    all(is.finite(out$ess_bulk)) &&
    max(out$rhat) < rhat_max && min(out$ess_bulk) > ess_min
  out
}

# ---------------------------------------------------------------------------
# Starting values
# ---------------------------------------------------------------------------
#
# Chain 1 starts at #187's per-trait REML lambda split of each trait's
# observed variance (.mvn_resolve_lambda(), R/joint_mvn_solver.R, read only),
# with the pairwise-complete trait correlation on Sigma_P. That estimator
# needs the dense n x n R, so above `reml_max_n` tips lambda = 0.9 is used
# instead. Chains 2+ perturb lambda on the logit scale, the variances on the
# log scale, and shrink the correlation towards 0 by a random factor.
.mip_reml_lambda <- function(prob, tree, reml_max_n = 2000L) {
  if (prob$n > reml_max_n) return(NULL)
  R <- stats::cov2cor(ape::vcv(tree))
  R <- R[rownames(prob$Y), rownames(prob$Y)]
  res <- tryCatch(
    .mvn_resolve_lambda("estimate", prob$Y, R, seq_len(prob$K)),
    error = function(e) NULL)
  if (is.null(res)) return(NULL)
  as.numeric(res$lambda_vec)
}

.mip_starts <- function(prob, lambda0, n_chains) {
  K <- prob$K
  s2 <- prob$obs_var
  s2[!is.finite(s2) | s2 <= 0] <- 1
  C <- if (K > 1L) {
    cc <- suppressWarnings(stats::cor(prob$Y, use = "pairwise.complete.obs"))
    cc[!is.finite(cc)] <- 0
    diag(cc) <- 1
    cc <- 0.9 * cc + 0.1 * diag(K)
    ev <- eigen(cc, symmetric = TRUE)
    stats::cov2cor(ev$vectors %*% diag(pmax(ev$values, 0.05), K) %*%
                     t(ev$vectors))
  } else matrix(1, 1, 1)
  make <- function(lam, s2c, Cc) {
    lam <- pmin(pmax(lam, 0.02), 0.98)
    sp <- sqrt(lam * s2c); se <- sqrt((1 - lam) * s2c)
    list(Sigma_P = (sp %o% sp) * Cc,
         Sigma_E = diag(se^2, K))
  }
  lapply(seq_len(n_chains), function(c) {
    if (c == 1L) return(make(lambda0, s2, C))
    lam <- stats::plogis(stats::qlogis(pmin(pmax(lambda0, 0.02), 0.98)) +
                           stats::rnorm(K, 0, 1.5))
    s2c <- s2 * exp(stats::rnorm(K, 0, 0.3))
    u <- stats::runif(1L)
    make(lam, s2c, u * C + (1 - u) * diag(K))
  })
}

# ---------------------------------------------------------------------------
# Control and driver
# ---------------------------------------------------------------------------

.mip_default_control <- function() {
  list(n_chains = 4L, n_iter = NULL, burnin = NULL, thin = NULL,
       keep_draws = 1000L, param_uncertainty = c("full", "none", "both"),
       seed = NULL)
}

.mip_resolve_control <- function(control, m) {
  def <- .mip_default_control()
  if (is.null(control)) control <- list()
  if (!is.list(control)) {
    stop("`posterior_control` must be a list.", call. = FALSE)
  }
  bad <- setdiff(names(control), names(def))
  if (length(bad)) {
    stop("Unknown `posterior_control` element(s): ",
         paste(bad, collapse = ", "), ". Allowed: ",
         paste(names(def), collapse = ", "), ".", call. = FALSE)
  }
  ctl <- utils::modifyList(def, control)
  ctl$param_uncertainty <- match.arg(ctl$param_uncertainty,
                                     c("full", "none", "both"))
  ctl$n_chains <- as.integer(ctl$n_chains)
  if (!is.finite(ctl$n_chains) || ctl$n_chains < 1L) {
    stop("`posterior_control$n_chains` must be a positive integer.",
         call. = FALSE)
  }
  ctl$keep_draws <- as.integer(ctl$keep_draws)
  if (!is.finite(ctl$keep_draws) || ctl$keep_draws < m) {
    stop("`posterior_control$keep_draws` must be at least m (", m, ").",
         call. = FALSE)
  }
  # Defaults, set from the n = 1000 smoke (K = 2, settings of G3): with the
  # collapsed moves, 4 chains x 3,000 kept-phase sweeps gave min bulk ESS
  # 473 in the worst setting (lambda = 0.05 / 0.95, Sigma_P[1,2]) and
  # >= 1,495 at lambda = 1. Defaults: 1,000 burn-in (adaptation) and 5,000
  # kept-phase sweeps per chain; thin so that the chains together keep
  # keep_draws.
  ctl$burnin <- as.integer(ctl$burnin %||% 1000L)
  ctl$n_iter <- as.integer(ctl$n_iter %||% 5000L)
  per_chain <- ceiling(ctl$keep_draws / ctl$n_chains)
  if (is.null(ctl$thin)) {
    ctl$thin <- max(1L, as.integer(ctl$n_iter %/% per_chain))
  }
  ctl$thin <- as.integer(ctl$thin)
  if (ctl$burnin < 0L || ctl$n_iter < 1L || ctl$thin < 1L) {
    stop("`posterior_control`: burnin >= 0, n_iter >= 1 and thin >= 1 are ",
         "required.", call. = FALSE)
  }
  n_kept <- ctl$n_chains * (ctl$n_iter %/% ctl$thin)
  if (n_kept < m) {
    stop("`posterior_control`: n_chains * floor(n_iter / thin) = ", n_kept,
         " kept sweeps, fewer than m = ", m, ".", call. = FALSE)
  }
  if (n_kept < ctl$keep_draws) {
    warning("posterior: only ", n_kept, " kept sweeps (keep_draws = ",
            ctl$keep_draws, "); cell intervals use ", n_kept, " draws.",
            call. = FALSE)
  }
  ctl
}

# Fit the posterior on a latent matrix Y (tip order) and a tree.
# Returns ymis draws (n_mis x kept), kept parameter draws, diagnostics.
.mip_fit <- function(Y, tree, ctl, verbose = FALSE) {
  prob <- .mip_problem(Y, tree)
  K <- prob$K
  hyper <- list(nu_W = K + 1, S_W = diag(K), V_alpha = 1000,
                nu_E = K + 1, S_E = diag(0.01 * prob$obs_var, K))
  if (!is.null(ctl$seed)) set.seed(as.integer(ctl$seed))
  chain_seeds <- sample.int(.Machine$integer.max, ctl$n_chains)
  lam_reml <- .mip_reml_lambda(prob, tree)
  lambda0 <- if (is.null(lam_reml)) rep(0.9, K) else lam_reml
  starts <- .mip_starts(prob, lambda0, ctl$n_chains)
  if (verbose) {
    message(sprintf("posterior: %d chains x (%d burn-in + %d sweeps, thin %d)",
                    ctl$n_chains, ctl$burnin, ctl$n_iter, ctl$thin))
  }
  t0 <- proc.time()[["elapsed"]]
  run1 <- function(c) {
    .mip_run_chain(prob, starts[[c]], hyper, ctl$burnin, ctl$n_iter,
                   ctl$thin, seed = chain_seeds[c])
  }
  # Chains run serially (each seeds itself, so results do not depend on
  # scheduling); parallelise across fits in the caller if needed.
  chains <- lapply(seq_len(ctl$n_chains), run1)
  wall <- proc.time()[["elapsed"]] - t0
  diagnostics <- .mip_diagnostics(lapply(chains, `[[`, "trace"), K)
  SP <- do.call(.mip_abind3, lapply(chains, `[[`, "Sigma_P"))
  SE <- do.call(.mip_abind3, lapply(chains, `[[`, "Sigma_E"))
  mu <- do.call(rbind, lapply(chains, `[[`, "mu"))
  ymis <- do.call(cbind, lapply(chains, `[[`, "ymis"))
  nm <- colnames(Y)
  improper <- NULL
  if (ctl$param_uncertainty %in% c("none", "both")) {
    # Plug-in draws at the posterior-mean covariances. Drawn after the
    # chains, so the RNG stream of the proper results is untouched.
    SP_hat <- apply(SP, c(1L, 2L), mean)
    SE_hat <- apply(SE, c(1L, 2L), mean)
    n_draws <- ncol(ymis)
    fx <- .mip_fixed_draws(prob, SP_hat, SE_hat, n_draws)
    if (identical(ctl$param_uncertainty, "none")) {
      ymis <- fx$ymis
      mu <- fx$mu
      SP <- array(SP_hat, c(K, K, n_draws))
      SE <- array(SE_hat, c(K, K, n_draws))
    } else {
      dimnames(SP_hat) <- dimnames(SE_hat) <- list(nm, nm)
      improper <- list(ymis = fx$ymis, Sigma_P = SP_hat, Sigma_E = SE_hat)
    }
  }
  # Index the diagonals directly: with K = 1, SP[, , d] drops to a scalar
  # and diag(scalar) would build an identity matrix instead.
  lambda <- t(vapply(seq_len(dim(SP)[3L]), function(d) {
    ii <- cbind(seq_len(K), seq_len(K), d)
    SP[ii] / (SP[ii] + SE[ii])
  }, numeric(K)))
  if (K == 1L) lambda <- matrix(lambda, ncol = 1L)
  dimnames(SP) <- dimnames(SE) <- list(nm, nm, NULL)
  colnames(lambda) <- colnames(mu) <- nm
  list(ymis = ymis, miss = prob$miss,
       params = list(Sigma_P = SP, Sigma_E = SE, lambda = lambda, mu = mu),
       diagnostics = diagnostics,
       start = list(lambda_reml = if (is.null(lam_reml)) NULL else
                      stats::setNames(lam_reml, nm),
                    lambda_start = stats::setNames(lambda0, nm)),
       sweeps = ctl$n_chains * (ctl$burnin + ctl$n_iter),
       wall_s = wall, hyper = hyper, improper = improper)
}

.mip_abind3 <- function(...) {
  xs <- list(...)
  K1 <- dim(xs[[1L]])[1L]; K2 <- dim(xs[[1L]])[2L]
  n <- sum(vapply(xs, function(x) dim(x)[3L], integer(1)))
  array(unlist(xs, use.names = FALSE), c(K1, K2, n))
}

# Refusal message for plug-in draws (param_uncertainty = "none"), shared by
# with_imputations() and pool_mi().
.mip_plugin_refusal <- function() {
  paste0("These are plug-in draws from multi_impute(draws_method = ",
         "\"posterior\") with posterior_control$param_uncertainty = \"none\" ",
         "(mi_workflow \"pigauto_posterior_plugin_diagnostic\"): the ",
         "covariance matrices are fixed at their posterior means, so the ",
         "draws leave out parameter uncertainty and are for validation only. ",
         "Downstream inference on them is unsupported. Rerun multi_impute() ",
         "with the default param_uncertainty = \"full\".")
}

# ---------------------------------------------------------------------------
# multi_impute(draws_method = "posterior") entry point
# ---------------------------------------------------------------------------
.multi_impute_posterior <- function(traits, tree, m, species_col, trait_types,
                                    multi_proportion_groups, log_transform,
                                    covariates, verbose, seed, control,
                                    ignored_args = character(0)) {
  if (!is.null(covariates)) {
    stop("draws_method = \"posterior\" does not support `covariates`: the ",
         "posterior imputation model contains only the imputed traits and ",
         "the tree. Drop `covariates`, or include a fully observed ",
         "continuous covariate as an extra trait column.", call. = FALSE)
  }
  if (length(ignored_args)) {
    message("draws_method = \"posterior\" uses the phylogenetic mixed model ",
            "only (no GNN is fitted); ignoring: ",
            paste(ignored_args, collapse = ", "), ".")
  }
  ctl <- .mip_resolve_control(control, m)
  if (is.null(ctl$seed) && !is.null(seed)) ctl$seed <- as.integer(seed)

  data <- preprocess_traits(traits, tree, species_col = species_col,
                            trait_types = trait_types,
                            multi_proportion_groups = multi_proportion_groups,
                            log_transform = log_transform)
  if (isTRUE(data$multi_obs)) {
    # preprocess_traits() sets multi_obs whenever species_col is given, so
    # check for repeated species before saying there are any.
    if (anyDuplicated(as.character(traits[[species_col]]))) {
      stop("draws_method = \"posterior\" supports one observation per ",
           "species only; `traits` has multiple rows for some species. ",
           "Summarise to one row per species first.", call. = FALSE)
    }
    stop("draws_method = \"posterior\" does not support `species_col` ",
         "(multi-observation input). Each species has one row here, so drop ",
         "`species_col` and give the species names as rownames(traits).",
         call. = FALSE)
  }
  trait_map <- data$trait_map
  types <- vapply(trait_map, `[[`, character(1), "type")
  names(types) <- vapply(trait_map, `[[`, character(1), "name")
  bad <- types[types != "continuous"]
  if (length(bad)) {
    stop("draws_method = \"posterior\" supports continuous traits only. ",
         "Non-continuous trait(s): ",
         paste(sprintf("%s (%s)", names(bad), bad), collapse = ", "),
         ". Drop them from `traits` to use posterior MI for the continuous ",
         "traits. pigauto has no analysis-aware MI path for missing ",
         "non-continuous traits: draws_method = \"conformal\" or ",
         "\"mc_dropout\" imputes them, but only as prediction-diagnostic ",
         "completions, which `with_imputations()` and `pool_mi()` refuse.",
         call. = FALSE)
  }

  X <- data$X_scaled
  rownames(X) <- data$species_names
  tree_use <- tree
  if (!identical(sort(tree$tip.label), sort(rownames(X)))) {
    tree_use <- ape::keep.tip(tree, rownames(X))
  }
  X <- X[tree_use$tip.label, , drop = FALSE]
  if (!anyNA(X)) {
    stop("draws_method = \"posterior\": no missing cells to impute; every ",
         "trait cell in `traits` is observed.", call. = FALSE)
  }
  fit <- .mip_fit(X, tree_use, ctl, verbose = verbose)

  miss <- fit$miss
  n_draws <- ncol(fit$ymis)
  idx <- unique(round(seq(1, n_draws, length.out = m)))
  input_row_order <- data$input_row_order
  input_row_order <- input_row_order[match(rownames(X), data$species_names)]

  decode_cells <- function(z, cols) {
    out <- z
    for (k in unique(cols)) {
      tm <- trait_map[[k]]
      sel <- cols == k
      v <- z[sel, , drop = FALSE] * tm$sd + tm$mean
      if (isTRUE(tm$log_transform)) v <- exp(v)
      out[sel, ] <- v
    }
    out
  }
  latent_cols <- vapply(trait_map, function(tm) tm$latent_cols[1L], integer(1))
  cell_trait <- match(miss[, 2L], latent_cols)
  in_row <- input_row_order[miss[, 1L]]

  # Latent y_mis draws (n_mis x draws) -> m completed datasets (original
  # scale), per-cell 95% intervals from all draws, per-cell SD.
  assemble <- function(ymis) {
    decoded <- decode_cells(ymis, cell_trait)
    datasets <- lapply(idx, function(d) {
      Xk <- X
      Xk[miss] <- ymis[, d]
      dec <- .dcb_decode_latent(Xk, trait_map)
      build_completed(traits, dec, species_col = NULL,
                      input_row_order = input_row_order)$completed
    })
    q <- t(apply(decoded, 1L, stats::quantile, probs = c(0.025, 0.5, 0.975),
                 names = FALSE))
    if (nrow(miss) == 1L) q <- matrix(q, nrow = 1L)
    cell_interval <- data.frame(
      row = as.integer(in_row),
      trait = names(types)[cell_trait],
      lower = q[, 1L], upper = q[, 3L], median = q[, 2L],
      stringsAsFactors = FALSE)
    cell_interval <- cell_interval[!is.na(cell_interval$row), , drop = FALSE]
    rownames(cell_interval) <- NULL
    list(datasets = datasets, cell_interval = cell_interval,
         sd_cell = apply(decoded, 1L, stats::sd))
  }

  proper <- assemble(fit$ymis)
  datasets <- proper$datasets
  cell_interval <- proper$cell_interval
  Xbar <- X
  Xbar[miss] <- rowMeans(fit$ymis)
  pooled <- build_completed(traits, .dcb_decode_latent(Xbar, trait_map),
                            species_col = NULL,
                            input_row_order = input_row_order)
  imputed_mask <- pooled$imputed_mask

  se <- matrix(0, nrow(imputed_mask), ncol(imputed_mask),
               dimnames = dimnames(imputed_mask))
  ok <- !is.na(in_row)
  se[cbind(in_row[ok], match(names(types)[cell_trait[ok]], colnames(se)))] <-
    proper$sd_cell[ok]

  posterior_improper <- NULL
  if (!is.null(fit$improper)) {
    imp <- assemble(fit$improper$ymis)
    posterior_improper <- list(datasets = imp$datasets,
                               cell_interval = imp$cell_interval,
                               Sigma_P = fit$improper$Sigma_P,
                               Sigma_E = fit$improper$Sigma_E)
  }

  diagnostics <- fit$diagnostics
  if (verbose || !isTRUE(attr(diagnostics, "converged"))) {
    if (!isTRUE(attr(diagnostics, "converged"))) {
      warning(sprintf(paste0(
        "draws_method = \"posterior\": chains did not meet the convergence ",
        "rule (max R-hat %.3f, needs < 1.05; min bulk ESS %.0f, needs > 400). ",
        "Increase posterior_control$n_iter / burnin."),
        max(diagnostics$rhat, na.rm = TRUE),
        min(diagnostics$ess_bulk, na.rm = TRUE)), call. = FALSE)
    }
  }

  out <- structure(
    list(
      datasets        = datasets,
      m               = length(datasets),
      draws_method    = "posterior",
      pooled_point    = pooled$completed,
      se              = se,
      imputed_mask    = imputed_mask,
      fit             = NULL,
      data            = data,
      tree            = tree,
      species_col     = species_col,
      # param_uncertainty = "none" returns plug-in draws (covariances fixed
      # at their posterior means, validation only); a distinct marker makes
      # with_imputations() and pool_mi() refuse them.
      mi_workflow     = if (identical(ctl$param_uncertainty, "none")) {
        "pigauto_posterior_plugin_diagnostic"
      } else {
        "pigauto_posterior_mi_v1"
      },
      posterior       = list(
        cell_interval = cell_interval,
        diagnostics   = diagnostics,
        params        = fit$params,
        converged     = isTRUE(attr(diagnostics, "converged")),
        control       = ctl,
        hyper         = fit$hyper,
        start         = fit$start,
        draw_index    = idx,
        wall_s        = fit$wall_s,
        sweeps        = fit$sweeps
      )
    ),
    class = c("pigauto_posterior_mi", "pigauto_mi", "list")
  )
  # param_uncertainty = "both" (validation only): plug-in draws from the same
  # chain run, beside (never instead of) the proper results.
  if (!is.null(posterior_improper)) out$posterior_improper <- posterior_improper
  out
}
