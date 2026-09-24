# script/rubin_lib.R
#
# Pooling + estimand library for the Rubin's-rules MI-arm comparison lane (arc/rubin-freq-bace).
# S1 owns this file only -- no edits to R/ or BACE/, no changes to the DGP in
# script/campaign_gnn_off_lib.R (read there for make_cell()/sim_latents() conventions).
#
# Defines:
#   rubin_pool()            -- Rubin's rules for a scalar estimate (mice::pool.scalar-compatible)
#   rubin_cell_intervals()  -- per-cell MI prediction intervals for missing values (U = 0)
#   est_pgls_slope()        -- gls(c2 ~ c1, corPagel) slope, its variance, and lambda_hat
#   est_phylo_cor()         -- GLS-whitened (Pagel's-lambda-weighted) correlation of c1, c2
#   pool_cor()               -- Rubin-pool a set of Fisher-z correlations
#
# rubin_pool()'s field names and arithmetic were checked field-for-field against mice 3.19.0's
# mice:::pool.scalar(rule = "rubin1987") and mice:::barnard.rubin() (2026-09-24) so it reproduces
# mice::pool() to numerical precision -- see script/tests-rubin/test-lib-pool.R part (b).

#' Barnard & Rubin (1999) adjusted degrees of freedom.
#' Reproduces mice:::barnard.rubin() exactly: dfcom = Inf collapses to the old Rubin (1987) df,
#' (m - 1) / lambda^2, where lambda here is the Barnard-Rubin "fraction of information missing due
#' to nonresponse" quantity (1 + 1/m) * b / t -- NOT the phylogenetic-signal lambda used elsewhere
#' in this file (est_pgls_slope()/est_phylo_cor()); the two are unrelated quantities that happen to
#' share a name in their respective literatures.
#'
#' @param m integer, number of imputations
#' @param b between-imputation variance
#' @param t total variance
#' @param dfcom complete-data degrees of freedom (Inf -> old Rubin df)
#' @return numeric degrees of freedom
.barnard_rubin_df <- function(m, b, t, dfcom = Inf) {
  lambda <- (1 + 1 / m) * b / t
  dfold <- (m - 1) / lambda^2
  if (is.infinite(dfcom)) return(dfold)
  tmp <- (1 - lambda) * (1 + dfcom) * dfcom
  (m - 1) * tmp / ((dfcom + 3) * (m - 1) + lambda^2 * tmp)
}

#' Rubin's rules pooling for a scalar estimate across M imputations.
#'
#' Matches mice:::pool.scalar(rule = "rubin1987") / mice::pool()$pooled field-for-field:
#' estimate = qbar = mean(q); W = ubar = mean(u); B = var(q); T = W + (1+1/M)*B;
#' riv = (1+1/M)*B/W; lambda = (1+1/M)*B/T (see .barnard_rubin_df() note on the name clash);
#' fmi = (riv + 2/(df+3)) / (riv+1); df via Barnard & Rubin (1999), df_com = Inf reproduces the
#' old Rubin (1987) df.
#'
#' @param q numeric length-M, point estimates
#' @param u numeric length-M, within-imputation variances
#' @param df_com complete-data df (default Inf)
#' @param conf confidence level for the CI (default 0.95)
#' @return list with estimate, W, B, T, riv, lambda, fmi, df, se, lower, upper
rubin_pool <- function(q, u, df_com = Inf, conf = 0.95) {
  stopifnot(length(q) == length(u), length(q) >= 2)
  m <- length(q)
  qbar <- mean(q)
  W <- mean(u)
  B <- stats::var(q)
  Tt <- W + (1 + 1 / m) * B
  riv <- (1 + 1 / m) * B / W
  lambda <- (1 + 1 / m) * B / Tt
  df <- .barnard_rubin_df(m, B, Tt, dfcom = df_com)
  fmi <- (riv + 2 / (df + 3)) / (riv + 1)
  se <- sqrt(Tt)
  alpha <- 1 - conf
  crit <- stats::qt(1 - alpha / 2, df)
  list(estimate = qbar, W = W, B = B, T = Tt, riv = riv, lambda = lambda, fmi = fmi, df = df,
       se = se, lower = qbar - crit * se, upper = qbar + crit * se)
}

#' Per-cell MI prediction intervals for missing values (within-imputation variance = 0).
#'
#' For each of K missing cells, given M imputed draws, returns the between-imputation-only Rubin
#' interval T = (1+1/M)*B with df = M-1 (fixed, per spec -- not Barnard-Rubin-adjusted, since
#' dfcom is undefined for a single missing cell with no model degrees of freedom). A cell with
#' B = 0 (all M draws identical) falls out of the same formula as a zero-width interval; it is
#' additionally flagged via zero_width rather than treated as an error.
#'
#' @param draws numeric M x K matrix, M draws (rows) for each of K missing cells (columns)
#' @param conf confidence level (default 0.95)
#' @return data.frame with one row per cell: cell, mean, B, T, df, lower, upper, zero_width
rubin_cell_intervals <- function(draws, conf = 0.95) {
  draws <- as.matrix(draws)
  m <- nrow(draws); k <- ncol(draws)
  stopifnot(m >= 2)
  mu <- colMeans(draws)
  B <- apply(draws, 2, stats::var)
  Tt <- (1 + 1 / m) * B
  df <- m - 1
  alpha <- 1 - conf
  crit <- stats::qt(1 - alpha / 2, df)
  se <- sqrt(Tt)
  cellnames <- colnames(draws); if (is.null(cellnames)) cellnames <- seq_len(k)
  data.frame(cell = cellnames, mean = mu, B = B, T = Tt, df = df,
             lower = mu - crit * se, upper = mu + crit * se, zero_width = B == 0,
             row.names = NULL)
}

#' PGLS slope of c2 ~ c1 under Pagel's lambda, with retry on convergence failure.
#'
#' nlme::gls(c2 ~ c1, correlation = ape::corPagel(start, phy = tree, form = ~sp, fixed = FALSE),
#' data = df, method = "REML"), where sp = rownames(df) (must match tree$tip.label). ape::corPagel's
#' internal optimiser is unconstrained -- it validates only the starting value, not the fitted one
#' -- so a fit can converge to a lambda outside [0,1] (observed: seed 89 of the G-S1c gate, start
#' 0.5 -> lambda = -0.194). Retries from alternate lambda starting values (0.1, then 0.9) whenever
#' the previous attempt either errors or returns a lambda outside [0,1]; if every start fails,
#' returns NA fields with converged = FALSE and a warning() -- never silently.
#'
#' @param df data.frame with columns c1, c2; rownames = species, matching tree$tip.label
#' @param tree phylo
#' @return list with estimate (slope), variance (its vcov diagonal entry), df_com (= n - 2),
#'   lambda_hat, converged
est_pgls_slope <- function(df, tree) {
  sp <- rownames(df)
  stopifnot(all(sp %in% tree$tip.label))
  d2 <- data.frame(df, sp = sp, stringsAsFactors = FALSE)
  n <- nrow(d2)
  starts <- c(0.5, 0.1, 0.9)
  fit <- NULL; lambda_hat <- NA_real_
  for (s in starts) {
    cand <- tryCatch(
      nlme::gls(c2 ~ c1,
                correlation = ape::corPagel(s, phy = tree, form = ~sp, fixed = FALSE),
                data = d2, method = "REML"),
      error = function(e) NULL)
    if (is.null(cand)) next
    lam <- unname(stats::coef(cand$modelStruct$corStruct, unconstrained = FALSE))
    if (is.finite(lam) && lam >= 0 && lam <= 1) { fit <- cand; lambda_hat <- lam; break }
  }
  # Bounded fallback: when every free start leaves [0,1] (typically the REML optimum sits on or
  # beyond a boundary; first seen at n = 60, seed 7001 of gate G-S4b), profile the REML log-likelihood
  # over lambda in [0,1] with lambda fixed, checking both endpoints, and refit at the maximiser.
  if (is.null(fit)) {
    fit_at <- function(l) tryCatch(
      nlme::gls(c2 ~ c1, correlation = ape::corPagel(l, phy = tree, form = ~sp, fixed = TRUE),
                data = d2, method = "REML"), error = function(e) NULL)
    ll_at <- function(l) { f <- fit_at(l); if (is.null(f)) -Inf else as.numeric(stats::logLik(f)) }
    opt <- stats::optimize(ll_at, c(0, 1), maximum = TRUE, tol = 1e-4)
    cand_l <- c(0, opt$maximum, 1)
    best <- cand_l[which.max(vapply(cand_l, ll_at, numeric(1)))]
    fit <- fit_at(best)
    if (!is.null(fit)) lambda_hat <- best
  }
  if (is.null(fit)) {
    warning("est_pgls_slope: gls() failed to converge to a valid lambda in [0,1] from starts ",
            paste(starts, collapse = ", "), " and the bounded profile")
    return(list(estimate = NA_real_, variance = NA_real_, df_com = n - 2,
                lambda_hat = NA_real_, converged = FALSE))
  }
  v <- stats::vcov(fit)
  list(estimate = unname(stats::coef(fit)["c1"]), variance = unname(v["c1", "c1"]),
       df_com = n - 2, lambda_hat = lambda_hat, converged = TRUE)
}

#' GLS-whitened (Pagel's-lambda-weighted) correlation of c1 and c2.
#'
#' Builds C_lambda = lambda * cov2cor(vcv(tree)) + (1 - lambda) * I, GLS-demeans c1, c2 under
#' C_lambda (weights w = C_lambda^-1 %*% 1), and returns the generalized correlation
#' r = (x-mx)'C^-1(y-my) / sqrt((x-mx)'C^-1(x-mx) * (y-my)'C^-1(y-my)). If lambda is NULL, it is
#' taken from est_pgls_slope(df, tree)$lambda_hat.
#'
#' @param df data.frame with columns c1, c2; rownames = species, matching tree$tip.label
#' @param tree phylo
#' @param lambda numeric in [0,1], or NULL to estimate via est_pgls_slope()
#' @return list with r, z (= atanh(r)), variance (= 1/(n-3)), df_com (= n-3), lambda_used
est_phylo_cor <- function(df, tree, lambda = NULL) {
  if (is.null(lambda)) lambda <- est_pgls_slope(df, tree)$lambda_hat
  sp <- rownames(df); n <- length(sp)
  # A failed slope fit leaves lambda NA; return NA fields rather than erroring in chol() (Meng B1).
  if (!is.finite(lambda)) return(list(r = NA_real_, z = NA_real_, variance = 1 / (n - 3), df_com = n - 3,
                                     lambda_used = NA_real_))
  V <- cov2cor(ape::vcv(tree))[sp, sp]
  C <- lambda * V + (1 - lambda) * diag(n)
  Cinv <- chol2inv(chol(C))
  w <- rowSums(Cinv); s <- sum(Cinv)
  x <- df$c1; y <- df$c2
  mx <- sum(w * x) / s; my <- sum(w * y) / s
  xc <- x - mx; yc <- y - my
  num <- as.numeric(t(xc) %*% Cinv %*% yc)
  den <- sqrt(as.numeric(t(xc) %*% Cinv %*% xc) * as.numeric(t(yc) %*% Cinv %*% yc))
  r <- num / den
  list(r = r, z = atanh(r), variance = 1 / (n - 3), df_com = n - 3, lambda_used = lambda)
}

#' Rubin-pool a set of Fisher-z-transformed correlations.
#'
#' u = 1/(n-3) per draw. Pools on the z scale via rubin_pool(), then back-transforms the point
#' estimate and CI via tanh().
#'
#' @param z_vec numeric length-M, Fisher-z correlations (e.g. from est_phylo_cor()$z per
#'   imputation)
#' @param n sample size used for each z (scalar, or length-M if it varies by draw)
#' @param conf confidence level (default 0.95)
#' @return list with r (pooled correlation), lower, upper (back-transformed CI), and z_pool (the
#'   full rubin_pool() result on the z scale)
pool_cor <- function(z_vec, n, conf = 0.95) {
  m <- length(z_vec)
  n <- rep_len(n, m)
  u <- 1 / (n - 3)
  # Fisher's z has a normal complete-data reference, so df_com = Inf (mice's default); a t_{n-3}
  # reference here disagreed with the qnorm complete-data check (Meng review N2).
  zp <- rubin_pool(z_vec, u, df_com = Inf, conf = conf)
  list(r = tanh(zp$estimate), lower = tanh(zp$lower), upper = tanh(zp$upper), z_pool = zp)
}
