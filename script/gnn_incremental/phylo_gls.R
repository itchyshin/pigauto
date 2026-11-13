# Observation-level phylogenetic random-intercept GLS.
# Y_i = X_i beta + u_{species(i)} + epsilon_i
# u ~ N(0, tau2 V_tree); epsilon ~ N(0, sigma2 I).

predict_phylo_gls_conditional <- function(y_obs, X_obs, X_new, V,
                                          species_obs, species_new,
                                          tau2, sigma2) {
  y_obs <- as.numeric(y_obs)
  X_obs <- as.matrix(X_obs)
  X_new <- as.matrix(X_new)
  V <- as.matrix(V)
  species_obs <- as.integer(species_obs)
  species_new <- as.integer(species_new)
  if (nrow(X_obs) != length(y_obs) || nrow(X_new) != length(species_new) ||
      ncol(X_obs) != ncol(X_new) ||
      length(species_obs) != length(y_obs) ||
      nrow(V) != ncol(V) || any(!is.finite(y_obs)) ||
      any(!is.finite(X_obs)) || any(!is.finite(X_new)) ||
      any(species_obs < 1L | species_obs > nrow(V)) ||
      any(species_new < 1L | species_new > nrow(V))) {
    stop("Inputs are dimensionally inconsistent or non-finite.", call. = FALSE)
  }
  if (length(tau2) != 1L || !is.finite(tau2) || tau2 < 0 ||
      length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0) {
    stop("tau2 must be non-negative and sigma2 must be positive.",
         call. = FALSE)
  }
  C_oo <- tau2 * V[species_obs, species_obs, drop = FALSE] +
    diag(sigma2, length(species_obs))
  factor <- tryCatch(chol(C_oo), error = function(e) NULL)
  if (is.null(factor)) stop("Observed covariance is not positive definite.",
                            call. = FALSE)
  solve_oo <- function(b) backsolve(factor, forwardsolve(t(factor), b))
  Cinv_X <- solve_oo(X_obs)
  Cinv_y <- solve_oo(y_obs)
  Xt_Cinv_X <- crossprod(X_obs, Cinv_X)
  beta <- tryCatch(solve(Xt_Cinv_X, crossprod(X_obs, Cinv_y)),
                   error = function(e) NULL)
  if (is.null(beta) || any(!is.finite(beta))) {
    stop("Generalized least-squares fixed effects are not identifiable.",
         call. = FALSE)
  }
  residual <- y_obs - drop(X_obs %*% beta)
  C_no <- tau2 * V[species_new, species_obs, drop = FALSE]
  pred <- drop(X_new %*% beta + C_no %*% solve_oo(residual))
  pred
}

fit_phylo_random_intercept <- function(y, X, species_index, V,
                                       maxit = 100L) {
  started <- proc.time()[["elapsed"]]
  y <- as.numeric(y)
  X <- as.matrix(X)
  V <- as.matrix(V)
  species_index <- as.integer(species_index)
  if (length(y) != nrow(X) || length(y) != length(species_index) ||
      any(!is.finite(y)) || any(!is.finite(X)) ||
      nrow(V) != ncol(V) ||
      any(species_index < 1L | species_index > nrow(V))) {
    stop("GLS inputs are dimensionally inconsistent or non-finite.",
         call. = FALSE)
  }
  if (qr(X)$rank < ncol(X)) {
    stop("GLS fixed-effect design is rank deficient.", call. = FALSE)
  }

  species_seen <- sort(unique(species_index))
  group <- match(species_index, species_seen)
  m <- length(species_seen)
  counts <- tabulate(group, nbins = m)
  V_obs <- V[species_seen, species_seen, drop = FALSE]
  chol_V <- tryCatch(chol(V_obs), error = function(e) NULL)
  if (is.null(chol_V)) stop("Observed tree covariance is not positive definite.",
                            call. = FALSE)
  V_inv <- chol2inv(chol_V)
  logdet_V <- 2 * sum(log(diag(chol_V)))
  n <- length(y)
  log2pi <- log(2 * pi)

  solve_C <- function(B, tau2, sigma2, chol_Q = NULL) {
    B <- as.matrix(B)
    if (is.null(chol_Q)) {
      Q <- V_inv / tau2 + diag(counts / sigma2, m)
      chol_Q <- chol(Q)
    }
    grouped <- rowsum(B, group = group, reorder = TRUE)
    qsolve <- backsolve(chol_Q,
                        forwardsolve(t(chol_Q), grouped))
    B / sigma2 - qsolve[group, , drop = FALSE] / sigma2^2
  }

  objective <- function(log_variances, return_fit = FALSE) {
    tau2 <- exp(log_variances[1L])
    sigma2 <- exp(log_variances[2L])
    ans <- tryCatch({
      Q <- V_inv / tau2 + diag(counts / sigma2, m)
      chol_Q <- chol(Q)
      logdet_C <- n * log(sigma2) + m * log(tau2) + logdet_V +
        2 * sum(log(diag(chol_Q)))
      Cinv_X <- solve_C(X, tau2, sigma2, chol_Q)
      Cinv_y <- solve_C(y, tau2, sigma2, chol_Q)
      Xt_Cinv_X <- crossprod(X, Cinv_X)
      beta <- solve(Xt_Cinv_X, crossprod(X, Cinv_y))
      residual <- y - drop(X %*% beta)
      Cinv_residual <- solve_C(residual, tau2, sigma2, chol_Q)
      nll <- 0.5 * (logdet_C + sum(residual * Cinv_residual) + n * log2pi)
      if (!is.finite(nll)) stop("Non-finite likelihood.")
      if (return_fit) {
        list(nll = nll, beta = drop(beta), residual = residual,
             Cinv_residual = drop(Cinv_residual), tau2 = tau2,
             sigma2 = sigma2, chol_Q = chol_Q)
      } else nll
    }, error = function(e) if (return_fit) NULL else Inf)
    ans
  }

  y_var <- max(stats::var(y), 1e-6)
  lower <- log(c(y_var * 1e-6, y_var * 1e-6))
  upper <- log(c(y_var * 100, y_var * 100))
  initial <- log(c(y_var * 0.5, y_var * 0.5))
  opt <- optim(initial, objective, method = "L-BFGS-B", lower = lower,
                upper = upper, control = list(maxit = as.integer(maxit)))
  fitted <- objective(opt$par, return_fit = TRUE)
  if (is.null(fitted)) stop("GLS optimizer returned an invalid fit.",
                            call. = FALSE)

  # Conditional mean of each species effect given all observed rows.
  zt_cinv_resid <- rowsum(fitted$Cinv_residual, group = group, reorder = TRUE)
  u_hat_seen <- fitted$tau2 * V_obs %*% zt_cinv_resid
  u_hat <- numeric(nrow(V))
  u_hat[species_seen] <- drop(u_hat_seen)
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    beta = fitted$beta,
    tau2 = fitted$tau2,
    sigma2 = fitted$sigma2,
    species_seen = species_seen,
    random_effect = u_hat,
    converged = identical(opt$convergence, 0L),
    convergence = opt$convergence,
    message = opt$message %||% "",
    logLik = -fitted$nll,
    elapsed_seconds = elapsed,
    predict = function(X_new, species_new) {
      X_new <- as.matrix(X_new)
      species_new <- as.integer(species_new)
      cross <- fitted$tau2 * V[species_new, species_seen, drop = FALSE]
      zt_resid <- rowsum(fitted$Cinv_residual, group = group, reorder = TRUE)
      beta_pred <- drop(X_new %*% fitted$beta)
      drop(beta_pred + cross %*% zt_resid)
    }
  )
}
