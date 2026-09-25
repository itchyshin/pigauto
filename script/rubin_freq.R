# script/rubin_freq.R
# Frequentist MI draw code for the continuous-family block (S2 lane). Reuses run_freq()'s
# (script/campaign_gnn_off_lib.R) block-fit call and transforms -- does not edit that file.
#
# run_freq() fits Rphylopars::phylopars(df_in, tree, model = phylo_model, phylo_correlated = TRUE,
# pheno_correlated = TRUE, REML = TRUE) jointly on the continuous-family columns, with proportion
# traits pre-transformed via qlogis (post-transformed via plogis) and, under the default
# count_method = "phyloglm_poisson_gee", count traits EXCLUDED from the joint block (they get their
# own marginal Poisson GEE model in run_freq -- not reproduced here). default_block_traits() below
# reproduces that exclusion.
#
# CHECK (2026-09-24, Rphylopars 0.3.10): for the campaign's data (one observation per species),
# Rphylopars' `pheno_error` argument defaults to FALSE ("unless <=1 observation per species is
# provided" -- ?phylopars), so `fit$pars` contains ONLY `phylocov`, no `phenocov`. There is no
# separately-estimated within-species/residual covariance to extract; Sigma_e is NULL in that
# regime. joint_cov() below still accepts a non-NULL Sigma_e (kron(Sigma_e, I_n) term) for
# completeness / multi-obs future use, but it will be NULL for every arm this lane actually runs.

# ---- block definition (mirrors run_freq's default continuous-family block) -------------------

default_block_traits <- function(cell) {
  ct <- cell$cont_traits
  is_cnt <- vapply(cell$truth[ct], is.integer, logical(1))
  setdiff(ct, ct[is_cnt])
}

transform_block <- function(df_raw, block_traits, is_prp) {
  Y <- df_raw[, block_traits, drop = FALSE]
  for (i in seq_along(block_traits)) {
    Y[[i]] <- if (is_prp[i]) stats::qlogis(as.numeric(Y[[i]])) else as.numeric(Y[[i]])
  }
  as.matrix(Y)
}

# ---- 1. fit_block ------------------------------------------------------------------------------
# Same call/transforms as run_freq's joint continuous-family fit (campaign_gnn_off_lib.R L446-465).
fit_block <- function(df_raw, tree, block_traits, trait_types = NULL, phylo_model = "lambda") {
  is_prp <- block_traits %in% names(trait_types)[trait_types == "proportion"]
  sp <- tree$tip.label
  Y <- transform_block(df_raw[sp, , drop = FALSE], block_traits, is_prp)
  df_in <- data.frame(species = sp, Y, stringsAsFactors = FALSE)
  fit <- Rphylopars::phylopars(df_in, tree = tree, model = phylo_model,
                                phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE)
  lambda <- if (!is.null(fit$model$lambda)) fit$model$lambda else NA_real_
  list(fit = fit, mu = fit$mu[block_traits], Sigma_p = fit$pars$phylocov,
       Sigma_e = fit$pars$phenocov,  # NULL for single-obs data (pheno_error defaults FALSE)
       lambda = lambda, block_traits = block_traits, is_prp = is_prp, species = sp,
       phylo_model = phylo_model)
}

# ---- 2. joint_cov -------------------------------------------------------------------------------
# Cov(vec(Y)), Y an n(species) x p(trait) matrix, vec = R's column-major as.vector (species fastest
# within each trait block) -- matches kronecker(Sigma_p, C_lambda). C_lambda verified (test-freq-
# condmean.R) against Rphylopars' own anc_recon: for a height-1 ultrametric tree,
# C_lambda = lambda * C + (1 - lambda) * diag(diag(C)), C = ape::vcv(tree).
joint_cov <- function(pars, tree) {
  sp <- pars$species
  n <- length(sp); p <- length(pars$block_traits)
  C <- ape::vcv(tree)[sp, sp]
  Clam <- pars$lambda * C + (1 - pars$lambda) * diag(diag(C))
  V <- kronecker(pars$Sigma_p, Clam)
  if (!is.null(pars$Sigma_e)) V <- V + kronecker(pars$Sigma_e, diag(n))
  mu_vec <- rep(pars$mu, each = n)
  list(V = V, mu_vec = mu_vec, species = sp, block_traits = pars$block_traits, n = n, p = p)
}

# ---- 3. cond_draw --------------------------------------------------------------------------------
# Y: n x p matrix (rows named by species, cols = pars$block_traits) on the FITTED/TRANSFORMED scale
# (i.e. already qlogis'd for proportion columns), NA for missing cells. Draws ALL missing cells
# jointly from N(cond_mean, cond_cov) via one Cholesky of the conditional covariance; observed cells
# pass through unchanged. Returns draws on the transformed scale plus a backtransform() helper
# (plogis for proportion columns, identity otherwise) matching run_freq's itf().
cond_draw <- function(Y, pars, tree, M) {
  jc <- joint_cov(pars, tree)
  sp <- jc$species; block_traits <- jc$block_traits; n <- jc$n; p <- jc$p
  Y <- as.matrix(Y[sp, block_traits, drop = FALSE])
  y_vec <- as.vector(Y)
  obs_idx <- which(!is.na(y_vec)); mis_idx <- which(is.na(y_vec))
  V <- jc$V; mu_vec <- jc$mu_vec

  backtransform <- function(mat) {
    out <- mat
    for (i in seq_along(block_traits)) if (pars$is_prp[i]) out[, i] <- stats::plogis(out[, i])
    out
  }

  if (!length(mis_idx)) {
    mat0 <- matrix(y_vec, n, p, dimnames = list(sp, block_traits))
    return(list(draws = replicate(M, mat0, simplify = FALSE), backtransform = backtransform,
                cond_mean = numeric(0), cond_cov = matrix(numeric(0), 0, 0),
                mis_idx = mis_idx, obs_idx = obs_idx))
  }

  V_oo <- V[obs_idx, obs_idx, drop = FALSE]
  V_mo <- V[mis_idx, obs_idx, drop = FALSE]
  V_mm <- V[mis_idx, mis_idx, drop = FALSE]
  mu_o <- mu_vec[obs_idx]; mu_m <- mu_vec[mis_idx]
  y_o <- y_vec[obs_idx]

  # Cholesky solves, then a hard validity check (2026-09-24 campaign): about 1 in 300 bootstrap refits lands
  # on a degenerate parameter set (lambda* = 1 with a near-singular phylogenetic covariance), and the old
  # solve() + silent ginv() fallback then returned conditional draws of 1e6 to 1e67. A conditional variance can
  # never exceed the marginal variance, and a conditional mean cannot sit 20 marginal SDs from the mean; a
  # violation (or a failed Cholesky) is a numerical failure, reported as an error the caller handles.
  R_oo <- tryCatch(chol(V_oo), error = function(e) NULL)
  if (is.null(R_oo)) stop("degenerate conditional distribution: chol(V_oo) failed")
  a <- backsolve(R_oo, forwardsolve(t(R_oo), y_o - mu_o))
  W <- backsolve(R_oo, forwardsolve(t(R_oo), t(V_mo)))
  cond_mean <- as.vector(mu_m + V_mo %*% a)
  cond_cov <- V_mm - V_mo %*% W
  cond_cov <- (cond_cov + t(cond_cov)) / 2
  marg_var <- diag(V_mm); cvar <- diag(cond_cov)
  if (any(!is.finite(cond_mean)) || any(!is.finite(cvar)) || any(cvar > 1.01 * marg_var + 1e-8) ||
      any(cvar < -1e-8 * max(marg_var)) || any(abs(cond_mean - mu_m) > 20 * sqrt(marg_var)))
    stop("degenerate conditional distribution: moments outside their bounds")
  L <- tryCatch(chol(cond_cov), error = function(e) NULL)
  if (is.null(L)) {                 # PSD projection: clip tiny negative eigenvalues from rounding
    e <- eigen(cond_cov, symmetric = TRUE)
    L <- t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values)))
  }

  n_mis <- length(mis_idx)
  draws <- vector("list", M)
  for (m in seq_len(M)) {
    z <- stats::rnorm(n_mis)
    y_full <- y_vec
    y_full[mis_idx] <- cond_mean + as.vector(crossprod(L, z))
    draws[[m]] <- matrix(y_full, n, p, dimnames = list(sp, block_traits))
  }
  list(draws = draws, backtransform = backtransform, cond_mean = cond_mean, cond_cov = cond_cov,
       mis_idx = mis_idx, obs_idx = obs_idx)
}

# apply_block_draw(): writes a backtransformed draw into ONLY the originally-missing block cells,
# leaving observed cells untouched (bit-identical to df_miss) -- matches run_freq's own convention
# (`comp[[v]] <- df_miss[[v]]; comp[mask[, v], v] <- itf(...)`), and avoids the ~1e-16 qlogis/plogis
# round-trip drift a full-column overwrite would introduce on observed proportion cells.
apply_block_draw <- function(df_miss, block_traits, raw) {
  out <- df_miss
  for (v in block_traits) {
    idx <- is.na(df_miss[[v]])
    if (any(idx)) out[[v]][idx] <- raw[rownames(df_miss)[idx], v]
  }
  out
}

# ---- 4. mi_freq_B ---------------------------------------------------------------------------
# One fit, M joint conditional draws (parameters fixed). Non-block columns are left as in df_miss.
mi_freq_B <- function(cell, M, block_traits = NULL, phylo_model = "lambda") {
  if (is.null(block_traits)) block_traits <- default_block_traits(cell)
  df_miss <- cell$df_miss; tree <- cell$tree; trait_types <- cell$trait_types

  pars <- fit_block(df_miss, tree, block_traits, trait_types, phylo_model)
  Yt <- transform_block(df_miss, block_traits, pars$is_prp)
  cd <- cond_draw(Yt, pars, tree, M)

  datasets <- lapply(cd$draws, function(mat) {
    raw <- cd$backtransform(mat)
    apply_block_draw(df_miss, block_traits, raw)
  })
  list(datasets = datasets, pars = pars)
}

# ---- 5. mi_freq_A -----------------------------------------------------------------------------
# Parametric bootstrap: fit once (pars0), then for m = 1..M simulate a full Y* ~ N(mu0, V0), apply
# the ORIGINAL missingness pattern, refit fit_block() on Y* (theta*_m), and draw the ORIGINAL data's
# missing cells jointly from the conditional normal under theta*_m. A failed refit (n_fail) or a
# degenerate conditional distribution (n_degenerate, see cond_draw) redraws the bootstrap sample, up to 5
# attempts per m; if all fail the draw m is NULL and the runner drops and counts it -- never silently.
mi_freq_A <- function(cell, M, block_traits = NULL, phylo_model = "lambda") {
  if (is.null(block_traits)) block_traits <- default_block_traits(cell)
  df_miss <- cell$df_miss; tree <- cell$tree; trait_types <- cell$trait_types; mask <- cell$mask

  pars0 <- fit_block(df_miss, tree, block_traits, trait_types, phylo_model)
  jc0 <- joint_cov(pars0, tree)
  sp <- jc0$species; n <- jc0$n; p <- jc0$p
  L0 <- chol(jc0$V)
  mask_block <- mask[sp, block_traits, drop = FALSE]

  simulate_boot <- function() {
    z <- stats::rnorm(n * p)
    y_star <- jc0$mu_vec + as.vector(crossprod(L0, z))
    Ystar <- matrix(y_star, n, p, dimnames = list(sp, block_traits))
    for (i in seq_along(block_traits)) if (pars0$is_prp[i]) Ystar[, i] <- stats::plogis(Ystar[, i])
    Ystar[mask_block] <- NA
    df_star <- df_miss
    df_star[sp, block_traits] <- Ystar[sp, block_traits]
    df_star
  }
  refit <- function(df_star) tryCatch(fit_block(df_star, tree, block_traits, trait_types, phylo_model),
                                       error = function(e) NULL)

  datasets <- vector("list", M)
  pars_star <- vector("list", M)
  n_fail <- 0L; n_degenerate <- 0L; max_attempts <- 5L

  for (m in seq_len(M)) {
    for (attempt in seq_len(max_attempts)) {
      fit_star <- refit(simulate_boot())
      if (is.null(fit_star)) { n_fail <- n_fail + 1L; next }
      Yt_orig <- transform_block(df_miss, block_traits, fit_star$is_prp)
      cd <- tryCatch(cond_draw(Yt_orig, fit_star, tree, 1L), error = function(e) NULL)
      if (is.null(cd)) { n_degenerate <- n_degenerate + 1L; next }   # redraw the bootstrap sample
      raw <- cd$backtransform(cd$draws[[1L]])
      datasets[[m]] <- apply_block_draw(df_miss, block_traits, raw)
      pars_star[[m]] <- list(lambda = fit_star$lambda, Sigma_p = fit_star$Sigma_p, Sigma_e = fit_star$Sigma_e)
      break
    }
  }
  list(datasets = datasets, pars_star = pars_star, n_fail = n_fail, n_degenerate = n_degenerate)
}
