#!/usr/bin/env Rscript
# script/bench_beta_sweep.R
#
# Beta-strength sweep at lambda=0.20, n=500.
# Sweep beta in {0.3, 0.5, 0.7, 1.0, 1.5} for nonlinear and interactive.
# 30 cells, ~1 hr.

suppressPackageStartupMessages({
  library(ape); library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
})
options(warn = -1L)

HERE <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script", "bench_beta_sweep.rds")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0), ..., "\n", sep = "")

N_SPECIES <- 500L; OBS_PER <- 5L; N_COVS <- 10L
LAMBDA <- 0.20
BETAS <- c(0.3, 0.5, 0.7, 1.0, 1.5)
F_TYPES <- c("nonlinear", "interactive")
N_REPS <- 3L; EPOCHS <- 250L
SP_MISSING_FRAC <- 0.5; WITHIN_MISS_FRAC <- 0.2; SIGMA_RES <- 0.5

log_line("BETA SWEEP: lambda=0.20, n=500, ncov=10, betas=", paste(BETAS, collapse=","))

sim_dgp <- function(n_species, obs_per, n_covs, lambda, f_type, beta_strength, seed) {
  set.seed(seed)
  tree <- ape::rtree(n_species); sp <- tree$tip.label
  phylo_vals <- if (lambda > 0) {
    as.numeric(ape::rTraitCont(tree, model = "BM", sigma = sqrt(lambda), root.value = 0)[sp])
  } else rep(0, n_species)
  names(phylo_vals) <- sp
  n_per_sp <- pmin(pmax(rpois(n_species, obs_per), 1L), 15L)
  species_vec <- rep(sp, n_per_sp); n_total <- length(species_vec)
  Z <- matrix(rnorm(n_total * n_covs), nrow = n_total, ncol = n_covs)
  colnames(Z) <- paste0("z", seq_len(n_covs))
  contrib <- switch(f_type,
    nonlinear = sin(Z[,1])*exp(0.3*Z[,2]) + 0.5*tanh(Z[,3]) +
                  0.3*sin(2*Z[,4]) - 0.2*cos(Z[,5]) +
                  0.4*tanh(Z[,6]+Z[,7]) - 0.3*sin(Z[,8]),
    interactive = Z[,1]*Z[,2] + 0.5*Z[,3]^2 + 0.3*Z[,1]*Z[,4] -
                    0.2*Z[,2]*Z[,5] + 0.25*Z[,6]*Z[,7] - 0.15*Z[,8]*Z[,9]^2
  )
  contrib <- as.numeric(scale(contrib))
  delta <- beta_strength * contrib
  epsilon <- rnorm(n_total, 0, SIGMA_RES)
  y <- phylo_vals[species_vec] + delta + epsilon
  miss_sp <- sample(sp, floor(n_species * SP_MISSING_FRAC))
  y_obs <- y; y_obs[species_vec %in% miss_sp] <- NA
  obs_idx <- which(!is.na(y_obs))
  n_within <- floor(WITHIN_MISS_FRAC * length(obs_idx))
  if (n_within > 0L) y_obs[sample(obs_idx, n_within)] <- NA
  mask <- is.na(y_obs)
  df_obs <- data.frame(species = species_vec, y = y_obs, Z, stringsAsFactors = FALSE)
  df_full <- data.frame(species = species_vec, y = y, Z, stringsAsFactors = FALSE)
  list(tree = tree, df_complete = df_full, df_observed = df_obs,
        mask = mask, predictor_names = colnames(Z))
}

score <- function(pred, truth, mask) {
  ok <- mask & is.finite(truth) & is.finite(pred)
  if (sum(ok) < 5L) return(c(rmse = NA_real_, r = NA_real_))
  c(rmse = sqrt(mean((pred[ok] - truth[ok])^2)),
    r = if (stats::sd(pred[ok]) > 1e-10) suppressWarnings(stats::cor(pred[ok], truth[ok])) else NA_real_)
}

m_lm_nonlinear <- function(d) {
  ddf <- d$df_observed; pred_n <- d$predictor_names
  fit_df <- ddf[!is.na(ddf$y), c("y", pred_n), drop = FALSE]
  poly_terms <- sapply(pred_n, function(p) sprintf("poly(%s, 2, raw=TRUE)", p))
  inter_terms <- combn(pred_n, 2L, function(p) sprintf("%s:%s", p[1], p[2]))
  fmla <- stats::as.formula(paste("y ~", paste(c(poly_terms, inter_terms), collapse=" + ")))
  fit <- tryCatch(stats::lm(fmla, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))
  tryCatch(predict(fit, newdata = ddf), error = function(e) rep(NA_real_, nrow(ddf)))
}

m_phylolm_blup <- function(d) {
  ddf <- d$df_observed; pred_n <- d$predictor_names
  sp_means_y <- tapply(ddf$y, ddf$species,
                        function(v) { v <- v[!is.na(v)]; if (length(v)) mean(v) else NA_real_ })
  sp_means_X <- as.data.frame(lapply(pred_n, function(p)
    tapply(ddf[[p]], ddf$species, mean, na.rm=TRUE)))
  names(sp_means_X) <- pred_n
  sp_df <- data.frame(species = names(sp_means_y), y = as.numeric(sp_means_y),
                       sp_means_X, stringsAsFactors = FALSE)
  sp_df <- sp_df[sp_df$species %in% d$tree$tip.label, , drop = FALSE]
  sp_df_obs <- sp_df[!is.na(sp_df$y), , drop = FALSE]
  if (nrow(sp_df_obs) < length(pred_n) + 5L) return(rep(NA_real_, nrow(ddf)))
  rownames(sp_df_obs) <- sp_df_obs$species
  tree_obs <- ape::keep.tip(d$tree, sp_df_obs$species)
  fmla <- stats::as.formula(paste("y ~", paste(pred_n, collapse=" + ")))
  fit <- tryCatch(phylolm::phylolm(fmla, data = sp_df_obs, phy = tree_obs, model = "lambda"),
                   error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))
  beta_hat <- coef(fit)
  Xmat <- model.matrix(stats::as.formula(paste("~", paste(pred_n, collapse=" + "))), sp_df)
  fixed <- as.numeric(Xmat %*% beta_hat)
  R <- stats::cov2cor(ape::vcv(d$tree)); R <- R[sp_df$species, sp_df$species]
  oi <- which(!is.na(sp_df$y)); mi <- which(is.na(sp_df$y))
  if (length(mi) == 0L) sp_pred <- sp_df$y else {
    e_obs <- sp_df$y[oi] - fixed[oi]
    blup <- as.numeric(R[mi, oi, drop = FALSE] %*%
                          solve(R[oi, oi, drop = FALSE] + diag(1e-6, length(oi)), e_obs))
    sp_pred <- numeric(nrow(sp_df))
    sp_pred[oi] <- sp_df$y[oi]; sp_pred[mi] <- fixed[mi] + blup
  }
  names(sp_pred) <- sp_df$species
  sp_pred[as.character(ddf$species)]
}

m_pigauto <- function(d, seed) {
  ddf <- d$df_observed; pred_n <- d$predictor_names
  cov_df <- as.data.frame(lapply(pred_n, function(p) as.numeric(ddf[[p]])))
  names(cov_df) <- pred_n
  res <- pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                          tree = d$tree, species_col = "species", covariates = cov_df,
                          missing_frac = 0.0, epochs = EPOCHS, verbose = FALSE,
                          seed = seed, safety_floor = TRUE)
  res$completed$y
}

cells <- expand.grid(beta = BETAS, f_type = F_TYPES, rep_id = seq_len(N_REPS),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
total <- nrow(cells); log_line("Total cells: ", total)

results <- list()
for (i in seq_len(total)) {
  row <- cells[i, ]; seed <- 27000L + i
  log_line(sprintf("Cell %d/%d: f=%s, beta=%.1f, rep=%d", i, total, row$f_type, row$beta, row$rep_id))
  d <- tryCatch(sim_dgp(N_SPECIES, OBS_PER, N_COVS, LAMBDA, row$f_type, row$beta, seed),
                 error = function(e) { log_line("DGP ERR: ", e$message); NULL })
  if (is.null(d)) next
  truth <- d$df_complete$y
  for (m_name in c("column_mean", "lm_nonlinear", "phylolm_lambda_blup", "pigauto_sfT")) {
    fn <- switch(m_name,
      column_mean = function() rep(mean(d$df_observed$y, na.rm=TRUE), nrow(d$df_observed)),
      lm_nonlinear = function() m_lm_nonlinear(d),
      phylolm_lambda_blup = function() m_phylolm_blup(d),
      pigauto_sfT = function() m_pigauto(d, seed))
    pred <- tryCatch(fn(), error = function(e) { log_line("  ", m_name, " ERR: ", e$message); rep(NA_real_, nrow(d$df_observed)) })
    s <- score(pred, truth, d$mask)
    results[[length(results) + 1L]] <- data.frame(
      lambda = LAMBDA, beta = row$beta, f_type = row$f_type, rep = row$rep_id,
      method = m_name, rmse = s["rmse"], r = s["r"], stringsAsFactors = FALSE)
  }
  if (i %% 4L == 0L || i == total) {
    df <- do.call(rbind, results); saveRDS(df, out_rds)
    log_line(sprintf("  [checkpoint] %d rows", nrow(df)))
  }
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results); saveRDS(df, out_rds)
log_line(sprintf("=== DONE: %d cells, %d rows ===", total, nrow(df)))
