#!/usr/bin/env Rscript
# script/bench_multi_obs_nonlinear.R
#
# THE bench that should showcase pigauto's GNN/transformer value-add.
#
# Design rationale (2026-04-27):
# Earlier benches (bench_sim_bace_pigauto, bench_sim_bace_nonlinear,
# bench_gnn_capacity) used BACE-sim's DGP, where predictors carry their
# own phylogenetic signal.  That makes the "nonlinear" addition partly
# recoverable by phylolm-lambda's linear+phylo regression -- pigauto's
# GNN had no clean room to add value.  Result: pigauto matched or
# slightly lost to phylolm-lambda on those benches.
#
# Conversely, bench_multi_obs.R (the bundled CTmax-style bench) uses
# acclim_temp ~ N(20, 5) i.i.d. per observation -- NO phylogenetic
# structure on the predictor.  In that DGP pigauto cleanly beats
# species_mean by 28-36% at strong covariate signal.
#
# This bench extends bench_multi_obs.R with the architectural
# ingredients we want to test:
#   - many predictors (5-10) all i.i.d. per observation
#   - nonlinear response (sin/exp/interaction in covariates)
#   - low phylogenetic signal (lambda in {0.1, 0.3}) -- minimises the
#     phylo channel so the covariate channel dominates
#   - multi-obs (5 obs per species avg)
#
# DGP:
#   For each species s in tree (n_species = 300):
#     phylo_s ~ BM(tree, sigma = sqrt(lambda))
#   For each observation i of species s (5 obs each on average):
#     z_i ~ N(0, 1) (n_covs covariates, all i.i.d. across observations)
#     contrib_i = beta_strength * f(z_i)   where f depends on f_type:
#       linear      : f = 0.5*z[1] + 0.3*z[2] + 0.2*z[3] + (linear sum)
#       nonlinear   : f = sin(z[1]) * exp(0.3*z[2]) + 0.5*tanh(z[3])
#       interactive : f = z[1]*z[2] + 0.5*z[3]^2 + 0.3*z[1]*z[4]
#     y_i = phylo_s + contrib_i + epsilon_i,  epsilon_i ~ N(0, 0.5)
#
# Mask:
#   - 50% of species: ALL observations missing (the "no-data species")
#   - For observed species: 20% MCAR within-species observation mask
#
# Methods (7):
#   1. column_mean              -- loss floor
#   2. species_mean             -- per-species observed mean (grand mean for unobserved)
#   3. lm                       -- lm(y ~ z1+z2+...+zk), no phylo
#   4. lm_nonlinear             -- lm(y ~ poly(z, 2) + z1:z2 + z1:z3 ...), no phylo
#   5. phylolm_lambda_blup      -- linear cov + lambda-BM-BLUP, the phylo-aware baseline
#   6. pigauto_sfT              -- pigauto with safety_floor TRUE (current default)
#   7. pigauto_sfF              -- pigauto with safety_floor FALSE (raw GNN)
#
# Decisive predictions (registered before harvest):
#   linear:      pigauto >= phylolm-lambda (linear DGP, both should match;
#                                            pigauto wins or ties)
#   nonlinear:   pigauto BEATS phylolm-lambda by >= 15% at beta=1.0
#                (the GNN should capture sin/exp; phylolm cannot)
#   interactive: pigauto BEATS phylolm-lambda by >= 10% at beta=1.0
#                (the GNN should capture interactions; lm_nonlinear comes
#                 close but cannot exploit phylo)
#
# Tier (env var PIGAUTO_TIER):
#   smoke:  18 cells = 2 lambda * 3 f_type * 3 beta * 1 reps   ~25-40 min
#   medium: 144 cells = 2 lambda * 3 f_type * 3 beta * 4 reps * 2 n_covs ~5-8 hr

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
})
options(warn = -1L)

TIER  <- toupper(Sys.getenv("PIGAUTO_TIER", "smoke"))
HERE  <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script",
                      sprintf("bench_multi_obs_nonlinear_%s.rds", tolower(TIER)))
out_log <- file.path(HERE, "script",
                      sprintf("bench_multi_obs_nonlinear_%s.log", tolower(TIER)))

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

if (TIER == "SMOKE") {
  N_SPECIES <- 300L
  OBS_PER   <- 5L
  N_COVS    <- 5L
  LAMBDAS   <- c(0.1, 0.3)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  BETAS     <- c(0.0, 0.5, 1.0)
  N_REPS    <- 1L
  EPOCHS    <- 200L
} else if (TIER == "MEDIUM") {
  N_SPECIES <- 300L
  OBS_PER   <- 5L
  N_COVS    <- c(5L, 10L)
  LAMBDAS   <- c(0.1, 0.3)
  F_TYPES   <- c("linear", "nonlinear", "interactive")
  BETAS     <- c(0.0, 0.5, 1.0)
  N_REPS    <- 4L
  EPOCHS    <- 300L
} else {
  stop("PIGAUTO_TIER must be smoke or medium")
}

EPOCHS <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", EPOCHS))
SP_MISSING_FRAC <- 0.5    # half species fully missing
WITHIN_MISS_FRAC <- 0.2   # 20% MCAR within observed species
SIGMA_RES <- 0.5

log_line("=========================================================")
log_line("Tier: ", TIER, "   (MULTI-OBS NONLINEAR -- low phylo, many covs, i.i.d.)")
log_line("Sweep: n=", N_SPECIES, " obs/sp=", OBS_PER,
          " | n_covs=", paste(N_COVS, collapse=","),
          " | lambda=", paste(LAMBDAS, collapse=","),
          " | f=", paste(F_TYPES, collapse=","),
          " | beta=", paste(BETAS, collapse=","),
          " | reps=", N_REPS,
          " | epochs=", EPOCHS)
log_line("Output: ", out_rds)

# ---------------- DGP ----------------------------------------------------
sim_multi_nonlinear_dgp <- function(n_species, obs_per, n_covs,
                                      lambda, f_type, beta_strength,
                                      seed,
                                      sp_miss_frac = SP_MISSING_FRAC,
                                      within_miss = WITHIN_MISS_FRAC,
                                      sigma_res = SIGMA_RES) {
  set.seed(seed)
  tree <- ape::rtree(n_species)
  sp <- tree$tip.label
  # Species-level phylogenetic effect
  phylo_vals <- ape::rTraitCont(tree, model = "BM",
                                  sigma = sqrt(lambda), root.value = 0)
  phylo_vals <- as.numeric(phylo_vals[sp])
  names(phylo_vals) <- sp
  # Variable obs per species (Poisson around obs_per, clamped 1..15)
  n_per_sp <- pmin(pmax(rpois(n_species, obs_per), 1L), 15L)
  names(n_per_sp) <- sp
  species_vec <- rep(sp, n_per_sp)
  n_total <- length(species_vec)
  # Generate i.i.d. covariates (no phylo structure)
  Z <- matrix(rnorm(n_total * n_covs), nrow = n_total, ncol = n_covs)
  colnames(Z) <- paste0("z", seq_len(n_covs))
  # Compute contribution
  contrib <- switch(f_type,
    linear = {
      coefs <- c(0.5, 0.3, 0.2, 0.1, -0.1, 0.15, -0.2, 0.05, 0.1, -0.15)[seq_len(n_covs)]
      as.numeric(Z %*% coefs)
    },
    nonlinear = {
      base <- sin(Z[, 1]) * exp(0.3 * Z[, 2]) + 0.5 * tanh(Z[, 3])
      if (n_covs >= 5L) base <- base + 0.3 * sin(2 * Z[, 4]) - 0.2 * cos(Z[, 5])
      base
    },
    interactive = {
      base <- Z[, 1] * Z[, 2] + 0.5 * Z[, 3]^2
      if (n_covs >= 4L) base <- base + 0.3 * Z[, 1] * Z[, 4]
      if (n_covs >= 5L) base <- base - 0.2 * Z[, 2] * Z[, 5]
      base
    }
  )
  contrib <- as.numeric(scale(contrib))   # standardise
  delta <- beta_strength * contrib
  # Build response
  epsilon <- rnorm(n_total, 0, sigma_res)
  y <- phylo_vals[species_vec] + delta + epsilon
  # Apply mask
  miss_sp <- sample(sp, floor(n_species * sp_miss_frac))
  y_obs <- y
  y_obs[species_vec %in% miss_sp] <- NA
  obs_idx <- which(!is.na(y_obs))
  n_within_miss <- floor(within_miss * length(obs_idx))
  if (n_within_miss > 0L) {
    within_miss_idx <- sample(obs_idx, n_within_miss)
    y_obs[within_miss_idx] <- NA
  }
  mask <- is.na(y_obs)

  df_observed <- data.frame(species = species_vec, y = y_obs, Z,
                              stringsAsFactors = FALSE)
  df_complete <- data.frame(species = species_vec, y = y, Z,
                              stringsAsFactors = FALSE)

  list(tree = tree, df_complete = df_complete, df_observed = df_observed,
        mask = mask, predictor_names = colnames(Z),
        meta = list(n_species = n_species, obs_per = obs_per,
                     n_covs = n_covs, lambda = lambda,
                     f_type = f_type, beta_strength = beta_strength,
                     missing_sp = miss_sp))
}

# ---------------- methods -------------------------------------------------
score <- function(pred, truth, mask) {
  ok <- mask & is.finite(truth) & is.finite(pred)
  if (sum(ok) < 5L) return(c(rmse = NA_real_, r = NA_real_))
  p <- pred[ok]; t <- truth[ok]
  rmse <- sqrt(mean((p - t)^2))
  r <- if (stats::sd(p) > 1e-10) suppressWarnings(stats::cor(p, t)) else NA_real_
  c(rmse = rmse, r = r)
}

method_column_mean <- function(d) rep(mean(d$df_observed$y, na.rm=TRUE),
                                       nrow(d$df_observed))

method_species_mean <- function(d) {
  ddf <- d$df_observed
  sp_means <- tapply(ddf$y, ddf$species, function(v) {
    v <- v[!is.na(v)]; if (length(v)) mean(v) else NA_real_
  })
  grand_mean <- mean(ddf$y, na.rm=TRUE)
  pred <- sp_means[as.character(ddf$species)]
  pred[is.na(pred)] <- grand_mean
  as.numeric(pred)
}

method_lm <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  fit_df <- ddf[!is.na(ddf$y), c("y", pred_n), drop = FALSE]
  fmla <- stats::as.formula(paste("y ~", paste(pred_n, collapse=" + ")))
  fit <- tryCatch(stats::lm(fmla, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))
  predict(fit, newdata = ddf)
}

method_lm_nonlinear <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  fit_df <- ddf[!is.na(ddf$y), c("y", pred_n), drop = FALSE]
  poly_terms <- sapply(pred_n, function(p) sprintf("poly(%s, 2, raw=TRUE)", p))
  inter_terms <- combn(pred_n, 2L, function(p) sprintf("%s:%s", p[1], p[2]))
  rhs <- paste(c(poly_terms, inter_terms), collapse=" + ")
  fmla <- stats::as.formula(paste("y ~", rhs))
  fit <- tryCatch(stats::lm(fmla, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(method_lm(d))
  tryCatch(predict(fit, newdata = ddf), error = function(e) rep(NA_real_, nrow(ddf)))
}

method_phylolm_blup <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  # Aggregate to species level
  sp_means_y <- tapply(ddf$y, ddf$species,
                        function(v) { v <- v[!is.na(v)]
                                       if (length(v)) mean(v) else NA_real_ })
  sp_means_X <- as.data.frame(lapply(pred_n, function(p)
    tapply(ddf[[p]], ddf$species, mean, na.rm=TRUE)))
  names(sp_means_X) <- pred_n
  sp_df <- data.frame(species = names(sp_means_y),
                       y = as.numeric(sp_means_y),
                       sp_means_X, stringsAsFactors = FALSE)
  sp_df <- sp_df[sp_df$species %in% d$tree$tip.label, , drop = FALSE]
  sp_df_obs <- sp_df[!is.na(sp_df$y), , drop = FALSE]
  if (nrow(sp_df_obs) < length(pred_n) + 5L) return(rep(NA_real_, nrow(ddf)))
  rownames(sp_df_obs) <- sp_df_obs$species
  tree_obs <- ape::keep.tip(d$tree, sp_df_obs$species)
  fmla <- stats::as.formula(paste("y ~", paste(pred_n, collapse=" + ")))
  fit <- tryCatch(
    phylolm::phylolm(fmla, data = sp_df_obs, phy = tree_obs, model = "lambda"),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(ddf)))
  beta_hat <- coef(fit)
  rhs_fmla <- stats::as.formula(paste("~", paste(pred_n, collapse=" + ")))
  Xmat <- model.matrix(rhs_fmla, sp_df)
  fixed <- as.numeric(Xmat %*% beta_hat)
  R <- ape::vcv(d$tree); R <- stats::cov2cor(R); R <- R[sp_df$species, sp_df$species]
  obs_idx <- which(!is.na(sp_df$y))
  miss_idx <- which(is.na(sp_df$y))
  if (length(miss_idx) == 0L) {
    sp_pred <- sp_df$y
  } else {
    e_obs <- sp_df$y[obs_idx] - fixed[obs_idx]
    R_oo <- R[obs_idx, obs_idx, drop = FALSE]
    R_mo <- R[miss_idx, obs_idx, drop = FALSE]
    blup <- as.numeric(R_mo %*% solve(R_oo + diag(1e-6, nrow(R_oo)), e_obs))
    sp_pred <- numeric(nrow(sp_df))
    sp_pred[obs_idx] <- sp_df$y[obs_idx]
    sp_pred[miss_idx] <- fixed[miss_idx] + blup
  }
  names(sp_pred) <- sp_df$species
  sp_pred[as.character(ddf$species)]
}

fit_pigauto_sf <- function(d, safety_floor, seed) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  cov_df <- as.data.frame(lapply(pred_n, function(p) as.numeric(ddf[[p]])))
  names(cov_df) <- pred_n
  res <- pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                          tree = d$tree, species_col = "species",
                          covariates = cov_df, missing_frac = 0.0,
                          epochs = EPOCHS, verbose = FALSE, seed = seed,
                          safety_floor = safety_floor)
  res$completed$y
}

method_pigauto_sfT <- function(d) fit_pigauto_sf(d, TRUE,  sample.int(1e7L,1L))
method_pigauto_sfF <- function(d) fit_pigauto_sf(d, FALSE, sample.int(1e7L,1L))

# ---------------- main loop -----------------------------------------------
cells <- expand.grid(
  n_covs   = N_COVS,
  lambda   = LAMBDAS,
  f_type   = F_TYPES,
  beta_str = BETAS,
  rep_id   = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_cells <- nrow(cells)
log_line("Total cells: ", total_cells)

method_funs <- list(
  column_mean         = method_column_mean,
  species_mean        = method_species_mean,
  lm                  = method_lm,
  lm_nonlinear        = method_lm_nonlinear,
  phylolm_lambda_blup = method_phylolm_blup,
  pigauto_sfT         = method_pigauto_sfT,
  pigauto_sfF         = method_pigauto_sfF
)

results <- list()
for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 17000L + i
  log_line(sprintf("Cell %d/%d: f=%s, lam=%.1f, beta=%.1f, n_covs=%d, rep=%d (seed=%d)",
                    i, total_cells, row$f_type, row$lambda, row$beta_str,
                    row$n_covs, row$rep_id, seed))
  d <- tryCatch(
    sim_multi_nonlinear_dgp(N_SPECIES, OBS_PER, row$n_covs,
                              row$lambda, row$f_type, row$beta_str, seed),
    error = function(e) { log_line("  DGP ERROR: ", e$message); NULL })
  if (is.null(d)) next
  truth <- d$df_complete$y
  for (m_name in names(method_funs)) {
    pred <- tryCatch(method_funs[[m_name]](d),
                      error = function(e) {
                        log_line("  ", m_name, " ERROR: ", e$message)
                        rep(NA_real_, nrow(d$df_observed))
                      })
    s <- score(pred, truth, d$mask)
    results[[length(results) + 1L]] <- data.frame(
      n_species = N_SPECIES, n_covs = row$n_covs,
      lambda = row$lambda, beta_strength = row$beta_str,
      f_type = row$f_type, rep = row$rep_id,
      method = m_name,
      rmse = s["rmse"], r = s["r"],
      stringsAsFactors = FALSE
    )
  }
  if (i %% 2L == 0L || i == total_cells) {
    df <- do.call(rbind, results)
    saveRDS(df, out_rds)
    log_line(sprintf("  [checkpoint] saved %d rows", nrow(df)))
  }
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results)
saveRDS(df, out_rds)
log_line(sprintf("=== DONE: %d cells, %d rows ===", total_cells, nrow(df)))
log_line(sprintf("  rds: %s", out_rds))
