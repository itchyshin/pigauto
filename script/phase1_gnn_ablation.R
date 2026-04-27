#!/usr/bin/env Rscript
# script/phase1_gnn_ablation.R
#
# Phase 1 of the 2026-04-27 exploration plan:
# Force the GNN delta to zero by overriding calibrated gates after
# fitting, then re-predict.  Compare to full pigauto and to phylolm-
# lambda on the multi-obs nonlinear cells where pigauto currently wins.
#
# Decisive question: when r_cal_gnn = 0 and r_cal_bm = 1 (i.e., the
# GNN delta is fully gated out), does pigauto still beat phylolm-lambda?
#   - If yes: pigauto's win comes from the multi-obs aggregation /
#     GLS baseline / cov_linear / obs_refine, NOT from the GNN AE.
#   - If no: the GNN AE earns the win and the architecture story
#     (though incidental) is at least about the AE concept rather
#     than the multi-head transformer specifically.
#
# Also: inspect the calibrated gate values across all cells to see
# whether r_cal_gnn is meaningfully > 0 on cells where pigauto wins.

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH",
                          unset = "/Users/z3437171/Dropbox/Github Local/pigauto")
  devtools::load_all(pkg_path, quiet = TRUE)
})
options(warn = -1L)

HERE  <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(HERE, "script", "phase1_gnn_ablation.rds")
out_log <- file.path(HERE, "script", "phase1_gnn_ablation.log")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.0fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = "")

EPOCHS    <- 200L
MISS_FRAC <- 0.30
N_SPECIES <- 300L
OBS_PER   <- 5L
N_COVS    <- 5L
SP_MISSING_FRAC  <- 0.5
WITHIN_MISS_FRAC <- 0.2
SIGMA_RES <- 0.5

# Use the cells where pigauto currently wins decisively
LAMBDAS  <- c(0.1, 0.3)
F_TYPES  <- c("linear", "nonlinear", "interactive")
BETAS    <- c(0.5, 1.0)
N_REPS   <- 2L

log_line("=========================================================")
log_line("PHASE 1 GNN ABLATION")
log_line("Sweep: lambda=", paste(LAMBDAS, collapse=","),
          " | f=", paste(F_TYPES, collapse=","),
          " | beta=", paste(BETAS, collapse=","),
          " | reps=", N_REPS, " | epochs=", EPOCHS)
log_line("Methods: column_mean, species_mean, phylolm_lambda_blup,",
          " pigauto_full, pigauto_GNN_off, pigauto_baseline_only")
log_line("Output: ", out_rds)

# ----------------- DGP (same as bench_multi_obs_nonlinear) -----------------
sim_dgp <- function(n_species, obs_per, n_covs, lambda, f_type,
                     beta_strength, seed) {
  set.seed(seed)
  tree <- ape::rtree(n_species)
  sp <- tree$tip.label
  phylo_vals <- ape::rTraitCont(tree, model = "BM",
                                  sigma = sqrt(lambda), root.value = 0)
  phylo_vals <- as.numeric(phylo_vals[sp])
  names(phylo_vals) <- sp
  n_per_sp <- pmin(pmax(rpois(n_species, obs_per), 1L), 15L)
  names(n_per_sp) <- sp
  species_vec <- rep(sp, n_per_sp)
  n_total <- length(species_vec)
  Z <- matrix(rnorm(n_total * n_covs), nrow = n_total, ncol = n_covs)
  colnames(Z) <- paste0("z", seq_len(n_covs))
  contrib <- switch(f_type,
    linear      = as.numeric(Z %*% c(0.5, 0.3, 0.2, 0.1, -0.1)[seq_len(n_covs)]),
    nonlinear   = sin(Z[,1]) * exp(0.3 * Z[,2]) + 0.5 * tanh(Z[,3]) +
                    if (n_covs >= 5L) 0.3 * sin(2*Z[,4]) - 0.2 * cos(Z[,5]) else 0,
    interactive = Z[,1]*Z[,2] + 0.5*Z[,3]^2 +
                    if (n_covs >= 4L) 0.3*Z[,1]*Z[,4] else 0 +
                    if (n_covs >= 5L) -0.2*Z[,2]*Z[,5] else 0
  )
  contrib <- as.numeric(scale(contrib))
  delta <- beta_strength * contrib
  epsilon <- rnorm(n_total, 0, SIGMA_RES)
  y <- phylo_vals[species_vec] + delta + epsilon

  miss_sp <- sample(sp, floor(n_species * SP_MISSING_FRAC))
  y_obs <- y; y_obs[species_vec %in% miss_sp] <- NA
  obs_idx <- which(!is.na(y_obs))
  n_within <- floor(WITHIN_MISS_FRAC * length(obs_idx))
  if (n_within > 0L)
    y_obs[sample(obs_idx, n_within)] <- NA
  mask <- is.na(y_obs)
  df_obs <- data.frame(species = species_vec, y = y_obs, Z, stringsAsFactors = FALSE)
  df_full <- data.frame(species = species_vec, y = y, Z, stringsAsFactors = FALSE)
  list(tree = tree, df_complete = df_full, df_observed = df_obs,
        mask = mask, predictor_names = colnames(Z))
}

score <- function(pred, truth, mask) {
  ok <- mask & is.finite(truth) & is.finite(pred)
  if (sum(ok) < 5L) return(c(rmse = NA_real_, r = NA_real_))
  sqrt(mean((pred[ok] - truth[ok])^2)) -> rmse
  if (stats::sd(pred[ok]) > 1e-10) suppressWarnings(stats::cor(pred[ok], truth[ok])) -> r
  else r <- NA_real_
  c(rmse = rmse, r = r)
}

# ----------------- Methods -------------------------------------------------
m_column_mean <- function(d) rep(mean(d$df_observed$y, na.rm=TRUE),
                                  nrow(d$df_observed))

m_species_mean <- function(d) {
  ddf <- d$df_observed
  sp_means <- tapply(ddf$y, ddf$species,
                      function(v) { v <- v[!is.na(v)]; if (length(v)) mean(v) else NA_real_ })
  grand_mean <- mean(ddf$y, na.rm=TRUE)
  pred <- sp_means[as.character(ddf$species)]
  pred[is.na(pred)] <- grand_mean
  as.numeric(pred)
}

m_phylolm_blup <- function(d) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
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
  rhs <- stats::as.formula(paste("~", paste(pred_n, collapse=" + ")))
  Xmat <- model.matrix(rhs, sp_df)
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

# Pigauto with three configurations: full / GNN_off / baseline_only
fit_pigauto_with_ablation <- function(d, seed) {
  ddf <- d$df_observed
  pred_n <- d$predictor_names
  cov_df <- as.data.frame(lapply(pred_n, function(p) as.numeric(ddf[[p]])))
  names(cov_df) <- pred_n

  # Run impute() once -- gives us the fit + full prediction
  res <- pigauto::impute(traits = ddf[, c("species", "y"), drop = FALSE],
                          tree = d$tree, species_col = "species",
                          covariates = cov_df, missing_frac = 0.0,
                          epochs = EPOCHS, verbose = FALSE, seed = seed,
                          safety_floor = TRUE)

  # Full pigauto prediction: pigauto's res$completed$y
  pred_full <- res$completed$y

  # Inspect calibrated gate values
  fit <- res$fit
  cg <- fit$calibrated_gates
  rg <- fit$r_cal_gnn %||% cg
  rb <- fit$r_cal_bm  %||% (1 - cg)
  rm <- fit$r_cal_mean %||% rep(0, length(cg))

  # Save gate values for inspection
  gates <- list(r_cal_gnn = rg, r_cal_bm = rb, r_cal_mean = rm)

  # GNN-OFF prediction: override gates to put 100% on BM baseline.
  fit_zero <- fit
  fit_zero$calibrated_gates <- rep(0, length(cg))
  fit_zero$r_cal_gnn  <- rep(0, length(rg))
  fit_zero$r_cal_bm   <- rep(1, length(rb))
  fit_zero$r_cal_mean <- rep(0, length(rm))

  pred_obj_zero <- predict(fit_zero, return_se = FALSE)
  pred_zero_imp <- pred_obj_zero$imputed[, "y"]
  # The imputed prediction is in internal (tree-tip) order; remap to user-input
  # row order using input_row_order
  if (!is.null(res$data$input_row_order)) {
    imp_row <- match(seq_len(nrow(ddf)), res$data$input_row_order)
    pred_gnn_off <- numeric(nrow(ddf))
    obs_y <- ddf$y
    for (i in seq_len(nrow(ddf))) {
      if (!is.na(obs_y[i])) {
        pred_gnn_off[i] <- obs_y[i]
      } else {
        pred_gnn_off[i] <- pred_zero_imp[imp_row[i]]
      }
    }
  } else {
    pred_gnn_off <- pred_zero_imp
  }

  list(pred_full = pred_full, pred_gnn_off = pred_gnn_off, gates = gates)
}

# ----------------- Main loop -----------------------------------------------
cells <- expand.grid(
  lambda  = LAMBDAS,
  f_type  = F_TYPES,
  beta    = BETAS,
  rep_id  = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_cells <- nrow(cells)
log_line("Total cells: ", total_cells)

results <- list()
gate_log <- list()

for (i in seq_len(total_cells)) {
  row <- cells[i, ]
  seed <- 19000L + i
  log_line(sprintf("Cell %d/%d: f=%s, lam=%.1f, beta=%.1f, rep=%d (seed=%d)",
                    i, total_cells, row$f_type, row$lambda, row$beta, row$rep_id, seed))
  d <- tryCatch(sim_dgp(N_SPECIES, OBS_PER, N_COVS, row$lambda,
                          row$f_type, row$beta, seed),
                 error = function(e) { log_line("  DGP ERROR: ", e$message); NULL })
  if (is.null(d)) next

  truth <- d$df_complete$y

  # Baseline methods
  for (m_name in c("column_mean", "species_mean", "phylolm_lambda_blup")) {
    fn <- list(column_mean = m_column_mean,
                species_mean = m_species_mean,
                phylolm_lambda_blup = m_phylolm_blup)[[m_name]]
    pred <- tryCatch(fn(d), error = function(e) { log_line("  ", m_name, " ERROR: ", e$message); rep(NA_real_, nrow(d$df_observed)) })
    s <- score(pred, truth, d$mask)
    results[[length(results) + 1L]] <- data.frame(
      lambda = row$lambda, f_type = row$f_type, beta = row$beta,
      rep = row$rep_id, method = m_name,
      rmse = s["rmse"], r = s["r"], stringsAsFactors = FALSE)
  }

  # Pigauto with ablation
  ab <- tryCatch(fit_pigauto_with_ablation(d, seed),
                  error = function(e) { log_line("  pigauto ERROR: ", e$message); NULL })
  if (!is.null(ab)) {
    s_full <- score(ab$pred_full, truth, d$mask)
    s_off  <- score(ab$pred_gnn_off, truth, d$mask)
    results[[length(results) + 1L]] <- data.frame(
      lambda = row$lambda, f_type = row$f_type, beta = row$beta,
      rep = row$rep_id, method = "pigauto_full",
      rmse = s_full["rmse"], r = s_full["r"], stringsAsFactors = FALSE)
    results[[length(results) + 1L]] <- data.frame(
      lambda = row$lambda, f_type = row$f_type, beta = row$beta,
      rep = row$rep_id, method = "pigauto_GNN_off",
      rmse = s_off["rmse"], r = s_off["r"], stringsAsFactors = FALSE)
    gate_log[[length(gate_log) + 1L]] <- list(
      lambda = row$lambda, f_type = row$f_type, beta = row$beta, rep = row$rep_id,
      r_cal_gnn = ab$gates$r_cal_gnn,
      r_cal_bm  = ab$gates$r_cal_bm,
      r_cal_mean= ab$gates$r_cal_mean
    )
    log_line(sprintf("  full RMSE=%.4f  GNN-off RMSE=%.4f  diff=%.4f  r_cal_gnn=%s",
                      s_full["rmse"], s_off["rmse"],
                      s_full["rmse"] - s_off["rmse"],
                      paste(round(ab$gates$r_cal_gnn, 3), collapse = ",")))
  }

  if (i %% 2L == 0L || i == total_cells) {
    df <- do.call(rbind, results)
    saveRDS(list(results = df, gates = gate_log), out_rds)
  }
  invisible(gc(verbose = FALSE))
}

df <- do.call(rbind, results)
saveRDS(list(results = df, gates = gate_log), out_rds)
log_line(sprintf("=== DONE: %d cells, %d rows, %d gate logs ===",
                  total_cells, nrow(df), length(gate_log)))
