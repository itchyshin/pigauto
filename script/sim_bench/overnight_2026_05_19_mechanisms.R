#!/usr/bin/env Rscript
# =============================================================================
# script/sim_bench/overnight_2026_05_19_mechanisms.R
# =============================================================================
#
# Sibling of overnight_2026_05_19.R, adding the MISSINGNESS-MECHANISM axis.
# Designed to merge cleanly with Dan's BACE simulation results when they
# arrive: same 3 mechanism names (phylo_MAR, trait_MAR, trait_MNAR),
# same parameter (DEP_STRENGTH = 1.5), single-trait adaptation of Dan's
# multi-trait injector.
#
# Design:
#   4 scenarios (bm_strong, ou_strong, nonlinear_cov, weak_signal)
#   x 2 N (200, 500)
#   x 3 mechanisms (phylo_MAR, trait_MAR, trait_MNAR)
#   x 10 sims
#   = 240 replicates
#
# Compared to overnight_2026_05_19.R (which used MCAR), this sim
# probes whether pigauto's results survive realistic non-random
# missingness patterns.
#
# Compute budget: ~10 hours wall on this laptop (serial, MPS).

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ape)
})

# ---- Config ---------------------------------------------------------------
N_SIMS       <- 10L
N_SPECIES    <- c(200L, 500L)
SCENARIOS    <- c("bm_strong", "ou_strong", "nonlinear_cov", "weak_signal")
MECHANISMS   <- c("phylo_MAR", "trait_MAR", "trait_MNAR")
MISS_RATE    <- 0.30
DEP_STRENGTH <- 1.5  # matches Dan's design (Penone 2014 ICCs)
SEED_BASE    <- 20260519L

OUT_DIR <- file.path("dev", "simulation_results_overnight_2026_05_19_mechanisms")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Overnight sim 2026-05-19 (mechanism axis) ===\n")
cat("scenarios     :", paste(SCENARIOS, collapse = ", "), "\n")
cat("mechanisms    :", paste(MECHANISMS, collapse = ", "), "\n")
cat("n_species     :", paste(N_SPECIES, collapse = ", "), "\n")
cat("n_sims/cell   :", N_SIMS, "\n")
cat("total reps    :", length(SCENARIOS) * length(N_SPECIES) *
                        length(MECHANISMS) * N_SIMS, "\n")
cat("output        :", OUT_DIR, "\n\n")

# ---- DGP helpers (same as overnight_2026_05_19.R) -------------------------
.sim_complete <- function(scenario, n, seed) {
  set.seed(seed)
  tree <- rcoal(n)
  R <- stats::cov2cor(vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  L_R <- t(chol(R))
  bm_signal <- as.numeric(L_R %*% rnorm(n))
  noise <- rnorm(n)
  cov_df <- NULL
  y <- switch(scenario,
    bm_strong = sqrt(0.9) * bm_signal + sqrt(0.1) * noise,
    ou_strong = {
      tree_ou <- pigauto:::transform_tree_pagel(tree, lambda = 0.3)
      R_ou <- stats::cov2cor(vcv.phylo(tree_ou))
      R_ou <- R_ou[tree$tip.label, tree$tip.label]
      L_ou <- t(chol(R_ou))
      ou_signal <- as.numeric(L_ou %*% rnorm(n))
      sqrt(0.9) * ou_signal + sqrt(0.1) * noise
    },
    nonlinear_cov = {
      env <- rnorm(n)
      cov_df <<- data.frame(env = env, row.names = tree$tip.label)
      sqrt(0.4) * bm_signal + sqrt(0.4) * sin(env * 2) + sqrt(0.2) * noise
    },
    weak_signal = sqrt(0.1) * bm_signal + sqrt(0.9) * noise
  )
  list(y_true = y, tree = tree, R = R, cov_df = cov_df)
}

# ---- Mechanism injector (single-trait adaptation of Dan's) ----------------
# For a single trait y on n tips, generate a missingness mask of expected
# rate `miss_rate` whose probability follows the chosen mechanism:
#   phylo_MAR  : P(miss) depends on a phylo-correlated latent z(tree)
#   trait_MAR  : P(miss) depends on the covariate (only available for
#                nonlinear_cov scenario; for others, falls back to MCAR
#                because there's no covariate to depend on)
#   trait_MNAR : P(miss) depends on y itself (the value being measured)
.inject_missingness_single <- function(y, tree, R, cov_df, mechanism,
                                         miss_rate, seed) {
  set.seed(seed)
  n <- length(y)
  if (identical(mechanism, "phylo_MAR")) {
    z <- as.numeric(MASS::mvrnorm(1L, mu = rep(0, n), Sigma = R))
    linpred <- -DEP_STRENGTH * as.numeric(scale(z))
  } else if (identical(mechanism, "trait_MAR")) {
    if (!is.null(cov_df) && "env" %in% colnames(cov_df)) {
      linpred <- -DEP_STRENGTH * as.numeric(scale(cov_df$env))
    } else {
      # No covariate available -> MCAR fallback
      linpred <- rep(0, n)
    }
  } else {  # trait_MNAR
    linpred <- -DEP_STRENGTH * as.numeric(scale(y))
  }
  # Calibrate intercept so expected miss rate matches target
  obj <- function(c) (mean(plogis(c + linpred)) - miss_rate)^2
  c_hat <- stats::optimize(obj, interval = c(-10, 10))$minimum
  p <- plogis(c_hat + linpred)
  miss <- stats::rbinom(n, 1L, p) == 1L
  miss_idx <- which(miss)
  list(y_obs = ifelse(miss, NA_real_, y),
       miss_idx = miss_idx,
       truth = y[miss_idx])
}

# ---- Method runners (same as overnight_2026_05_19.R) ----------------------
.run_column_mean <- function(y_obs, miss_idx, tree, R, cov_df) {
  t0 <- Sys.time()
  pred <- rep(mean(y_obs, na.rm = TRUE), length(miss_idx))
  list(pred = pred,
       wall = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       gate = NA_real_)
}
.run_bm_kriging <- function(y_obs, miss_idx, tree, R, cov_df) {
  t0 <- Sys.time()
  res <- pigauto:::bm_impute_col(y_obs, R)
  list(pred = res$mu[miss_idx],
       wall = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       gate = NA_real_)
}
.run_pigauto <- function(y_obs, miss_idx, tree, R, cov_df, lambda_mode) {
  df <- data.frame(y = y_obs, row.names = tree$tip.label)
  t0 <- Sys.time()
  args <- list(df, tree, epochs = 300L, eval_every = 50L, patience = 5L,
               verbose = FALSE, seed = 1L, lambda_mode = lambda_mode)
  if (!is.null(cov_df)) args$covariates <- cov_df
  res <- tryCatch(do.call(impute, args), error = function(e) NULL)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(res)) return(list(pred = NA_real_, wall = wall, gate = NA_real_))
  pred <- res$completed[, "y"][miss_idx]
  gate <- if (!is.null(res$fit$calibrated_gates))
    mean(res$fit$calibrated_gates) else NA_real_
  list(pred = pred, wall = wall, gate = gate)
}

# ---- One replicate runs all 4 methods --------------------------------------
run_one_rep <- function(scenario, n, mechanism, sim_id) {
  sci <- as.integer(which(SCENARIOS == scenario))
  ni  <- as.integer(which(N_SPECIES == n))
  mi  <- as.integer(which(MECHANISMS == mechanism))
  seed <- as.integer(SEED_BASE + 10000L * sci + 1000L * ni +
                       100L * mi + sim_id)
  sim <- .sim_complete(scenario, n, seed = seed)
  m   <- .inject_missingness_single(sim$y_true, sim$tree, sim$R, sim$cov_df,
                                      mechanism, MISS_RATE, seed = seed + 1L)
  if (length(m$miss_idx) == 0L) return(NULL)

  methods <- list(
    column_mean     = .run_column_mean(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df),
    bm_kriging      = .run_bm_kriging(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df),
    pigauto_default = .run_pigauto(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df, "fixed_1"),
    pigauto_bayes   = .run_pigauto(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df, "bayes")
  )
  rmse <- function(p, t) sqrt(mean((p - t)^2, na.rm = TRUE))
  z_sd <- stats::sd(m$y_obs, na.rm = TRUE)
  rows <- lapply(names(methods), function(nm) {
    out <- methods[[nm]]
    data.frame(
      scenario = scenario, n = n, mechanism = mechanism, sim_id = sim_id,
      method = nm,
      rmse = rmse(out$pred, m$truth),
      z_rmse = rmse(out$pred, m$truth) / z_sd,
      wall_sec = out$wall,
      gate_mean = out$gate,
      n_miss = length(m$miss_idx)
    )
  })
  do.call(rbind, rows)
}

# ---- Design + execution ----------------------------------------------------
design <- expand.grid(scenario = SCENARIOS, n = N_SPECIES,
                      mechanism = MECHANISMS,
                      sim_id = seq_len(N_SIMS),
                      stringsAsFactors = FALSE)
cat("Total replicates:", nrow(design), "\n")
cat("Started at:", format(Sys.time()), "\n\n")

CHECKPOINT_PATH <- file.path(OUT_DIR, "results_partial.rds")
results_list <- vector("list", nrow(design))
for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  cat(sprintf("[%s] rep %d/%d: %s mech=%s n=%d sim=%d\n",
              format(Sys.time(), "%H:%M:%S"),
              i, nrow(design), row$scenario, row$mechanism, row$n, row$sim_id))
  results_list[[i]] <- tryCatch(
    run_one_rep(row$scenario, row$n, row$mechanism, row$sim_id),
    error = function(e) {
      cat(sprintf("  REP FAILED: %s\n", e$message))
      NULL
    })
  partial <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  if (!is.null(partial) && nrow(partial) > 0L) {
    saveRDS(partial, CHECKPOINT_PATH)
  }
}

results <- do.call(rbind, results_list[!sapply(results_list, is.null)])
cat("\nFinished at:", format(Sys.time()), "\n")
cat("Saved", nrow(results), "rows of", nrow(design) * 4L, "expected.\n")

saveRDS(results, file.path(OUT_DIR, "results.rds"))
cat("Output written to:", file.path(OUT_DIR, "results.rds"), "\n\n")

# Per-(scenario, n, mechanism, method) summary for quick eyeballing
agg <- aggregate(z_rmse ~ scenario + n + mechanism + method, data = results,
                  FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x)))
print(agg)
