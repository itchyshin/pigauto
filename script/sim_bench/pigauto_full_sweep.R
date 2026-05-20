#!/usr/bin/env Rscript
# =============================================================================
# script/sim_bench/pigauto_full_sweep.R
# =============================================================================
#
# Dan-comparable single-trait pigauto sweep (2026-05-19 rewrite).
#
# **What this is.** The high-rep version of overnight_2026_05_19_mechanisms.R.
# Same single-trait DGP (`.sim_complete()`), same mechanism injector
# (`.inject_missingness_single()`), same four scenarios, same three
# mechanisms, same evaluator schema (rmse / z_rmse / wall_sec / gate_mean).
# Only the design knobs change:
#
#   mech sim      : 4 scen x 3 mech x {200,500} x 10 sims =  240 reps, 4 methods
#   full sweep    : 4 scen x 3 mech x {500}     x 50 sims =  600 reps, 3 methods
#
# Dropping pigauto_bayes (mech sim showed bayes ~= default within 1 MC SE on
# every cell) halves pigauto wall and brings total wall to ~12h. Keeping n=500
# only matches Dan's BACE sim (BACE/dev/02_benchmark_simulated_full.R), so the
# two RDSs can be merged for the paper figure.
#
# **What this is NOT.** This is not a replication of Dan's multi-trait
# sim_bace() DGP. The earlier attempt at that (the BACE-side companion) tripped
# on multi_impute's handling of sim_bace's species-column layout and produced
# zero successful reps. Single-trait gives a cleaner contrast and is directly
# comparable to the mech sim, so cross-N + cross-mechanism + cross-mode
# conclusions all land in one schema.
#
# Output: dev/simulation_results_pigauto/results.rds (long format, same
# columns as the mech sim) + a per-rep results_partial.rds checkpoint.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ape)
})

# ---- Config ---------------------------------------------------------------
N_SIMS       <- 50L
N_SPECIES    <- c(500L)
SCENARIOS    <- c("bm_strong", "ou_strong", "nonlinear_cov", "weak_signal")
MECHANISMS   <- c("phylo_MAR", "trait_MAR", "trait_MNAR")
MISS_RATE    <- 0.30
DEP_STRENGTH <- 1.5
SEED_BASE    <- 20260519L

OUT_DIR <- file.path("dev", "simulation_results_pigauto")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== pigauto_full_sweep (Dan-comparable single-trait) ===\n")
cat("scenarios     :", paste(SCENARIOS, collapse = ", "), "\n")
cat("mechanisms    :", paste(MECHANISMS, collapse = ", "), "\n")
cat("n_species     :", paste(N_SPECIES, collapse = ", "), "\n")
cat("n_sims/cell   :", N_SIMS, "\n")
cat("total reps    :", length(SCENARIOS) * length(N_SPECIES) *
                        length(MECHANISMS) * N_SIMS, "\n")
cat("output        :", OUT_DIR, "\n\n")

# ---- DGP helpers (identical to overnight_2026_05_19_mechanisms.R) ---------
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
      linpred <- rep(0, n)
    }
  } else {
    linpred <- -DEP_STRENGTH * as.numeric(scale(y))
  }
  obj <- function(c) (mean(plogis(c + linpred)) - miss_rate)^2
  c_hat <- stats::optimize(obj, interval = c(-10, 10))$minimum
  p <- plogis(c_hat + linpred)
  miss <- stats::rbinom(n, 1L, p) == 1L
  miss_idx <- which(miss)
  list(y_obs = ifelse(miss, NA_real_, y),
       miss_idx = miss_idx,
       truth = y[miss_idx])
}

# ---- Method runners -------------------------------------------------------
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
.run_pigauto <- function(y_obs, miss_idx, tree, R, cov_df) {
  df <- data.frame(y = y_obs, row.names = tree$tip.label)
  t0 <- Sys.time()
  args <- list(df, tree, epochs = 300L, eval_every = 50L, patience = 5L,
               verbose = FALSE, seed = 1L, lambda_mode = "fixed_1")
  if (!is.null(cov_df)) args$covariates <- cov_df
  res <- tryCatch(do.call(impute, args), error = function(e) NULL)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(res)) return(list(pred = NA_real_, wall = wall, gate = NA_real_))
  pred <- res$completed[, "y"][miss_idx]
  gate <- if (!is.null(res$fit$calibrated_gates))
    mean(res$fit$calibrated_gates) else NA_real_
  list(pred = pred, wall = wall, gate = gate)
}

# ---- One replicate runs all 3 methods --------------------------------------
run_one_rep <- function(scenario, n, mechanism, sim_id) {
  sci <- as.integer(which(SCENARIOS == scenario))
  ni  <- as.integer(which(N_SPECIES == n))
  mi  <- as.integer(which(MECHANISMS == mechanism))
  # NOTE: shift seed offset away from the mech sim's seed space (10000*sci +
  # 1000*ni + 100*mi + sim_id, sim_id <= 10). full_sweep uses sim_id 1..50, so
  # adding +200000 makes the rep seeds disjoint from any mech sim cell.
  seed <- as.integer(SEED_BASE + 200000L + 10000L * sci + 1000L * ni +
                       100L * mi + sim_id)
  sim <- .sim_complete(scenario, n, seed = seed)
  m   <- .inject_missingness_single(sim$y_true, sim$tree, sim$R, sim$cov_df,
                                      mechanism, MISS_RATE, seed = seed + 1L)
  if (length(m$miss_idx) == 0L) return(NULL)

  methods <- list(
    column_mean     = .run_column_mean(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df),
    bm_kriging      = .run_bm_kriging(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df),
    pigauto_default = .run_pigauto(m$y_obs, m$miss_idx, sim$tree, sim$R, sim$cov_df)
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
cat("Saved", nrow(results), "rows of", nrow(design) * 3L, "expected.\n")

saveRDS(results, file.path(OUT_DIR, "results.rds"))
cat("Output written to:", file.path(OUT_DIR, "results.rds"), "\n\n")

agg <- aggregate(z_rmse ~ scenario + n + mechanism + method, data = results,
                  FUN = function(x) c(mean = mean(x), sd = stats::sd(x),
                                       n = length(x)))
print(agg)
