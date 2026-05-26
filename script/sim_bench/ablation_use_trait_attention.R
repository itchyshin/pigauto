#!/usr/bin/env Rscript
# =============================================================================
# script/sim_bench/ablation_use_trait_attention.R
# =============================================================================
# B3 LAZY-OPTIMIZER ABLATION for use_trait_attention (N=2000, Multi-Trait DGP)
#
# Tests whether the within-row cross-trait attention head added in v0.9.3
# is genuinely redundant against the joint-MVN baseline at BIEN scale, or
# whether a "lazy optimiser" is hiding useful signal behind the mu shortcut.
#
# Three arms per cell:
#   pigauto_OFF    use_trait_attention = FALSE  (the v0.9.3 default)
#   pigauto_ON     use_trait_attention = TRUE   (full B3 path)
#   pigauto_ON_L0  use_trait_attention = TRUE + baseline_mu masked +
#                  lambda_shrink = 0   (lazy-optimizer trap disarmed)
#
# Headline result (60 reps, 5 seeds per cell — see NEWS.md and
# useful/MEMO_2026-05-26_use_trait_attention_ablation.md):
#   z_rmse  pigauto_OFF 1.038 < pigauto_ON 1.056 < pigauto_ON_L0 1.057
# i.e. the attention path regresses against the joint-MVN baseline at scale.

# Load required libraries
library(ape)
library(methods)
library(MASS)
library(stats)
# Manually source all R files in the local 'R' folder to load functions
for (f in list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)) {
	source(f)
}

# ---- Config ----------------------
N_SIMS       <- 5L      
N_SPECIES    <- c(2000L) 
SCENARIOS    <- c("bm_strong", "ou_strong", "nonlinear_cov", "weak_signal")
MECHANISMS   <- c("phylo_MAR", "trait_MAR", "trait_MNAR")
MISS_RATE    <- 0.30
DEP_STRENGTH <- 1.5
SEED_BASE    <- 20260519L

OUT_DIR <- file.path("dev", "simulation_results_lazy_opt")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== pigauto_full_sweep (Lazy-Optimizer Ablation) ===\n")
cat("scenarios    :", paste(SCENARIOS, collapse = ", "), "\n")
cat("mechanisms   :", paste(MECHANISMS, collapse = ", "), "\n")
cat("n_species    :", paste(N_SPECIES, collapse = ", "), "\n")
cat("n_sims/cell  :", N_SIMS, "\n")
cat("total reps   :", length(SCENARIOS) * length(N_SPECIES) * length(MECHANISMS) * N_SIMS, "\n")
cat("output       :", OUT_DIR, "\n\n")

# ---- Multi-Trait DGP -------------------
.sim_complete <- function(scenario, n, seed) {
  set.seed(seed)
  tree <- rcoal(n)
  R <- stats::cov2cor(vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  L_R <- t(chol(R))
  bm_signal <- as.numeric(L_R %*% rnorm(n))
  noise <- rnorm(n)
  cov_df <- NULL
  
  y1 <- switch(scenario,
    bm_strong = sqrt(0.9) * bm_signal + sqrt(0.1) * noise,
    ou_strong = {
      tree_ou <- transform_tree_pagel(tree, lambda = 0.3)
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
  
  # Generate Secondary Trait (y2) highly correlated with y1 to engage Joint MVN
  y2 <- 0.7 * y1 + sqrt(1 - 0.7^2) * as.numeric(L_R %*% rnorm(n))
  y_df <- data.frame(y1 = y1, y2 = y2, row.names = tree$tip.label)
  
  list(y_true = y_df, tree = tree, R = R, cov_df = cov_df)
}

.inject_missingness_single <- function(y_df, tree, R, cov_df, mechanism, miss_rate, seed) {
  set.seed(seed)
  n <- nrow(y_df)
  y1 <- y_df$y1
  
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
    linpred <- -DEP_STRENGTH * as.numeric(scale(y1))
  }
  
  obj <- function(c) (mean(plogis(c + linpred)) - miss_rate)^2
  c_hat <- stats::optimize(obj, interval = c(-10, 10))$minimum
  p <- plogis(c_hat + linpred)
  miss <- stats::rbinom(n, 1L, p) == 1L
  miss_idx <- which(miss)
  
  y_obs <- y_df
  y_obs$y1[miss_idx] <- NA_real_ 
  
  list(y_obs = y_obs, miss_idx = miss_idx, truth = y1[miss_idx])
}

# ---- Method runners -------------------------------------------------------
.run_column_mean <- function(y_obs, miss_idx, truth, tree, R, cov_df) {
  t0 <- Sys.time()
  pred <- rep(mean(y_obs$y1, na.rm = TRUE), length(miss_idx))
  list(pred = pred, wall = as.numeric(difftime(Sys.time(), t0, units = "secs")), gate = NA_real_, coverage = NA_real_, width = NA_real_)
}

.run_bm_kriging <- function(y_obs, miss_idx, truth, tree, R, cov_df) {
  t0 <- Sys.time()
  res <- bm_impute_col(y_obs$y1, R)
  list(pred = res$mu[miss_idx], wall = as.numeric(difftime(Sys.time(), t0, units = "secs")), gate = NA_real_, coverage = NA_real_, width = NA_real_)
}

.run_pigauto <- function(y_obs, miss_idx, truth, tree, R, cov_df, ...) {
  t0 <- Sys.time()
  args <- list(y_obs, tree, epochs = 300L, eval_every = 50L, patience = 5L, verbose = FALSE, seed = 1L, lambda_mode = "fixed_1", ...)
  if (!is.null(cov_df)) args$covariates <- cov_df
  res <- tryCatch(do.call(impute, args), error = function(e) NULL)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  if (is.null(res)) return(list(pred = NA_real_, wall = wall, gate = NA_real_, coverage = NA_real_, width = NA_real_))
  
  pred <- res$completed[, "y1"][miss_idx]
  gate <- if (!is.null(res$fit$calibrated_gates)) mean(res$fit$calibrated_gates) else NA_real_
  
  coverage <- NA_real_
  width <- NA_real_
  if ("conformal_lower" %in% names(res$prediction) && "conformal_upper" %in% names(res$prediction)) {
    lower <- res$prediction$conformal_lower[miss_idx]
    upper <- res$prediction$conformal_upper[miss_idx]
    coverage <- mean(truth >= lower & truth <= upper, na.rm = TRUE)
    width <- mean(upper - lower, na.rm = TRUE)
  }
  
  list(pred = pred, wall = wall, gate = gate, coverage = coverage, width = width)
}

# ---- One replicate runs all methods --------------------------------------
run_one_rep <- function(scenario, n, mechanism, sim_id) {
  sci <- as.integer(which(SCENARIOS == scenario))
  ni  <- as.integer(which(N_SPECIES == n))
  mi  <- as.integer(which(MECHANISMS == mechanism))
  seed <- as.integer(SEED_BASE + 200000L + 10000L * sci + 1000L * ni + 100L * mi + sim_id)
  
  sim <- .sim_complete(scenario, n, seed = seed)
  m   <- .inject_missingness_single(sim$y_true, sim$tree, sim$R, sim$cov_df, mechanism, MISS_RATE, seed = seed + 1L)
  if (length(m$miss_idx) == 0L) return(NULL)

  methods <- list(
    column_mean   = .run_column_mean(m$y_obs, m$miss_idx, m$truth, sim$tree, sim$R, sim$cov_df),
    bm_kriging    = .run_bm_kriging(m$y_obs, m$miss_idx, m$truth, sim$tree, sim$R, sim$cov_df),
    pigauto_OFF   = .run_pigauto(m$y_obs, m$miss_idx, m$truth, sim$tree, sim$R, sim$cov_df, use_trait_attention = FALSE),
    pigauto_ON    = .run_pigauto(m$y_obs, m$miss_idx, m$truth, sim$tree, sim$R, sim$cov_df, use_trait_attention = TRUE),
    pigauto_ON_L0 = .run_pigauto(m$y_obs, m$miss_idx, m$truth, sim$tree, sim$R, sim$cov_df, use_trait_attention = TRUE, lambda_shrink = 0.0)
  )
  
  rmse <- function(p, t) sqrt(mean((p - t)^2, na.rm = TRUE))
  z_sd <- stats::sd(m$y_obs$y1, na.rm = TRUE)
  
  rows <- lapply(names(methods), function(nm) {
    out <- methods[[nm]]
    data.frame(
      scenario = scenario, n = n, mechanism = mechanism, sim_id = sim_id,
      method = nm, rmse = rmse(out$pred, m$truth), z_rmse = rmse(out$pred, m$truth) / z_sd,
      wall_sec = out$wall, gate_mean = out$gate, coverage = out$coverage, width = out$width, n_miss = length(m$miss_idx)
    )
  })
  do.call(rbind, rows)
}

# ---- Design + execution ----------------------------------------------------
design <- expand.grid(scenario = SCENARIOS, n = N_SPECIES,
                      mechanism = MECHANISMS, sim_id = seq_len(N_SIMS),
                      stringsAsFactors = FALSE)

CHECKPOINT_PATH <- file.path(OUT_DIR, "results_partial.rds")
results_list <- vector("list", nrow(design))

# --- RESUME BLOCK ---
if (file.exists(CHECKPOINT_PATH)) {
  done_res <- readRDS(CHECKPOINT_PATH)
  done_combos <- unique(paste(done_res$scenario, done_res$mechanism, done_res$sim_id))
  curr_combos <- paste(design$scenario, design$mechanism, design$sim_id)
  design <- design[!curr_combos %in% done_combos, ]
  cat(sprintf("Resuming from checkpoint: %d reps left to run.\n", nrow(design)))
}
# --------------------

for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  cat(sprintf("[%s] rep %d/%d: %s mech=%s n=%d sim=%d\n",
              format(Sys.time(), "%H:%M:%S"), i, nrow(design), 
              row$scenario, row$mechanism, row$n, row$sim_id))
  results_list[[i]] <- tryCatch(
    run_one_rep(row$scenario, row$n, row$mechanism, row$sim_id),
    error = function(e) { cat(sprintf("  REP FAILED: %s\n", e$message)); NULL }
  )
  partial <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  if (!is.null(partial) && nrow(partial) > 0L) saveRDS(partial, CHECKPOINT_PATH)
}

results <- do.call(rbind, results_list[!sapply(results_list, is.null)])
cat("\nFinished at:", format(Sys.time()), "\n")
saveRDS(results, file.path(OUT_DIR, "results.rds"))

agg <- aggregate(cbind(z_rmse, coverage, width) ~ scenario + mechanism + method, data = results, FUN = function(x) mean(x, na.rm=TRUE), na.action = na.pass)
print(agg)
