#!/usr/bin/env Rscript
# =============================================================================
# script/sim_bench/overnight_2026_05_19.R
# =============================================================================
#
# Focused overnight simulation: does pigauto's GNN earn its keep?
#
# Design choice: 4 scenarios chosen to STRESS-TEST pigauto's claim that
# its gated GNN provides value over a pure BM-kriging baseline AND over
# a column-mean floor. Two scenarios are "easy" (BM is the true model
# and a strong phylo baseline should win); two are "hard" for BM (the
# true model isn't BM, so a method that can flex beyond BM should win).
#
# Design:
#   - 4 scenarios  x  2 N_SPECIES (200, 500)  x  10 sims  =  80 replicates
#   - 4 imputation methods compared per replicate:
#         column_mean       grand-mean of observed (trivial floor)
#         bm_kriging        pure BM conditional MVN (the "right" model when
#                           DGP is BM; analogue of phylolm)
#         pigauto_default   pigauto with lambda_mode='fixed_1' (v0.10 default)
#         pigauto_bayes     pigauto with lambda_mode='bayes' (v0.11 opt-in)
#   - 30% MCAR mask on y
#   - Metric: per-cell RMSE on masked entries (z-scored)
#
# Scenarios:
#   bm_strong        y = sqrt(0.9)*BM + sqrt(0.1)*noise   (BM is right; easy)
#   ou_strong        y = sqrt(0.9)*OU(a=2) + sqrt(0.1)*noise  (BM is wrong)
#   nonlinear_cov    y = sqrt(0.4)*BM + sqrt(0.4)*sin(env) + sqrt(0.2)*noise
#                    plus 'env' as a covariate column. Linear BM baseline
#                    can't capture sin(env); pigauto's GNN with cov should.
#   weak_signal      y = sqrt(0.1)*BM + sqrt(0.9)*noise   (no phylo signal)
#                    Tests pigauto's safety floor: no method should win much.
#
# Compute budget at default config (m_imp=20, epochs=300):
#   n=200: ~1 min per pigauto fit
#   n=500: ~4 min per pigauto fit
#   80 reps total; pigauto x 2 modes per rep = 160 pigauto fits
#   sequential: ~5-7 hours; with 4 cores: ~1.5-2 hours
#
# Output: dev/simulation_results_overnight_2026_05_19/results.rds (one row
# per (scenario, n, method, sim_id) cell with RMSE + gate stats).
#
# Run:
#   nohup Rscript script/sim_bench/overnight_2026_05_19.R \
#     > /tmp/overnight_sim.log 2>&1 &

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ape)
  library(parallel)
})

# ---- Config ---------------------------------------------------------------
# SIM_ID_START lets the cascade launcher run "more sims" by shifting the
# seed range: stage 1 uses sim_id 1..10, stage 2 uses 11..20, stage 3
# uses 21..30. Each stage writes to a separate RDS so partial results
# from any stage are never overwritten.
N_SIMS       <- 10L
SIM_ID_START <- as.integer(Sys.getenv("PIGAUTO_SIM_ID_START", unset = "1"))
SIM_ID_TAG   <- Sys.getenv("PIGAUTO_SIM_TAG", unset = "stage1")
N_SPECIES    <- c(200L, 500L)
SCENARIOS    <- c("bm_strong", "ou_strong", "nonlinear_cov", "weak_signal")
# Serial only -- torch + MPS does not fork-safe; mclapply with mc.cores > 1
# segfaults on macOS when multiple workers try to share the MPS device.
# Serial MPS is faster than parallel CPU at this n, so the right tradeoff
# is to keep MPS and run one rep at a time.
N_CORES      <- 1L
MISS_RATE    <- 0.30
SEED_BASE    <- 20260519L

OUT_DIR      <- file.path("dev", "simulation_results_overnight_2026_05_19")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Overnight sim 2026-05-19 ===\n")
cat("scenarios     :", paste(SCENARIOS, collapse = ", "), "\n")
cat("n_species     :", paste(N_SPECIES, collapse = ", "), "\n")
cat("n_sims/cell   :", N_SIMS, "\n")
cat("total reps    :", length(SCENARIOS) * length(N_SPECIES) * N_SIMS, "\n")
cat("cores         :", N_CORES, "\n")
cat("output        :", OUT_DIR, "\n\n")

# ---- DGP helpers ---------------------------------------------------------
# Simulate a single replicate's complete data + tree + covariate (if any).
# Returns list(y_true, tree, R, cov_df). cov_df is NULL for scenarios w/o cov.
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
      # Approximate OU(alpha=2) via tree-transform: shrink the off-diagonal
      # of R toward zero by a Pagel-lambda-like factor. Lambda ~ 0.3 maps to
      # strong selection (most variance is non-phylo).
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

.apply_mask <- function(y, miss_rate, seed) {
  set.seed(seed)
  n <- length(y)
  miss_idx <- sample(n, floor(miss_rate * n))
  y_obs <- y
  y_obs[miss_idx] <- NA
  list(y_obs = y_obs, miss_idx = miss_idx, truth = y[miss_idx])
}

# ---- Methods: each returns (mu_pred[miss_idx], wall_sec, optional gate)
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

# ---- One replicate runs all 4 methods ------------------------------------
run_one_rep <- function(scenario, n, sim_id) {
  # Deterministic per-cell seed. Coerce everything to integer for set.seed.
  sci <- as.integer(which(SCENARIOS == scenario))
  ni  <- as.integer(which(N_SPECIES == n))
  seed <- as.integer(SEED_BASE + 1000L * sci + 100L * ni + sim_id)
  sim <- .sim_complete(scenario, n, seed = seed)
  m   <- .apply_mask(sim$y_true, MISS_RATE, seed = seed + 1L)

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
      scenario = scenario, n = n, sim_id = sim_id, method = nm,
      rmse = rmse(out$pred, m$truth),
      z_rmse = rmse(out$pred, m$truth) / z_sd,
      wall_sec = out$wall,
      gate_mean = out$gate
    )
  })
  do.call(rbind, rows)
}

# ---- Design + execution --------------------------------------------------
sim_id_seq <- SIM_ID_START + seq_len(N_SIMS) - 1L
design <- expand.grid(scenario = SCENARIOS, n = N_SPECIES,
                      sim_id = sim_id_seq,
                      stringsAsFactors = FALSE)
cat("  sim_id range  :", min(sim_id_seq), "..", max(sim_id_seq), "\n")
cat("  tag           :", SIM_ID_TAG, "\n")
cat("Total replicates:", nrow(design), "\n")
cat("Started at:", format(Sys.time()), "\n\n")

# Parallel over replicates. mclapply forks N_CORES R sessions.
# Each fork uses the already-loaded pigauto. Peak RAM is roughly
# N_CORES * (1.5 GB per pigauto fit at n=500) = 6 GB; well within
# laptop budget.
CHECKPOINT_PATH <- file.path(OUT_DIR, paste0("results_partial_", SIM_ID_TAG, ".rds"))
results_list <- vector("list", nrow(design))
for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  cat(sprintf("[%s] rep %d/%d: %s n=%d sim=%d\n",
              format(Sys.time(), "%H:%M:%S"),
              i, nrow(design), row$scenario, row$n, row$sim_id))
  results_list[[i]] <- tryCatch(
    run_one_rep(row$scenario, row$n, row$sim_id),
    error = function(e) {
      cat(sprintf("  REP FAILED: %s\n", e$message))
      NULL
    })
  # Checkpoint after each rep so an overnight crash doesn't lose work.
  partial <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  if (!is.null(partial) && nrow(partial) > 0L) {
    saveRDS(partial, CHECKPOINT_PATH)
  }
}

results <- do.call(rbind, results_list[!sapply(results_list, is.null)])
cat("\nFinished at:", format(Sys.time()), "\n")
cat("Saved", nrow(results), "rows of", nrow(design) * 4L, "expected.\n")

final_path <- file.path(OUT_DIR, paste0("results_", SIM_ID_TAG, ".rds"))
saveRDS(results, final_path)
cat("Output written to:", final_path, "\n\n")

# ---- Summary table for quick eyeballing ---------------------------------
agg <- aggregate(z_rmse ~ scenario + n + method, data = results,
                  FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x)))
print(agg)
