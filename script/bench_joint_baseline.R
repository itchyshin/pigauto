#!/usr/bin/env Rscript
# script/bench_joint_baseline.R
#
# A/B comparison: per-column BM vs joint MVN BM baseline on correlated traits.
#
# Purpose
#   Demonstrate that the Level-C Phase-2 joint multivariate BM baseline
#   (Rphylopars::phylopars) improves imputation accuracy over the legacy
#   per-column BM path when traits are correlated under a shared Brownian
#   motion process.
#
# Design
#   - Tree:     ape::rcoal(100) per rep
#   - Traits:   3 continuous traits; t2 and t3 correlated with t1 (0.7, 0.5)
#   - Holes:    25 random cells of t2 masked to NA (25% missingness)
#   - Metric:   RMSE on the held-out cells in latent (z-scored) space
#   - Reps:     10 independent trees + datasets
#   - Methods:  per_col_bm  — per-column BM via fit_baseline() with
#                              joint_mvn_available mocked to FALSE
#               joint_mvn_bm — fit_joint_mvn_baseline() directly
#
# Output
#   script/bench_joint_baseline.rds   list(results, agg, n_reps, commit)
#
# Run with
#   Rscript script/bench_joint_baseline.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-2",
    quiet = TRUE
  )
  library(ape)
})

WDIR    <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-2"
out_rds <- file.path(WDIR, "script", "bench_joint_baseline.rds")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%5.1fs] ", proc.time()[["elapsed"]] - t0), ..., "\n", sep = "")
  flush.console()
}

set.seed(2026)
n_reps  <- 10
results <- data.frame(rep = integer(), method = character(),
                      rmse = numeric(), stringsAsFactors = FALSE)

ns <- asNamespace("pigauto")

for (rep in seq_len(n_reps)) {

  log_line(sprintf("Rep %d/%d ...", rep, n_reps))

  # ---- Simulate correlated BM traits on a random tree ----------------------
  tree <- ape::rcoal(100)
  C    <- ape::vcv(tree)
  L    <- t(chol(C))                    # lower Cholesky: C = L %*% t(L)

  u1 <- stats::rnorm(100)
  u2 <- stats::rnorm(100)
  u3 <- stats::rnorm(100)

  # Traits share a common BM process; t2, t3 correlated with t1
  x1 <- as.numeric(L %*% u1)
  x2 <- 0.7 * x1 + as.numeric(L %*% u2) * sqrt(1 - 0.7^2)
  x3 <- 0.5 * x1 + as.numeric(L %*% u3) * sqrt(1 - 0.5^2)

  df <- data.frame(t1 = x1, t2 = x2, t3 = x3,
                   row.names = tree$tip.label)

  # ---- Mask 25 cells of t2 at random --------------------------------------
  miss <- sample(100, 25)
  truth_t2 <- df$t2          # full ground truth, same order as tip.label
  df$t2[miss] <- NA

  # ---- Preprocess + splits -------------------------------------------------
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                val_frac = 0.5, seed = rep,
                                trait_map = pd$trait_map)
  graph  <- build_phylo_graph(tree, k_eigen = "auto")

  # Latent-scale truth for the 25 masked cells of t2
  tm_t2 <- pd$trait_map[["t2"]]
  truth_latent <- (truth_t2 - tm_t2$mean) / tm_t2$sd

  # ---- Path A: joint MVN baseline (new) ------------------------------------
  bl_joint  <- fit_joint_mvn_baseline(pd, tree, splits = splits, graph = graph)
  # mu is (n_species x p_latent); rows are in tip.label order
  pred_joint <- bl_joint$mu[miss, "t2"]
  rmse_joint <- sqrt(mean((pred_joint - truth_latent[miss])^2))

  # ---- Path B: per-column BM (legacy) via namespace mock -------------------
  # fit_baseline() dispatches to joint when joint_mvn_available() == TRUE.
  # We replace the function in the pigauto namespace temporarily so the
  # existing dispatch logic falls through to the per-column loop.
  orig_fn  <- get("joint_mvn_available", envir = ns)
  mock_fn  <- function() FALSE

  tryCatch(unlockBinding("joint_mvn_available", ns),
           error = function(e) NULL)   # already unlocked under load_all
  assign("joint_mvn_available", mock_fn, envir = ns)

  bl_old   <- fit_baseline(pd, tree, splits = splits, graph = graph)

  assign("joint_mvn_available", orig_fn, envir = ns)
  tryCatch(lockBinding("joint_mvn_available", ns),
           error = function(e) NULL)

  pred_old <- bl_old$mu[miss, "t2"]
  rmse_old <- sqrt(mean((pred_old - truth_latent[miss])^2))

  results <- rbind(
    results,
    data.frame(rep = rep, method = "per_col_bm",   rmse = rmse_old),
    data.frame(rep = rep, method = "joint_mvn_bm", rmse = rmse_joint)
  )
}

# ---- Summary ---------------------------------------------------------------
agg <- aggregate(rmse ~ method, data = results, FUN = mean)
agg <- agg[order(agg$rmse), ]

cat("\nMean RMSE (lower = better):\n")
print(agg, row.names = FALSE)

rmse_old   <- agg$rmse[agg$method == "per_col_bm"]
rmse_joint <- agg$rmse[agg$method == "joint_mvn_bm"]
lift_pct   <- 100 * (rmse_old - rmse_joint) / rmse_old

cat(sprintf("\nLift (joint vs per-column): %.1f%%\n", lift_pct))
cat(sprintf("Wall time: %.1f s\n", proc.time()[["elapsed"]] - t0))

# ---- Save RDS --------------------------------------------------------------
commit <- tryCatch(
  system("git rev-parse HEAD", intern = TRUE),
  error = function(e) "unknown"
)

saveRDS(
  list(
    results  = results,
    agg      = agg,
    n_reps   = n_reps,
    lift_pct = lift_pct,
    commit   = commit
  ),
  out_rds
)
log_line("Wrote ", out_rds)
