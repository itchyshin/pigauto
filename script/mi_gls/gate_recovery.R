#!/usr/bin/env Rscript
# script/mi_gls/gate_recovery.R
#
# G3 (docs/dev-log/mi-posterior/design.md section 5b): parameter recovery of
# the posterior sampler behind multi_impute(draws_method = "posterior").
#
# DGP: vec(Y) ~ N(0, Sigma_P %x% R + Sigma_E %x% I_n), R = cov2cor(vcv(tree)),
# tree = ape::rtree(n) (the same generator as script/mi_gls/dgp_v2.R),
# n = 1000, K = 2, unit total variance per trait
# (Sigma_P[k,k] = lambda_k, Sigma_E[k,k] = 1 - lambda_k, Sigma_E diagonal),
# 25% MCAR missingness per trait. Seeds 1, 2, 3 per setting, fixed in advance.
#
# Settings (lambda_1, lambda_2; phylogenetic correlation rho_P):
#   A (1, 1; 0.7)   B (0.5, 0.5; 0.7)   C (0.3, 0.9; 0)   D (0.05, 0.95; 0.7)
#
# Pass rules (all fits must pass):
#   1. |posterior-mean lambda_k - truth| < 0.1 for every trait;
#   2. |posterior-mean cov2cor(Sigma_P)[1,2] - rho_P| < 0.1 where both
#      lambda_k >= 0.3 (settings A, B, C);
#   3. |posterior-mean lambda_k - REML lambda_k| < 2 posterior SD, where the
#      REML lambda is #187's per-trait estimate (.mvn_resolve_lambda()).
#      #187's REML search is confined to [0.005, 0.995]
#      (.pagel_lambda_from_cache()), so the posterior mean is clamped to that
#      interval before the comparison; otherwise a posterior mean of 0.997
#      with SD 0.001 "disagrees" with a REML estimate that cannot exceed
#      0.995. The unclamped values are printed.
# Boundary bias (posterior mean - truth at lambda = 1 and lambda = 0.05) is
# printed even when it passes. Prints RECOVERY_OK only if every check passes.
#
# Usage (cwd = worktree root): Rscript script/mi_gls/gate_recovery.R
# Env: MIP_G3_CORES    fits run in parallel (default 4; parallel::mclapply)
#      MIP_G3_SETTINGS subset of settings, e.g. "B" or "A,D" (default all)
#      MIP_G3_SEEDS    subset of seeds, e.g. "1" (default "1,2,3")
# RECOVERY_OK is printed only for the full run (all settings, all seeds);
# a subset run prints its table and RECOVERY_SUBSET_OK / G3 FAILED.
# Cost: 4 chains x 6,000 sweeps per fit, chains serial within a fit.
# Measured on the Mac Studio (M1 Ultra), n = 1000, K = 2: see the report;
# the full 12-fit gate is sized for Totoro / DRAC with MIP_G3_CORES fits in
# parallel.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
suppressMessages(devtools::load_all(quiet = TRUE))

n <- 1000L
miss_frac <- 0.25
cores <- as.integer(Sys.getenv("MIP_G3_CORES", "4"))
settings <- data.frame(
  id = c("A", "B", "C", "D"),
  lambda1 = c(1, 0.5, 0.3, 0.05),
  lambda2 = c(1, 0.5, 0.9, 0.95),
  rho = c(0.7, 0.7, 0, 0.7))
seeds <- as.integer(strsplit(Sys.getenv("MIP_G3_SEEDS", "1,2,3"), ",")[[1]])
sel <- strsplit(Sys.getenv("MIP_G3_SETTINGS", "A,B,C,D"), ",")[[1]]
full_run <- setequal(sel, settings$id) && setequal(seeds, 1:3)
settings <- settings[settings$id %in% sel, , drop = FALSE]

simulate_g3 <- function(lambda, rho, seed) {
  set.seed(seed)
  tree <- ape::rtree(n)
  R <- stats::cov2cor(ape::vcv(tree))
  Cr <- matrix(c(1, rho, rho, 1), 2)
  SP <- diag(sqrt(lambda)) %*% Cr %*% diag(sqrt(lambda))
  Lr <- t(chol(R + diag(1e-10, n)))
  Y <- Lr %*% matrix(stats::rnorm(2 * n), n) %*% chol(SP)
  se <- sqrt(1 - lambda)
  Y <- Y + sweep(matrix(stats::rnorm(2 * n), n), 2L, se, "*")
  dimnames(Y) <- list(tree$tip.label, c("t1", "t2"))
  for (k in 1:2) Y[sample.int(n, round(miss_frac * n)), k] <- NA
  list(Y = Y, tree = tree)
}

ctl <- .mip_resolve_control(list(), m = 20L)
cat(sprintf(paste0("G3: n = %d, K = 2, %d chains x (%d burn-in + %d sweeps, ",
                   "thin %d); %d fits in parallel\n"),
            n, ctl$n_chains, ctl$burnin, ctl$n_iter, ctl$thin, cores))

jobs <- expand.grid(s = seq_len(nrow(settings)), seed = seeds)
one_fit <- function(j) {
  st <- settings[jobs$s[j], ]
  seed <- jobs$seed[j]
  lam <- c(st$lambda1, st$lambda2)
  d <- simulate_g3(lam, st$rho, seed)
  ctl_j <- ctl
  ctl_j$seed <- seed
  t0 <- proc.time()[["elapsed"]]
  f <- .mip_fit(d$Y, d$tree, ctl_j)
  wall <- proc.time()[["elapsed"]] - t0
  lam_mean <- colMeans(f$params$lambda)
  lam_sd <- apply(f$params$lambda, 2L, stats::sd)
  rho_hat <- mean(apply(f$params$Sigma_P, 3L,
                        function(S) stats::cov2cor(S)[1L, 2L]))
  reml <- f$start$lambda_reml
  chk_lam <- all(abs(lam_mean - lam) < 0.1)
  chk_rho <- if (all(lam >= 0.3)) abs(rho_hat - st$rho) < 0.1 else NA
  chk_reml <- all(abs(pmin(pmax(lam_mean, 0.005), 0.995) - reml) <
                    2 * lam_sd)
  dg <- f$diagnostics
  sweeps <- ctl$n_chains * (ctl$burnin + ctl$n_iter)
  data.frame(
    setting = st$id, seed = seed,
    lam1_true = lam[1], lam1_post = lam_mean[1], lam1_sd = lam_sd[1],
    lam1_reml = reml[1],
    lam2_true = lam[2], lam2_post = lam_mean[2], lam2_sd = lam_sd[2],
    lam2_reml = reml[2],
    rho_true = st$rho, rho_post = rho_hat,
    max_rhat = max(dg$rhat), min_ess = min(dg$ess_bulk),
    converged = isTRUE(attr(dg, "converged")),
    ok_lambda = chk_lam, ok_rho = chk_rho, ok_reml = chk_reml,
    wall_s = wall, ms_per_sweep = 1000 * f$wall_s / sweeps)
}
t_all <- proc.time()[["elapsed"]]
rows <- parallel::mclapply(seq_len(nrow(jobs)), one_fit,
                           mc.cores = min(cores, nrow(jobs)))
bad <- vapply(rows, inherits, logical(1), "try-error")
if (any(bad)) {
  cat("G3 ERROR:", as.character(rows[[which(bad)[1L]]]), "\n")
  quit(save = "no", status = 1L)
}
res <- do.call(rbind, rows)
cat(sprintf("\nTotal wall time: %.1f s\n\n", proc.time()[["elapsed"]] - t_all))

num <- vapply(res, is.numeric, logical(1))
out <- res
out[num] <- lapply(res[num], function(x) round(x, 3))
print(out, row.names = FALSE)

cat("\nBoundary bias (posterior mean - truth):\n")
mk <- function(w, x) if (length(x)) data.frame(where = w, bias = x) else NULL
bb <- rbind(
  mk("lambda = 1 (A, trait 1)", res$lam1_post[res$setting == "A"] - 1),
  mk("lambda = 1 (A, trait 2)", res$lam2_post[res$setting == "A"] - 1),
  mk("lambda = 0.05 (D, trait 1)", res$lam1_post[res$setting == "D"] - 0.05),
  mk("lambda = 0.95 (D, trait 2)", res$lam2_post[res$setting == "D"] - 0.95))
if (!is.null(bb) && nrow(bb)) {
  bb <- stats::aggregate(bias ~ where, data = bb,
                         FUN = function(x) c(mean = mean(x), min = min(x),
                                             max = max(x)))
  print(bb)
} else {
  cat("  (settings A and D not in this run)\n")
}
cat(sprintf("\nConvergence: %d of %d fits converged (max R-hat < 1.05, min ESS > 400)\n",
            sum(res$converged), nrow(res)))

fail <- !res$ok_lambda | (!is.na(res$ok_rho) & !res$ok_rho) | !res$ok_reml
if (any(fail)) {
  cat("\nG3 FAILED in:\n")
  print(out[fail, c("setting", "seed", "lam1_true", "lam1_post", "lam1_reml",
                    "lam2_true", "lam2_post", "lam2_reml", "rho_true",
                    "rho_post", "ok_lambda", "ok_rho", "ok_reml")],
        row.names = FALSE)
  quit(save = "no", status = 1L)
}
if (full_run) cat("RECOVERY_OK\n") else cat("RECOVERY_SUBSET_OK\n")
