#!/usr/bin/env Rscript
# script/mi_gls/gate_calibration.R
#
# G3 (revised 2026-09-24, Shinichi at CP1): calibration of the posterior
# sampler behind multi_impute(draws_method = "posterior"). Replaces the
# per-fit 0.1 tolerance of gate_recovery.R, which is about 1.7 posterior SD
# at lambda = 0.5, so honest fits fail it (evidence/README.md).
#
# DGP: identical to gate_recovery.R (n = 1000, K = 2, 25% MCAR per trait,
# ape::rtree, unit total variance per trait, Sigma_E diagonal).
# Settings (lambda_1, lambda_2; rho_P):
#   A (1, 1; 0.7)  B (0.5, 0.5; 0.7)  C (0.3, 0.9; 0)  D (0.05, 0.95; 0.7)
#
# Two modes.
#   run   (cwd = the code tree to test; slow, run on Totoro):
#     Rscript <path>/gate_calibration.R run <out.rds> [settings=A,B,C,D]
#             [seed_from=101] [seed_to=150] [cores=4]
#     Export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 at
#     launch (setting them inside R is too late for BLAS). The code SHA is
#     read from ./SHA (a git archive) or git HEAD, and stored in every row.
#   check (cwd = worktree root; fast; the GATES.md G3 CHECK):
#     Rscript script/mi_gls/gate_calibration.R check
#             [rds=docs/dev-log/mi-posterior/evidence/g3_calibration.rds]
#
# Pass rules (all must hold; RECOVERY_OK printed only then):
#   1. every setting A-D present with >= 50 fits, all rows one code SHA;
#   2. converged fraction >= 0.98 per setting (max R-hat < 1.05, min ESS > 400);
#   3. 95% credible-interval coverage of the true lambda_k in [0.88, 1.00]
#      wherever the truth is inside (0, 1) (B, C, D; the lambda = 1 truths of
#      A sit on the boundary, where an interval of a (0, 1) quantity cannot
#      contain them, so A is judged by rule 5 only);
#   4. 95% coverage of the true rho_P in [0.88, 1.00] where both lambda_k >=
#      0.3 (A, B, C);
#   5. |mean over fits of (posterior mean - truth)| <= 0.03 for every lambda_k
#      in every setting, and for rho_P where rule 4 applies;
#   6. agreement with #187's REML lambda: in >= 90% of fits per setting,
#      |posterior mean (clamped to #187's search interval [0.005, 0.995]) -
#      REML| < 2 posterior SD for both traits.
# With 50 fits a coverage estimate has MCSE about 0.03, so [0.88, 1.00] is
# about 2.3 MCSE below 0.95.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "check"

settings_all <- data.frame(
  id = c("A", "B", "C", "D"),
  lambda1 = c(1, 0.5, 0.3, 0.05),
  lambda2 = c(1, 0.5, 0.9, 0.95),
  rho = c(0.7, 0.7, 0, 0.7))

if (identical(mode, "run")) {
  if (length(args) < 2L) stop("run mode needs <out.rds>", call. = FALSE)
  out <- args[[2L]]
  sel <- strsplit(if (length(args) >= 3L) args[[3L]] else "A,B,C,D", ",")[[1]]
  seed_from <- as.integer(if (length(args) >= 4L) args[[4L]] else 101L)
  seed_to <- as.integer(if (length(args) >= 5L) args[[5L]] else 150L)
  cores <- as.integer(if (length(args) >= 6L) args[[6L]] else 4L)
  suppressMessages(devtools::load_all(quiet = TRUE))
  code_sha <- if (file.exists("SHA")) readLines("SHA", n = 1L) else
    tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
  n <- 1000L
  miss_frac <- 0.25
  settings <- settings_all[settings_all$id %in% sel, , drop = FALSE]
  sim <- function(lambda, rho, seed) {
    set.seed(seed)
    tree <- ape::rtree(n)
    R <- stats::cov2cor(ape::vcv(tree))
    Cr <- matrix(c(1, rho, rho, 1), 2)
    SP <- diag(sqrt(lambda)) %*% Cr %*% diag(sqrt(lambda))
    Lr <- t(chol(R + diag(1e-10, n)))
    Y <- Lr %*% matrix(stats::rnorm(2 * n), n) %*% chol(SP)
    Y <- Y + sweep(matrix(stats::rnorm(2 * n), n), 2L, sqrt(1 - lambda), "*")
    dimnames(Y) <- list(tree$tip.label, c("t1", "t2"))
    for (k in 1:2) Y[sample.int(n, round(miss_frac * n)), k] <- NA
    list(Y = Y, tree = tree)
  }
  ctl <- .mip_resolve_control(list(), m = 20L)
  jobs <- expand.grid(s = seq_len(nrow(settings)), seed = seed_from:seed_to)
  q <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE)
  one <- function(j) {
    st <- settings[jobs$s[j], ]
    seed <- jobs$seed[j]
    lam <- c(st$lambda1, st$lambda2)
    d <- sim(lam, st$rho, seed)
    c1 <- ctl
    c1$seed <- seed
    f <- .mip_fit(d$Y, d$tree, c1)
    L <- f$params$lambda
    rho <- apply(f$params$Sigma_P, 3L, function(S) stats::cov2cor(S)[1L, 2L])
    data.frame(
      setting = st$id, seed = seed, code_sha = code_sha,
      lam1_true = lam[1], lam1_mean = mean(L[, 1]), lam1_sd = stats::sd(L[, 1]),
      lam1_lo = q(L[, 1])[1], lam1_hi = q(L[, 1])[2],
      lam1_reml = f$start$lambda_reml[1],
      lam2_true = lam[2], lam2_mean = mean(L[, 2]), lam2_sd = stats::sd(L[, 2]),
      lam2_lo = q(L[, 2])[1], lam2_hi = q(L[, 2])[2],
      lam2_reml = f$start$lambda_reml[2],
      rho_true = st$rho, rho_mean = mean(rho), rho_sd = stats::sd(rho),
      rho_lo = q(rho)[1], rho_hi = q(rho)[2],
      converged = isTRUE(attr(f$diagnostics, "converged")),
      max_rhat = max(f$diagnostics$rhat), min_ess = min(f$diagnostics$ess_bulk),
      wall_s = f$wall_s)
  }
  t0 <- proc.time()[["elapsed"]]
  rows <- parallel::mclapply(seq_len(nrow(jobs)), function(j)
    tryCatch(one(j), error = function(e) {
      message("fit ", j, " failed: ", conditionMessage(e)); NULL }),
    mc.cores = cores)
  res <- do.call(rbind, rows)
  saveRDS(res, out)
  cat(sprintf("CALIBRATION_RUN_DONE fits=%d of %d wall=%.0fs sha=%s out=%s\n",
              nrow(res), nrow(jobs), proc.time()[["elapsed"]] - t0, code_sha, out))
  quit(save = "no", status = 0L)
}

if (!identical(mode, "check")) stop("mode must be 'run' or 'check'", call. = FALSE)
rds <- if (length(args) >= 2L) args[[2L]] else
  file.path("docs", "dev-log", "mi-posterior", "evidence", "g3_calibration.rds")
if (!file.exists(rds)) {
  cat("G3 FAILED: calibration file not found:", rds, "\n")
  quit(save = "no", status = 1L)
}
res <- readRDS(rds)
fail <- character(0)
shas <- unique(res$code_sha)
if (length(shas) != 1L || is.na(shas)) {
  fail <- c(fail, sprintf("code SHA not unique or missing: %s", paste(shas, collapse = ",")))
}
cover <- function(t, lo, hi) mean(t >= lo & t <= hi)
in_band <- function(x) is.finite(x) && x >= 0.88 && x <= 1.00
rows <- list()
for (s in settings_all$id) {
  r <- res[res$setting == s, , drop = FALSE]
  st <- settings_all[settings_all$id == s, ]
  if (nrow(r) < 50L) {
    fail <- c(fail, sprintf("setting %s: %d fits (< 50)", s, nrow(r)))
    if (!nrow(r)) next
  }
  conv <- mean(r$converged)
  if (!(conv >= 0.98)) fail <- c(fail, sprintf("setting %s: converged fraction %.3f < 0.98", s, conv))
  lam_true <- c(st$lambda1, st$lambda2)
  cov_l <- c(cover(r$lam1_true, r$lam1_lo, r$lam1_hi), cover(r$lam2_true, r$lam2_lo, r$lam2_hi))
  bias_l <- c(mean(r$lam1_mean - r$lam1_true), mean(r$lam2_mean - r$lam2_true))
  for (k in 1:2) {
    if (lam_true[k] > 0 && lam_true[k] < 1 && !in_band(cov_l[k])) {
      fail <- c(fail, sprintf("setting %s: lambda_%d coverage %.3f outside [0.88, 1.00]", s, k, cov_l[k]))
    }
    if (!(abs(bias_l[k]) <= 0.03)) {
      fail <- c(fail, sprintf("setting %s: lambda_%d mean bias %+.4f (|.| > 0.03)", s, k, bias_l[k]))
    }
  }
  rho_rule <- all(lam_true >= 0.3)
  cov_r <- cover(r$rho_true, r$rho_lo, r$rho_hi)
  bias_r <- mean(r$rho_mean - r$rho_true)
  if (rho_rule) {
    if (!in_band(cov_r)) fail <- c(fail, sprintf("setting %s: rho_P coverage %.3f outside [0.88, 1.00]", s, cov_r))
    if (!(abs(bias_r) <= 0.03)) fail <- c(fail, sprintf("setting %s: rho_P mean bias %+.4f (|.| > 0.03)", s, bias_r))
  }
  cl <- function(x) pmin(pmax(x, 0.005), 0.995)
  agree <- abs(cl(r$lam1_mean) - r$lam1_reml) < 2 * r$lam1_sd &
    abs(cl(r$lam2_mean) - r$lam2_reml) < 2 * r$lam2_sd
  agree_frac <- mean(agree)
  if (!(agree_frac >= 0.90)) fail <- c(fail, sprintf("setting %s: REML agreement %.3f < 0.90", s, agree_frac))
  rows[[s]] <- data.frame(setting = s, fits = nrow(r), converged = conv,
    lam1_cover = cov_l[1], lam2_cover = cov_l[2], rho_cover = cov_r,
    lam1_bias = bias_l[1], lam2_bias = bias_l[2], rho_bias = bias_r,
    lam1_sd_of_means = stats::sd(r$lam1_mean), lam1_mean_post_sd = mean(r$lam1_sd),
    lam2_sd_of_means = stats::sd(r$lam2_mean), lam2_mean_post_sd = mean(r$lam2_sd),
    reml_agree = agree_frac, rho_gated = rho_rule)
}
tab <- do.call(rbind, rows)
cat(sprintf("G3 calibration check: %s (code SHA %s)\n\n", rds, paste(shas, collapse = ",")))
num <- vapply(tab, is.numeric, logical(1))
tab[num] <- lapply(tab[num], function(x) round(x, 3))
print(tab, row.names = FALSE)
cat("\n(lambda coverage is not gated where the truth is 1 (setting A); rho_P is gated where both lambda >= 0.3)\n")
if (length(fail)) {
  cat("\nG3 FAILED:\n", paste0("  ", fail, "\n"), sep = "")
  quit(save = "no", status = 1L)
}
cat("RECOVERY_OK\n")
