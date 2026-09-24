#!/usr/bin/env Rscript
# script/mi_gls/diag/oracle_cell.R
#
# Diagnosis of the posterior-MI campaign (docs/dev-log/mi-posterior/diagnosis.md),
# failure 1 (negative paired slope bias in the lambda = 1 regimes).
#
# For ONE (regime, rep) of regimes 1-16 (same seed, tree, truth and masks as
# the campaign, via dgp_v2.R), imputes the masked cells m = 20 times with
# FIXED, known parameters (no sampler) and runs the campaign's downstream
# fits and pooling (verbatim from 01_cell_v2.R). Three arms:
#
#   oracle_true     exact draws of y_mis | y_obs under the DGP's OWN
#                   covariance, vec(Y) ~ N(0, Sig %x% V_sim), with
#                   V_sim = vcv(sim_tree) / max(vcv(sim_tree)) exactly as
#                   dgp_v2.R builds it (dense conditional; mu = 0 known).
#                   Tests H2 (a finite-n property of the paired estimand).
#   oracle_nominal  the package's fixed-parameter draws (.mip_fixed_draws(),
#                   the function behind G2) under the SAMPLER's covariance
#                   family Sigma_P %x% R + Sigma_E %x% I, R = cov2cor(vcv(tree)),
#                   at the nominal DGP values: lambda = 1: Sigma_P = Sig,
#                   Sigma_E = 1e-6 I; lambda < 1: Sigma_P = lambda Sig,
#                   Sigma_E = (1 - lambda) Sig; mu = 0 fixed.
#   oracle_kl       the same package draws at the in-model pseudo-true values
#                   Sigma_P = a* Sig, Sigma_E = b* Sig (+ 1e-6 I), where
#                   (a*, b*) minimise KL(N(0, V_sim) || N(0, a R + b I)).
#                   Because Sig %x% V_sim is separable and the model family
#                   is closed under linear maps of the traits, this pair is
#                   the KL projection of the DGP onto the sampler's model: the
#                   best parameters the sampler could possibly find.
#
# Sig = [[1, rho], [rho, 1]], rho = 0.7. The data are used on their raw
# scale (the sampler z-scores first; the z-scoring is linear, so draws at
# fixed parameters are equivalent).
#
# Usage (cwd = package source tree, e.g. frozen .../69670d44f9/code):
#   Rscript <path>/oracle_cell.R <regime_id> <rep> <outdir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("expected: regime_id rep outdir", call. = FALSE)
regime_id <- as.integer(args[[1L]])
rep_i     <- as.integer(args[[2L]])
outdir    <- args[[3L]]
if (regime_id > 16L) stop("oracle_cell.R covers regimes 1-16 only", call. = FALSE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({ devtools::load_all(quiet = TRUE); library(ape); library(nlme) })
source(file.path("script", "mi_gls", "dgp_v2.R"))
m <- 20L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(outdir, sprintf("oracle_regime_%d_rep_%d.rds", regime_id, rep_i))
if (file.exists(out_f)) { cat("SKIP", out_f, "\n"); quit(save = "no") }

cell <- simulate_regime_cell(regime_id, rep_i)
reg <- cell$regime
if (!is.null(reg$dgp) && !identical(reg$dgp, "tree_raw")) {
  stop("oracle_cell.R builds the raw-tree covariance of regimes 1-16 only; regime ",
       reg$regime_id, " has dgp '", reg$dgp, "'", call. = FALSE)
}
tree <- cell$tree
tips <- tree$tip.label
stopifnot(identical(rownames(cell$truth), tips))

# ---- the DGP's covariance (as in dgp_v2.R) ------------------------------------
sim_tree <- if (reg$lambda == 1) tree else transform_tree_pagel(tree, reg$lambda)
V_sim <- ape::vcv(sim_tree)[tips, tips]
V_sim <- V_sim / max(V_sim)
Sig <- matrix(reg$rho, 2, 2); diag(Sig) <- 1
n <- length(tips)

# ---- downstream fit + pooling (verbatim from 01_cell_v2.R) -------------------
fit_downstream <- function(dat) {
  dat$species <- rownames(dat)
  gls_fit <- tryCatch(
    nlme::gls(y ~ x, correlation = ape::corBrownian(1, cell$tree, form = ~species),
              data = dat, method = "ML"),
    error = function(e) NULL)
  phylolm_fit <- tryCatch(
    phylolm::phylolm(y ~ x, data = dat, phy = cell$tree, model = "lambda"),
    error = function(e) NULL)
  list(gls = gls_fit, phylolm = phylolm_fit)
}
single_est <- function(model_fit, kind) {
  if (is.null(model_fit)) return(c(estimate = NA_real_, se = NA_real_, df = NA_real_))
  if (kind == "phylolm") {
    co <- summary(model_fit)$coefficients
    c(estimate = unname(co["x", "Estimate"]), se = unname(co["x", "StdErr"]),
      df = model_fit$n - model_fit$d)
  } else {
    tt <- summary(model_fit)$tTable
    c(estimate = unname(tt["x", "Value"]), se = unname(tt["x", "Std.Error"]),
      df = model_fit$dims$N - model_fit$dims$p)
  }
}
pool_est <- function(fits_list, kind) {
  fits <- lapply(fits_list, `[[`, kind)
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) return(c(estimate = NA_real_, se = NA_real_, df = NA_real_))
  pooled <- tryCatch(suppressWarnings(pigauto::pool_mi(fits[ok])), error = function(e) e)
  if (inherits(pooled, "error")) return(c(estimate = NA_real_, se = NA_real_, df = NA_real_))
  row <- pooled[pooled$term == "x", , drop = FALSE]
  c(estimate = row$estimate, se = row$std.error, df = row$df)
}
arm_result <- function(arm, ymis_mat, miss_idx) {
  # ymis_mat: n_mis x m draws; miss_idx: n_mis x 2 (row, col) into truth
  datasets <- lapply(seq_len(ncol(ymis_mat)), function(d) {
    dat <- cell$df
    Z <- as.matrix(dat[, c("x", "y")])
    Z[miss_idx] <- ymis_mat[, d]
    dat$x <- Z[, 1L]; dat$y <- Z[, 2L]
    dat
  })
  fl <- lapply(datasets, fit_downstream)
  rbind(data.frame(arm = arm, downstream = "gls", t(pool_est(fl, "gls"))),
        data.frame(arm = arm, downstream = "phylolm", t(pool_est(fl, "phylolm"))))
}

# ---- complete data -------------------------------------------------------------
cf <- fit_downstream(cell$truth)
complete <- rbind(data.frame(arm = "complete", downstream = "gls", t(single_est(cf$gls, "gls"))),
                  data.frame(arm = "complete", downstream = "phylolm", t(single_est(cf$phylolm, "phylolm"))))

Y <- as.matrix(cell$df[, c("x", "y")])
miss_idx <- which(is.na(Y), arr.ind = TRUE)
miss_idx <- miss_idx[order(miss_idx[, 2L], miss_idx[, 1L]), , drop = FALSE]
draw_seed <- cell$seed + 555555L

# ---- arm 1: oracle_true (dense exact conditional under Sig %x% V_sim) --------
t0 <- proc.time()[["elapsed"]]
C <- kronecker(Sig, V_sim + 1e-8 * diag(n))
v <- as.vector(Y)
mi <- which(is.na(v)); o <- which(!is.na(v))
stopifnot(identical(mi, as.integer((miss_idx[, 2L] - 1L) * n + miss_idx[, 1L])))
Lo <- chol(C[o, o])
A <- backsolve(Lo, forwardsolve(t(Lo), C[o, mi]))       # C_oo^{-1} C_om
cmean <- as.vector(crossprod(A, v[o]))
ccov <- C[mi, mi] - crossprod(C[o, mi], A)
ccov <- (ccov + t(ccov)) / 2
Lc <- chol(ccov + 1e-10 * diag(length(mi)))
set.seed(draw_seed)
ymis_true <- cmean + crossprod(Lc, matrix(stats::rnorm(length(mi) * m), length(mi), m))
res_true <- arm_result("oracle_true", ymis_true, miss_idx)
wall_true <- proc.time()[["elapsed"]] - t0
rm(C, A, ccov, Lc, Lo)

# ---- in-model arms via the package's fixed-parameter draws --------------------
prob <- .mip_problem(Y[tree$tip.label, , drop = FALSE], tree)
stopifnot(identical(prob$miss, unname(miss_idx)))

# KL projection of N(0, V_sim) onto {a R + b I}, R = cov2cor(vcv(tree)).
R <- stats::cov2cor(ape::vcv(tree)[tips, tips])
ev <- eigen(R, symmetric = TRUE)
lam <- ev$values
w <- colSums(ev$vectors * (V_sim %*% ev$vectors))          # diag(U' V U)
kl <- function(p) {
  s <- exp(p[1L]) * lam + exp(p[2L])
  sum(w / s) + sum(log(s))
}
opt <- stats::optim(c(log(mean(diag(V_sim))), log(0.01)), kl, method = "L-BFGS-B",
                    lower = c(-20, -20), upper = c(5, 5))
a_star <- exp(opt$par[1L]); b_star <- exp(opt$par[2L])
c_star <- sum(w / lam) / n                                  # b fixed at 0

SP_nom <- if (reg$lambda == 1) Sig else reg$lambda * Sig
SE_nom <- if (reg$lambda == 1) diag(1e-6, 2L) else (1 - reg$lambda) * Sig
SP_kl <- a_star * Sig
SE_kl <- b_star * Sig + diag(1e-6, 2L)

t0 <- proc.time()[["elapsed"]]
set.seed(draw_seed)
fx_nom <- .mip_fixed_draws(prob, SP_nom, SE_nom, m, mu = c(0, 0))
res_nom <- arm_result("oracle_nominal", fx_nom$ymis, miss_idx)
wall_nom <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
set.seed(draw_seed)
fx_kl <- .mip_fixed_draws(prob, SP_kl, SE_kl, m, mu = c(0, 0))
res_kl <- arm_result("oracle_kl", fx_kl$ymis, miss_idx)
wall_kl <- proc.time()[["elapsed"]] - t0

res <- rbind(complete, res_true, res_nom, res_kl)
res$regime_id <- regime_id; res$rep <- rep_i
out <- list(regime_id = regime_id, rep = rep_i, seed = cell$seed, regime = reg,
            m = m, results = res,
            kl = c(a_star = a_star, b_star = b_star, lambda_star = a_star / (a_star + b_star),
                   c_star_b0 = c_star, mean_diag_V = mean(diag(V_sim)),
                   min_diag_V = min(diag(V_sim)), kl_conv = opt$convergence),
            n_missing = nrow(miss_idx),
            wall = c(true = wall_true, nominal = wall_nom, kl = wall_kl),
            pkg_sha = Sys.getenv("MI_POST_SHA", NA_character_),
            diag_sha = Sys.getenv("DIAG_SHA", NA_character_),
            host = unname(Sys.info()[["nodename"]]))
saveRDS(out, out_f)
cat(sprintf("OK oracle regime=%d rep=%d a*=%.3f b*=%.4f walls=%.0f/%.0f/%.0f\n",
            regime_id, rep_i, a_star, b_star, wall_true, wall_nom, wall_kl))
