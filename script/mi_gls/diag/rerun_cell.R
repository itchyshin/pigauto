#!/usr/bin/env Rscript
# script/mi_gls/diag/rerun_cell.R
#
# Diagnosis of the posterior-MI campaign (docs/dev-log/mi-posterior/diagnosis.md).
# Re-runs ONE (regime, rep) cell of script/mi_gls/01_cell_v2.R with the
# posterior sampler and saves what the campaign output did not keep:
# per-parameter convergence diagnostics and posterior summaries of Sigma_P,
# Sigma_E, lambda and the phylogenetic and residual correlations. The
# downstream fits and pooling are copied verbatim from 01_cell_v2.R, so with
# the default chain settings and the frozen campaign code the pooled slopes
# must reproduce the campaign output exactly (checked by the collector).
#
# Usage (cwd = the package source tree to load, e.g. the frozen campaign
# copy .../69670d44f9/code, or a scratch copy for the prior-sensitivity arm):
#   Rscript <path>/rerun_cell.R <regime_id> <rep> <outdir> <tag>
#
# Env:
#   MI_POST_NITER, MI_POST_BURNIN  chain length overrides (unset = package
#                                  defaults, as in the campaign)
#   DIAG_TWIN=1                    replace the regime 1-16 truth by its
#                                  in-model twin: each tip's row divided by
#                                  sqrt(diag(V_sim)[i]), so the data are
#                                  exactly N(0, Sig %x% cov2cor(V_sim)), the
#                                  sampler's own covariance family (same
#                                  seed, tree, noise and masks)
#   MI_POST_SHA                    SHA of the package tree loaded
#   DIAG_SHA                       SHA of this diagnosis script tree

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("expected: regime_id rep outdir tag", call. = FALSE)
regime_id <- as.integer(args[[1L]])
rep_i     <- as.integer(args[[2L]])
outdir    <- args[[3L]]
tag       <- args[[4L]]

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({ devtools::load_all(quiet = TRUE); library(ape); library(nlme) })
source(file.path("script", "mi_gls", "dgp_v2.R"))

# DIAG_PAR_CHAINS=k runs the k chains of ONE fit in parallel (forked). This
# only shortens the wall time of the long 2x / 4x re-runs. Each chain seeds
# itself (set.seed(chain_seed) at the top of .mip_run_chain()), so the chains,
# the proper draws, the diagnostics and the posterior_full slopes are
# identical to the serial run. The plug-in draws (posterior_none) are NOT:
# they are drawn after the chains from the global RNG stream, which the
# serial run leaves at chain 4's end state. The package source is not
# edited; the loaded namespace copy of .mip_fit() is swapped in memory.
par_chains <- as.integer(Sys.getenv("DIAG_PAR_CHAINS", "1"))
if (par_chains > 1L) {
  ns <- asNamespace("pigauto")
  src <- deparse(get(".mip_fit", envir = ns))
  hit <- grep("chains <- lapply(seq_len(ctl$n_chains), run1)", src, fixed = TRUE)
  if (length(hit) != 1L) stop("DIAG_PAR_CHAINS: .mip_fit() chain loop not found", call. = FALSE)
  src[hit] <- sub("lapply(seq_len(ctl$n_chains), run1)",
                  sprintf("parallel::mclapply(seq_len(ctl$n_chains), run1, mc.cores = %dL)", par_chains),
                  src[hit], fixed = TRUE)
  f <- eval(parse(text = src))
  environment(f) <- ns
  assignInNamespace(".mip_fit", f, ns = "pigauto")
}

env_int_or_null <- function(name) {
  v <- Sys.getenv(name, "")
  if (nzchar(v)) as.integer(v) else NULL
}
n_iter_override <- env_int_or_null("MI_POST_NITER")
burnin_override <- env_int_or_null("MI_POST_BURNIN")
twin <- identical(Sys.getenv("DIAG_TWIN", "0"), "1")
m <- 20L; keep_draws <- 1000L; n_chains <- 4L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(outdir, sprintf("%s_regime_%d_rep_%d.rds", tag, regime_id, rep_i))
if (file.exists(out_f)) { cat("SKIP", out_f, "\n"); quit(save = "no") }

cell <- simulate_regime_cell(regime_id, rep_i)
if (twin) {
  if (regime_id > 16L) stop("DIAG_TWIN applies to regimes 1-16 only", call. = FALSE)
  reg <- cell$regime
  sim_tree <- if (reg$lambda == 1) cell$tree else transform_tree_pagel(cell$tree, reg$lambda)
  V_sim <- ape::vcv(sim_tree); V_sim <- V_sim / max(V_sim)
  s <- sqrt(diag(V_sim))[rownames(cell$truth)]
  cell$truth$x <- cell$truth$x / s
  cell$truth$y <- cell$truth$y / s
  cell$df <- cell$truth
  cell$df$x[cell$mask_x] <- NA
  cell$df$y[cell$mask_y] <- NA
}

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
  if (is.null(model_fit)) return(c(estimate = NA_real_, se = NA_real_))
  if (kind == "phylolm") {
    co <- summary(model_fit)$coefficients
    c(estimate = unname(co["x", "Estimate"]), se = unname(co["x", "StdErr"]))
  } else {
    tt <- summary(model_fit)$tTable
    c(estimate = unname(tt["x", "Value"]), se = unname(tt["x", "Std.Error"]))
  }
}
pool_est <- function(fits_list, kind) {
  fits <- lapply(fits_list, `[[`, kind)
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) return(c(estimate = NA_real_, se = NA_real_, df = NA_real_))
  pooled <- tryCatch(pigauto::pool_mi(fits[ok]), error = function(e) e)
  if (inherits(pooled, "error")) return(c(estimate = NA_real_, se = NA_real_, df = NA_real_))
  row <- pooled[pooled$term == "x", , drop = FALSE]
  c(estimate = row$estimate, se = row$std.error, df = row$df)
}

cfit <- fit_downstream(cell$truth)
complete <- rbind(gls = single_est(cfit$gls, "gls"),
                  phylolm = single_est(cfit$phylolm, "phylolm"))

pu_mode <- if (identical(cell$regime$missing, "both")) "both" else "full"
ctl <- list(n_chains = n_chains, keep_draws = keep_draws,
            param_uncertainty = pu_mode, seed = cell$seed)
if (!is.null(n_iter_override)) ctl$n_iter <- n_iter_override
if (!is.null(burnin_override)) ctl$burnin <- burnin_override
t0 <- proc.time()[["elapsed"]]
mi <- suppressWarnings(multi_impute(cell$df, cell$tree, m = m,
                                    draws_method = "posterior",
                                    posterior_control = ctl))
wall <- proc.time()[["elapsed"]] - t0

fl <- lapply(mi$datasets, fit_downstream)
post <- rbind(gls = pool_est(fl, "gls"), phylolm = pool_est(fl, "phylolm"))
none <- NULL
if (!is.null(mi$posterior_improper)) {
  fl0 <- lapply(mi$posterior_improper$datasets, fit_downstream)
  none <- rbind(gls = pool_est(fl0, "gls"), phylolm = pool_est(fl0, "phylolm"))
}

# ---- posterior summaries (latent = z-scored scale) ----------------------------
pp <- mi$posterior$params
SP <- pp$Sigma_P; SE <- pp$Sigma_E
corr <- function(A) A[1L, 2L, ] / sqrt(A[1L, 1L, ] * A[2L, 2L, ])
psum <- data.frame(
  quantity = c("Sigma_P[1,1]", "Sigma_P[1,2]", "Sigma_P[2,2]",
               "Sigma_E[1,1]", "Sigma_E[1,2]", "Sigma_E[2,2]",
               "lambda_x", "lambda_y", "corr_P", "corr_E",
               "beta_P (SP12/SP22)", "beta_E (SE12/SE22)"),
  mean = c(mean(SP[1, 1, ]), mean(SP[1, 2, ]), mean(SP[2, 2, ]),
           mean(SE[1, 1, ]), mean(SE[1, 2, ]), mean(SE[2, 2, ]),
           mean(pp$lambda[, 1L]), mean(pp$lambda[, 2L]),
           mean(corr(SP)), mean(corr(SE)),
           mean(SP[1, 2, ] / SP[2, 2, ]), mean(SE[1, 2, ] / SE[2, 2, ])),
  median = c(stats::median(SP[1, 1, ]), stats::median(SP[1, 2, ]), stats::median(SP[2, 2, ]),
             stats::median(SE[1, 1, ]), stats::median(SE[1, 2, ]), stats::median(SE[2, 2, ]),
             stats::median(pp$lambda[, 1L]), stats::median(pp$lambda[, 2L]),
             stats::median(corr(SP)), stats::median(corr(SE)),
             stats::median(SP[1, 2, ] / SP[2, 2, ]), stats::median(SE[1, 2, ] / SE[2, 2, ])),
  stringsAsFactors = FALSE)

# ---- H3 checks: decode / back-transform -------------------------------------
tm <- mi$data$trait_map
obs_x <- !cell$mask_x; obs_y <- !cell$mask_y
max_obs_change <- max(vapply(mi$datasets, function(d) {
  max(abs(d[rownames(cell$df), "x"][obs_x] - cell$df$x[obs_x]),
      abs(d[rownames(cell$df), "y"][obs_y] - cell$df$y[obs_y]))
}, numeric(1)))
h3 <- list(
  log_transform = vapply(tm, function(t) isTRUE(t$log_transform), logical(1)),
  latent_mean = vapply(tm, function(t) t$mean, numeric(1)),
  latent_sd = vapply(tm, function(t) t$sd, numeric(1)),
  max_obs_change = max_obs_change,
  any_na_completed = any(vapply(mi$datasets, anyNA, logical(1))))

out <- list(
  regime_id = regime_id, rep = rep_i, seed = cell$seed, tag = tag, twin = twin,
  control = mi$posterior$control[c("n_chains", "burnin", "n_iter", "thin", "keep_draws",
                                   "param_uncertainty")],
  complete = complete, posterior_full = post, posterior_none = none,
  diagnostics = as.data.frame(mi$posterior$diagnostics),
  converged = isTRUE(mi$posterior$converged),
  posterior_summary = psum,
  plugin_Sigma_P = if (!is.null(mi$posterior_improper)) mi$posterior_improper$Sigma_P else NULL,
  plugin_Sigma_E = if (!is.null(mi$posterior_improper)) mi$posterior_improper$Sigma_E else NULL,
  hyper = mi$posterior$hyper,
  h3 = h3, wall_s = wall,
  pkg_sha = Sys.getenv("MI_POST_SHA", NA_character_),
  diag_sha = Sys.getenv("DIAG_SHA", NA_character_),
  host = unname(Sys.info()[["nodename"]]))
saveRDS(out, out_f)
cat(sprintf("OK %s regime=%d rep=%d wall=%.0fs converged=%s gls=%.5f phylolm=%.5f\n",
            tag, regime_id, rep_i, wall, out$converged,
            post["gls", "estimate"], post["phylolm", "estimate"]))
