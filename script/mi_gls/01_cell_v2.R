#!/usr/bin/env Rscript
# script/mi_gls/01_cell_v2.R
#
# One (regime, rep) job for the posterior-MI validation sweep
# (arc/mi-posterior). Extends script/mi_gls/01_cell.R to
# multi_impute(draws_method = "posterior") across ALL 24 regimes (see
# script/mi_gls/regimes.R). Usage:
#   OPENBLAS_NUM_THREADS=1 Rscript script/mi_gls/01_cell_v2.R <regime_id> <rep> <outdir>
#
# Methods:
#   complete        -- reference: fit on the (unobserved-in-practice) truth
#   posterior_full  -- multi_impute(draws_method = "posterior",
#                       posterior_control = list(param_uncertainty = "full"))
#   posterior_none  -- same, param_uncertainty = "none" (parameters fixed at
#                       their posterior mean; comparison only, per
#                       docs/dev-log/mi-posterior/design.md section 3)
#
# Downstream, per completed dataset:
#   gls     -- nlme::gls(y ~ x, corBrownian(1, tree, form = ~species), method = "ML")
#   phylolm -- phylolm::phylolm(y ~ x, phy = tree, model = "lambda")
# Pooled via pigauto::pool_mi(); "complete" (1 dataset) reports the raw
# model SE/df.
#
# Per-cell: coverage of the masked truth (and interval width) by
# mi$posterior$cell_interval's 95% predictive interval, split by missing
# trait (x, y). Diagnostics: max split R-hat and min bulk ESS over Sigma_P,
# Sigma_E, lambda (mi$posterior$diagnostics), plus the converged flag
# (attr(diagnostics, "converged")), saved for BOTH posterior_full and
# posterior_none so 04_acceptance.R / 05_cell_coverage.R can gate only on
# converged fits and flag regimes with excess non-convergence (design
# review item B3).
#
# Env vars: MI_POST_M (default 20), MI_POST_KEEP (default 1000),
# MI_POST_CHAINS (default 4).
#
# You cannot run this end to end until multi_impute() accepts
# draws_method = "posterior" (R/mi_posterior.R, built concurrently on this
# branch by another agent against the frozen API in
# docs/dev-log/mi-posterior/design.md section 4) -- this script checks for
# that up front and fails with a clear message if it is missing.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("expected: regime_id rep outdir", call. = FALSE)
regime_id <- as.integer(args[[1L]])
rep_i     <- as.integer(args[[2L]])
outdir    <- args[[3L]]
if (!is.finite(regime_id) || !is.finite(rep_i)) {
  stop("regime_id and rep must be integers", call. = FALSE)
}

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({ devtools::load_all(quiet = TRUE); library(ape); library(nlme) })
if (!requireNamespace("phylolm", quietly = TRUE)) {
  stop("package 'phylolm' is required by 01_cell_v2.R (the lambda downstream ",
      "model) but is not installed.", call. = FALSE)
}

source(file.path("script", "mi_gls", "dgp_v2.R"))

# ---- fail fast and clearly if the frozen API has not landed yet -------------
draws_methods <- tryCatch(eval(formals(multi_impute)$draws_method),
                          error = function(e) character(0))
if (!("posterior" %in% draws_methods)) {
  stop("multi_impute() does not yet support draws_method = \"posterior\" ",
      "(R/mi_posterior.R and the multi_impute() wiring have not landed on ",
      "this branch/worktree yet). This script codes against the frozen API ",
      "in docs/dev-log/mi-posterior/design.md section 4 and cannot run end ",
      "to end until that API exists. Currently supported draws_method ",
      "values: ", paste(draws_methods, collapse = ", "), ".", call. = FALSE)
}

m          <- as.integer(Sys.getenv("MI_POST_M", "20"))
keep_draws <- as.integer(Sys.getenv("MI_POST_KEEP", "1000"))
n_chains   <- as.integer(Sys.getenv("MI_POST_CHAINS", "4"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(outdir, sprintf("regime_%d_rep_%d.rds", regime_id, rep_i))
if (file.exists(out_f)) { cat("SKIP", out_f, "\n"); quit(save = "no") }

cell <- simulate_regime_cell(regime_id, rep_i)
true_beta <- cell$true_beta

# ---- downstream fit + pooling helpers ----------------------------------------
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

covers <- function(est, se, dfr, truth_val = true_beta, conf = 0.95) {
  if (!is.finite(est) || !is.finite(se) || !is.finite(dfr) || dfr <= 0) return(NA)
  tcrit <- stats::qt(1 - (1 - conf) / 2, dfr)
  truth_val >= est - tcrit * se & truth_val <= est + tcrit * se
}

na_row <- function(kind) {
  data.frame(downstream = kind, estimate = NA_real_, se = NA_real_,
            df = NA_real_, covered = NA)
}

single_row <- function(model_fit, kind) {
  if (is.null(model_fit)) return(na_row(kind))
  if (kind == "phylolm") {
    co  <- summary(model_fit)$coefficients
    est <- unname(co["x", "Estimate"]); se <- unname(co["x", "StdErr"])
    dfr <- model_fit$n - model_fit$d
  } else {
    tt  <- summary(model_fit)$tTable
    est <- unname(tt["x", "Value"]); se <- unname(tt["x", "Std.Error"])
    dfr <- model_fit$dims$N - model_fit$dims$p
  }
  data.frame(downstream = kind, estimate = est, se = se, df = dfr,
            covered = covers(est, se, dfr))
}

pool_row <- function(fits_list, kind) {
  fits <- lapply(fits_list, `[[`, kind)
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) return(na_row(kind))
  pooled <- tryCatch(pigauto::pool_mi(fits[ok]), error = function(e) e)
  if (inherits(pooled, "error")) return(na_row(kind))
  row <- pooled[pooled$term == "x", , drop = FALSE]
  data.frame(downstream = kind, estimate = row$estimate, se = row$std.error,
            df = row$df, covered = covers(row$estimate, row$std.error, row$df))
}

single_row_or_pool <- function(fits_list, kind) {
  if (length(fits_list) == 1L) single_row(fits_list[[1L]][[kind]], kind)
  else pool_row(fits_list, kind)
}

method_result <- function(method_name, datasets, wall_s) {
  fits_list <- lapply(datasets, fit_downstream)
  rbind(
    cbind(method = method_name, wall_s = wall_s,
         single_row_or_pool(fits_list, "gls")),
    cbind(method = method_name, wall_s = wall_s,
         single_row_or_pool(fits_list, "phylolm"))
  )
}

error_result <- function(method_name, wall_s) {
  cbind(method = method_name, wall_s = wall_s,
       rbind(na_row("gls"), na_row("phylolm")))
}

# ---- per-cell predictive-interval coverage helpers ---------------------------
cell_true_values <- function(ci, truth) {
  row_id <- ci$row
  idx <- if (is.character(row_id) || is.factor(row_id)) {
    match(as.character(row_id), rownames(truth))
  } else {
    as.integer(row_id)
  }
  vals <- rep(NA_real_, nrow(ci))
  for (tr in unique(ci$trait)) {
    sel <- ci$trait == tr
    vals[sel] <- truth[[tr]][idx[sel]]
  }
  vals
}

summarise_cell_coverage <- function(cell_df) {
  if (is.null(cell_df) || !nrow(cell_df)) {
    return(data.frame(trait = character(0), n = integer(0),
                      coverage = numeric(0), mean_width = numeric(0)))
  }
  do.call(rbind, lapply(split(cell_df, cell_df$trait), function(s) {
    data.frame(trait = s$trait[1L], n = nrow(s),
              coverage = mean(s$covered), mean_width = mean(s$width))
  }))
}

extract_diag <- function(mi) {
  d <- mi$posterior$diagnostics
  if (is.null(d) || !nrow(d)) {
    return(list(max_rhat = NA_real_, min_ess = NA_real_, converged = NA))
  }
  list(max_rhat  = max(d$rhat, na.rm = TRUE),
      min_ess   = min(d$ess_bulk, na.rm = TRUE),
      converged = isTRUE(attr(d, "converged")))
}

run_posterior <- function(param_uncertainty, method_name) {
  t0 <- proc.time()[["elapsed"]]
  mi <- tryCatch(
    multi_impute(cell$df, cell$tree, m = m, draws_method = "posterior",
                posterior_control = list(n_chains = n_chains, keep_draws = keep_draws,
                                         param_uncertainty = param_uncertainty,
                                         seed = cell$seed)),
    error = function(e) e)
  wall <- proc.time()[["elapsed"]] - t0
  if (inherits(mi, "error")) {
    return(list(method_res = error_result(method_name, wall),
               cell_detail = NULL, cell_summary = NULL,
               diag = list(max_rhat = NA_real_, min_ess = NA_real_, converged = NA),
               err = conditionMessage(mi)))
  }
  ci <- mi$posterior$cell_interval
  cell_detail <- if (!is.null(ci) && nrow(ci)) {
    true_vals <- cell_true_values(ci, cell$truth)
    data.frame(method = method_name, trait = ci$trait,
              covered = as.integer(true_vals >= ci$lower & true_vals <= ci$upper),
              width = ci$upper - ci$lower)
  } else NULL
  diag <- extract_diag(mi)
  list(method_res = method_result(method_name, mi$datasets, wall),
      cell_detail = cell_detail,
      cell_summary = if (!is.null(cell_detail)) {
        cbind(method = method_name, summarise_cell_coverage(cell_detail))
      } else NULL,
      diag = diag, err = NA_character_)
}

# ---- run methods --------------------------------------------------------------
results <- list()
cell_details   <- list()
cell_summaries <- list()
diagnostics    <- list()

t0 <- proc.time()[["elapsed"]]
complete_result <- method_result("complete", list(cell$truth), NA_real_)
complete_result$wall_s <- proc.time()[["elapsed"]] - t0
results$complete <- complete_result

posterior_errors <- list()
for (pu in c("full", "none")) {
  method_name <- paste0("posterior_", pu)
  out <- run_posterior(pu, method_name)
  results[[method_name]] <- out$method_res
  if (!is.null(out$cell_detail))   cell_details[[method_name]]   <- out$cell_detail
  if (!is.null(out$cell_summary))  cell_summaries[[method_name]] <- out$cell_summary
  diagnostics[[method_name]] <- data.frame(
    method = method_name, max_rhat = out$diag$max_rhat,
    min_ess = out$diag$min_ess, converged = out$diag$converged)
  posterior_errors[[method_name]] <- out$err
}

res_df <- do.call(rbind, results); rownames(res_df) <- NULL
cell_detail_df  <- if (length(cell_details))   do.call(rbind, cell_details)   else NULL
cell_summary_df <- if (length(cell_summaries)) do.call(rbind, cell_summaries) else NULL
diag_df <- do.call(rbind, diagnostics); rownames(diag_df) <- NULL

# ---- assemble + save -----------------------------------------------------------
out <- list(
  regime_id = regime_id, rep = rep_i, seed = cell$seed, regime = cell$regime,
  m = m, keep_draws = keep_draws, n_chains = n_chains, true_beta = true_beta,
  n_missing_x = sum(cell$mask_x), n_missing_y = sum(cell$mask_y),
  posterior_errors = posterior_errors,
  results = res_df,
  cell_detail = cell_detail_df,
  cell_coverage = cell_summary_df,
  diagnostics = diag_df
)
saveRDS(out, out_f)

cat(sprintf("OK regime=%d rep=%d n=%d mech=%s missing=%s\n",
           regime_id, rep_i, cell$regime$n, cell$regime$mechanism, cell$regime$missing))
for (nm in names(diagnostics)) {
  d <- diagnostics[[nm]]
  if (!isTRUE(d$converged)) {
    cat(sprintf("NONCONVERGED regime=%d rep=%d method=%s max_rhat=%.4f min_ess=%.1f\n",
               regime_id, rep_i, nm, d$max_rhat, d$min_ess))
  }
}
print(res_df)
