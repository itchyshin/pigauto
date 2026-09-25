#!/usr/bin/env Rscript
# script/mi_gls/01_cell_v2.R
#
# One (regime, rep) job for the posterior-MI validation sweep
# (arc/mi-posterior). Extends script/mi_gls/01_cell.R to
# multi_impute(draws_method = "posterior") across ALL 40 regimes (see
# script/mi_gls/regimes.R; 25-40 are the in-model twins of 1-16, simulated
# in script/mi_gls/dgp_v2.R with the source regime's seed). Usage:
#   OPENBLAS_NUM_THREADS=1 Rscript script/mi_gls/01_cell_v2.R <regime_id> <rep> <outdir>
#
# Methods:
#   complete        -- reference: fit on the (unobserved-in-practice) truth
#   posterior_full  -- multi_impute(draws_method = "posterior",
#                       posterior_control = list(param_uncertainty = "full"))
#   posterior_none  -- improper plug-in draws (covariances fixed at their
#                       posterior mean; comparison only, per
#                       docs/dev-log/mi-posterior/design.md section 3). Runs
#                       ONLY in regimes with missing == "both" (it exists for
#                       G6 rule 4; design.md section 5c, D6). There the
#                       sampler runs ONCE with param_uncertainty = "both":
#                       posterior_full is mi$datasets / mi$posterior (exactly
#                       as "full" with the same seed) and posterior_none is
#                       mi$posterior_improper, from the same chains, so the
#                       two share diagnostics (review estimands#4: no second
#                       4-chain run). Other regimes use "full".
#   conformal       -- per-cell comparator only (no downstream rows):
#                       impute(<same masked df>, tree, gnn = FALSE,
#                       seed = <cell seed>), pigauto defaults otherwise;
#                       prediction$conformal_lower/upper scored on exactly
#                       the masked cells (design.md section 5c, D7).
#
# Downstream, per completed dataset:
#   gls     -- nlme::gls(y ~ x, corBrownian(1, tree, form = ~species), method = "ML")
#   phylolm -- phylolm::phylolm(y ~ x, phy = tree, model = "lambda")
# Pooled via pigauto::pool_mi(); "complete" (1 dataset) reports the raw
# model SE/df. Each row keeps estimate, se, df and n_ok (downstream fits
# pooled). `covered_pop` is coverage of regimes$true_beta_pop and is
# descriptive in regimes 17-24 (exact, 0.7, in 1-16 and 25-40);
# 03_summarise_v2.R recomputes the gated coverage from estimate/se/df
# against the right truth (D1).
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
# MI_POST_CHAINS (default 4), MI_POST_NITER and MI_POST_BURNIN (default
# unset, meaning the package defaults; set them only for cheap smoke runs,
# never for the campaign), MI_POST_SHA (git SHA of the frozen code tree,
# recorded in the output; 12_totoro_campaign.sh and 11_fir_array.sbatch
# set it).
#
# Provenance (D8): the output records code_sha, pigauto_version, R version
# and host. The package is loaded with devtools::load_all() from the
# working directory; the campaign runs from a git archive of one SHA, so
# that directory is frozen.
#
# The script checks up front that multi_impute() accepts
# draws_method = "posterior" and fails with a clear message if it does not.

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
pu_modes <- tryCatch(.mip_default_control()$param_uncertainty, error = function(e) character(0))
if (!("both" %in% pu_modes)) {
  stop("posterior_control$param_uncertainty = \"both\" is not available in this ",
      "code tree (needed for posterior_none in the both-missing regimes).", call. = FALSE)
}

m          <- as.integer(Sys.getenv("MI_POST_M", "20"))
keep_draws <- as.integer(Sys.getenv("MI_POST_KEEP", "1000"))
n_chains   <- as.integer(Sys.getenv("MI_POST_CHAINS", "4"))
env_int_or_null <- function(name) {
  v <- Sys.getenv(name, "")
  if (nzchar(v)) as.integer(v) else NULL
}
n_iter_override <- env_int_or_null("MI_POST_NITER")    # NULL = package default
burnin_override <- env_int_or_null("MI_POST_BURNIN")   # NULL = package default

# ---- provenance (D8) -----------------------------------------------------------
code_sha <- Sys.getenv("MI_POST_SHA", "")
code_sha_source <- "MI_POST_SHA"
if (!nzchar(code_sha)) {
  code_sha <- tryCatch(suppressWarnings(system2("git", c("rev-parse", "HEAD"),
                                                stdout = TRUE, stderr = FALSE)),
                       error = function(e) character(0))
  if (length(code_sha) == 1L && grepl("^[0-9a-f]{7,40}$", code_sha)) {
    code_sha_source <- "git HEAD of the working tree (uncommitted edits not captured)"
  } else {
    code_sha <- NA_character_
    code_sha_source <- "unknown"
  }
}
pigauto_version <- tryCatch(as.character(utils::packageVersion("pigauto")),
                            error = function(e) {
                              tryCatch(unname(read.dcf("DESCRIPTION", fields = "Version")[1, 1]),
                                       error = function(e2) NA_character_)
                            })

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

# Coverage of regimes$true_beta_pop (column covered_pop). Exact truth in
# regimes 1-16 and 25-40, descriptive in 17-24 (see header and D1).
covers <- function(est, se, dfr, truth_val = true_beta, conf = 0.95) {
  if (!is.finite(est) || !is.finite(se) || !is.finite(dfr) || dfr <= 0) return(NA)
  tcrit <- stats::qt(1 - (1 - conf) / 2, dfr)
  truth_val >= est - tcrit * se & truth_val <= est + tcrit * se
}

na_row <- function(kind) {
  data.frame(downstream = kind, estimate = NA_real_, se = NA_real_,
            df = NA_real_, n_ok = NA_integer_, covered_pop = NA)
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
            n_ok = 1L, covered_pop = covers(est, se, dfr))
}

pool_row <- function(fits_list, kind) {
  fits <- lapply(fits_list, `[[`, kind)
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) return(na_row(kind))
  pooled <- tryCatch(pigauto::pool_mi(fits[ok]), error = function(e) e)
  if (inherits(pooled, "error")) return(na_row(kind))
  row <- pooled[pooled$term == "x", , drop = FALSE]
  data.frame(downstream = kind, estimate = row$estimate, se = row$std.error,
            df = row$df, n_ok = sum(ok),
            covered_pop = covers(row$estimate, row$std.error, row$df))
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
    return(list(max_rhat = NA_real_, min_ess = NA_real_, converged = NA,
                n_extensions = NA_integer_, sweeps_per_chain = NA_integer_))
  }
  # n_extensions / sweeps_per_chain: automatic chain extension (design.md 5e);
  # NULL on fits from code before that feature, recorded as NA.
  ne <- attr(d, "n_extensions"); sp <- attr(d, "sweeps_per_chain")
  list(max_rhat  = max(d$rhat, na.rm = TRUE),
      min_ess   = min(d$ess_bulk, na.rm = TRUE),
      converged = isTRUE(attr(d, "converged")),
      n_extensions = if (is.null(ne)) NA_integer_ else as.integer(ne),
      sweeps_per_chain = if (is.null(sp)) NA_integer_ else as.integer(sp))
}

na_diag <- list(max_rhat = NA_real_, min_ess = NA_real_, converged = NA,
                n_extensions = NA_integer_, sweeps_per_chain = NA_integer_)

# Score one posterior result set (datasets + cell_interval) as one method.
score_posterior <- function(method_name, datasets, ci, diag, wall) {
  cell_detail <- if (!is.null(ci) && nrow(ci)) {
    true_vals <- cell_true_values(ci, cell$truth)
    data.frame(method = method_name, trait = ci$trait,
              covered = as.integer(true_vals >= ci$lower & true_vals <= ci$upper),
              width = ci$upper - ci$lower)
  } else NULL
  list(method_res = method_result(method_name, datasets, wall),
      cell_detail = cell_detail,
      cell_summary = if (!is.null(cell_detail)) {
        cbind(method = method_name, summarise_cell_coverage(cell_detail))
      } else NULL,
      diag = diag, err = NA_character_)
}

# One sampler run. param_uncertainty = "full" returns posterior_full;
# "both" returns posterior_full and posterior_none from the same chains.
# Returns a named list of per-method results plus the resolved control.
run_posterior <- function(param_uncertainty) {
  methods <- if (identical(param_uncertainty, "both")) {
    c("posterior_full", "posterior_none")
  } else "posterior_full"
  ctl <- list(n_chains = n_chains, keep_draws = keep_draws,
              param_uncertainty = param_uncertainty, seed = cell$seed)
  if (!is.null(n_iter_override)) ctl$n_iter <- n_iter_override
  if (!is.null(burnin_override)) ctl$burnin <- burnin_override
  t0 <- proc.time()[["elapsed"]]
  mi <- tryCatch(
    multi_impute(cell$df, cell$tree, m = m, draws_method = "posterior",
                posterior_control = ctl),
    error = function(e) e)
  wall <- proc.time()[["elapsed"]] - t0   # one run, shared by both methods
  failed <- function(nm, msg) {
    list(method_res = error_result(nm, wall), cell_detail = NULL,
         cell_summary = NULL, diag = na_diag, err = msg)
  }
  if (inherits(mi, "error")) {
    res <- lapply(methods, failed, msg = conditionMessage(mi))
    names(res) <- methods
    return(list(res = res, ctl = NULL))
  }
  diag <- extract_diag(mi)
  res <- list(posterior_full = score_posterior("posterior_full", mi$datasets,
                                               mi$posterior$cell_interval, diag, wall))
  if ("posterior_none" %in% methods) {
    imp <- mi$posterior_improper
    res$posterior_none <- if (is.null(imp) || !length(imp$datasets)) {
      failed("posterior_none", "multi_impute() returned no posterior_improper result")
    } else {
      # Same chains as posterior_full, so the same convergence diagnostics.
      score_posterior("posterior_none", imp$datasets, imp$cell_interval, diag, wall)
    }
  }
  pc <- mi$posterior$control
  list(res = res,
       ctl = if (!is.null(pc)) pc[intersect(c("n_chains", "burnin", "n_iter", "thin", "keep_draws"),
                                            names(pc))] else NULL)
}

# Conformal comparator (D7): impute() on the SAME masked df, gnn = FALSE,
# pigauto defaults otherwise (verbose off only). prediction$conformal_*
# rows are species (tip) names; they are matched to the truth by name, on
# exactly the masked cells (mask_x for x, mask_y for y), on the original
# scale, the same cells and truth that cell_true_values() scores for the
# posterior intervals. Any failure is caught and recorded; it never kills
# the cell.
run_conformal <- function() {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(impute(cell$df, cell$tree, gnn = FALSE, seed = cell$seed, verbose = FALSE),
                  error = function(e) e)
  wall <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error")) {
    return(list(detail = NULL, wall = wall, err = conditionMessage(res)))
  }
  det <- tryCatch({
    lo <- res$prediction$conformal_lower
    hi <- res$prediction$conformal_upper
    if (is.null(lo) || is.null(hi)) stop("impute() returned no conformal intervals")
    masks <- list(x = cell$mask_x, y = cell$mask_y)
    do.call(rbind, lapply(names(masks), function(tr) {
      sp <- rownames(cell$truth)[masks[[tr]]]
      if (!length(sp)) return(NULL)
      if (!(tr %in% colnames(lo))) stop("no conformal column for trait ", tr)
      idx <- match(sp, rownames(lo))
      if (anyNA(idx)) stop(sum(is.na(idx)), " masked ", tr, " cells missing from the conformal rows")
      true_vals <- cell$truth[sp, tr]
      data.frame(method = "conformal", trait = tr,
                 covered = as.integer(true_vals >= lo[idx, tr] & true_vals <= hi[idx, tr]),
                 width = hi[idx, tr] - lo[idx, tr])
    }))
  }, error = function(e) e)
  if (inherits(det, "error")) {
    return(list(detail = NULL, wall = wall, err = conditionMessage(det)))
  }
  list(detail = det, wall = wall, err = NA_character_)
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
posterior_methods <- mi_gls_v2_methods(cell$regime$missing)   # D6
pu_mode <- if ("posterior_none" %in% posterior_methods) "both" else "full"
post <- run_posterior(pu_mode)
stopifnot(identical(names(post$res), posterior_methods))
sampler_control <- post$ctl
for (method_name in posterior_methods) {
  out <- post$res[[method_name]]
  results[[method_name]] <- out$method_res
  if (!is.null(out$cell_detail))   cell_details[[method_name]]   <- out$cell_detail
  if (!is.null(out$cell_summary))  cell_summaries[[method_name]] <- out$cell_summary
  diagnostics[[method_name]] <- data.frame(
    method = method_name, max_rhat = out$diag$max_rhat,
    min_ess = out$diag$min_ess, converged = out$diag$converged,
    n_extensions = out$diag$n_extensions, sweeps_per_chain = out$diag$sweeps_per_chain)
  posterior_errors[[method_name]] <- out$err
}

# Conformal arm last, so the posterior arms (seeded by cell$seed) are
# unaffected by anything impute() does to the RNG.
conf <- run_conformal()
if (!is.null(conf$detail)) {
  cell_details$conformal   <- conf$detail
  cell_summaries$conformal <- cbind(method = "conformal", summarise_cell_coverage(conf$detail))
}

res_df <- do.call(rbind, results); rownames(res_df) <- NULL
cell_detail_df  <- if (length(cell_details))   do.call(rbind, cell_details)   else NULL
cell_summary_df <- if (length(cell_summaries)) do.call(rbind, cell_summaries) else NULL
diag_df <- do.call(rbind, diagnostics); rownames(diag_df) <- NULL

# ---- assemble + save -----------------------------------------------------------
out <- list(
  regime_id = regime_id, rep = rep_i, seed = cell$seed, regime = cell$regime,
  m = m, keep_draws = keep_draws, n_chains = n_chains, true_beta = true_beta,
  n_iter_override = n_iter_override, burnin_override = burnin_override,
  sampler_control = sampler_control,
  n_missing_x = sum(cell$mask_x), n_missing_y = sum(cell$mask_y),
  posterior_methods = posterior_methods, param_uncertainty_mode = pu_mode,
  posterior_errors = posterior_errors,
  conformal_error = conf$err, conformal_wall_s = conf$wall,
  code_sha = code_sha, code_sha_source = code_sha_source,
  pigauto_version = pigauto_version, r_version = R.version.string,
  pkg_versions = vapply(c("Matrix", "ape", "nlme", "phylolm"), function(p)
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), ""),
  host = unname(Sys.info()[["nodename"]]),
  results = res_df,
  cell_detail = cell_detail_df,
  cell_coverage = cell_summary_df,
  diagnostics = diag_df
)
saveRDS(out, out_f)

cat(sprintf("OK regime=%d rep=%d n=%d mech=%s missing=%s code_sha=%s pigauto=%s\n",
           regime_id, rep_i, cell$regime$n, cell$regime$mechanism, cell$regime$missing,
           code_sha, pigauto_version))
if (!is.na(conf$err)) cat(sprintf("CONFORMAL_FAILED regime=%d rep=%d: %s\n",
                                  regime_id, rep_i, conf$err))
for (nm in names(diagnostics)) {
  d <- diagnostics[[nm]]
  if (!isTRUE(d$converged)) {
    cat(sprintf("NONCONVERGED regime=%d rep=%d method=%s max_rhat=%.4f min_ess=%.1f\n",
               regime_id, rep_i, nm, d$max_rhat, d$min_ess))
  }
}
print(res_df)
