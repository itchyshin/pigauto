#!/usr/bin/env Rscript
# script/mi_gls/gate_convergence.R
#
# G4 smoke (.unlazy/mi-posterior/GATES.md): split R-hat < 1.05 and bulk
# ESS > 400 for every element of Sigma_P, Sigma_E and lambda_k, on
# regime 1 rep 1 (single-shared-lambda tree-transform DGP, missing =
# "x_only") and regime 23 rep 1 (two-lambda Sigma_P/Sigma_E DGP,
# MAR_phylo, both traits missing) from multi_impute(...,
# draws_method = "posterior", posterior_control = list(param_uncertainty =
# "full")). Reuses the cell's own data generator
# (script/mi_gls/dgp_v2.R::simulate_regime_cell()) so this checks the exact
# regimes the acceptance sweep (04_acceptance.R / G6) will run.
#
# Usage: Rscript script/mi_gls/gate_convergence.R
#
# Env vars: MI_POST_CHAINS (default 4), MI_POST_KEEP (default 1000).
#
# The printed wall_s is ONE posterior_full fit at m = 5 on two n = 300
# regimes. It is NOT a campaign sizing source: a campaign cell also runs
# the complete-data fits, posterior_none (both-missing regimes) and the
# conformal arm, and half the regimes are n = 1000. Size campaigns from a
# timed 01_cell_v2.R cell (see 11_fir_array.sbatch).

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({ devtools::load_all(quiet = TRUE); library(ape) })

source(file.path("script", "mi_gls", "dgp_v2.R"))

draws_methods <- tryCatch(eval(formals(multi_impute)$draws_method),
                          error = function(e) character(0))
if (!("posterior" %in% draws_methods)) {
  stop("multi_impute() does not yet support draws_method = \"posterior\" ",
      "(R/mi_posterior.R and the multi_impute() wiring have not landed on ",
      "this branch/worktree yet). See docs/dev-log/mi-posterior/design.md ",
      "section 4 for the frozen API this gate codes against. G4 cannot run ",
      "until that lands.", call. = FALSE)
}

n_chains <- as.integer(Sys.getenv("MI_POST_CHAINS", "4"))
keep     <- as.integer(Sys.getenv("MI_POST_KEEP", "1000"))

check_one <- function(regime_id, rep_i) {
  cell <- simulate_regime_cell(regime_id, rep_i)
  t0 <- proc.time()[["elapsed"]]
  mi <- multi_impute(cell$df, cell$tree, m = 5L, draws_method = "posterior",
                     posterior_control = list(n_chains = n_chains, keep_draws = keep,
                                              param_uncertainty = "full",
                                              seed = cell$seed))
  wall_s <- proc.time()[["elapsed"]] - t0
  diag <- mi$posterior$diagnostics
  if (is.null(diag) || !nrow(diag)) {
    stop("regime ", regime_id, " rep ", rep_i,
        ": mi$posterior$diagnostics is empty or NULL", call. = FALSE)
  }
  list(regime_id = regime_id, rep = rep_i,
      max_rhat  = max(diag$rhat, na.rm = TRUE),
      min_ess   = min(diag$ess_bulk, na.rm = TRUE),
      converged = isTRUE(attr(diag, "converged")), wall_s = wall_s,
      n_extensions = as.integer(attr(diag, "n_extensions") %||% NA_integer_))
}

results <- tryCatch(
  list(check_one(1L, 1L), check_one(23L, 1L)),
  error = function(e) {
    cat("G4 ERROR:", conditionMessage(e), "\n")
    NULL
  })

if (is.null(results)) quit(save = "no", status = 1L)

ok <- TRUE
for (r in results) {
  cat(sprintf(
    "regime %d rep %d: n_chains=%d max_rhat=%.4f min_ess=%.1f converged=%s n_extensions=%s wall_s=%.1f\n",
    r$regime_id, r$rep, n_chains, r$max_rhat, r$min_ess, r$converged, format(r$n_extensions), r$wall_s))
  if (!isTRUE(r$converged)) {
    cat(sprintf("NONCONVERGED regime %d: attr(diagnostics, \"converged\") is not TRUE\n",
               r$regime_id))
  }
  if (!(is.finite(r$max_rhat) && r$max_rhat < 1.05 &&
       is.finite(r$min_ess) && r$min_ess > 400 && isTRUE(r$converged))) {
    ok <- FALSE
  }
}

if (ok) cat("CONVERGENCE_OK\n") else cat("G4 FAILED: rhat/ess/converged thresholds not met\n")
quit(save = "no", status = if (ok) 0L else 1L)
