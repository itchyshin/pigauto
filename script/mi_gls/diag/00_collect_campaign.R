#!/usr/bin/env Rscript
# script/mi_gls/diag/00_collect_campaign.R
#
# Diagnosis of the posterior-MI campaign (docs/dev-log/mi-posterior/diagnosis.md).
# Reads every per-cell output of the frozen campaign (code SHA 69670d44f9)
# and writes one row per (regime, rep): complete-data and MI slopes for both
# analysis models, the plug-in (posterior_none) slopes where present, the
# saved convergence summaries (max split R-hat, min bulk ESS, converged
# flag) and the sampler wall time. No model is refitted.
#
# Usage: Rscript 00_collect_campaign.R <campaign_sim_dir> <out_csv>

args <- commandArgs(trailingOnly = TRUE)
sim_dir <- args[[1L]]
out_csv <- args[[2L]]

files <- list.files(sim_dir, pattern = "^regime_[0-9]+_rep_[0-9]+\\.rds$",
                    full.names = TRUE)
get_est <- function(res, meth, ds, col = "estimate") {
  v <- res[[col]][res$method == meth & res$downstream == ds]
  if (length(v)) v[[1L]] else NA_real_
}
rows <- lapply(files, function(f) {
  x <- readRDS(f)
  r <- x$results
  d <- x$diagnostics
  dfull <- d[d$method == "posterior_full", , drop = FALSE]
  data.frame(
    regime_id = x$regime_id, rep = x$rep, seed = x$seed,
    n_missing_x = x$n_missing_x, n_missing_y = x$n_missing_y,
    complete_gls = get_est(r, "complete", "gls"),
    complete_phylolm = get_est(r, "complete", "phylolm"),
    full_gls = get_est(r, "posterior_full", "gls"),
    full_phylolm = get_est(r, "posterior_full", "phylolm"),
    none_gls = get_est(r, "posterior_none", "gls"),
    none_phylolm = get_est(r, "posterior_none", "phylolm"),
    full_wall_s = get_est(r, "posterior_full", "gls", "wall_s"),
    max_rhat = if (nrow(dfull)) dfull$max_rhat[1L] else NA_real_,
    min_ess = if (nrow(dfull)) dfull$min_ess[1L] else NA_real_,
    converged = if (nrow(dfull)) dfull$converged[1L] else NA,
    code_sha = x$code_sha,
    stringsAsFactors = FALSE)
})
out <- do.call(rbind, rows)
out <- out[order(out$regime_id, out$rep), ]
utils::write.csv(out, out_csv, row.names = FALSE)
cat("wrote", nrow(out), "rows to", out_csv, "\n")
