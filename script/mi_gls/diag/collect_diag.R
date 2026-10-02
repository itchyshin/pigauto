#!/usr/bin/env Rscript
# script/mi_gls/diag/collect_diag.R
#
# Collects the diagnosis outputs (oracle_cell.R, rerun_cell.R) into small
# CSVs for docs/dev-log/mi-posterior/evidence/diagnosis/.
#
# Usage: Rscript collect_diag.R <out_dir> <dest_dir>
#   out_dir  .../diag_<sha>/out (holds oracle/ and rerun/)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1L]]; dest <- args[[2L]]
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

# ---- oracle arms ----------------------------------------------------------------
of <- list.files(file.path(out_dir, "oracle"), pattern = "\\.rds$", full.names = TRUE)
if (length(of)) {
  orows <- lapply(of, function(f) {
    x <- readRDS(f)
    r <- x$results
    r$a_star <- x$kl[["a_star"]]; r$b_star <- x$kl[["b_star"]]
    r$lambda_star <- x$kl[["lambda_star"]]; r$mean_diag_V <- x$kl[["mean_diag_V"]]
    r$min_diag_V <- x$kl[["min_diag_V"]]; r$kl_conv <- x$kl[["kl_conv"]]
    r$n_missing <- x$n_missing; r$pkg_sha <- x$pkg_sha; r$diag_sha <- x$diag_sha
    r
  })
  o <- do.call(rbind, orows)
  o <- o[order(o$regime_id, o$rep, o$arm, o$downstream), ]
  utils::write.csv(o, file.path(dest, "oracle_long.csv"), row.names = FALSE)
  cat("oracle rows:", nrow(o), "from", length(of), "files\n")
}

# ---- sampler re-runs ------------------------------------------------------------
rf <- list.files(file.path(out_dir, "rerun"), pattern = "\\.rds$", full.names = TRUE)
if (length(rf)) {
  dl <- list(); sl <- list()
  for (f in rf) {
    x <- readRDS(f)
    d <- x$diagnostics
    dl[[f]] <- data.frame(tag = x$tag, regime_id = x$regime_id, rep = x$rep,
                          parameter = d$parameter, rhat = d$rhat, ess_bulk = d$ess_bulk)
    ps <- x$posterior_summary
    pm <- stats::setNames(ps$mean, paste0("pmean_", make.names(ps$quantity)))
    pmed <- stats::setNames(ps$median, paste0("pmed_", make.names(ps$quantity)))
    nn <- function(v) if (is.null(v)) NA_real_ else v
    sl[[f]] <- data.frame(
      tag = x$tag, regime_id = x$regime_id, rep = x$rep, seed = x$seed, twin = x$twin,
      n_iter = x$control$n_iter, burnin = x$control$burnin, thin = x$control$thin,
      converged = x$converged, max_rhat = max(d$rhat), min_ess = min(d$ess_bulk),
      rhat_param = d$parameter[which.max(d$rhat)],
      ess_param = d$parameter[which.min(d$ess_bulk)],
      complete_gls = x$complete["gls", "estimate"],
      complete_phylolm = x$complete["phylolm", "estimate"],
      full_gls = x$posterior_full["gls", "estimate"],
      full_phylolm = x$posterior_full["phylolm", "estimate"],
      full_gls_se = x$posterior_full["gls", "se"],
      full_phylolm_se = x$posterior_full["phylolm", "se"],
      none_gls = nn(x$posterior_none["gls", "estimate"]),
      none_phylolm = nn(x$posterior_none["phylolm", "estimate"]),
      t(pm), t(pmed),
      h3_log_transform_any = any(x$h3$log_transform),
      h3_max_obs_change = x$h3$max_obs_change,
      h3_any_na = x$h3$any_na_completed,
      latent_sd_x = x$h3$latent_sd[["x"]], latent_sd_y = x$h3$latent_sd[["y"]],
      wall_s = x$wall_s, pkg_sha = x$pkg_sha, diag_sha = x$diag_sha,
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  s <- do.call(rbind, sl); s <- s[order(s$tag, s$regime_id, s$rep), ]
  d <- do.call(rbind, dl); d <- d[order(d$tag, d$regime_id, d$rep, d$parameter), ]
  utils::write.csv(s, file.path(dest, "rerun_summary.csv"), row.names = FALSE)
  utils::write.csv(d, file.path(dest, "rerun_diag_long.csv"), row.names = FALSE)
  cat("rerun rows:", nrow(s), "\n")
  print(table(s$tag, s$regime_id))
}
