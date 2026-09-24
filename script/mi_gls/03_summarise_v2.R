#!/usr/bin/env Rscript
# script/mi_gls/03_summarise_v2.R
#
# Summarise script/mi_gls/01_cell_v2.R output into:
#   (a) a regime x method x downstream table (paired bias vs complete and
#       its MCSE, SE ratio, coverage and MCSE, fit-failure rate, plus
#       per-regime x method convergence summary: n_fits, n_converged,
#       median max R-hat, median min ESS) -- read by 04_acceptance.R (G6).
#   (b) a regime x method x trait x mechanism per-cell coverage table
#       (n, covered_sum, plus the same convergence summary) -- read by
#       05_cell_coverage.R (G7).
#
# Design-review item B3 (2026-09-24): the statistical accuracy rules
# (paired bias, SE ratio, coverage) are computed from CONVERGED reps only
# for posterior_full / posterior_none; "complete" has no convergence
# concept. fit_failure_rate is computed over ALL attempted reps (a fit
# that errors is a different failure mode than one that converges poorly).
# The gate scripts separately flag regimes with > 2% non-converged fits.
#
# Usage:
#   Rscript script/mi_gls/03_summarise_v2.R <indir> <out_summary.csv> [out.md] [out_cell_coverage.csv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected: indir out_summary.csv [out.md] [out_cell_coverage.csv]", call. = FALSE)
indir       <- args[[1L]]
out_csv     <- args[[2L]]
out_md      <- if (length(args) >= 3L) args[[3L]] else NA_character_
out_cell_csv <- if (length(args) >= 4L) args[[4L]] else NA_character_

source(file.path("script", "mi_gls", "regimes.R"))

files <- list.files(indir, pattern = "^regime_[0-9]+_rep_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0L) stop("no .rds files found in ", indir, call. = FALSE)

cells <- list(); diags <- list(); cell_dets <- list()
for (f in files) {
  r <- readRDS(f)
  d <- r$results; d$regime_id <- r$regime_id; d$rep <- r$rep
  cells[[length(cells) + 1L]] <- d
  if (!is.null(r$diagnostics) && nrow(r$diagnostics)) {
    dg <- r$diagnostics; dg$regime_id <- r$regime_id; dg$rep <- r$rep
    diags[[length(diags) + 1L]] <- dg
  }
  if (!is.null(r$cell_detail) && nrow(r$cell_detail)) {
    cd <- r$cell_detail; cd$regime_id <- r$regime_id; cd$rep <- r$rep
    cell_dets[[length(cell_dets) + 1L]] <- cd
  }
}
all_cells <- do.call(rbind, cells)
all_diag  <- if (length(diags))     do.call(rbind, diags)     else NULL
all_cellD <- if (length(cell_dets)) do.call(rbind, cell_dets) else NULL

converged_reps <- function(regime_id, method) {
  if (is.null(all_diag)) return(integer(0))
  sub <- all_diag[all_diag$regime_id == regime_id & all_diag$method == method &
                  (all_diag$converged %in% TRUE), ]
  sub$rep
}

# ---- (a) regime x method x downstream summary --------------------------------
pair_and_summarise <- function(sub, complete_sub) {
  m <- merge(sub, complete_sub[, c("rep", "estimate", "se")],
            by = "rep", suffixes = c("", "_complete"))
  ok <- is.finite(m$estimate) & is.finite(m$estimate_complete)
  R  <- sum(ok)
  paired <- m$estimate[ok] - m$estimate_complete[ok]
  bias      <- if (R > 0) mean(paired) else NA_real_
  bias_mcse <- if (R > 1) stats::sd(paired) / sqrt(R) else NA_real_
  emp_sd  <- if (sum(is.finite(m$estimate)) > 1) stats::sd(m$estimate[is.finite(m$estimate)]) else NA_real_
  mean_se <- if (any(is.finite(m$se))) mean(m$se, na.rm = TRUE) else NA_real_
  se_ratio <- if (is.finite(mean_se) && is.finite(emp_sd) && emp_sd > 0) mean_se / emp_sd else NA_real_
  cov_ok <- is.finite(m$covered)
  coverage <- if (any(cov_ok)) mean(m$covered[cov_ok]) else NA_real_
  coverage_mcse <- if (any(cov_ok) && is.finite(coverage)) {
    sqrt(coverage * (1 - coverage) / sum(cov_ok))
  } else NA_real_
  complete_cov_ok <- is.finite(complete_sub$covered)
  complete_coverage <- if (any(complete_cov_ok)) mean(complete_sub$covered[complete_cov_ok]) else NA_real_
  data.frame(R = R, paired_bias = bias, paired_bias_mcse = bias_mcse,
            emp_sd = emp_sd, mean_se = mean_se, se_ratio = se_ratio,
            coverage = coverage, coverage_mcse = coverage_mcse,
            complete_coverage = complete_coverage)
}

keys <- unique(all_cells[all_cells$method != "complete", c("regime_id", "method", "downstream")])
rows <- lapply(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]
  sub_all <- all_cells[all_cells$regime_id == k$regime_id & all_cells$method == k$method &
                       all_cells$downstream == k$downstream, ]
  n_total <- nrow(sub_all)
  fit_failure_rate <- if (n_total > 0) mean(!is.finite(sub_all$estimate)) else NA_real_

  sub <- sub_all
  # Convergence filtering applies to BOTH posterior_full and
  # posterior_none: R/mi_posterior.R computes diagnostics from the SAME
  # MCMC chains before branching on param_uncertainty (posterior_none just
  # plugs in those chains' posterior mean instead of a fresh draw per
  # completed dataset), so a chain that failed to converge invalidates the
  # posterior_none "fixed" estimate exactly as much as posterior_full's.
  if (k$method %in% c("posterior_full", "posterior_none")) {
    keep <- converged_reps(k$regime_id, k$method)
    sub <- sub[sub$rep %in% keep, , drop = FALSE]
  }
  complete_sub <- all_cells[all_cells$regime_id == k$regime_id & all_cells$method == "complete" &
                            all_cells$downstream == k$downstream, ]
  s <- pair_and_summarise(sub, complete_sub)
  cbind(k, n_total = n_total, fit_failure_rate = fit_failure_rate, s)
})
summary_df <- do.call(rbind, rows)

# convergence summary per regime x method
conv_summary <- function(regime_id, method) {
  if (is.null(all_diag)) {
    return(data.frame(n_fits = NA_integer_, n_converged = NA_integer_,
                      median_max_rhat = NA_real_, median_min_ess = NA_real_))
  }
  sub <- all_diag[all_diag$regime_id == regime_id & all_diag$method == method, ]
  data.frame(n_fits = nrow(sub), n_converged = sum(sub$converged %in% TRUE),
            median_max_rhat = stats::median(sub$max_rhat, na.rm = TRUE),
            median_min_ess  = stats::median(sub$min_ess, na.rm = TRUE))
}
conv_rows <- lapply(seq_len(nrow(summary_df)), function(i) {
  if (!(summary_df$method[i] %in% c("posterior_full", "posterior_none"))) {
    return(data.frame(n_fits = NA_integer_, n_converged = NA_integer_,
                      median_max_rhat = NA_real_, median_min_ess = NA_real_))
  }
  conv_summary(summary_df$regime_id[i], summary_df$method[i])
})
summary_df <- cbind(summary_df, do.call(rbind, conv_rows))
summary_df <- summary_df[order(summary_df$regime_id, summary_df$method, summary_df$downstream), ]
rownames(summary_df) <- NULL

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary_df, out_csv, row.names = FALSE)
cat("wrote", out_csv, "(", nrow(summary_df), "rows )\n")

if (!is.na(out_md)) {
  fmt <- function(x, digits = 4) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  md <- c("# Posterior MI simulation summary", "",
         sprintf("Generated from %d rep file(s) in `%s`.", length(files), indir), "",
         paste0("| ", paste(names(summary_df), collapse = " | "), " |"),
         paste0("|", paste(rep("---", ncol(summary_df)), collapse = "|"), "|"))
  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    md <- c(md, paste0("| ", paste(vapply(r, function(v)
      if (is.numeric(v)) fmt(v) else as.character(v), ""), collapse = " | "), " |"))
  }
  dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)
  writeLines(md, out_md)
  cat("wrote", out_md, "\n")
}

# ---- (b) per-cell coverage table (regime x method x trait x mechanism) -------
if (!is.na(out_cell_csv)) {
  if (is.null(all_cellD)) {
    cell_cov_df <- data.frame(regime_id = integer(0), method = character(0),
                              trait = character(0), mechanism = character(0),
                              n = integer(0), covered_sum = integer(0),
                              mean_width = numeric(0),
                              n_fits = integer(0), n_converged = integer(0),
                              median_max_rhat = numeric(0), median_min_ess = numeric(0))
  } else {
    ck <- unique(all_cellD[, c("regime_id", "method", "trait")])
    cell_rows <- lapply(seq_len(nrow(ck)), function(i) {
      k <- ck[i, ]
      # See the pairing loop above: convergence filtering applies to both
      # posterior_full and posterior_none (shared diagnostics from the same
      # MCMC chains; R/mi_posterior.R computes them before branching on
      # param_uncertainty).
      keep <- if (k$method %in% c("posterior_full", "posterior_none")) {
        converged_reps(k$regime_id, k$method)
      } else {
        unique(all_cellD$rep[all_cellD$regime_id == k$regime_id])
      }
      sub <- all_cellD[all_cellD$regime_id == k$regime_id & all_cellD$method == k$method &
                       all_cellD$trait == k$trait & all_cellD$rep %in% keep, ]
      reg_row <- regimes[regimes$regime_id == k$regime_id, ]
      cv <- if (k$method %in% c("posterior_full", "posterior_none")) {
        conv_summary(k$regime_id, k$method)
      } else {
        data.frame(n_fits = NA_integer_, n_converged = NA_integer_,
                  median_max_rhat = NA_real_, median_min_ess = NA_real_)
      }
      cbind(k, mechanism = reg_row$mechanism, n = nrow(sub),
           covered_sum = sum(sub$covered), mean_width = mean(sub$width), cv)
    })
    cell_cov_df <- do.call(rbind, cell_rows)
    cell_cov_df <- cell_cov_df[order(cell_cov_df$regime_id, cell_cov_df$method, cell_cov_df$trait), ]
    rownames(cell_cov_df) <- NULL
  }
  dir.create(dirname(out_cell_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cell_cov_df, out_cell_csv, row.names = FALSE)
  cat("wrote", out_cell_csv, "(", nrow(cell_cov_df), "rows )\n")
}
