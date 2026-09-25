#!/usr/bin/env Rscript
# script/mi_gls/06_results_tables.R
#
# Results tables for docs/dev-log/mi-posterior/results.md, built only from the
# summary CSVs written by 03_summarise_v2.R (so every number in the results
# document traces to a file). Prints markdown.
#
# Rows are grouped by DGP (regimes.R column `dgp`; CP2 follow-up, Shinichi
# 2026-09-24): the in-model regimes first (twins 25-40, then Kronecker
# 17-24; both gated), then the stress test (regimes 1-16, raw tree
# covariance outside the sampler's model; reported, not gated). In the
# downstream table the twin rows carry the source regime's paired bias and
# coverage beside their own (source_paired_bias, source_coverage).
#
# Usage:
#   Rscript script/mi_gls/06_results_tables.R <sim_summary.csv> <cell_coverage.csv> [out.md]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected: sim_summary.csv cell_coverage.csv [out.md]", call. = FALSE)
s <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE)
cc <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE)
out <- if (length(args) >= 3L) args[[3L]] else ""
source(file.path("script", "mi_gls", "regimes.R"))

f3 <- function(x) ifelse(is.finite(x), formatC(x, digits = 3, format = "f"), "NA")
f2 <- function(x) ifelse(is.finite(x), formatC(x, digits = 2, format = "f"), "NA")
md_table <- function(d) {
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
reg_label <- function(id) {
  r <- regimes[regimes$regime_id == id, ]
  lam <- if (is.na(r$lambda)) sprintf("%.1f/%.1f", r$lambda_x, r$lambda_y) else sprintf("%.1f", r$lambda)
  cor <- if (is.na(r$corr_phylo)) "" else sprintf(" rP %.1f rE %.1f", r$corr_phylo, r$corr_resid)
  twin <- if (identical(r$dgp, "twin")) sprintf("twin of %d; ", r$twin_of) else ""
  sprintf("%d (%slambda %s%s, n %d, %s, %s)", id, twin, lam, cor, r$n, r$mechanism, r$missing)
}
dgp_of <- function(id) regimes$dgp[match(id, regimes$regime_id)]

# In-model regimes first, then the stress test.
groups <- list(
  twin      = "In-model twins of regimes 1-16 (regimes 25-40; gated)",
  kronecker = "Two-lambda Kronecker DGP (regimes 17-24; gated)",
  tree_raw  = "Stress test: raw tree covariance, outside the sampler's model (regimes 1-16; reported, not gated)")

lines <- c("## Downstream slope (posterior_full vs complete data, same analysis model)", "",
           "Paired bias = mean over reps of (MI pooled slope - complete-data slope in the same rep); MCSE in brackets.",
           "SE ratio = mean pooled SE / empirical SD of the pooled slope; 'rel' = MI ratio / complete ratio (gated under phylolm, [0.90, 1.15]).",
           "rel_mcse: approximate Monte Carlo SE of 'rel', allowing for the pairing of the two empirical SDs over the same reps: rel * sqrt((1 - rho^2) / (R - 1)), where rho (column rho) is the per-rep correlation of the MI and complete-data slopes, recovered from emp_sd of both rows and paired_bias_mcse (see 07_se_ratio_noise.R); reported only, not part of any gate.",
           "Coverage truth: 0.7 in regimes 1-16 and 25-40; mean complete-data slope in 17-24.",
           "Twin rows: source_paired_bias and source_coverage are the source regime's (twin_of) values under the same analysis model.", "")
pf <- s[s$method == "posterior_full", ]
pn <- s[s$method == "posterior_none", ]
pc <- s[s$method == "complete", ]
pf <- pf[order(pf$regime_id, pf$downstream), ]
down_row <- function(x) {
  n <- pn[pn$regime_id == x$regime_id & pn$downstream == x$downstream, ]
  cm <- pc[pc$regime_id == x$regime_id & pc$downstream == x$downstream, ]
  sd_c <- if (nrow(cm) == 1L) cm$emp_sd else NA_real_
  rho <- min(1, (x$emp_sd^2 + sd_c^2 - x$paired_bias_mcse^2 * x$R) / (2 * x$emp_sd * sd_c))
  rel <- x$se_ratio / x$complete_se_ratio
  d <- data.frame(
    regime = reg_label(x$regime_id), analysis = x$downstream,
    paired_bias = sprintf("%s (%s)", f3(x$paired_bias), f3(x$paired_bias_mcse)),
    se_ratio = f2(x$se_ratio), complete_se_ratio = f2(x$complete_se_ratio),
    rel = f2(rel), rho = f2(rho),
    rel_mcse = f3(rel * sqrt((1 - rho^2) / (x$R - 1))),
    coverage = f3(x$coverage), complete_coverage = f3(x$complete_coverage),
    stringsAsFactors = FALSE)
  if (identical(dgp_of(x$regime_id), "twin")) {
    src <- pf[pf$regime_id == regimes$twin_of[regimes$regime_id == x$regime_id] &
              pf$downstream == x$downstream, ]
    d$source_paired_bias <- if (nrow(src)) sprintf("%s (%s)", f3(src$paired_bias), f3(src$paired_bias_mcse)) else "not in CSV"
    d$source_coverage <- if (nrow(src)) f3(src$coverage) else "not in CSV"
  }
  d$plugin_se_ratio <- if (nrow(n)) f2(n$se_ratio) else ""
  d$plugin_coverage <- if (nrow(n)) f3(n$coverage) else ""
  d$converged <- sprintf("%d/%d", x$n_converged, x$n_expected)
  d
}
for (g in names(groups)) {
  sub <- pf[dgp_of(pf$regime_id) == g, ]
  if (!nrow(sub)) next
  rows <- lapply(seq_len(nrow(sub)), function(i) down_row(sub[i, ]))
  lines <- c(lines, sprintf("### %s", groups[[g]]), "", md_table(do.call(rbind, rows)), "")
}

lines <- c(lines, "## Per-cell 95% predictive coverage (masked truth)", "",
           "Coverage = covered cells / scored cells over all converged reps; width = mean interval width on the latent scale.",
           "posterior_full vs conformal (impute(gnn = FALSE), same masked cells). MCAR is gated in [0.92, 0.98] in regimes 17-40; regimes 1-16 are the stress test (not gated); MAR_phylo (clade-biased) is reported only.", "")
cc$coverage <- cc$covered_sum / cc$n
pc <- cc[cc$method == "posterior_full", ]
pc <- pc[order(pc$regime_id, pc$trait), ]
for (g in names(groups)) {
  sub <- pc[dgp_of(pc$regime_id) == g, ]
  if (!nrow(sub)) next
  rows <- lapply(seq_len(nrow(sub)), function(i) {
    x <- sub[i, ]
    k <- cc[cc$method == "conformal" & cc$regime_id == x$regime_id & cc$trait == x$trait, ]
    data.frame(regime = reg_label(x$regime_id), trait = x$trait, mask = x$mechanism,
               cells = x$n, posterior_cov = f3(x$coverage), posterior_width = f2(x$mean_width),
               conformal_cov = if (nrow(k)) f3(k$coverage) else "",
               conformal_width = if (nrow(k)) f2(k$mean_width) else "",
               width_ratio = if (nrow(k)) f2(x$mean_width / k$mean_width) else "",
               stringsAsFactors = FALSE)
  })
  lines <- c(lines, sprintf("### %s", groups[[g]]), "", md_table(do.call(rbind, rows)), "")
}

conv <- pf[pf$downstream == "gls", ]
lines <- c(lines, "## Convergence (posterior_full)", "")
rng <- function(x, f) { x <- x[is.finite(x)]; if (length(x)) c(f(min(x)), f(max(x))) else c("NA", "NA") }
fmt_ess <- function(x) formatC(x, format = "d", big.mark = ",")
for (g in names(groups)) {
  cg <- conv[dgp_of(conv$regime_id) == g, ]
  if (!nrow(cg)) next
  rh <- rng(cg$median_max_rhat, f3); es <- rng(round(cg$median_min_ess), fmt_ess)
  lines <- c(lines, sprintf("%s. Fits: %d; converged: %d (%.2f%%). Median of max R-hat per regime: %s to %s; median of min bulk ESS: %s to %s.",
                            groups[[g]], sum(cg$n_fits), sum(cg$n_converged),
                            100 * sum(cg$n_converged) / sum(cg$n_fits), rh[1], rh[2], es[1], es[2]), "")
}
sha <- unique(s$code_sha)
lines <- c(lines, sprintf("Source: `%s`, `%s`; code SHA %s.", args[[1L]], args[[2L]], paste(sha, collapse = ",")))
if (nzchar(out)) writeLines(lines, out) else cat(lines, sep = "\n")
