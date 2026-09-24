#!/usr/bin/env Rscript
# script/mi_gls/06_results_tables.R
#
# Results tables for docs/dev-log/mi-posterior/results.md, built only from the
# summary CSVs written by 03_summarise_v2.R (so every number in the results
# document traces to a file). Prints markdown.
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
  sprintf("%d (lambda %s%s, n %d, %s, %s)", id, lam, cor, r$n, r$mechanism, r$missing)
}

lines <- c("## Downstream slope (posterior_full vs complete data, same analysis model)", "",
           "Paired bias = mean over reps of (MI pooled slope - complete-data slope in the same rep); MCSE in brackets.",
           "SE ratio = mean pooled SE / empirical SD of the pooled slope; 'rel' = MI ratio / complete ratio (gated under phylolm, [0.90, 1.15]).",
           "Coverage truth: 0.7 in regimes 1-16; mean complete-data slope in 17-24.", "")
pf <- s[s$method == "posterior_full", ]
pn <- s[s$method == "posterior_none", ]
rows <- lapply(seq_len(nrow(pf)), function(i) {
  x <- pf[i, ]
  n <- pn[pn$regime_id == x$regime_id & pn$downstream == x$downstream, ]
  data.frame(
    regime = reg_label(x$regime_id), analysis = x$downstream,
    paired_bias = sprintf("%s (%s)", f3(x$paired_bias), f3(x$paired_bias_mcse)),
    se_ratio = f2(x$se_ratio), complete_se_ratio = f2(x$complete_se_ratio),
    rel = f2(x$se_ratio / x$complete_se_ratio),
    coverage = f3(x$coverage), complete_coverage = f3(x$complete_coverage),
    plugin_se_ratio = if (nrow(n)) f2(n$se_ratio) else "",
    plugin_coverage = if (nrow(n)) f3(n$coverage) else "",
    converged = sprintf("%d/%d", x$n_converged, x$n_expected),
    stringsAsFactors = FALSE)
})
lines <- c(lines, md_table(do.call(rbind, rows)), "")

lines <- c(lines, "## Per-cell 95% predictive coverage (masked truth)", "",
           "Coverage = covered cells / scored cells over all converged reps; width = mean interval width on the latent scale.",
           "posterior_full vs conformal (impute(gnn = FALSE), same masked cells). MCAR is gated in [0.92, 0.98]; MAR_phylo (clade-biased) is reported only.", "")
cc$coverage <- cc$covered_sum / cc$n
pc <- cc[cc$method == "posterior_full", ]
rows <- lapply(seq_len(nrow(pc)), function(i) {
  x <- pc[i, ]
  k <- cc[cc$method == "conformal" & cc$regime_id == x$regime_id & cc$trait == x$trait, ]
  data.frame(regime = reg_label(x$regime_id), trait = x$trait, mask = x$mechanism,
             cells = x$n, posterior_cov = f3(x$coverage), posterior_width = f2(x$mean_width),
             conformal_cov = if (nrow(k)) f3(k$coverage) else "",
             conformal_width = if (nrow(k)) f2(k$mean_width) else "",
             width_ratio = if (nrow(k)) f2(x$mean_width / k$mean_width) else "",
             stringsAsFactors = FALSE)
})
lines <- c(lines, md_table(do.call(rbind, rows)), "")

conv <- pf[pf$downstream == "gls", ]
lines <- c(lines, "## Convergence (posterior_full)", "",
           sprintf("Fits: %d; converged: %d (%.2f%%). Median of max R-hat per regime: %s to %s; median of min bulk ESS: %s to %s.",
                   sum(conv$n_fits), sum(conv$n_converged), 100 * sum(conv$n_converged) / sum(conv$n_fits),
                   f3(min(conv$median_max_rhat)), f3(max(conv$median_max_rhat)),
                   formatC(min(conv$median_min_ess), format = "d", big.mark = ","),
                   formatC(max(conv$median_min_ess), format = "d", big.mark = ",")), "")
sha <- unique(s$code_sha)
lines <- c(lines, sprintf("Source: `%s`, `%s`; code SHA %s.", args[[1L]], args[[2L]], paste(sha, collapse = ",")))
if (nzchar(out)) writeLines(lines, out) else cat(lines, sep = "\n")
