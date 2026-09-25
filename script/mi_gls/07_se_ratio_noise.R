#!/usr/bin/env Rscript
# script/mi_gls/07_se_ratio_noise.R
#
# Is the G6 relative SE-ratio rule ([0.90, 1.15] on MI SE ratio / complete-data
# SE ratio, phylolm, gated regimes 17-40) able to separate miscalibration from
# Monte Carlo noise at 200 reps? Reported only; changes no gate.
#
# Each SE ratio is mean SE / empirical SD over R reps; the empirical SD has
# relative SE about 1/sqrt(2(R-1)), and the relative ratio divides two such
# ratios, so its MCSE is about rel * sqrt(2 / (2(R-1))) (treating the two SDs
# as independent, which overstates the noise a little because they share data).
#
# Usage: Rscript script/mi_gls/07_se_ratio_noise.R docs/dev-log/mi-posterior/sim_summary.csv
a <- commandArgs(trailingOnly = TRUE)
s <- utils::read.csv(a[[1L]], stringsAsFactors = FALSE)
g <- s[s$regime_id >= 17 & s$method == "posterior_full" & s$downstream == "phylolm", ]
rel <- g$se_ratio / g$complete_se_ratio
R <- stats::median(g$R)
mcse <- mean(rel) * sqrt(1 / (2 * (R - 1)) + 1 / (2 * (R - 1)))
p_fail <- stats::pnorm(0.90, mean(rel), mcse) + 1 - stats::pnorm(1.15, mean(rel), mcse)
cat(sprintf("gated phylolm rows: %d; reps per row: %d\n", nrow(g), as.integer(R)))
cat(sprintf("relative SE ratio: mean %.3f, SD across regimes %.3f, approx MCSE per row %.3f\n",
            mean(rel), stats::sd(rel), mcse))
cat(sprintf("rows outside [0.90, 1.15]: %d observed; %.1f expected from noise alone (p = %.3f per row)\n",
            sum(rel < 0.90 | rel > 1.15), nrow(g) * p_fail, p_fail))
cat(sprintf("pooled mean relative ratio %.3f, SE %.3f (t = %.1f against 1)\n",
            mean(rel), mcse / sqrt(nrow(g)), (mean(rel) - 1) / (mcse / sqrt(nrow(g)))))
cat(sprintf("absolute MI SE ratio: range %.3f to %.3f; outside [0.90, 1.15]: %d (complete-data ratio range %.3f to %.3f)\n",
            min(g$se_ratio), max(g$se_ratio), sum(g$se_ratio < 0.90 | g$se_ratio > 1.15),
            min(g$complete_se_ratio), max(g$complete_se_ratio)))
o <- data.frame(regime = g$regime_id, mi = round(g$se_ratio, 3), complete = round(g$complete_se_ratio, 3),
                rel = round(rel, 3))
print(o[order(-o$rel), ], row.names = FALSE)
