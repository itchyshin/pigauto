#!/usr/bin/env Rscript
# script/mi_gls/07_se_ratio_noise.R
#
# Can Monte Carlo noise at 200 reps explain the G6 relative SE-ratio
# failures ([0.90, 1.15] on MI SE ratio / complete-data SE ratio, phylolm,
# gated regimes 17-40)? Reported only; changes no gate.
#
# MCSE of rel = (mean se_MI / sd_MI) / (mean se_c / sd_c). The two empirical
# SDs come from the same R replicates, and within a replicate the MI and
# complete-data slopes are strongly correlated, so the MCSE must allow for
# the pairing (M2 review, 2026-09-24). Under normal theory
# corr(log s_MI, log s_c) is about rho^2, so
#   var(log rel) ~ (1 - rho^2) / (R - 1),  MCSE(rel) ~ rel * sqrt((1 - rho^2) / (R - 1)).
# rho is recovered per row from committed sim_summary.csv columns (no rep
# files needed): sd_diff = paired_bias_mcse * sqrt(R) is the SD of the
# per-rep difference (MI minus complete), so
#   rho = (emp_sd_MI^2 + emp_sd_c^2 - sd_diff^2) / (2 * emp_sd_MI * emp_sd_c).
# Left out: (1) noise in the two mean-SE numerators (small: the per-rep SEs
# vary little and the MI and complete SEs move together); (2) the
# complete-data SD uses every finite rep, while the MI SD and sd_diff use the
# converged reps (they differ only in a row with fewer than 200 converged).
#
# Check of the formula (M2 skeptic, 2026-09-24), on the only committed
# per-rep slopes, evidence/diagnosis/campaign_cells.csv (first campaign,
# 69670d4, regimes 17-24, phylolm, 200 reps each): a paired bootstrap
# (4,000 resamples of reps, seed 1) of sd(complete)/sd(MI) gave a mean
# relative SE of 0.0494, against 0.0500 from this rho formula (rho =
# cor(complete_phylolm, full_phylolm)) and 0.0709 from the independence
# formula sqrt(1/(R-1)). This script used the independence formula until
# 2026-09-24; it overstated the noise and made the three failures look like
# noise. Per-rep files of the final campaign (regimes 17-40) are on Totoro,
# not in the repo, so the committed method is the rho formula above.
#
# Usage: Rscript script/mi_gls/07_se_ratio_noise.R docs/dev-log/mi-posterior/sim_summary.csv
a <- commandArgs(trailingOnly = TRUE)
s <- utils::read.csv(a[[1L]], stringsAsFactors = FALSE)
band <- c(0.90, 1.15)
g <- s[s$regime_id >= 17 & s$method == "posterior_full" & s$downstream == "phylolm", ]
cr <- s[s$regime_id >= 17 & s$method == "complete" & s$downstream == "phylolm", ]
g <- g[order(g$regime_id), ]
cr <- cr[match(g$regime_id, cr$regime_id), ]
stopifnot(identical(cr$regime_id, g$regime_id), all(is.finite(cr$emp_sd)))
k <- nrow(g)
R <- g$R
rel <- g$se_ratio / g$complete_se_ratio
sd_diff <- g$paired_bias_mcse * sqrt(R)
rho <- pmin(1, (g$emp_sd^2 + cr$emp_sd^2 - sd_diff^2) / (2 * g$emp_sd * cr$emp_sd))
mcse <- rel * sqrt((1 - rho^2) / (R - 1))
mcse_ind <- rel * sqrt(1 / (2 * (R - 1)) + 1 / (2 * (R - 1)))   # the pre-2026-09-24 formula
out <- rel < band[1] | rel > band[2]
m <- mean(rel)
p_fail <- function(centre) stats::pnorm(band[1], centre, mcse) + 1 - stats::pnorm(band[2], centre, mcse)
p_m <- p_fail(m); p_1 <- p_fail(1)
# Poisson-binomial P(X >= x) for independent rows with probabilities p.
pb_upper <- function(p, x) {
  d <- 1
  for (pi in p) d <- c(d * (1 - pi), 0) + c(0, d * pi)
  sum(d[(x + 1L):length(d)])
}
w <- 1 / mcse^2
m_w <- sum(w * rel) / sum(w)
Q <- sum(w * (rel - m_w)^2)
Q_unw <- sum(((rel - m) / mcse)^2)
tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
se_mc <- sqrt(sum(mcse^2)) / k
se_between <- stats::sd(rel) / sqrt(k)

cat(sprintf("gated phylolm rows: %d; reps per row: %d to %d\n", k, min(R), max(R)))
cat(sprintf("MI vs complete-data slope correlation rho (recovered): %.2f to %.2f (median %.2f)\n",
            min(rho), max(rho), stats::median(rho)))
cat(sprintf("MCSE of rel, correlation-aware: %.3f to %.3f (mean %.4f); independence formula: %.3f to %.3f (mean %.4f), %.2f to %.2f times larger (mean %.2f)\n",
            min(mcse), max(mcse), mean(mcse), min(mcse_ind), max(mcse_ind), mean(mcse_ind),
            min(mcse_ind / mcse), max(mcse_ind / mcse), mean(mcse_ind / mcse)))
cat(sprintf("relative SE ratio: mean %.3f, SD across regimes %.3f (mean MCSE %.3f)\n",
            m, stats::sd(rel), mean(mcse)))
cat(sprintf("rows outside [%.2f, %.2f]: %d observed; expected from noise alone: %.2f if every row's true value is the mean %.3f (P(>= %d) = %.3f, Poisson-binomial), %.2f if it is 1 (P(>= %d) = %.3f)\n",
            band[1], band[2], sum(out), sum(p_m), m, sum(out), pb_upper(p_m, sum(out)),
            sum(p_1), sum(out), pb_upper(p_1, sum(out))))
cat(sprintf("heterogeneity: Cochran Q = %.1f on %d df (p = %.3f; inverse-variance mean %.3f); around the unweighted mean %.1f (p = %.3f); DL between-regime SD %.3f\n",
            Q, k - 1L, stats::pchisq(Q, k - 1L, lower.tail = FALSE), m_w,
            Q_unw, stats::pchisq(Q_unw, k - 1L, lower.tail = FALSE), sqrt(tau2)))
blocks <- list(`17-24 (two lambdas)` = 17:24, `25-32 (twins, x only)` = 25:32, `33-40 (twins, both missing)` = 33:40)
cat("heterogeneity by block (Cochran Q, inverse-variance):",
    paste(vapply(names(blocks), function(b) {
      i <- g$regime_id %in% blocks[[b]]
      if (sum(i) < 2L) return(sprintf("%s: n/a", b))
      mb <- sum(w[i] * rel[i]) / sum(w[i]); qb <- sum(w[i] * (rel[i] - mb)^2)
      sprintf("%s: mean %.3f, Q = %.1f on %d df (p = %.3f)", b, mean(rel[i]), qb, sum(i) - 1L,
              stats::pchisq(qb, sum(i) - 1L, lower.tail = FALSE))
    }, character(1)), collapse = "; "), "\n")
cat(sprintf("pooled mean relative ratio %.3f: SE %.4f from the per-row MCSEs (t = %.1f against 1); SE %.4f from the between-regime SD (t = %.1f; 95%% CI %.3f to %.3f)\n",
            m, se_mc, (m - 1) / se_mc, se_between, (m - 1) / se_between,
            m - stats::qt(0.975, k - 1L) * se_between, m + stats::qt(0.975, k - 1L) * se_between))
abs_out <- g$se_ratio < band[1] | g$se_ratio > band[2]
cat(sprintf("rows failing the relative rule: %s; absolute MI SE ratio there %.3f to %.3f; complete-data ratio there %.3f to %.3f\n",
            paste(g$regime_id[out], collapse = ", "), min(g$se_ratio[out]), max(g$se_ratio[out]),
            min(g$complete_se_ratio[out]), max(g$complete_se_ratio[out])))
cat(sprintf("absolute MI SE ratio: range %.3f to %.3f over all rows; outside [%.2f, %.2f]: %d (regimes %s), where the complete-data ratio is %.3f to %.3f (all rows: %.3f to %.3f)\n",
            min(g$se_ratio), max(g$se_ratio), band[1], band[2], sum(abs_out),
            paste(g$regime_id[abs_out], collapse = ", "),
            min(g$complete_se_ratio[abs_out]), max(g$complete_se_ratio[abs_out]),
            min(g$complete_se_ratio), max(g$complete_se_ratio)))

# Consequences of candidate decisions on the current numbers (reported; no
# gate is changed here).
c_bad <- g$complete_se_ratio < band[1] | g$complete_se_ratio > band[2]
opt_b <- ifelse(c_bad, abs_out, out)
opt_either <- out & abs_out
opt_mcse <- rel < 1 - 2.5 * mcse | rel > 1 + 2.5 * mcse
ids <- function(x) if (any(x)) paste(g$regime_id[x], collapse = ", ") else "none"
cat("OPTIONS (current numbers; no gate changed):\n")
cat(sprintf("  keep the relative rule: %d row(s) fail (%s)\n", sum(out), ids(out)))
cat(sprintf("  absolute MI ratio in the %d rows whose complete-data ratio is outside the band (%s), relative elsewhere: %d row(s) fail (%s)\n",
            sum(c_bad), ids(c_bad), sum(opt_b), ids(opt_b)))
cat(sprintf("  pass if either the relative or the absolute ratio is in the band: %d row(s) fail (%s)\n",
            sum(opt_either), ids(opt_either)))
cat(sprintf("  relative ratio within 1 +/- 2.5 correlation-aware MCSE per row: %d row(s) fail (%s)\n",
            sum(opt_mcse), ids(opt_mcse)))
cat(sprintf("  pooled mean relative ratio in [0.95, 1.10]: %.3f, %s\n", m,
            if (m >= 0.95 && m <= 1.10) "passes" else "fails"))

o <- data.frame(regime = g$regime_id, mi = round(g$se_ratio, 3), complete = round(g$complete_se_ratio, 3),
                rel = round(rel, 3), rho = round(rho, 2), mcse = round(mcse, 3),
                mcse_indep = round(mcse_ind, 3), z_mean = round((rel - m) / mcse, 2),
                z_one = round((rel - 1) / mcse, 2), p_out = round(p_m, 3),
                out = ifelse(out, "FAIL", ""))
print(o[order(-o$rel), ], row.names = FALSE)
