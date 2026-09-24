# Compare a new 18-core-cell aggregate (pigauto arms, lambda estimated by default) with the committed
# yardstick from arc/imputation-sim. Prints a table and BENCH_OK or BENCH_FAIL.
#
# Usage: Rscript docs/dev-log/lambda-default/compare_benchmark.R [new_summary.csv] [committed_summary.csv]
# Yardstick definition (fixed here so every report cites the same numbers): per (n, lambda, arm), the mean
# zRMSE over the continuous-family traits c1, c2, cnt, prp and over rho in {0, 0.5}. Coverage: c1, c2.
args <- commandArgs(trailingOnly = TRUE)
new_path <- if (length(args) >= 1) args[1] else "docs/dev-log/lambda-default/core_lambda_agg_summary.csv"
old_path <- if (length(args) >= 2) args[2] else
  "/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/script/campaign_sim_results/summary.csv"
old <- read.csv(old_path); new <- read.csv(new_path)
cont <- c("c1", "c2", "cnt", "prp")
agg <- function(s, m, traits) {
  z <- s[s$metric == m & s$trait %in% traits & s$dgp == "types_mixed" & s$miss == "mcar", ]
  aggregate(mean ~ n + lambda + arm, z, mean)
}
zo <- agg(old, "zRMSE", cont); zn <- agg(new, "zRMSE", cont)
co <- agg(old, "coverage", c("c1", "c2")); cn <- agg(new, "coverage", c("c1", "c2"))
arms <- intersect(unique(zn$arm), c("gnn_off", "gnn_off_rphylopars", "gnn_on"))
tab <- merge(zo[zo$arm %in% arms, ], zn, by = c("n", "lambda", "arm"), suffixes = c("_old", "_new"))
ref <- zo[zo$arm == "freq_lambda", c("n", "lambda", "mean")]; names(ref)[3] <- "freq_lambda"
tab <- merge(tab, ref, by = c("n", "lambda"), all.x = TRUE)
tab$delta <- tab$mean_new - tab$mean_old
tab$gap_closed <- ifelse(tab$lambda < 1, (tab$mean_old - tab$mean_new) / (tab$mean_old - tab$freq_lambda), NA)
tab <- tab[order(tab$arm, tab$lambda, tab$n), ]
cat("\n== zRMSE (mean over c1,c2,cnt,prp and rho) ==\n"); print(tab, digits = 3, row.names = FALSE)
ctab <- merge(co[co$arm %in% arms, ], cn, by = c("n", "lambda", "arm"), suffixes = c("_old", "_new"))
ctab$delta <- ctab$mean_new - ctab$mean_old
cat("\n== conformal coverage (c1, c2) ==\n"); print(ctab[order(ctab$arm, ctab$lambda, ctab$n), ], digits = 3, row.names = FALSE)

g <- tab[tab$arm == "gnn_off", ]
guard_l1 <- all(abs(g$delta[g$lambda == 1]) < 0.01)
row03 <- g[g$lambda == 0.3 & g$n == 100, ]
half_gap <- nrow(row03) == 1 && is.finite(row03$gap_closed) && row03$gap_closed >= 0.5
cg <- ctab[ctab$arm == "gnn_off" & ctab$n >= 300, ]
cov_ok <- nrow(cg) > 0 && all(cg$mean_new >= 0.94)   # one-sided: split conformal guarantees >= 0.95 nominal; over-coverage is not a defect
complete <- all(table(zn$arm[zn$arm %in% arms]) == 9)   # 3 lambda x 3 n per arm
cat(sprintf("\nGUARD lambda=1 |delta| < 0.01 : %s\n", guard_l1))
cat(sprintf("GAP lambda 0.3 n 100 closed   : %.2f (>= 0.50 required)\n", row03$gap_closed))
cat(sprintf("COVERAGE n>=300 >= 0.94       : %s\n", cov_ok))
cat(sprintf("COMPLETE 9 cells per arm      : %s\n", complete))
cat(if (guard_l1 && half_gap && cov_ok && complete) "BENCH_OK\n" else "BENCH_FAIL\n")
