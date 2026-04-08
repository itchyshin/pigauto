# Post-benchmark analysis and reporting for pigauto v2
suppressPackageStartupMessages({ devtools::load_all(quiet = TRUE) })

bench <- readRDS("script/bench_v2.rds")

cat("\n===================================================\n")
cat("pigauto benchmark v2: direct delta training results\n")
cat("===================================================\n\n")

# Summary table
res <- bench$results

# --- Continuous (RMSE) ---------------------------------------------------
cont <- res[res$type == "continuous" & res$metric == "rmse", ]
cat("CONTINUOUS TRAITS (RMSE, lower = better)\n")
cat("----------------------------------------\n")
cat(sprintf("%-15s %-10s %8s %8s %8s\n",
            "scenario", "trait", "baseline", "pigauto", "delta%"))
for (scen in unique(cont$scenario)) {
  for (tr in unique(cont$trait[cont$scenario == scen])) {
    bl <- mean(cont$value[cont$scenario==scen & cont$trait==tr & cont$method=="baseline"])
    pa <- mean(cont$value[cont$scenario==scen & cont$trait==tr & cont$method=="pigauto"])
    delta <- 100 * (bl - pa) / bl
    cat(sprintf("%-15s %-10s %8.4f %8.4f %+7.2f%%\n", scen, tr, bl, pa, delta))
  }
}

# --- Discrete (accuracy) -------------------------------------------------
disc <- res[res$type %in% c("binary", "categorical") & res$metric == "accuracy", ]
if (nrow(disc) > 0) {
  cat("\nDISCRETE TRAITS (accuracy, higher = better)\n")
  cat("-------------------------------------------\n")
  cat(sprintf("%-15s %-15s %-10s %8s %8s %8s\n",
              "scenario", "trait", "type", "baseline", "pigauto", "delta"))
  for (scen in unique(disc$scenario)) {
    for (tr in unique(disc$trait[disc$scenario == scen])) {
      type <- disc$type[disc$scenario==scen & disc$trait==tr][1]
      bl <- mean(disc$value[disc$scenario==scen & disc$trait==tr & disc$method=="baseline"])
      pa <- mean(disc$value[disc$scenario==scen & disc$trait==tr & disc$method=="pigauto"])
      delta <- 100 * (pa - bl)
      cat(sprintf("%-15s %-15s %-10s %8.3f %8.3f %+7.1fpp\n",
                  scen, tr, type, bl, pa, delta))
    }
  }
}

# --- Overall improvement -------------------------------------------------
cat("\n\nOVERALL IMPROVEMENT SUMMARY\n")
cat("===========================\n")
for (scen in unique(res$scenario)) {
  rmse_sub <- res[res$scenario == scen & res$metric == "rmse", ]
  if (nrow(rmse_sub) > 0) {
    bl <- mean(rmse_sub$value[rmse_sub$method == "baseline"])
    pa <- mean(rmse_sub$value[rmse_sub$method == "pigauto"])
    delta <- 100 * (bl - pa) / bl
    cat(sprintf("  %-15s  RMSE  %+6.2f%%   (%.4f -> %.4f)\n",
                scen, delta, bl, pa))
  }

  acc_sub <- res[res$scenario == scen & res$metric == "accuracy", ]
  if (nrow(acc_sub) > 0) {
    bl <- mean(acc_sub$value[acc_sub$method == "baseline"])
    pa <- mean(acc_sub$value[acc_sub$method == "pigauto"])
    delta <- 100 * (pa - bl)
    cat(sprintf("  %-15s  ACC   %+6.2fpp  (%.3f -> %.3f)\n",
                scen, delta, bl, pa))
  }
}

# Save compact CSV
write.csv(bench$results, "script/bench_v2_results.csv", row.names = FALSE)
write.csv(bench$summary, "script/bench_v2_summary.csv", row.names = FALSE)

# Make a plot
pdf("script/bench_v2.pdf", width = 12, height = 5)
plot(bench, metric = "rmse")
dev.off()
cat("\nWrote script/bench_v2.pdf, bench_v2_results.csv, bench_v2_summary.csv\n")
