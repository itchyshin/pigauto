#!/usr/bin/env Rscript
# script/analyze_sim_bace.R
#
# Aggregate sim_bace bench results, generate verdict + figures.
# Usage:
#   Rscript script/analyze_sim_bace.R          # default tier=smoke
#   PIGAUTO_TIER=medium Rscript script/analyze_sim_bace.R

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required for figures")
  library(ggplot2)
})

TIER <- tolower(Sys.getenv("PIGAUTO_TIER", "smoke"))
HERE <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(HERE, "script", sprintf("bench_sim_bace_pigauto_%s.rds", TIER))
out_md  <- file.path(HERE, "useful", sprintf("sim_bace_verdict_%s.md", TIER))
out_png <- file.path(HERE, "useful", sprintf("sim_bace_heatmap_%s.png", TIER))
out_pdf <- file.path(HERE, "useful", sprintf("sim_bace_heatmap_%s.pdf", TIER))

if (!file.exists(in_rds)) {
  stop("results .rds not found: ", in_rds, "\nRun bench_sim_bace_pigauto.R first.")
}

df <- readRDS(in_rds)
df$method <- factor(df$method,
  levels = c("column_mean", "lm", "phylolm_lambda_blup", "pigauto_cov_sfT"))

# --- Per-cell aggregation: mean ± SE per (cell-config, method) ----------------
keys <- c("response_type", "n_predictors", "n_species", "multi_obs",
           "phylo_signal", "beta_strength", "method")
agg_one <- function(d) {
  if (d$response_type[1] == "gaussian") {
    list(mean_score = mean(d$rmse, na.rm = TRUE),
          se_score   = stats::sd(d$rmse, na.rm = TRUE) / sqrt(sum(!is.na(d$rmse))),
          score_metric = "rmse")
  } else {
    list(mean_score = mean(d$acc, na.rm = TRUE),
          se_score   = stats::sd(d$acc, na.rm = TRUE) / sqrt(sum(!is.na(d$acc))),
          score_metric = "accuracy")
  }
}
agg <- do.call(rbind, by(df, df[, keys], function(d) {
  row <- d[1, keys, drop = FALSE]
  cbind(row, as.data.frame(agg_one(d), stringsAsFactors = FALSE))
}))
rownames(agg) <- NULL
agg <- agg[order(agg$response_type, agg$n_species, agg$multi_obs,
                   agg$phylo_signal, agg$beta_strength, agg$method), ]

# --- Wide format: one row per cell, methods as cols ---------------------------
keys_no_method <- setdiff(keys, "method")
wide_rmse <- reshape(agg[agg$score_metric == "rmse",
                          c(keys, "mean_score")],
                       idvar = keys_no_method,
                       timevar = "method", direction = "wide",
                       v.names = "mean_score")
# Ratio: pigauto / phylolm — defined only for gaussian
if ("mean_score.pigauto_cov_sfT" %in% names(wide_rmse) &&
    "mean_score.phylolm_lambda_blup" %in% names(wide_rmse)) {
  wide_rmse$ratio_pig_phy <- with(wide_rmse,
    mean_score.pigauto_cov_sfT / mean_score.phylolm_lambda_blup)
  wide_rmse$ratio_pig_mean <- with(wide_rmse,
    mean_score.pigauto_cov_sfT / mean_score.column_mean)
}

# --- Verdict markdown ---------------------------------------------------------
v <- c(
  sprintf("# sim_bace verdict (%s tier)", TIER),
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Source: %s", in_rds),
  "",
  "## Per-cell mean scores (RMSE for gaussian, accuracy for binary/threshold)",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE, max = 1500)),
  "```",
  "")

if (exists("wide_rmse") && "ratio_pig_phy" %in% names(wide_rmse)) {
  v <- c(v,
    "## Decisive comparison: pigauto / phylolm-lambda BLUP (gaussian only)",
    "",
    "Below 1.0 = pigauto wins; above 1.0 = phylolm wins.",
    "",
    "```",
    capture.output(print(wide_rmse[, c(setdiff(keys_no_method, "score_metric"),
                                         "mean_score.column_mean",
                                         "mean_score.lm",
                                         "mean_score.phylolm_lambda_blup",
                                         "mean_score.pigauto_cov_sfT",
                                         "ratio_pig_phy", "ratio_pig_mean")],
                          row.names = FALSE)),
    "```",
    "",
    "### Headline summary",
    sprintf("- pigauto / phylolm-lambda median ratio: **%.3f**", median(wide_rmse$ratio_pig_phy, na.rm = TRUE)),
    sprintf("- pigauto / column-mean median ratio:   **%.3f**", median(wide_rmse$ratio_pig_mean, na.rm = TRUE)),
    sprintf("- cells where pigauto beats phylolm:    %d / %d",
              sum(wide_rmse$ratio_pig_phy < 1.0, na.rm = TRUE),
              sum(!is.na(wide_rmse$ratio_pig_phy))),
    sprintf("- cells where pigauto beats column-mean: %d / %d",
              sum(wide_rmse$ratio_pig_mean < 1.0, na.rm = TRUE),
              sum(!is.na(wide_rmse$ratio_pig_mean))),
    "")
}
writeLines(v, out_md)
cat(sprintf("Wrote %s\n", out_md))

# --- Figure: pigauto/phylolm ratio heatmap (gaussian cells) -------------------
if (exists("wide_rmse") && "ratio_pig_phy" %in% names(wide_rmse)) {
  hm <- wide_rmse
  hm$ratio_cap <- pmin(pmax(hm$ratio_pig_phy, 0.5), 2.0)
  hm$beta_lab <- factor(sprintf("β=%.1f", hm$beta_strength))
  hm$ps_lab   <- factor(sprintf("phy=%.1f", hm$phylo_signal))
  hm$multi_lab <- factor(ifelse(hm$multi_obs == 1L, "single-obs", "multi-obs"))
  hm$n_lab     <- factor(sprintf("n=%d", hm$n_species))

  p <- ggplot(hm, aes(x = beta_lab, y = ps_lab, fill = ratio_cap)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", ratio_pig_phy)), size = 2.5) +
    scale_fill_gradient2(low = "#3a7d44", mid = "white", high = "#a93a3a",
                           midpoint = 1.0, limits = c(0.5, 2.0),
                           name = "pigauto / phylolm-λ\n(< 1 = pigauto wins)") +
    facet_grid(n_lab ~ multi_lab) +
    labs(x = NULL, y = NULL,
         title = sprintf("pigauto vs phylolm-λ on BACE-simulated gaussian (%s tier)", TIER),
         subtitle = "Green = pigauto wins.  Red = phylolm-λ wins.") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave(out_png, p, width = 9, height = 6, dpi = 110)
  ggsave(out_pdf, p, width = 9, height = 6)
  cat(sprintf("Wrote %s\n", out_png))
}

cat("Done.\n")
