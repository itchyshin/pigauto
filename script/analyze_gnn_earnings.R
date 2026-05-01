#!/usr/bin/env Rscript
# script/analyze_gnn_earnings.R
#
# Harvest the GNN earnings simulation results (v1 and v2) and
# produce the decisive verdict.
#
# Pre-registered decision rule (set BEFORE looking at v2 data):
#   pigauto+cov earns its keep if:
#     mean(RMSE_pigauto_cov_sfT) < 0.90 * mean(RMSE_phylolm_lambda_blup)
#   AND
#     |mean_RMSE_pigauto - mean_RMSE_phylolm| > 1.96 * pooled_SE
#   on f in {nonlinear, interactive} and beta = 0.4 cells.
#
# Output:
#   useful/gnn_earnings_verdict.md
#   useful/gnn_earnings_per_cell.md  (raw table)
#   useful/gnn_earnings_figure.png   (heatmaps)

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_md_verdict <- file.path(here, "useful", "gnn_earnings_verdict.md")
out_md_table   <- file.path(here, "useful", "gnn_earnings_per_cell.md")
out_png        <- file.path(here, "useful", "gnn_earnings_figure.png")
out_pdf        <- file.path(here, "useful", "gnn_earnings_figure.pdf")

cat("Looking for results files ...\n")

# v1 phases (informative but biased -- BM-misspecified phylolm + iid cov)
v1_tags <- c("tree300", "mammal_n1000", "amphibio_n1500")
v1_files <- file.path(here, "script",
                       sprintf("bench_gnn_earnings_%s.rds", v1_tags))
v1_avail <- v1_files[file.exists(v1_files)]
cat(sprintf("  v1 phases: %d/%d available\n", length(v1_avail), length(v1_files)))

# v2 phases (decisive -- lambda phylolm, alpha sweep, phylo-structured cov)
v2_tags <- c("v2_tree300", "v2_mammal_n1000", "v2_amphibio_n1500")
v2_files <- file.path(here, "script",
                       sprintf("bench_gnn_earnings_%s.rds", v2_tags))
v2_avail <- v2_files[file.exists(v2_files)]
cat(sprintf("  v2 phases: %d/%d available\n", length(v2_avail), length(v2_files)))

if (length(v2_avail) == 0L) {
  cat("\nNo v2 results yet -- analysis can't produce decisive verdict.\n")
  cat("Re-run after at least one v2 phase finishes.\n")
  quit(save = "no")
}

# --------- Load v2 data ---------
load_v2 <- function(f) {
  x <- readRDS(f)
  d <- x$results
  d$tag <- x$meta$tag
  d$n_species <- x$meta$n
  d
}
v2_all <- do.call(rbind, lapply(v2_avail, load_v2))
cat(sprintf("Loaded v2: %d rows, tags = [%s]\n",
              nrow(v2_all), paste(unique(v2_all$tag), collapse=", ")))

# --------- Compute per-(alpha, beta, f, cov, method, tag) summary ---------
agg_fun <- function(d) {
  data.frame(
    mean_rmse = mean(d$rmse, na.rm = TRUE),
    se_rmse   = stats::sd(d$rmse, na.rm = TRUE) / sqrt(sum(!is.na(d$rmse))),
    mean_r    = mean(d$pearson_r, na.rm = TRUE),
    n_reps_ok = sum(!is.na(d$rmse))
  )
}
keys <- c("tag", "alpha", "beta", "f_type", "cov_type", "method")
agg <- do.call(rbind, by(v2_all, lapply(keys, function(k) v2_all[[k]]),
                            function(d) {
                              row <- d[1, keys]
                              row <- cbind(row, agg_fun(d))
                              row
                            }))
rownames(agg) <- NULL
agg <- agg[order(agg$tag, agg$cov_type, agg$alpha, agg$beta,
                   agg$f_type, agg$method), ]

# --------- Apply decision rule ---------
decisive_cells <- subset(agg,
  f_type %in% c("nonlinear", "interactive") &
  abs(beta - 0.4) < 1e-9 &
  method %in% c("pigauto_cov_sfT", "phylolm_lambda_blup"))

# Reshape: one row per (tag, alpha, f_type, cov_type), columns for the two methods
wide <- list()
for (tag in unique(decisive_cells$tag)) {
  for (a in unique(decisive_cells$alpha[decisive_cells$tag == tag])) {
    for (ft in c("nonlinear", "interactive")) {
      for (ct in c("iid", "phylo_structured")) {
        sub <- subset(decisive_cells, tag == tag & alpha == a &
                                          f_type == ft & cov_type == ct)
        if (nrow(sub) < 2L) next
        pig <- sub[sub$method == "pigauto_cov_sfT", ]
        phy <- sub[sub$method == "phylolm_lambda_blup", ]
        if (nrow(pig) == 0L || nrow(phy) == 0L) next
        ratio    <- pig$mean_rmse / phy$mean_rmse
        se_pool  <- sqrt(pig$se_rmse^2 + phy$se_rmse^2)
        gap_z    <- (phy$mean_rmse - pig$mean_rmse) / se_pool
        wins     <- (ratio < 0.90) & (gap_z > 1.96)
        wide[[length(wide) + 1L]] <- data.frame(
          tag = tag, alpha = a, f_type = ft, cov_type = ct,
          phylolm_rmse = phy$mean_rmse, phylolm_se = phy$se_rmse,
          pigauto_rmse = pig$mean_rmse, pigauto_se = pig$se_rmse,
          ratio = ratio, gap_z = gap_z, pigauto_wins = wins,
          stringsAsFactors = FALSE)
      }
    }
  }
}
verdict_df <- do.call(rbind, wide)

# --------- Write verdict markdown ---------
verdict_lines <- c(
  "# GNN earnings — decisive verdict",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Pre-registered decision rule",
  "",
  "pigauto + cov **earns its keep** vs phylolm-lambda BLUP (the smart linear baseline) if BOTH:",
  "",
  "1. mean(RMSE_pigauto_cov_sfT) < 0.90 * mean(RMSE_phylolm_lambda_blup), and",
  "2. |mean RMSE gap| > 1.96 * pooled SE.",
  "",
  "Evaluated only on cells with f in {nonlinear, interactive} and beta = 0.4 (the regime where covariates carry real signal that linear methods cannot fully exploit).",
  "",
  "## Per-cell verdict (v2 sims, smart baseline + phylo-structured cov + lambda-fitted phylolm)",
  "",
  "```",
  sprintf("%-20s %5s %-12s %-17s %10s %10s %6s %5s %s",
            "tag", "alpha", "f_type", "cov_type",
            "phylolm", "pigauto", "ratio", "z", "verdict"),
  capture.output({
    for (i in seq_len(nrow(verdict_df))) {
      r <- verdict_df[i, ]
      verdict_str <- if (is.na(r$pigauto_wins)) "?"
                     else if (r$pigauto_wins) ">>> WIN"
                     else if (r$ratio < 1) "marginal"
                     else "loses"
      cat(sprintf("%-20s %5.1f %-12s %-17s %10.3f %10.3f %6.3f %+5.2f %s\n",
                    r$tag, r$alpha, r$f_type, r$cov_type,
                    r$phylolm_rmse, r$pigauto_rmse, r$ratio, r$gap_z,
                    verdict_str))
    }
  }),
  "```",
  "",
  "## Summary by cov_type",
  "")

for (ct in c("iid", "phylo_structured")) {
  sub <- subset(verdict_df, cov_type == ct)
  if (nrow(sub) == 0L) next
  n_win <- sum(sub$pigauto_wins, na.rm = TRUE)
  median_ratio <- median(sub$ratio, na.rm = TRUE)
  verdict_lines <- c(verdict_lines,
    sprintf("- **cov_type = %s**: pigauto wins %d/%d cells (median ratio = %.3f)",
              ct, n_win, nrow(sub), median_ratio))
}

# Headline conclusion
all_wins <- sum(verdict_df$pigauto_wins, na.rm = TRUE)
n_total  <- nrow(verdict_df)
pct_win  <- 100 * all_wins / n_total

verdict_lines <- c(verdict_lines,
  "",
  "## Headline",
  "",
  sprintf("Across all decisive cells (n=%d), pigauto+cov beat the smart linear baseline by the pre-registered rule on **%d (%.1f%%)** of cells.",
            n_total, all_wins, pct_win),
  "")
if (pct_win >= 50) {
  verdict_lines <- c(verdict_lines,
    "**The GNN earns its keep.**  pigauto's neural architecture provides a meaningful improvement over phylolm-lambda BLUP on a majority of nonlinear/interactive covariate-effect regimes.")
} else if (pct_win >= 25) {
  verdict_lines <- c(verdict_lines,
    "**Mixed result.**  pigauto helps in some regimes but not others.  The paper claim should be scoped to specific (alpha, cov_type) regimes where pigauto wins.")
} else {
  verdict_lines <- c(verdict_lines,
    "**The GNN does NOT beat the smart linear baseline by a defensible margin.**  Scope the paper claim down: pigauto's value is the safety floor + multi-obs handling + uncertainty quantification, not nonlinear signal extraction beyond what phylolm-lambda already captures.")
}

writeLines(verdict_lines, out_md_verdict)
cat(sprintf("Wrote %s\n", out_md_verdict))

# --------- Per-cell raw table ---------
table_lines <- c(
  "# GNN earnings v2 — per-cell raw results",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Per (tag, alpha, beta, f_type, cov_type, method): mean RMSE, SE RMSE, mean Pearson r, n_reps_ok.",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE, max = 5000)),
  "```")
writeLines(table_lines, out_md_table)
cat(sprintf("Wrote %s\n", out_md_table))

# --------- Figure: ratio heatmap per (tag, cov_type) ---------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))

  # Compute ratio for ALL (alpha, beta, f, cov, tag) cells
  ratio_df <- list()
  for (tag in unique(agg$tag)) for (a in unique(agg$alpha)) {
    for (b in unique(agg$beta)) for (ft in unique(agg$f_type)) {
      for (ct in unique(agg$cov_type)) {
        pig <- subset(agg, tag == tag & alpha == a & beta == b &
                              f_type == ft & cov_type == ct &
                              method == "pigauto_cov_sfT")
        phy <- subset(agg, tag == tag & alpha == a & beta == b &
                              f_type == ft & cov_type == ct &
                              method == "phylolm_lambda_blup")
        if (nrow(pig) == 0L || nrow(phy) == 0L) next
        ratio_df[[length(ratio_df) + 1L]] <- data.frame(
          tag = tag, alpha = a, beta = b, f_type = ft, cov_type = ct,
          ratio = pig$mean_rmse / phy$mean_rmse,
          stringsAsFactors = FALSE)
      }
    }
  }
  rdf <- do.call(rbind, ratio_df)
  rdf$alpha_lab <- factor(sprintf("alpha = %.1f", rdf$alpha))
  rdf$beta_lab  <- factor(sprintf("beta = %.1f", rdf$beta))
  rdf$f_type    <- factor(rdf$f_type,
                            levels = c("linear", "nonlinear", "interactive"))
  rdf$cov_type  <- factor(rdf$cov_type,
                            levels = c("iid", "phylo_structured"))
  rdf$tag       <- factor(rdf$tag)

  # Cap ratio for color scale (avoid extreme outliers blowing scale)
  rdf$ratio_cap <- pmin(pmax(rdf$ratio, 0.5), 1.5)

  p <- ggplot(rdf, aes(x = beta_lab, y = f_type, fill = ratio_cap)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", ratio)),
                size = 2.5) +
    scale_fill_gradient2(low = "#3a7d44", mid = "white", high = "#a93a3a",
                           midpoint = 1.0,
                           limits = c(0.5, 1.5),
                           name = "pigauto / phylolm\n(< 1 = pigauto wins)") +
    facet_grid(tag + cov_type ~ alpha_lab,
                 labeller = label_value) +
    labs(x = NULL, y = NULL,
         title = "pigauto / phylolm-lambda RMSE ratio (v2 sim)",
         subtitle = "Green: pigauto better.  Red: phylolm better.  Threshold (pre-reg): pigauto wins iff ratio < 0.90 AND gap > 1.96 SE.") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey40"),
          strip.text = element_text(size = 8))

  n_facets_y <- length(unique(rdf$tag)) * length(unique(rdf$cov_type))
  n_facets_x <- length(unique(rdf$alpha))
  ggsave(out_png, p,
          width = 3 + 2.2 * n_facets_x,
          height = 2 + 1.5 * n_facets_y, dpi = 110)
  ggsave(out_pdf, p,
          width = 3 + 2.2 * n_facets_x,
          height = 2 + 1.5 * n_facets_y)
  cat(sprintf("Wrote %s\n", out_png))
  cat(sprintf("Wrote %s\n", out_pdf))
}

cat("Done.\n")
