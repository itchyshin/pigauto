#!/usr/bin/env Rscript
# script/make_paper_figures.R
#
# Generates the 3 paper-ready figures from the overnight Phase 1 sweeps:
#   Fig 1: lambda-threshold sweep
#   Fig 2: n-scaling sweep
#   Fig 3: beta-strength sweep
#
# Output: useful/paper_figures/{fig1_lambda.png, fig2_n.png, fig3_beta.png}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    cat("dplyr not installed; using base R aggregation\n")
  }
})

HERE <- "/Users/z3437171/Dropbox/Github Local/pigauto"
fig_dir <- file.path(HERE, "useful", "paper_figures")
dir.create(fig_dir, showWarnings = FALSE)

# Common helpers
make_ratio_df <- function(rds_path, x_col) {
  df <- readRDS(rds_path)
  agg <- aggregate(rmse ~ method + .data, data = df,
                    FUN = function(x) mean(x, na.rm = TRUE))
  invisible(NULL)  # not used
}

# -----------------------------------------------------------------
# FIGURE 1: lambda-threshold
# -----------------------------------------------------------------
df1 <- readRDS(file.path(HERE, "script", "bench_lambda_sweep.rds"))
agg1 <- aggregate(rmse ~ method + lambda + f_type, data = df1,
                   FUN = function(x) mean(x, na.rm = TRUE))
agg1_se <- aggregate(rmse ~ method + lambda + f_type, data = df1,
                      FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
names(agg1_se)[4] <- "se"
agg1 <- merge(agg1, agg1_se)

# Compute pigauto / lm_nonlinear ratio (with delta-method SE approximation)
w1 <- reshape(agg1[, c("method","lambda","f_type","rmse","se")],
               idvar = c("lambda","f_type"),
               timevar = "method", direction = "wide")
w1$ratio <- w1$rmse.pigauto_sfT / w1$rmse.lm_nonlinear
# Approximate ratio SE (delta-method)
w1$ratio_se <- w1$ratio * sqrt(
  (w1$se.pigauto_sfT / w1$rmse.pigauto_sfT)^2 +
  (w1$se.lm_nonlinear / w1$rmse.lm_nonlinear)^2
)

p1 <- ggplot(w1, aes(x = lambda, y = ratio, color = f_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = ratio - 1.96*ratio_se, ymax = ratio + 1.96*ratio_se,
                    fill = f_type), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_fill_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)) +
  labs(x = expression(paste("Phylogenetic signal (Pagel's ", lambda, ")")),
        y = "RMSE ratio: pigauto / lm_nonlinear",
        color = "DGP",
        fill = "DGP",
        title = "Pigauto wins above a phylogenetic-signal threshold",
        subtitle = "Lower ratio = pigauto wins more (n=500, β=1.0, ncov=10, 3 reps)") +
  annotate("text", x = 0.0, y = 1.6, label = "lm_nonlinear wins",
            hjust = 0, color = "grey20", fontface = "italic") +
  annotate("text", x = 0.4, y = 0.65, label = "pigauto wins",
            hjust = 0.5, color = "grey20", fontface = "italic") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig1_lambda_threshold.png"), p1,
        width = 7, height = 5, dpi = 150)
cat("Saved fig1_lambda_threshold.png\n")

# -----------------------------------------------------------------
# FIGURE 2: n-scaling
# -----------------------------------------------------------------
df2 <- readRDS(file.path(HERE, "script", "bench_n_scaling.rds"))
agg2 <- aggregate(rmse ~ method + n_species + f_type, data = df2,
                   FUN = function(x) mean(x, na.rm = TRUE))
agg2_se <- aggregate(rmse ~ method + n_species + f_type, data = df2,
                      FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
names(agg2_se)[4] <- "se"
agg2 <- merge(agg2, agg2_se)
w2 <- reshape(agg2[, c("method","n_species","f_type","rmse","se")],
               idvar = c("n_species","f_type"),
               timevar = "method", direction = "wide")
w2$ratio <- w2$rmse.pigauto_sfT / w2$rmse.lm_nonlinear
w2$ratio_se <- w2$ratio * sqrt(
  (w2$se.pigauto_sfT / w2$rmse.pigauto_sfT)^2 +
  (w2$se.lm_nonlinear / w2$rmse.lm_nonlinear)^2
)

p2 <- ggplot(w2, aes(x = n_species, y = ratio, color = f_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = ratio - 1.96*ratio_se, ymax = ratio + 1.96*ratio_se,
                    fill = f_type), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_fill_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_x_log10(breaks = c(200, 500, 1000, 2000)) +
  labs(x = "Number of species (log scale)",
        y = "RMSE ratio: pigauto / lm_nonlinear",
        color = "DGP",
        fill = "DGP",
        title = "Pigauto's autoencoder advantage grows with sample size",
        subtitle = expression(paste("At ", lambda, " = 0.20, β = 1.0, ncov = 10 (3 reps)"))) +
  annotate("text", x = 200, y = 1.20, label = "lm_nonlinear wins",
            hjust = 0, color = "grey20", fontface = "italic") +
  annotate("text", x = 1500, y = 0.85, label = "pigauto wins",
            hjust = 0.5, color = "grey20", fontface = "italic") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig2_n_scaling.png"), p2,
        width = 7, height = 5, dpi = 150)
cat("Saved fig2_n_scaling.png\n")

# -----------------------------------------------------------------
# FIGURE 3: beta-strength
# -----------------------------------------------------------------
df3 <- readRDS(file.path(HERE, "script", "bench_beta_sweep.rds"))
agg3 <- aggregate(rmse ~ method + beta + f_type, data = df3,
                   FUN = function(x) mean(x, na.rm = TRUE))
agg3_se <- aggregate(rmse ~ method + beta + f_type, data = df3,
                      FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
names(agg3_se)[4] <- "se"
agg3 <- merge(agg3, agg3_se)
w3 <- reshape(agg3[, c("method","beta","f_type","rmse","se")],
               idvar = c("beta","f_type"),
               timevar = "method", direction = "wide")
w3$ratio <- w3$rmse.pigauto_sfT / w3$rmse.lm_nonlinear
w3$ratio_se <- w3$ratio * sqrt(
  (w3$se.pigauto_sfT / w3$rmse.pigauto_sfT)^2 +
  (w3$se.lm_nonlinear / w3$rmse.lm_nonlinear)^2
)

p3 <- ggplot(w3, aes(x = beta, y = ratio, color = f_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = ratio - 1.96*ratio_se, ymax = ratio + 1.96*ratio_se,
                    fill = f_type), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_fill_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_x_continuous(breaks = c(0.3, 0.5, 0.7, 1.0, 1.5)) +
  labs(x = "Covariate signal strength (β)",
        y = "RMSE ratio: pigauto / lm_nonlinear",
        color = "DGP",
        fill = "DGP",
        title = "Pigauto wins at moderate, not overwhelming, signal",
        subtitle = expression(paste("At ", lambda, " = 0.20, n = 500, ncov = 10 (3 reps)"))) +
  annotate("text", x = 0.4, y = 0.6, label = "pigauto wins",
            hjust = 0.5, color = "grey20", fontface = "italic") +
  annotate("text", x = 1.4, y = 1.25, label = "lm_nonlinear wins",
            hjust = 1, color = "grey20", fontface = "italic") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig3_beta_strength.png"), p3,
        width = 7, height = 5, dpi = 150)
cat("Saved fig3_beta_strength.png\n")

# -----------------------------------------------------------------
# COMBINED: 3-panel figure
# -----------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_combined <- p1 / p2 / p3
  ggsave(file.path(fig_dir, "fig_combined_3panel.png"), p_combined,
          width = 7, height = 14, dpi = 150)
  cat("Saved fig_combined_3panel.png\n")
}

# -----------------------------------------------------------------
# Also: full RMSE comparison (not just ratio) for Fig 1
# -----------------------------------------------------------------
methods_to_plot <- c("column_mean", "lm_nonlinear", "phylolm_lambda_blup", "pigauto_sfT")
agg1_long <- agg1[agg1$method %in% methods_to_plot, ]
agg1_long$method <- factor(agg1_long$method,
  levels = c("column_mean", "lm_nonlinear", "phylolm_lambda_blup", "pigauto_sfT"),
  labels = c("column mean", "lm_nonlinear", "phylolm-λ", "pigauto"))

p1_abs <- ggplot(agg1_long, aes(x = lambda, y = rmse, color = method)) +
  geom_ribbon(aes(ymin = rmse - 1.96*se, ymax = rmse + 1.96*se, fill = method),
                alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_wrap(~ f_type, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("column mean" = "grey60",
                                  "lm_nonlinear" = "#33a02c",
                                  "phylolm-λ" = "#ff7f00",
                                  "pigauto" = "#1f78b4")) +
  scale_fill_manual(values = c("column mean" = "grey60",
                                  "lm_nonlinear" = "#33a02c",
                                  "phylolm-λ" = "#ff7f00",
                                  "pigauto" = "#1f78b4")) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)) +
  labs(x = expression(paste("Phylogenetic signal (Pagel's ", lambda, ")")),
        y = "RMSE",
        color = "Method", fill = "Method",
        title = "Per-method RMSE across phylogenetic signal",
        subtitle = "n=500, β=1.0, ncov=10, 3 reps") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig1b_per_method_rmse.png"), p1_abs,
        width = 7, height = 8, dpi = 150)
cat("Saved fig1b_per_method_rmse.png\n")

cat("\nAll figures saved in", fig_dir, "\n")
