#!/usr/bin/env Rscript
# script/make_paper_figures_v3.R
#
# Combined 2x2 panel for the paper -- pigauto's regime of advantage.

suppressPackageStartupMessages(library(ggplot2))
HERE <- "/Users/z3437171/Dropbox/Github Local/pigauto"
fig_dir <- file.path(HERE, "useful", "paper_figures")

# Helper to make a ratio plot
make_ratio_plot <- function(x_var, x_lab, x_breaks, df_path, title_text, subtitle_text) {
  df <- readRDS(df_path)
  agg <- aggregate(rmse ~ method + .data, data = df,
                    FUN = function(x) mean(x, na.rm = TRUE))
  invisible(NULL)
}

# Lambda
df1 <- readRDS(file.path(HERE, "script", "bench_lambda_sweep.rds"))
agg1 <- aggregate(rmse ~ method + lambda + f_type, data = df1,
                   FUN = function(x) mean(x, na.rm = TRUE))
ag_se <- aggregate(rmse ~ method + lambda + f_type, data = df1,
                    FUN = function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
names(ag_se)[4] <- "se"
agg1 <- merge(agg1, ag_se)
w1 <- reshape(agg1, idvar = c("lambda","f_type"), timevar = "method",
              direction = "wide")
w1$ratio <- w1$rmse.pigauto_sfT / w1$rmse.lm_nonlinear
w1$ratio_se <- w1$ratio * sqrt((w1$se.pigauto_sfT / w1$rmse.pigauto_sfT)^2 +
                                  (w1$se.lm_nonlinear / w1$rmse.lm_nonlinear)^2)

# n-scaling
df2 <- readRDS(file.path(HERE, "script", "bench_n_scaling.rds"))
agg2 <- aggregate(rmse ~ method + n_species + f_type, data = df2,
                   FUN = function(x) mean(x, na.rm = TRUE))
ag_se <- aggregate(rmse ~ method + n_species + f_type, data = df2,
                    FUN = function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
names(ag_se)[4] <- "se"
agg2 <- merge(agg2, ag_se)
w2 <- reshape(agg2, idvar = c("n_species","f_type"), timevar = "method",
              direction = "wide")
w2$ratio <- w2$rmse.pigauto_sfT / w2$rmse.lm_nonlinear
w2$ratio_se <- w2$ratio * sqrt((w2$se.pigauto_sfT / w2$rmse.pigauto_sfT)^2 +
                                  (w2$se.lm_nonlinear / w2$rmse.lm_nonlinear)^2)

# beta
df3 <- readRDS(file.path(HERE, "script", "bench_beta_sweep.rds"))
agg3 <- aggregate(rmse ~ method + beta + f_type, data = df3,
                   FUN = function(x) mean(x, na.rm = TRUE))
ag_se <- aggregate(rmse ~ method + beta + f_type, data = df3,
                    FUN = function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
names(ag_se)[4] <- "se"
agg3 <- merge(agg3, ag_se)
w3 <- reshape(agg3, idvar = c("beta","f_type"), timevar = "method",
              direction = "wide")
w3$ratio <- w3$rmse.pigauto_sfT / w3$rmse.lm_nonlinear
w3$ratio_se <- w3$ratio * sqrt((w3$se.pigauto_sfT / w3$rmse.pigauto_sfT)^2 +
                                  (w3$se.lm_nonlinear / w3$rmse.lm_nonlinear)^2)

# missingness
df4 <- readRDS(file.path(HERE, "script", "bench_missingness_sweep.rds"))
agg4 <- aggregate(rmse ~ method + sp_miss + within_miss + f_type, data = df4,
                   FUN = function(x) mean(x, na.rm = TRUE))
ag_se <- aggregate(rmse ~ method + sp_miss + within_miss + f_type, data = df4,
                    FUN = function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
names(ag_se)[5] <- "se"
agg4 <- merge(agg4, ag_se)
miss_total <- aggregate(miss_total ~ sp_miss + within_miss, data = df4, FUN = mean)
miss_total$miss_total <- round(miss_total$miss_total, 3)
agg4 <- merge(agg4, miss_total)
w4 <- reshape(agg4[, c("method","sp_miss","within_miss","f_type","miss_total","rmse","se")],
              idvar = c("sp_miss","within_miss","f_type","miss_total"),
              timevar = "method", direction = "wide")
w4$ratio <- w4$rmse.pigauto_sfT / w4$rmse.lm_nonlinear
w4$ratio_se <- w4$ratio * sqrt((w4$se.pigauto_sfT / w4$rmse.pigauto_sfT)^2 +
                                  (w4$se.lm_nonlinear / w4$rmse.lm_nonlinear)^2)

# Build a unified plot with shared aesthetics
common_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot")

clr <- c(nonlinear = "#1f78b4", interactive = "#e31a1c")

plot_panel <- function(data, x_var, x_breaks, x_lab, title) {
  ggplot(data, aes(x = .data[[x_var]], y = ratio, color = f_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = ratio - 1.96*ratio_se,
                      ymax = ratio + 1.96*ratio_se,
                      fill = f_type), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = clr) +
    scale_fill_manual(values = clr) +
    labs(x = x_lab, y = "pigauto / lm_nonlinear",
          color = "DGP", fill = "DGP", title = title) +
    common_theme
}

p1 <- plot_panel(w1, "lambda", c(0,0.05,0.1,0.15,0.2,0.3,0.5),
                  expression(paste("Phylogenetic signal (Pagel's ", lambda, ")")),
                  "(a) λ-threshold")
p2 <- plot_panel(w2, "n_species", c(200,500,1000,2000),
                  "Number of species (log)", "(b) n-scaling") +
  scale_x_log10(breaks = c(200, 500, 1000, 2000))
p3 <- plot_panel(w3, "beta", c(0.3,0.5,0.7,1.0,1.5),
                  "Covariate signal strength (β)", "(c) β-strength")
p4 <- plot_panel(w4, "miss_total", seq(0.3, 0.9, 0.1),
                  "Total missingness fraction", "(d) missingness")

# Combine
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_all <- (p1 + p2) / (p3 + p4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "fig_combined_4axes.png"),
          p_all, width = 11, height = 9, dpi = 150)
  cat("Saved fig_combined_4axes.png (paper main figure)\n")
} else {
  cat("patchwork not installed; saving panels separately\n")
  ggsave(file.path(fig_dir, "fig_panel_a_lambda.png"), p1, width = 5, height = 4, dpi = 150)
  ggsave(file.path(fig_dir, "fig_panel_b_n.png"), p2, width = 5, height = 4, dpi = 150)
  ggsave(file.path(fig_dir, "fig_panel_c_beta.png"), p3, width = 5, height = 4, dpi = 150)
  ggsave(file.path(fig_dir, "fig_panel_d_miss.png"), p4, width = 5, height = 4, dpi = 150)
}
cat("\nDone.\n")
