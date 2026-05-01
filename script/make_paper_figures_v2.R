#!/usr/bin/env Rscript
# script/make_paper_figures_v2.R
#
# Adds Figure 4 (missingness sweep) and a 4-panel combined plot.

suppressPackageStartupMessages(library(ggplot2))
HERE <- "/Users/z3437171/Dropbox/Github Local/pigauto"
fig_dir <- file.path(HERE, "useful", "paper_figures")

# ----------- Fig 4: missingness -----------
df <- readRDS(file.path(HERE, "script", "bench_missingness_sweep.rds"))
agg <- aggregate(rmse ~ method + sp_miss + within_miss + f_type, data = df,
                  FUN = function(x) mean(x, na.rm=TRUE))
agg_se <- aggregate(rmse ~ method + sp_miss + within_miss + f_type, data = df,
                     FUN = function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
names(agg_se)[5] <- "se"
agg <- merge(agg, agg_se)
miss_total <- aggregate(miss_total ~ sp_miss + within_miss, data = df,
                          FUN = function(x) round(mean(x), 3))
agg <- merge(agg, miss_total, by = c("sp_miss", "within_miss"))

w <- reshape(agg[, c("method","sp_miss","within_miss","f_type","miss_total","rmse","se")],
              idvar = c("sp_miss","within_miss","f_type","miss_total"),
              timevar = "method", direction = "wide")
w$ratio <- w$rmse.pigauto_sfT / w$rmse.lm_nonlinear
w$ratio_se <- w$ratio * sqrt(
  (w$se.pigauto_sfT / w$rmse.pigauto_sfT)^2 +
  (w$se.lm_nonlinear / w$rmse.lm_nonlinear)^2
)

p4 <- ggplot(w, aes(x = miss_total, y = ratio, color = f_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.0,
               span = 0.9, formula = y ~ x) +
  geom_errorbar(aes(ymin = ratio - 1.96*ratio_se, ymax = ratio + 1.96*ratio_se),
                  width = 0.01, alpha = 0.4) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_x_continuous(breaks = seq(0.3, 0.9, 0.1)) +
  labs(x = "Total missingness fraction",
        y = "RMSE ratio: pigauto / lm_nonlinear",
        color = "DGP",
        title = "On nonlinear DGPs, pigauto's advantage GROWS with missingness",
        subtitle = expression(paste("At ", lambda, "=0.20, n=500, ", beta, "=1.0, ncov=10 (3 reps each, 18 cells total)"))) +
  annotate("text", x = 0.4, y = 1.18, label = "lm_nonlinear wins",
            hjust = 0, color = "grey20", fontface = "italic") +
  annotate("text", x = 0.85, y = 0.83, label = "pigauto wins",
            hjust = 1, color = "grey20", fontface = "italic") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig4_missingness.png"), p4, width = 7, height = 5, dpi = 150)
cat("Saved fig4_missingness.png\n")

# ----------- vs phylolm-lambda -----------
w$ratio_phylolm <- w$rmse.pigauto_sfT / w$rmse.phylolm_lambda_blup
w$ratio_phylolm_se <- w$ratio_phylolm * sqrt(
  (w$se.pigauto_sfT / w$rmse.pigauto_sfT)^2 +
  (w$se.phylolm_lambda_blup / w$rmse.phylolm_lambda_blup)^2
)

p4b <- ggplot(w, aes(x = miss_total, y = ratio_phylolm, color = f_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.0,
               span = 0.9, formula = y ~ x) +
  geom_errorbar(aes(ymin = ratio_phylolm - 1.96*ratio_phylolm_se,
                      ymax = ratio_phylolm + 1.96*ratio_phylolm_se),
                  width = 0.01, alpha = 0.4) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(nonlinear = "#1f78b4", interactive = "#e31a1c")) +
  scale_x_continuous(breaks = seq(0.3, 0.9, 0.1)) +
  labs(x = "Total missingness fraction",
        y = expression(paste("RMSE ratio: pigauto / phylolm-", lambda)),
        color = "DGP",
        title = "Pigauto consistently beats phylolm-λ across all missingness",
        subtitle = "9/9 wins on each f_type; the AE's value vs linear+phylo is robust") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig4b_missingness_vs_phylolm.png"), p4b,
        width = 7, height = 5, dpi = 150)
cat("Saved fig4b_missingness_vs_phylolm.png\n")

cat("\nDone. All paper figures in", fig_dir, "\n")
