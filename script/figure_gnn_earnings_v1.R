#!/usr/bin/env Rscript
# Quick figure for v1 GNN earnings results.
# pigauto_cov_sfT vs phylolm_blup ratio across (beta, f_type, tree).

suppressPackageStartupMessages({ library(ggplot2) })
here <- "/Users/z3437171/Dropbox/Github Local/pigauto"

rows <- list()
for (tg in c("tree300", "mammal_n1000", "amphibio_n1500")) {
  rds <- file.path(here, "script", sprintf("bench_gnn_earnings_%s.rds", tg))
  if (!file.exists(rds)) next
  x <- readRDS(rds)
  d <- x$results
  agg <- aggregate(rmse ~ method + beta + f_type, data = d, FUN = mean)
  w <- reshape(agg, idvar = c("beta", "f_type"), timevar = "method", direction = "wide")
  w$tree <- sprintf("%s (n=%d)", tg, x$meta$n)
  w$ratio <- w$rmse.pigauto_cov_sfT / w$rmse.phylolm_blup
  rows[[length(rows)+1L]] <- w[w$beta > 0,
                                  c("tree", "f_type", "beta", "ratio")]
}
rdf <- do.call(rbind, rows)
rdf$tree   <- factor(rdf$tree,
                       levels = c("tree300 (n=300)",
                                    "mammal_n1000 (n=1000)",
                                    "amphibio_n1500 (n=1500)"))
rdf$f_type <- factor(rdf$f_type,
                       levels = c("linear", "nonlinear", "interactive"))
rdf$beta_lab <- factor(sprintf("beta=%.1f", rdf$beta))
rdf$ratio_cap <- pmin(pmax(rdf$ratio, 0.5), 2.0)

p <- ggplot(rdf, aes(x = beta_lab, y = f_type, fill = ratio_cap)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", ratio)), size = 3) +
  scale_fill_gradient2(low = "#3a7d44", mid = "white", high = "#a93a3a",
                         midpoint = 1.0, limits = c(0.5, 2.0),
                         name = "pigauto_cov_sfT / phylolm_BM_BLUP\n(<1 pigauto wins)",
                         oob = scales::squish) +
  facet_wrap(~ tree, ncol = 3) +
  labs(x = NULL, y = NULL,
       title = "v1 GNN earnings: pigauto+cov vs phylolm BLUP (BM-misspecified)",
       subtitle = paste("Each cell: 3 reps, alpha=0.4 fixed, 30% MCAR.",
                          "phylolm here uses model='BM' (handicap;",
                          "v2 with model='lambda' is decisive).",
                          sep = "  ")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))

ggsave("useful/gnn_earnings_v1_figure.png", p, width = 11, height = 4, dpi = 110)
ggsave("useful/gnn_earnings_v1_figure.pdf", p, width = 11, height = 4)
cat("wrote useful/gnn_earnings_v1_figure.{png,pdf}\n")
