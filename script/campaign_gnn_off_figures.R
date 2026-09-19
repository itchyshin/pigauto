# script/campaign_gnn_off_figures.R
# Figures for the GNN-off campaign, built on the aggregate produced by
# script/campaign_gnn_off_aggregate.R (see that script for object structure:
# summary_arm, per_trait, walls, cell_arm). Handles any subset of dgp x n
# present in the aggregate (partial or full campaign).
# Usage: Rscript script/campaign_gnn_off_figures.R <agg_rds> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
agg_rds <- if (length(args) >= 1) args[1] else "/tmp/campaign_agg.rds"
out_dir <- if (length(args) >= 2) args[2] else "/tmp/campaign_fig"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (!have_ggplot) {
  cat("ggplot2 not installed - falling back to base graphics (layout will be simpler)\n")
}

x <- readRDS(agg_rds)
summary_arm <- x$summary_arm
per_trait <- x$per_trait
walls <- x$walls
cell_arm <- x$cell_arm

mcse <- function(v) { v <- v[!is.na(v)]; if (length(v) < 2) return(NA_real_); stats::sd(v) / sqrt(length(v)) }

# ---- arm ordering / labels (shared across figures) --------------------------
arm_levels <- c("gnn_off_pure", "gnn_off", "rphylopars", "gnn_on_full", "gnn_on", "bace", "floor")
arm_labels <- c(
  gnn_off_pure = "pigauto, GNN off (pure baseline)",
  gnn_off      = "pigauto, GNN off",
  rphylopars   = "raw Rphylopars",
  gnn_on_full  = "pigauto, GNN on, tax-free baseline",
  gnn_on       = "pigauto, GNN on",
  bace         = "BACE",
  floor        = "mean/mode floor"
)
# Okabe-Ito colour-blind-safe palette, one colour per arm (fixed mapping).
arm_colors <- c(
  gnn_off_pure = "#0072B2",
  gnn_off      = "#56B4E9",
  rphylopars   = "#009E73",
  gnn_on_full  = "#E69F00",
  gnn_on       = "#D55E00",
  bace         = "#CC79A7",
  floor        = "#999999"
)

order_arms <- function(d) {
  d$arm <- factor(d$arm, levels = arm_levels[arm_levels %in% unique(d$arm)])
  d
}

n_seeds_caption <- function(d) {
  ns <- sort(unique(d$n_seeds))
  if (length(ns) == 1) paste0("seeds = ", ns) else paste0("seeds = ", paste(ns, collapse = "/"))
}

regime_title <- function(main, d) {
  bquote(.(main) ~ "(30% MCAR," ~ .(n_seeds_caption(d)) ~ ")")
}

facet_dgp_n <- if (have_ggplot) ggplot2::facet_grid(dgp ~ n, labeller = ggplot2::labeller(n = function(x) paste0("n = ", x))) else NULL

theme_fig <- if (have_ggplot) {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    )
} else NULL

save_png <- function(p, path, width = 7.5, height = 5.5) {
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 300)
}

# ==============================================================================
# 1. fig_zrmse.png - mean z-RMSE (continuous traits) by arm, dot + MCSE bar
# ==============================================================================
make_dotplot <- function(d, y, y_mcse, ylab, main, floor_note = TRUE) {
  d <- order_arms(d)
  d$lo <- d[[y]] - d[[y_mcse]]
  d$hi <- d[[y]] + d[[y_mcse]]

  non_floor <- d[[y]][d$arm != "floor" & !is.na(d[[y]])]
  y_clip <- if (length(non_floor)) stats::quantile(non_floor, 0.99, na.rm = TRUE) * 1.1 else NA_real_
  floor_present <- "floor" %in% d$arm && any(d[[y]][d$arm == "floor"] > y_clip, na.rm = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = arm, y = .data[[y]], colour = arm)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.2, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    facet_dgp_n +
    ggplot2::scale_colour_manual(values = arm_colors, breaks = arm_levels, labels = arm_labels) +
    ggplot2::labs(x = NULL, y = ylab, title = regime_title(main, d)) +
    theme_fig +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))

  if (floor_present && !is.na(y_clip)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, y_clip)) +
      ggplot2::labs(caption = "floor arm is off-scale (clipped at ~99th percentile of the other arms); see fig caption / summary table for its value")
  }
  p
}

d_zrmse <- summary_arm[!is.na(summary_arm$zrmse), ]
p_zrmse <- make_dotplot(d_zrmse, "zrmse", "zrmse_mcse", "mean z-RMSE (continuous traits)", "Continuous-trait z-RMSE by arm")
if (have_ggplot) save_png(p_zrmse, file.path(out_dir, "fig_zrmse.png"))

# ==============================================================================
# 2. fig_accuracy.png - discrete accuracy by arm, same layout
# ==============================================================================
d_acc <- summary_arm[!is.na(summary_arm$acc) & !is.nan(summary_arm$acc), ]
p_acc <- make_dotplot(d_acc, "acc", "acc_mcse", "mean accuracy (discrete traits)", "Discrete-trait accuracy by arm")
if (have_ggplot) save_png(p_acc, file.path(out_dir, "fig_accuracy.png"))

# ==============================================================================
# 3. fig_coverage.png - production-interval coverage for pigauto arms, 0.95 line
# ==============================================================================
# summary_arm has no coverage MCSE column; recompute mean/MCSE across seeds
# from cell_arm (already averaged across continuous traits per seed/arm).
pigauto_arms <- c("gnn_off_pure", "gnn_off", "gnn_on_full", "gnn_on")
ca_cov <- cell_arm[cell_arm$arm %in% pigauto_arms & !is.na(cell_arm$coverage), ]
cov_summary <- do.call(rbind, lapply(
  split(ca_cov, list(ca_cov$dgp, ca_cov$n, ca_cov$arm), drop = TRUE),
  function(d) data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_seeds = nrow(d),
                          coverage = mean(d$coverage), coverage_mcse = mcse(d$coverage))
))

if (!is.null(cov_summary) && nrow(cov_summary)) {
  cov_summary <- order_arms(cov_summary)
  cov_summary$lo <- cov_summary$coverage - cov_summary$coverage_mcse
  cov_summary$hi <- cov_summary$coverage + cov_summary$coverage_mcse

  p_cov <- ggplot2::ggplot(cov_summary, ggplot2::aes(x = arm, y = coverage, colour = arm)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.2, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    facet_dgp_n +
    ggplot2::scale_colour_manual(values = arm_colors, breaks = arm_levels, labels = arm_labels) +
    ggplot2::labs(x = NULL, y = "production-interval coverage",
                  title = regime_title("Prediction-interval coverage (pigauto arms)", cov_summary)) +
    theme_fig +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
  if (have_ggplot) save_png(p_cov, file.path(out_dir, "fig_coverage.png"))
} else {
  cat("no coverage data for pigauto arms - fig_coverage.png skipped\n")
}

# ==============================================================================
# 4. fig_wall.png - wall time per arm, log10 seconds, facet by n, colour by dgp
#    gnn_on_full has wall = 0 (derived) -> drop it.
# ==============================================================================
w <- walls[walls$arm != "gnn_on_full", ]
w <- w[w$wall > 0, ]
w <- order_arms(w)
w_summary <- do.call(rbind, lapply(
  split(w, list(w$dgp, w$n, w$arm), drop = TRUE),
  function(d) data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_seeds = nrow(d),
                          log10_wall = mean(log10(d$wall)), log10_wall_mcse = mcse(log10(d$wall)))
))
w_summary <- order_arms(w_summary)
w_summary$lo <- w_summary$log10_wall - w_summary$log10_wall_mcse
w_summary$hi <- w_summary$log10_wall + w_summary$log10_wall_mcse

p_wall <- ggplot2::ggplot(w_summary, ggplot2::aes(x = arm, y = log10_wall, colour = dgp, group = dgp)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.2,
                          position = ggplot2::position_dodge(width = 0.4), na.rm = TRUE) +
  ggplot2::geom_point(size = 2.2, position = ggplot2::position_dodge(width = 0.4), na.rm = TRUE) +
  ggplot2::facet_wrap(~ n, labeller = ggplot2::labeller(n = function(x) paste0("n = ", x))) +
  ggplot2::scale_colour_manual(values = c(bm_mixed = "#0072B2", ou_mixed = "#D55E00",
                                           bace_dgp = "#009E73", avonet = "#CC79A7")) +
  ggplot2::scale_x_discrete(labels = arm_labels) +
  ggplot2::labs(x = NULL, y = expression(log[10] ~ "wall time (s)"), colour = "dgp",
                title = regime_title("Wall time by arm (gnn_on_full omitted - wall is derived, = 0)", w_summary)) +
  theme_fig +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
if (have_ggplot) save_png(p_wall, file.path(out_dir, "fig_wall.png"), width = 8, height = 5.5)

# ==============================================================================
# 5. fig_tax.png - tax decomposition: gnn_off -> gnn_on_full -> gnn_on per seed
# ==============================================================================
tax_arms <- c("gnn_off", "gnn_on_full", "gnn_on")
d_tax <- cell_arm[cell_arm$arm %in% tax_arms, ]
d_tax$arm <- factor(d_tax$arm, levels = tax_arms)

tax_mean <- do.call(rbind, lapply(
  split(d_tax, list(d_tax$dgp, d_tax$n, d_tax$arm), drop = TRUE),
  function(d) data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_seeds = nrow(d),
                          zrmse = mean(d$zrmse, na.rm = TRUE), zrmse_mcse = mcse(d$zrmse))
))
tax_mean$arm <- factor(tax_mean$arm, levels = tax_arms)

p_tax <- ggplot2::ggplot(d_tax, ggplot2::aes(x = arm, y = zrmse, group = seed)) +
  ggplot2::geom_line(colour = "grey70", alpha = 0.5, na.rm = TRUE) +
  ggplot2::geom_point(colour = "grey70", alpha = 0.5, size = 1, na.rm = TRUE) +
  ggplot2::geom_line(data = tax_mean, ggplot2::aes(group = 1), colour = "#D55E00", linewidth = 1, na.rm = TRUE) +
  ggplot2::geom_errorbar(data = tax_mean,
                          ggplot2::aes(ymin = zrmse - zrmse_mcse, ymax = zrmse + zrmse_mcse, group = 1),
                          colour = "#D55E00", width = 0.15, na.rm = TRUE) +
  ggplot2::geom_point(data = tax_mean, ggplot2::aes(group = 1), colour = "#D55E00", size = 2.5, na.rm = TRUE) +
  facet_dgp_n +
  ggplot2::scale_x_discrete(labels = arm_labels[tax_arms]) +
  ggplot2::labs(x = NULL, y = "mean z-RMSE (continuous traits)",
                title = regime_title("Tax decomposition: baseline tax paid by GNN-on", tax_mean),
                caption = "grey lines: per-seed mean; orange: seed mean +/- MCSE") +
  theme_fig +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
if (have_ggplot) save_png(p_tax, file.path(out_dir, "fig_tax.png"))

# ==============================================================================
# Confirm outputs
# ==============================================================================
pngs <- c("fig_zrmse.png", "fig_accuracy.png", "fig_coverage.png", "fig_wall.png", "fig_tax.png")
for (f in pngs) {
  fp <- file.path(out_dir, f)
  if (file.exists(fp)) {
    cat(sprintf("%s: %d bytes\n", fp, file.info(fp)$size))
  } else {
    cat(sprintf("%s: MISSING\n", fp))
  }
}
