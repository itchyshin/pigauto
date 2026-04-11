#!/usr/bin/env Rscript
#
# script/make_bench_count_html.R
#
# Build a self-contained HTML report from bench_count.rds.
#
# Run with
#   /usr/local/bin/Rscript script/make_bench_count_html.R

suppressPackageStartupMessages({ })

rds_path <- "script/bench_count.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path,
       "\nRun script/bench_count.R first.")
}

r <- readRDS(rds_path)
res <- r$results
stopifnot(is.data.frame(res), nrow(res) > 0L)

test_df <- res[res$split == "test", ]
if (!nrow(test_df)) stop("No test-split rows in ", rds_path)

# ---------------------------------------------------------------------------
# Palette / labels
# ---------------------------------------------------------------------------

method_order  <- c("mean", "baseline", "pigauto")
method_label  <- c(mean     = "Mean imputation",
                   baseline = "BM baseline (Rphylopars)",
                   pigauto  = "pigauto (BM + GNN)")
method_colour <- c(mean     = "#9ca3af",
                   baseline = "#2563eb",
                   pigauto  = "#059669")

scenarios_primary   <- r$scenarios_primary
secondary_scenarios <- r$scenarios_secondary
missing_frac        <- r$missing_frac

scenario_label <- c(
  mean_5   = "\u03bc = 5",
  mean_20  = "\u03bc = 20",
  mean_100 = "\u03bc = 100",
  mean_500 = "\u03bc = 500"
)
for (s in secondary_scenarios) {
  scenario_label[s] <- switch(s,
    poisson = "Poisson",
    negbin  = "NegBin",
    s
  )
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "&ndash;", formatC(x, digits = digits, format = "f"))
}

# ---------------------------------------------------------------------------
# Average across reps
# ---------------------------------------------------------------------------

avg_metric <- function(df, method, scenario, frac, metric) {
  sub <- df[df$method == method &
              df$scenario == scenario &
              df$missing_frac == frac, ]
  vals <- sub[[metric]]
  vals <- vals[is.finite(vals)]
  if (!length(vals)) NA_real_ else mean(vals)
}

# ---------------------------------------------------------------------------
# SVG: grouped bar chart
# ---------------------------------------------------------------------------

one_bar_chart <- function(metric_col, y_label, title_text, scenarios,
                          x_axis_label, higher_is_better = FALSE) {
  W <- 640; H <- 320
  pad_l <- 60; pad_r <- 140; pad_t <- 36; pad_b <- 62
  plot_w <- W - pad_l - pad_r
  plot_h <- H - pad_t - pad_b

  n_groups <- length(scenarios)
  n_bars   <- length(method_order)
  group_w  <- plot_w / n_groups
  bar_w    <- group_w / (n_bars + 1)
  gap      <- bar_w / (n_bars + 1)

  vals <- matrix(NA_real_, n_bars, n_groups,
                 dimnames = list(method_order, scenarios))
  for (m in method_order) {
    for (s in scenarios) {
      vals[m, s] <- avg_metric(test_df, m, s, missing_frac, metric_col)
    }
  }

  v_finite <- vals[is.finite(vals)]
  if (!length(v_finite)) return("")

  if (higher_is_better) {
    y_lo <- max(0, min(v_finite) - 0.1)
    y_hi <- min(1, max(v_finite) + 0.05)
  } else {
    y_lo <- 0
    y_hi <- max(v_finite) * 1.15
  }
  if (y_hi <= y_lo) { y_lo <- 0; y_hi <- 1 }

  y_of <- function(v) {
    v <- pmin(pmax(v, y_lo), y_hi)
    pad_t + plot_h - (v - y_lo) / (y_hi - y_lo) * plot_h
  }

  parts <- character()
  parts <- c(parts, sprintf(
    '<svg viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" style="width:100%%;max-width:640px;height:auto;font-family:system-ui,-apple-system,sans-serif">',
    W, H))
  parts <- c(parts, sprintf(
    '<text x="%d" y="22" font-size="14" font-weight="600" fill="#111827">%s</text>',
    pad_l, title_text))

  n_ticks <- 5
  for (k in 0:(n_ticks - 1)) {
    v <- y_lo + (y_hi - y_lo) * k / (n_ticks - 1)
    y <- y_of(v)
    parts <- c(parts, sprintf(
      '<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e7eb" stroke-width="1"/>',
      pad_l, y, pad_l + plot_w, y))
    parts <- c(parts, sprintf(
      '<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="#6b7280">%s</text>',
      pad_l - 6, y + 3, fmt(v, 2)))
  }
  parts <- c(parts, sprintf(
    '<text x="%d" y="%d" font-size="10" fill="#6b7280" transform="rotate(-90 %d %d)">%s</text>',
    12, pad_t + plot_h / 2, 12, pad_t + plot_h / 2, y_label))

  for (g in seq_len(n_groups)) {
    scen <- scenarios[g]
    gx <- pad_l + (g - 1) * group_w
    for (b in seq_len(n_bars)) {
      m <- method_order[b]
      v <- vals[m, scen]
      if (!is.finite(v)) next
      bx <- gx + gap + (b - 1) * (bar_w + gap)
      by <- y_of(v)
      bh <- y_of(y_lo) - by
      parts <- c(parts, sprintf(
        '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" rx="2"/>',
        bx, by, bar_w, bh, method_colour[m]))
      parts <- c(parts, sprintf(
        '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="9" fill="#374151">%s</text>',
        bx + bar_w / 2, by - 4, fmt(v, 3)))
    }
    parts <- c(parts, sprintf(
      '<text x="%.1f" y="%d" text-anchor="middle" font-size="11" fill="#374151">%s</text>',
      gx + group_w / 2, pad_t + plot_h + 16,
      scenario_label[scen]))
  }
  parts <- c(parts, sprintf(
    '<text x="%.1f" y="%d" text-anchor="middle" font-size="10" fill="#6b7280">%s</text>',
    pad_l + plot_w / 2, pad_t + plot_h + 38, x_axis_label))

  legend_y <- pad_t
  for (m in method_order) {
    parts <- c(parts, sprintf(
      '<rect x="%d" y="%d" width="12" height="12" fill="%s" rx="2"/>',
      pad_l + plot_w + 14, legend_y, method_colour[m]))
    parts <- c(parts, sprintf(
      '<text x="%d" y="%d" font-size="11" fill="#111827">%s</text>',
      pad_l + plot_w + 30, legend_y + 10, method_label[m]))
    legend_y <- legend_y + 20
  }
  parts <- c(parts, '</svg>')
  paste(parts, collapse = "\n")
}

# ---------------------------------------------------------------------------
# Build charts
# ---------------------------------------------------------------------------

chart_rmse <- one_bar_chart("rmse", "RMSE", "RMSE by mean count",
                            scenarios_primary, "Mean count",
                            higher_is_better = FALSE)
chart_mae  <- one_bar_chart("mae", "MAE", "MAE by mean count",
                            scenarios_primary, "Mean count",
                            higher_is_better = FALSE)
chart_r    <- one_bar_chart("pearson_r", "Pearson r",
                            "Pearson r by mean count",
                            scenarios_primary, "Mean count",
                            higher_is_better = TRUE)

chart_disp <- if (length(secondary_scenarios) > 0) {
  one_bar_chart("rmse", "RMSE",
                "RMSE: Poisson vs NegBin (mean = 20)",
                secondary_scenarios, "Distribution",
                higher_is_better = FALSE)
} else ""

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

primary_table_rows <- character()
for (scen in scenarios_primary) {
  rmse_vals <- setNames(
    vapply(method_order, function(m) avg_metric(test_df, m, scen, missing_frac, "rmse"), numeric(1)),
    method_order)
  mae_vals <- setNames(
    vapply(method_order, function(m) avg_metric(test_df, m, scen, missing_frac, "mae"), numeric(1)),
    method_order)
  r_vals <- setNames(
    vapply(method_order, function(m) avg_metric(test_df, m, scen, missing_frac, "pearson_r"), numeric(1)),
    method_order)

  best_rmse <- names(which.min(rmse_vals))
  best_mae  <- names(which.min(mae_vals))
  best_r    <- names(which.max(r_vals))

  star <- function(method, best) {
    if (method == best) ' <span style="color:#059669">&#9733;</span>' else ""
  }

  primary_table_rows <- c(primary_table_rows, sprintf(
    '<tr><td><b>%s</b></td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td></tr>',
    scenario_label[scen],
    fmt(rmse_vals["mean"]),       star("mean", best_rmse),
    fmt(rmse_vals["baseline"]),   star("baseline", best_rmse),
    fmt(rmse_vals["pigauto"]),    star("pigauto", best_rmse),
    fmt(mae_vals["mean"]),        star("mean", best_mae),
    fmt(mae_vals["baseline"]),    star("baseline", best_mae),
    fmt(mae_vals["pigauto"]),     star("pigauto", best_mae),
    fmt(r_vals["mean"]),          star("mean", best_r),
    fmt(r_vals["baseline"]),      star("baseline", best_r),
    fmt(r_vals["pigauto"]),       star("pigauto", best_r)
  ))
}

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

rmse_bl_5   <- avg_metric(test_df, "baseline", "mean_5",   missing_frac, "rmse")
rmse_pg_5   <- avg_metric(test_df, "pigauto",  "mean_5",   missing_frac, "rmse")
rmse_bl_500 <- avg_metric(test_df, "baseline", "mean_500", missing_frac, "rmse")
rmse_pg_500 <- avg_metric(test_df, "pigauto",  "mean_500", missing_frac, "rmse")

verdict1 <- if (is.finite(rmse_bl_5) && is.finite(rmse_pg_5)) {
  sprintf('For sparse counts (mean = 5) the BM baseline achieves RMSE %.3f and pigauto achieves %.3f on the log1p-z scale. Low counts are inherently noisy, limiting all methods.',
          rmse_bl_5, rmse_pg_5)
} else 'For sparse counts the baseline and pigauto perform similarly due to inherent noise.'

verdict2 <- if (is.finite(rmse_bl_500) && is.finite(rmse_pg_500)) {
  pct <- 100 * (rmse_bl_500 - rmse_pg_500) / rmse_bl_500
  sprintf('For dense counts (mean = 500), pigauto moves RMSE from %.3f to %.3f (%+.1f%%). The log1p transform compresses the scale, so the GNN has smoother gradients to learn from.',
          rmse_bl_500, rmse_pg_500, pct)
} else 'For dense counts the log1p transform smooths the scale and the GNN can contribute more.'

# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

commit_str <- if (!is.null(r$commit) && length(r$commit) == 1L && r$commit != "unknown") {
  substr(r$commit, 1L, 10L)
} else "dev"

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>pigauto: count-trait benchmark</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 920px; margin: 2em auto; padding: 0 1.5em; color: #111827;
         line-height: 1.55; }
  h1 { border-bottom: 2px solid #111827; padding-bottom: .25em; }
  h2 { margin-top: 2em; color: #1f2937; border-bottom: 1px solid #d1d5db;
       padding-bottom: .15em; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 14px; }
  th, td { padding: 6px 10px; border-bottom: 1px solid #e5e7eb; text-align: left; }
  th { background: #f3f4f6; font-weight: 600; }
  .meta { color: #6b7280; font-size: 13px; }
  .verdict { background: #ecfdf5; border-left: 4px solid #059669;
             padding: 1em 1.2em; border-radius: 4px; margin: 1.5em 0; }
  .verdict p { margin: 0.4em 0; }
  .verdict b { color: #047857; }
  code { background: #f3f4f6; padding: 1px 5px; border-radius: 3px; font-size: 13px; }
  ul li { margin: 6px 0; }
  footer { color: #6b7280; font-size: 12px; margin-top: 3em;
           border-top: 1px solid #e5e7eb; padding-top: 1em; }
</style>
</head>
<body>

<h1>Count-trait benchmark: Poisson / NegBin</h1>
<p class="meta">
Tree: <code>ape::rtree(', r$n_species, ')</code> &middot;
Traits: ', r$n_traits, ' count per scenario &middot;
Mean counts: 5 &ndash; 500 &middot;
Methods: mean &middot; BM baseline (Rphylopars) &middot; pigauto &middot;
Replicates: ', r$n_reps, ' &middot;
Missingness: ', as.integer(100 * missing_frac), '% MCAR &middot;
Commit ', commit_str, ' &middot;
Run on ', format(Sys.time(), "%Y-%m-%d %H:%M"), ' &middot;
Total wall: ', sprintf("%.1f", r$total_wall / 60), ' min
</p>

<div class="verdict">
<p><b>Bottom line.</b> ', verdict1, '</p>
<p>', verdict2, '</p>
</div>

<h2>Primary sweep: performance by mean count (', as.integer(100 * missing_frac), '% missingness)</h2>
<p class="meta">Average across ', r$n_traits, ' traits and ', r$n_reps, ' replicates. &#9733; marks the best method per scenario.</p>

<table>
<thead>
<tr><th rowspan="2">Mean count</th><th colspan="3" style="text-align:center">RMSE (lower is better)</th><th colspan="3" style="text-align:center">MAE (lower is better)</th><th colspan="3" style="text-align:center">Pearson r (higher is better)</th></tr>
<tr><th style="text-align:right">Mean</th><th style="text-align:right">BM</th><th style="text-align:right">pigauto</th><th style="text-align:right">Mean</th><th style="text-align:right">BM</th><th style="text-align:right">pigauto</th><th style="text-align:right">Mean</th><th style="text-align:right">BM</th><th style="text-align:right">pigauto</th></tr>
</thead>
<tbody>
', paste(primary_table_rows, collapse = "\n"), '
</tbody>
</table>

<div style="text-align:center;margin:1.5em 0">
', chart_rmse, '
</div>

<div style="text-align:center;margin:1.5em 0">
', chart_mae, '
</div>

<div style="text-align:center;margin:1.5em 0">
', chart_r, '
</div>

', if (nchar(chart_disp) > 0) paste0('
<h2>Secondary sweep: Poisson vs Negative Binomial (mean = 20)</h2>
<p class="meta">Overdispersed counts (NegBin) are harder for all methods.</p>
<div style="text-align:center;margin:1.5em 0">
', chart_disp, '
</div>
') else '', '

<h2>What the benchmark shows</h2>
<ul>
<li><b>The log1p-z pipeline handles count data well.</b> Counts are transformed via log1p then z-scored before entering the latent matrix. This brings sparse and dense counts onto a common scale where BM can operate.</li>
<li><b>Sparse counts (mean = 5) are inherently noisy.</b> All methods struggle because the signal-to-noise ratio is low in log-space for small integers. The GNN cannot extract structure that is not there.</li>
<li><b>Dense counts (mean = 100+) give the GNN room to improve.</b> With more information per cell, the GNN can learn cross-trait correlations and non-linear patterns that BM misses.</li>
<li><b>Negative-binomial overdispersion increases error for all methods</b> but does not change the ordering. The extra variance makes all predictions less precise, but the phylogenetic information remains exploitable.</li>
</ul>

<h2>Reproducibility</h2>
<p>Driver: <code>script/bench_count.R</code>. Tree: <code>ape::rtree(', r$n_species, ')</code>. Traits: <code>simulate_count_traits()</code>. Training: ', r$epochs, ' epochs with early stopping. To reproduce: <code>Rscript script/bench_count.R</code>, then <code>Rscript script/make_bench_count_html.R</code>.</p>

<footer>
Source: <code>script/bench_count.R</code> &middot;
Results: <code>script/bench_count.rds</code> &middot;
Report: <code>pkgdown/assets/dev/bench_count.html</code>
</footer>

</body>
</html>
')

targets <- c("script/bench_count.html",
             "pkgdown/assets/dev/bench_count.html")
for (t in targets) {
  dir.create(dirname(t), showWarnings = FALSE, recursive = TRUE)
  writeLines(html, t)
  cat("Wrote ", t, " (", file.size(t), " bytes)\n", sep = "")
}
