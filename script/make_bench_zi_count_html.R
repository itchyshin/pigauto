#!/usr/bin/env Rscript
#
# script/make_bench_zi_count_html.R
#
# Build a self-contained HTML report from bench_zi_count.rds.
#
# Run with
#   /usr/local/bin/Rscript script/make_bench_zi_count_html.R

suppressPackageStartupMessages({ })

rds_path <- "script/bench_zi_count.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path,
       "\nRun script/bench_zi_count.R first.")
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
                   baseline = "LP + BM baseline",
                   pigauto  = "pigauto (LP + BM + GNN)")
method_colour <- c(mean     = "#9ca3af",
                   baseline = "#2563eb",
                   pigauto  = "#059669")

scenarios_primary   <- r$scenarios_primary
secondary_scenarios <- r$secondary_scenarios
primary_frac        <- r$primary_frac

scenario_label <- c(
  zf_0.2 = "20% zeros",
  zf_0.4 = "40% zeros",
  zf_0.6 = "60% zeros",
  zf_0.8 = "80% zeros"
)
for (s in secondary_scenarios) {
  scenario_label[s] <- switch(s,
    mean_nz_5   = "\u03bc(nz) = 5",
    mean_nz_20  = "\u03bc(nz) = 20",
    mean_nz_100 = "\u03bc(nz) = 100",
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
  pad_l <- 60; pad_r <- 160; pad_t <- 36; pad_b <- 62
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
      vals[m, s] <- avg_metric(test_df, m, s, primary_frac, metric_col)
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

chart_rmse <- one_bar_chart("rmse", "RMSE",
                            "RMSE by zero fraction",
                            scenarios_primary, "Zero fraction",
                            higher_is_better = FALSE)
chart_zacc <- one_bar_chart("zero_accuracy", "Zero accuracy",
                            "Zero-class accuracy by zero fraction",
                            scenarios_primary, "Zero fraction",
                            higher_is_better = TRUE)

chart_mean_nz <- if (length(secondary_scenarios) > 0) {
  one_bar_chart("rmse", "RMSE",
                "RMSE by non-zero mean (zero frac = 0.5)",
                secondary_scenarios, "Mean of non-zero counts",
                higher_is_better = FALSE)
} else ""

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

primary_table_rows <- character()
for (scen in scenarios_primary) {
  rmse_vals <- setNames(
    vapply(method_order, function(m) avg_metric(test_df, m, scen, primary_frac, "rmse"), numeric(1)),
    method_order)
  zacc_vals <- setNames(
    vapply(method_order, function(m) avg_metric(test_df, m, scen, primary_frac, "zero_accuracy"), numeric(1)),
    method_order)

  best_rmse <- names(which.min(rmse_vals))
  best_zacc <- names(which.max(zacc_vals))

  star <- function(method, best) {
    if (method == best) ' <span style="color:#059669">&#9733;</span>' else ""
  }

  primary_table_rows <- c(primary_table_rows, sprintf(
    '<tr><td><b>%s</b></td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td><td style="text-align:right">%s%s</td></tr>',
    scenario_label[scen],
    fmt(rmse_vals["mean"]),      star("mean", best_rmse),
    fmt(rmse_vals["baseline"]),  star("baseline", best_rmse),
    fmt(rmse_vals["pigauto"]),   star("pigauto", best_rmse),
    fmt(zacc_vals["mean"]),      star("mean", best_zacc),
    fmt(zacc_vals["baseline"]),  star("baseline", best_zacc),
    fmt(zacc_vals["pigauto"]),   star("pigauto", best_zacc)
  ))
}

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

rmse_bl_20 <- avg_metric(test_df, "baseline", "zf_0.2", primary_frac, "rmse")
rmse_pg_20 <- avg_metric(test_df, "pigauto",  "zf_0.2", primary_frac, "rmse")
rmse_bl_80 <- avg_metric(test_df, "baseline", "zf_0.8", primary_frac, "rmse")
rmse_pg_80 <- avg_metric(test_df, "pigauto",  "zf_0.8", primary_frac, "rmse")
zacc_bl_80 <- avg_metric(test_df, "baseline", "zf_0.8", primary_frac, "zero_accuracy")
zacc_pg_80 <- avg_metric(test_df, "pigauto",  "zf_0.8", primary_frac, "zero_accuracy")

verdict1 <- if (is.finite(rmse_bl_20) && is.finite(rmse_pg_20)) {
  sprintf('At 20%% zero inflation, the baseline achieves RMSE %.3f and pigauto achieves %.3f. With few structural zeros the count component dominates and BM on the log1p-z scale performs well.',
          rmse_bl_20, rmse_pg_20)
} else 'At low zero inflation the count component dominates and baseline performs well.'

verdict2 <- if (is.finite(rmse_bl_80) && is.finite(rmse_pg_80) &&
                is.finite(zacc_bl_80) && is.finite(zacc_pg_80)) {
  sprintf('At 80%% zero inflation, RMSE rises to %.3f (baseline) vs %.3f (pigauto), and zero-class accuracy is %.1f%% vs %.1f%%. Heavy zero inflation creates a two-component mixture that challenges all methods, but the GNN can learn the zero/non-zero boundary from cross-trait patterns.',
          rmse_bl_80, rmse_pg_80, 100 * zacc_bl_80, 100 * zacc_pg_80)
} else 'At high zero inflation the two-component mixture challenges all methods.'

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
<title>pigauto: zero-inflated count benchmark</title>
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

<h1>Zero-inflated count benchmark: zero fraction sweep</h1>
<p class="meta">
Tree: <code>ape::rtree(', r$n_species, ')</code> &middot;
Traits: ', r$n_traits, ' zero-inflated count per scenario &middot;
Zero fraction: 0.2 &ndash; 0.8 &middot;
Methods: mean &middot; LP + BM baseline &middot; pigauto &middot;
Replicates: ', r$n_reps, ' &middot;
Missingness: ', as.integer(100 * primary_frac), '% MCAR &middot;
Commit ', commit_str, ' &middot;
Run on ', format(Sys.time(), "%Y-%m-%d %H:%M"), ' &middot;
Total wall: ', sprintf("%.1f", r$total_wall / 60), ' min
</p>

<div class="verdict">
<p><b>Bottom line.</b> ', verdict1, '</p>
<p>', verdict2, '</p>
</div>

<h2>Primary sweep: performance by zero fraction (', as.integer(100 * primary_frac), '% missingness)</h2>
<p class="meta">Average across ', r$n_traits, ' traits and ', r$n_reps, ' replicates. &#9733; marks the best method per scenario.</p>

<table>
<thead>
<tr><th rowspan="2">Zero fraction</th><th colspan="3" style="text-align:center">RMSE (lower is better)</th><th colspan="3" style="text-align:center">Zero accuracy (higher is better)</th></tr>
<tr><th style="text-align:right">Mean</th><th style="text-align:right">Baseline</th><th style="text-align:right">pigauto</th><th style="text-align:right">Mean</th><th style="text-align:right">Baseline</th><th style="text-align:right">pigauto</th></tr>
</thead>
<tbody>
', paste(primary_table_rows, collapse = "\n"), '
</tbody>
</table>

<div style="text-align:center;margin:1.5em 0">
', chart_rmse, '
</div>

<div style="text-align:center;margin:1.5em 0">
', chart_zacc, '
</div>

', if (nchar(chart_mean_nz) > 0) paste0('
<h2>Secondary sweep: non-zero mean (zero frac = 0.5)</h2>
<p class="meta">How the mean of non-zero counts affects imputation quality at fixed 50% zero inflation.</p>
<div style="text-align:center;margin:1.5em 0">
', chart_mean_nz, '
</div>
') else '', '

<h2>What the benchmark shows</h2>
<ul>
<li><b>Zero-inflated counts are a two-component mixture.</b> The data has structural zeros (drawn from a zero-generating process) and positive counts (drawn from a count process). Imputation must handle both components.</li>
<li><b>Low zero inflation behaves like regular counts.</b> At 20% zeros the count component dominates, and the log1p-z BM baseline captures most of the variation. The GNN gate stays near zero.</li>
<li><b>High zero inflation is challenging for all methods.</b> At 80% zeros, distinguishing structural zeros from missing data becomes difficult. The GNN can learn the zero/non-zero boundary from cross-trait correlations, but the fundamental identifiability challenge remains.</li>
<li><b>Zero-class accuracy is a key secondary metric.</b> Overall RMSE may be misleading when most true values are zero. Tracking whether each cell correctly identifies as zero or non-zero gives a clearer picture of structural accuracy.</li>
</ul>

<h2>Reproducibility</h2>
<p>Driver: <code>script/bench_zi_count.R</code>. Tree: <code>ape::rtree(', r$n_species, ')</code>. Traits: <code>simulate_zi_count_traits()</code>. Training: ', r$epochs, ' epochs with early stopping. To reproduce: <code>Rscript script/bench_zi_count.R</code>, then <code>Rscript script/make_bench_zi_count_html.R</code>.</p>

<footer>
Source: <code>script/bench_zi_count.R</code> &middot;
Results: <code>script/bench_zi_count.rds</code> &middot;
Report: <code>pkgdown/assets/dev/bench_zi_count.html</code>
</footer>

</body>
</html>
')

targets <- c("script/bench_zi_count.html",
             "pkgdown/assets/dev/bench_zi_count.html")
for (t in targets) {
  dir.create(dirname(t), showWarnings = FALSE, recursive = TRUE)
  writeLines(html, t)
  cat("Wrote ", t, " (", file.size(t), " bytes)\n", sep = "")
}
