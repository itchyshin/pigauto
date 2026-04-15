#!/usr/bin/env Rscript
#
# script/make_bench_multi_proportion_html.R
#
# Build a self-contained HTML report from bench_multi_proportion.rds.

suppressPackageStartupMessages({ })

rds_path <- "script/bench_multi_proportion.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path, "\nRun script/bench_multi_proportion.R first.")
}

r   <- readRDS(rds_path)
res <- r$results
stopifnot(is.data.frame(res), nrow(res) > 0L)

test_df <- res[res$split == "test", ]
if (!nrow(test_df)) stop("No test-split rows in ", rds_path)

method_order  <- c("mean", "baseline", "pigauto")
method_label  <- c(mean = "Mean imputation (no phylo)",
                   baseline = "BM baseline (phylogeny only)",
                   pigauto = "pigauto (BM + GNN)")
method_colour <- c(mean = "#9ca3af", baseline = "#2563eb",
                   pigauto = "#059669")

scenarios_primary   <- r$scenarios_primary
secondary_scenarios <- r$secondary_scenarios
primary_frac        <- r$primary_frac

scen_label <- c(
  setNames(gsub("signal_", "Signal = ", scenarios_primary), scenarios_primary),
  setNames(gsub("K_",      "K = ",      secondary_scenarios), secondary_scenarios)
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "&ndash;", formatC(x, digits = digits, format = "f"))
}

avg_metric <- function(method, scenario, frac, metric) {
  sub <- test_df[test_df$method == method &
                   test_df$scenario == scenario &
                   test_df$missing_frac == frac, ]
  v <- sub[[metric]]
  v <- v[is.finite(v)]
  if (!length(v)) NA_real_ else mean(v)
}

bar_chart <- function(metric_col, y_label, title_text, scenarios) {
  W <- 680; H <- 320
  pad_l <- 68; pad_r <- 160; pad_t <- 40; pad_b <- 64
  plot_w <- W - pad_l - pad_r
  plot_h <- H - pad_t - pad_b
  n_groups <- length(scenarios)
  n_bars   <- length(method_order)
  group_w  <- plot_w / n_groups
  bar_w    <- group_w / (n_bars + 1)
  gap      <- bar_w / (n_bars + 1)

  vals <- matrix(NA_real_, n_bars, n_groups,
                 dimnames = list(method_order, scenarios))
  for (m in method_order) for (s in scenarios) {
    vals[m, s] <- avg_metric(m, s, primary_frac, metric_col)
  }
  v_fin <- vals[is.finite(vals)]
  if (!length(v_fin)) return("")

  y_lo <- 0
  y_hi <- max(v_fin) * 1.18
  if (y_hi <= y_lo) { y_lo <- 0; y_hi <- 1 }

  y_of <- function(v) {
    v <- pmin(pmax(v, y_lo), y_hi)
    pad_t + plot_h - (v - y_lo) / (y_hi - y_lo) * plot_h
  }

  p <- c()
  p <- c(p, sprintf(
    '<svg viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" style="width:100%%;max-width:680px;height:auto;font-family:system-ui,-apple-system,sans-serif">',
    W, H))
  p <- c(p, sprintf(
    '<text x="%d" y="24" font-size="14" font-weight="600" fill="#111827">%s</text>',
    pad_l, title_text))

  for (k in 0:4) {
    v <- y_lo + (y_hi - y_lo) * k / 4
    y <- y_of(v)
    p <- c(p, sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e7eb"/>',
                       pad_l, y, pad_l + plot_w, y))
    p <- c(p, sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="#6b7280">%s</text>',
                       pad_l - 6, y + 3, fmt(v, 2)))
  }
  p <- c(p, sprintf('<text x="12" y="%d" font-size="10" fill="#6b7280" transform="rotate(-90 12 %d)">%s</text>',
                     pad_t + plot_h / 2, pad_t + plot_h / 2, y_label))

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
      p <- c(p, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" rx="2"/>',
                         bx, by, bar_w, bh, method_colour[m]))
      p <- c(p, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="9" fill="#374151">%s</text>',
                         bx + bar_w / 2, by - 4, fmt(v, 3)))
    }
    p <- c(p, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="11" fill="#374151">%s</text>',
                       gx + group_w / 2, pad_t + plot_h + 22, scen_label[scen]))
  }

  lx <- pad_l + plot_w + 12
  for (b in seq_len(n_bars)) {
    m <- method_order[b]
    ly <- pad_t + 20 + (b - 1) * 20
    p <- c(p, sprintf('<rect x="%d" y="%d" width="14" height="10" fill="%s" rx="2"/>',
                       lx, ly, method_colour[m]))
    p <- c(p, sprintf('<text x="%d" y="%d" font-size="11" fill="#374151">%s</text>',
                       lx + 18, ly + 9, method_label[m]))
  }

  p <- c(p, '</svg>')
  paste(p, collapse = "\n")
}

chart_aitch <- bar_chart("aitchison",   "Aitchison distance",
                          "Aitchison distance by signal strength",
                          scenarios_primary)
chart_mae   <- bar_chart("simplex_mae", "Simplex MAE",
                          "Simplex MAE by signal strength",
                          scenarios_primary)
chart_K     <- if (length(secondary_scenarios) > 0L) {
  bar_chart("aitchison", "Aitchison distance",
            "Aitchison distance by number of components (signal = 0.6)",
            secondary_scenarios)
} else ""

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "&ndash;", paste0(formatC(100 * x, digits = digits, format = "f"), "%"))
}

rows <- character()
for (s in scenarios_primary) {
  ait      <- setNames(vapply(method_order, function(m)
    avg_metric(m, s, primary_frac, "aitchison"), numeric(1)), method_order)
  rmse_clr <- setNames(vapply(method_order, function(m)
    avg_metric(m, s, primary_frac, "rmse_clr"), numeric(1)), method_order)
  smae     <- setNames(vapply(method_order, function(m)
    avg_metric(m, s, primary_frac, "simplex_mae"), numeric(1)), method_order)
  argmax   <- setNames(vapply(method_order, function(m)
    avg_metric(m, s, primary_frac, "accuracy"), numeric(1)), method_order)

  # R-squared on z-scored CLR: 1 - rmse_clr^2 (null baseline has RMSE = 1 on z-scored data).
  # Negative R^2 means the method is WORSE than predicting the column mean.
  r2 <- setNames(1 - rmse_clr^2, method_order)

  # Lift over mean baseline on Aitchison distance: positive = better than mean.
  lift_bl <- (ait["mean"] - ait["baseline"]) / ait["mean"]
  lift_pg <- (ait["mean"] - ait["pigauto"])  / ait["mean"]

  best <- names(which.min(ait))
  star <- function(m, b) if (!is.null(b) && m == b)
    '<span style="color:#059669">&#9733;</span>' else ''
  rows <- c(rows, sprintf(
    paste0('<tr><td>%s</td>',
           '<td>%s</td><td>%s</td><td>%s</td>',                 # argmax acc
           '<td>%s%s</td><td>%s%s</td><td>%s%s</td>',             # Aitchison + winner
           '<td>%s</td><td>%s</td><td>%s</td>',                  # R^2 (CLR)
           '<td>%s</td><td>%s</td>',                              # lift (BM/pigauto vs mean)
           '<td>%s</td><td>%s</td><td>%s</td></tr>'),             # simplex MAE
    scen_label[s],
    fmt_pct(argmax["mean"]), fmt_pct(argmax["baseline"]), fmt_pct(argmax["pigauto"]),
    fmt(ait["mean"]),     star("mean", best),
    fmt(ait["baseline"]), star("baseline", best),
    fmt(ait["pigauto"]),  star("pigauto", best),
    fmt(r2["mean"], 2), fmt(r2["baseline"], 2), fmt(r2["pigauto"], 2),
    fmt_pct(lift_bl), fmt_pct(lift_pg),
    fmt(smae["mean"]), fmt(smae["baseline"]), fmt(smae["pigauto"])
  ))
}

ait_bl_hi <- avg_metric("baseline", "signal_1.0", primary_frac, "aitchison")
ait_pg_hi <- avg_metric("pigauto",  "signal_1.0", primary_frac, "aitchison")
ait_bl_lo <- avg_metric("baseline", "signal_0.2", primary_frac, "aitchison")
ait_pg_lo <- avg_metric("pigauto",  "signal_0.2", primary_frac, "aitchison")
v_hi <- if (is.finite(ait_bl_hi) && is.finite(ait_pg_hi))
  sprintf("At high signal (1.0), baseline Aitchison = %s and pigauto = %s (gate closes when baseline is near-optimal).",
          fmt(ait_bl_hi), fmt(ait_pg_hi))  else "At high signal: data not available."
v_lo <- if (is.finite(ait_bl_lo) && is.finite(ait_pg_lo))
  sprintf("At low signal (0.2), baseline Aitchison = %s and pigauto = %s.",
          fmt(ait_bl_lo), fmt(ait_pg_lo))  else "At low signal: data not available."

html <- paste0('<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Multi-proportion trait benchmark &mdash; pigauto</title>
<style>
  body{font-family:system-ui,-apple-system,sans-serif;line-height:1.5;color:#111827;max-width:900px;margin:2rem auto;padding:0 1rem}
  h1,h2{color:#111827}
  h1{border-bottom:2px solid #059669;padding-bottom:.3rem}
  h2{margin-top:2rem;color:#059669}
  table{border-collapse:collapse;margin:1rem 0;font-size:14px;width:100%}
  th{background:#f3f4f6;text-align:left;padding:.5rem .75rem;border-bottom:2px solid #d1d5db;font-weight:600}
  td{padding:.45rem .75rem;border-bottom:1px solid #e5e7eb}
  .summary-card{background:#f9fafb;border-left:4px solid #059669;padding:1rem 1.2rem;margin:1.2rem 0;border-radius:0 4px 4px 0}
  code{background:#f3f4f6;padding:.1rem .35rem;border-radius:3px;font-size:90%}
</style>
</head>
<body>

<h1>Multi-proportion trait benchmark</h1>

<div class="summary-card">
  <strong>Trait type:</strong> K-category compositional data (rows sum to 1).<br>
  <strong>Encoding:</strong> Centred log-ratio (CLR) + per-component z-score.<br>
  <strong>Baseline:</strong> Brownian motion on CLR space (K independent BM fits).<br>
  <strong>GNN output:</strong> K logits &rarr; sum-to-zero projection &rarr; softmax &rarr; simplex.<br>
  <strong>Species:</strong> ', r$n_species, '.
  <strong>Reps:</strong> ', r$n_reps, '.
  <strong>Missing fraction:</strong> ', primary_frac, '.
  <strong>Epochs:</strong> ', r$epochs, '.
  <strong>Wall time:</strong> ', sprintf("%.1f", r$total_wall / 60), ' min.
</div>

<h2>Metrics</h2>
<p>Each held-out row has its entire composition (all K values) masked together, because
compositional data is either complete for a row or missing &mdash; you cannot observe
3 of 5 components without the rest.</p>
<ul>
  <li><strong>Aitchison distance</strong>: Euclidean distance in CLR space; the natural metric
      for compositional data. Lower is better.</li>
  <li><strong>RMSE on CLR</strong>: root-mean-square error on z-scored CLR latent values;
      directly comparable to continuous-trait RMSE.</li>
  <li><strong>Simplex MAE</strong>: mean absolute error on the simplex after softmax decode;
      interpretable in original proportion units.</li>
</ul>

<h2>Primary sweep: phylogenetic signal (K = 5)</h2>
', chart_aitch, '

<h2>Summary table (test split, averaged across reps)</h2>
<table>
<thead>
<tr>
  <th rowspan="2">Scenario</th>
  <th colspan="3">Argmax accuracy</th>
  <th colspan="3">Aitchison distance</th>
  <th colspan="3">R&sup2; (CLR)</th>
  <th colspan="2">Lift vs mean</th>
  <th colspan="3">Simplex MAE</th>
</tr>
<tr>
  <th>mean</th><th>BM</th><th>pigauto</th>
  <th>mean</th><th>BM</th><th>pigauto</th>
  <th>mean</th><th>BM</th><th>pigauto</th>
  <th>BM</th><th>pigauto</th>
  <th>mean</th><th>BM</th><th>pigauto</th>
</tr>
</thead>
<tbody>
', paste(rows, collapse = "\n"), '
</tbody></table>
<p><small>
&#9733; = best method on Aitchison distance for that scenario.
<strong>Argmax accuracy</strong>: fraction of held-out rows where predicted dominant component matches truth.
<strong>R&sup2; (CLR)</strong>: 1 &minus; RMSE&sup2; on z-scored CLR; negative means worse than column-mean null.
<strong>Lift vs mean</strong>: % reduction in Aitchison distance relative to naive column-mean.
</small></p>

', chart_mae, '

<h2>Secondary sweep: K (signal = 0.6)</h2>
', chart_K, '

<h2>Commentary</h2>
<p>', v_hi, '</p>
<p>', v_lo, '</p>
<p>The gated safety in pigauto ensures that when the BM-on-CLR baseline is already optimal
(high phylogenetic signal), the calibrated gate closes and pigauto returns the baseline
&mdash; so it never degrades accuracy. When signal is low, the GNN can learn cross-component
correlations that the independent-per-component BM misses.</p>

<h2>How to reproduce</h2>
<pre><code>cd pigauto &amp;&amp; Rscript script/bench_multi_proportion.R
Rscript script/make_bench_multi_proportion_html.R</code></pre>

<p><small>Generated from <code>script/bench_multi_proportion.rds</code> on ',
format(Sys.time()), '. Git commit: <code>',
substr(r$commit, 1, 7), '</code>.</small></p>

</body></html>
')

out1 <- "script/bench_multi_proportion.html"
writeLines(html, out1)
cat("Wrote ", out1, " (", nchar(html), " bytes)\n", sep = "")

out2 <- "pkgdown/assets/dev/bench_multi_proportion.html"
if (dir.exists(dirname(out2))) {
  writeLines(html, out2)
  cat("Wrote ", out2, " (", nchar(html), " bytes)\n", sep = "")
}
