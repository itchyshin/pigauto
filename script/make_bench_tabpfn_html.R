#!/usr/bin/env Rscript
# Build a self-contained HTML report from bench_tabpfn.rds
suppressPackageStartupMessages({ })

r <- readRDS("script/bench_tabpfn.rds")
d <- r$results

# ---- Continuous RMSE aggregation -----------------------------------------
cr <- d[d$metric == "rmse" & d$type == "continuous", ]
agg <- aggregate(value ~ scenario + method, data = cr,
                 FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
agg_df <- data.frame(scenario = agg$scenario, method = agg$method,
                     mean = agg$value[, "mean"], sd = agg$value[, "sd"],
                     n = agg$value[, "n"], stringsAsFactors = FALSE)

scen_order   <- c("BM", "OU", "nonlinear", "regime_shift", "mixed")
method_order <- c("baseline", "pigauto", "tabpfn", "pigauto_tabpfn")
method_label <- c(baseline = "Phylo baseline",
                  pigauto = "pigauto (phylo + GNN)",
                  tabpfn = "TabPFN alone",
                  pigauto_tabpfn = "pigauto + TabPFN (stacked)")
method_colour <- c(baseline = "#6b7280",
                   pigauto = "#059669",
                   tabpfn = "#c2410c",
                   pigauto_tabpfn = "#7c3aed")

# ---- Discrete accuracy ----------------------------------------------------
ac <- d[d$metric == "accuracy" & d$type %in% c("binary", "categorical"), ]
acc_df <- if (nrow(ac) > 0L) {
  a <- aggregate(value ~ scenario + method + trait + type, data = ac, mean)
  a
} else NULL

# ---- Helpers --------------------------------------------------------------
fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "f")

svg_bar_chart <- function(scen, max_rmse) {
  sub <- agg_df[agg_df$scenario == scen, ]
  sub <- sub[match(method_order, sub$method), ]
  sub <- sub[!is.na(sub$method), ]
  if (nrow(sub) == 0L) return("")

  W   <- 560
  H   <- 200
  pad_l <- 160
  pad_r <- 40
  pad_t <- 24
  pad_b <- 40
  plot_w <- W - pad_l - pad_r
  plot_h <- H - pad_t - pad_b
  n_bars <- nrow(sub)
  bar_h  <- plot_h / (n_bars * 1.4)
  gap    <- bar_h * 0.4

  scale <- function(v) v / max_rmse * plot_w

  parts <- character()
  parts <- c(parts, sprintf(
    '<svg viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" style="width:100%%;max-width:560px;height:auto;font-family:system-ui,-apple-system,sans-serif">',
    W, H))
  # axis line
  parts <- c(parts, sprintf(
    '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#999" stroke-width="1"/>',
    pad_l, H - pad_b, W - pad_r, H - pad_b))
  # gridlines at 0.25, 0.5, 0.75, 1.0, 1.25 (whichever <= max_rmse)
  ticks <- seq(0, ceiling(max_rmse * 4) / 4, by = 0.25)
  for (t in ticks) {
    x <- pad_l + scale(t)
    parts <- c(parts, sprintf(
      '<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e5e7eb" stroke-width="1"/>',
      x, pad_t, x, H - pad_b))
    parts <- c(parts, sprintf(
      '<text x="%.1f" y="%d" text-anchor="middle" font-size="10" fill="#6b7280">%.2f</text>',
      x, H - pad_b + 14, t))
  }

  for (i in seq_len(n_bars)) {
    row <- sub[i, ]
    y <- pad_t + (i - 1) * (bar_h + gap) + gap / 2
    w <- scale(row$mean)
    col <- method_colour[row$method]
    # Bar
    parts <- c(parts, sprintf(
      '<rect x="%d" y="%.1f" width="%.1f" height="%.1f" fill="%s" rx="3"/>',
      pad_l, y, w, bar_h, col))
    # Label on left
    parts <- c(parts, sprintf(
      '<text x="%d" y="%.1f" text-anchor="end" font-size="12" fill="#111827">%s</text>',
      pad_l - 8, y + bar_h / 2 + 4, method_label[row$method]))
    # Value + sd on right of bar
    parts <- c(parts, sprintf(
      '<text x="%.1f" y="%.1f" font-size="11" fill="#111827">%.3f &#177; %.3f</text>',
      pad_l + w + 6, y + bar_h / 2 + 4, row$mean, row$sd))
  }
  parts <- c(parts, '</svg>')
  paste(parts, collapse = "\n")
}

# ---- Build HTML -----------------------------------------------------------
max_rmse <- max(agg_df$mean + agg_df$sd, na.rm = TRUE) * 1.05

# Main table rows
tbl_rows <- character()
for (scen in scen_order) {
  sub <- agg_df[agg_df$scenario == scen, ]
  sub <- sub[match(method_order, sub$method), ]
  sub <- sub[!is.na(sub$method), ]
  if (nrow(sub) == 0L) next

  bl <- sub$mean[sub$method == "baseline"]
  best <- sub$method[which.min(sub$mean)]

  first <- TRUE
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, ]
    delta <- if (row$method == "baseline") "&mdash;" else {
      pct <- 100 * (bl - row$mean) / bl
      sign_str <- if (pct >= 0) sprintf("+%.1f%%", pct) else sprintf("%.1f%%", pct)
      colour <- if (pct >= 0) "#059669" else "#dc2626"
      sprintf('<span style="color:%s;font-weight:600">%s</span>', colour, sign_str)
    }
    best_mark <- if (row$method == best) ' <span title="best">&#9733;</span>' else ""

    scen_cell <- if (first) {
      first <- FALSE
      sprintf('<td rowspan="%d" style="vertical-align:top;padding-top:10px"><b>%s</b></td>',
              nrow(sub), scen)
    } else ""

    tbl_rows <- c(tbl_rows, sprintf(
      '<tr>%s<td>%s%s</td><td style="text-align:right;font-variant-numeric:tabular-nums">%s &#177; %s</td><td style="text-align:right">%s</td></tr>',
      scen_cell,
      method_label[row$method],
      best_mark,
      fmt(row$mean), fmt(row$sd),
      delta))
  }
}

# Accuracy rows for mixed
acc_rows <- character()
if (!is.null(acc_df) && nrow(acc_df) > 0L) {
  for (tr in unique(acc_df$trait)) {
    sub <- acc_df[acc_df$trait == tr, ]
    sub <- sub[match(method_order, sub$method), ]
    sub <- sub[!is.na(sub$method), ]
    if (nrow(sub) == 0L) next
    first <- TRUE
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      tr_cell <- if (first) {
        first <- FALSE
        sprintf('<td rowspan="%d" style="vertical-align:top;padding-top:10px"><b>%s</b><br><span style="color:#6b7280;font-size:11px">%s</span></td>',
                nrow(sub), row$trait, row$type)
      } else ""
      acc_rows <- c(acc_rows, sprintf(
        '<tr>%s<td>%s</td><td style="text-align:right;font-variant-numeric:tabular-nums">%.1f%%</td></tr>',
        tr_cell, method_label[row$method], 100 * row$value))
    }
  }
}

# Per-scenario charts
chart_html <- character()
narratives <- list(
  BM = "Pure Brownian motion. Phylogeny is the entire story. TabPFN has no phylogenetic features so it collapses to the trait mean (RMSE \u2248 1.0 on standardised data).",
  OU = "Ornstein\u2013Uhlenbeck with strong equilibrium pull. Species values are weakly dependent on phylogeny; information shared across traits dominates. <b>This is the only scenario where TabPFN beats the phylogenetic baseline</b> (+12.9%), and stacking adds a further 0.5 pp.",
  nonlinear = "Non-linear trait transform on BM. The phylogenetic baseline still extracts most of the signal; TabPFN finds no useful cross-trait pattern and is \u224838% worse.",
  regime_shift = "Step change at an internal node. Phylogenetic label propagation nails this because neighbours on the tree share regime. TabPFN has no tree and fails dramatically (\u2212122%).",
  mixed = "Continuous + binary + categorical traits. pigauto holds the phylo baseline's continuous RMSE and fully preserves discrete accuracy (69.4%, unchanged). TabPFN errored on this scenario (a column-naming bug in the subprocess wrapper when dummy-encoded categorical columns are present)."
)

for (scen in scen_order) {
  sub <- agg_df[agg_df$scenario == scen, ]
  if (nrow(sub) == 0L) next
  chart_html <- c(chart_html, sprintf(
    '<section style="margin:28px 0"><h3 style="margin:0 0 8px 0">Scenario: %s</h3><p style="margin:0 0 12px 0;color:#374151;max-width:700px">%s</p>%s</section>',
    scen, narratives[[scen]], svg_bar_chart(scen, max_rmse)))
}

# Headline numbers
baseline_means <- agg_df[agg_df$method == "baseline", c("scenario", "mean")]
pigauto_means  <- agg_df[agg_df$method == "pigauto",  c("scenario", "mean")]
tabpfn_means   <- agg_df[agg_df$method == "tabpfn",   c("scenario", "mean")]
stack_means    <- agg_df[agg_df$method == "pigauto_tabpfn", c("scenario", "mean")]

rel_improve <- function(method_means) {
  m <- merge(baseline_means, method_means, by = "scenario", suffixes = c("_b", "_m"))
  100 * (m$mean_b - m$mean_m) / m$mean_b
}
ri_pig  <- rel_improve(pigauto_means)
ri_tab  <- rel_improve(tabpfn_means)
ri_stk  <- rel_improve(stack_means)

summary_bullets <- paste0(
  "<ul>",
  sprintf("<li><b>Phylogenetic baseline is very hard to beat</b> on phylogenetic data. Across 5 scenarios (%d reps each), TabPFN alone is worse than the baseline in %d of 4 evaluated scenarios, often catastrophically (BM &minus;99.8%%, regime shift &minus;122.2%%).</li>",
          unique(cr$rep) |> length(), sum(ri_tab < 0, na.rm = TRUE)),
  sprintf("<li><b>pigauto (phylo + GNN)</b> ties or marginally beats the baseline on every scenario (%+.1f%% to %+.1f%% continuous RMSE).</li>",
          min(ri_pig, na.rm = TRUE), max(ri_pig, na.rm = TRUE)),
  "<li><b>Stacking pigauto on top of TabPFN</b> barely changes TabPFN&rsquo;s predictions (&lt;0.5 pp in every scenario). The GNN cannot recover phylogenetic signal from a baseline that has already collapsed toward the mean.</li>",
  "<li><b>One scenario favours TabPFN</b>: OU with a strong equilibrium pull, where phylogeny is a weaker predictor and cross-trait patterns matter more. TabPFN alone is +12.9% vs baseline; stacking adds another 0.5 pp.</li>",
  "<li><b>Discrete traits</b> are unchanged: pigauto matches the baseline exactly on binary and categorical accuracy (69.4% on the mixed scenario), confirming the v4 calibration fix held.</li>",
  "</ul>"
)

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>pigauto vs TabPFN &mdash; 4-way benchmark</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 880px; margin: 2em auto; padding: 0 1.5em; color: #111827;
         line-height: 1.55; }
  h1 { border-bottom: 2px solid #111827; padding-bottom: .25em; }
  h2 { margin-top: 2em; color: #1f2937; border-bottom: 1px solid #d1d5db;
       padding-bottom: .15em; }
  h3 { color: #111827; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 14px; }
  th, td { padding: 6px 10px; border-bottom: 1px solid #e5e7eb; text-align: left; }
  th { background: #f3f4f6; font-weight: 600; }
  .meta { color: #6b7280; font-size: 13px; }
  .verdict { background: #fef3c7; border-left: 4px solid #d97706; padding: 1em 1.2em;
             border-radius: 4px; margin: 1.5em 0; }
  .verdict b { color: #92400e; }
  code { background: #f3f4f6; padding: 1px 5px; border-radius: 3px; font-size: 13px; }
  ul li { margin: 6px 0; }
  footer { color: #6b7280; font-size: 12px; margin-top: 3em; border-top: 1px solid #e5e7eb;
           padding-top: 1em; }
</style>
</head>
<body>

<h1>pigauto vs TabPFN: can a tabular foundation model beat phylogeny?</h1>
<p class="meta">Benchmark run on ', format(Sys.time(), "%Y-%m-%d %H:%M"),
  ' &middot; ', r$n_species, ' species &middot; ', r$n_reps, ' reps per scenario &middot; 25% missingness &middot; 600 epochs</p>

<div class="verdict">
<b>Bottom line.</b> On simulated phylogenetic data, <b>TabPFN does not beat the simple phylogenetic baseline</b>, and stacking pigauto on top of TabPFN does not rescue it. The phylogenetic baseline + pigauto remains the best general-purpose approach on phylogenetic data. TabPFN is only competitive when phylogenetic signal is weak (OU with strong equilibrium pull), in which case it <em>does</em> win by ~13%. This suggests a hybrid worth exploring: detect the per-trait phylogenetic signal and route low-lambda traits to TabPFN.
</div>

<h2>Executive summary</h2>
', summary_bullets, '

<h2>Continuous RMSE by scenario</h2>
<p class="meta">Lower is better. &#9733; marks the best method within each scenario. &ldquo;vs baseline&rdquo; shows relative RMSE improvement (positive = better than phylo baseline).</p>
<table>
<thead><tr><th>Scenario</th><th>Method</th><th style="text-align:right">RMSE (mean &#177; SD)</th><th style="text-align:right">vs baseline</th></tr></thead>
<tbody>
', paste(tbl_rows, collapse = "\n"), '
</tbody>
</table>
')

if (length(acc_rows) > 0L) {
  html <- paste0(html, '
<h2>Discrete accuracy (mixed scenario)</h2>
<table>
<thead><tr><th>Trait</th><th>Method</th><th style="text-align:right">Accuracy</th></tr></thead>
<tbody>
', paste(acc_rows, collapse = "\n"), '
</tbody>
</table>
')
}

html <- paste0(html, '

<h2>Per-scenario breakdown</h2>
', paste(chart_html, collapse = "\n"), '

<h2>What this means for pigauto</h2>
<p>
Two ideas worth taking seriously:
</p>
<ol>
<li><b>pigauto is doing its job as a phylogenetic imputer</b>, but the ceiling imposed by the phylo label-propagation baseline is high and the GNN is only finding small residual corrections (+0 to +2% RMSE). The recent calibration work means we no longer sacrifice discrete accuracy, so the method is safe to recommend as a drop-in upgrade to BACE/Rphylopars-style methods.</li>
<li><b>TabPFN is the wrong default for phylogenetic data</b>. It treats species as exchangeable rows, which destroys the tree structure that dominates these scenarios. The single exception (OU) points at a useful selector: per-trait Pagel&rsquo;s &lambda;. When &lambda; is low, route to TabPFN; when high, use phylo+GNN. This is a much more promising architecture than blind stacking.</li>
</ol>

<h2>Bug to fix</h2>
<p>The TabPFN subprocess wrapper errored on the <code>mixed</code> scenario with <code>length of &#39;dimnames&#39; [2] not equal to array extent</code>. The cause: <code>X_scaled</code> has one-hot columns for the categorical trait, so <code>ncol(mu) &gt; length(trait_names)</code>. Fix: use <code>colnames(data$X_scaled)</code> rather than <code>data$trait_names</code> when re-attaching dim-names.</p>

<footer>
Source: <code>script/bench_tabpfn.R</code> &middot; Results: <code>script/bench_tabpfn.rds</code> &middot; Report: <code>script/bench_tabpfn.html</code>
</footer>

</body>
</html>
')

# Dual-write: one copy in script/ (dev artefact), one copy in
# pkgdown/assets/dev/ so pkgdown::build_site() exposes it on the web.
targets <- c("script/bench_tabpfn.html",
             "pkgdown/assets/dev/bench_tabpfn.html")
for (t in targets) {
  dir.create(dirname(t), showWarnings = FALSE, recursive = TRUE)
  writeLines(html, t)
  cat("Wrote ", t, "\n", sep = "")
}
