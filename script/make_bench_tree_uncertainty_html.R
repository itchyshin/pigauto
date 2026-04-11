#!/usr/bin/env Rscript
#
# script/make_bench_tree_uncertainty_html.R
#
# Build a self-contained HTML report from bench_tree_uncertainty.rds.
#
# Source data:  script/bench_tree_uncertainty.rds
# Writes two identical copies:
#   script/bench_tree_uncertainty.html
#   pkgdown/assets/dev/bench_tree_uncertainty.html
#
# Run with
#   /usr/local/bin/Rscript script/make_bench_tree_uncertainty_html.R

suppressPackageStartupMessages({ })

rds_path <- "script/bench_tree_uncertainty.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path,
       "\nRun script/bench_tree_uncertainty.R first.")
}

r <- readRDS(rds_path)
res <- r$results
meta <- r$meta
stopifnot(is.data.frame(res), nrow(res) > 0L)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), "&ndash;", formatC(x, digits = digits, format = "f"))
}
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "&ndash;", sprintf("%.*f%%", digits, 100 * x))
}

method_label <- c(
  single_tree = "Single-tree MI",
  multi_tree  = "Multi-tree MI (50 trees)"
)
method_colour <- c(
  single_tree = "#2563eb",
  multi_tree  = "#059669"
)

# ---------------------------------------------------------------------------
# Build HTML
# ---------------------------------------------------------------------------

html <- character()
h <- function(...) html <<- c(html, paste0(...))

h('<!DOCTYPE html>')
h('<html lang="en"><head>')
h('<meta charset="utf-8">')
h('<meta name="viewport" content="width=device-width, initial-scale=1">')
h('<title>pigauto benchmark: tree uncertainty</title>')
h('<style>')
h('body{font-family:Inter,system-ui,sans-serif;max-width:900px;margin:2rem auto;',
  'padding:0 1rem;color:#1f2937;line-height:1.6}')
h('h1{color:#059669;border-bottom:2px solid #059669;padding-bottom:.3rem}')
h('h2{color:#374151;margin-top:2rem}')
h('table{border-collapse:collapse;width:100%;margin:1rem 0}')
h('th,td{border:1px solid #d1d5db;padding:6px 10px;text-align:right}')
h('th{background:#f3f4f6;text-align:center}')
h('td:first-child,th:first-child{text-align:left}')
h('.best{font-weight:700;color:#059669}')
h('.note{font-size:.85rem;color:#6b7280;margin-top:.5rem}')
h('</style>')
h('</head><body>')

h('<h1>Tree uncertainty benchmark</h1>')
h('<p>Comparison of single-tree MI vs multi-tree MI on <strong>avonet300</strong> ',
  '(', meta$n_species, ' species, mixed-type traits).</p>')
h('<p>Downstream model: <code>', meta$model_formula, '</code></p>')
h('<p>Single-tree: MCC tree, m = ', meta$m_single, ' MC-dropout draws.<br>')
h('Multi-tree: ', meta$n_trees, ' posterior trees &times; ',
  meta$m_per_tree, ' draws = ', meta$n_trees * meta$m_per_tree, ' total datasets.</p>')

# ---- Main results table ----------------------------------------------------
h('<h2>Pooled estimates (Rubin\'s rules)</h2>')

for (frac in meta$miss_fracs) {
  h('<h3>Missing fraction: ', 100 * frac, '%</h3>')

  sub <- res[res$missing_frac == frac, ]
  terms <- unique(sub$term)

  h('<table>')
  h('<tr><th>Term</th><th>Method</th><th>Estimate</th><th>SE</th>',
    '<th>df</th><th>FMI</th><th>RIV</th></tr>')

  for (tm in terms) {
    for (meth in c("single_tree", "multi_tree")) {
      rows <- sub[sub$term == tm & sub$method == meth, ]
      if (nrow(rows) == 0) next
      est <- mean(rows$estimate, na.rm = TRUE)
      se  <- mean(rows$std.error, na.rm = TRUE)
      df  <- mean(rows$df, na.rm = TRUE)
      fmi <- mean(rows$fmi, na.rm = TRUE)
      riv <- mean(rows$riv, na.rm = TRUE)

      h('<tr>')
      h('<td>', tm, '</td>')
      h('<td style="color:', method_colour[meth], '">', method_label[meth], '</td>')
      h('<td>', fmt(est), '</td>')
      h('<td>', fmt(se), '</td>')
      h('<td>', fmt(df, 1), '</td>')
      h('<td>', fmt(fmi, 3), '</td>')
      h('<td>', fmt(riv, 3), '</td>')
      h('</tr>')
    }
  }
  h('</table>')
}

# ---- SE ratio table --------------------------------------------------------
h('<h2>SE ratio (multi-tree / single-tree)</h2>')
h('<p>Values > 1 indicate that phylogenetic uncertainty widens standard errors ',
  'beyond what single-tree MI provides.</p>')

h('<table>')
h('<tr><th>Missing %</th><th>Term</th><th>SE (single)</th>',
  '<th>SE (multi)</th><th>Ratio</th></tr>')

for (frac in meta$miss_fracs) {
  sub <- res[res$missing_frac == frac, ]
  terms <- unique(sub$term)
  for (tm in terms) {
    se_s <- mean(sub$std.error[sub$term == tm & sub$method == "single_tree"],
                 na.rm = TRUE)
    se_m <- mean(sub$std.error[sub$term == tm & sub$method == "multi_tree"],
                 na.rm = TRUE)
    ratio <- se_m / se_s

    cls <- if (ratio > 1.05) ' class="best"' else ''
    h('<tr>')
    h('<td>', 100 * frac, '%</td>')
    h('<td>', tm, '</td>')
    h('<td>', fmt(se_s), '</td>')
    h('<td>', fmt(se_m), '</td>')
    h('<td', cls, '>', fmt(ratio, 2), '</td>')
    h('</tr>')
  }
}
h('</table>')

# ---- Footer ----------------------------------------------------------------
h('<p class="note">Wall time: ', round(meta$wall_time / 60, 1), ' min. ',
  'Generated: ', format(meta$timestamp, "%Y-%m-%d %H:%M"), '</p>')
h('<p class="note">Reference: Nakagawa S, de Villemereuil P (2019). ',
  '<em>Systematic Biology</em> 68(4): 632-641.</p>')
h('</body></html>')

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------

html_text <- paste(html, collapse = "\n")

out_script  <- sub("\\.rds$", ".html", rds_path)
out_pkgdown <- file.path("pkgdown", "assets", "dev",
                         basename(out_script))

writeLines(html_text, out_script)
dir.create(dirname(out_pkgdown), recursive = TRUE, showWarnings = FALSE)
writeLines(html_text, out_pkgdown)

cat("Wrote", out_script, "(", file.size(out_script), "bytes)\n")
cat("Wrote", out_pkgdown, "(", file.size(out_pkgdown), "bytes)\n")
