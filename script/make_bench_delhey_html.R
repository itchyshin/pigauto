#!/usr/bin/env Rscript
#
# script/make_bench_delhey_html.R
#
# Build a self-contained HTML report from bench_delhey.rds.
#
# Source data:  script/bench_delhey.rds
# Writes two identical copies:
#   script/bench_delhey.html
#   pkgdown/assets/dev/bench_delhey.html
#
# Run with
#   /usr/local/bin/Rscript script/make_bench_delhey_html.R

suppressPackageStartupMessages({ })

rds_path <- "script/bench_delhey.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path,
       "\nRun script/bench_delhey.R first.")
}

r <- readRDS(rds_path)
res <- r$results
meta <- r$meta
stopifnot(is.data.frame(res), nrow(res) > 0L)

# Keep test-split rows only
test_df <- res[res$split == "test", ]
if (!nrow(test_df)) stop("No test-split rows in ", rds_path)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), "&ndash;", formatC(x, digits = digits, format = "f"))
}
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "&ndash;", sprintf("%.*f%%", digits, 100 * x))
}

method_order  <- c("mean", "BM_only", "pigauto", "pigauto_covs")
method_label  <- c(
  mean          = "Mean imputation",
  BM_only       = "BM baseline",
  pigauto       = "pigauto (no covariates)",
  pigauto_covs  = "pigauto + covariates"
)
method_colour <- c(
  mean          = "#9ca3af",
  BM_only       = "#2563eb",
  pigauto       = "#7c3aed",
  pigauto_covs  = "#059669"
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
h('<title>pigauto benchmark: Delhey plumage lightness (covariates)</title>')
h('<style>')
h('body{font-family:Inter,system-ui,sans-serif;max-width:960px;margin:2rem auto;',
  'padding:0 1rem;color:#1f2937;line-height:1.6}')
h('h1{color:#059669;border-bottom:2px solid #059669;padding-bottom:.3rem}')
h('h2{color:#374151;margin-top:2rem}')
h('h3{color:#6b7280}')
h('table{border-collapse:collapse;width:100%;margin:1rem 0}')
h('th,td{border:1px solid #d1d5db;padding:6px 10px;text-align:right}')
h('th{background:#f3f4f6;text-align:center}')
h('td:first-child,th:first-child{text-align:left}')
h('.best{font-weight:700;color:#059669}')
h('.note{font-size:.85rem;color:#6b7280;margin-top:.5rem}')
h('</style>')
h('</head><body>')

h('<h1>Delhey plumage lightness benchmark</h1>')
h('<p>Environmental covariates for trait imputation on <strong>',
  meta$n_species, '</strong> passerine species (Delhey et al. 2019).</p>')
h('<p>Traits: <code>', paste(meta$trait_cols, collapse = ', '), '</code></p>')
h('<p>Covariates: <code>', paste(meta$cov_cols, collapse = ', '), '</code></p>')

# ---- RMSE table per missingness level --------------------------------------
h('<h2>Test-set RMSE</h2>')

for (frac in meta$miss_fracs) {
  h('<h3>Missing: ', 100 * frac, '%</h3>')

  sub <- test_df[test_df$missing_frac == frac, ]
  traits <- unique(sub$trait)

  h('<table>')
  h('<tr><th>Trait</th>')
  for (m in method_order) {
    if (m %in% sub$method) {
      h('<th style="color:', method_colour[m], '">', method_label[m], '</th>')
    }
  }
  h('</tr>')

  for (tr in traits) {
    h('<tr><td>', tr, '</td>')
    vals <- numeric(0)
    for (m in method_order) {
      rows <- sub[sub$trait == tr & sub$method == m, ]
      if (nrow(rows) > 0) {
        v <- mean(rows$rmse, na.rm = TRUE)
        vals <- c(vals, v)
      }
    }
    best_val <- if (length(vals) > 0 && any(is.finite(vals))) min(vals[is.finite(vals)]) else Inf

    for (m in method_order) {
      rows <- sub[sub$trait == tr & sub$method == m, ]
      if (nrow(rows) > 0) {
        v <- mean(rows$rmse, na.rm = TRUE)
        cls <- if (is.finite(v) && is.finite(best_val) && abs(v - best_val) < 1e-8) ' class="best"' else ''
        h('<td', cls, '>', fmt(v), '</td>')
      }
    }
    h('</tr>')
  }
  h('</table>')
}

# ---- Pearson r table -------------------------------------------------------
h('<h2>Test-set Pearson r</h2>')

for (frac in meta$miss_fracs) {
  h('<h3>Missing: ', 100 * frac, '%</h3>')

  sub <- test_df[test_df$missing_frac == frac, ]
  traits <- unique(sub$trait)

  h('<table>')
  h('<tr><th>Trait</th>')
  for (m in method_order) {
    if (m %in% sub$method) {
      h('<th style="color:', method_colour[m], '">', method_label[m], '</th>')
    }
  }
  h('</tr>')

  for (tr in traits) {
    h('<tr><td>', tr, '</td>')
    vals <- numeric(0)
    for (m in method_order) {
      rows <- sub[sub$trait == tr & sub$method == m, ]
      if (nrow(rows) > 0) {
        v <- mean(rows$pearson_r, na.rm = TRUE)
        vals <- c(vals, v)
      }
    }
    best_val <- if (length(vals) > 0 && any(is.finite(vals))) max(vals[is.finite(vals)]) else -Inf

    for (m in method_order) {
      rows <- sub[sub$trait == tr & sub$method == m, ]
      if (nrow(rows) > 0) {
        v <- mean(rows$pearson_r, na.rm = TRUE)
        cls <- if (is.finite(v) && is.finite(best_val) && abs(v - best_val) < 1e-8) ' class="best"' else ''
        h('<td', cls, '>', fmt(v), '</td>')
      }
    }
    h('</tr>')
  }
  h('</table>')
}

# ---- Covariate lift summary ------------------------------------------------
h('<h2>Covariate lift summary</h2>')
h('<p>RMSE ratio of <span style="color:#059669">pigauto + covariates</span> ',
  'relative to <span style="color:#7c3aed">pigauto (no covariates)</span>. ',
  'Values &lt; 1 indicate covariates help.</p>')

rmse_pg <- aggregate(rmse ~ missing_frac + trait,
                     data = test_df[test_df$method == "pigauto", ],
                     FUN = mean, na.rm = TRUE)
rmse_pc <- aggregate(rmse ~ missing_frac + trait,
                     data = test_df[test_df$method == "pigauto_covs", ],
                     FUN = mean, na.rm = TRUE)

if (nrow(rmse_pg) > 0 && nrow(rmse_pc) > 0) {
  merged <- merge(rmse_pg, rmse_pc, by = c("missing_frac", "trait"),
                  suffixes = c("_nocov", "_cov"))
  merged$ratio <- merged$rmse_cov / merged$rmse_nocov

  h('<table>')
  h('<tr><th>Missing %</th><th>Trait</th><th>RMSE (no cov)</th>',
    '<th>RMSE (with cov)</th><th>Ratio</th></tr>')

  for (i in seq_len(nrow(merged))) {
    cls <- if (merged$ratio[i] < 0.98) ' class="best"' else ''
    h('<tr>')
    h('<td>', 100 * merged$missing_frac[i], '%</td>')
    h('<td>', merged$trait[i], '</td>')
    h('<td>', fmt(merged$rmse_nocov[i]), '</td>')
    h('<td>', fmt(merged$rmse_cov[i]), '</td>')
    h('<td', cls, '>', fmt(merged$ratio[i], 3), '</td>')
    h('</tr>')
  }
  h('</table>')
}

# ---- Footer ----------------------------------------------------------------
h('<p class="note">Wall time: ', round(meta$wall_time / 60, 1), ' min. ',
  'Generated: ', format(meta$timestamp, "%Y-%m-%d %H:%M"), '</p>')
h('<p class="note">Reference: Delhey K, Dale J, Valcu M, Kempenaers B (2019). ',
  '<em>Ecology Letters</em>, 22(5): 726-736.</p>')
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
