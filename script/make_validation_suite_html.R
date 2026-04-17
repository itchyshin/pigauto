#!/usr/bin/env Rscript
#
# script/make_validation_suite_html.R
#
# Generate pkgdown/assets/validation_suite.html — an index of all 15
# benchmark results for pigauto v0.9.0.
#
# Run with:
#   Rscript script/make_validation_suite_html.R
#
# Reads: script/bench_*.rds
# Writes: pkgdown/assets/validation_suite.html
#
# If a .rds file is missing (script not yet run), that row shows "pending".

suppressPackageStartupMessages({})

here <- normalizePath(dirname(dirname(
  if (nchar(Sys.getenv("R_FILE")) > 0) Sys.getenv("R_FILE")
  else {
    args <- commandArgs(trailingOnly = FALSE)
    m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
    if (length(m) > 0) m[1] else "script/make_validation_suite_html.R"
  }
)), mustWork = FALSE)
if (!file.exists(file.path(here, "DESCRIPTION")))
  here <- normalizePath(".", mustWork = FALSE)

out_html_1 <- file.path(here, "pkgdown", "assets", "validation_suite.html")
out_html_2 <- file.path(here, "script", "validation_suite.html")

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

load_rds <- function(name) {
  path <- file.path(here, "script", paste0(name, ".rds"))
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

fmt <- function(x, d = 3) {
  if (is.null(x) || is.na(x) || !is.finite(x)) return("&ndash;")
  formatC(x, digits = d, format = "f")
}

pct <- function(x) {
  if (is.null(x) || is.na(x) || !is.finite(x)) return("&ndash;")
  sprintf("%.1f%%", x * 100)
}

lift_pct <- function(baseline_rmse, pigauto_rmse) {
  if (any(is.na(c(baseline_rmse, pigauto_rmse)))) return(NA_real_)
  (baseline_rmse - pigauto_rmse) / baseline_rmse
}

best_rmse <- function(df, method_val, scenario_val = NULL, frac_val = NULL,
                      split_val = "test", metric = "rmse") {
  sub <- df[df$method == method_val, ]
  if (!is.null(split_val) && "split" %in% names(sub))
    sub <- sub[sub$split == split_val, ]
  if (!is.null(scenario_val) && "scenario" %in% names(sub))
    sub <- sub[sub$scenario == scenario_val, ]
  if (!is.null(frac_val) && "missing_frac" %in% names(sub))
    sub <- sub[sub$missing_frac == frac_val, ]
  vals <- sub[[metric]]
  vals <- vals[is.finite(vals)]
  if (!length(vals)) NA_real_ else mean(vals)
}

best_acc <- function(df, method_val, scenario_val = NULL, frac_val = NULL,
                     split_val = "test") {
  best_rmse(df, method_val, scenario_val, frac_val, split_val, "accuracy")
}

# ---------------------------------------------------------------------------
# Extract headline numbers from each benchmark
# ---------------------------------------------------------------------------

extract_pertype <- function(name, baseline_method = "BM",
                             metric = "rmse", higher_is_better = FALSE,
                             frac = 0.25) {
  r <- load_rds(paste0("bench_", name))
  if (is.null(r) || is.null(r$results)) return(NULL)
  test_df <- r$results[r$results$split == "test", ]

  # Pick the "primary" scenario (first in scenarios_primary, or highest signal)
  scenarios <- r$scenarios_primary %||% sort(unique(test_df$scenario))
  scen <- scenarios[length(scenarios)]  # usually highest signal / most informative

  b <- best_rmse(test_df, baseline_method, scen, frac, metric = metric)
  p <- best_rmse(test_df, "pigauto",       scen, frac, metric = metric)

  if (higher_is_better) {
    improvement <- if (all(is.finite(c(b, p)))) p - b else NA_real_
  } else {
    improvement <- if (all(is.finite(c(b, p)))) lift_pct(b, p) else NA_real_
  }
  list(baseline = b, pigauto = p, improvement = improvement,
       nrow = nrow(r$results), epochs = r$epochs %||% NA, scen = scen)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Git commit
# ---------------------------------------------------------------------------

commit_short <- tryCatch(
  substr(system("git -C \".\" rev-parse HEAD 2>/dev/null", intern = TRUE)[1], 1, 10),
  error = function(e) "unknown"
)

# ---------------------------------------------------------------------------
# Extract metrics for all 14 benchmarks
# ---------------------------------------------------------------------------

cat("Loading benchmark results...\n")

b_binary      <- extract_pertype("binary",   "baseline", "accuracy", TRUE,  0.25)
b_categorical <- extract_pertype("categorical", "baseline", "accuracy", TRUE, 0.25)
b_continuous  <- extract_pertype("continuous", "BM", "rmse", FALSE, 0.25)
b_count       <- extract_pertype("count",    "baseline", "rmse",     FALSE, 0.25)
b_ordinal     <- extract_pertype("ordinal",  "baseline", "rmse",     FALSE, 0.25)
b_proportion  <- extract_pertype("proportion", "baseline","rmse",    FALSE, 0.25)
b_zi_count    <- extract_pertype("zi_count", "baseline", "rmse",     FALSE, 0.25)
b_mechanism         <- extract_pertype("missingness_mechanism", "baseline", "rmse", FALSE, 0.25)
b_multi_proportion  <- extract_pertype("multi_proportion", "baseline", "aitchison", FALSE, 0.25)

# Covariate simulation: pigauto_covs vs pigauto at strong_env scenario
b_covsim <- local({
  r <- load_rds("bench_covariate_sim")
  if (is.null(r) || is.null(r$results)) return(NULL)
  res <- r$results
  # Use the scenario with strongest env signal for main headline
  scen <- "mod_phylo_strong_env"
  frac <- 0.40
  no_cov <- mean(res$rmse[res$method == "pigauto" &
                           res$scenario == scen &
                           abs(res$missing_frac - frac) < 0.01], na.rm = TRUE)
  with_cov <- mean(res$rmse[res$method == "pigauto_covs" &
                              res$scenario == scen &
                              abs(res$missing_frac - frac) < 0.01], na.rm = TRUE)
  list(no_cov = no_cov, with_cov = with_cov,
       improvement = lift_pct(no_cov, with_cov),
       nrow = nrow(res))
})

# Tree uncertainty: report FMI range across coefficients at 20% missingness
b_tree_unc <- local({
  r <- load_rds("bench_tree_uncertainty")
  if (is.null(r) || is.null(r$results)) return(NULL)
  res <- r$results
  sub <- res[abs(res$missing_frac - 0.20) < 0.01 & res$method == "single_tree", ]
  fmi_vals <- sub$fmi[is.finite(sub$fmi)]
  list(fmi_mean = if (length(fmi_vals)) mean(fmi_vals) else NA,
       fmi_max  = if (length(fmi_vals)) max(fmi_vals)  else NA,
       nrow = nrow(res))
})

# Multi-obs: obs_rmse for pigauto_cov vs pigauto_no_cov
b_multiobs <- local({
  r <- load_rds("bench_multi_obs")
  if (is.null(r) || is.null(r$results)) return(NULL)
  res <- r$results
  no_cov   <- mean(res$obs_rmse[res$method == "pigauto_no_cov"], na.rm = TRUE)
  with_cov <- mean(res$obs_rmse[res$method == "pigauto_cov"],    na.rm = TRUE)
  list(no_cov = no_cov, with_cov = with_cov,
       improvement = lift_pct(no_cov, with_cov),
       nrow = nrow(res))
})

# Delhey: RMSE for BM_only vs pigauto vs pigauto_covs at 40% missingness
b_delhey <- local({
  r <- load_rds("bench_delhey")
  if (is.null(r) || is.null(r$results)) return(NULL)
  res <- r$results
  test_df <- if ("split" %in% names(res)) res[res$split == "test", ] else res
  frac <- 0.40
  bm  <- mean(test_df$rmse[test_df$method == "BM_only" &
                              abs(test_df$missing_frac - frac) < 0.01], na.rm = TRUE)
  pig <- mean(test_df$rmse[test_df$method == "pigauto" &
                              abs(test_df$missing_frac - frac) < 0.01], na.rm = TRUE)
  cov <- mean(test_df$rmse[test_df$method == "pigauto_covs" &
                              abs(test_df$missing_frac - frac) < 0.01], na.rm = TRUE)
  list(bm = bm, pigauto = pig, pigauto_covs = cov,
       improvement_bm = lift_pct(bm, pig),
       improvement_cov = lift_pct(pig, cov),
       nrow = nrow(res))
})

# Avonet missingness: RMSE at 50% missingness for BM vs pigauto
b_avonet <- local({
  r <- load_rds("bench_avonet_missingness")
  if (is.null(r) || is.null(r$results)) return(NULL)
  res <- r$results
  test_df <- if ("split" %in% names(res)) res[res$split == "test", ] else res
  frac <- 0.50
  bm  <- mean(test_df$rmse[test_df$method == "BM" &
                              abs(test_df$missing_frac - frac) < 0.01], na.rm = TRUE)
  pig <- mean(test_df$rmse[test_df$method == "pigauto" &
                              abs(test_df$missing_frac - frac) < 0.01], na.rm = TRUE)
  list(bm = bm, pigauto = pig, n_species = r$n_species,
       improvement = lift_pct(bm, pig),
       nrow = nrow(res))
})

# Scaling: max N with status "ok" across all stages
b_scaling <- local({
  r <- load_rds("bench_scaling_v090")
  if (is.null(r)) return(NULL)
  df <- if (is.data.frame(r)) r else r$results
  if (is.null(df)) return(NULL)
  ok_rows <- df[df$status == "ok", ]
  list(max_n = if (nrow(ok_rows)) max(ok_rows$n) else NA,
       n_rows = nrow(df),
       n_ok   = nrow(ok_rows))
})

# ---------------------------------------------------------------------------
# HTML generation helpers
# ---------------------------------------------------------------------------

lines <- character(0)
h <- function(...) lines <<- c(lines, paste0(...))

# Create one row for the summary table
#   status: "ok" | "pending" | "error"
#   benchmark_file: relative link to bench HTML
make_row <- function(name, trait_type, dataset, method_cmp,
                     baseline_val, pigauto_val, improvement,
                     metric_label, bench_file, status = "ok") {
  if (status == "pending") {
    h('<tr>',
      '<td>', name, '</td>',
      '<td><em>', trait_type, '</em></td>',
      '<td>', dataset, '</td>',
      '<td colspan="4" style="color:#9ca3af;font-style:italic">',
      'pending &mdash; script not yet run</td>',
      '<td>&mdash;</td>',
      '</tr>')
    return(invisible(NULL))
  }

  impr_str <- if (!is.null(improvement) && is.finite(improvement)) {
    col <- if (improvement > 0) "#059669" else "#dc2626"
    arrow <- if (improvement > 0) "&#9650;" else "&#9660;"
    sprintf('<span style="color:%s">%s %.1f%%</span>', col, arrow, abs(improvement * 100))
  } else "&ndash;"

  h('<tr>',
    '<td><a href="', bench_file, '">', name, '</a></td>',
    '<td><em>', trait_type, '</em></td>',
    '<td>', dataset, '</td>',
    '<td>', method_cmp, '</td>',
    '<td>', fmt(baseline_val), '</td>',
    '<td>', fmt(pigauto_val), '</td>',
    '<td>', impr_str, '</td>',
    '<td>', metric_label, '</td>',
    '</tr>')
}

# ---------------------------------------------------------------------------
# Build HTML
# ---------------------------------------------------------------------------

h('<!doctype html>')
h('<html lang="en">')
h('<head>')
h('<meta charset="utf-8">')
h('<title>pigauto v0.9.0 &mdash; Validation Suite</title>')
h('<style>')
h('  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;')
h('         max-width: 1040px; margin: 2em auto; padding: 0 1.5em; color: #111827;')
h('         line-height: 1.55; }')
h('  h1 { border-bottom: 2px solid #111827; padding-bottom: .25em; }')
h('  h2 { margin-top: 2em; color: #1f2937; border-bottom: 1px solid #d1d5db;')
h('       padding-bottom: .15em; }')
h('  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 13px; }')
h('  th, td { padding: 6px 10px; border-bottom: 1px solid #e5e7eb; text-align: left; vertical-align: top; }')
h('  th { background: #f3f4f6; font-weight: 600; }')
h('  .meta { color: #6b7280; font-size: 13px; }')
h('  .note { background: #fef3c7; border-left: 4px solid #d97706;')
h('         padding: 1em 1.2em; border-radius: 4px; margin: 1.5em 0; }')
h('  .note b { color: #92400e; }')
h('  .badge { display: inline-block; background: #eef2ff; color: #3730a3;')
h('           padding: 2px 10px; border-radius: 999px; font-size: 13px;')
h('           font-weight: 600; margin-left: 8px; vertical-align: middle; }')
h('  code { background: #f3f4f6; padding: 1px 5px; border-radius: 3px; font-size: 12px; }')
h('  a { color: #2563eb; }')
h('  footer { color: #6b7280; font-size: 12px; margin-top: 3em;')
h('           border-top: 1px solid #e5e7eb; padding-top: 1em; }')
h('</style>')
h('</head>')
h('<body>')
h()
h('<h1>pigauto <span class="badge">v0.9.0</span> &mdash; Validation Suite</h1>')
h('<p class="meta">Generated ', timestamp, ' &middot; Commit ', commit_short, '</p>')
h()
h('<p>')
h('This page summarises 15 benchmark experiments validating pigauto v0.9.0 across all')
h('supported trait types, multiple real datasets, and a range of missingness scenarios.')
h('Each row links to a full benchmark report with per-scenario tables, bar charts, and')
h('methodology notes. Benchmarks were run with <code>devtools::load_all()</code> against')
h('the current working tree.')
h('</p>')
h()
h('<h2>Per-trait-type benchmarks</h2>')
h('<p class="meta">Simulated data: <code>ape::rtree(300)</code>, 25% MCAR, 5 replicates,')
h('500 epochs. Improvement = (baseline &minus; pigauto) / baseline &times; 100%.')
h('Positive values (green) indicate pigauto reduces error.</p>')
h()
h('<table>')
h('<thead><tr>')
h('<th>Benchmark</th><th>Type</th><th>Dataset</th><th>Baseline vs</th>')
h('<th>Baseline</th><th>pigauto</th><th>Improvement</th><th>Metric</th>')
h('</tr></thead>')
h('<tbody>')

# Binary
if (!is.null(b_binary)) {
  make_row("Binary traits", "binary", "Simulated (n=300)",
           "Phylo LP", b_binary$baseline, b_binary$pigauto, b_binary$improvement,
           "Accuracy (test)", "dev/bench_binary.html")
} else {
  make_row("Binary traits", "binary", "Simulated", "", NA, NA, NA, "Accuracy", "dev/bench_binary.html", "pending")
}

# Categorical
if (!is.null(b_categorical)) {
  make_row("Categorical traits", "categorical", "Simulated (n=300)",
           "Phylo LP", b_categorical$baseline, b_categorical$pigauto, b_categorical$improvement,
           "Accuracy (test)", "dev/bench_categorical.html")
} else {
  make_row("Categorical traits", "categorical", "Simulated", "", NA, NA, NA, "Accuracy", "dev/bench_categorical.html", "pending")
}

# Continuous
if (!is.null(b_continuous)) {
  make_row("Continuous traits", "continuous", "Simulated (n=300)",
           "BM baseline", b_continuous$baseline, b_continuous$pigauto, b_continuous$improvement,
           "RMSE (test, z-score)", "dev/bench_continuous.html")
} else {
  make_row("Continuous traits", "continuous", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_continuous.html", "pending")
}

# Count
if (!is.null(b_count)) {
  make_row("Count traits", "count", "Simulated (n=300)",
           "BM baseline", b_count$baseline, b_count$pigauto, b_count$improvement,
           "RMSE (test, log1p-z)", "dev/bench_count.html")
} else {
  make_row("Count traits", "count", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_count.html", "pending")
}

# Ordinal
if (!is.null(b_ordinal)) {
  make_row("Ordinal traits", "ordinal", "Simulated (n=300)",
           "BM baseline", b_ordinal$baseline, b_ordinal$pigauto, b_ordinal$improvement,
           "RMSE (test, z-score)", "dev/bench_ordinal.html")
} else {
  make_row("Ordinal traits", "ordinal", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_ordinal.html", "pending")
}

# Proportion
if (!is.null(b_proportion)) {
  make_row("Proportion traits", "proportion", "Simulated (n=300)",
           "BM baseline", b_proportion$baseline, b_proportion$pigauto, b_proportion$improvement,
           "RMSE (test, logit-z)", "dev/bench_proportion.html")
} else {
  make_row("Proportion traits", "proportion", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_proportion.html", "pending")
}

# Multi-proportion
if (!is.null(b_multi_proportion)) {
  make_row("Multi-proportion traits", "multi_proportion", "Simulated (n=300)",
           "BM baseline (CLR)", b_multi_proportion$baseline, b_multi_proportion$pigauto,
           b_multi_proportion$improvement,
           "Aitchison distance (test)", "dev/bench_multi_proportion.html")
} else {
  make_row("Multi-proportion traits", "multi_proportion", "Simulated", "", NA, NA, NA,
           "Aitchison distance", "dev/bench_multi_proportion.html", "pending")
}

# ZI count
if (!is.null(b_zi_count)) {
  make_row("ZI count traits", "zi_count", "Simulated (n=300)",
           "BM/LP baseline", b_zi_count$baseline, b_zi_count$pigauto, b_zi_count$improvement,
           "RMSE magnitude (test)", "dev/bench_zi_count.html")
} else {
  make_row("ZI count traits", "zi_count", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_zi_count.html", "pending")
}

# Missingness mechanism
if (!is.null(b_mechanism)) {
  make_row("Missingness mechanisms", "mixed", "Simulated (n=300)",
           "BM baseline", b_mechanism$baseline, b_mechanism$pigauto, b_mechanism$improvement,
           "RMSE (test, MCAR)", "dev/bench_missingness_mechanism.html")
} else {
  make_row("Missingness mechanisms", "mixed", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_missingness_mechanism.html", "pending")
}

h('</tbody></table>')

# ---------------------------------------------------------------------------
h('<h2>Real-data and feature benchmarks</h2>')
h('<p class="meta">Uses bundled datasets. See individual reports for full details.</p>')
h()
h('<table>')
h('<thead><tr>')
h('<th>Benchmark</th><th>Type</th><th>Dataset</th><th>Question</th>')
h('<th>Method A</th><th>Method B</th><th>Improvement</th><th>Metric</th>')
h('</tr></thead>')
h('<tbody>')

# Avonet missingness
if (!is.null(b_avonet)) {
  make_row(sprintf("Avonet full (%s sp.)", format(b_avonet$n_species, big.mark=",")),
           "mixed (7 traits)", "AVONET + BirdTree (n=9,993)",
           "BM vs pigauto", b_avonet$bm, b_avonet$pigauto, b_avonet$improvement,
           "RMSE at 50% miss.", "dev/bench_avonet_missingness.html")
} else {
  make_row("Avonet full missingness", "mixed", "AVONET (9,993 sp)", "", NA, NA, NA, "RMSE", "dev/bench_avonet_missingness.html", "pending")
}

# Delhey covariate study
if (!is.null(b_delhey)) {
  make_row("Delhey covariate study", "continuous", "Delhey 2019 (n=5,809)",
           "BM vs pigauto", b_delhey$bm, b_delhey$pigauto, b_delhey$improvement_bm,
           "RMSE at 40% miss.", "dev/bench_delhey.html")
} else {
  make_row("Delhey covariate study", "continuous", "Delhey 2019 (5,809 sp)", "", NA, NA, NA, "RMSE", "dev/bench_delhey.html", "pending")
}

# Covariate simulation
if (!is.null(b_covsim)) {
  make_row("Covariate benefit (sim)", "continuous", "Simulated (n=300)",
           "pigauto vs +covariates", b_covsim$no_cov, b_covsim$with_cov, b_covsim$improvement,
           "RMSE (mod_env, 40% miss)", "dev/bench_covariate_sim.html")
} else {
  make_row("Covariate benefit (sim)", "continuous", "Simulated", "", NA, NA, NA, "RMSE", "dev/bench_covariate_sim.html", "pending")
}

# Multi-obs
if (!is.null(b_multiobs)) {
  make_row("Multi-obs + covariates", "continuous", "Simulated CTmax (multi-obs)",
           "no-cov vs cov", b_multiobs$no_cov, b_multiobs$with_cov, b_multiobs$improvement,
           "Obs-level RMSE (mean)", "dev/bench_multi_obs.html")
} else {
  make_row("Multi-obs + covariates", "continuous", "Simulated CTmax", "", NA, NA, NA, "Obs RMSE", "dev/bench_multi_obs.html", "pending")
}

# Tree uncertainty
if (!is.null(b_tree_unc)) {
  fmi_str <- if (is.finite(b_tree_unc$fmi_mean)) fmt(b_tree_unc$fmi_mean) else "&ndash;"
  h('<tr>',
    '<td><a href="dev/bench_tree_uncertainty.html">Tree uncertainty</a></td>',
    '<td><em>continuous</em></td>',
    '<td>AVONET300 + BirdTree (50 trees)</td>',
    '<td>single-tree vs multi-tree MI</td>',
    '<td colspan="2">Mean FMI = ', fmi_str, '</td>',
    '<td>&mdash;</td>',
    '<td>FMI (pooled coefs)</td>',
    '</tr>')
} else {
  make_row("Tree uncertainty", "continuous", "AVONET300 + 50 trees", "", NA, NA, NA, "FMI", "dev/bench_tree_uncertainty.html", "pending")
}

# Scaling
if (!is.null(b_scaling)) {
  max_n_str <- if (!is.na(b_scaling$max_n)) format(b_scaling$max_n, big.mark=",") else "&ndash;"
  h('<tr>',
    '<td><a href="dev/bench_scaling_v090.html">Scaling (v0.9.0)</a></td>',
    '<td><em>mixed</em></td>',
    '<td>Simulated (up to ', max_n_str, ' species)</td>',
    '<td>wall time per N</td>',
    '<td colspan="2">Max N = ', max_n_str, ' (', b_scaling$n_ok, '/', b_scaling$n_rows, ' stages ok)</td>',
    '<td>&mdash;</td>',
    '<td>Wall time (s)</td>',
    '</tr>')
} else {
  make_row("Scaling (v0.9.0)", "mixed", "Simulated (up to 10k sp)", "", NA, NA, NA, "Wall time", "dev/bench_scaling_v090.html", "pending")
}

h('</tbody></table>')

# ---------------------------------------------------------------------------
h()
h('<div class="note">')
h('<b>How to read the Improvement column.</b> For RMSE metrics, positive values (')
h('<span style="color:#059669">&#9650;</span>) mean pigauto achieved lower error than the baseline.')
h('For accuracy metrics, positive values mean pigauto achieved higher accuracy.')
h('A value near zero is expected when phylogenetic signal dominates and the calibrated gate')
h('closes to zero &mdash; in that case pigauto falls back to the baseline by design.')
h('</div>')

h()
h('<h2>Reproducibility</h2>')
h('<ul>')
h('<li>Package version: <code>pigauto 0.9.0</code></li>')
h('<li>Commit: <code>', commit_short, '</code></li>')
h('<li>Run on: ', timestamp, '</li>')
h('<li>All scripts in <code>script/bench_*.R</code>; generators in <code>script/make_bench_*_html.R</code></li>')
h('<li>All .rds result files checked into <code>script/</code></li>')
h('</ul>')

h()
h('<footer>')
h('Generator: <code>script/make_validation_suite_html.R</code>')
h('</footer>')
h()
h('</body>')
h('</html>')

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------

writeLines(lines, out_html_1)
cat("Wrote", out_html_1, "\n")
writeLines(lines, out_html_2)
cat("Wrote", out_html_2, "\n")
