#!/usr/bin/env Rscript
#
# script/make_bench_phase8_aggregate_html.R
#
# Phase 8 aggregate report: reads the 5 Phase 8 bench RDS files and emits
# a single landing page linking them plus a TL;DR table.

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
rds_sweep <- file.path(here, "script", "bench_signal_sweep.rds")
rds_head  <- file.path(here, "script", "bench_bace_avonet_head_to_head.rds")
rds_corr  <- file.path(here, "script", "bench_correlation_sweep.rds")
rds_evo   <- file.path(here, "script", "bench_evo_model_sweep.rds")
rds_clade <- file.path(here, "script", "bench_clade_missingness.rds")
out_html_script  <- file.path(here, "script", "phase8_summary.html")
out_html_pkgdown <- file.path(here, "pkgdown", "assets", "dev",
                                "phase8_summary.html")
dir.create(dirname(out_html_pkgdown), showWarnings = FALSE, recursive = TRUE)

# ---- Load RDS files (tolerate missing) --------------------------------

load_or_null <- function(p) if (file.exists(p)) readRDS(p) else NULL
sweep <- load_or_null(rds_sweep)
head2 <- load_or_null(rds_head)
corr  <- load_or_null(rds_corr)
evo   <- load_or_null(rds_evo)
clade <- load_or_null(rds_clade)

# ---- TL;DR summary: mean metric per method across all benches ---------

summary_rows <- list()
append_summary <- function(bench_obj, tag) {
  if (is.null(bench_obj)) return(NULL)
  res <- bench_obj$results
  res <- res[res$metric %in% c("accuracy", "rmse", "pearson_r"), ]
  m <- aggregate(value ~ method + metric, data = res, FUN = mean)
  m$source <- tag
  summary_rows[[length(summary_rows) + 1L]] <<- m
}
append_summary(sweep, "signal_sweep")
append_summary(head2, "bace_head_to_head")
append_summary(corr,  "correlation_sweep")
append_summary(evo,   "evo_model_sweep")
append_summary(clade, "clade_missingness")

summary_tbl <- if (length(summary_rows)) do.call(rbind, summary_rows) else NULL

render_summary <- function(tbl) {
  if (is.null(tbl)) return("<p><em>No bench outputs found yet.</em></p>")
  rows <- apply(tbl, 1L, function(r) {
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td></tr>",
             r[["source"]], r[["method"]], r[["metric"]],
             as.numeric(r[["value"]]))
  })
  paste0(
    "<table><tr>",
    "<th>source</th><th>method</th><th>metric</th><th>mean value</th></tr>",
    paste(rows, collapse = "\n"),
    "</table>"
  )
}

# ---- Sub-report list --------------------------------------------------

sub_reports <- c(
  signal_sweep       = "Pagel's λ signal-strength sweep (MVP)",
  bace_head_to_head  = "AVONET 300 head-to-head vs BACE (MVP)",
  correlation_sweep  = "Cross-trait-correlation (ρ) sweep (Phase 8.1)",
  evo_model_sweep    = "Evolutionary-model sweep — BM / OU / regime_shift / nonlinear (Phase 8.2)",
  clade_missingness  = "Clade-correlated missingness (realistic MAR) (Phase 8.3)"
)
bench_rds_map <- list(
  signal_sweep       = sweep,
  bace_head_to_head  = head2,
  correlation_sweep  = corr,
  evo_model_sweep    = evo,
  clade_missingness  = clade
)
bench_html_map <- list(
  signal_sweep       = "bench_signal_sweep.html",
  bace_head_to_head  = "bench_bace_avonet_head_to_head.html",
  correlation_sweep  = "bench_correlation_sweep.html",
  evo_model_sweep    = "bench_evo_model_sweep.html",
  clade_missingness  = "bench_clade_missingness.html"
)
sub_list <- paste0(
  "<ul>",
  paste0(
    vapply(names(sub_reports), function(k) {
      if (!is.null(bench_rds_map[[k]])) {
        sprintf("<li><a href='%s'>%s</a></li>",
                bench_html_map[[k]], sub_reports[[k]])
      } else {
        sprintf("<li><em>%s — not yet run.</em></li>", sub_reports[[k]])
      }
    }, character(1)),
    collapse = ""
  ),
  "</ul>"
)

# ---- Reproducibility block --------------------------------------------

sess_info <- paste(capture.output(utils::sessionInfo()), collapse = "\n")
pigauto_version <- tryCatch(as.character(utils::packageVersion("pigauto")),
                              error = function(e) "?")
bace_version <- tryCatch(as.character(utils::packageVersion("BACE")),
                           error = function(e) "not installed")
torch_version <- tryCatch(as.character(utils::packageVersion("torch")),
                            error = function(e) "?")

# ---- Render -----------------------------------------------------------

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Phase 8 benchmark summary</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222} h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px} pre{background:#f6f6f6;padding:1em;overflow:auto}</style>",
  "</head><body>",
  "<h1>Phase 8: discriminative benchmark suite</h1>",
  "<p>pigauto v", pigauto_version,
  " · BACE ", bace_version,
  " · torch ", torch_version, "</p>",
  "<h2>TL;DR — mean metric per method across all benches</h2>",
  render_summary(summary_tbl),
  "<p><em>Lower is better for RMSE; higher is better for accuracy and Pearson r.</em></p>",
  "<h2>Sub-reports</h2>", sub_list,
  "<h2>Reproducibility</h2><pre>", sess_info, "</pre>",
  "</body></html>"
)

writeLines(html, out_html_script)
writeLines(html, out_html_pkgdown)
cat("Wrote", out_html_script, "\n")
cat("Wrote", out_html_pkgdown, "\n")
