#!/usr/bin/env Rscript
#
# script/make_bench_phase8_aggregate_html.R
#
# Phase 8 MVP aggregate report: reads both bench RDS files and emits a
# single landing page linking them plus a TL;DR table.

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
rds_sweep  <- file.path(here, "script", "bench_signal_sweep.rds")
rds_head   <- file.path(here, "script", "bench_bace_avonet_head_to_head.rds")
out_html_script  <- file.path(here, "script", "phase8_summary.html")
out_html_pkgdown <- file.path(here, "pkgdown", "assets", "dev",
                                "phase8_summary.html")
dir.create(dirname(out_html_pkgdown), showWarnings = FALSE, recursive = TRUE)

# ---- Load RDS files (tolerate missing) --------------------------------

sweep <- if (file.exists(rds_sweep)) readRDS(rds_sweep) else NULL
head2 <- if (file.exists(rds_head))  readRDS(rds_head)  else NULL

# ---- TL;DR summary: accuracy / RMSE means per method across all tests --

summary_rows <- list()

if (!is.null(sweep)) {
  res <- sweep$results
  sw_mean <- aggregate(value ~ method + metric,
                        data = res[res$metric %in% c("accuracy", "rmse", "pearson_r"), ],
                        FUN = mean)
  sw_mean$source <- "signal_sweep"
  summary_rows[[length(summary_rows) + 1L]] <- sw_mean
}

if (!is.null(head2)) {
  res2 <- head2$results
  hh_mean <- aggregate(value ~ method + metric,
                        data = res2[res2$metric %in% c("accuracy", "rmse", "pearson_r"), ],
                        FUN = mean)
  hh_mean$source <- "bace_head_to_head"
  summary_rows[[length(summary_rows) + 1L]] <- hh_mean
}

summary_tbl <- if (length(summary_rows)) do.call(rbind, summary_rows) else NULL

render_summary <- function(tbl) {
  if (is.null(tbl)) return("<p><em>No bench outputs found yet.</em></p>")
  rows <- apply(tbl, 1L, function(r) {
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td></tr>",
             r[["source"]], r[["method"]], r[["metric"]],
             as.numeric(r[["value"]]))
  })
  paste0(
    "<table style='border-collapse:collapse'><tr>",
    "<th>source</th><th>method</th><th>metric</th><th>mean value</th></tr>",
    paste(rows, collapse = "\n"),
    "</table>"
  )
}

# ---- Reproducibility block --------------------------------------------

sess_info <- paste(capture.output(utils::sessionInfo()), collapse = "\n")
pigauto_version <- tryCatch(as.character(utils::packageVersion("pigauto")),
                              error = function(e) "?")
bace_version <- tryCatch(as.character(utils::packageVersion("BACE")),
                           error = function(e) "not installed")
torch_version <- tryCatch(as.character(utils::packageVersion("torch")),
                            error = function(e) "?")

# ---- Render page -------------------------------------------------------

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

  "<h2>TL;DR — mean metric per method across both benches</h2>",
  render_summary(summary_tbl),
  "<p><em>Lower is better for RMSE; higher is better for accuracy and Pearson r.</em></p>",

  "<h2>Sub-reports</h2>",
  "<ul>",
  if (!is.null(sweep))
    "<li><a href='bench_signal_sweep.html'>Pagel's λ signal-strength sweep</a></li>"
  else
    "<li><em>Signal sweep not yet run.</em></li>",
  if (!is.null(head2))
    "<li><a href='bench_bace_avonet_head_to_head.html'>AVONET 300 head-to-head vs BACE</a></li>"
  else
    "<li><em>AVONET head-to-head not yet run.</em></li>",
  "</ul>",

  "<h2>Reproducibility</h2>",
  "<pre>", sess_info, "</pre>",

  "<p><em>Follow-ups deferred to Phase 8.x: ρ sweep, evolutionary-model sweep, clade-correlated missingness.</em></p>",
  "</body></html>"
)

writeLines(html, out_html_script)
writeLines(html, out_html_pkgdown)
cat("Wrote", out_html_script, "\n")
cat("Wrote", out_html_pkgdown, "\n")
