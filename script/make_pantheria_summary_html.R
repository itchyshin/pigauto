#!/usr/bin/env Rscript
# PanTHERIA aggregate summary landing page.

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
rds_full  <- file.path(here, "script", "bench_pantheria_full.rds")
rds_bace  <- file.path(here, "script", "bench_pantheria_bace_head_to_head.rds")
out_s <- file.path(here, "script", "pantheria_summary.html")
out_p <- file.path(here, "pkgdown", "assets", "dev", "pantheria_summary.html")
dir.create(dirname(out_p), showWarnings = FALSE, recursive = TRUE)

load_or_null <- function(p) if (file.exists(p)) readRDS(p) else NULL
b_full <- load_or_null(rds_full)
b_bace <- load_or_null(rds_bace)

summary_rows <- list()
append_summary <- function(bench_obj, tag) {
  if (is.null(bench_obj)) return(NULL)
  res <- bench_obj$results
  res <- res[res$metric %in% c("accuracy","rmse","pearson_r"), ]
  m <- aggregate(value ~ method + metric, data = res, FUN = mean)
  m$source <- tag
  summary_rows[[length(summary_rows) + 1L]] <<- m
}
append_summary(b_full, "full_scale")
append_summary(b_bace, "bace_head_to_head")
summary_tbl <- if (length(summary_rows)) do.call(rbind, summary_rows) else NULL

render_summary <- function(tbl) {
  if (is.null(tbl)) return("<p><em>No bench outputs found yet.</em></p>")
  rows <- apply(tbl, 1L, function(r)
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td></tr>",
             r[["source"]], r[["method"]], r[["metric"]],
             as.numeric(r[["value"]])))
  paste0(
    "<table><tr><th>source</th><th>method</th><th>metric</th><th>mean value</th></tr>",
    paste(rows, collapse = "\n"), "</table>")
}

sub_list <- paste0(
  "<ul>",
  if (!is.null(b_full))
    "<li><a href='bench_pantheria_full.html'>Full-scale PanTHERIA (pigauto paths)</a></li>"
  else "<li><em>Full-scale bench — not yet run.</em></li>",
  if (!is.null(b_bace))
    "<li><a href='bench_pantheria_bace_head_to_head.html'>PanTHERIA vs BACE head-to-head (n=500 subset)</a></li>"
  else "<li><em>BACE head-to-head — not yet run.</em></li>",
  "</ul>")

sess <- paste(capture.output(utils::sessionInfo()), collapse = "\n")
pigv <- tryCatch(as.character(utils::packageVersion("pigauto")),
                  error = function(e) "?")
bav  <- tryCatch(as.character(utils::packageVersion("BACE")),
                  error = function(e) "not installed")
tv   <- tryCatch(as.character(utils::packageVersion("torch")),
                  error = function(e) "?")

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>PanTHERIA benchmark summary</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px}pre{background:#f6f6f6;padding:1em;overflow:auto}</style>",
  "</head><body>",
  "<h1>PanTHERIA mammal benchmark — summary</h1>",
  "<p>pigauto v", pigv, " · BACE ", bav, " · torch ", tv, "</p>",
  "<p><em>Second real-data benchmark for pigauto after AVONET (birds). ",
  "PanTHERIA (Jones et al. 2009, <em>Ecology</em>) provides 4,629 mammals × ",
  "~55 traits; this bench uses 8 canonical traits matched against a ",
  "taxonomy-derived cladogram.</em></p>",
  "<h2>TL;DR — mean metric per method across both benches</h2>",
  render_summary(summary_tbl),
  "<p><em>Lower is better for RMSE; higher is better for accuracy and Pearson r.</em></p>",
  "<h2>Sub-reports</h2>", sub_list,
  "<h2>Reproducibility</h2><pre>", sess, "</pre>",
  "</body></html>"
)
writeLines(html, out_s); writeLines(html, out_p)
cat("Wrote", out_s, "+", out_p, "\n")
