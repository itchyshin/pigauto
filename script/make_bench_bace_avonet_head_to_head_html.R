#!/usr/bin/env Rscript
#
# Render the pigauto-vs-BACE head-to-head RDS into HTML.

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_bace_avonet_head_to_head.rds")
out_html_script  <- file.path(here, "script", "bench_bace_avonet_head_to_head.html")
out_html_pkgdown <- file.path(here, "pkgdown", "assets", "dev",
                                "bench_bace_avonet_head_to_head.html")
dir.create(dirname(out_html_pkgdown), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) {
  stop("Run script/bench_bace_avonet_head_to_head.R first.")
}

bench    <- readRDS(in_rds)
results  <- bench$results
bace_ran <- isTRUE(bench$bace_ran)

# Pivot: one row per (trait, metric), one column per method
pivot <- reshape(results[, c("method", "trait", "metric", "value")],
                  timevar = "method", idvar = c("trait", "metric"),
                  direction = "wide")
names(pivot) <- sub("^value\\.", "", names(pivot))

# Winner highlighting: for 'rmse', lower is better; others higher is better.
method_cols <- setdiff(names(pivot), c("trait", "metric"))
highlight_row <- function(row) {
  vals <- as.numeric(row[method_cols])
  if (all(is.na(vals))) return(rep(FALSE, length(method_cols)))
  if (identical(row[["metric"]], "rmse")) {
    vals == min(vals, na.rm = TRUE)
  } else {
    vals == max(vals, na.rm = TRUE)
  }
}

render_cell <- function(v, is_winner) {
  if (is.na(v)) return("<td>-</td>")
  cls <- if (is_winner) " style='background:#def6d0;font-weight:bold'" else ""
  sprintf("<td%s>%.3f</td>", cls, v)
}

rows_html <- apply(pivot, 1L, function(row) {
  wn <- highlight_row(row)
  cells <- mapply(
    render_cell,
    as.numeric(row[method_cols]),
    wn,
    SIMPLIFY = FALSE
  )
  paste0("<tr><td>", row[["trait"]], "</td><td>", row[["metric"]], "</td>",
          paste(cells, collapse = ""), "</tr>")
})

wall_by_method <- aggregate(wall_s ~ method, data = results,
                             FUN = function(x) unique(x)[1])

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Phase 8: AVONET head-to-head</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 8px}",
  "</style></head><body>",
  "<h1>Phase 8 MVP: AVONET 300 — pigauto vs BACE</h1>",
  sprintf("<p>Single-seed comparison on avonet300/tree300 (n = 300). Seed = %d, miss_frac = %.2f.</p>",
          bench$seed, bench$miss_frac),
  if (!bace_ran)
    "<p><strong>BACE did not run</strong> (not installed, or output shape unrecognised). Comparison below is pigauto paths only.</p>"
  else "",
  "<h2>Per-trait metrics (winner highlighted)</h2>",
  "<table><tr><th>trait</th><th>metric</th>",
  paste0("<th>", method_cols, "</th>", collapse = ""), "</tr>",
  paste(rows_html, collapse = "\n"),
  "</table>",
  "<h2>Wall time per method</h2>",
  "<table><tr><th>method</th><th>wall (s)</th></tr>",
  paste0(
    apply(wall_by_method, 1L, function(r) {
      sprintf("<tr><td>%s</td><td>%.1f</td></tr>",
               r[["method"]], as.numeric(r[["wall_s"]]))
    }),
    collapse = "\n"
  ),
  "</table>",
  "</body></html>"
)

writeLines(html, out_html_script)
writeLines(html, out_html_pkgdown)
cat("Wrote", out_html_script, "\n")
cat("Wrote", out_html_pkgdown, "\n")
