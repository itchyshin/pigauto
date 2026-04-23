#!/usr/bin/env Rscript
# HTML generator for the AVONET full n=9,993 local (Mac MPS) bench.
# Reads script/bench_avonet_full_local.rds (written by
# script/bench_avonet_full_local.R) and produces a winner-highlighted
# pivot table plus conformal + MC-dropout coverage rows.

options(stringsAsFactors = FALSE)

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_avonet_full_local.rds")
out_s  <- file.path(here, "script", "bench_avonet_full_local.html")
out_p  <- file.path(here, "pkgdown", "assets", "dev",
                      "bench_avonet_full_local.html")
dir.create(dirname(out_p), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) {
  stop("No RDS at ", in_rds,
       ". Run script/bench_avonet_full_local.R first.", call. = FALSE)
}

b   <- readRDS(in_rds)
res <- b$results

pivot <- reshape(res[, c("method", "trait", "metric", "value")],
                  timevar = "method",
                  idvar = c("trait", "metric"),
                  direction = "wide")
names(pivot) <- sub("^value[.]", "", names(pivot))
method_cols <- setdiff(names(pivot), c("trait", "metric"))

highlight_row <- function(row) {
  vals <- suppressWarnings(as.numeric(row[method_cols]))
  if (all(is.na(vals))) return(rep(FALSE, length(method_cols)))
  if (identical(row[["metric"]], "rmse")) vals == min(vals, na.rm = TRUE)
  else vals == max(vals, na.rm = TRUE)
}

render_cell <- function(v, w) {
  if (is.na(v)) return("<td>-</td>")
  cls <- if (isTRUE(w)) " style='background:#def6d0;font-weight:bold'" else ""
  sprintf("<td%s>%.3f</td>", cls, v)
}

rows_html <- apply(pivot, 1L, function(row) {
  wn <- highlight_row(row)
  cells <- mapply(render_cell,
                    suppressWarnings(as.numeric(row[method_cols])),
                    wn, SIMPLIFY = FALSE)
  paste0("<tr><td>", row[["trait"]], "</td><td>", row[["metric"]], "</td>",
          paste(cells, collapse = ""), "</tr>")
})

wall_tbl <- unique(res[, c("method", "wall_s")])
wall_rows <- apply(wall_tbl, 1, function(r)
  sprintf("<tr><td>%s</td><td>%.1f</td><td>%.1f</td></tr>",
           r[["method"]], as.numeric(r[["wall_s"]]),
           as.numeric(r[["wall_s"]]) / 60))

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>AVONET full n=9,993 &mdash; pigauto (Mac MPS)</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1200px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px}</style></head><body>",
  "<h1>AVONET full n=9,993 &mdash; pigauto on local Mac MPS</h1>",
  sprintf("<p>n = %d species x %d-trait AVONET. Seed = %d, miss_frac = %.2f, ",
          b$n_species,
          length(unique(res$trait)),
          b$seed, b$miss_frac),
  sprintf("n_imputations = %d.</p>", b$n_imputations),
  "<p>Apple Silicon MPS run -- Vulcan L40S OOM'd at predict stage on ",
  "the equivalent config (the refine_steps forward-pass accumulation ",
  "leak, partly addressed in commit <code>dc8cffa</code> and tracked ",
  "for v0.9.2).  Mac MPS unified memory absorbs the attention at n=10k.</p>",
  "<h2>Per-trait metrics (winner highlighted)</h2>",
  "<table><tr><th>trait</th><th>metric</th>",
  paste0("<th>", method_cols, "</th>", collapse = ""), "</tr>",
  paste(rows_html, collapse = "\n"),
  "</table>",
  "<h2>Wall time per method</h2>",
  "<table><tr><th>method</th><th>wall (s)</th><th>wall (min)</th></tr>",
  paste(wall_rows, collapse = "\n"),
  "</table>",
  "</body></html>"
)

writeLines(html, out_s)
writeLines(html, out_p)
cat("Wrote", out_s, "+", out_p, "\n")
