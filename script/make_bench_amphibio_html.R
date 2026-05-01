#!/usr/bin/env Rscript
# HTML generator for the AmphiBIO amphibian bench.  Reads
# script/bench_amphibio.rds, writes a winner-highlighted pivot table.
options(stringsAsFactors = FALSE)

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds  <- file.path(here, "script", "bench_amphibio.rds")
out_s   <- file.path(here, "script", "bench_amphibio.html")
out_p   <- file.path(here, "pkgdown", "assets", "dev", "bench_amphibio.html")
dir.create(dirname(out_p), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) stop("No RDS at ", in_rds)

b   <- readRDS(in_rds)
res <- b$results

pivot <- reshape(res[, c("method", "trait", "metric", "value")],
                  timevar = "method",
                  idvar = c("trait", "metric"),
                  direction = "wide")
names(pivot) <- sub("^value[.]", "", names(pivot))
method_cols <- setdiff(names(pivot), c("trait", "metric"))

render_cell <- function(v, is_winner) {
  if (is.na(v)) return("<td>-</td>")
  cls <- if (isTRUE(is_winner))
    " style='background:#def6d0;font-weight:bold'" else ""
  sprintf("<td%s>%.3f</td>", cls, v)
}

highlight_row <- function(row) {
  vals <- suppressWarnings(as.numeric(row[method_cols]))
  if (all(is.na(vals))) return(rep(FALSE, length(method_cols)))
  if (identical(row[["metric"]], "rmse")) vals == min(vals, na.rm = TRUE)
  else vals == max(vals, na.rm = TRUE)
}

rows_html <- apply(pivot, 1L, function(row) {
  wn <- highlight_row(row)
  cells <- mapply(render_cell,
                    suppressWarnings(as.numeric(row[method_cols])),
                    wn, SIMPLIFY = FALSE)
  paste0("<tr><td>", row[["trait"]], "</td><td>", row[["metric"]], "</td>",
          paste(cells, collapse = ""), "</tr>")
})

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>AmphiBIO x pigauto</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1040px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px}</style></head><body>",
  "<h1>AmphiBIO + taxonomic tree &mdash; pigauto</h1>",
  sprintf("<p>n = %d amphibian species, seed = %d, miss_frac = %.2f, n_imputations = %d.</p>",
          b$n_species, b$seed, b$miss_frac, b$n_imputations),
  "<p>Completes the tetrapod story (birds, mammals, fish, amphibians). ",
  "Continuous-only v1; binary/categorical AmphiBIO columns currently ",
  "hit a pigauto threshold-joint baseline type issue, parked for v0.9.2.</p>",
  "<h2>Per-trait metrics (winner highlighted)</h2>",
  "<table><tr><th>trait</th><th>metric</th>",
  paste0("<th>", method_cols, "</th>", collapse = ""), "</tr>",
  paste(rows_html, collapse = "\n"),
  "</table>",
  "</body></html>"
)

writeLines(html, out_s)
writeLines(html, out_p)
cat("Wrote", out_s, "+", out_p, "\n")
