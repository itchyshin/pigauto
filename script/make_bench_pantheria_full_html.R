#!/usr/bin/env Rscript
# HTML generator for PanTHERIA full-scale bench.

options(stringsAsFactors = FALSE)

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_pantheria_full.rds")
out_s  <- file.path(here, "script", "bench_pantheria_full.html")
out_p  <- file.path(here, "pkgdown", "assets", "dev",
                      "bench_pantheria_full.html")
dir.create(dirname(out_p), showWarnings = FALSE, recursive = TRUE)
if (!file.exists(in_rds)) stop("Run script/bench_pantheria_full.R first.")

b   <- readRDS(in_rds)
res <- b$results

# Pivot: trait x metric x method
pivot <- reshape(res[, c("method","trait","metric","value")],
                  timevar = "method", idvar = c("trait","metric"),
                  direction = "wide")
names(pivot) <- sub("^value\\.", "", names(pivot))
method_cols <- setdiff(names(pivot), c("trait","metric"))

highlight_row <- function(row) {
  vals <- as.numeric(row[method_cols])
  if (all(is.na(vals))) return(rep(FALSE, length(method_cols)))
  if (identical(row[["metric"]], "rmse")) vals == min(vals, na.rm = TRUE)
  else vals == max(vals, na.rm = TRUE)
}
render_cell <- function(v, w) {
  if (is.na(v)) return("<td>-</td>")
  cls <- if (w) " style='background:#def6d0;font-weight:bold'" else ""
  sprintf("<td%s>%.3f</td>", cls, v)
}
rows_html <- apply(pivot, 1L, function(row) {
  wn <- highlight_row(row)
  cells <- mapply(render_cell, as.numeric(row[method_cols]), wn, SIMPLIFY = FALSE)
  paste0("<tr><td>", row[["trait"]], "</td><td>", row[["metric"]], "</td>",
          paste(cells, collapse = ""), "</tr>")
})
wall_tbl <- unique(res[, c("method","wall_s")])
wall_rows <- apply(wall_tbl, 1, function(r)
  sprintf("<tr><td>%s</td><td>%.1f</td></tr>", r[["method"]],
           as.numeric(r[["wall_s"]])))

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>PanTHERIA full-scale bench</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px}</style></head><body>",
  "<h1>PanTHERIA mammal benchmark (full scale)</h1>",
  sprintf("<p>n = %d mammals (PanTHERIA ∩ taxonomy tree), %d traits, seed = %d, miss_frac = %.2f.</p>",
          b$n_species, b$n_traits, b$seed, b$miss_frac),
  "<p><em>Tree: taxonomy-derived cladogram (Order/Family/Genus/Species via ",
  "<code>ape::as.phylo.formula()</code>) with Grafen branch lengths — not ",
  "a time-calibrated phylogeny. Adequate for testing pigauto's hierarchical ",
  "phylo prior.</em></p>",
  "<h2>Per-trait metrics (winner highlighted)</h2>",
  "<table><tr><th>trait</th><th>metric</th>",
  paste0("<th>", method_cols, "</th>", collapse = ""), "</tr>",
  paste(rows_html, collapse = "\n"),
  "</table>",
  "<h2>Wall time per method</h2>",
  "<table><tr><th>method</th><th>wall (s)</th></tr>",
  paste(wall_rows, collapse = "\n"),
  "</table>",
  "</body></html>"
)
writeLines(html, out_s); writeLines(html, out_p)
cat("Wrote", out_s, "+", out_p, "\n")
