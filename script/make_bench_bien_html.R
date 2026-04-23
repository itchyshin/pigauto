#!/usr/bin/env Rscript
# HTML generator for the BIEN + V.PhyloMaker2 plant bench.  Reads
# script/bench_bien.rds (written by script/bench_bien.R) and produces
# a winner-highlighted pivot table plus conformal + MC-dropout
# coverage rows.  This is the "kingdom jump" in the pigauto
# taxonomic-breadth series -- birds, mammals, fish, amphibians (all
# animal kingdom) + plants.

options(stringsAsFactors = FALSE)

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_bien.rds")
out_s  <- file.path(here, "script", "bench_bien.html")
out_p  <- file.path(here, "pkgdown", "assets", "dev", "bench_bien.html")
dir.create(dirname(out_p), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) {
  stop("No RDS at ", in_rds,
       ". Run script/bench_bien.R first.", call. = FALSE)
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
  "<title>BIEN + V.PhyloMaker2 plants &mdash; pigauto</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1200px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px}</style></head><body>",
  "<h1>BIEN + V.PhyloMaker2 plants &mdash; pigauto (kingdom jump)</h1>",
  sprintf("<p>n = %d species x %d traits. Seed = %d, miss_frac = %.2f, ",
          b$n_species,
          length(unique(res$trait)),
          b$seed, b$miss_frac),
  sprintf("n_imputations = %d.</p>", b$n_imputations),
  "<p>Plant traits from BIEN (Botanical Information and Ecology Network) ",
  "via <code>BIEN_trait_trait()</code> aggregated to species means; ",
  "tree from V.PhyloMaker2 on the Smith &amp; Brown 2018 seed-plant ",
  "megaphylogeny backbone (scenario S3 &mdash; random within-genus ",
  "placement of species not in the backbone).</p>",
  "<h2>Interpretation &mdash; a scope-of-phylogenetic-imputation result</h2>",
  "<p>Unlike the four vertebrate classes (birds r=0.85&ndash;0.97, mammals ",
  "r=0.83&ndash;0.98, fish r=0.42&ndash;0.83, amphibians r=0.87&ndash;0.98), ",
  "these BIEN plant traits show Pearson r in the range <b>&ndash;0.02 to 0.43</b> ",
  "on the held-out test cells. Only <b>wood_density</b> (r=0.43) shows ",
  "a clean lift (&ndash;6 % RMSE vs grand mean); on the other four traits ",
  "pigauto is 15&ndash;101 % worse than the grand-mean baseline on RMSE. ",
  "Two compounding causes:</p><ul>",
  "<li><b>Weak phylogenetic signal in pooled BIEN traits.</b> Species-level ",
  "means are aggregated from heterogeneous individual observations (often ",
  "few per species), which dilutes any genuine phylogenetic structure.</li>",
  "<li><b>Random within-genus polytomy resolution.</b> V.PhyloMaker2 ",
  "scenario S3 places species without backbone entries randomly within ",
  "their genus, adding noise to the tree input for tips far from the ",
  "Smith &amp; Brown 2018 ~70-genera backbone.</li></ul>",
  "<p>pigauto's gate safety contains the damage: the calibrated gate closes ",
  "to 0 on weak-signal traits, so the pipeline degrades gracefully toward ",
  "the (also weak) BM prior rather than blowing up. <b>Conformal coverage ",
  "still lands at 0.90&ndash;0.98 on all five traits</b> &mdash; the ",
  "distribution-free uncertainty quantification remains valid even when ",
  "point estimates have low informativeness.</p>",
  "<p>This is the honest paper-ready <i>boundary case</i> for phylogenetic ",
  "trait imputation: across-kingdom generalisation requires either ",
  "(a) a tree with resolved species-level branching (not a backbone + random ",
  "placement), or (b) traits with demonstrably strong phylogenetic signal ",
  "at the scale studied. Wood density passes both; height / leaf area / ",
  "SLA / seed mass do not at this species pool.</p>",
  "<h3>Jensen back-transform note</h3>",
  "<p>An earlier version of this bench ran at <code>n_imputations = 1</code> ",
  "(to avoid a known predict-stage GPU memory leak in v0.9.1). On that ",
  "run, <code>height_m</code> RMSE was 47.7 &mdash; 4&times; the grand-mean ",
  "baseline. The re-run at <code>n_imputations = 20</code> (this page) ",
  "dropped it to 15.15 (3&times; better) by activating the median-pool ",
  "MI correction for log-decoded continuous traits (commit <code>dc8cffa</code>) ",
  "and stripping Jensen back-transform bias. The underlying weak-signal ",
  "problem above is what remains after that fix.</p>",
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
