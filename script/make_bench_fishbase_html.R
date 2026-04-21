#!/usr/bin/env Rscript
# HTML generator for the FishBase + fishtree bench.
# Reads script/bench_fishbase.rds, writes a winner-highlighted pivot
# with per-method RMSE/accuracy/coverage.

options(stringsAsFactors = FALSE)

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_fishbase.rds")
out_s  <- file.path(here, "script", "bench_fishbase.html")
out_p  <- file.path(here, "pkgdown", "assets", "dev", "bench_fishbase.html")
dir.create(dirname(out_p), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) {
  stop("No RDS at ", in_rds, ". Run bench_fishbase.R first.")
}

b   <- readRDS(in_rds)
res <- b$results
bace_ran <- isTRUE(b$bace_ran)

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
  cls <- if (w) " style='background:#def6d0;font-weight:bold'" else ""
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

banner <- if (!is.null(b$n_subset) && !is.na(b$n_subset)) {
  sprintf("<p>Random subset from the full matched intersection of fishtree + FishBase (seed %d, N_SUBSET=%d).</p>",
          b$seed, b$n_subset)
} else ""
banner_bace <- if (!bace_ran) "<p><em>BACE skipped</em> (not installed or PIGAUTO_RUN_BACE != 1).</p>" else ""

# Known-issue caveat for the log-transform + MI-pooling + expm1 decode
# amplification that blows up Length / Weight RMSE at n >= 5000.
# Applies only when n_species is large enough for the outlier
# pattern to dominate (empirically > ~3000).
caveat <- if (!is.null(b$n_species) && b$n_species >= 5000L) {
  paste0(
    "<div style='background:#fef3c7;border-left:4px solid #d97706;",
    "padding:1em 1.2em;border-radius:4px;margin:1em 0'>",
    "<b>Known issue:</b> pigauto's Length / Weight RMSE are inflated on ",
    "this run by an MI-pooling + expm1 decode bug. With ",
    "n_imputations &gt; 1 and log-transformed continuous traits, ",
    "a few dropout-noisy draws with high latent values decode to absurd ",
    "magnitudes, and the mean pool across M imputations is dominated by ",
    "those outliers. Pattern mirrors Issue #40 (fixed for count / ",
    "proportion via PR #41 median pooling) but has not yet been ",
    "extended to log-transformed continuous traits. Calibrated gate ",
    "on Length was 0 on this run &mdash; pigauto falls back to the ",
    "baseline mean anyway, so the decoded Length RMSE here is NOT a ",
    "signal about pigauto&apos;s modelling quality, just about the ",
    "MI-pooling decode path. The BodyShapeI / Troph / Vulnerability / ",
    "DepthRangeDeep lifts are unaffected and show pigauto&apos;s real ",
    "performance at this scale.</div>")
} else ""

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>FishBase + fishtree - pigauto + BACE</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1200px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;",
  "padding:4px 8px}</style></head><body>",
  "<h1>FishBase + fishtree - pigauto vs BACE</h1>",
  sprintf("<p>n = %d species, seed = %d, miss_frac = %.2f. Third real-data ",
          b$n_species, b$seed, b$miss_frac),
  "benchmark completing the vertebrate breadth triad ",
  "(birds via AVONET, mammals via PanTHERIA, fish via FishBase + fishtree).</p>",
  banner, banner_bace, caveat,
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
