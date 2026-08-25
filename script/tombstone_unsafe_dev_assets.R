#!/usr/bin/env Rscript

targets <- c(
  "bench_avonet9993_bace.html",
  "bench_avonet9993_bace_index.html",
  "bench_avonet9993_bace_n3000.html",
  "bench_bace_avonet_head_to_head.html",
  "bench_bien.html",
  "bench_fishbase.html",
  "bench_pantheria_bace_head_to_head.html",
  "calibration_grid.html",
  "pantheria_summary.html",
  "phase8_summary.html",
  "tests_overview.html"
)
root <- file.path("pkgdown", "assets", "dev")
missing <- targets[!file.exists(file.path(root, targets))]
if (length(missing)) stop("Missing tombstone target(s): ", paste(missing, collapse = ", "), call. = FALSE)

for (name in targets) {
  title <- sub("[.]html$", "", name)
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    sprintf("<title>%s — historical development page withdrawn</title>", title),
    "<meta name=\"robots\" content=\"noindex\"></head><body>",
    "<main><h1>Historical development page withdrawn</h1>",
    sprintf("<p><code>%s</code> is retained only as a stable tombstone for an old link.</p>", name),
    "<p>Its former contents are not current evidence for pigauto 0.11 and must not be read as a performance, external-competitiveness, interval-coverage, or downstream-inference claim.</p>",
    "<p>The external-comparator gate remains open. Use the current package articles and reference pages for supported workflows and claim boundaries.</p>",
    "</main></body></html>"
  )
  writeLines(html, file.path(root, name))
}

cat(sprintf("Tombstoned %d unsafe historical development pages\n", length(targets)))
