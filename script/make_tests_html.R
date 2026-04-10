#!/usr/bin/env Rscript
# Build a self-contained HTML index of every test_that() block in tests/testthat.
# No package state required; reads test files line by line and extracts
# descriptions via regex. Does NOT run the tests — this is a static index.

test_dir <- "tests/testthat"
if (!dir.exists(test_dir)) {
  stop("Must be run from the pigauto project root (couldn't find ", test_dir, ")")
}

files <- list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE)
files <- sort(files)

# ---- Human-readable area for each file ------------------------------------
area_map <- list(
  "test-preprocess.R"    = "Trait preprocessing, type detection, multi-obs alignment",
  "test-graph.R"         = "Phylogenetic graph, adjacency, spectral coordinates",
  "test-masking.R"       = "Missing-data splits and observed-cell masks",
  "test-fit-predict.R"   = "Training, prediction, attention, conformal, impute()",
  "test-mixed-types.R"   = "5-type detection, encode/decode round-trip",
  "test-multi-impute.R"  = "multi_impute(), with_imputations(), pool_mi() (Rubin's rules)",
  "test-new-features.R"  = "evaluate, summary, plot, report, CV, benchmark, save/load"
)

# ---- Extract test descriptions --------------------------------------------
extract_tests <- function(path) {
  lines <- readLines(path, warn = FALSE)
  # Match test_that("desc", ... or test_that('desc', ...
  m <- regmatches(
    lines,
    regexec("^test_that\\([\"'](.+?)[\"'],", lines)
  )
  hits <- which(vapply(m, length, integer(1)) == 2L)
  if (length(hits) == 0L) return(data.frame(line = integer(), desc = character()))
  data.frame(
    line = hits,
    desc = vapply(m[hits], `[`, character(1), 2L),
    stringsAsFactors = FALSE
  )
}

file_data <- lapply(files, function(f) {
  list(
    name  = basename(f),
    path  = f,
    tests = extract_tests(f)
  )
})

total_tests <- sum(vapply(file_data, function(x) nrow(x$tests), integer(1)))
n_files     <- length(file_data)

# ---- HTML escaping --------------------------------------------------------
esc <- function(s) {
  s <- gsub("&",  "&amp;",  s, fixed = TRUE)
  s <- gsub("<",  "&lt;",   s, fixed = TRUE)
  s <- gsub(">",  "&gt;",   s, fixed = TRUE)
  s <- gsub("\"", "&quot;", s, fixed = TRUE)
  s
}

# ---- Build sections -------------------------------------------------------
sections <- character()
for (fd in file_data) {
  area <- area_map[[fd$name]] %||% ""
  # Simple %||% in case user doesn't load rlang
  if (is.null(area) || !nzchar(area)) area <- "&mdash;"

  n <- nrow(fd$tests)
  items <- character()
  if (n > 0L) {
    for (i in seq_len(n)) {
      items <- c(items, sprintf(
        '    <li>%s <span class="lineref">(line %d)</span></li>',
        esc(fd$tests$desc[i]), fd$tests$line[i]
      ))
    }
  } else {
    items <- '    <li class="empty">(no tests found)</li>'
  }
  sections <- c(sections, sprintf(
'<section>
  <div class="filehead">
    <h3><code>%s</code></h3>
    <span class="count">%d test%s</span>
  </div>
  <p class="area">%s</p>
  <ol>
%s
  </ol>
</section>',
    esc(fd$name), n, if (n == 1L) "" else "s", esc(area),
    paste(items, collapse = "\n")))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Assemble HTML --------------------------------------------------------
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>pigauto &mdash; test index</title>
<style>
  :root { --fg:#111827; --muted:#6b7280; --line:#e5e7eb; --soft:#f3f4f6;
          --accent:#059669; --badge:#eef2ff; --badge-fg:#3730a3; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 880px; margin: 2em auto; padding: 0 1.5em; color: var(--fg);
         line-height: 1.55; }
  h1 { border-bottom: 2px solid var(--fg); padding-bottom: .25em; }
  h3 { margin: 0; font-size: 16px; }
  h3 code { background: var(--soft); padding: 2px 8px; border-radius: 4px;
            font-size: 15px; }
  section { margin: 1.5em 0; padding: 1em 1.2em; border: 1px solid var(--line);
            border-radius: 6px; background: #fff; }
  .filehead { display: flex; align-items: baseline; justify-content: space-between;
              gap: 1em; }
  .count { background: var(--badge); color: var(--badge-fg); padding: 2px 10px;
           border-radius: 999px; font-size: 13px; font-weight: 600; }
  .area { color: var(--muted); margin: .4em 0 .8em 0; font-size: 13px; }
  ol { margin: 0; padding-left: 1.5em; font-size: 14px; }
  li { margin: 4px 0; }
  .lineref { color: var(--muted); font-size: 12px; margin-left: 6px; }
  .empty { color: var(--muted); font-style: italic; list-style: none;
           margin-left: -1.5em; }
  .summary { background: #ecfdf5; border-left: 4px solid var(--accent);
             padding: 1em 1.2em; border-radius: 4px; margin: 1.5em 0; }
  .summary b { color: #065f46; }
  .meta { color: var(--muted); font-size: 13px; }
  footer { color: var(--muted); font-size: 12px; margin-top: 3em;
           border-top: 1px solid var(--line); padding-top: 1em; }
</style>
</head>
<body>

<h1>pigauto &mdash; test index</h1>
<p class="meta">Static index of every <code>test_that()</code> block in <code>tests/testthat/</code>.
Generated ', timestamp, '. Does not run the tests &mdash; use <code>devtools::test()</code> for pass/fail status.</p>

<div class="summary">
<b>', total_tests, ' tests</b> across <b>', n_files, ' files</b>.
</div>

', paste(sections, collapse = "\n\n"), '

<footer>
Generator: <code>script/make_tests_html.R</code> &middot;
To run the tests: <code>devtools::test()</code> from the pigauto root.
</footer>

</body>
</html>
')

# Dual-write: one copy in script/ (dev artefact), one copy in
# pkgdown/assets/dev/ so pkgdown::build_site() exposes it on the web.
targets <- c("script/tests_overview.html",
             "pkgdown/assets/dev/tests_overview.html")
for (t in targets) {
  dir.create(dirname(t), showWarnings = FALSE, recursive = TRUE)
  writeLines(html, t)
  cat("Wrote ", t, " (", total_tests, " tests across ", n_files,
      " files)\n", sep = "")
}
