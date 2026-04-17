#!/usr/bin/env Rscript
#
# script/make_bench_scaling_v090_html.R
#
# Minimal wrapper: reads bench_scaling_v090.md + embeds bench_scaling_v090.png
# as a base64 data URI, then writes self-contained HTML to:
#   script/bench_scaling_v090.html
#   pkgdown/assets/dev/bench_scaling_v090.html
#
# Run with:
#   Rscript script/make_bench_scaling_v090_html.R

md_path  <- "script/bench_scaling_v090.md"
png_path <- "script/bench_scaling_v090.png"

if (!file.exists(md_path))  stop("Missing ", md_path)
if (!file.exists(png_path)) stop("Missing ", png_path)

# Read markdown text
md_lines <- readLines(md_path, warn = FALSE)

# Encode PNG as base64 data URI
png_bytes  <- readBin(png_path, what = "raw", n = file.info(png_path)$size)
png_b64    <- paste0("data:image/png;base64,",
                     paste(base64enc::base64encode(png_bytes), collapse = ""))

# Convert markdown to very simple HTML:
#   - H1 headings  -> <h1>
#   - H2 headings  -> <h2>
#   - H3 headings  -> <h3>
#   - Fenced code  -> <pre>
#   - Pipe tables  -> <pre> (preserve as monospace)
#   - Blank lines  -> paragraph breaks
#   - Other lines  -> text lines

convert_md <- function(lines) {
  out <- character()
  i   <- 1L
  n   <- length(lines)

  while (i <= n) {
    line <- lines[[i]]

    # Fenced code block
    if (grepl("^```", line)) {
      out <- c(out, "<pre><code>")
      i <- i + 1L
      while (i <= n && !grepl("^```", lines[[i]])) {
        out <- c(out, htmlEscape(lines[[i]]))
        i <- i + 1L
      }
      out <- c(out, "</code></pre>")
      i <- i + 1L
      next
    }

    # Table rows (pipe table)
    if (grepl("^\\|", line)) {
      out <- c(out, "<pre>")
      while (i <= n && grepl("^\\|", lines[[i]])) {
        out <- c(out, htmlEscape(lines[[i]]))
        i <- i + 1L
      }
      out <- c(out, "</pre>")
      next
    }

    # ATX headings
    if (grepl("^# ", line)) {
      out <- c(out, paste0("<h1>", htmlEscape(sub("^# +", "", line)), "</h1>"))
      i <- i + 1L; next
    }
    if (grepl("^## ", line)) {
      out <- c(out, paste0("<h2>", htmlEscape(sub("^## +", "", line)), "</h2>"))
      i <- i + 1L; next
    }
    if (grepl("^### ", line)) {
      out <- c(out, paste0("<h3>", htmlEscape(sub("^### +", "", line)), "</h3>"))
      i <- i + 1L; next
    }

    # Blank line
    if (grepl("^\\s*$", line)) {
      i <- i + 1L; next
    }

    # List items
    if (grepl("^- ", line)) {
      out <- c(out, paste0("<li>", htmlEscape(sub("^- +", "", line)), "</li>"))
      i <- i + 1L; next
    }

    # Plain text -> paragraph
    out <- c(out, paste0("<p>", htmlEscape(line), "</p>"))
    i <- i + 1L
  }
  paste(out, collapse = "\n")
}

htmlEscape <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  s <- gsub(">", "&gt;",  s, fixed = TRUE)
  s
}

body_html <- convert_md(md_lines)

html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>pigauto v0.9.0 scaling curve</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 960px; margin: 2rem auto; padding: 0 1.5rem;
         color: #111827; line-height: 1.6; }
  h1 { border-bottom: 2px solid #059669; padding-bottom: .25em; color: #059669; }
  h2 { margin-top: 2em; border-bottom: 1px solid #d1d5db; padding-bottom: .15em; }
  h3 { margin-top: 1.5em; color: #374151; }
  pre { background: #f3f4f6; padding: 1em; border-radius: 4px;
        overflow-x: auto; font-size: 13px; line-height: 1.45;
        white-space: pre; }
  code { font-size: 13px; }
  p { margin: .5em 0; }
  li { margin: 4px 0; }
  ul { padding-left: 1.5em; }
  img { max-width: 100%; height: auto; display: block; margin: 1.5em auto;
        border: 1px solid #e5e7eb; border-radius: 4px; }
  footer { margin-top: 3em; padding-top: 1em; border-top: 1px solid #e5e7eb;
           font-size: 12px; color: #6b7280; }
</style>
</head>
<body>

', body_html, '

<h2>Scaling curve</h2>
<img src="', png_b64, '" alt="pigauto v0.9.0 scaling curve chart">

<footer>
Source: <code>script/bench_scaling_v090.R</code> &middot;
Results: <code>script/bench_scaling_v090.rds</code> &middot;
Report: <code>pkgdown/assets/dev/bench_scaling_v090.html</code> &middot;
Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M"), '
</footer>
</body>
</html>
')

targets <- c("script/bench_scaling_v090.html",
             "pkgdown/assets/dev/bench_scaling_v090.html")
for (t in targets) {
  dir.create(dirname(t), showWarnings = FALSE, recursive = TRUE)
  writeLines(html, t)
  cat("Wrote", t, "(", file.size(t), "bytes)\n")
}
