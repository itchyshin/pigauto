#!/usr/bin/env Rscript

# pkgdown intentionally renders every top-level Markdown file. Remove internal
# coordination and validation pages from the deployable site, sitemap, and
# search index.
`%||%` <- function(x, y) if (is.null(x)) y else x

# The site is built into _site/ (see `destination` in _pkgdown.yml). The
# environment override lets the same cleanup run against a temporary build.
site <- Sys.getenv("PIGAUTO_SITE_DIR", unset = "_site")

stems <- c("AGENTS", "CLAUDE", "goodagents", "VALIDATION_LEDGER")
unlink(file.path(site, c(paste0(stems, ".html"), paste0(stems, ".md"))))

# pkgdown renders this build-excluded historical article from the source tree
# even though it is not listed in the article menu. Keep the source and its
# campaign evidence in Git, but retire the rendered HTML and Markdown, plus
# pkgdown's nested HTML redirect, from the public site.
retired_routes <- c(
  "articles/simulation-study.html",
  "articles/simulation-study.md",
  "articles/articles/simulation-study.html",
  "articles/articles/simulation-study.md"
)
unlink(file.path(site, retired_routes))

article_index <- file.path(site, "articles", "index.html")
if (file.exists(article_index)) {
  index_html <- paste(readLines(article_index, warn = FALSE), collapse = "\n")
  index_pattern <- paste0(
    "<dt><a href=\"simulation-study\\.html\">[^<]*</a></dt>",
    "\\s*<dd>\\s*</dd>"
  )
  index_html <- gsub(index_pattern, "", index_html, perl = TRUE)
  if (grepl("href=\"simulation-study\\.html\"", index_html)) {
    stop("Could not remove the retired simulation-study link from articles/index.html")
  }
  writeLines(index_html, article_index)
}

sitemap <- file.path(site, "sitemap.xml")
if (file.exists(sitemap)) {
  lines <- readLines(sitemap, warn = FALSE)
  pattern <- paste0("/(", paste(stems, collapse = "|"), ")\\.html")
  keep_internal <- !grepl(pattern, lines)
  keep_retired <- !vapply(lines, function(line) {
    any(vapply(retired_routes, grepl, logical(1), x = line, fixed = TRUE))
  }, logical(1))
  writeLines(lines[keep_internal & keep_retired], sitemap)
}

search <- file.path(site, "search.json")
if (file.exists(search)) {
  entries <- jsonlite::fromJSON(search, simplifyVector = FALSE)
  pattern <- paste0("/(", paste(stems, collapse = "|"), ")\\.html$")
  keep <- !vapply(entries, function(entry) {
    path <- entry$path %||% ""
    length(path) == 1L && (
      grepl(pattern, path) ||
        any(vapply(retired_routes, grepl, logical(1), x = path, fixed = TRUE))
    )
  }, logical(1))
  jsonlite::write_json(entries[keep], search, auto_unbox = TRUE)
}
