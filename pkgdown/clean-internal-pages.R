#!/usr/bin/env Rscript

# pkgdown intentionally renders every top-level Markdown file. Remove internal
# agent coordination pages from the deployable site, sitemap, and search index.
`%||%` <- function(x, y) if (is.null(x)) y else x

# The site is built into _site/ (see `destination` in _pkgdown.yml).
site <- "_site"

stems <- c("AGENTS", "CLAUDE", "goodagents")
unlink(file.path(site, c(paste0(stems, ".html"), paste0(stems, ".md"))))

sitemap <- file.path(site, "sitemap.xml")
if (file.exists(sitemap)) {
  lines <- readLines(sitemap, warn = FALSE)
  pattern <- paste0("/(", paste(stems, collapse = "|"), ")\\.html")
  writeLines(lines[!grepl(pattern, lines)], sitemap)
}

search <- file.path(site, "search.json")
if (file.exists(search)) {
  entries <- jsonlite::fromJSON(search, simplifyVector = FALSE)
  pattern <- paste0("/(", paste(stems, collapse = "|"), ")\\.html$")
  keep <- !vapply(entries, function(entry) {
    path <- entry$path %||% ""
    length(path) == 1L && grepl(pattern, path)
  }, logical(1))
  jsonlite::write_json(entries[keep], search, auto_unbox = TRUE)
}
