#!/usr/bin/env Rscript

# pkgdown intentionally renders every top-level Markdown file. Remove internal
# coordination and validation pages from the deployable site, sitemap, and
# search index.
`%||%` <- function(x, y) if (is.null(x)) y else x

# The site is built into _site/ (see `destination` in _pkgdown.yml). The
# environment override lets the same cleanup run against a temporary build.
site <- Sys.getenv("PIGAUTO_SITE_DIR", unset = "_site")

# pkgdown currently places its inline MathJax configuration after the loader.
# MathJax needs configuration first; otherwise lazily loaded TeX extensions
# such as boldsymbol can fail during typesetting. Move only the adjacent
# generated config/loader pair, then assert the order for every rendered page.
mathjax_pages <- list.files(
  site, pattern = "\\.html$", recursive = TRUE, full.names = TRUE
)
mathjax_pair <- paste0(
  "(?s)(<script\\b[^>]*\\bsrc=\"[^\"]*mathjax[^\"]*\"[^>]*>\\s*</script>)",
  "(<script\\s*>\\s*window\\.MathJax\\s*=\\s*\\{.*?\\};\\s*</script>)"
)
mathjax_loader <- "<script\\b[^>]*\\bsrc=\"[^\"]*mathjax[^\"]*\"[^>]*>"
mathjax_config <- "<script\\s*>\\s*window\\.MathJax\\s*="
for (path in mathjax_pages) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  reordered <- sub(mathjax_pair, "\\2\\1", html, perl = TRUE)
  # The Pandoc renderer used by pkgdown can serialize the HTML void element
  # <source> with an invalid closing tag. Remove that tag while retaining the
  # responsive source and its following img inside <picture>.
  cleaned <- gsub("</source\\s*>", "", reordered, perl = TRUE)
  # pkgdown's Bootstrap template currently emits curly quotes around this
  # favicon MIME type. Keep the generated link, but normalize its attribute.
  cleaned <- gsub(
    'type="”image/svg+xml”"', 'type="image/svg+xml"', cleaned, fixed = TRUE
  )
  if (grepl("</source\\s*>", cleaned, perl = TRUE)) {
    stop("Invalid closing source element remains in ", path)
  }
  if (grepl('type="”image/svg+xml”"', cleaned, fixed = TRUE)) {
    stop("Malformed SVG favicon MIME type remains in ", path)
  }
  if (!identical(cleaned, html)) writeLines(cleaned, path)

  loader_at <- regexpr(mathjax_loader, cleaned, perl = TRUE)[[1]]
  config_at <- regexpr(mathjax_config, cleaned, perl = TRUE)[[1]]
  if (loader_at > 0L && config_at > loader_at) {
    stop("MathJax configuration still follows its loader in ", path)
  }
}

# Give pkgdown's immediate reference redirects a document language. These
# noindex pages have no body content, but remain reachable through old URLs;
# declaring English keeps their document metadata complete for assistive
# technology and translation tools.
redirect_pages <- list.files(
  file.path(site, "reference"), pattern = "\\.html$",
  recursive = TRUE, full.names = TRUE
)
source_config <- yaml::read_yaml("_pkgdown.yml")
source_url <- sub("/+$", "", as.character(source_config$url))
if (length(source_url) != 1L || is.na(source_url) || !nzchar(source_url)) {
  stop("Could not identify the configured pkgdown site URL")
}
reference_url_prefix <- paste0(source_url, "/reference/")
for (path in redirect_pages) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  is_redirect <- grepl('<meta http-equiv="refresh"', html, fixed = TRUE) &&
    grepl('<meta name="robots" content="noindex"', html, fixed = TRUE)
  changed <- FALSE
  if (is_redirect) {
    refresh_tag_pattern <- '<meta http-equiv="refresh"[^>]*>'
    refresh_match <- regexpr(refresh_tag_pattern, html, perl = TRUE)
    if (refresh_match[[1]] < 0L) {
      stop("Could not identify the generated reference refresh tag in ", path)
    }
    refresh_tag <- regmatches(html, refresh_match)
    if (length(regmatches(
      html, gregexpr(refresh_tag_pattern, html, perl = TRUE)
    )[[1]]) != 1L) {
      stop("Expected one generated reference refresh tag in ", path)
    }
    destination <- sub(
      '^.*URL=([^";[:space:]]+).*$', "\\1", refresh_tag, perl = TRUE
    )
    if (startsWith(destination, reference_url_prefix)) {
      target <- substring(destination, nchar(reference_url_prefix) + 1L)
      target_path <- file.path(dirname(path), target)
      if (!grepl("^[[:alnum:]_.-]+\\.html$", target) ||
          !file.exists(target_path)) {
        stop("Generated reference redirect has an invalid local target in ", path)
      }
      local_refresh_tag <- sub(
        destination, paste0("./", target), refresh_tag, fixed = TRUE
      )
      html <- sub(refresh_tag, local_refresh_tag, html, fixed = TRUE)
      changed <- !identical(local_refresh_tag, refresh_tag)
    }
  }
  has_root_language <- grepl("<html[^>]*\\blang=", html, perl = TRUE)
  if (is_redirect && !has_root_language) {
    html <- sub("<html", '<html lang="en"', html, fixed = TRUE)
    changed <- TRUE
  }
  if (changed) {
    writeLines(html, path)
  }
}

# Harden the reviewed scientific tables without changing article source prose.
# These tables use column headers only; scope makes that relationship explicit
# to assistive technology. The wrapper works with extra.js to add a cue and a
# keyboard stop only when a table actually overflows.
responsive_table_pages <- c(
  "index.html",
  "articles/mixed-types.html",
  "articles/gnn-architecture.html",
  "articles/multiple-imputation.html",
  "reference/multi_impute_trees.html"
)
for (relative_path in responsive_table_pages) {
  path <- file.path(site, relative_path)
  if (!file.exists(path)) next
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  html <- gsub("<th>", "<th scope=\"col\">", html, fixed = TRUE)
  if (!grepl("class=\"table-responsive\"", html, fixed = TRUE)) {
    html <- gsub(
      "<table class=\"table\">",
      "<div class=\"table-responsive\">\n<table class=\"table\">",
      html, fixed = TRUE
    )
    html <- gsub("</table>", "</table>\n</div>", html, fixed = TRUE)
  }
  writeLines(html, path)
}

# Give the 11 reviewed example plots concise, output-level text alternatives.
# These descriptions report what the generated figures show; the axes-only
# gate rendering is described as rendered, without inferring why it is empty.
figure_alt <- c(
  "plot.pigauto_benchmark-1.png" = paste(
    "BM benchmark RMSE by trait: for trait1, the baseline bar is about 1.0",
    "and the pigauto bar about 0.8."
  ),
  "simulate_benchmark-1.png" = paste(
    "Grouped RMSE bars compare baseline and pigauto for two traits.",
    "The plotted values are similar for trait1 (about 0.64) and trait2 (about 0.32)."
  ),
  "plot.pigauto_pred-1.png" = paste(
    "Ranked predictions for Mass and Wing.Length across species; both rise",
    "toward higher ranks, with a steep final increase."
  ),
  "plot.pigauto_pred-2.png" = paste(
    "Observed versus predicted Mass and Wing.Length. Points are shown with",
    "a dashed equality line and a broad shaded interval; the displayed",
    "correlation is 1 for both traits."
  ),
  "plot.pigauto_pred-3.png" = paste(
    "Mass predictions by species rank with a shaded uncertainty band.",
    "The line rises gradually, then sharply for the highest-ranked species;",
    "the band also widens at the upper end."
  ),
  "plot.pigauto_pred-4.png" = paste(
    "Observed versus predicted Mass and Wing.Length for held-out cells.",
    "Each panel has a dashed equality line and shaded interval; Mass has",
    "three visible points and Wing.Length includes one point well below equality."
  ),
  "plot.pigauto_fit-1.png" = paste(
    "Blend gate per trait for Mass and Wing.Length. The fit uses gnn = FALSE,",
    "so both GNN gate values are zero. The chart states that there is no GNN",
    "contribution; this is the expected baseline-only result."
  ),
  "plot.pigauto_fit-2.png" = paste(
    "Conformal prediction scores for Mass and Wing.Length on the latent",
    "scale. The dashed horizontal line marks their median; its value is",
    "labelled above the bars."
  ),
  "plot_uncertainty-1.png" = paste(
    "Mass predictions for three species sorted by prediction, with a shaded",
    "uncertainty band around the line; both prediction and band rise across",
    "the displayed species."
  ),
  "plot_history_gg-1.png" = paste(
    "Training history showing validation loss decreasing from 1.0 at epoch 1",
    "to 0.7 at epoch 2 and 0.5 at epoch 3."
  ),
  "plot_comparison-1.png" = paste(
    "RMSE comparison for Mass: the baseline point is about 1.0 and the",
    "pigauto triangle about 0.8. The legend marks whether the GNN is better",
    "or worse than baseline."
  )
)
reference_pages <- list.files(
  file.path(site, "reference"), pattern = "\\.html$",
  recursive = TRUE, full.names = TRUE
)
for (path in reference_pages) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  original <- html
  for (filename in names(figure_alt)) {
    needle <- paste0('(<img[^>]*src="', filename,
                     '"[^>]*alt=")[^"]*(")')
    replacement <- paste0("\\1", figure_alt[[filename]], "\\2")
    html <- gsub(needle, replacement, html, perl = TRUE)
  }
  if (!identical(html, original)) writeLines(html, path)
}

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

# Published pkgdown metadata feeds cross-package article links. Remove only
# this retired article key, preserving every other metadata byte and field.
metadata_path <- file.path(site, "pkgdown.yml")
if (!file.exists(metadata_path)) stop("Candidate site is missing pkgdown.yml")
metadata_bytes <- readBin(
  metadata_path, what = "raw", n = file.info(metadata_path)$size
)
metadata <- yaml::yaml.load(rawToChar(metadata_bytes))
retired_article_key <- "articles/simulation-study"
if (retired_article_key %in% names(metadata$articles)) {
  expected_metadata <- metadata
  expected_metadata$articles[[retired_article_key]] <- NULL
  # Locate complete lines in the original bytes, including their LF/CRLF
  # terminators when present. No surviving byte is reserialized.
  line_breaks <- which(metadata_bytes == as.raw(10L))
  line_starts <- c(1L, line_breaks + 1L)
  line_ends <- c(line_breaks, length(metadata_bytes))
  nonempty <- line_starts <= line_ends
  line_starts <- line_starts[nonempty]
  line_ends <- line_ends[nonempty]
  metadata_lines <- vapply(seq_along(line_starts), function(i) {
    rawToChar(metadata_bytes[seq.int(line_starts[i], line_ends[i])])
  }, character(1))
  retired_line <- grep(
    "^  articles/simulation-study:[[:space:]]", metadata_lines
  )
  if (length(retired_line) != 1L) {
    stop("Could not identify the exact retired article metadata line")
  }
  retired_bytes <- seq.int(line_starts[retired_line], line_ends[retired_line])
  cleaned_bytes <- metadata_bytes[-retired_bytes]
  cleaned_metadata <- yaml::yaml.load(rawToChar(cleaned_bytes))
  if (!identical(cleaned_metadata, expected_metadata)) {
    stop("Retired article cleanup would change unrelated pkgdown metadata")
  }
  writeBin(cleaned_bytes, metadata_path)
}

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
  # pkgdown's generated article index starts its in-content list at h3 after
  # the page h1. Raise it one level so the outline has no skipped heading.
  index_html <- sub(
    "<h3>All vignettes</h3>", "<h2>All vignettes</h2>", index_html,
    fixed = TRUE
  )
  writeLines(index_html, article_index)
}

# The generated Markdown search aids have their own article lists. Remove the
# retired entry there too, preserving all surrounding bytes and line endings.
remove_markdown_entry <- function(path, target) {
  if (!file.exists(path)) return(invisible(FALSE))
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  breaks <- which(bytes == as.raw(10L))
  starts <- c(1L, breaks + 1L)
  ends <- c(breaks, length(bytes))
  nonempty <- starts <= ends
  starts <- starts[nonempty]
  ends <- ends[nonempty]
  lines <- vapply(seq_along(starts), function(i) {
    line <- bytes[seq.int(starts[i], ends[i])]
    if (length(line) && identical(line[[length(line)]], as.raw(10L))) {
      line <- line[-length(line)]
    }
    if (length(line) && identical(line[[length(line)]], as.raw(13L))) {
      line <- line[-length(line)]
    }
    rawToChar(line)
  }, character(1))
  hit <- which(grepl(target, lines, fixed = TRUE))
  if (!length(hit)) return(invisible(FALSE))
  if (length(hit) != 1L) stop("Expected one retired link in ", path)
  first <- hit[[1]]
  while (first > 1L && !grepl("^- ", lines[[first]], perl = TRUE)) {
    first <- first - 1L
  }
  if (!grepl("^- ", lines[[first]], perl = TRUE)) {
    stop("Could not identify the retired Markdown list item in ", path)
  }
  last <- hit[[1]]
  while (last < length(lines) && !grepl("^- ", lines[[last + 1L]], perl = TRUE)) {
    last <- last + 1L
  }
  if (!grepl("\\):[[:space:]]*$", lines[[last]], perl = TRUE)) {
    stop("Could not identify the end of the retired Markdown link in ", path)
  }
  cleaned_lines <- lines[-seq.int(first, last)]
  if (any(grepl(target, cleaned_lines, fixed = TRUE))) {
    stop("Retired link remains in ", path)
  }
  remove_start <- starts[first]
  remove_end <- ends[last]
  cleaned_bytes <- bytes[-seq.int(remove_start, remove_end)]
  writeBin(cleaned_bytes, path)
  invisible(TRUE)
}
remove_markdown_entry(
  file.path(site, "articles", "index.md"),
  "articles/simulation-study.md"
)
remove_markdown_entry(
  file.path(site, "llms.txt"),
  "articles/simulation-study.md"
)

# The generated not-found page is served at the requested local route. Keep its
# skip link on that page instead of sending keyboard users to the hosted site.
not_found_page <- file.path(site, "404.html")
if (file.exists(not_found_page)) {
  not_found_html <- paste(readLines(not_found_page, warn = FALSE), collapse = "\n")
  skip_link_pattern <- '<a href="[^"]*" class="visually-hidden-focusable">Skip to contents</a>'
  skip_link_count <- lengths(regmatches(
    not_found_html, gregexpr(skip_link_pattern, not_found_html, perl = TRUE)
  ))
  if (skip_link_count != 1L) {
    stop("Could not identify exactly one generated 404 skip link")
  }
  not_found_html <- sub(
    skip_link_pattern,
    '<a href="#main" class="visually-hidden-focusable">Skip to contents</a>',
    not_found_html, perl = TRUE
  )
  # A custom 404 can be returned for an arbitrary nested path. Its search
  # input therefore cannot use pkgdown's relative search.json URL: browsers
  # would request the index beneath the missing path. Point only this page's
  # search input at the configured site-root index.
  search_index_pattern <- 'data-search-index="[^"]*"'
  search_index_locations <- gregexpr(
    search_index_pattern, not_found_html, perl = TRUE
  )
  search_index_attributes <- regmatches(
    not_found_html, search_index_locations
  )[[1]]
  if (length(search_index_attributes) != 1L) {
    stop("Could not identify exactly one generated 404 search index")
  }
  search_index_url <- paste0(source_url, "/search.json")
  not_found_html <- sub(
    search_index_pattern,
    paste0('data-search-index="', search_index_url, '"'),
    not_found_html, perl = TRUE
  )
  if (!grepl(
    paste0('data-search-index="', search_index_url, '"'),
    not_found_html, fixed = TRUE
  )) {
    stop("Could not set the generated 404 search index to the site root")
  }
  malformed_mathjax_source <- paste0(source_url, "/NA")
  script_tag_locations <- gregexpr("<script\\b[^>]*>", not_found_html, perl = TRUE)
  script_tags <- regmatches(not_found_html, script_tag_locations)[[1]]
  if (length(script_tags) == 1L && identical(script_tags, character(0))) {
    script_tags <- character()
  }
  malformed_tags <- which(vapply(
    script_tags, grepl, logical(1), pattern = malformed_mathjax_source,
    fixed = TRUE
  ))
  if (length(malformed_tags) > 1L) {
    stop("Found multiple same-site /NA script sources on the 404 page")
  }
  if (length(malformed_tags) == 1L) {
    script_index <- malformed_tags[[1]]
    script_start <- script_tag_locations[[1]][[script_index]]
    script_tag <- script_tags[[script_index]]
    body_start <- script_start + nchar(script_tag)
    script_body <- substring(not_found_html, body_start)
    script_end <- regexpr("</script\\s*>", script_body, perl = TRUE)[1]
    if (script_end < 0L) {
      stop("The malformed 404 script tag has no closing script element")
    }
    script_body <- substring(script_body, 1L, script_end - 1L)
    if (!grepl("^\\s*window\\.MathJax\\s*=", script_body, perl = TRUE)) {
      stop("The malformed 404 script body is not a MathJax configuration")
    }
    malformed_src_attribute <- paste0(' src="', malformed_mathjax_source, '"')
    repaired_tag <- sub(
      malformed_src_attribute, "", script_tag, fixed = TRUE
    )
    if (identical(repaired_tag, script_tag)) {
      stop("Could not remove the exact generated same-site /NA src attribute")
    }
    not_found_html <- sub(script_tag, repaired_tag, not_found_html, fixed = TRUE)
  }
  remaining_script_locations <- gregexpr(
    "<script\\b[^>]*>", not_found_html, perl = TRUE
  )
  remaining_script_tags <- regmatches(not_found_html, remaining_script_locations)[[1]]
  if (length(remaining_script_tags) == 1L && identical(remaining_script_tags, character(0))) {
    remaining_script_tags <- character()
  }
  if (any(vapply(
    remaining_script_tags, grepl, logical(1), pattern = malformed_mathjax_source,
    fixed = TRUE
  ))) {
    stop("Same-site /NA script source remains on the 404 page")
  }
  # The earlier all-page ordering pass cannot recognize this 404 config until
  # its malformed same-site src is removed. Reorder this one page now and
  # assert the final output before writing it.
  not_found_html <- sub(mathjax_pair, "\\2\\1", not_found_html, perl = TRUE)
  not_found_loader_at <- regexpr(mathjax_loader, not_found_html, perl = TRUE)[[1]]
  not_found_config_at <- regexpr(mathjax_config, not_found_html, perl = TRUE)[[1]]
  if (not_found_loader_at > 0L && not_found_config_at > not_found_loader_at) {
    stop("MathJax configuration still follows its loader in 404.html")
  }
  writeLines(not_found_html, not_found_page)
}

sitemap <- file.path(site, "sitemap.xml")
if (file.exists(sitemap)) {
  lines <- readLines(sitemap, warn = FALSE)
  pattern <- paste0("/(", paste(stems, collapse = "|"), ")\\.html")
  keep_internal <- !grepl(pattern, lines)
  keep_retired <- !vapply(lines, function(line) {
    any(vapply(retired_routes, grepl, logical(1), x = line, fixed = TRUE))
  }, logical(1))
  keep_not_found <- !grepl("/404.html", lines, fixed = TRUE)
  writeLines(lines[keep_internal & keep_retired & keep_not_found], sitemap)
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
  entries <- lapply(entries[keep], function(entry) {
    path <- entry$path
    navigable <- is.character(path) && length(path) == 1L &&
      !is.na(path) && nzchar(path)
    if (!navigable) return(entry)
    if (!("id" %in% names(entry))) {
      stop("Navigable search record is missing its id key")
    }
    id <- entry[["id"]]
    empty_object <- is.list(id) && length(id) == 0L && !is.null(names(id))
    if (is.null(id) || empty_object) {
      entry[["id"]] <- ""
    } else if (!(is.character(id) && length(id) == 1L && !is.na(id))) {
      stop("Navigable search record has an invalid id type or length")
    }
    entry
  })
  jsonlite::write_json(entries, search, auto_unbox = TRUE, null = "null")
}
