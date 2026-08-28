#!/usr/bin/env Rscript

fail <- function(...) stop(sprintf(...), call. = FALSE)
root <- normalizePath(Sys.getenv("PIGAUTO_SOURCE", unset = "."), mustWork = TRUE)
paths <- c(
  "DESCRIPTION", "README.md", "NEWS.md", "_pkgdown.yml",
  list.files(file.path(root, "R"), pattern = "[.]R$", full.names = FALSE),
  list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = FALSE),
  list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$", full.names = FALSE)
)
paths[5:length(paths)] <- c(
  file.path("R", list.files(file.path(root, "R"), pattern = "[.]R$")),
  file.path("man", list.files(file.path(root, "man"), pattern = "[.]Rd$")),
  file.path("vignettes", list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$"))
)
paths <- unique(paths[file.exists(file.path(root, paths))])
asset_paths <- list.files(file.path(root, "pkgdown", "assets"), pattern = "[.]html$",
                          recursive = TRUE, full.names = TRUE)
paths <- c(file.path(root, paths), asset_paths)
text <- paste(vapply(paths, function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}, character(1)), collapse = "\n")
flat <- gsub("[[:space:]]+", " ", text)
news <- paste(readLines(file.path(root, "NEWS.md"), warn = FALSE), collapse = "\n")
news_flat <- gsub("[[:space:]]+", " ", news)

forbidden <- c(
  "competitive with BACE",
  "beating BACE",
  "accuracy switch",
  "properly calibrated (standard errors|SEs)",
  "honest (standard errors|SEs)",
  "guaranteed 95% coverage",
  "distribution-free 95% coverage",
  "(guarantees?|provides?) (an? )?optimal sampling",
  "optimal sampling (strategy|design)",
  "provably optimal"
)
hits <- forbidden[vapply(forbidden, grepl, logical(1), x = flat, ignore.case = TRUE, perl = TRUE)]
if (length(hits)) fail("Forbidden public claim pattern(s): %s", paste(hits, collapse = ", "))

proximity <- c(
  "bace_final_imp.{0,240}(proper multiple-imputation|Rubin [(]1987[)])",
  "PMM.{0,240}(pool_mi|Rubin|honest standard errors|properly calibrated SEs)",
  "tree uncertainty.{0,120}downstream standard errors.{0,120}Rubin",
  "pooled via Rubin's rules.{0,120}default path|Rubin's rules stable"
)
hits <- proximity[vapply(proximity, grepl, logical(1), x = news_flat, ignore.case = TRUE, perl = TRUE)]
if (length(hits)) fail("Unsupported inferential proximity pattern(s): %s", paste(hits, collapse = ", "))

required <- c(
  "current supported downstream route begins with `multi_impute_analysis()`",
  "prediction diagnostics rather than validated inferential imputations",
  "not an optimal sampling guarantee",
  "does not certify package-wide coverage"
)
missing <- required[!vapply(required, function(pattern) {
  grepl(tolower(pattern), tolower(flat), fixed = TRUE)
}, logical(1))]
if (length(missing)) fail("Required claim boundary missing: %s", paste(missing, collapse = "; "))

positive_controls <- c(
  "bace_final_imp() produces\nproper multiple-imputation draws",
  "PMM draws can be passed to\npool_mi() for Rubin inference",
  "tree uncertainty can propagate into\ndownstream standard errors through Rubin pooling",
  "datasets are pooled via Rubin's rules\nand become the default path"
)
detected <- vapply(positive_controls, function(x) {
  probe <- gsub("[[:space:]]+", " ", x)
  any(vapply(proximity, grepl, logical(1), x = probe, ignore.case = TRUE, perl = TRUE))
}, logical(1))
if (!all(detected)) fail("Claim verifier failed its planted wrapped-text negative controls")

withdrawn_names <- c(
  "bench_avonet9993_bace.html", "bench_avonet9993_bace_index.html",
  "bench_avonet9993_bace_n3000.html", "bench_bace_avonet_head_to_head.html",
  "bench_bien.html", "bench_fishbase.html", "bench_pantheria_bace_head_to_head.html",
  "calibration_grid.html", "pantheria_summary.html", "phase8_summary.html",
  "tests_overview.html"
)
for (name in withdrawn_names) {
  path <- file.path(root, "pkgdown", "assets", "dev", name)
  page <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!grepl("Historical development page withdrawn", page, fixed = TRUE)) {
    fail("Unsafe source dev page is not tombstoned: %s", path)
  }
}

site <- Sys.getenv("PIGAUTO_SITE", unset = "")
if (!nzchar(site)) {
  pointer <- file.path(root, "checkpoints", "pigauto-011", "CURRENT")
  if (file.exists(pointer)) {
    current <- normalizePath(readLines(pointer, n = 1L, warn = FALSE), mustWork = TRUE)
    site <- file.path(current, "site")
  }
}
if (nzchar(site)) {
  site <- normalizePath(site, mustWork = TRUE)
  site_paths <- list.files(site, pattern = "[.]html$", recursive = TRUE, full.names = TRUE)
  if (!length(site_paths)) fail("Rendered site contains no HTML files")
  site_text <- paste(vapply(site_paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1)), collapse = "\n")
  site_flat <- gsub("[[:space:]]+", " ", site_text)
  site_hits <- forbidden[vapply(forbidden, grepl, logical(1), x = site_flat,
                                ignore.case = TRUE, perl = TRUE)]
  if (length(site_hits)) fail("Forbidden rendered-site claim pattern(s): %s", paste(site_hits, collapse = ", "))
  for (name in withdrawn_names) {
    path <- file.path(site, "dev", name)
    if (!file.exists(path)) fail("Required historical tombstone missing from rendered site: %s", path)
    page <- paste(readLines(path, warn = FALSE), collapse = "\n")
    if (!grepl("Historical development page withdrawn", page, fixed = TRUE)) {
      fail("Unsafe historical dev page is not tombstoned in rendered site: %s", path)
    }
  }
}

cat("G9 claim boundary verified\n")
