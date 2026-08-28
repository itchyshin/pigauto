#!/usr/bin/env Rscript

fail <- function(...) stop(sprintf(...), call. = FALSE)
candidate_current <- function() {
  pointer <- file.path("checkpoints", "pigauto-011", "CURRENT")
  if (!file.exists(pointer)) return("")
  normalizePath(readLines(pointer, n = 1L, warn = FALSE), mustWork = TRUE)
}
need_path <- function(name, fallback) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) value <- fallback
  if (!nzchar(value)) fail("%s is missing and no frozen CURRENT candidate exists", name)
  normalizePath(value, mustWork = TRUE)
}

current <- candidate_current()
manifest <- if (nzchar(current)) readRDS(file.path(current, "manifest.rds")) else NULL
tarball <- need_path("PIGAUTO_TARBALL", if (is.null(manifest)) "" else file.path(current, manifest$tarball))
library_path <- need_path("PIGAUTO_LIBRARY", if (nzchar(current)) file.path(current, "candidate-library") else "")
source_root <- normalizePath(Sys.getenv("PIGAUTO_SOURCE", unset = "."), mustWork = TRUE)
entries <- utils::untar(tarball, list = TRUE)
if (!length(entries)) fail("Candidate tarball inventory is empty")
if (any(grepl("(^|/)BACE(/|$)", entries))) fail("Candidate tarball contains BACE/")
if (any(grepl("(^|/)_pkgdown[.]yml$", entries))) fail("Candidate tarball contains _pkgdown.yml")

desc_path <- file.path(library_path, "pigauto", "DESCRIPTION")
desc <- read.dcf(desc_path)
dependency_fields <- intersect(c("Depends", "Imports", "Suggests", "Enhances", "LinkingTo"), colnames(desc))
dependencies <- paste(desc[1, dependency_fields, drop = TRUE], collapse = "\n")
has_bace_dependency <- function(x) grepl("(^|[^[:alnum:].])BACE([^[:alnum:].]|$)", x, perl = TRUE)
if (!all(vapply(c("BACE", "BACE, torch", "foo, BACE (>= 1)"), has_bace_dependency, logical(1)))) {
  fail("BACE dependency detector failed its planted positive controls")
}
if (has_bace_dependency(dependencies)) {
  fail("Installed DESCRIPTION retains a BACE dependency")
}

pkgdown <- paste(readLines(file.path(source_root, "_pkgdown.yml"), warn = FALSE), collapse = "\n")
if (grepl("fit_baseline_bace|BACE::|bace_final_imp", pkgdown, ignore.case = TRUE, perl = TRUE)) {
  fail("Source pkgdown configuration retains an installed BACE entry")
}

ns <- loadNamespace("pigauto", lib.loc = library_path)
exports <- getNamespaceExports(ns)
if ("fit_baseline_bace" %in% exports) fail("fit_baseline_bace() remains exported")
if (exists("fit_baseline_bace", envir = ns, inherits = FALSE)) {
  fail("fit_baseline_bace() remains bound in the installed namespace")
}

rd <- tools::Rd_db("pigauto", lib.loc = library_path)
rd_text <- vapply(rd, function(topic) paste(as.character(topic), collapse = "\n"), character(1))
if (any(grepl("fit_baseline_bace", rd_text, fixed = TRUE))) {
  fail("Installed help retains fit_baseline_bace()")
}

cat("G4 BACE surface verified\n")
