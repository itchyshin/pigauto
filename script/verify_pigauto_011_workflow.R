#!/usr/bin/env Rscript

fail <- function(...) stop(sprintf(...), call. = FALSE)
library_path <- Sys.getenv("PIGAUTO_LIBRARY", unset = "")
if (!nzchar(library_path)) {
  pointer <- file.path("checkpoints", "pigauto-011", "CURRENT")
  if (file.exists(pointer)) {
    current <- normalizePath(readLines(pointer, n = 1L, warn = FALSE), mustWork = TRUE)
    library_path <- file.path(current, "candidate-library")
  }
}
if (!nzchar(library_path)) fail("PIGAUTO_LIBRARY is missing and no frozen CURRENT candidate exists")
library_path <- normalizePath(library_path, mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))

installed <- normalizePath(find.package("pigauto", lib.loc = library_path), mustWork = TRUE)
if (!identical(installed, normalizePath(file.path(library_path, "pigauto"), mustWork = TRUE))) {
  fail("pigauto was not loaded from the candidate library")
}
script <- system.file("examples", "novice-workflow.R", package = "pigauto", lib.loc = library_path)
if (!nzchar(script) || !file.exists(script)) fail("Installed novice workflow is missing")

expr <- parse(script)
if (length(expr) != 6L) fail("Installed novice workflow must contain exactly six expressions")

ns <- asNamespace("pigauto")
original_get_device <- get("get_device", envir = ns, inherits = FALSE)
assignInNamespace("get_device", function() torch::torch_device("cpu"), ns = "pigauto")
on.exit(assignInNamespace("get_device", original_get_device, ns = "pigauto"), add = TRUE)
if (!identical(as.character(get("get_device", envir = ns, inherits = FALSE)()), "cpu")) {
  fail("Novice verifier could not force the candidate workflow to CPU")
}

run_dir <- tempfile("pigauto-011-workflow-")
dir.create(run_dir)
old_dir <- setwd(run_dir)
on.exit(setwd(old_dir), add = TRUE)
old_browser <- getOption("browser")
options(browser = function(...) invisible(TRUE))
on.exit(options(browser = old_browser), add = TRUE)

env <- new.env(parent = globalenv())
sys.source(script, envir = env)
if (!inherits(env$result, "pigauto_result")) fail("Workflow did not return pigauto_result")
if (!inherits(env$result$check, "pigauto_check")) fail("Result does not retain pigauto_check")
if (!identical(env$completed, pigauto::completed_data(env$result))) fail("completed_data() disagrees with result")
if (!identical(rownames(env$traits), rownames(env$completed))) fail("Species row identity changed")
if (!identical(names(env$traits), names(env$completed))) fail("Trait column order changed")
if (!identical(vapply(env$traits, class, character(1)), vapply(env$completed, class, character(1)))) {
  fail("Trait column classes changed")
}
observed_unchanged <- vapply(names(env$traits), function(name) {
  observed <- !is.na(env$traits[[name]])
  identical(env$completed[[name]][observed], env$traits[[name]][observed])
}, logical(1))
if (!all(observed_unchanged)) {
  fail("Observed cells changed in: %s", paste(names(observed_unchanged)[!observed_unchanged], collapse = ", "))
}
if (anyNA(env$completed)) fail("Novice fixture retains unresolved missing cells")
if (!file.exists("pigauto_report.html")) fail("Workflow did not write pigauto_report.html")
report <- paste(readLines("pigauto_report.html", warn = FALSE), collapse = "\n")
if (!grepl("Completed data", report, fixed = TRUE)) fail("Report omits completed-data role")
required_report_boundaries <- c(
  "How to use outputs",
  "A pigauto_result does not authorize inference",
  "multi_impute_analysis() only in its documented narrow route"
)
missing_report_boundaries <- required_report_boundaries[
  !vapply(required_report_boundaries, grepl, logical(1), x = report, fixed = TRUE)
]
if (length(missing_report_boundaries)) {
  fail("Report omits boundary text: %s", paste(missing_report_boundaries, collapse = "; "))
}

cat("G8 installed novice workflow verified\n")
