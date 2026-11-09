# Package check for a frozen development candidate, independent of the CRAN archive.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
sha <- args[3]
stopifnot(identical(readLines(file.path(root, "SOURCE_SHA"), warn = FALSE), sha))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
result <- rcmdcheck::rcmdcheck(path = root, args = "--no-manual",
  check_dir = normalizePath(out), error_on = "never", quiet = FALSE)
saveRDS(result, file.path(out, "check-result.rds"))
jsonlite::write_json(list(source_sha = sha, errors = result$errors,
  warnings = result$warnings, notes = result$notes,
  manual = FALSE, r_version = R.version.string), file.path(out, "check-summary.json"),
  auto_unbox = TRUE, pretty = TRUE)
stopifnot(length(result$errors) == 0L, length(result$warnings) == 0L)
cat("PACKAGE_CHECK_PASS", sha, "notes", length(result$notes), "\n")
