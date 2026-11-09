# Verify assertions that intentionally cannot run from an installed tarball.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
sha <- args[3]
stopifnot(identical(readLines(file.path(root, "SOURCE_SHA")), sha))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
setwd(root)
Sys.setenv(PIGAUTO_REQUIRE_SOURCE = "true", NOT_CRAN = "true")
r <- devtools::test(
  filter = "^(community-surface|coverage-workflow|gha-helpers|docs-phytools-warning|bace-compat-eval)$",
  reporter = "summary", stop_on_failure = FALSE)
saveRDS(r, file.path(out, "results.rds"))
d <- as.data.frame(r)
d$result <- NULL
write.csv(d, file.path(out, "results.csv"), row.names = FALSE)
ok <- nrow(d) > 0L && !any(d$failed) && !any(d$error) &&
  !any(d$warning) && !any(d$skipped)
jsonlite::write_json(list(source_sha = sha, status = if (ok) "ok" else "failed",
                         assertions = sum(d$nb), failures = sum(d$failed),
                         warnings = sum(d$warning), skips = sum(d$skipped)),
                     file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
stopifnot(ok)
cat("SOURCE_GATE_PASS", sha, "assertions", sum(d$nb), "\n")
