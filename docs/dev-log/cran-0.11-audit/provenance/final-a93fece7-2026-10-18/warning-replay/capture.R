args <- commandArgs(TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
lib <- normalizePath(file.path(root, "library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
stopifnot(identical(Sys.getenv("NOT_CRAN", unset = NA_character_), NA_character_))
stopifnot(identical(Sys.getenv("_R_CHECK_FORCE_SUGGESTS_"), "true"))
stopifnot(startsWith(find.package("pigauto"), lib))
stopifnot(as.character(packageVersion("testthat")) == "3.3.2")
torch::torch_set_num_threads(1L)
stopifnot(torch::torch_tensor(2)$item() == 2)
setwd(file.path(root, "warning-replay.Rcheck", "tests"))
stopifnot(!file.exists("../DESCRIPTION"))
options(useFancyQuotes = FALSE)
reporter <- testthat::CheckReporter$new()
start <- proc.time()[["elapsed"]]
result <- tryCatch(testthat::test_check("pigauto", reporter = reporter), error = identity)
elapsed <- proc.time()[["elapsed"]] - start
warnings <- reporter$warnings$as_list()
skips <- reporter$skips$as_list()
problems <- reporter$problems$as_list()
saveRDS(list(warnings = warnings, skips = skips, problems = problems),
        file.path(root, "warning-objects.rds"), version = 2)
warning_text <- capture.output(lapply(seq_along(warnings), function(i) {
  cat("WARNING", i, "\n")
  str(warnings[[i]], max.level = 2)
}))
writeLines(warning_text, file.path(root, "warning-objects.txt"))
messages <- vapply(warnings, function(x) x$message, character(1))
write.csv(data.frame(id = seq_along(messages), message = messages),
          file.path(root, "warning-messages.csv"), row.names = FALSE)
summary <- list(package = as.character(packageVersion("pigauto")),
                package_path = find.package("pigauto"),
                testthat = as.character(packageVersion("testthat")),
                not_cran = "UNSET", force_suggests = Sys.getenv("_R_CHECK_FORCE_SUGGESTS_"),
                pass = reporter$n_ok, fail = length(problems),
                warnings = length(warnings), skips = length(skips),
                elapsed_seconds = elapsed,
                error = if (inherits(result, "error")) conditionMessage(result) else NULL)
jsonlite::write_json(summary, file.path(root, "warning-replay.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(root, "warning-session.txt"))
stopifnot(!inherits(result, "error"), length(problems) == 0L)
cat("EXACT_ARCHIVE_WARNING_CAPTURE_COMPLETE\n")
