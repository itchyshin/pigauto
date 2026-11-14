# Separate diagnostic replay after the native exact-archive check.
# The original check logs remain the check verdict and timing evidence.
parse_native_summary <- function(lines) {
  pattern <- "\\[ FAIL ([0-9]+) \\| WARN ([0-9]+) \\| SKIP ([0-9]+) \\| PASS ([0-9]+) \\]"
  matches <- regmatches(lines, regexec(pattern, lines))
  matches <- Filter(function(x) length(x) == 5L, matches)
  if (!length(matches)) stop("Native test summary is missing")
  setNames(as.integer(tail(matches, 1L)[[1]][-1]),
           c("fail", "warnings", "skips", "pass"))
}

if (identical(commandArgs(TRUE), "--selftest")) {
  expected <- c(fail = 0L, warnings = 2L, skips = 3L, pass = 4L)
  stopifnot(identical(parse_native_summary(c("noise", "[ FAIL 0 | WARN 2 | SKIP 3 | PASS 4 ]")), expected))
  refused <- tryCatch({parse_native_summary("unfinished"); FALSE}, error = function(e) TRUE)
  stopifnot(refused)
  stopifnot(!identical(parse_native_summary("[ FAIL 1 | WARN 2 | SKIP 3 | PASS 4 ]"), expected))
  cat("WARNING_CAPTURE_SUMMARY_CONTROLS_PASS\n")
} else {
  evidence <- normalizePath("platform-evidence", mustWork = TRUE)
  check_root <- normalizePath("frozen-check/pigauto.Rcheck", mustWork = TRUE)
  native <- parse_native_summary(readLines(file.path(check_root, "tests", "testthat.Rout"), warn = FALSE))
  .libPaths(c(check_root, .libPaths()))
  loaded <- normalizePath(find.package("pigauto"), mustWork = TRUE)
  stopifnot(identical(loaded, normalizePath(file.path(check_root, "pigauto"))))
  stopifnot(identical(Sys.getenv("NOT_CRAN", unset = NA_character_), NA_character_))
  stopifnot(identical(Sys.getenv("_R_CHECK_FORCE_SUGGESTS_"), "true"))
  torch::torch_set_num_threads(1L)
  stopifnot(torch::torch_tensor(2)$item() == 2)
  setwd(file.path(check_root, "tests"))
  stopifnot(!file.exists("../DESCRIPTION"))
  options(useFancyQuotes = FALSE)
  reporter <- testthat::CheckReporter$new()
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(testthat::test_check("pigauto", reporter = reporter), error = identity)
  elapsed <- proc.time()[["elapsed"]] - started
  warnings <- reporter$warnings$as_list()
  skips <- reporter$skips$as_list()
  problems <- reporter$problems$as_list()
  saveRDS(list(warnings = warnings, skips = skips, problems = problems),
          file.path(evidence, "warning-objects.rds"), version = 2)
  details <- lapply(seq_along(warnings), function(i) {
    w <- warnings[[i]]
    src <- attr(w$srcref, "srcfile")
    calls <- w$trace$call
    data.frame(id = i, message = w$message, test = w$test,
               file = basename(src$filename), line = as.integer(w$srcref[[1]]),
               last_call = paste(deparse(tail(calls, 1L)[[1]]), collapse = " "))
  })
  write.csv(do.call(rbind, details), file.path(evidence, "warning-details.csv"), row.names = FALSE)
  observed <- c(fail = length(problems), warnings = length(warnings),
                skips = length(skips), pass = reporter$n_ok)
  summary <- list(package = as.character(packageVersion("pigauto")), package_path = loaded,
                  testthat = as.character(packageVersion("testthat")), not_cran = "UNSET",
                  force_suggests = Sys.getenv("_R_CHECK_FORCE_SUGGESTS_"),
                  native_counts = as.list(native), replay_counts = as.list(observed),
                  diagnostic_elapsed_seconds = elapsed,
                  scope = "Separate native-layout replay; not part of R CMD check timing",
                  error = if (inherits(result, "error")) conditionMessage(result) else NULL)
  jsonlite::write_json(summary, file.path(evidence, "warning-replay.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(evidence, "warning-session.txt"))
  stopifnot(!inherits(result, "error"), length(problems) == 0L,
            identical(unname(as.integer(native)), unname(as.integer(observed))))
  cat("EXACT_ARCHIVE_WARNING_CAPTURE_COMPLETE\n")
}
