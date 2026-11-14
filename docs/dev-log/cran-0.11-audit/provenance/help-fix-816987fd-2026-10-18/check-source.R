repo <- "/Users/z3437171/.codex/worktrees/cran-011-final-artifact/pigauto"
out <- "/private/tmp/pigauto-final-closure-20261017/help-fix"
git <- "/Library/Developer/CommandLineTools/usr/bin/git"
base <- "8e275ca738129c8d608f06ffc407b8141c1a2636"
setwd(repo)
files <- system2(git, c("ls-tree", "-r", "--name-only", base, "--", "R"), stdout = TRUE)
files <- files[endsWith(files, ".R")]
same <- vapply(files, function(path) {
  old <- system2(git, c("show", shQuote(paste0(base, ":", path))), stdout = TRUE)
  identical(parse(text = old, keep.source = FALSE), parse(path, keep.source = FALSE))
}, logical(1))
stopifnot(length(files) == 49L, all(same))
write.csv(data.frame(file = files, executable_AST_identical = same),
          file.path(out, "runtime-AST-comparison.csv"), row.names = FALSE)
for (path in c("man/impute.Rd", "man/pigauto_report.Rd")) {
  rd <- tools::parse_Rd(path)
  writeLines(capture.output(tools::Rd2txt(rd)),
             file.path(out, paste0(basename(path), ".txt")))
}
cat("ALL_49_R_EXECUTABLE_ASTS_UNCHANGED\n")
torch::torch_set_num_threads(1L)
result <- devtools::test(filter = "report-213", reporter = "summary")
df <- as.data.frame(result)
flat <- df[, c("file", "test", "failed", "skipped", "error", "warning", "passed")]
write.csv(flat, file.path(out, "report-tests.csv"), row.names = FALSE)
stopifnot(sum(df$failed) == 0L, !any(df$error))
cat("EXISTING_REPORT_TESTS_PASS passed=", sum(df$passed),
    " warnings=", sum(df$warning), " skipped=", sum(df$skipped), "\n", sep = "")
