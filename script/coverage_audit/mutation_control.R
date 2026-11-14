# Require the specific existing gate assertion to detect an in-memory defect.
# No package source file is changed.
main <- function() {
  devtools::load_all(quiet = TRUE)
  original <- get(".plot_gates", envir = asNamespace("pigauto"))
  on.exit(assignInNamespace(".plot_gates", original, ns = "pigauto"), add = TRUE)
  baseline <- testthat::test_file("tests/testthat/test-plot-contract.R",
    reporter = "silent", stop_on_failure = FALSE)
  clean <- as.data.frame(baseline)
  stopifnot(nrow(clean) > 0L, !any(clean$failed), !any(clean$error),
            !any(clean$warning), !any(clean$skipped))
  mutated <- original
  code <- paste(deparse(body(original)), collapse = "\n")
  stopifnot(grepl("mean(gates[lc])", code, fixed = TRUE))
  body(mutated) <- parse(text = sub("mean(gates[lc])", "sum(gates[lc])", code, fixed = TRUE))[[1]]
  assignInNamespace(".plot_gates", mutated, ns = "pigauto")
  result <- testthat::test_file("tests/testthat/test-plot-contract.R",
    reporter = "silent", stop_on_failure = FALSE)
  counts <- as.data.frame(result)
  target <- grepl("gate plots average", counts$test, fixed = TRUE)
  stopifnot(sum(target) == 1L, sum(counts$failed[target]) > 0L,
            !any(counts$error), !any(counts$failed[!target]))
  cat("MUTATION_CONTROL_PASS gate mean-to-sum was detected by its gate assertion\n")
}
main()
