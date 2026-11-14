lib <- "/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/mi-matrix/neither"
.libPaths(c(lib, .libPaths()))
stopifnot(startsWith(find.package("pigauto"), lib))
messages <- character()
run_example <- function() {
  messages <<- character()
  capture.output(withCallingHandlers(
    utils::example("compare_methods", package = "pigauto", lib.loc = lib,
                   run.donttest = TRUE, ask = FALSE),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  stopifnot(any(grepl("Skipping GNN comparison: install torch and its Lantern runtime.", messages, fixed = TRUE)))
}
# The actual installed example must skip before calling torch_tensor.
testthat::with_mocked_bindings(run_example(),
  torch_is_installed = function(...) FALSE,
  torch_tensor = function(...) stop("GUARD_FAILURE: tensor called while runtime absent"),
  .package = "torch")
cat("ABSENT_RUNTIME_INDICATOR_SKIP_PASS\n")
# A broken runtime must be caught by the actual tensor probe guard.
testthat::with_mocked_bindings(run_example(),
  torch_is_installed = function(...) TRUE,
  torch_tensor = function(...) stop("Lantern is not loaded. simulated exact failure"),
  .package = "torch")
cat("BROKEN_LANTERN_PROBE_SKIP_PASS\n")
cat("EXACT_INSTALLED_COMPARE_METHODS_DONTTEST_GUARD_PASS\n")
