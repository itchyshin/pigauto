test_that("all failed benchmark replicates retain the first cause and attempt count", {
  for (verbose in c(FALSE, TRUE)) {
    messages <- character()
    error <- withCallingHandlers(
      tryCatch(simulate_benchmark(n_species = 8L, n_traits = 1L,
        scenarios = c("BM", "OU"), missing_frac = 1, n_reps = 2L,
        seed = 3L, gnn = FALSE, verbose = verbose), error = identity),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error),
                 "No benchmark replicates succeeded \\(0 of 4\\)")
    expect_match(conditionMessage(error), "BM replicate 1.*missing_frac")
    expect_false(grepl("first argument must be a vector", conditionMessage(error),
                        fixed = TRUE))
    if (verbose) {
      expect_true(any(grepl("ERROR:.*missing_frac", messages)))
    } else {
      expect_length(messages, 0L)
    }
  }
})

test_that("a benchmark with no attempted replicates has a direct diagnostic", {
  expect_error(simulate_benchmark(n_species = 8L, n_traits = 1L,
    scenarios = "BM", n_reps = 0L, gnn = FALSE, verbose = FALSE),
    "No benchmark replicates succeeded \\(0 of 0\\).*No fits were run")
})

test_that("successful benchmark aggregation still contains its evaluation rows", {
  expect_warning(result <- simulate_benchmark(n_species = 16L, n_traits = 1L,
    scenarios = "BM", missing_frac = .25, n_reps = 1L, seed = 3L,
    gnn = FALSE, verbose = FALSE), "Small validation set")
  expect_s3_class(result, "pigauto_benchmark")
  expect_gt(nrow(result$results), 0L)
  expect_identical(result$n_reps, 1L)
  expect_identical(unique(result$results$scenario), "BM")
  rmse <- result$summary[result$summary$metric == "rmse", ]
  expect_equal(nrow(rmse), 2L)
  expect_length(unique(rmse$method), 2L)
  expect_true(all(is.finite(rmse$mean)))
  expect_equal(rmse$n_reps, rep(1L, 2L))
})
