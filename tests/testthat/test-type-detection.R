# P1-13: integer -> count auto-detection is lossy and, for negative integers,
# actively wrong (count encoding is log1p, so negatives produce NaN).

test_that("[P1-13] negative integer column is NOT auto-detected as count", {
  df <- data.frame(a = c(-3L, 4L, 7L, -1L))
  types <- suppressMessages(pigauto:::detect_trait_types(df))
  expect_equal(unname(types[["a"]]), "continuous")
})

test_that("[P1-13] the negative-integer correction prevents a NaN encoding", {
  # Demonstrates the defect the guard closes: count encoding is log1p, and
  # log1p() of a negative count is NaN.
  expect_true(is.nan(suppressWarnings(log1p(-3))))
  df <- data.frame(a = c(-3L, 4L, 7L, -1L))
  types <- suppressMessages(pigauto:::detect_trait_types(df))
  expect_false(identical(unname(types[["a"]]), "count"))
})

test_that("[P1-13] non-negative integers still auto-detect as count", {
  df <- data.frame(b = c(0L, 2L, 5L, 9L))
  types <- suppressMessages(pigauto:::detect_trait_types(df))
  expect_equal(unname(types[["b"]]), "count")
})

test_that("[P1-13] count auto-detection announces itself (message, not warning)", {
  df <- data.frame(b = c(0L, 2L, 5L))
  expect_message(pigauto:::detect_trait_types(df), "Auto-detected trait type")
  # Must not be a warning: genuine counts are common and would cry wolf.
  expect_silent(suppressMessages(pigauto:::detect_trait_types(df)))
})

test_that("[P1-13] an explicit override wins and silences the notice", {
  df <- data.frame(b = c(0L, 2L, 5L))
  types <- pigauto:::detect_trait_types(df, overrides = c(b = "continuous"))
  expect_equal(unname(types[["b"]]), "continuous")
  expect_silent(pigauto:::detect_trait_types(df, overrides = c(b = "continuous")))
})

test_that("[P1-13] NA-only and mixed-NA integer columns are handled", {
  # any(x < 0, na.rm = TRUE) must not error or return NA on all-NA input.
  expect_equal(
    unname(suppressMessages(
      pigauto:::detect_trait_types(data.frame(c = c(NA_integer_, NA_integer_)))
    )[["c"]]),
    "count")
  expect_equal(
    unname(suppressMessages(
      pigauto:::detect_trait_types(data.frame(d = c(NA_integer_, -2L, 5L)))
    )[["d"]]),
    "continuous")
})
