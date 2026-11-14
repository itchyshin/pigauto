# P1-13 (auto-detection) and P1-10 (encode guards): data that contradicts its
# declared or detected type is now refused rather than silently mangled.

make_tree <- function(n, labels) {
  set.seed(1L); tr <- ape::rtree(n); tr$tip.label <- labels; tr
}

# ---- P1-10: fail-closed encode guards ---------------------------------------

test_that("[P1-10] explicit count type errors on negative values", {
  sp <- paste0("s", 1:6)
  df <- data.frame(row.names = sp, a = c(-2L, 1L, 3L, 4L, 5L, 6L))
  tr <- make_tree(6L, sp)
  expect_error(
    preprocess_traits(df, tr, trait_types = c(a = "count")),
    "negative")
})

test_that("[P1-10] proportion out of [0,1] errors instead of silent clamping", {
  sp <- paste0("s", 1:6)
  tr <- make_tree(6L, sp)
  # Percentages mistyped as a proportion — previously squashed to 0.999.
  df <- data.frame(row.names = sp, p = c(10, 25, 50, 75, 90, 99))
  expect_error(
    preprocess_traits(df, tr, trait_types = c(p = "proportion")),
    "outside \\[0, 1\\]")
  # Legitimate boundary values 0 and 1 must still be accepted (clamped).
  df_ok <- data.frame(row.names = sp, p = c(0, 0.25, 0.5, 0.75, 1, 0.5))
  expect_no_error(preprocess_traits(df_ok, tr, trait_types = c(p = "proportion")))
})

test_that("[P1-10] multi_proportion rows must sum to 1", {
  sp <- paste0("s", 1:6)
  tr <- make_tree(6L, sp)
  ok <- data.frame(row.names = sp,
                   a = rep(0.5, 6), b = rep(0.3, 6), c = rep(0.2, 6))
  expect_no_error(preprocess_traits(ok, tr,
    multi_proportion_groups = list(g = c("a", "b", "c"))))

  # Percentage scale (sums to 100) — previously renormalised silently.
  pct <- data.frame(row.names = sp,
                    a = rep(50, 6), b = rep(30, 6), c = rep(20, 6))
  expect_error(preprocess_traits(pct, tr,
    multi_proportion_groups = list(g = c("a", "b", "c"))), "sum to 1")

  # An all-zero row previously became an invented uniform composition.
  z <- ok; z[3, ] <- 0
  expect_error(preprocess_traits(z, tr,
    multi_proportion_groups = list(g = c("a", "b", "c"))), "sum to 1")

  # Negative component.
  neg <- ok; neg[2, "a"] <- -0.1; neg[2, "b"] <- 0.9
  expect_error(preprocess_traits(neg, tr,
    multi_proportion_groups = list(g = c("a", "b", "c"))), "non-negative")
})

test_that("[P1-10] NA rows in a multi_proportion group are still allowed", {
  sp <- paste0("s", 1:6)
  tr <- make_tree(6L, sp)
  df <- data.frame(row.names = sp,
                   a = rep(0.5, 6), b = rep(0.3, 6), c = rep(0.2, 6))
  df[4, ] <- NA_real_   # a genuinely unobserved composition
  expect_no_error(preprocess_traits(df, tr,
    multi_proportion_groups = list(g = c("a", "b", "c"))))
})

# ---- P1-13: integer -> count auto-detection ---------------------------------
# Lossy in general and, for negative integers, actively wrong (count encoding
# is log1p, so negatives produce NaN).

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
