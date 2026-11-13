test_that("read_traits accepts an explicit ecological level order", {
  raw <- data.frame(species = c("a", "b", "c"),
                    status = c("VU", NA, "LC"))
  out <- read_traits(raw, ordinal_levels = list(status = c("LC", "NT", "VU")))
  expect_true(is.ordered(out$status))
  expect_identical(levels(out$status), c("LC", "NT", "VU"))
  expect_identical(as.character(out$status), raw$status)
  expect_error(read_traits(raw, ordinal_levels = list(status = c("LC", "NT"))),
               "not included")
  expect_error(read_traits(raw, ordinal_levels = list(absent = c("a", "b"))),
               "not found")
  expect_error(read_traits(raw, ordinal_levels = list(status = c("LC", "LC"))),
               "unique")
})

test_that("predict rejects unused arguments rather than silently ignoring them", {
  # Argument validation precedes fitted-state validation so typos are clear.
  expect_error(predict.pigauto_fit(list(), evaluaton = TRUE), "Unused.*evaluaton")
})
