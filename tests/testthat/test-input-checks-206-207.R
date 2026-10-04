# Issues #206 / #207: bad scalars and contradictory trait_types fail before a fit.

tiny_input <- function() {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:1);")
  df <- data.frame(
    SVL = c(10.1, 12.2, NA_real_),
    row.names = tree$tip.label
  )
  list(tree = tree, df = df)
}

test_that("n_imputations and epochs must be positive integers before a fit", {
  x <- tiny_input()
  expect_error(impute(x$df, x$tree, n_imputations = 0, verbose = FALSE),
               "positive integer")
  expect_error(impute(x$df, x$tree, n_imputations = -2, verbose = FALSE),
               "positive integer")
  expect_error(impute(x$df, x$tree, n_imputations = "a", verbose = FALSE),
               "positive integer")
  expect_error(impute(x$df, x$tree, epochs = 0, verbose = FALSE),
               "positive integer")
  expect_error(impute(x$df, x$tree, epochs = -2, verbose = FALSE),
               "positive integer")
  expect_error(impute(x$df, x$tree, epochs = "a", verbose = FALSE),
               "positive integer")
})

test_that("missing_frac must be in (0, 1) before a fit", {
  x <- tiny_input()
  expect_error(impute(x$df, x$tree, missing_frac = -0.1, verbose = FALSE),
               "\\(0, 1\\)")
  expect_error(impute(x$df, x$tree, missing_frac = 0, verbose = FALSE),
               "\\(0, 1\\)")
})

test_that("unknown trait_types names error before a fit", {
  x <- tiny_input()
  expect_error(
    impute(x$df, x$tree, trait_types = c(not_a_column = "continuous"),
           verbose = FALSE),
    "not found"
  )
})

test_that("declared type that contradicts the column errors before a fit", {
  x <- tiny_input()
  expect_error(
    impute(x$df, x$tree, trait_types = c(SVL = "binary"), verbose = FALSE),
    "contradict"
  )
})
