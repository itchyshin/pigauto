# Issues #210 / #214: Inf covariates and ordinal-on-non-ordered columns
# fail before a fit. Trees stay tiny; no GNN training.

tiny_cov_input <- function() {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  df <- data.frame(
    y = c(1.1, NA_real_, 3.3, 4.4),
    row.names = tree$tip.label
  )
  cv <- data.frame(
    body_height = c(0.2, 0.3, 0.4, 0.5),
    row.names = tree$tip.label
  )
  list(tree = tree, df = df, cv = cv)
}

tiny_ordinal_input <- function() {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  df <- data.frame(
    cont = c(1.1, 2.2, NA_real_, 4.4),
    lik = c(1L, 2L, 3L, 5L),
    row.names = tree$tip.label
  )
  df$lik[3L] <- NA_integer_
  list(tree = tree, df = df)
}

test_that("Inf in a covariate errors before a fit and names the column", {
  x <- tiny_cov_input()
  x$cv$body_height[2L] <- Inf
  expect_error(
    impute(x$df, x$tree, covariates = x$cv, epochs = 2L, verbose = FALSE),
    "body_height"
  )
  expect_error(
    impute(x$df, x$tree, covariates = x$cv, epochs = 2L, verbose = FALSE),
    "Inf|non-finite"
  )
})

test_that("check_pigauto() does not pass Inf in a covariate", {
  x <- tiny_cov_input()
  x$cv$body_height[2L] <- Inf
  chk <- check_pigauto(x$df, x$tree, covariates = x$cv)
  expect_false(isTRUE(chk$ready))
  expect_identical(chk$status, "error")
  expect_true(any(grepl("body_height", chk$messages$text)))
  expect_true(any(grepl("Inf|non-finite", chk$messages$text)))
})

test_that("integer declared ordinal errors before a fit", {
  x <- tiny_ordinal_input()
  expect_error(
    impute(x$df, x$tree, trait_types = c(lik = "ordinal"),
           epochs = 2L, verbose = FALSE),
    "ordinal"
  )
  chk <- check_pigauto(x$df, x$tree, trait_types = c(lik = "ordinal"))
  expect_false(isTRUE(chk$ready))
  expect_identical(chk$status, "error")
})

test_that("numeric declared ordinal errors before a fit", {
  x <- tiny_ordinal_input()
  x$df$lik <- as.numeric(x$df$lik)
  expect_error(
    impute(x$df, x$tree, trait_types = c(lik = "ordinal"),
           epochs = 2L, verbose = FALSE),
    "ordinal"
  )
  chk <- check_pigauto(x$df, x$tree, trait_types = c(lik = "ordinal"))
  expect_false(isTRUE(chk$ready))
  expect_identical(chk$status, "error")
})

test_that("2-level integer declared ordinal errors before a fit", {
  x <- tiny_ordinal_input()
  x$df$lik <- ifelse(is.na(x$df$lik), NA_integer_, ifelse(x$df$lik > 2L, 2L, 1L))
  expect_error(
    impute(x$df, x$tree, trait_types = c(lik = "ordinal"),
           epochs = 2L, verbose = FALSE),
    "ordinal"
  )
  chk <- check_pigauto(x$df, x$tree, trait_types = c(lik = "ordinal"))
  expect_false(isTRUE(chk$ready))
  expect_identical(chk$status, "error")
})

test_that("ordered factor declared ordinal still preprocesses", {
  x <- tiny_ordinal_input()
  x$df$lik <- factor(x$df$lik, levels = 1:5, ordered = TRUE)
  expect_no_error(
    preprocess_traits(x$df, x$tree, trait_types = c(lik = "ordinal"))
  )
  chk <- check_pigauto(x$df, x$tree, trait_types = c(lik = "ordinal"))
  expect_true(chk$status %in% c("ready", "ready_with_warnings"))
})
