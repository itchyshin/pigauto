# Issue #205: species_col survives data-only rows and tree-only tips.

issue205_data <- function() {
  set.seed(20261004)
  tr <- ape::rcoal(20)
  tr$tip.label <- paste0("sp", 1:20)
  x <- as.numeric(ape::rTraitCont(tr))
  df <- data.frame(
    species = tr$tip.label,
    x = x,
    y = 0.5 * x + stats::rnorm(20),
    stringsAsFactors = FALSE
  )
  df$x[5] <- NA
  dfa <- df
  dfa$species[1:3] <- paste0("zz", 1:3)
  dfa$x[1] <- NA_real_
  dfb <- df[1:12, ]
  list(tree = tr, df = df, dfa = dfa, dfb = dfb)
}

test_that("species_col imputes when 3 of 20 data species are missing from the tree", {
  d <- issue205_data()
  res <- impute(d$dfa, d$tree, species_col = "species",
                epochs = 8L, verbose = FALSE, seed = 20261004)
  expect_s3_class(res, "pigauto_result")
  completed <- completed_data(res)
  expect_equal(nrow(completed), 20L)
  expect_true(is.na(completed$x[completed$species == "zz1"]))
  expect_output(print(summary(res)), "not-modeled")
  expect_false(grepl("unresolved cells: 0",
                     paste(capture.output(print(summary(res))), collapse = "\n")) &&
                 sum(is.na(completed$x)) > 0 &&
                 !grepl("not-modeled",
                        paste(capture.output(print(summary(res))), collapse = "\n")))
})

test_that("species_col imputes when 8 of 20 tree tips have no data", {
  d <- issue205_data()
  res <- impute(d$dfb, d$tree, species_col = "species",
                epochs = 8L, verbose = FALSE, seed = 20261004)
  expect_s3_class(res, "pigauto_result")
  expect_equal(nrow(completed_data(res)), 12L)
})

test_that("summary() keeps r_cal = 0 as a legal fallback", {
  d <- issue205_data()
  res <- impute(d$df, d$tree, species_col = "species",
                epochs = 8L, verbose = FALSE, seed = 20261004)
  expect_true(all(res$fit$r_cal >= 0))
  expect_true(isTRUE(any(res$fit$r_cal == 0) || all(res$fit$r_cal > 0)))
})
