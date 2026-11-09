test_that("discrete uncertainty scores are refused as Gaussian SE bands", {
  for (type in c("binary", "categorical")) {
    levels <- if (type == "binary") c("no", "yes") else c("a", "b", "c")
    pred <- structure(list(imputed = data.frame(x = factor(levels)),
      se = matrix(rep(.2, length(levels)), dimnames = list(NULL, "x")),
      trait_map = list(list(name = "x", type = type, levels = levels))),
      class = "pigauto_pred")
    expect_error(plot_uncertainty(pred, trait_name = "x"),
                 "not Gaussian standard errors", fixed = TRUE)
    pred$trait_map <- NULL
    expect_error(plot_uncertainty(pred, trait_name = "x"),
                 "not Gaussian standard errors", fixed = TRUE)
  }
})

test_that("missing analytic SE is diagnosed before matrix indexing", {
  for (se in list(NULL, matrix(1, 3, dimnames = list(NULL, "other")), 1:3)) {
    pred <- list(imputed = data.frame(x = 1:3), se = se)
    expect_error(plot_uncertainty(pred, trait_name = "x"), "No SE available", fixed = TRUE)
  }
  pred <- list(imputed = data.frame(x = c(1, 2, 3)),
    conformal_lower = matrix(c(0, 0, 1), dimnames = list(NULL, "x")),
    conformal_upper = matrix(c(3, 4, 5), dimnames = list(NULL, "x")))
  p <- plot_uncertainty(pred, trait_name = "x")
  expect_equal(p$data$lower, c(0, 0, 1))
  expect_equal(p$data$upper, c(3, 4, 5))
})

test_that("ordered ordinal predictions retain their supported interval plot", {
  pred <- list(imputed = data.frame(x = ordered(c("low", "high", "mid"),
    levels = c("low", "mid", "high"))),
    trait_map = list(list(name = "x", type = "ordinal")),
    conformal_lower = matrix(c(.5, 2.5, 1.5), dimnames = list(NULL, "x")),
    conformal_upper = matrix(c(1.5, 3.5, 2.5), dimnames = list(NULL, "x")))
  p <- plot_uncertainty(pred, trait_name = "x")
  expect_equal(p$data$pred, c(1, 2, 3))
  expect_equal(p$data$lower, c(.5, 1.5, 2.5))
})
