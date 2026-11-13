test_that("uncertainty plots use prediction ranks for every aligned layer", {
  pred <- structure(list(imputed = data.frame(x = c(3, 1, 2)),
    conformal_lower = matrix(c(2.5, .5, 1.5), dimnames = list(NULL, "x")),
    conformal_upper = matrix(c(3.5, 1.5, 2.5), dimnames = list(NULL, "x"))),
    class = "pigauto_pred")
  p <- plot_uncertainty(pred, trait_name = "x", truth = data.frame(x = c(30, 10, 20)))
  expect_equal(p$data$i, 1:3)
  expect_equal(p$data$pred, c(1, 2, 3))
  expect_equal(p$data$lower, c(.5, 1.5, 2.5))
  expect_equal(p$data$upper, c(1.5, 2.5, 3.5))
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[2]]$y, c(1, 2, 3))
  expect_equal(built$data[[3]]$x, c(1, 2, 3))
  expect_equal(built$data[[3]]$y, c(10, 20, 30))
})
