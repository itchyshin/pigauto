# plot_comparison() must default to the method labels that
# evaluate() and compare_methods() actually emit (#211).

test_that("plot_comparison default methods match evaluate()/compare_methods()", {
  results <- data.frame(
    trait  = c("Mass", "Mass"),
    type   = c("continuous", "continuous"),
    metric = c("rmse", "rmse"),
    method = c("baseline", "pigauto"),
    value  = c(1.0, 0.8),
    stringsAsFactors = FALSE
  )
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot_comparison(results))
})
