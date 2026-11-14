test_that("comparison improvement labels identify the supplied second method", {
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  labels <- list()
  testthat::local_mocked_bindings(legend = function(x, legend, ...) {
    labels[[length(labels) + 1L]] <<- legend
    graphics::legend(x, legend = legend, ...)
  }, .package = "pigauto")
  for (methods in list(c("baseline", "pigauto"), c("reference", "candidate"))) {
    results <- data.frame(trait = rep(c("Mass", "Class"), each = 2),
                          type = rep(c("continuous", "binary"), each = 2),
                          metric = rep(c("rmse", "accuracy"), each = 2),
                          method = rep(methods, 2), value = c(1, .8, .6, .7))
    plot_comparison(results, methods = methods)
    for (legend in tail(labels, 2)) {
      expect_identical(legend, c(methods, paste(methods[2], c("better", "worse"))))
      expect_false(any(grepl("GNN", legend, fixed = TRUE)))
    }
  }
})
