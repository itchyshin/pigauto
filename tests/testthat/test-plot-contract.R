# Synthetic fit objects isolate the graphics contract from model fitting.
.pc_fit <- function() structure(list(
  calibrated_gates = c(.1, .2, .4),
  trait_map = list(list(name = "mass", type = "continuous", latent_cols = 1L),
                   list(name = "class", type = "categorical", latent_cols = 2:3)),
  conformal_scores = c(mass = .5, class = .8),
  model_config = list(gnn = TRUE, gate_cap = .8),
  history = data.frame(epoch = 1:3, loss_rec = c(3, 2, 1),
    loss_shrink = c(.3, .2, .1), loss_gate = c(.2, .1, 0),
    val_loss = c(4, 3, 2))), class = "pigauto_fit")

.pc_pred <- function() structure(list(
  trait_map = list(list(name = "mass", type = "continuous"),
                   list(name = "bin", type = "binary", levels = c("no", "yes")),
                   list(name = "cat", type = "categorical", levels = c("a", "b", "c"))),
  trait_names = c("mass", "bin", "cat"),
  imputed = data.frame(mass = c(3, 1, 2), bin = c("yes", "no", "yes"),
                       cat = c("a", "b", "c")),
  conformal_lower = matrix(c(2, 0, 1), dimnames = list(NULL, "mass")),
  conformal_upper = matrix(c(4, 2, 3), dimnames = list(NULL, "mass")),
  probabilities = list(bin = c(.8, .1, .6),
    cat = matrix(c(.8, .1, .1, .2, .7, .1, .1, .2, .7), 3, byrow = TRUE,
      dimnames = list(NULL, c("a", "b", "c"))))), class = "pigauto_pred")

test_that("gate plots average latent gates and restore graphics parameters", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  before <- graphics::par(c("mar", "las"))
  heights <- list()
  testthat::local_mocked_bindings(barplot = function(height, ...) {
    heights[[length(heights) + 1L]] <<- height
    graphics::barplot(height, ...)
  }, .package = "pigauto")
  expect_invisible(plot(.pc_fit(), type = "gates"))
  expect_equal(heights[[1]], c(.1, .3))
  expect_equal(graphics::par(c("mar", "las")), before)
  fit <- .pc_fit()
  fit$trait_map <- NULL
  fit$latent_names <- c("x", "y", "z")
  plot(fit, type = "gates")
  expect_equal(heights[[2]], fit$calibrated_gates)
  fit$history <- NULL
  fit$model_config$gnn <- FALSE
  plot(fit)
  expect_equal(heights[[3]], fit$calibrated_gates)
  fit$calibrated_gates <- NULL
  expect_error(plot(fit, type = "gates"), "No calibrated gates")
})

test_that("history plots use finite epochs and support the legacy validation column", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  calls <- list()
  testthat::local_mocked_bindings(lines = function(x, y, ...) {
    calls[[length(calls) + 1L]] <<- list(x = x, y = y)
    graphics::lines(x, y, ...)
  }, .package = "pigauto")
  fit <- .pc_fit()
  fit$history$loss_shrink[] <- NA_real_
  fit$history$val_loss[2] <- NA_real_
  plot(fit)
  expect_length(calls, 3L)
  expect_equal(calls[[3]], list(x = c(1L, 3L), y = c(4, 2)))
  names(fit$history)[names(fit$history) == "val_loss"] <- "val_rmse"
  plot(structure(list(fit = fit), class = "pigauto_result"))
  expect_equal(calls[[6]], calls[[3]])
  fit$history <- NULL
  expect_error(plot(fit, type = "history"), "No training history")
  fit$model_config$gnn <- FALSE
  expect_error(plot(fit, type = "history"), "gnn = FALSE")
  expect_error(plot.pigauto_result(list()), "must be a pigauto_result")
})

test_that("conformal plots omit missing scores and reject absent calibration", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  heights <- NULL
  testthat::local_mocked_bindings(barplot = function(height, ...) {
    heights <<- height
    graphics::barplot(height, ...)
  }, .package = "pigauto")
  fit <- .pc_fit()
  fit$conformal_scores <- c(mass = .5, unknown = .7, class = NA_real_)
  plot(fit, type = "conformal")
  expect_equal(heights, c(mass = .5, unknown = .7))
  fit$trait_map <- NULL
  plot(fit, type = "conformal")
  expect_equal(heights, c(mass = .5, unknown = .7))
  fit$conformal_scores[] <- NA_real_
  expect_error(plot(fit, type = "conformal"), "All conformal scores are NA")
  fit$conformal_scores <- NULL
  expect_error(plot(fit, type = "conformal"), "No conformal scores")
})

test_that("scatter plots pair finite truth by position and filter trait types", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  paired <- list()
  testthat::local_mocked_bindings(points = function(x, y, ...) {
    paired[[length(paired) + 1L]] <<- list(x = x, y = y)
    graphics::points(x, y, ...)
  }, .package = "pigauto")
  pred <- .pc_pred()
  truth <- data.frame(mass = c(4, NA, 1))
  plot(pred, data = truth)
  expect_equal(paired[[1]], list(x = c(4, 1), y = c(3, 2)))
  plot(pred, data = matrix(1:3, dimnames = list(NULL, "mass")))
  expect_equal(paired[[2]], list(x = as.numeric(1:3), y = c(3, 1, 2)))
  plot(pred, data = data.frame(mass = rep(NA_real_, 3)))
  expect_length(paired, 2L)
  expect_error(plot(pred, trait = "bin"), "No continuous, count, or proportion")
  pred$trait_map <- NULL
  pred$trait_names <- "mass"
  pred$imputed <- as.matrix(pred$imputed["mass"])
  plot(pred, trait = "mass")
  expect_error(plot(pred, trait = "absent"), "No traits available")
})

test_that("interval plots sort bounds and observed points together", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ribbon <- overlay <- NULL
  testthat::local_mocked_bindings(
    polygon = function(x, y, ...) {
      ribbon <<- list(x = x, y = y)
      graphics::polygon(x, y, ...)
    },
    points = function(x, y, ...) {
      overlay <<- list(x = x, y = y)
      graphics::points(x, y, ...)
    }, .package = "pigauto")
  pred <- .pc_pred()
  plot(pred, type = "intervals", data = data.frame(mass = c(9, 7, NA)))
  expect_equal(ribbon, list(x = c(1L, 2L, 3L, 3L, 2L, 1L),
                           y = c(0, 1, 2, 4, 3, 2)))
  expect_equal(overlay, list(x = c(1L, 3L), y = c(7, 9)))
  expect_error(plot(pred, type = "intervals", trait = "bin"), "No traits with interval")
  pred$trait_map <- NULL
  pred$trait_names <- "mass"
  pred$conformal_lower <- pred$conformal_upper <- NULL
  pred$imputed <- matrix(2, 3, dimnames = list(NULL, "mass"))
  plot(pred, type = "intervals", data = matrix(1:3, dimnames = list(NULL, "mass")))
  expect_equal(overlay, list(x = 1:3, y = as.numeric(1:3)))
  plot(pred, type = "intervals", data = matrix(1:3, dimnames = list(NULL, "other")))
})

test_that("categorical plots average usable labels and order unlabelled predictions", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  heights <- list()
  testthat::local_mocked_bindings(barplot = function(height, ...) {
    heights[[length(heights) + 1L]] <<- height
    graphics::barplot(height, ...)
  }, .package = "pigauto")
  pred <- .pc_pred()
  truth <- data.frame(cat = c("a", "a", "unknown"))
  plot(pred, type = "probabilities", trait = "cat", data = truth)
  expected <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  expected[, "a"] <- colMeans(pred$probabilities$cat[1:2, ])
  expect_equal(heights[[1]], expected)
  plot(pred, type = "probabilities", trait = "cat", data = as.matrix(truth))
  expect_equal(heights[[2]], expected)
  plot(pred, type = "probabilities", trait = "cat")
  expect_equal(heights[[3]], t(pred$probabilities$cat))
  p1 <- .6 + seq_len(120) / 1000
  pred$probabilities$cat <- cbind(p1, 1 - p1, 0, 0, 0, 0)
  plot(pred, type = "probabilities", trait = "cat")
  expect_equal(dim(heights[[4]]), c(6L, 100L))
  chosen <- 120L - as.integer(round(seq(0, 119, length.out = 100)))
  expect_equal(unname(heights[[4]]), unname(t(pred$probabilities$cat[chosen, ])))
  pred$probabilities$cat <- matrix(1, 3, 1)
  plot(pred, type = "probabilities", trait = "cat")
  expect_length(heights, 4L)
})

test_that("binary plots omit unusable labels and fall back to probability histograms", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  histograms <- list()
  boxes <- list()
  testthat::local_mocked_bindings(hist = function(x, ...) {
    histograms[[length(histograms) + 1L]] <<- x
    graphics::hist(x, ...)
  }, boxplot = function(formula, ...) {
    boxes[[length(boxes) + 1L]] <<- stats::model.frame(formula)
    graphics::boxplot(formula, ...)
  }, .package = "pigauto")
  pred <- .pc_pred()
  plot(pred, type = "probabilities", trait = "bin",
       data = data.frame(bin = c("yes", "no", "invalid")))
  expect_length(histograms, 0L)
  expect_equal(unname(boxes[[1]][[1]]), c(.8, .1))
  expect_equal(as.character(boxes[[1]][[2]]), c("yes", "no"))
  plot(pred, type = "probabilities", trait = "bin",
       data = matrix(c("yes", "no", NA), dimnames = list(NULL, "bin")))
  expect_length(histograms, 0L)
  expect_equal(unname(boxes[[2]][[1]]), c(.8, .1))
  expect_equal(as.character(boxes[[2]][[2]]), c("yes", "no"))
  plot(pred, type = "probabilities", trait = "bin", data = data.frame(bin = rep(NA, 3)))
  expect_equal(histograms[[1]], pred$probabilities$bin)
  pred$trait_map <- NULL
  plot(pred, type = "probabilities", trait = "bin")
  expect_equal(histograms[[2]], pred$probabilities$bin)
  expect_error(plot(pred, type = "probabilities", trait = "absent"), "No matching probability")
  pred$probabilities <- NULL
  expect_error(plot(pred, type = "probabilities"), "No probability predictions")
})
