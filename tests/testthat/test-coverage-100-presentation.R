# Synthetic fit objects isolate the graphics contract from model fitting.
.p100_fit <- function() structure(list(
  calibrated_gates = c(.1, .2, .4),
  trait_map = list(list(name = "mass", type = "continuous", latent_cols = 1L),
                   list(name = "class", type = "categorical", latent_cols = 2:3)),
  conformal_scores = c(mass = .5, class = .8),
  model_config = list(gnn = TRUE, gate_cap = .8),
  history = data.frame(epoch = 1:3, loss_rec = c(3, 2, 1),
    loss_shrink = c(.3, .2, .1), loss_gate = c(.2, .1, 0),
    val_loss = c(4, 3, 2))), class = "pigauto_fit")

.p100_pred <- function() structure(list(
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


test_that("comparison plots aggregate replicates and sort traits by baseline", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  axes <- list(); pts <- list(); seg <- list(); legends <- list()
  local_mocked_bindings(
    axis = function(side, at, labels, ...) {
      axes[[length(axes) + 1L]] <<- list(at = at, labels = labels)
      graphics::axis(side, at = at, labels = labels, ...)
    },
    points = function(x, y, ...) {
      pts[[length(pts) + 1L]] <<- list(x = x, y = y)
      graphics::points(x, y, ...)
    },
    segments = function(x0, y0, x1, y1, ...) {
      seg[[length(seg) + 1L]] <<- c(x0, y0, x1, y1)
      graphics::segments(x0, y0, x1, y1, ...)
    },
    legend = function(x, legend, ...) {
      legends[[length(legends) + 1L]] <<- legend
      graphics::legend(x, legend = legend, ...)
    }, .package = "pigauto")
  d <- data.frame(trait = c("b", "b", "a", "a", "a", "a", "c"),
    type = "continuous", metric = "rmse",
    method = c("baseline", "pigauto", "baseline", "baseline", "pigauto", "pigauto", "baseline"),
    value = c(3, 4, 1, 3, .5, 1.5, 5))
  expect_invisible(plot_comparison(d))
  expect_equal(axes[[1]]$labels, c("a", "b", "c"))
  expect_equal(pts[[1]]$x, c(2, 3, 5))
  expect_equal(pts[[2]]$x, c(1, 4, NA))
  expect_equal(seg[[1]], c(2, 1, 1, 1))
  expect_equal(legends[[1]][1:2], c("baseline", "pigauto"))
  d2 <- transform(d, metric = "accuracy")
  expect_invisible(plot_comparison(rbind(d, d2)))
  expect_invisible(plot_comparison(d2))
  expect_invisible(plot_comparison(d, metric = "rmse"))
  expect_error(plot_comparison(1), "data.frame")
  expect_error(plot_comparison(d["trait"]), "Missing columns")
  expect_error(plot_comparison(d, methods = c("x", "y")), "No rows match")
  expect_error(plot_comparison(transform(d, metric = "other")), "No rmse or accuracy")
  expect_error(plot_comparison(d, metric = "mae"), "No data for metric")
})

test_that("ggplot history and uncertainty retain sorted values and layers", {
  skip_if_not_installed("ggplot2")
  fit <- .p100_fit()
  fit$history$val_loss[2] <- NA_real_
  p <- plot_history_gg(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$epoch, c(1L, 3L))
  expect_length(p$layers, 2L)
  names(fit$history)[names(fit$history) == "val_loss"] <- "val_rmse"
  expect_equal(plot_history_gg(fit)$data$val_rmse, c(4, 2))
  expect_error(plot_history_gg(list()), "pigauto_fit")
  fit$history <- fit$history[FALSE, ]
  expect_error(plot_history_gg(fit), "No training history")
  fit$model_config$gnn <- FALSE
  expect_error(plot_history_gg(fit), "gnn = FALSE")

  pred <- .p100_pred()
  pred$se <- matrix(c(.1, .2, .3, .1, .2, .3), ncol = 2,
                    dimnames = list(NULL, c("mass", "bin")))
  p <- plot_uncertainty(structure(list(prediction = pred), class = "pigauto_result"),
                        truth = data.frame(mass = c(9, 7, NA)), trait_name = "mass")
  expect_equal(p$data$pred, c(1, 2, 3))
  expect_equal(p$data$lower, c(0, 1, 2))
  expect_equal(p$data$upper, c(2, 3, 4))
  expect_equal(p$layers[[3]]$data$truth, c(7, 9))
  expect_error(plot_uncertainty(pred, trait_name = "absent"), "not found")
  pred$imputed <- as.matrix(pred$imputed["mass"])
  expect_error(plot_uncertainty(pred, trait_name = "absent"), "not found")
  expect_equal(plot_uncertainty(pred, trait_name = "mass")$data$pred, c(1, 2, 3))
})

test_that("missingness respects latent groups and phylogenetic order", {
  skip_if_not_installed("ape")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- structure(list(X_scaled = matrix(c(NA, 1, NA, 1, NA, NA, 1, 1, NA), 3,
    dimnames = list(c("s2", "s1", "s3"), NULL)),
    trait_map = list(list(name = "a", type = "continuous", latent_cols = 1L),
                     list(name = "b", type = "categorical", latent_cols = 2:3)),
    species_names = c("s2", "s1", "s3")), class = "pigauto_data")
  expected <- matrix(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE), 3,
    dimnames = list(c("s2", "s1", "s3"), c("a", "b")))
  expect_equal(plot_missingness(x, "summary"), expected)
  expect_equal(plot_missingness(x, "matrix"), expected)
  expect_equal(plot_missingness(x, "co_missing"), expected)
  tree <- ape::read.tree(text = "((s1:1,s2:1):1,s3:2);")
  expect_equal(plot_missingness(x, "phylo", tree), expected[c(2L, 1L, 3L), , drop = FALSE])
  expect_error(plot_missingness(list()), "pigauto_data")
  expect_error(plot_missingness(x, "phylo"), "phylo object")
  x$species_names[1] <- "absent"
  expect_error(plot_missingness(x, "phylo", tree), "All data species")
})

test_that("report rejects invalid inputs and escapes HTML and JavaScript", {
  expect_error(pigauto_report(list(), title = NA), "title")
  expect_error(pigauto_report(list(), output_path = ""), "output_path")
  expect_error(pigauto_report(list(), open = NA), "open")
  expect_error(pigauto_report(list(), data = list()), "data")
  expect_error(pigauto_report(list(), splits = list()), "splits")
  expect_error(pigauto_report(list(), open = FALSE), "pigauto_fit or pigauto_result")
  expect_equal(pigauto:::.escape_report_html("<&\"'>"), "&lt;&amp;&quot;&#39;&gt;")
  expect_equal(pigauto:::.report_js_string("</script>😀"), '"\\u003c\\u002f\\u0073\\u0063\\u0072\\u0069\\u0070\\u0074\\u003e\\ud83d\\ude00"')
})

test_that("report compares paired discrete accuracies and labels the result", {
  tm <- list(name = "class", type = "binary", latent_cols = 1L)
  fit <- list(species_names = c("s1", "s2"), trait_map = list(tm),
              model_config = list(gnn = TRUE, k_eigen = 1L), val_rmse = 0, test_rmse = 0)
  m <- list(name = "class", type = "binary", n_test = 2L, n_bm = 2L,
            n_gnn = 2L, n_paired = 2L, accuracy_bm = .5, accuracy_gnn = 1,
            accuracy_paired_bm = .5, accuracy_paired_gnn = 1)
  html <- pigauto:::build_report_html("report", fit, list(), list(m), fit$trait_map)
  expect_match(html, "+50.0 pp", fixed = TRUE)
  expect_match(html, "Blend better", fixed = TRUE)
  expect_match(html, "paired n=2", fixed = TRUE)
  m$accuracy_paired_gnn <- 0
  expect_match(pigauto:::build_report_html("report", fit, list(), list(m), fit$trait_map), "Baseline better", fixed = TRUE)
  m$accuracy_paired_gnn <- .5
  expect_match(pigauto:::build_report_html("report", fit, list(), list(m), fit$trait_map), "tied", fixed = TRUE)
  fit$model_config$gnn <- FALSE
  expect_match(pigauto:::build_report_html("report", fit, list(), list(m), fit$trait_map), "GNN off", fixed = TRUE)
})

test_that("report validates multi-observation dimensions and split metadata", {
  fit <- list(multi_obs = TRUE, obs_to_species = c(1L, 1L, 2L),
              species_names = c("a", "b"), trait_names = "mass",
              trait_map = list(list(name = "mass", type = "continuous", latent_cols = 1L)))
  expect_equal(pigauto:::.report_latent_dims(fit), c(n = 3L, p = 1L))
  data <- structure(list(species_names = fit$species_names, trait_names = fit$trait_names,
                         X_scaled = matrix(1:3, 3)), class = "pigauto_data")
  expect_invisible(pigauto:::.check_report_data_matches_fit(data, fit))
  valid <- list(val_idx = 1L, test_idx = 3L, n = 3L, p = 1L)
  expect_invisible(pigauto:::.check_report_splits_match_fit(valid, fit, data))
  invalid <- valid; invalid$n <- 2L
  expect_error(pigauto:::.check_report_splits_match_fit(invalid, fit, data), "do not match")
  invalid <- valid; invalid$p <- 2L
  expect_error(pigauto:::.check_report_splits_match_fit(invalid, fit, data), "do not match")
  invalid <- valid; invalid$test_idx <- 4L
  expect_error(pigauto:::.check_report_splits_match_fit(invalid, fit, data), "do not match")
})

test_that("report computes proportion coverage and binary Brier score", {
  tm <- list(name = "trait", type = "proportion", latent_cols = 1L,
             mean = 0, sd = 1, log_transform = FALSE)
  fit <- list(trait_map = list(tm), baseline = list(mu = matrix(c(0, 1), 2)),
              species_names = c("a", "b"), multi_obs = FALSE)
  data <- list(X_scaled = matrix(c(0, 1), 2))
  pred <- list(imputed_latent = matrix(c(0, 1), 2),
               conformal_lower = matrix(c(.4, .7), 2, dimnames = list(NULL, "trait")),
               conformal_upper = matrix(c(.6, .8), 2, dimnames = list(NULL, "trait")))
  splits <- list(test_idx = 1:2)
  m <- pigauto:::collect_report_metrics(fit, pred, data, splits)[[1]]
  expect_equal(m$n_coverage, 2L)
  expect_equal(m$conformal_coverage, 1)
  fit$trait_map[[1]]$type <- "binary"
  data$X_scaled <- matrix(c(0, 1), 2)
  pred <- list(probabilities = list(trait = c(.2, .6)))
  m <- pigauto:::collect_report_metrics(fit, pred, data, splits)[[1]]
  expect_equal(m$brier_gnn, .1)
  expect_equal(m$n_paired, 2L)
})

test_that("plots handle absent matrix truth and tied method values", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  pred <- .p100_pred()
  truth <- matrix(c("a", "b", "c"), dimnames = list(NULL, "other"))
  expect_invisible(plot(pred, type = "probabilities", trait = "bin", data = truth))
  expect_invisible(plot(pred, type = "probabilities", trait = "cat", data = truth))
  d <- data.frame(trait = "a", type = "continuous", metric = "rmse",
                  method = c("baseline", "pigauto"), value = c(1, 1))
  expect_invisible(plot_comparison(d))
  expect_error(plot(structure(list(), class = "pigauto_result")), "fit component")
})

test_that("report requires trait metadata on legacy results", {
  result <- structure(list(fit = list(), prediction = list()), class = "pigauto_result")
  expect_error(pigauto_report(result, open = FALSE, output_path = tempfile()), "requires a trait_map")
  expect_null(pigauto:::.check_report_data_matches_fit(NULL, list()))
  expect_null(pigauto:::.check_report_splits_match_fit(NULL, list()))
})

test_that("report open requests the written report once", {
  opened <- character()
  local_mocked_bindings(browseURL = function(url, ...) {
    opened <<- c(opened, url)
    invisible(NULL)
  }, .package = "utils")
  tm <- list(name = "mass", type = "continuous", latent_cols = 1L)
  fit <- list(species_names = c("s1", "s2"), trait_map = list(tm),
              model_config = list(gnn = FALSE, k_eigen = 1L),
              val_rmse = 0, test_rmse = 0)
  result <- structure(list(fit = fit, prediction = list()), class = "pigauto_result")
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  expect_invisible(pigauto_report(result, output_path = path, title = "Open test", open = TRUE))
  expect_true(file.exists(path))
  expect_match(paste(readLines(path, warn = FALSE), collapse = "\n"), "Open test", fixed = TRUE)
  expect_identical(opened, path)
})

test_that("unnamed category probabilities receive stable labels and ordering", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  heights <- NULL
  labels <- NULL
  testthat::local_mocked_bindings(
    barplot = function(height, ...) {
      heights <<- height
      graphics::barplot(height, ...)
    },
    legend = function(x, legend, ...) {
      labels <<- legend
      graphics::legend(x, legend = legend, ...)
    }, .package = "pigauto")
  probabilities <- rbind(c(.8, .2), c(.1, .9), c(.6, .4))
  .plot_categorical_probs("class", probabilities, NULL, NULL)
  expect_identical(labels, c("C2", "C1"))
  expect_equal(heights, t(probabilities[c(1, 3, 2), ]))
  expect_equal(colSums(heights), rep(1, 3))
})

test_that("uncertainty plots reject absent dataframe traits", {
  pred <- .p100_pred()
  expect_error(plot_uncertainty(pred, trait_name = "absent"),
               "'absent' not found in pred_result\\$imputed")
})

test_that("continuous uncertainty fallback uses matching SEs and matrix truth", {
  pred <- .p100_pred()
  pred$conformal_lower <- pred$conformal_upper <- NULL
  pred$se <- matrix(c(.1, .2, .3), 3, dimnames = list(NULL, "mass"))
  truth <- matrix(c(30, 10, NA), 3, dimnames = list(NULL, "mass"))
  p <- plot_uncertainty(pred, truth = truth, trait_name = "mass")
  expect_equal(p$data$pred, c(1, 2, 3))
  expect_equal(p$data$lower, c(1 - 1.96 * .2, 2 - 1.96 * .3, 3 - 1.96 * .1))
  expect_equal(p$data$upper, c(1 + 1.96 * .2, 2 + 1.96 * .3, 3 + 1.96 * .1))
  expect_equal(p$layers[[3]]$data$truth, c(10, 30))
})

test_that("uncertainty plots omit truth overlays for an unrelated matrix column", {
  pred <- .p100_pred()
  truth <- matrix(c(30, 10, 20), 3, dimnames = list(NULL, "other_trait"))
  p <- plot_uncertainty(pred, truth = truth, trait_name = "mass")
  expect_length(p$layers, 2L)
  expect_equal(p$data$pred, c(1, 2, 3))
  expect_equal(p$data$lower, c(0, 1, 2))
  expect_equal(p$data$upper, c(2, 3, 4))
})
