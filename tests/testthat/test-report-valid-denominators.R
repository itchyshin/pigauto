report_fixture <- function(type, truth, predicted, baseline, lower = NULL, upper = NULL) {
  n <- nrow(truth)
  tm <- list(name = "trait", type = type, latent_cols = seq_len(ncol(truth)),
             mean = 0, sd = 1, log_transform = FALSE)
  fit <- list(trait_map = list(tm), baseline = list(mu = baseline), multi_obs = FALSE,
              species_names = paste0("s", seq_len(n)), model_config = list(gnn = FALSE, k_eigen = 1L),
              val_rmse = 0, test_rmse = 0)
  if (type %in% c("continuous", "ordinal", "count", "proportion")) {
    pred <- list(imputed_latent = predicted, conformal_lower = lower, conformal_upper = upper)
  } else {
    pred <- list(probabilities = list(trait = predicted))
  }
  list(fit = fit, pred = pred, data = list(X_scaled = truth),
       splits = list(test_idx = seq_len(n)))
}
report_metrics <- function(x) pigauto:::collect_report_metrics(x$fit, x$pred, x$data, x$splits)[[1]]

test_that("report scores finite cells and compares matched cells", {
  x <- report_fixture("continuous", matrix(c(0, 1, 2), ncol = 1),
    matrix(c(1, NA, 12), ncol = 1), matrix(c(2, 1, NA), ncol = 1),
    matrix(c(-1, NA, 1), ncol = 1, dimnames = list(NULL, "trait")),
    matrix(c(1, NA, 3), ncol = 1, dimnames = list(NULL, "trait")))
  m <- report_metrics(x)
  expect_equal(c(m$n_test, m$n_bm, m$n_gnn, m$n_paired, m$n_coverage), c(3, 2, 2, 1, 2))
  expect_equal(c(m$rmse_paired_bm, m$rmse_paired_gnn), c(2, 1))
  html <- pigauto:::build_report_html("report", x$fit, x$pred, list(trait = m), x$fit$trait_map)
  expect_match(html, "+50.0%", fixed = TRUE)
  expect_match(html, "paired n=1", fixed = TRUE)
  expect_match(html, "coverage n=2", fixed = TRUE)
  expect_match(html, "root mean squared", fixed = TRUE)
})

test_that("missing continuous and binary scores do not break the report", {
  for (type in c("continuous", "binary")) {
    x <- report_fixture(type, matrix(c(0, 1), ncol = 1),
      if (type == "continuous") matrix(NA_real_, 2, 1) else c(NA_real_, NA_real_),
      matrix(NA_real_, 2, 1))
    m <- report_metrics(x)
    expect_equal(c(m$n_bm, m$n_gnn, m$n_paired), c(0, 0, 0))
    expect_silent(pigauto:::build_report_html("report", x$fit, x$pred, list(trait = m), x$fit$trait_map))
  }
})

test_that("categorical scoring preserves rows when scores are missing", {
  truth <- rbind(c(1, 0), c(0, 1), c(1, 0))
  p <- rbind(c(0.8, 0.2), c(NA, NA), c(0.3, 0.7))
  b <- rbind(c(0, -1), c(-1, 0), c(NA, NA))
  x <- report_fixture("categorical", truth, p, b)
  m <- report_metrics(x)
  expect_equal(c(m$n_bm, m$n_gnn, m$n_paired), c(2, 2, 1))
  expect_equal(c(m$accuracy_paired_bm, m$accuracy_paired_gnn), c(1, 1))
  x$pred$probabilities$trait[,] <- NA_real_; x$fit$baseline$mu[,] <- NA_real_
  m <- report_metrics(x)
  expect_equal(c(m$n_bm, m$n_gnn, m$n_paired), c(0, 0, 0))
  expect_silent(pigauto:::build_report_html("report", x$fit, x$pred, list(trait = m), x$fit$trait_map))
})


test_that("nonfinite and reversed intervals have no scored coverage", {
  x <- report_fixture("continuous", matrix(c(0, 1), ncol = 1),
    matrix(c(0, Inf), ncol = 1), matrix(c(-Inf, 1), ncol = 1),
    matrix(c(2, -Inf), ncol = 1, dimnames = list(NULL, "trait")),
    matrix(c(1, Inf), ncol = 1, dimnames = list(NULL, "trait")))
  m <- report_metrics(x)
  expect_equal(c(m$n_bm, m$n_gnn, m$n_paired, m$n_coverage), c(1, 1, 0, 0))
  expect_null(m$conformal_coverage)
  x$pred$imputed_latent <- NULL
  expect_equal(report_metrics(x)$n_gnn, 0)
})

test_that("report escapes trait names in HTML and chart JavaScript", {
  x <- report_fixture("continuous", matrix(c(0, 1), ncol = 1),
                      matrix(c(0, 1), ncol = 1), matrix(c(0, 1), ncol = 1))
  x$fit$trait_map[[1]]$name <- "a'</script><script>bad()</script>"
  x$fit$model_config$gnn <- TRUE
  m <- report_metrics(x)
  html <- pigauto:::build_report_html("report", x$fit, x$pred, list(m), x$fit$trait_map)
  expect_false(grepl("<script>bad()</script>", html, fixed = TRUE))
  expect_match(html, "&lt;script&gt;bad()&lt;/script&gt;", fixed = TRUE)
  expect_match(html, "\\u003c", fixed = TRUE)
})

test_that("ordinal and count coverage handles exact class endpoints", {
  for (case in list(list(type = "ordinal", latent = 1.0000000000000002, mean = -1, sd = 1, bound = 0),
                    list(type = "count", latent = log1p(2), mean = 0, sd = 1, bound = 2))) {
    x <- report_fixture(case$type, matrix(case$latent, ncol = 1), matrix(case$latent, ncol = 1),
                        matrix(case$latent, ncol = 1),
                        matrix(case$bound, ncol = 1, dimnames = list(NULL, "trait")),
                        matrix(case$bound, ncol = 1, dimnames = list(NULL, "trait")))
    x$fit$trait_map[[1]]$mean <- case$mean
    x$fit$trait_map[[1]]$sd <- case$sd
    m <- report_metrics(x)
    expect_equal(m$n_coverage, 1)
    expect_equal(m$conformal_coverage, 1)
  }
})

test_that("report names the calibrated blend and does not attribute its gain to GNN", {
  x <- report_fixture("continuous", matrix(c(0, 1), ncol = 1),
                      matrix(c(0.5, 1.5), ncol = 1), matrix(c(1, 2), ncol = 1))
  for (enabled in c(FALSE, TRUE)) {
    x$fit$model_config$gnn <- enabled
    m <- report_metrics(x)
    html <- pigauto:::build_report_html("report", x$fit, x$pred, list(m), x$fit$trait_map)
    expect_match(html, "Calibrated blend", fixed = TRUE)
    expect_false(grepl("GNN better", html, fixed = TRUE))
    expect_match(html, "does not isolate", fixed = TRUE)
    if (!enabled) expect_match(html, "GNN off", fixed = TRUE)
  }
})


test_that("baseline scores expand to repeated observations before masking", {
  x <- report_fixture("continuous", matrix(c(0, 1, 2), ncol = 1),
    matrix(c(0, 1, 2), ncol = 1), matrix(c(0, 2), ncol = 1))
  x$fit$multi_obs <- TRUE
  x$fit$obs_to_species <- c(1L, 1L, 2L)
  x$splits$test_idx <- c(3L, 1L, 2L)
  m <- report_metrics(x)
  expect_equal(c(m$n_bm, m$n_gnn, m$n_paired), c(3, 3, 3))
  expect_equal(m$rmse_bm, sqrt(1 / 3))
})
