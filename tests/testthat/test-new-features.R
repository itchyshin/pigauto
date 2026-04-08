# Tests for evaluate(), summary.pigauto_fit, plot.pigauto_fit,
# pigauto_report(), build_phylo_graph auto k_eigen, and cross_validate()

# ---- Shared helpers --------------------------------------------------------

make_test_data_new <- function(n = 40, p = 2, seed = 100) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df   <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5
  )
  list(tree = tree, df = df)
}

# Build a fit object once for reuse across several tests.
# Kept minimal (20 epochs, 40 tips) for speed.
build_quick_fit <- function(seed = 100) {
  td  <- make_test_data_new(n = 40, seed = seed)
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = seed, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = seed)
  list(fit = fit, pd = pd, spl = spl, tree = td$tree)
}


# ---- 1. evaluate() --------------------------------------------------------

test_that("evaluate returns a data.frame with expected columns", {
  obj <- build_quick_fit(seed = 101)
  eval_df <- evaluate(obj$fit, data = obj$pd, splits = obj$spl)

  expect_s3_class(eval_df, "data.frame")
  expect_true("trait"  %in% names(eval_df))
  expect_true("metric" %in% names(eval_df))
  expect_true("value"  %in% names(eval_df))
  expect_true("method" %in% names(eval_df))
  expect_true("type"   %in% names(eval_df))
  expect_true("n_test" %in% names(eval_df))
  # Should contain pigauto results
  expect_true("pigauto" %in% eval_df$method)
  # Traits should match input
  expect_true(all(c("tr1", "tr2") %in% eval_df$trait))
  # Values should be finite where present
  expect_true(all(is.finite(eval_df$value) | is.na(eval_df$value)))
})

test_that("evaluate errors when data is missing", {
  obj <- build_quick_fit(seed = 102)
  expect_error(evaluate(obj$fit), "data.*required")
})


# ---- 2. summary.pigauto_fit -----------------------------------------------

test_that("summary.pigauto_fit runs without error and prints header", {
  obj <- build_quick_fit(seed = 103)
  # Without data: prints header, skips per-trait metrics
  expect_output(summary(obj$fit), "pigauto_fit")
})

test_that("summary.pigauto_fit with data prints trait metrics", {
  obj <- build_quick_fit(seed = 104)
  # With data: should print trait performance table
  out <- capture.output(summary(obj$fit, data = obj$pd))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("pigauto_fit", combined))
  expect_true(grepl("Species", combined))
  expect_true(grepl("Traits", combined))
})

test_that("summary.pigauto_fit returns evaluation invisibly", {
  obj <- build_quick_fit(seed = 105)
  result <- invisible(capture.output(
    ret <- summary(obj$fit, data = obj$pd)
  ))
  # Should return the evaluation data.frame (or NULL if no test data)
  if (!is.null(ret)) {
    expect_s3_class(ret, "data.frame")
  }
})


# ---- 3. plot.pigauto_fit ---------------------------------------------------

test_that("plot.pigauto_fit history produces plots without error", {
  obj <- build_quick_fit(seed = 106)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_no_error(plot(obj$fit, type = "history"))
  dev.off()
  unlink(tmp)
})

test_that("plot.pigauto_fit gates produces plots without error", {
  obj <- build_quick_fit(seed = 107)
  # gates requires calibrated_gates
  skip_if(is.null(obj$fit$calibrated_gates),
          "No calibrated gates in this fit")
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_no_error(plot(obj$fit, type = "gates"))
  dev.off()
  unlink(tmp)
})

test_that("plot.pigauto_fit conformal produces plots without error", {
  obj <- build_quick_fit(seed = 108)
  skip_if(is.null(obj$fit$conformal_scores) ||
            all(is.na(obj$fit$conformal_scores)),
          "No conformal scores in this fit")
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_no_error(plot(obj$fit, type = "conformal"))
  dev.off()
  unlink(tmp)
})


# ---- 4. pigauto_report generates HTML --------------------------------------

test_that("pigauto_report creates an HTML file from pigauto_result", {
  td <- make_test_data_new(n = 20, seed = 109)
  result <- impute(td$df, td$tree, epochs = 20L, verbose = FALSE, seed = 109)

  expect_s3_class(result, "pigauto_result")

  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)

  pigauto_report(result, output_path = tmp, open = FALSE)
  expect_true(file.exists(tmp))

  html <- readLines(tmp)
  # Should contain pigauto branding
  expect_true(any(grepl("pigauto", html, ignore.case = TRUE)))
  # Should be valid HTML
  expect_true(any(grepl("<html", html)))
  expect_true(any(grepl("</html>", html)))
})

test_that("pigauto_report works with a pigauto_fit object and data", {
  obj <- build_quick_fit(seed = 110)

  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)

  pigauto_report(obj$fit, data = obj$pd, splits = obj$spl,
                 output_path = tmp, open = FALSE)
  expect_true(file.exists(tmp))

  html <- readLines(tmp)
  expect_true(any(grepl("pigauto", html, ignore.case = TRUE)))
})


# ---- 5. adaptive k_eigen in build_phylo_graph ------------------------------

test_that("build_phylo_graph auto k_eigen scales with tree size", {
  set.seed(200)
  g30  <- build_phylo_graph(ape::rtree(30))
  g200 <- build_phylo_graph(ape::rtree(200))

  # Larger tree should get more spectral features

  expect_true(ncol(g200$coords) > ncol(g30$coords))

  # Verify auto formula: min(max(ceiling(n/20), 4), 32)
  expect_equal(ncol(g30$coords),  max(ceiling(30 / 20), 4L))   # 4
  expect_equal(ncol(g200$coords), max(ceiling(200 / 20), 4L))  # 10
})

test_that("build_phylo_graph auto k_eigen is capped at 32", {
  set.seed(201)
  # 700 tips: ceiling(700/20)=35, should cap at 32
  g700 <- build_phylo_graph(ape::rtree(700))
  expect_equal(ncol(g700$coords), 32L)
})

test_that("build_phylo_graph auto k_eigen floors at 4 for small trees", {
  set.seed(202)
  # 10 tips: ceiling(10/20)=1, should floor at 4
  g10 <- build_phylo_graph(ape::rtree(10))
  expect_equal(ncol(g10$coords), 4L)
})


# ---- 6. cross_validate basic test ------------------------------------------

test_that("cross_validate returns pigauto_cv object", {
  # k must be >= 3 so that at least one fold remains as training data
  # (with k=2 the test + val folds mask all cells, leaving nothing for
  # the Rphylopars baseline).  Use 50 species and 3 traits.
  set.seed(301)
  tree <- ape::rtree(50)
  df <- data.frame(
    row.names = tree$tip.label,
    tr1 = abs(stats::rnorm(50)) + 0.5,
    tr2 = abs(stats::rnorm(50)) + 0.5,
    tr3 = abs(stats::rnorm(50)) + 0.5
  )
  pd <- preprocess_traits(df, tree)

  cv <- cross_validate(pd, tree, k = 3L, seeds = 1L,
                       epochs = 20L, eval_every = 10L,
                       patience = 5L, verbose = FALSE)

  expect_s3_class(cv, "pigauto_cv")
  expect_true("results" %in% names(cv))
  expect_true("summary" %in% names(cv))
  expect_true("k" %in% names(cv))
  expect_true("n_reps" %in% names(cv))
  expect_equal(cv$k, 3L)
  expect_equal(cv$n_reps, 1L)
})

test_that("cross_validate results have expected structure", {
  set.seed(302)
  tree <- ape::rtree(50)
  df <- data.frame(
    row.names = tree$tip.label,
    tr1 = abs(stats::rnorm(50)) + 0.5,
    tr2 = abs(stats::rnorm(50)) + 0.5,
    tr3 = abs(stats::rnorm(50)) + 0.5
  )
  pd <- preprocess_traits(df, tree)

  cv <- cross_validate(pd, tree, k = 3L, seeds = 1L,
                       epochs = 20L, eval_every = 10L,
                       patience = 5L, verbose = FALSE)

  # Results data.frame
  res <- cv$results
  expect_s3_class(res, "data.frame")
  if (nrow(res) > 0L) {
    expect_true(all(c("fold", "rep", "trait", "type", "metric", "value")
                    %in% names(res)))
    # All three traits should appear
    expect_true(all(c("tr1", "tr2", "tr3") %in% res$trait))
  }

  # Summary data.frame
  summ <- cv$summary
  expect_s3_class(summ, "data.frame")
  if (nrow(summ) > 0L) {
    expect_true(all(c("trait", "type", "metric", "mean", "sd", "n_folds")
                    %in% names(summ)))
  }
})

test_that("cross_validate print and summary methods work", {
  set.seed(303)
  tree <- ape::rtree(50)
  df <- data.frame(
    row.names = tree$tip.label,
    tr1 = abs(stats::rnorm(50)) + 0.5,
    tr2 = abs(stats::rnorm(50)) + 0.5,
    tr3 = abs(stats::rnorm(50)) + 0.5
  )
  pd <- preprocess_traits(df, tree)

  cv <- cross_validate(pd, tree, k = 3L, seeds = 1L,
                       epochs = 20L, eval_every = 10L,
                       patience = 5L, verbose = FALSE)

  expect_output(print(cv), "fold cross-validation")
  expect_output(summary(cv), "fold cross-validation")
})


# ---- 7. simulate_benchmark --------------------------------------------------

test_that("simulate_benchmark runs and returns expected structure", {
  bench <- simulate_benchmark(
    n_species = 20, n_traits = 2, scenarios = "BM",
    n_reps = 1, epochs = 20L, verbose = FALSE
  )

  expect_s3_class(bench, "pigauto_benchmark")
  expect_true("results" %in% names(bench))
  expect_true("summary" %in% names(bench))
  expect_s3_class(bench$results, "data.frame")
  expect_true(all(c("scenario", "method", "trait", "metric", "value")
                  %in% names(bench$results)))
  expect_equal(bench$n_species, 20L)
  expect_equal(bench$n_reps, 1L)

  # Print and summary work
  expect_output(print(bench), "pigauto simulation benchmark")
  expect_output(summary(bench), "pigauto simulation benchmark")
})

# ---- 8. save_pigauto / load_pigauto ----------------------------------------

test_that("save_pigauto and load_pigauto round-trip correctly", {
  obj <- build_quick_fit(seed = 400)

  tmp <- tempfile(fileext = ".pigauto")
  on.exit(unlink(tmp), add = TRUE)

  save_pigauto(obj$fit, tmp)
  expect_true(file.exists(tmp))

  fit2 <- load_pigauto(tmp)
  expect_s3_class(fit2, "pigauto_fit")
  expect_equal(length(fit2$species_names), length(obj$fit$species_names))
  expect_equal(fit2$model_config$hidden_dim, obj$fit$model_config$hidden_dim)

  # Predictions from loaded model should match
  pred1 <- predict(obj$fit, return_se = FALSE)
  pred2 <- predict(fit2, return_se = FALSE)
  expect_equal(pred1$imputed_latent, pred2$imputed_latent,
               tolerance = 1e-5)
})


test_that("simulate_benchmark mixed scenario handles discrete traits", {
  bench <- simulate_benchmark(
    n_species = 25, scenarios = "mixed",
    n_reps = 1, epochs = 20L, verbose = FALSE
  )

  expect_s3_class(bench, "pigauto_benchmark")
  # Should have categorical and binary traits
  types <- unique(bench$results$type)
  expect_true("binary" %in% types || "categorical" %in% types)
})
