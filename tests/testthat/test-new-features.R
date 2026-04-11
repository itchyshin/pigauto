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
  # the BM baseline).  Use 50 species and 3 traits.
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


# ---- Environmental covariates -----------------------------------------------

test_that("impute() accepts covariates and threads them through the pipeline", {
  n <- 40
  set.seed(500)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label

  df <- data.frame(
    row.names = sp,
    y  = abs(rnorm(n)) + 0.5,
    x1 = abs(rnorm(n)) + 0.5
  )
  df$y[sample(n, 8)] <- NA

  covs <- data.frame(
    temperature   = rnorm(n, 20, 5),
    precipitation = rnorm(n, 1000, 200)
  )

  res <- impute(df, tree, covariates = covs, epochs = 20L,
                verbose = FALSE, seed = 500L,
                eval_every = 10L, patience = 5L)

  expect_s3_class(res, "pigauto_result")
  expect_equal(sum(res$imputed_mask), 8L)

  # Covariates stored in data object
  expect_true(!is.null(res$data$covariates))
  expect_equal(ncol(res$data$covariates), 2L)
  expect_equal(res$data$cov_names, c("temperature", "precipitation"))

  # Covariates stored in fit object (for predict)
  expect_true(!is.null(res$fit$covariates))
  expect_equal(ncol(res$fit$covariates), 2L)

  # Model config has larger cov_dim
  base_cov_dim <- res$data$p_latent + 1L  # baseline + mask_ind
  expect_equal(res$fit$model_config$cov_dim, base_cov_dim + 2L)
})

test_that("covariates with NAs are rejected", {
  set.seed(501)
  tree <- ape::rtree(20)
  sp <- tree$tip.label
  df <- data.frame(row.names = sp, y = rnorm(20))

  covs_bad <- data.frame(temp = c(NA, rnorm(19)))
  expect_error(
    impute(df, tree, covariates = covs_bad, epochs = 10L, verbose = FALSE),
    regexp = "fully observed"
  )
})

test_that("covariates with wrong row count are rejected", {
  set.seed(502)
  tree <- ape::rtree(20)
  sp <- tree$tip.label
  df <- data.frame(row.names = sp, y = rnorm(20))

  covs_bad <- data.frame(temp = rnorm(10))  # wrong n
  expect_error(
    impute(df, tree, covariates = covs_bad, epochs = 10L, verbose = FALSE),
    regexp = "rows"
  )
})


# ---- Multi-obs + covariates (Phase 1 validation) ----------------------------
# These tests validate the current behaviour of species_col + covariates
# used together, establishing a baseline before the Phase 2 architecture
# improvement (observation-level refinement MLP).

test_that("impute() with species_col + covariates runs end-to-end", {
  set.seed(600)
  tree <- ape::rtree(20)
  n_obs <- 60  # 3 obs per species

  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    y1 = abs(rnorm(n_obs)) + 0.5,
    y2 = abs(rnorm(n_obs)) + 0.5
  )
  # Introduce some missingness
  df$y1[sample(n_obs, 12)] <- NA
  df$y2[sample(n_obs, 12)] <- NA

  covs <- data.frame(
    temperature   = rnorm(n_obs, 20, 5),
    precipitation = rnorm(n_obs, 1000, 200)
  )

  result <- impute(df, tree, species_col = "species",
                   covariates = covs,
                   epochs = 20L, verbose = FALSE, seed = 600L,
                   eval_every = 10L, patience = 5L)

  expect_s3_class(result, "pigauto_result")

  # Predictions at observation level (n_obs, not n_species)
  expect_equal(nrow(result$prediction$imputed), n_obs)
  expect_true(result$prediction$multi_obs)

  # Covariates stored correctly
  expect_true(!is.null(result$data$covariates))
  expect_equal(ncol(result$data$covariates), 2L)
  expect_equal(nrow(result$data$covariates), n_obs)
  expect_equal(result$data$cov_names, c("temperature", "precipitation"))

  # Covariates stored in fit for predict()
  expect_true(!is.null(result$fit$covariates))
  expect_equal(ncol(result$fit$covariates), 2L)
  expect_equal(nrow(result$fit$covariates), n_obs)

  # Model config has expanded cov_dim: p_latent + 1 (mask_ind) + 2 (user covs)
  base_cov_dim <- result$data$p_latent + 1L
  expect_equal(result$fit$model_config$cov_dim, base_cov_dim + 2L)

  # All imputed values are finite
  imp <- result$prediction$imputed
  expect_true(all(sapply(imp, function(col) all(is.finite(col)))))
})

test_that("multi_impute() with species_col + covariates runs end-to-end", {
  set.seed(601)
  tree <- ape::rtree(15)
  n_obs <- 45  # 3 obs per species

  df <- data.frame(
    sp = rep(tree$tip.label, each = 3),
    trait1 = abs(rnorm(n_obs)) + 0.5
  )
  df$trait1[sample(n_obs, 9)] <- NA

  covs <- data.frame(acclim_temp = rnorm(n_obs, 25, 5))

  mi <- multi_impute(df, tree, m = 3L, species_col = "sp",
                     covariates = covs,
                     epochs = 20L, verbose = FALSE, seed = 601L,
                     eval_every = 10L, patience = 5L)

  expect_s3_class(mi, "pigauto_mi")
  expect_equal(mi$m, 3L)
  expect_equal(length(mi$datasets), 3L)

  # Each imputed dataset has observation-level rows
  for (d in mi$datasets) {
    expect_equal(nrow(d), n_obs)
    expect_true(all(is.finite(d$trait1)))
  }

  # Pooled point estimate also at obs level
  expect_equal(nrow(mi$pooled_point), n_obs)
})

test_that("multi-obs + covariates: obs-level refinement produces within-species variation", {
  # After the Phase 2 architecture change (obs_refine MLP), observations
  # of the same species with different covariate values should receive
  # different predictions.  The refinement layer re-injects user
  # covariates after the species-level GNN broadcast, enabling
  # covariate-conditional predictions within species.

  set.seed(602)
  tree <- ape::rtree(10)

  # 4 observations per species — varied covariate values within species
  df <- data.frame(
    species = rep(tree$tip.label, each = 4),
    y = abs(rnorm(40)) + 0.5
  )
  covs <- data.frame(
    temp = rep(c(10, 20, 30, 40), times = 10)  # systematic within-species variation
  )

  result <- impute(df, tree, species_col = "species",
                   covariates = covs,
                   epochs = 30L, verbose = FALSE, seed = 602L,
                   eval_every = 10L, patience = 5L,
                   missing_frac = 0.25)

  pred <- result$prediction$imputed
  # Check a few species: obs with different covariates should now get
  # different predictions (non-zero within-species range)
  any_varies <- FALSE
  for (sp in tree$tip.label[1:3]) {
    sp_rows <- which(df$species == sp)
    sp_preds <- pred$y[sp_rows]
    rng <- max(sp_preds) - min(sp_preds)
    if (rng > 1e-6) any_varies <- TRUE
  }
  expect_true(any_varies,
              label = "At least one species should have within-species prediction variation")
})

test_that("multi-obs WITHOUT covariates: within-species predictions are uniform", {
  # Without user covariates, no obs_refine MLP is created, so the
  # species-level broadcast produces identical predictions per species.
  set.seed(604)
  tree <- ape::rtree(10)

  df <- data.frame(
    species = rep(tree$tip.label, each = 4),
    y = abs(rnorm(40)) + 0.5
  )
  # No covariates supplied

  result <- impute(df, tree, species_col = "species",
                   epochs = 30L, verbose = FALSE, seed = 604L,
                   eval_every = 10L, patience = 5L,
                   missing_frac = 0.25)

  pred <- result$prediction$imputed
  for (sp in tree$tip.label[1:3]) {
    sp_rows <- which(df$species == sp)
    sp_preds <- pred$y[sp_rows]
    expect_equal(max(sp_preds) - min(sp_preds), 0,
                 tolerance = 1e-6,
                 label = paste("Within-species range for", sp, "(no covariates)"))
  }
})

test_that("multi-obs + covariates: fit stores multi_obs flag and obs metadata", {
  set.seed(603)
  tree <- ape::rtree(15)
  n_obs <- 30

  df <- data.frame(
    sp = rep(tree$tip.label, each = 2),
    val = abs(rnorm(n_obs)) + 0.5
  )
  covs <- data.frame(env = rnorm(n_obs))

  result <- impute(df, tree, species_col = "sp",
                   covariates = covs,
                   epochs = 20L, verbose = FALSE, seed = 603L,
                   eval_every = 10L, patience = 5L)

  fit <- result$fit
  expect_true(fit$multi_obs)
  expect_equal(length(fit$obs_species), n_obs)
  expect_equal(length(unique(fit$obs_species)), 15L)  # 15 species
  expect_equal(length(fit$species_names), 15L)

  # predict() from the fit also returns obs-level output
  pred2 <- predict(fit, return_se = TRUE)
  expect_equal(nrow(pred2$imputed), n_obs)
  expect_equal(nrow(pred2$se), n_obs)
  expect_true(pred2$multi_obs)
})
