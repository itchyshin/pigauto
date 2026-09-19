# Tests for gnn = FALSE (pure phylogenetic-baseline mode), per the frozen
# contract docs/dev-log/arc/2026-09-18-gnn-off-contract.md.
#
# gnn = FALSE tests deliberately do NOT call skip_if_no_libtorch(): the whole
# point of gnn = FALSE is that fit_pigauto()/predict() make zero torch::
# calls, so these tests must run even when the torch backend cannot execute
# a real tensor op. torch is a hard Import of pigauto, so torch::torch_tensor
# and the pigauto::ResidualPhyloDAE generator are always resolvable for
# trace()-based call counting even in that state. The one exception is test
# 10 (the gnn = TRUE default path), which actually trains a GNN and so needs
# a working torch backend.

make_gnn_off_fixture <- function(n = 40L, seed = 11L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df <- data.frame(
    mass      = exp(ape::rTraitCont(tree, sigma = 0.6) + 3),
    clutch    = as.integer(pmax(1L, round(exp(ape::rTraitCont(tree, sigma = 0.4) + 1)))),
    nocturnal = factor(ifelse(ape::rTraitCont(tree) > 0, "yes", "no")),
    diet      = factor(sample(c("herb", "omni", "carn"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  rownames(df) <- sp
  set.seed(seed + 1L)
  for (j in seq_along(df)) df[sample(n, 6L), j] <- NA
  list(tree = tree, df = df)
}

# The gnn = FALSE three-way blend collapses to two terms (r_gnn is always 0):
#   pred = r_bm * mu + r_mean * mean_baseline_per_col
# ("predict.pigauto_fit under gnn = FALSE" in the frozen contract). `mu` is
# either baseline_full$mu (production mode) or baseline$mu (evaluation mode);
# r_cal_bm / r_cal_mean / mean_baseline_per_col are shared across both.
.gnn_off_blend <- function(fit, mu) {
  latent_names <- names(fit$r_cal_bm)
  p <- length(latent_names)

  r_mean <- fit$r_cal_mean
  if (is.null(r_mean)) {
    r_mean <- stats::setNames(rep(0, p), latent_names)
  } else {
    r_mean <- r_mean[latent_names]
  }
  mean_bl <- fit$mean_baseline_per_col
  if (is.null(mean_bl)) {
    mean_bl <- stats::setNames(rep(0, p), latent_names)
  } else {
    mean_bl <- mean_bl[latent_names]
  }

  sweep(mu, 2, fit$r_cal_bm[latent_names], `*`) +
    matrix(r_mean * mean_bl, nrow = nrow(mu), ncol = p, byrow = TRUE)
}


test_that("fit_pigauto(gnn = FALSE) returns a contract-conformant fit with zero torch calls", {
  skip_on_cran()
  fx  <- make_gnn_off_fixture()
  pd  <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 3L,
                              trait_map = pd$trait_map)

  n_tensor <- 0L
  trace(torch::torch_tensor, tracer = function() n_tensor <<- n_tensor + 1L,
        print = FALSE, where = asNamespace("torch"))
  on.exit(try(untrace(torch::torch_tensor, where = asNamespace("torch")),
              silent = TRUE), add = TRUE)

  fit <- fit_pigauto(pd, fx$tree, splits = spl, gnn = FALSE, verbose = FALSE,
                      seed = 1L)

  expect_s3_class(fit, "pigauto_fit")
  expect_true(isFALSE(fit$model_config$gnn))

  p <- ncol(pd$X_scaled)
  expect_named(fit$r_cal_gnn)
  expect_length(fit$r_cal_gnn, p)
  expect_true(all(fit$r_cal_gnn == 0))

  expect_type(fit$model_state, "list")
  expect_length(fit$model_state, 0L)

  expect_s3_class(fit$history, "data.frame")
  expect_equal(nrow(fit$history), 0L)

  expect_true(is.numeric(fit$val_rmse))
  expect_length(fit$val_rmse, 1L)

  expect_false(is.null(fit$conformal_scores))
  expect_true(is.finite(fit$conformal_scores[["mass"]]))

  expect_false(is.null(fit$baseline$path))
  expect_true(is.character(fit$baseline$path))

  expect_false(is.null(fit$baseline_full$mu))
  expect_equal(dim(fit$baseline_full$mu), dim(fit$baseline$mu))

  expect_equal(n_tensor, 0L)
})


test_that("predict() on a GNN-off fit builds no model and reproduces the blended baseline", {
  skip_on_cran()
  fx  <- make_gnn_off_fixture()
  pd  <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 3L,
                              trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, fx$tree, splits = spl, gnn = FALSE, verbose = FALSE,
                      seed = 1L)

  # base::trace() cannot instrument ResidualPhyloDAE directly: it is a
  # torch::nn_module generator, whose class attribute is
  # c("ResidualPhyloDAE", "nn_module", "nn_module_generator") rather than
  # plain "function" -- trace() sees the extra classes and routes through
  # methods:::.TraceWithMethods()/.classEnv(), which then fails looking for
  # an S4 class definition named "ResidualPhyloDAE" that does not exist.
  # This is a limitation of base::trace() on multiply-classed closures, not
  # a pigauto bug. torch::nn_parameter() is the plain-function call every
  # ResidualPhyloDAE$initialize() makes to register a learnable weight (92
  # calls were observed tracing it through a real GNN-on fit), so it is used
  # here as the "no model was constructed" proxy instead.
  n_model  <- 0L
  n_tensor <- 0L
  trace(torch::nn_parameter, tracer = function() n_model <<- n_model + 1L,
        print = FALSE, where = asNamespace("torch"))
  on.exit(try(untrace(torch::nn_parameter, where = asNamespace("torch")),
              silent = TRUE), add = TRUE)
  trace(torch::torch_tensor, tracer = function() n_tensor <<- n_tensor + 1L,
        print = FALSE, where = asNamespace("torch"))
  on.exit(try(untrace(torch::torch_tensor, where = asNamespace("torch")),
              silent = TRUE), add = TRUE)

  pred  <- predict(fit, return_se = TRUE)                     # production mode
  predm <- predict(fit, return_se = TRUE,                     # evaluation mode
                    .mask_observed_idx = c(spl$val_idx, spl$test_idx))

  expect_equal(n_model, 0L)
  expect_equal(n_tensor, 0L)

  miss <- is.na(pd$X_scaled)
  expected_full <- .gnn_off_blend(fit, fit$baseline_full$mu)
  expect_equal(pred$imputed_latent[miss], expected_full[miss], tolerance = 1e-6)

  hold <- matrix(FALSE, nrow(pd$X_scaled), ncol(pd$X_scaled))
  hold[c(spl$val_idx, spl$test_idx)] <- TRUE
  expected_held <- .gnn_off_blend(fit, fit$baseline$mu)
  expect_equal(predm$imputed_latent[hold], expected_held[hold], tolerance = 1e-6)

  expect_true(all(is.finite(pred$conformal_lower[, "mass"])))
  expect_true(all(is.finite(pred$conformal_upper[, "mass"])))
  expect_true(all(pred$conformal_lower[, "mass"] <= pred$imputed$mass + 1e-8))
  expect_true(all(pred$imputed$mass <= pred$conformal_upper[, "mass"] + 1e-8))
})


test_that("predict(n_imputations = 5) on a GNN-off fit yields BM-posterior draws", {
  skip_on_cran()
  fx  <- make_gnn_off_fixture()
  pd  <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 3L,
                              trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, fx$tree, splits = spl, gnn = FALSE, verbose = FALSE,
                      seed = 1L)

  pred5 <- predict(fit, n_imputations = 5L)
  expect_length(pred5$imputed_datasets, 5L)

  miss_mass <- is.na(pd$X_scaled[, "mass"])
  ds_mass <- vapply(pred5$imputed_datasets, function(d) d$mass,
                     numeric(nrow(pd$X_scaled)))

  expect_true(any(apply(ds_mass[miss_mass, , drop = FALSE], 1,
                         function(row) length(unique(row)) > 1L)),
              info = "draws should differ across datasets at missing cells")
  expect_true(all(apply(ds_mass[!miss_mass, , drop = FALSE], 1,
                         function(row) length(unique(row)) == 1L)),
              info = "draws must be identical across datasets at observed cells")

  pred5_again <- predict(fit, n_imputations = 5L)
  expect_equal(pred5_again$imputed_datasets, pred5$imputed_datasets)

  set.seed(1L)
  rng_before <- .Random.seed
  invisible(predict(fit, n_imputations = 5L))
  rng_after <- .Random.seed
  expect_identical(rng_before, rng_after)
})


test_that("impute(gnn = FALSE) runs the whole pipeline fast", {
  skip_on_cran()
  fx <- make_gnn_off_fixture()

  t0  <- Sys.time()
  res <- impute(fx$df, fx$tree, gnn = FALSE, verbose = FALSE, seed = 5L)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  expect_s3_class(res, "pigauto_result")
  for (col in c("mass", "clutch", "nocturnal", "diet")) {
    expect_false(anyNA(res$completed[[col]]), info = col)
  }
  expect_lt(wall, 60)
  expect_false(res$fit$model_config$gnn)
})


test_that("evaluate() on a GNN-off fit scores the held-out baseline (no leakage)", {
  skip_on_cran()
  fx <- make_gnn_off_fixture()

  # Leakage-freedom (this test's namesake claim) is a property of the two
  # baseline fits (held-out vs full) and does not depend on the calibrated
  # blend, so it is checked on the default gnn = FALSE arm.
  res <- impute(fx$df, fx$tree, gnn = FALSE, verbose = FALSE, seed = 5L)
  held <- res$fit$baseline$mu
  full <- res$fit$baseline_full$mu
  test_mask <- matrix(FALSE, nrow(held), ncol(held))
  test_mask[res$splits$test_idx] <- TRUE
  expect_false(isTRUE(all.equal(held[test_mask], full[test_mask])),
               label = "held-out and full baselines must differ at test cells")

  # The default gnn = FALSE arm keeps GNN-on's safety-floor blend (r_cal_mean
  # can be > 0 even with no GNN), so pigauto's blended prediction need not
  # equal the raw baseline mu there. The pigauto-rmse == baseline-rmse
  # equality only holds on the PURE traditional-stats arm (safety_floor =
  # FALSE, phylo_signal_gate = FALSE), where r_cal_bm is pinned to 1 and
  # evaluate()'s "baseline" row (fit$baseline$mu) is exactly what predict()
  # blends to.
  res_pure <- impute(fx$df, fx$tree, gnn = FALSE, safety_floor = FALSE,
                      phylo_signal_gate = FALSE, verbose = FALSE, seed = 5L)
  expect_equal(unname(res_pure$fit$r_cal_bm),
               rep(1, length(res_pure$fit$r_cal_bm)), tolerance = 1e-8)

  ev <- evaluate(res_pure$fit, data = res_pure$data, splits = res_pure$splits)
  pg <- ev[ev$method == "pigauto" & ev$metric == "rmse", ]
  bl <- ev[ev$method == "baseline" & ev$metric == "rmse", ]
  m  <- merge(pg, bl, by = "trait", suffixes = c("_pg", "_bl"))
  m  <- m[m$trait %in% c("mass", "clutch"), ]

  expect_equal(nrow(m), 2L)
  expect_equal(m$value_pg, m$value_bl, tolerance = 1e-6)
})


test_that("multi_impute(gnn = FALSE) works with both draw methods", {
  skip_on_cran()
  fx <- make_gnn_off_fixture()

  mi <- multi_impute(fx$df, fx$tree, m = 3L, gnn = FALSE, verbose = FALSE,
                      seed = 5L)
  expect_s3_class(mi, "pigauto_mi")
  expect_length(mi$datasets, 3L)

  expect_message(
    mi2 <- multi_impute(fx$df, fx$tree, m = 2L, gnn = FALSE,
                         draws_method = "mc_dropout", verbose = FALSE,
                         seed = 5L),
    "BM-posterior"
  )
  expect_s3_class(mi2, "pigauto_mi")
  expect_length(mi2$datasets, 2L)
})


test_that("multi_impute_trees(gnn = FALSE) works and calls no torch", {
  skip_on_cran()
  fx <- make_gnn_off_fixture()
  trees <- list(fx$tree, fx$tree)
  class(trees) <- "multiPhylo"

  n_tensor <- 0L
  trace(torch::torch_tensor, tracer = function() n_tensor <<- n_tensor + 1L,
        print = FALSE, where = asNamespace("torch"))
  on.exit(try(untrace(torch::torch_tensor, where = asNamespace("torch")),
              silent = TRUE), add = TRUE)

  mi_t <- suppressWarnings(multi_impute_trees(
    fx$df, trees, m_per_tree = 2L, gnn = FALSE, verbose = FALSE, seed = 5L
  ))

  expect_s3_class(mi_t, "pigauto_mi_trees")
  expect_length(mi_t$datasets, 4L)
  expect_equal(n_tensor, 0L)
})


test_that("cross_validate(gnn = FALSE) runs", {
  skip_on_cran()
  fx <- make_gnn_off_fixture()
  pd <- preprocess_traits(fx$df, fx$tree)

  cv <- cross_validate(pd, fx$tree, k = 3L, seeds = 1L, gnn = FALSE,
                        verbose = FALSE)

  expect_s3_class(cv, "pigauto_cv")
  expect_true(is.data.frame(cv$results))
})


test_that("save_pigauto / load_pigauto round-trip a GNN-off fit", {
  skip_on_cran()
  fx  <- make_gnn_off_fixture()
  pd  <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 3L,
                              trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, fx$tree, splits = spl, gnn = FALSE, verbose = FALSE,
                      seed = 1L)

  tf <- withr::local_tempfile(fileext = ".pigauto")
  suppressMessages(save_pigauto(fit, tf))
  fit2 <- load_pigauto(tf)

  expect_true(isFALSE(fit2$model_config$gnn))
  expect_equal(predict(fit2)$imputed, predict(fit)$imputed)
})


test_that("gnn = TRUE default path is unchanged", {
  skip_on_cran()
  skip_if_no_libtorch()
  fx <- make_gnn_off_fixture()

  res <- impute(fx$df, fx$tree, epochs = 5L, verbose = FALSE)

  expect_true(res$fit$model_config$gnn)
  bf <- res$fit$baseline_full
  if (!is.null(bf)) {
    expect_equal(dim(bf$mu), dim(res$fit$baseline$mu))
  }
})


test_that("pure traditional-stats arm: gnn = FALSE with safety_floor = FALSE and phylo_signal_gate = FALSE", {
  skip_on_cran()
  fx  <- make_gnn_off_fixture()
  pd  <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 3L,
                              trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, fx$tree, splits = spl, gnn = FALSE,
                      safety_floor = FALSE, phylo_signal_gate = FALSE,
                      verbose = FALSE, seed = 1L)

  expect_equal(unname(fit$r_cal_bm), rep(1, length(fit$r_cal_bm)),
               tolerance = 1e-8)
  if (!is.null(fit$r_cal_mean)) {
    expect_equal(unname(fit$r_cal_mean), rep(0, length(fit$r_cal_mean)),
                 tolerance = 1e-8)
  }
})
