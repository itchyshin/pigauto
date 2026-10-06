# Bounded S1 effective-path checks. Keep every real fit baseline-only and tiny;
# longer posterior-chain and GNN benchmark paths belong to other evidence.

test_that("exported entry points resolve and policy defaults stay aligned", {
  exports <- getNamespaceExports("pigauto")
  expect_length(exports, 33L)
  for (name in exports) {
    expect_true(is.function(getExportedValue("pigauto", name)), info = name)
  }
  expect_false(formals(impute)$gnn)
  expect_false(formals(fit_pigauto)$gnn)
  expect_false(formals(multi_impute)$gnn)
  expect_false(formals(multi_impute_trees)$gnn)
  expect_true(formals(compare_methods)$gnn)
  expect_true(formals(simulate_benchmark)$gnn)
  expect_identical(eval(formals(multi_impute)$draws_method)[1L], "auto")
  expect_identical(eval(formals(fit_baseline)$predict_method)[1L], "auto")
  expect_identical(eval(formals(fit_baseline)$lambda_mode)[1L], "estimate")
  expect_identical(eval(formals(fit_baseline)$discrete_lambda)[1L], "estimate")
  expect_identical(eval(formals(fit_baseline)$joint_solver)[1L], "inhouse")
  expect_false(formals(impute)$safety_floor)
  expect_false(formals(impute)$phylo_signal_gate)
})

cran_defaults_fixture <- function(n = 24L, seed = 1106L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  comp <- matrix(stats::rgamma(n * 3L, shape = 2), ncol = 3L,
                 dimnames = list(sp, c("comp_a", "comp_b", "comp_c")))
  comp <- comp / rowSums(comp)
  df <- data.frame(
    continuous = exp(stats::rnorm(n, mean = 2, sd = 0.3)),
    count = as.integer(stats::rpois(n, 3)),
    binary = factor(rep(c("no", "yes"), length.out = n)),
    categorical = factor(rep(c("a", "b", "c"), length.out = n)),
    ordinal = ordered(rep(c("low", "mid", "high"), length.out = n),
                      levels = c("low", "mid", "high")),
    proportion = stats::runif(n, 0.1, 0.9),
    zi = as.integer(ifelse(seq_len(n) %% 3L == 0L, 0L,
                           stats::rpois(n, 2) + 1L)),
    row.names = sp
  )
  df <- cbind(df, as.data.frame(comp))
  # Keep observed levels plentiful while exercising cell and row masking.
  for (nm in c("continuous", "count", "binary", "categorical",
               "ordinal", "proportion", "zi")) {
    df[[nm]][c(2L, 8L)] <- NA
  }
  df[c(3L, 11L), c("comp_a", "comp_b", "comp_c")] <- NA
  list(tree = tree, traits = df)
}

cran_defaults_single <- function(df, tree, col, type = NULL) {
  x <- df[, col, drop = FALSE]
  override <- if (is.null(type)) NULL else stats::setNames(type, col)
  preprocess_traits(x, tree, trait_types = override)
}

test_that("all eight type-specific baseline paths accept one trait at a time", {
  fx <- cran_defaults_fixture()
  cases <- list(
    continuous = list(col = "continuous", type = NULL),
    count = list(col = "count", type = NULL),
    binary = list(col = "binary", type = NULL),
    categorical = list(col = "categorical", type = NULL),
    ordinal = list(col = "ordinal", type = NULL),
    proportion = list(col = "proportion", type = "proportion"),
    zi_count = list(col = "zi", type = "zi_count"),
    multi_proportion = list(col = c("comp_a", "comp_b", "comp_c"), type = NULL)
  )

  for (expected in names(cases)) {
    spec <- cases[[expected]]
    x <- fx$traits[, spec$col, drop = FALSE]
    if (expected == "multi_proportion") {
      pd <- preprocess_traits(x, fx$tree,
        multi_proportion_groups = list(comp = names(x)))
    } else {
      override <- if (is.null(spec$type)) NULL else
        stats::setNames(spec$type, spec$col)
      pd <- preprocess_traits(x, fx$tree, trait_types = override)
    }
    expect_identical(pd$trait_map[[1L]]$type, expected, info = expected)
    bl <- fit_baseline(pd, fx$tree)
    expect_true(all(is.finite(bl$mu)), info = expected)
    expect_equal(dim(bl$mu), dim(pd$X_scaled), info = expected)
  }
})

test_that("default mixed imputation preserves observations and masks compositions by row", {
  fx <- cran_defaults_fixture()
  observed <- !is.na(fx$traits)
  result <- suppressWarnings(impute(
    fx$traits, fx$tree,
    multi_proportion_groups = list(comp = c("comp_a", "comp_b", "comp_c")),
    verbose = FALSE, seed = 1106L
  ))

  expect_false(isTRUE(result$fit$model_config$gnn))
  expect_false(result$fit$safety_floor)
  expect_false(any(result$fit$phylo_gate_triggered))
  expect_identical(result$fit$model_config$lambda_mode, "estimate")
  expect_identical(result$fit$model_config$discrete_lambda, "estimate")
  expect_identical(result$fit$model_config$joint_solver, "inhouse")
  expect_true(result$data$trait_map$continuous$log_transform)
  expect_true("continuous" %in% names(result$fit$baseline$path))
  for (nm in names(fx$traits)) {
    expect_equal(result$completed[[nm]][observed[, nm]],
                 fx$traits[[nm]][observed[, nm]],
                 tolerance = 1e-10, info = nm)
  }
  group <- c("comp_a", "comp_b", "comp_c")
  miss_rows <- which(rowSums(is.na(fx$traits[, group, drop = FALSE])) > 0L)
  expect_true(all(is.na(fx$traits[miss_rows, group, drop = FALSE])))
  expect_true(all(is.finite(as.matrix(result$completed[miss_rows, group]))))
})

test_that("transformation and baseline route overrides are recorded", {
  fx <- cran_defaults_fixture()
  result <- suppressWarnings(impute(
    fx$traits[, c("continuous", "proportion"), drop = FALSE], fx$tree,
    trait_types = c(proportion = "proportion"),
    log_transform = FALSE, lambda_mode = "fixed_1",
    discrete_lambda = "fixed_1", joint_solver = "inhouse",
    predict_method = "per_column", verbose = FALSE, seed = 1107L
  ))
  expect_false(result$data$trait_map$continuous$log_transform)
  expect_identical(result$fit$model_config$lambda_mode, "fixed_1")
  expect_identical(result$fit$model_config$discrete_lambda, "fixed_1")
  expect_identical(result$fit$model_config$joint_solver, "inhouse")
  expect_true(all(result$fit$model_config$predict_method_by_trait == "per_column"))

  # With no validation cells, auto records the actual per-column fallback.
  pd <- preprocess_traits(fx$traits[, "continuous", drop = FALSE], fx$tree)
  fallback <- fit_baseline(pd, fx$tree, predict_method = "auto")
  expect_identical(unname(fallback$predict_method_by_trait), "per_column")
})

test_that("multi-observation covariates warn under the default GNN-off path", {
  fx <- cran_defaults_fixture(n = 12L, seed = 1108L)
  dat <- data.frame(
    species = rep(fx$tree$tip.label, each = 2L),
    value = rep(fx$traits$continuous, each = 2L),
    row.names = paste0("obs", seq_len(24L))
  )
  dat$value[c(2L, 9L)] <- NA_real_
  covs <- matrix(seq_len(nrow(dat)), ncol = 1L,
                 dimnames = list(rownames(dat), "temperature"))
  warnings <- character()
  res <- withCallingHandlers(
    impute(dat, fx$tree, species_col = "species",
      covariates = covs, verbose = FALSE, seed = 1108L),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("covariates are used only by the GNN", warnings)))
  expect_false(isTRUE(res$fit$model_config$gnn))
  expect_true(res$data$multi_obs)
})

test_that("two-tree diagnostic runs with the default GNN-off fit", {
  fx <- cran_defaults_fixture(n = 12L, seed = 1109L)
  tree2 <- ape::rtree(length(fx$tree$tip.label))
  tree2$tip.label <- fx$tree$tip.label
  trees <- list(fx$tree, tree2)
  out <- suppressWarnings(multi_impute_trees(
    fx$traits[, c("continuous", "count"), drop = FALSE], trees,
    m_per_tree = 1L, verbose = FALSE, seed = 1109L
  ))
  expect_s3_class(out, "pigauto_mi_trees")
  expect_equal(out$n_trees, 2L)
  expect_false(isTRUE(out$fit$model_config$gnn))
  expect_identical(out$mi_workflow, "pigauto_tree_sensitivity_diagnostic")
  expect_error(with_imputations(out, function(data) stats::lm(continuous ~ count,
                                                              data = data)))
})

test_that("saved default GNN-off fit reloads and predicts", {
  fx <- cran_defaults_fixture()
  res <- suppressWarnings(impute(
    fx$traits[, c("continuous", "count"), drop = FALSE], fx$tree,
    verbose = FALSE, seed = 1110L
  ))
  path <- withr::local_tempfile(fileext = ".pigauto")
  suppressMessages(save_pigauto(res$fit, path))
  loaded <- load_pigauto(path)
  expect_false(isTRUE(loaded$model_config$gnn))
  expect_equal(predict(loaded)$imputed, predict(res$fit)$imputed)
})

test_that("wrapper gnn defaults and overrides reach fit_pigauto", {
  fx <- cran_defaults_fixture(n = 12L, seed = 1111L)
  pd <- preprocess_traits(fx$traits[, "continuous", drop = FALSE], fx$tree)
  seen <- NULL
  testthat::local_mocked_bindings(
    fit_pigauto = function(...) {
      seen <<- list(...)
      stop("forwarding captured", call. = FALSE)
    }, .package = "pigauto"
  )

  expect_error(cross_validate(pd, fx$tree, k = 2L, seeds = 1L),
               "forwarding captured")
  expect_null(seen$gnn)
  expect_error(cross_validate(pd, fx$tree, k = 2L, seeds = 1L,
                              gnn = TRUE), "forwarding captured")
  expect_true(isTRUE(seen$gnn))
  expect_error(compare_methods(pd, fx$tree, seeds = 1L),
               "forwarding captured")
  expect_true(isTRUE(seen$gnn))
  expect_error(compare_methods(pd, fx$tree, seeds = 1L, gnn = FALSE),
               "forwarding captured")
  expect_false(seen$gnn)
})

test_that("benchmark helper forwards its default GNN-on arm without fitting it", {
  fx <- cran_defaults_fixture(n = 10L, seed = 1112L)
  seen <- NULL
  testthat::local_mocked_bindings(
    simulate_bm_traits = function(tree, n_traits, seed = NULL) {
      fx$traits[tree$tip.label, "continuous", drop = FALSE]
    },
    fit_pigauto = function(...) {
      seen <<- list(...)
      stop("forwarding captured", call. = FALSE)
    }, .package = "pigauto"
  )
  invisible(try(simulate_benchmark(
    n_species = 10L, n_traits = 1L, scenarios = "BM", n_reps = 1L,
    seed = 1112L, verbose = FALSE
  ), silent = TRUE))
  expect_false(is.null(seen))
  expect_true(isTRUE(seen$gnn))
})

test_that("automatic multiple-imputation route resolves without sampling", {
  fx <- cran_defaults_fixture(n = 12L, seed = 1113L)
  route <- function(cols, species_col = NULL, covariates = NULL,
                    multi_proportion_groups = NULL) {
    suppressMessages(.mi_resolve_draws_auto(
      traits = fx$traits[, cols, drop = FALSE], tree = fx$tree,
      species_col = species_col, trait_types = NULL,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = TRUE, covariates = covariates, verbose = FALSE))
  }
  expect_identical(route("continuous"), "posterior")
  expect_identical(route(c("continuous", "binary")), "conformal")
  expect_identical(route("continuous", species_col = "species"), "conformal")
  expect_identical(route("continuous", covariates = matrix(1)), "conformal")
  expect_identical(route(c("comp_a", "comp_b", "comp_c"),
                         multi_proportion_groups = list(comp =
                           c("comp_a", "comp_b", "comp_c"))), "conformal")
})

test_that("multi_impute actually dispatches the automatic route", {
  fx <- cran_defaults_fixture(n = 12L, seed = 1114L)
  reached <- NULL
  testthat::local_mocked_bindings(
    .multi_impute_posterior = function(...) {
      reached <<- "posterior"
      list(route = "posterior")
    },
    impute = function(...) {
      reached <<- "conformal"
      stop("conformal route captured", call. = FALSE)
    }, .package = "pigauto"
  )
  out <- multi_impute(fx$traits[, "continuous", drop = FALSE], fx$tree,
                      m = 2L, verbose = FALSE)
  expect_identical(out$route, "posterior")
  expect_identical(reached, "posterior")
  expect_error(multi_impute(fx$traits[, c("continuous", "binary"), drop = FALSE],
                            fx$tree, m = 2L, verbose = FALSE),
               "conformal route captured")
  expect_identical(reached, "conformal")
})
