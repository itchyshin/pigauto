# API of multi_impute(draws_method = "posterior") (design.md sections 4 and
# 5a items 7 to 11). Short chains keep this fast; non-convergence warnings
# are expected and muffled.

post_data <- function(n = 30L, seed = 21L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  R <- stats::cov2cor(ape::vcv(tree))
  Z <- t(chol(R + diag(1e-10, n))) %*% matrix(stats::rnorm(2L * n), n) %*%
    chol(matrix(c(1, 0.6, 0.6, 1), 2))
  df <- data.frame(mass = exp(Z[, 1] + 3), wing = Z[, 2] * 2 + 10,
                   row.names = tree$tip.label)
  df$mass[c(2, 5, 9, 14)] <- NA
  df$wing[c(5, 7, 20, 25, 28)] <- NA
  # Shuffle rows so input order differs from tree tip order.
  df <- df[sample.int(n), , drop = FALSE]
  list(df = df, tree = tree)
}

fast_ctl <- list(n_chains = 2L, burnin = 60L, n_iter = 120L, keep_draws = 40L)

run_post <- function(df, tree, m = 4L, seed = 3L, ctl = fast_ctl, ...) {
  suppressWarnings(multi_impute(df, tree, m = m, draws_method = "posterior",
                                posterior_control = ctl, seed = seed,
                                verbose = FALSE, ...))
}

test_that("posterior output has the frozen section-4 shape", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree)
  expect_s3_class(mi, "pigauto_posterior_mi")
  expect_s3_class(mi, "pigauto_mi")
  expect_identical(mi$draws_method, "posterior")
  expect_identical(mi$mi_workflow, "pigauto_posterior_mi_v1")
  expect_length(mi$datasets, 4L)
  expect_identical(mi$m, 4L)
  for (d in mi$datasets) {
    expect_identical(dim(d), dim(pd$df))
    expect_identical(rownames(d), rownames(pd$df))
    expect_false(anyNA(d))
  }
  ci <- mi$posterior$cell_interval
  expect_identical(names(ci), c("row", "trait", "lower", "upper", "median"))
  expect_identical(nrow(ci), sum(is.na(pd$df)))
  expect_true(all(ci$lower <= ci$median & ci$median <= ci$upper))
  expect_true(all(is.na(pd$df[cbind(ci$row, match(ci$trait, names(pd$df)))])))
  # Back-transformed to the original scale (not left on the log scale).
  expect_true(all(ci$median[ci$trait == "mass"] > min(pd$df$mass, na.rm = TRUE) / 5))
  dg <- mi$posterior$diagnostics
  expect_identical(names(dg), c("parameter", "rhat", "ess_bulk"))
  expect_true(is.logical(attr(dg, "converged")))
  expect_setequal(dg$parameter,
                  c("Sigma_P[1,1]", "Sigma_P[1,2]", "Sigma_P[2,2]",
                    "Sigma_E[1,1]", "Sigma_E[1,2]", "Sigma_E[2,2]",
                    "lambda[1]", "lambda[2]"))
  pr <- mi$posterior$params
  expect_identical(names(pr)[1:4], c("Sigma_P", "Sigma_E", "lambda", "mu"))
  expect_equal(dim(pr$Sigma_P), c(2L, 2L, 40L))
  expect_equal(dim(pr$Sigma_E), c(2L, 2L, 40L))
  expect_equal(dim(pr$lambda), c(40L, 2L))
  expect_equal(dim(pr$mu), c(40L, 2L))
  expect_identical(colnames(pr$lambda), c("mass", "wing"))
})

test_that("observed cells are never altered and imputed cells vary", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree, m = 5L)
  obs <- !is.na(as.matrix(pd$df))
  for (d in mi$datasets) {
    expect_identical(as.matrix(d)[obs], as.matrix(pd$df)[obs])
  }
  stack <- vapply(mi$datasets, function(d) as.matrix(d)[!obs], numeric(sum(!obs)))
  expect_true(all(apply(stack, 1L, stats::sd) > 0))
  expect_identical(mi$imputed_mask, is.na(as.matrix(pd$df)))
  # Imputed values lie inside or near their own cell's interval (row
  # alignment): at least most of them are within the 95% interval.
  ci <- mi$posterior$cell_interval
  v <- as.matrix(mi$datasets[[1]])[cbind(ci$row, match(ci$trait, names(pd$df)))]
  expect_gte(mean(v >= ci$lower & v <= ci$upper), 0.6)
})

test_that("the same seed reproduces the draws; another seed does not", {
  pd <- post_data()
  a <- run_post(pd$df, pd$tree, seed = 5L)
  b <- run_post(pd$df, pd$tree, seed = 5L)
  c <- run_post(pd$df, pd$tree, seed = 6L)
  expect_identical(a$datasets, b$datasets)
  expect_identical(a$posterior$params, b$posterior$params)
  expect_false(identical(a$datasets, c$datasets))
  # posterior_control$seed overrides seed.
  d <- run_post(pd$df, pd$tree, seed = 99L, ctl = c(fast_ctl, list(seed = 5L)))
  expect_identical(a$datasets, d$datasets)
})

test_that("non-continuous traits, multi-obs data and covariates are clear errors", {
  pd <- post_data()
  df <- pd$df
  df$diet <- factor(sample(c("a", "b"), nrow(df), replace = TRUE))
  expect_error(run_post(df, pd$tree), "continuous traits only.*diet \\(binary\\)")
  # The message does not send users to a path that with_imputations() refuses.
  expect_error(run_post(df, pd$tree),
               "no analysis-aware MI path.*`with_imputations\\(\\)` and `pool_mi\\(\\)` refuse")
  df2 <- pd$df
  df2$clutch <- sample(1:5, nrow(df2), replace = TRUE)
  expect_error(run_post(df2, pd$tree), "clutch \\(count\\)")
  df3 <- pd$df
  df3$size <- factor(sample(c("s", "m", "l"), nrow(df3), replace = TRUE),
                     levels = c("s", "m", "l"), ordered = TRUE)
  expect_error(run_post(df3, pd$tree), "size \\(ordinal\\)")
  multi <- rbind(pd$df, pd$df[1:3, ])
  multi$species <- c(rownames(pd$df), rownames(pd$df)[1:3])
  rownames(multi) <- NULL
  expect_error(run_post(multi, pd$tree, species_col = "species"),
               "one observation per species")
  # species_col with one row per species: no false claim of duplicates.
  single <- pd$df
  single$species <- rownames(single)
  rownames(single) <- NULL
  expect_error(run_post(single, pd$tree, species_col = "species"),
               "does not support `species_col`.*one row here")
  cov <- data.frame(temp = stats::rnorm(nrow(pd$df)))
  expect_error(run_post(pd$df, pd$tree, covariates = cov),
               "does not support `covariates`")
  expect_error(run_post(pd$df, pd$tree, ctl = list(chains = 2L)),
               "Unknown `posterior_control`")
  expect_error(run_post(pd$df, pd$tree, m = 50L), "keep_draws.*at least m")
})

test_that("GNN and fitting arguments are ignored with a message", {
  pd <- post_data()
  expect_message(run_post(pd$df, pd$tree, gnn = FALSE, epochs = 3L),
                 "ignoring: gnn, epochs")
  expect_no_message(run_post(pd$df, pd$tree))
})

test_that("with_imputations() and pool_mi() run on a posterior object", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree, m = 4L)
  fits <- with_imputations(mi, function(d) stats::lm(wing ~ log(mass), data = d),
                           .progress = FALSE)
  expect_s3_class(fits, "pigauto_mi_fits")
  expect_identical(attr(fits, "mi_workflow"), "pigauto_posterior_mi_v1")
  pooled <- pool_mi(fits)
  expect_s3_class(pooled, "pigauto_pooled")
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(pooled$std.error > 0))
  # A hand-built list stamped with the posterior marker is inconsistent.
  forged <- lapply(mi$datasets, function(d) stats::lm(wing ~ mass, data = d))
  attr(forged, "mi_workflow") <- "pigauto_posterior_mi_v1"
  expect_error(pool_mi(forged), "inconsistent provenance")
  # A posterior-looking object without the class is still refused.
  fake <- unclass(mi)
  class(fake) <- c("pigauto_mi", "list")
  expect_error(with_imputations(fake, stats::lm, .progress = FALSE),
               "Legacy `pigauto_mi`")
  # pool_mi() on the MI object itself says what to do instead.
  expect_error(pool_mi(mi),
               "takes the list of fits returned by `with_imputations\\(\\)`")
  # Diagnostic completions are refused with a pointer to the posterior path.
  diag_obj <- structure(list(datasets = mi$datasets,
                             mi_workflow = "pigauto_diagnostic_mi"),
                        class = c("pigauto_mi", "list"))
  expect_error(with_imputations(diag_obj, stats::lm, .progress = FALSE),
               "prediction-diagnostic.*multi_impute\\(draws_method = \"posterior\"\\)")
})

test_that("print() reports the method and convergence", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree)
  out <- utils::capture.output(print(mi))
  expect_true(any(grepl("posterior multiple imputation", out)))
  expect_true(any(grepl("Converged", out)))
  expect_true(any(grepl("with_imputations", out)))
})

test_that("a non-converged run warns", {
  pd <- post_data()
  expect_warning(
    multi_impute(pd$df, pd$tree, m = 2L, draws_method = "posterior",
                 posterior_control = list(n_chains = 2L, burnin = 5L,
                                          n_iter = 20L, keep_draws = 10L),
                 seed = 1L, verbose = FALSE),
    "did not meet the convergence rule")
})

test_that("improper mode returns fixed parameters through the API", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree,
                 ctl = c(fast_ctl, list(param_uncertainty = "none")))
  sp <- mi$posterior$params$Sigma_P
  expect_equal(max(abs(sweep(sp, c(1L, 2L), sp[, , 1L]))), 0)
  expect_identical(mi$posterior$control$param_uncertainty, "none")
})

test_that("plug-in draws (param_uncertainty = 'none') cannot be pooled", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree,
                 ctl = c(fast_ctl, list(param_uncertainty = "none")))
  expect_s3_class(mi, "pigauto_posterior_mi")
  expect_identical(mi$mi_workflow, "pigauto_posterior_plugin_diagnostic")
  f <- function(d) stats::lm(wing ~ log(mass), data = d)
  expect_error(with_imputations(mi, f, .progress = FALSE),
               "plug-in draws.*param_uncertainty = \"none\"")
  expect_error(pool_mi(mi), "plug-in draws")
  # Fits stamped with the plug-in marker are refused, with or without the
  # pigauto_mi_fits class.
  stamped <- lapply(mi$datasets, f)
  attr(stamped, "mi_workflow") <- "pigauto_posterior_plugin_diagnostic"
  expect_error(pool_mi(stamped), "plug-in draws")
  class(stamped) <- c("pigauto_mi_fits", "list")
  expect_error(pool_mi(stamped), "plug-in draws")
  out <- utils::capture.output(print(mi))
  expect_true(any(grepl("Downstream inference: +unsupported for these draws",
                        out)))
  expect_false(any(grepl("with_imputations\\(mi, f\\) then pool_mi", out)))
})

test_that("a single continuous trait (K = 1) works through multi_impute()", {
  pd <- post_data()
  df1 <- pd$df[, "mass", drop = FALSE]
  for (pu in c("full", "none")) {
    mi <- run_post(df1, pd$tree,
                   ctl = c(fast_ctl, list(param_uncertainty = pu)))
    pr <- mi$posterior$params
    expect_equal(dim(pr$lambda), c(40L, 1L))
    expect_identical(colnames(pr$lambda), "mass")
    expect_equal(as.numeric(pr$lambda),
                 pr$Sigma_P[1, 1, ] / (pr$Sigma_P[1, 1, ] + pr$Sigma_E[1, 1, ]))
    expect_true(all(pr$lambda > 0 & pr$lambda < 1))
    expect_identical(nrow(mi$posterior$cell_interval), sum(is.na(df1$mass)))
    for (d in mi$datasets) expect_false(anyNA(d))
  }
})

test_that("a fully observed input stops before any MCMC is run", {
  pd <- post_data()
  full <- pd$df
  full$mass[is.na(full$mass)] <- 50
  full$wing[is.na(full$wing)] <- 10
  local_mocked_bindings(.mip_fit = function(...) stop("the MCMC ran"))
  expect_error(run_post(full, pd$tree), "no missing cells to impute")
})

test_that("the m datasets are spaced across chains; intervals use every kept draw", {
  pd <- post_data()
  mi <- run_post(pd$df, pd$tree, m = 4L)
  ctl <- mi$posterior$control
  n_kept <- ctl$n_chains * (ctl$n_iter %/% ctl$thin)       # 2 chains x 20
  expect_identical(n_kept, 40L)
  idx <- mi$posterior$draw_index
  expect_equal(idx, round(seq(1, n_kept, length.out = 4L)))
  per_chain <- n_kept / ctl$n_chains
  expect_true(any(idx <= per_chain) && any(idx > per_chain))  # both chains
  # Rebuild the same run (.mip_fit reseeds from ctl$seed) and decode every
  # kept draw to the original scale.
  X <- mi$data$X_scaled
  rownames(X) <- mi$data$species_names
  X <- X[pd$tree$tip.label, , drop = FALSE]
  f <- .mip_fit(X, pd$tree, ctl)
  expect_equal(ncol(f$ymis), n_kept)
  dec <- f$ymis
  for (tm in mi$data$trait_map) {
    sel <- f$miss[, 2L] == tm$latent_cols[1L]
    v <- dec[sel, , drop = FALSE] * tm$sd + tm$mean
    if (isTRUE(tm$log_transform)) v <- exp(v)
    dec[sel, ] <- v
  }
  q <- t(apply(dec, 1L, stats::quantile, probs = c(0.025, 0.5, 0.975),
               names = FALSE))
  ci <- mi$posterior$cell_interval
  expect_identical(nrow(ci), nrow(f$miss))
  # The intervals are the quantiles of all kept draws...
  expect_equal(ci$lower, q[, 1L], tolerance = 1e-10)
  expect_equal(ci$median, q[, 2L], tolerance = 1e-10)
  expect_equal(ci$upper, q[, 3L], tolerance = 1e-10)
  # ...the m datasets are the kept draws at draw_index...
  j <- match(ci$trait, names(pd$df))
  stack <- vapply(mi$datasets, function(d) as.matrix(d)[cbind(ci$row, j)],
                  numeric(nrow(ci)))
  expect_equal(stack, dec[, idx], tolerance = 1e-10, ignore_attr = TRUE)
  # ...and the intervals are not the spread of those m draws.
  expect_false(isTRUE(all.equal(
    ci$upper - ci$lower, apply(stack, 1L, function(v) diff(range(v))))))
  expect_false(isTRUE(all.equal(
    ci$lower, apply(stack, 1L, stats::quantile, 0.025, names = FALSE))))
})

test_that("conformal multi_impute() objects are still refused by with_imputations()", {
  skip_if_no_libtorch()
  pd <- post_data()
  mi <- suppressWarnings(
    multi_impute(pd$df, pd$tree, m = 2L, draws_method = "conformal",
                 epochs = 5L, verbose = FALSE, seed = 1L))
  expect_error(with_imputations(mi, stats::lm, .progress = FALSE),
               "prediction-diagnostic completions.*draws_method = \"posterior\"")
})

test_that("param_uncertainty = 'both' adds plug-in draws from the same run", {
  pd <- post_data()
  full <- run_post(pd$df, pd$tree, m = 4L, seed = 8L)
  both <- run_post(pd$df, pd$tree, m = 4L, seed = 8L,
                   ctl = c(fast_ctl, list(param_uncertainty = "both")))
  # 1. The proper results are exactly those of "full" with the same seed.
  expect_identical(both$datasets, full$datasets)
  expect_identical(both$posterior$params, full$posterior$params)
  expect_identical(both$posterior$cell_interval, full$posterior$cell_interval)
  expect_null(full$posterior_improper)
  # The proper draws keep the poolable marker in "full" and "both" modes.
  expect_identical(full$mi_workflow, "pigauto_posterior_mi_v1")
  expect_identical(both$mi_workflow, "pigauto_posterior_mi_v1")
  # 2. m plug-in datasets; observed cells unchanged; imputed cells filled.
  imp <- both$posterior_improper
  expect_length(imp$datasets, 4L)
  obs <- !is.na(as.matrix(pd$df))
  for (d in imp$datasets) {
    expect_identical(as.matrix(d)[obs], as.matrix(pd$df)[obs])
    expect_false(anyNA(d))
  }
  expect_false(identical(imp$datasets, both$datasets))
  # The plug-in covariances are the posterior means of the kept draws, and
  # the plug-in intervals cover the same cells with the same columns.
  expect_equal(imp$Sigma_P, apply(both$posterior$params$Sigma_P, c(1L, 2L), mean),
               ignore_attr = TRUE)
  expect_equal(imp$Sigma_E, apply(both$posterior$params$Sigma_E, c(1L, 2L), mean),
               ignore_attr = TRUE)
  expect_identical(names(imp$cell_interval), names(both$posterior$cell_interval))
  expect_identical(imp$cell_interval[, c("row", "trait")],
                   both$posterior$cell_interval[, c("row", "trait")])
  expect_true(all(imp$cell_interval$lower <= imp$cell_interval$upper))
  # Item 3 of the request (plug-in variance <= proper variance) is NOT
  # asserted: computed exactly from the kept parameter draws on this fixture
  # (converged chains), the plug-in predictive variance at the posterior-mean
  # covariances is 1.4-1.6% LARGER than the proper one. The conditional
  # variance is concave in the covariance, so the Jensen gap of the
  # posterior-mean plug-in (about 0.010 on the latent scale) exceeds the
  # parameter-uncertainty term Var(E[y | theta]) (about 0.005).
})
