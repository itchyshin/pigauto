test_that("legacy decoding retains scale, names and both sources of uncertainty", {
  obj <- list(species_names = c("a", "b"), trait_names = "mass",
              norm = list(sds = 2, means = 1, log_transform = FALSE))
  base <- list(se = matrix(c(.5, 1), ncol = 1))
  a <- matrix(c(0, 1), ncol = 1)
  b <- matrix(c(2, 3), ncol = 1)
  one <- decode_continuous_legacy(list(a), obj, 0, TRUE, 1L, base)
  expect_equal(unname(one$imputed[, 1]), c(1, 3))
  expect_equal(unname(one$se[, 1]), c(1, 2))
  expect_equal(dimnames(one$imputed), list(c("a", "b"), "mass"))
  two <- decode_continuous_legacy(list(a, b), obj, 0, TRUE, 2L, base)
  expect_equal(unname(two$imputed[, 1]), c(3, 5))
  expect_equal(unname(two$se[, 1]), sqrt(c(1, 4) + 8))
  obj$norm$log_transform <- TRUE
  log_one <- decode_continuous_legacy(list(a), obj, 0, TRUE, 1L, base)
  expect_equal(unname(log_one$imputed[, 1]), exp(c(1, 3)))
  expect_equal(unname(log_one$se[, 1]), exp(c(1, 3)) * c(1, 2))
  log_two <- decode_continuous_legacy(list(a, b), obj, 0, TRUE, 2L, base)
  expect_equal(unname(log_two$se[, 1]),
               sqrt((exp(c(3, 5)) * c(1, 2))^2 +
                      (exp(c(5, 7)) - exp(c(1, 3)))^2 / 2))
  expect_null(decode_continuous_legacy(list(a), obj, 0, FALSE, 1L, base)$se)
})

test_that("PMM reconstructs donors on their original scales and preserves observations", {
  tm <- list(list(name = "n", type = "count", latent_cols = 1, mean = 0, sd = 1),
             list(name = "p", type = "proportion", latent_cols = 2, mean = 0, sd = 1),
             list(name = "z", type = "zi_count", latent_cols = 3:4, mean = 0, sd = 1))
  x <- cbind(log1p(c(2, NA, 8)), qlogis(c(.2, NA, .8)),
             c(0, NA, 1), log1p(c(NA, NA, 4)))
  rownames(x) <- c("a", "b", "c")
  expect_null(recover_X_orig(NULL, tm))
  orig <- recover_X_orig(x, tm)
  expect_equal(orig$n, c(2L, NA_integer_, 8L))
  expect_equal(orig$p, c(.2, NA, .8))
  expect_equal(orig$z, c(0L, NA_integer_, 4L))
  expect_equal(rownames(orig), c("a", "b", "c"))
  decoded <- list(list(imputed = data.frame(n = c(2, 7, 8), p = c(.2, .7, .8), z = c(0, 3, 4))))
  got <- apply_pmm_to_decoded(decoded, tm, orig, K = 1L)
  expect_identical(got[[1]]$imputed$n, c(2L, 8L, 8L))
  expect_equal(got[[1]]$imputed$p, c(.2, .8, .8))
  expect_identical(got[[1]]$imputed$z, c(0L, 4L, 4L))
  mismatch <- apply_pmm_to_decoded(decoded, tm, orig[1:2, ], K = 1L)
  expect_equal(mismatch, decoded)
})

test_that("conditional Gaussian prior agrees with hand calculations", {
  S <- matrix(c(4, 1, 1, 9), 2)
  L <- rbind(c(2, 3), c(NA, NA))
  p <- build_conditional_prior(S, L, eps = 0)
  expect_equal(p$mu_prior[1, ], c(1 / 3, .5))
  expect_equal(p$sd_prior[1, ], sqrt(c(4 - 1 / 9, 9 - 1 / 4)))
  expect_equal(p$mu_prior[2, ], c(0, 0))
  expect_equal(p$sd_prior[2, ], c(2, 3))
  expect_error(build_conditional_prior(matrix(1, 2, 3), L), "square")
  expect_error(build_conditional_prior(diag(3), L), "ncol")
})

test_that("lambda validation and fixed values retain latent-column alignment", {
  L <- cbind(a = 1:3, b = 4:6)
  expect_identical(.mvn_resolve_lambda_cols(c(TRUE, FALSE), 2), 1L)
  expect_error(.mvn_resolve_lambda_cols(TRUE, 2), "length")
  expect_error(.mvn_resolve_lambda_cols(c(1, 3), 2), "out of range")
  expect_error(.mvn_resolve_lambda(-.1, L, diag(3), 1:2), "scalar")
  expect_error(.mvn_resolve_lambda(c(.3, 2), L, diag(3), 1:2), "entries")
  expect_error(.mvn_resolve_lambda("unknown", L, diag(3), 1:2), "fixed_1")
  v <- structure(c(.2, .8), lambda_block = .4)
  out <- .mvn_resolve_lambda(v, L, diag(3), 1:2)
  expect_equal(out$lambda_vec, c(a = .2, b = .8))
  expect_equal(out$lambda_block, .4)
  expect_equal(.mvn_resolve_lambda("estimate", L, diag(3), integer())$lambda_block, 1)
  expect_equal(.mvn_gls_mean_at_lambda(c(NA, NA), diag(2), 1), 0)
  b <- bayes_lambda_for_col(1:5, diag(5))
  expect_equal(b$lambda_grid, 1)
  expect_equal(b$weights, 1)
  expect_equal(b$lambda_post_entropy, 0)
})

test_that("completing logical values respects original row alignment", {
  original <- data.frame(species = c("b", "a", "c"),
                         flag = c(TRUE, NA, NA), extra = c(1, 2, 3))
  pred <- data.frame(flag = factor(c("yes", "1", "FALSE")))
  got <- build_completed(original, pred, species_col = "species",
                         input_row_order = c(2, 1, 3))
  expect_identical(got$completed$flag, c(TRUE, TRUE, FALSE))
  expect_identical(got$completed$extra, original$extra)
  expect_identical(unname(got$imputed_mask[, "flag"]), c(FALSE, TRUE, TRUE))
  pred$flag <- c(FALSE, TRUE, TRUE)
  got <- build_completed(original, pred, species_col = "species")
  expect_identical(got$completed$flag, c(TRUE, TRUE, TRUE))
})

test_that("active sampling handles empty observations and failed factorization", {
  sim <- matrix(1, 3, 3); diag(sim) <- 0
  expect_length(lp_entropy_reduction_binary(c(NA, NA, NA), sim), 0)
  expect_length(lp_entropy_reduction_binary(c(0, 1, 1), sim), 0)
  expect_error(suggest_next_observation(list()), "pigauto_result")
  y <- setNames(c(1:5, NA), letters[1:6])
  R <- diag(c(rep(-1, 5), 1))
  expect_warning(v <- bm_variance_reduction(y, R), "Cholesky failed")
  expect_equal(v, c(f = 0))
})

test_that("latent utilities enforce shapes and preserve probability normalization", {
  expect_error(build_phylo_graph(list()), "phylo")
  expect_error(spectral_features_from_D(as.matrix(dist(1:4)), 4, .5), "too large")
  expect_equal(estep_liability(list(type = "zi_count"), c(1, 2), 0, 1),
               list(mean = c(1, 2), var = c(0, 0)))
  expect_error(estep_liability(list(type = "unknown"), 1, 0, 1), "unsupported")
  expect_error(.phylo_signal_training_data(list(X_scaled = 1:3), list(val_idx = 1)), "matrix")
  expect_error(.phylo_signal_training_data(list(X_scaled = matrix(1, 2, 2)),
                                           list(val_idx = 5)), "valid linear")
  expect_null(mondrian_cell_scores(list(), NULL, list(), 2))
  p <- matrix(c(.9, .9, NA, .1, .2, .7), 2, byrow = TRUE)
  decoded <- exp(decode_ovr_categorical(p))
  expect_equal(rowSums(decoded), c(1, 1))
  expect_true(all(is.finite(decoded)))
})

test_that("GNN-off prediction uses the correct baseline and rejects malformed fits", {
  set.seed(123)
  tree <- ape::rtree(20)
  traits <- data.frame(x = seq_len(20) + sin(seq_len(20)), row.names = tree$tip.label)
  traits$x[c(2, 7)] <- NA_real_
  expect_warning(
    result <- impute(traits, tree, gnn = FALSE, verbose = FALSE, seed = 42,
                     safety_floor = FALSE, phylo_signal_gate = FALSE),
    "Small validation set"
  )
  fit <- result$fit
  prod <- predict(fit)
  held <- c(fit$splits$val_idx, fit$splits$test_idx)
  expect_gt(length(held), 0)
  eval <- predict(fit, evaluation = TRUE)
  expect_equal(unname(prod$imputed_latent), unname(fit$baseline_full$mu))
  expect_equal(unname(eval$imputed_latent), unname(fit$baseline$mu))
  expect_true(all(held %in% eval$evaluation_safe_idx))
  test_idx <- fit$splits$test_idx
  expect_gt(length(test_idx), 0L)
  expect_true(any(abs(fit$baseline_full$mu[test_idx] - fit$baseline$mu[test_idx]) > 1e-8))
  truth <- fit$X_scaled[test_idx]
  masked_rmse <- sqrt(mean((truth - fit$baseline$mu[test_idx])^2))
  full_rmse <- sqrt(mean((truth - fit$baseline_full$mu[test_idx])^2))
  expect_gt(abs(masked_rmse - full_rmse), 1e-8)
  scored <- evaluate(fit)
  baseline_rmse <- scored$value[scored$method == "baseline" & scored$trait == "x" & scored$metric == "rmse"]
  expect_length(baseline_rmse, 1L)
  expect_equal(baseline_rmse, masked_rmse, tolerance = 1e-10)
  expect_equal(result$completed$x[!is.na(traits$x)], traits$x[!is.na(traits$x)])
  expect_true(all(fit$r_cal_gnn == 0))
  expect_error(predict(fit, clamp_factor = .5), "clamp_factor")
  expect_error(predict(fit, pmm_K = 0), "pmm_K")
  expect_error(predict(fit, newdata = traits), "newdata support")
  expect_error(predict(fit, .mask_observed_idx = 0), "valid linear")
  bad <- fit; bad$r_cal_bm <- c(1, 1)
  expect_error(predict(bad), "r_cal_bm.*length")
  bad <- fit; bad$r_cal_bm <- NA_real_
  expect_error(predict(bad), "r_cal_bm.*finite")
  bad <- fit; bad$baseline_full$se <- matrix(0, 1, 1)
  expect_error(predict(bad), "same dimensions")
  bad <- fit; bad$baseline_full$mu <- NULL
  expect_error(predict(bad), "must contain")
  bad <- fit; bad$calibrated_gates <- NULL
  expect_error(predict(bad), "require calibrated gates")
  expect_output(print(prod), "Imputations:")
  legacy <- fit; legacy$X_scaled <- NULL
  expect_error(predict(legacy, match_observed = "pmm"), "older pigauto version")
  legacy$trait_map <- NULL
  legacy$norm <- list(sds = 1, means = 0, log_transform = FALSE)
  legacy_pred <- predict(legacy)
  expect_equal(unname(legacy_pred$imputed), unname(fit$baseline_full$mu))
  expect_error(.pigauto_assert_predictable_state(list(model_state = list(as.raw(1)))), "serialised")
  expect_error(.pigauto_assert_predictable_state(list(model_state = list(1))), "save_pigauto")
})

test_that("sparse conditional prediction agrees with dense Gaussian conditioning", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  h <- build_henderson_S_inv(tree)
  R <- stats::cov2cor(ape::vcv.phylo(tree))[tree$tip.label, tree$tip.label]
  y <- c(1, NA, -1, .5)
  p <- henderson_bm_predict(y, h, sigma2 = 2, cor_scale = TRUE)
  obs <- which(!is.na(y)); miss <- which(is.na(y))
  expected <- R[miss, obs, drop = FALSE] %*% solve(R[obs, obs], y[obs])
  variance <- 2 * (R[miss, miss] - R[miss, obs, drop = FALSE] %*%
                     solve(R[obs, obs], R[obs, miss, drop = FALSE]))
  expect_equal(p$mu[miss], as.numeric(expected), tolerance = 1e-6)
  expect_equal(p$se[miss], sqrt(as.numeric(variance)), tolerance = 1e-6)
  expect_equal(p$mu[obs], y[obs])
  expect_equal(p$se[obs], rep(0, length(obs)))
  A <- ape::vcv.phylo(tree)[tree$tip.label, tree$tip.label]
  raw <- henderson_bm_predict(y, h, sigma2 = -1, cor_scale = FALSE)
  raw_mean <- A[miss, obs, drop = FALSE] %*% solve(A[obs, obs], y[obs])
  raw_var <- A[miss, miss] - A[miss, obs, drop = FALSE] %*%
    solve(A[obs, obs], A[obs, miss, drop = FALSE])
  expect_equal(raw$mu[miss], as.numeric(raw_mean), tolerance = 1e-6)
  expect_equal(raw$se[miss], sqrt(as.numeric(raw_var)), tolerance = 1e-6)
  fallback <- henderson_bm_predict(rep(NA_real_, 4), h)
  expect_equal(fallback$mu, rep(0, 4))
  expect_equal(fallback$se, rep(1, 4))
  bad <- h; bad$tip_sqrt_d <- NULL
  expect_error(henderson_bm_predict(y, bad, cor_scale = TRUE), "tip_sqrt_d")
  expect_error(henderson_R_inv_apply(rep(1, 4), bad, cor_scale = TRUE), "tip_sqrt_d")
  L <- matrix(y, ncol = 1, dimnames = list(tree$tip.label, "x"))
  exact <- exact_conditional_mvn(L, matrix(2, 1), h, cor_scale = TRUE)
  expect_equal(exact$mu[miss, 1], as.numeric(expected), tolerance = 1e-6)
  expect_equal(exact$var[miss, 1], as.numeric(variance), tolerance = 1e-6)
  expect_null(exact_conditional_mvn(L, matrix(2, 1), NULL))
  expect_null(exact_conditional_mvn(matrix(NA_real_, 4, 1), matrix(2, 1), h))
  expect_length(.ec_inv_diag(NULL, 1, integer()), 0)
  expect_error(build_henderson_S_inv(tree, tip_depths = 1:3), "one value per tip")
})

test_that("BM fallback keeps observed cells and reports sample uncertainty", {
  y <- c(1:5, NA_real_)
  broken <- diag(c(rep(-1, 5), 1))
  expect_warning(out <- bm_impute_col(y, broken), "Cholesky.*failed")
  expect_equal(out$mu, c(1:5, 3))
  expect_equal(out$se, c(rep(0, 5), sd(1:5)))
  expect_error(bm_impute_col_with_cov(y, matrix(1, 6, 1), diag(6), lambda = -1),
               "lambda")
})

test_that("baseline and scoring validation reject unsupported inputs", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- data.frame(x = c(1, 2, 3, NA), row.names = tree$tip.label)
  dat <- preprocess_traits(x, tree)
  expect_error(fit_baseline(dat, tree, lambda_fixed = c(.5)), "named numeric")
  expect_error(fit_baseline(dat, tree, em_iterations = -1), "non-negative")
  expect_error(fit_baseline(list(), tree), "pigauto_data")
  expect_error(fit_baseline(dat, list()), "phylo")
  truth <- matrix(NA_real_, 2, 2)
  mu <- matrix(0, 2, 2)
  expect_equal(.pigauto_route_val_loss(list(type = "binary", latent_cols = 1),
                                       1:2, 1:2, truth, mu),
               list(n = 0L, loss = NA_real_))
  expect_equal(.pigauto_route_val_loss(list(type = "categorical", latent_cols = 1:2),
                                       1:2, 1:2, truth, mu),
               list(n = 0L, loss = NA_real_))
})

test_that("multiple-draw uncertainty combines ordinal and proportion terms", {
  tm <- list(list(name = "o", type = "ordinal", latent_cols = 1, sd = 2),
             list(name = "p", type = "proportion", latent_cols = 2, sd = 1, mean = 0))
  d1 <- list(imputed = data.frame(o = factor(c("low", "low"), levels = c("low", "high"), ordered = TRUE), p = c(.2, .3)))
  d2 <- list(imputed = data.frame(o = factor(c("high", "high"), levels = c("low", "high"), ordered = TRUE), p = c(.4, .5)))
  se <- compute_mc_se(list(d1, d2), tm, matrix(.5, 2, 2), matrix(0, 2, 2), c("a", "b"))
  expect_equal(unname(se[, "o"]), rep(sqrt(1 + .5), 2))
  expect_equal(unname(se[, "p"]), rep(sqrt((.5 * .25)^2 + .02), 2))
  expect_equal(dimnames(se), list(c("a", "b"), c("o", "p")))
})

test_that("species aggregation masks observations before combining ZI traits", {
  dat <- list(multi_obs = TRUE, n_obs = 4L, n_species = 2L,
              species_names = c("a", "b"), obs_to_species = NULL,
              obs_species = c("a", "a", "b", "b"),
              X_scaled = cbind(gate = c(0, 1, 1, 1), mag = c(NA, 2, 3, 5)),
              trait_map = list(list(type = "zi_count", latent_cols = 1:2)))
  s <- list(val_idx = c(2L, 6L), test_idx = integer())
  hard <- aggregate_to_species(dat, s)
  expect_equal(unname(hard$X_species[, 1]), c(0, 1))
  expect_equal(unname(hard$X_species[, 2]), c(NA_real_, 4))
  soft <- aggregate_to_species(dat, NULL, soft_aggregate = TRUE)
  expect_equal(unname(soft$X_species[, 1]), c(.5, 1))
  expect_identical(soft$is_proportion_col, c(TRUE, FALSE))
  expect_length(lp_entropy_reduction_categorical(matrix(NA_real_, 2, 2), diag(2)), 0)
})

test_that("label propagation handles all-held-out and partially observed discrete traits", {
  set.seed(41)
  tree <- ape::rtree(9)
  df <- data.frame(bin = factor(rep(c("no", "yes", "no"), 3)),
                   cat = factor(rep(c("a", "b", "c"), 3)),
                   row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree)
  bm <- pd$trait_map[[1]]$latent_cols
  cm <- pd$trait_map[[2]]$latent_cols
  allheld <- list(val_idx = seq_along(pd$X_scaled), test_idx = integer())
  out <- fit_baseline(pd, tree, splits = allheld, predict_method = "per_column")
  expect_equal(unname(out$mu[, bm]), rep(0, 9))
  expect_equal(unname(exp(out$mu[, cm])), matrix(1/3, 9, 3))
  expect_true(all(out$se[, c(bm, cm)] == 0))
  expect_identical(unname(out$path["cat"]), "label_propagation")
  holdcat <- which(row(pd$X_scaled) == 1 & col(pd$X_scaled) %in% cm)
  part <- fit_baseline(pd, tree, splits = list(val_idx = holdcat, test_idx = integer()),
                       predict_method = "per_column")
  expect_equal(unname(rowSums(exp(part$mu[, cm, drop = FALSE]))), rep(1, 9), tolerance = 1e-8)
  expect_error(fit_ovr_categorical_fits(pd, tree, "absent"), "no categorical trait")
})

test_that("low nonzero counts use the constant magnitude fallback", {
  set.seed(42)
  tree <- ape::rtree(10)
  df <- data.frame(z = as.integer(c(0,0,0,0,0,0,2,3,NA,NA)), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, trait_types = c(z = "zi_count"))
  lc <- pd$trait_map[[1]]$latent_cols
  out <- fit_baseline(pd, tree, predict_method = "per_column")
  vals <- pd$X_scaled[, lc[2]]
  expect_equal(unname(out$mu[, lc[2]]), rep(mean(vals, na.rm = TRUE), 10))
  expect_equal(unname(out$se[, lc[2]]), rep(sd(vals, na.rm = TRUE), 10))
  allheld <- list(val_idx = seq_along(pd$X_scaled), test_idx = integer())
  empty <- fit_baseline(pd, tree, splits = allheld, predict_method = "per_column")
  expect_equal(unname(empty$mu[, lc[1]]), rep(0, 10))
})

test_that("cache backfill and replacement recover graph quantities", {
  set.seed(43); tree <- ape::rtree(8)
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file))
  graph <- build_phylo_graph(tree, k_eigen = 2)
  old <- graph; old$D <- old$D_sq <- old$R_phy <- NULL
  saveRDS(old, file)
  loaded <- build_phylo_graph(tree, k_eigen = 2, cache_path = file)
  expect_equal(loaded$D, ape::cophenetic.phylo(tree))
  expect_equal(loaded$D_sq, loaded$D^2)
  expect_equal(loaded$R_phy, phylo_cor_matrix(tree))
  expect_message(other <- build_phylo_graph(tree, k_eigen = 3, cache_path = file), "dimensions mismatch")
  expect_equal(ncol(other$coords), 3)
})

test_that("active ranking distinguishes absent inputs, ineligible cells and ZI entropy", {
  set.seed(44); tree <- ape::rtree(8); spp <- tree$tip.label
  X <- cbind(cont = 1:8, binary = rep(0:1, 4),
             ca = rep(c(1, 0), 4), cb = rep(c(0, 1), 4),
             gate = rep(c(0, 1), 4), mag = c(NA, 1, NA, 2, NA, 3, NA, 4))
  rownames(X) <- spp
  maps <- list(list(name = "cont", type = "continuous", latent_cols = 1),
               list(name = "binary", type = "binary", latent_cols = 2),
               list(name = "cat", type = "categorical", latent_cols = 3:4),
               list(name = "z", type = "zi_count", latent_cols = 5:6))
  result <- structure(list(data = list(X_scaled = X, trait_map = maps, species_names = spp),
                           fit = list(graph = list(R_phy = phylo_cor_matrix(tree))), tree = tree),
                      class = "pigauto_result")
  expect_error(suggest_next_observation(result, top_n = 0), "top_n")
  absent <- result; absent$tree <- NULL; absent$fit$graph$R_phy <- NULL
  expect_error(suggest_next_observation(absent, types = "continuous"), "correlation matrix")
  absent$fit$graph$R_phy <- phylo_cor_matrix(tree)
  expect_error(suggest_next_observation(absent, types = "binary"), "tree is missing")
  no_candidates <- suggest_next_observation(result)
  expect_equal(nrow(no_candidates), 0L)
  expect_output(print(no_candidates), "No eligible candidates")
  result$data$X_scaled[8, 2:6] <- NA_real_
  z <- suggest_next_observation(result, types = "zi_count")
  expect_equal(z$species, spp[8])
  expect_equal(z$delta_var_total, 0)
  expect_true(is.finite(z$delta_entropy_total) && z$delta_entropy_total >= 0)
  expect_output(print(z), "proxy rankings")
  unclassed_by <- z; attr(unclassed_by, "by") <- NULL
  expect_output(print(unclassed_by), "by = cell")
  cat_only <- suggest_next_observation(result, types = "categorical")
  expect_equal(cat_only$species, spp[8])
  expect_true(is.finite(cat_only$delta) && cat_only$delta >= 0)
  empty_disc <- result; empty_disc$data$X_scaled[, 2:5] <- NA_real_
  empty <- suggest_next_observation(empty_disc, types = c("binary", "categorical", "zi_count"))
  expect_equal(nrow(empty), 0L)
})

test_that("phylogenetic signal dispatch respects data availability and ZI components", {
  set.seed(45); tree <- ape::rtree(12); spp <- tree$tip.label
  X <- cbind(gate = rep(c(0, 1), 6), mag = seq_len(12) + sin(seq_len(12)))
  rownames(X) <- spp
  dat <- list(X_scaled = X, species_names = spp)
  tm <- list(name = "z", type = "zi_count", latent_cols = 1:2)
  all_missing <- dat; all_missing$X_scaled[,] <- NA_real_
  expect_true(is.na(.phylo_signal_one_trait(tm, all_missing, tree, "lambda", 5)))
  multi <- list(X_scaled = matrix(NA_real_, 24, 1),
                obs_to_species = rep(1:12, each = 2), species_names = spp)
  expect_true(is.na(.phylo_signal_lambda_continuous(list(latent_cols = 1), multi, tree, "lambda", 5)))
  if (requireNamespace("phytools", quietly = TRUE)) {
    gate <- .phylo_signal_lambda_continuous(list(latent_cols = 1), dat, tree, "lambda", 5)
    mag <- .phylo_signal_lambda_continuous(list(latent_cols = 2), dat, tree, "lambda", 5)
    expect_equal(.phylo_signal_one_trait(tm, dat, tree, "lambda", 5), max(gate, mag))
    only_gate <- dat; only_gate$X_scaled[, 2] <- NA_real_
    only_mag <- dat; only_mag$X_scaled[, 1] <- NA_real_
    expect_equal(.phylo_signal_one_trait(tm, only_gate, tree, "lambda", 5), gate)
    expect_equal(.phylo_signal_one_trait(tm, only_mag, tree, "lambda", 5), mag)
    no_names <- dat; no_names$species_names <- NULL
    expect_equal(.phylo_signal_lambda_continuous(list(latent_cols = 2), no_names, tree, "lambda", 5), mag)
    mp <- list(type = "multi_proportion", latent_cols = 2)
    expect_equal(.phylo_signal_one_trait(mp, dat, tree, "lambda", 5), mag)
  }
})

test_that("joint solver validates inputs and uses dense covariance without a tree", {
  L <- cbind(a = c(1, 2, NA, 4, 5, 6), b = c(2, NA, 3, 5, 7, 9))
  rownames(L) <- paste0("s", 1:6)
  R <- diag(6); dimnames(R) <- list(rownames(L), rownames(L))
  expect_error(fit_mvn_bm_inhouse(L), "either tree or R")
  expect_error(fit_mvn_bm_inhouse(unname(L), R = R), "rownames")
  out <- fit_mvn_bm_inhouse(L, R = R, predict_method = "per_column")
  expect_equal(dim(out$anc_recon), dim(L))
  expect_equal(unname(out$anc_recon[!is.na(L)]), unname(L[!is.na(L)]))
  expect_true(all(is.finite(out$anc_recon)))
  expect_true(all(diag(out$pars$phylocov) > 0))
  notified_before <- .pigauto_exact_fallback_env$notified
  withr::defer(.pigauto_exact_fallback_env$notified <- notified_before)
  .pigauto_exact_fallback_reset()
  expect_message(fallback <- fit_mvn_bm_inhouse(L, R = R), "no Henderson")
  expect_identical(fallback$predict_method_used, "per_column")
  sigma <- .mvn_sigma_ml(matrix(c(1, 2, 3, 4), 2), matrix(0, 2, 2), diag(2))
  expect_equal(sigma, crossprod(matrix(c(1, 2, 3, 4), 2)) / 2)
  pd <- .mvn_ensure_pd(matrix(c(1, .1, .3, 1), 2))
  expect_true(isSymmetric(pd))
  expect_true(all(eigen(pd, symmetric = TRUE)$values > 0))
})

test_that("multi-observation fallback aggregates discrete traits by species", {
  set.seed(46); tree <- ape::rtree(8)
  df <- data.frame(species = rep(tree$tip.label, each = 2),
                   binary = factor(rep(c("yes", "no"), 8)),
                   cat = factor(rep(c("a", "b", "c", "a"), 4)),
                   z = as.integer(rep(c(0, 2), 8)))
  pd <- preprocess_traits(df, tree, species_col = "species", trait_types = c(z = "zi_count"))
  binary <- pd$trait_map[[1]]$latent_cols
  catcols <- pd$trait_map[[2]]$latent_cols
  gate <- pd$trait_map[[3]]$latent_cols[1]
  # Leave one categorical observation and all binary/ZI observations available;
  # the OVR fit sees too little category information and falls back to LP.
  cat_idx <- which(col(pd$X_scaled) %in% catcols & row(pd$X_scaled) > 1)
  out <- fit_baseline(pd, tree, splits = list(val_idx = cat_idx, test_idx = integer()),
                      predict_method = "per_column")
  expect_true(all(is.finite(out$mu[, c(binary, catcols, gate)])))
  expect_equal(unname(rowSums(exp(out$mu[, catcols, drop = FALSE]))), rep(1, 8))
  expect_true(all(out$se[, c(binary, catcols, gate)] == 0))
  expect_true(all(stats::plogis(out$mu[, gate]) >= .01))
})

test_that("result printing labels each held-out metric", {
  evaluation <- data.frame(split = rep("test", 4), trait = c("x", "n", "o", "b"),
                           type = c("continuous", "count", "ordinal", "binary"),
                           rmse = c(1, 2, NA, NA), pearson_r = c(.4, NA, NA, NA),
                           spearman_rho = c(NA, NA, .3, NA), accuracy = c(NA, NA, NA, .75))
  res <- structure(list(fit = list(model_config = list(gnn = FALSE)),
                        data = list(species_names = "a", trait_map = list(list(name = "x", type = "continuous"))),
                        imputed_mask = matrix(TRUE, 1, 1, dimnames = list(NULL, "x")),
                        evaluation = evaluation), class = "pigauto_result")
  printed <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(printed, "GNN: off")
  expect_match(printed, "Held-out test-set metrics")
  expect_match(printed, "RMSE")
  expect_match(printed, "0.4000")
  expect_match(printed, "rho")
  expect_match(printed, "accuracy")
  expect_match(printed, "does not authorize inference")
})

test_that("joint E-step pools uncertainty without changing observed cells", {
  L <- rbind(c(1, 2), c(3, 4))
  mask <- matrix(c(TRUE, FALSE, TRUE, TRUE), 2)
  variances <- matrix(c(0, 2, 0, 0), 2)
  out <- .mvn_estep_refine(mask, L, variances, diag(2), eps = 0, refine_variance = "pooled")
  expect_equal(out$L_hat[mask], L[mask])
  expect_equal(out$L_var[2, 1], 2/3)
  expect_equal(out$L_hat[2, 1], 1)
  sparse_L <- cbind(a = c(1, NA, NA), b = c(1, 2, 3))
  init <- .mvn_init_per_column(sparse_L, diag(3))
  expect_equal(unname(init$L_hat[, 1]), c(1, 0, 0))
  expect_equal(unname(init$L_var[, 1]), c(0, 1, 1))
})

test_that("Pagel likelihood fallback is explicit for insufficient observations", {
  c1 <- build_pagel_nll_cache(c(1, NA, NA), diag(3))
  expect_equal(c1$n_o, 1L)
  expect_equal(c1$nll(.5), .Machine$double.xmax)
  expect_true(is.na(c1$mu_at(.5)))
  c2 <- build_pagel_nll_cache(c(1, 2), diag(2), X = diag(2))
  expect_equal(c2$nll(.5), .Machine$double.xmax)
  expect_null(c2$mu_at)
  expect_error(transform_tree_pagel(list(), .5), "phylo")
})

test_that("multi-observation discrete baselines retain species-level evidence", {
  set.seed(47); tree <- ape::rtree(8); n <- 16
  df <- data.frame(species = rep(tree$tip.label, each = 2),
                   b = factor(rep(c("yes", "no", "yes", "yes"), 4)),
                   z = as.integer(rep(c(0, 2, 1, 0), 4)),
                   p = seq(.1, .8, length.out = n))
  df$q <- 1 - df$p
  pd <- preprocess_traits(df, tree, species_col = "species", trait_types = c(z = "zi_count"),
                          multi_proportion_groups = list(comp = c("p", "q")))
  out <- fit_baseline(pd, tree, predict_method = "per_column")
  b <- pd$trait_map[[which(vapply(pd$trait_map, function(tm) tm$name == "b", logical(1)))]]$latent_cols
  z <- pd$trait_map[[which(vapply(pd$trait_map, function(tm) tm$name == "z", logical(1)))]]$latent_cols[1]
  expect_equal(dim(out$mu)[1], 8L)
  expect_true(all(is.finite(out$mu[, c(b, z)])))
  expect_true(all(stats::plogis(out$mu[, c(b, z)]) >= .01))
  expect_true(all(stats::plogis(out$mu[, c(b, z)]) <= .99))
})

test_that("EM and joint ZI fits retain latent shape and finite covariance", {
  set.seed(48); tree <- ape::rtree(12)
  df <- data.frame(x = seq_len(12) + sin(seq_len(12)),
                   z = as.integer(rep(c(0, 1, 2), 4)),
                   b = factor(rep(c("no", "yes"), 6)),
                   cat = factor(rep(c("a", "b", "c"), 4)), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, trait_types = c(z = "zi_count"))
  joint <- fit_joint_mvn_baseline(pd, tree, splits = NULL, predict_method = "per_column")
  expect_equal(dim(joint$mu), dim(pd$X_scaled))
  expect_true(all(is.finite(joint$mu[, pd$trait_map[[2]]$latent_cols[2]])))
  empty <- list(val_idx = seq_along(pd$X_scaled), test_idx = integer())
  no_fit <- fit_joint_threshold_baseline_em(pd, tree, splits = empty,
                                             em_iterations = 2, predict_method = "per_column")
  expect_null(no_fit$phylopars_fit)
  expect_equal(no_fit$em_state$iterations_run, 1L)
  em <- fit_joint_threshold_baseline_em(pd, tree, splits = NULL, em_iterations = 2,
                                         em_tol = 1000, predict_method = "per_column")
  expect_true(isTRUE(em$em_state$converged))
  expect_equal(em$em_state$iterations_run, 2L)
  ovr <- fit_ovr_categorical_fits_em(pd, tree, "cat", em_iterations = 2,
                                     em_tol = 1000, predict_method = "per_column")
  expect_equal(dim(ovr), c(12L, 3L))
  expect_true(all(is.finite(ovr)))
  expect_true(isTRUE(attr(ovr, "em_state")$converged))
})

test_that("refinement with dense precision has a finite one-step outcome", {
  L <- cbind(a = c(1, 2, NA, 4, 5, 6), b = c(2, NA, 3, 5, 7, 9))
  rownames(L) <- paste0("s", 1:6)
  R <- diag(6); dimnames(R) <- list(rownames(L), rownames(L))
  out <- fit_mvn_bm_inhouse(L, R = R, max_iter = 1, tol = 1e6, predict_method = "per_column")
  expect_equal(out$n_iter, 1L)
  expect_true(out$converged)
  expect_true(all(is.finite(out$anc_recon)))
  expect_equal(unname(out$anc_recon[!is.na(L)]), unname(L[!is.na(L)]))
})

test_that("covariate BM fallback reports mean and standard deviation", {
  y <- c(1:10, NA_real_); X <- matrix(1, 11, 1)
  R <- diag(c(rep(-1, 10), 1))
  expect_warning(out <- bm_impute_col_with_cov(y, X, R), "Cholesky failed")
  expect_equal(out$mu, c(1:10, mean(1:10)))
  expect_equal(out$se, c(rep(0, 10), sd(1:10)))
  expect_true(all(is.na(out$beta_hat)))
})

test_that("OU fallback announces BM and partial lambda vectors default to one", {
  set.seed(49); tree <- ape::rtree(10)
  df <- data.frame(x = c(as.numeric(1:9), NA), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree)
  expect_message(ou <- fit_baseline(pd, tree, model = "OU", predict_method = "per_column"), "using BM")
  bm <- fit_baseline(pd, tree, lambda_fixed = c(other = .4), predict_method = "per_column")
  expect_equal(unname(ou$mu), unname(bm$mu))
  expect_equal(unname(bm$lambda_per_trait), 1)
})

test_that("legacy auto route defaults to exact when held-out evidence is absent", {
  dat <- list(trait_map = list(list(name = "x", type = "continuous", latent_cols = 1)),
              X_scaled = matrix(1:3, 3, 1))
  fit <- list(mu = matrix(0, 3, 1))
  route <- .pigauto_choose_predict_route(dat, list(val_idx = integer()), fit, fit)
  expect_identical(route, c(x = "exact"))
})

test_that("Henderson tip-only algebra matches direct precision calculations", {
  Q <- Matrix::Matrix(matrix(c(2, -.25, -.25, 3), 2), sparse = TRUE)
  h <- list(Q = Q, tip_idx = 1:2, int_idx = integer(), tip_sqrt_d = c(2, 3))
  b <- cbind(c(1, 4), c(-2, 3))
  expect_equal(henderson_R_inv_apply(b, h), as.matrix(Q %*% b))
  expected_cor <- diag(c(2, 3)) %*% as.matrix(Q) %*% diag(c(2, 3)) %*% b
  expect_equal(henderson_R_inv_apply(b, h, cor_scale = TRUE), expected_cor)
  expect_equal(henderson_log_det_R(h), -as.numeric(determinant(as.matrix(Q), logarithm = TRUE)$modulus))
})

test_that("observed likelihood omits fully missing rows and failed GLS uses the mean", {
  L <- rbind(c(1, 2), c(3, NA), c(NA, NA))
  par <- .mvn_chol_to_par(diag(2))
  expect_equal(.mvn_obs_nll(par, L, 2), .mvn_obs_nll(par, L[1:2, ], 2))
  expect_equal(.mvn_obs_nll(par, L[3, , drop = FALSE], 2), 0)
  bad_R <- matrix(c(1, 2, 2, 1), 2)
  expect_equal(.mvn_gls_mean_at_lambda(c(2, 6), bad_R, 1), 4)
  expect_equal(.mvn_gls_mean_at_lambda(c(NA, NA), diag(2), 1), 0)
})

test_that("legacy observation aggregation recovers species indices from row names", {
  X <- matrix(c(2, 4, 6, 8), 4, 1,
              dimnames = list(c("a", "a", "b", "b"), "x"))
  d <- list(X_scaled = X, n_obs = 4L, n_species = 2L, species_names = c("a", "b"),
            multi_obs = TRUE, trait_map = list(list(type = "continuous", latent_cols = 1L)))
  got <- aggregate_to_species(d)
  expect_equal(unname(got$X_species), matrix(c(3, 7), 2, 1))
  expect_identical(rownames(got$X_species), c("a", "b"))
  expect_identical(d$X_scaled, X)
  dimnames(d$X_scaled) <- NULL
  expect_error(aggregate_to_species(d), "Observation-to-species index length")
})

test_that("decoding respects missing caps and reports probability traits", {
  latent <- cbind(count = log1p(c(3, 20)), binary = c(-2, 2))
  count <- list(name = "count", type = "count", latent_cols = 1L,
                mean = 0, sd = 1)
  binary <- list(name = "binary", type = "binary", latent_cols = 2L,
                 levels = c("no", "yes"))
  decode <- function(tm) decode_from_latent(latent, list(tm, binary),
    c("a", "b"), clamp_outliers = TRUE, clamp_factor = 2)
  a <- decode(count)
  expect_equal(a$imputed$count, c(3L, 20L))
  count$obs_max <- Inf
  expect_equal(decode(count)$imputed$count, c(3L, 20L))
  count$obs_max <- 4
  capped <- decode(count)
  expect_equal(capped$imputed$count, c(3L, 8L))
  expect_equal(unname(capped$probabilities$binary), plogis(c(-2, 2)))
  pred <- structure(c(capped, list(species_names = c("a", "b"),
    trait_names = c("count", "binary"), n_imputations = 1L)), class = "pigauto_pred")
  expect_output(print(pred), "Probability traits: binary", fixed = TRUE)
})

test_that("fit print reports ungated phylogenetic signal", {
  fit <- structure(list(model_config = list(gnn = FALSE, hidden_dim = 8L, k_eigen = 2L),
    species_names = c("a", "b"), trait_names = "x", val_rmse = NA_real_,
    test_rmse = NA_real_, phylo_gate_triggered = c(x = FALSE),
    phylo_signal_per_trait = c(x = .8), phylo_signal_threshold = .2,
    phylo_signal_method = "lambda"), class = "pigauto_fit")
  expect_output(print(fit), "all 1 traits have lambda >= 0.20", fixed = TRUE)
  fit$phylo_signal_method <- "K"
  expect_output(print(fit), "all 1 traits have K >= 0.20", fixed = TRUE)
})

test_that("lambda likelihood refuses singular covariance and exact zero residuals", {
  R <- matrix(1, 3, 3)
  ordinary <- build_pagel_nll_cache(c(1, 2, 4), R, nugget = 0)
  expect_identical(ordinary$nll(1), .Machine$double.xmax)
  expect_true(is.na(ordinary$mu_at(1)))
  covariate <- build_pagel_nll_cache(c(1, 2, 4), R, nugget = 0,
                                    X = matrix(1, 3, 1))
  expect_identical(covariate$nll(1), .Machine$double.xmax)
  exact <- build_pagel_nll_cache(rep(0, 3), diag(3), X = matrix(1, 3, 1))
  expect_identical(exact$nll(.5), .Machine$double.xmax)
})

test_that("Fisher-ML refinement retains observed cells on a small dense correlation", {
  L <- cbind(a = c(1, 2, NA, 4, 5, 6), b = c(2, NA, 3, 5, 7, 9))
  rownames(L) <- paste0("s", seq_len(nrow(L)))
  R <- diag(6); dimnames(R) <- list(rownames(L), rownames(L))
  out <- fit_mvn_bm_inhouse(L, R = R, sigma_method = "fisher_ml",
                            max_iter = 1L, tol = 1e6, predict_method = "per_column")
  expect_true(out$converged)
  expect_true(all(is.finite(out$anc_recon)))
  expect_equal(unname(out$anc_recon[!is.na(L)]), unname(L[!is.na(L)]))
  expect_equal(out$pars$phylocov, t(out$pars$phylocov), tolerance = 1e-10)
  expect_true(all(eigen(out$pars$phylocov, symmetric = TRUE)$values > 0))
})

test_that("joint solver handles healthy, constant, and collinear continuous traits", {
  skip_if_not_installed("Rphylopars")
  set.seed(241)
  tree <- ape::rtree(12)
  i <- seq_len(12)
  x <- as.numeric(i / 10 + sin(i) / 4)
  for (kind in c("healthy", "constant", "collinear")) {
    y <- switch(kind, healthy = cos(i / 2) + i / 20,
                constant = rep(3, 12), collinear = 2 * x + 3)
    traits <- data.frame(x = x, y = y, row.names = tree$tip.label)
    traits$x[c(2, 7)] <- NA_real_
    traits$y[if (kind == "collinear") c(2, 7) else c(3, 9)] <- NA_real_
    d <- preprocess_traits(traits, tree, trait_types = c(x = "continuous", y = "continuous"),
                           log_transform = FALSE)
    expect_identical(vapply(d$trait_map, `[[`, character(1), "type"),
                     c(x = "continuous", y = "continuous"), info = kind)
    if (kind == "constant") {
      expect_true(all(d$X_scaled[!is.na(d$X_scaled[, "y"]), "y"] == 0))
    } else if (kind == "collinear") {
      observed_both <- complete.cases(d$X_scaled)
      expect_equal(unname(d$X_scaled[observed_both, "x"]),
                   unname(d$X_scaled[observed_both, "y"]), tolerance = 1e-12)
    }
    run_capturing_fallback <- function(expr) {
      fallback <- character()
      value <- withCallingHandlers(expr, warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl('fit_joint_solver: joint_solver = "rphylopars" failed', msg, fixed = TRUE) &&
            grepl("falling back to the in-house solver", msg, fixed = TRUE)) {
          fallback <<- c(fallback, msg)
          invokeRestart("muffleWarning")
        }
      })
      list(value = value, fallback = fallback)
    }
    joint <- run_capturing_fallback(fit_joint_solver(d$X_scaled, tree,
      joint_solver = "rphylopars", predict_method = "per_column", lambda = "fixed_1"))
    out <- joint$value
    tips <- out$anc_recon[match(tree$tip.label, rownames(out$anc_recon)), , drop = FALSE]
    expect_equal(dim(tips), c(12L, 2L), info = kind)
    expect_identical(rownames(tips), tree$tip.label, info = kind)
    expect_identical(colnames(tips), colnames(d$X_scaled), info = kind)
    expect_true(all(is.finite(tips)), info = kind)
    expect_equal(unname(tips[!is.na(d$X_scaled)]),
                 unname(d$X_scaled[!is.na(d$X_scaled)]), tolerance = 1e-6, info = kind)
    if (length(joint$fallback)) {
      expect_identical(out$predict_method_used, "per_column", info = kind)
    }
    dispatched <- run_capturing_fallback(fit_baseline(d, tree, joint_solver = "rphylopars",
      predict_method = "per_column", lambda_mode = "fixed_1"))
    base <- dispatched$value
    expect_equal(dim(base$mu), c(12L, 2L), info = kind)
    expect_identical(dimnames(base$mu), dimnames(d$X_scaled), info = kind)
    expect_identical(dimnames(base$se), dimnames(base$mu), info = kind)
    expect_true(all(is.finite(base$mu)) && all(is.finite(base$se)), info = kind)
    expect_equal(unname(base$mu[!is.na(d$X_scaled)]),
                 unname(d$X_scaled[!is.na(d$X_scaled)]), tolerance = 1e-6, info = kind)
    expect_identical(unname(base$path), c("joint_mvn", "joint_mvn"), info = kind)
    if (length(dispatched$fallback)) {
      expect_identical(base$predict_method_used, "per_column", info = kind)
    }
  }
})

test_that("inconsistent legacy GNN configuration reports a model loading error", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  traits$x[2] <- NA_real_
  d <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                         log_transform = FALSE)
  fit <- fit_pigauto(d, tree, gnn = TRUE, epochs = 2L, eval_every = 1L,
    hidden_dim = 8L, n_gnn_layers = 1L, corruption_rate = .8, corruption_ramp = 0L,
    dropout = 0, edge_dropout = 0, use_transformer_blocks = FALSE,
    use_attention = TRUE, verbose = FALSE, seed = 41, predict_method = "per_column",
    refine_steps = 1L)
  expect_true(all(is.finite(predict(fit, return_se = FALSE)$imputed_latent)))
  inconsistent <- fit
  inconsistent$model_config$hidden_dim <- 16L
  expect_error(predict(inconsistent, return_se = FALSE), "save_pigauto.*load_pigauto")
})

test_that("ordinal baseline keeps threshold route without ordinal validation cells", {
  set.seed(28)
  tree <- ape::rtree(20)
  i <- seq_len(20)
  traits <- data.frame(x = sin(i / 2),
    grade = ordered(rep(c("low", "mid", "high", "mid"), 5),
                    levels = c("low", "mid", "high")),
    row.names = tree$tip.label)
  d <- preprocess_traits(traits, tree, trait_types = c(x = "continuous", grade = "ordinal"),
                         log_transform = FALSE)
  ordinal_col <- d$trait_map$grade$latent_cols
  expect_identical(ordinal_col, 2L)
  # Only the continuous column has held-out cells; ordinal has no validation data.
  splits <- list(val_idx = 2L, test_idx = 3L)
  for (method in c("per_column", "auto")) {
    out <- fit_baseline(d, tree, splits = splits, predict_method = method,
                        lambda_mode = "fixed_1", seed = 42)
    expect_identical(unname(out$ordinal_path_chosen[as.character(ordinal_col)]),
                     "threshold_joint", info = method)
    expect_identical(unname(out$path["grade"]), "threshold_joint", info = method)
    expect_identical(dimnames(out$mu), dimnames(d$X_scaled), info = method)
    expect_true(all(is.finite(out$mu)) && all(is.finite(out$se)), info = method)
    expect_true(all(out$se >= 0), info = method)
  }
})


test_that("pooling one decoded row preserves count median and ordinal mode", {
  tm <- list(list(name = "count", type = "count", latent_cols = 1L, mean = 0, sd = 1),
             list(name = "grade", type = "ordinal", latent_cols = 2L, mean = 0, sd = 1,
                  levels = c("low", "mid", "high")))
  draws <- lapply(list(c(log1p(1), 0), c(log1p(3), 1), c(log1p(7), 1)),
    function(x) matrix(x, nrow = 1, dimnames = list("tip", c("count", "grade"))))
  decoded <- lapply(draws, decode_from_latent, trait_map = tm, species_names = "tip")
  pooled <- pool_imputations(decoded, draws, tm, pool_method = "mode")
  expect_equal(dim(pooled$imputed), c(1L, 2L))
  expect_identical(rownames(pooled$imputed), "tip")
  expect_identical(pooled$imputed$count, 3L)
  expect_true(is.ordered(pooled$imputed$grade))
  expect_identical(as.character(pooled$imputed$grade), "mid")
  expect_equal(pooled$latent_mean, Reduce(`+`, draws) / 3)
})

test_that("baseline draws restore an absent R random seed", {
  local({
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if
        (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(41)
    tree <- ape::rtree(12)
    traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
    traits$x[2] <- NA_real_
    data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                              log_transform = FALSE)
    fit <- fit_pigauto(data, tree, gnn = FALSE, seed = 41, verbose = FALSE,
                       predict_method = "per_column")
    rm(list = ".Random.seed", envir = .GlobalEnv)
    first <- predict(fit, n_imputations = 3L, return_se = TRUE)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    second <- predict(fit, n_imputations = 3L, return_se = TRUE)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    expect_equal(first$imputed_datasets, second$imputed_datasets)
    expect_true(all(is.finite(first$imputed_latent)))
    expect_gt(first$se_latent[2, 1], 0)
  })
})

test_that("missing optional phytools returns named unavailable signal values", {
  set.seed(52)
  tree <- ape::rtree(8)
  traits <- data.frame(x = seq_len(8) / 10, y = cos(seq_len(8)),
                       row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous", y = "continuous"),
                            log_transform = FALSE)
  fn <- compute_phylo_signal_per_trait
  env <- new.env(parent = environment(fn))
  env$requireNamespace <- function(package, ...) {
    if (identical(package, "phytools")) FALSE else base::requireNamespace(package, ...)
  }
  environment(fn) <- env
  expect_warning(out <- fn(data, tree), "requires the 'phytools' package")
  expect_identical(names(out), c("x", "y"))
  expect_identical(unname(out), c(NA_real_, NA_real_))
})

test_that("exact conditional refuses singular trait covariance without a ridge", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  h <- build_henderson_S_inv(tree)
  L <- cbind(x = c(1, NA, -1, .5), y = c(1, NA, -1, .5))
  rownames(L) <- tree$tip.label
  Sigma <- matrix(1, 2, 2)
  expect_equal(qr(Sigma)$rank, 1L)
  expect_null(exact_conditional_mvn(L, Sigma, h, eps = 0))
  regularized <- exact_conditional_mvn(L, Sigma, h, eps = 1e-6)
  expect_false(is.null(regularized))
  expect_true(all(is.finite(regularized$mu)) && all(is.finite(regularized$var)))
  expect_equal(regularized$mu[!is.na(L)], L[!is.na(L)])
  expect_true(all(regularized$var >= 0))
})

test_that("categorical OVR retains unavailable class fits when all cells are held out", {
  set.seed(51)
  tree <- ape::rtree(9)
  traits <- data.frame(cat = factor(rep(c("a", "b", "c"), 3)),
                       row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree)
  # All one-hot cells are held out, leaving no observed category to fit.
  splits <- list(val_idx = seq_along(data$X_scaled), test_idx = integer())
  out <- fit_ovr_categorical_fits_em(data, tree, "cat", splits = splits,
                                    predict_method = "per_column", em_iterations = 2L)
  expect_equal(dim(out), c(9L, 3L))
  expect_identical(rownames(out), tree$tip.label)
  expect_identical(colnames(out), colnames(data$X_scaled))
  expect_true(all(is.na(out)))
  fits <- attr(out, "per_class_fits")
  expect_length(fits, 3L)
  expect_true(all(vapply(fits, function(x) is.null(x$jt$phylopars_fit), logical(1))))
  expect_true(is.na(attr(out, "predict_method_used")))
  expect_equal(attr(out, "em_state")$sigma_cat, rep(1, 3))
  expect_identical(attr(out, "em_state")$iterations_run, 2L)
})

test_that("active composition ranking requires five observed compositions", {
  set.seed(60)
  tree <- ape::rtree(8)
  p <- seq(.15, .85, length.out = 8)
  for (n_observed in c(4L, 5L)) {
    traits <- data.frame(p = p, q = 1 - p, row.names = tree$tip.label)
    traits[seq.int(n_observed + 1L, 8L), c("p", "q")] <- NA_real_
    fit <- impute(traits, tree, multi_proportion_groups = list(comp = c("p", "q")),
                  gnn = FALSE, missing_frac = 0, verbose = FALSE, seed = 40,
                  predict_method = "per_column")
    ranked <- suggest_next_observation(fit, top_n = 20, types = "multi_proportion")
    expect_s3_class(ranked, "pigauto_active")
    expect_equal(nrow(ranked), if (n_observed == 4L) 0L else 3L)
    if (n_observed == 5L) {
      expect_setequal(ranked$species, tree$tip.label[6:8])
      expect_true(all(ranked$trait == "comp"))
      expect_true(all(is.finite(ranked$delta)) && all(ranked$delta >= 0))
      expect_equal(ranked$delta, ranked$delta_var_total)
    }
  }
})

test_that("a single-value Pagel grid retains its fixed lambda and GLS mean", {
  y <- c(1, 2, 4, 3, 5, NA_real_)
  cache <- build_pagel_nll_cache(y, diag(6))
  fit <- .pagel_lambda_from_cache(cache, lambda_grid = c(.5, .5))
  expect_identical(fit$lambda_hat, .5)
  expect_equal(fit$mu_hat, mean(y, na.rm = TRUE), tolerance = 1e-8)
  expect_equal(fit$mu_hat, cache$mu_at(.5))
})

test_that("Pagel CV skips empty folds without changing leave-one-out selection", {
  i <- seq_len(50)
  y <- sin(i / 7) + i / 20
  R <- .8^abs(outer(i, i, "-"))
  grid <- c(.1, .5, .9)
  nugget <- 1e-6
  # Independent GLS conditioning: recreate the fold assignment, then solve
  # the held-out mean directly rather than calling bm_impute_col().
  heldout_sse <- function(k) {
    set.seed(17)
    fold_id <- sample(rep(seq_len(k), length.out = length(y)))
    total <- numeric(length(grid))
    for (f in seq_len(k)) {
      held <- which(fold_id == f)
      if (!length(held)) next
      obs <- setdiff(i, held)
      for (j in seq_along(grid)) {
        lambda <- grid[j]
        R_lambda <- lambda * R
        diag(R_lambda) <- lambda * diag(R) + (1 - lambda)
        R_oo <- R_lambda[obs, obs, drop = FALSE] + diag(nugget, length(obs))
        weights <- solve(R_oo, cbind(1, y[obs]))
        mu <- sum(weights[, 2]) / sum(weights[, 1])
        prediction <- mu + as.vector(R_lambda[held, obs, drop = FALSE] %*%
                                      solve(R_oo, y[obs] - mu))
        total[j] <- total[j] + sum((prediction - y[held])^2)
      }
    }
    total
  }
  sse_50 <- heldout_sse(50L)
  sse_55 <- heldout_sse(55L)
  expect_equal(sse_55, sse_50, tolerance = 1e-10)
  expect_true(all(is.finite(sse_50)))
  expect_identical(which.min(sse_50), 3L)
  expect_gt(sort(sse_50)[2] - min(sse_50), 1)
  expected <- grid[which.min(sse_50)]
  expect_equal(cv_lambda_for_col(y, R, nugget = nugget, lambda_grid = grid,
                                  k = 50L, seed = 17), expected)
  expect_equal(cv_lambda_for_col(y, R, nugget = nugget, lambda_grid = grid,
                                  k = 55L, seed = 17), expected)
})

test_that("Pagel CV falls back when finite inputs overflow squared errors", {
  y <- (seq_len(50) - 25) * 1e155
  R <- diag(50)
  grid <- c(.1, .9)
  expect_true(all(is.finite(y)))
  expect_true(is.infinite(sum(y^2)))
  expect_identical(cv_lambda_for_col(y, R, lambda_grid = grid, k = 5L, seed = 17), 1)
  expect_true(cv_lambda_for_col(y / 1e155, R, lambda_grid = grid,
                                k = 5L, seed = 17) %in% grid)
})

test_that("joint baseline falls back when discrete traits are entirely held out", {
  set.seed(72)
  tree <- ape::rtree(12)
  i <- seq_len(12)
  traits <- data.frame(x = sin(i), b = factor(rep(c("no", "yes"), 6)),
    o = ordered(rep(c("lo", "mid", "hi"), 4), levels = c("lo", "mid", "hi")),
    row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                            log_transform = FALSE)
  held <- which(col(data$X_scaled) > 1L)
  splits <- list(val_idx = held, test_idx = integer())
  joint <- fit_joint_threshold_baseline(data, tree, splits = splits,
                                         predict_method = "per_column")
  b <- which(joint$liab_types == "binary")
  o <- which(joint$liab_types == "ordinal")
  expect_length(b, 1L)
  expect_length(o, 1L)
  expect_true(all(is.na(joint$mu_liab[, b])))
  expect_true(all(is.na(joint$mu_liab[, o])))
  base <- fit_baseline(data, tree, splits = splits, predict_method = "per_column")
  expect_identical(unname(base$path), c("threshold_joint", "label_propagation", "per_column_bm"))
  expect_identical(dimnames(base$mu), dimnames(data$X_scaled))
  expect_true(all(is.finite(base$mu)) && all(is.finite(base$se)))
  expect_equal(unname(base$mu[, "x"]), unname(data$X_scaled[, "x"]), tolerance = 1e-6)
  expect_equal(unname(base$mu[, "b"]), rep(0, 12))
  expect_equal(unname(base$mu[, "o"]), rep(0, 12))
})

test_that("Fisher covariance falls back safely on collinear real inputs", {
  L <- cbind(x = c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA),
             y = c(2, 4, 6, 8, 10, 12, 14, 16, NA, NA),
             z = c(rep(1, 8), NA, NA))
  rownames(L) <- paste0("s", seq_len(nrow(L)))
  R <- diag(10)
  dimnames(R) <- list(rownames(L), rownames(L))
  expect_equal(qr(L[complete.cases(L), ])$rank, 2L)
  for (iterations in c(0L, 1L)) {
    warnings <- character()
    fit <- withCallingHandlers(fit_mvn_bm_inhouse(L, R = R,
      sigma_method = "fisher_ml", predict_method = "per_column", max_iter = iterations),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl('sigma_method = "fisher_ml" optim() did not converge', msg, fixed = TRUE) &&
            grepl("falling back to the single_pass", msg, fixed = TRUE)) {
          warnings <<- c(warnings, msg)
          invokeRestart("muffleWarning")
        }
      })
    expect_true(all(is.finite(fit$anc_recon)), info = iterations)
    expect_true(all(is.finite(fit$pars$phylocov)), info = iterations)
    expect_equal(fit$pars$phylocov, t(fit$pars$phylocov), tolerance = 1e-8, info = iterations)
    expect_true(all(eigen(fit$pars$phylocov, symmetric = TRUE)$values > 0), info = iterations)
    expect_equal(unname(fit$anc_recon[!is.na(L)]), unname(L[!is.na(L)]), info = iterations)
    # If every attempted optimizer call falls back, compare against the real
    # single-pass solver. A future optimizer may converge on this input.
    if (length(warnings) == iterations + 1L) {
      reference <- fit_mvn_bm_inhouse(L, R = R, sigma_method = "single_pass",
                                       predict_method = "per_column", max_iter = iterations)
      expect_equal(fit$pars$phylocov, reference$pars$phylocov, tolerance = 1e-8,
                   info = iterations)
      expect_equal(fit$anc_recon, reference$anc_recon, tolerance = 1e-8,
                   info = iterations)
    }
  }
})

test_that("sparse helpers handle an unavailable Matrix namespace", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  h <- build_henderson_S_inv(tree)
  L <- matrix(c(1, NA, -1, .5), ncol = 1,
              dimnames = list(tree$tip.label, "x"))
  available <- exact_conditional_mvn(L, matrix(2, 1), h)
  expect_true(all(is.finite(available$mu)))
  # Fault injection is limited to the namespace availability check; tree,
  # precision and numerical solvers are real. Matrix is normally an Import.
  namespace_unavailable <- function(fun) {
    env <- new.env(parent = environment(fun))
    env$requireNamespace <- function(package, ...) {
      if (identical(package, "Matrix")) FALSE else base::requireNamespace(package, ...)
    }
    environment(fun) <- env
    fun
  }
  conditional <- namespace_unavailable(exact_conditional_mvn)
  henderson <- namespace_unavailable(build_henderson_S_inv)
  expect_null(conditional(L, matrix(2, 1), h))
  expect_error(henderson(tree), "requires the 'Matrix' package")
})

test_that("Pagel likelihood retains a better boundary grid point", {
  i <- seq_len(50)
  y <- sin(i / 7) + i / 20
  R <- .8^abs(outer(i, i, "-"))
  nugget <- 1e-6
  cache <- build_pagel_nll_cache(y, R, nugget = nugget)
  grid <- seq(.1, .9, length.out = 11L)
  # Dense restricted-likelihood calculation independent of the eigen cache.
  dense <- function(lambda) {
    S <- lambda * R
    diag(S) <- lambda * diag(R) + (1 - lambda) + nugget
    u <- solve(S, cbind(1, y))
    mu <- sum(u[, 2]) / sum(u[, 1])
    e <- y - mu
    variance <- as.numeric(crossprod(e, solve(S, e))) / (length(y) - 1L)
    nll <- .5 * ((length(y) - 1L) * log(variance) +
                  as.numeric(determinant(S, logarithm = TRUE)$modulus) + log(sum(u[, 1])))
    c(nll = nll, mu = mu)
  }
  direct <- vapply(grid, dense, numeric(2L))
  from_cache <- vapply(grid, cache$nll, numeric(1L))
  expect_equal(unname(from_cache), unname(direct["nll", ]), tolerance = 1e-7)
  expect_true(all(diff(direct["nll", ]) < 0))
  inside <- stats::optimize(cache$nll, interval = tail(grid, 2), tol = 1e-4)
  expect_gt(inside$objective - tail(from_cache, 1), 1e-8)
  fit <- .pagel_lambda_from_cache(cache, lambda_grid = c(.1, .9))
  expect_equal(fit$lambda_hat, .9)
  expect_equal(fit$mu_hat, unname(direct["mu", 11]), tolerance = 1e-7)
  # Alternating values favour the opposite boundary on the same correlation.
  y <- rep(c(-1, 1), 25)
  cache <- build_pagel_nll_cache(y, R, nugget = nugget)
  direct <- vapply(grid, dense, numeric(2L))
  from_cache <- vapply(grid, cache$nll, numeric(1L))
  expect_equal(unname(from_cache), unname(direct["nll", ]), tolerance = 1e-7)
  expect_true(all(diff(direct["nll", ]) > 0))
  inside <- stats::optimize(cache$nll, interval = head(grid, 2), tol = 1e-4)
  expect_gt(inside$objective - from_cache[1], 1e-8)
  fit <- .pagel_lambda_from_cache(cache, lambda_grid = c(.1, .9))
  expect_equal(fit$lambda_hat, .1)
  expect_equal(fit$mu_hat, unname(direct["mu", 1]), tolerance = 1e-7)
})

test_that("torch state from a separate raw RDS process reports the portable-save remedy", {
  skip_if_no_libtorch()
  dir <- tempfile("pigauto-external-state-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  script <- file.path(dir, "write-state.R")
  payload <- file.path(dir, "state.rds")
  writeLines(c(paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
               "saveRDS(torch::torch_tensor(c(1, 2)), commandArgs(TRUE)[1])"), script)
  rscript <- file.path(R.home("bin"), "Rscript")
  result <- system2(rscript, c("--vanilla", shQuote(script), shQuote(payload)),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  state <- readRDS(payload)
  expect_true(inherits(state, "torch_tensor"))
  expect_error(.pigauto_assert_predictable_state(list(model_state = list(state))),
               "save_pigauto.*load_pigauto")
})

test_that("tree-based Fisher refinement returns valid values on collinear traits", {
  set.seed(87)
  tree <- ape::rtree(12)
  x <- as.numeric(scale(seq_len(9)))
  L <- rbind(cbind(x = x, y = 2 * x, z = rep(1, 9)), matrix(NA_real_, 3, 3))
  rownames(L) <- tree$tip.label
  expect_equal(qr(L[complete.cases(L), ])$rank, 2L)
  warnings <- character()
  fit <- withCallingHandlers(fit_mvn_bm_inhouse(L, tree = tree, sigma_method = "fisher_ml",
    predict_method = "per_column", max_iter = 1L), warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl('sigma_method = "fisher_ml" optim() did not converge', msg, fixed = TRUE) &&
          grepl("falling back to the single_pass", msg, fixed = TRUE)) {
        warnings <<- c(warnings, msg)
        invokeRestart("muffleWarning")
      }
    })
  expect_equal(fit$n_iter, 1L)
  expect_identical(dimnames(fit$anc_recon), dimnames(L))
  expect_true(all(is.finite(fit$anc_recon)))
  expect_true(all(is.finite(fit$anc_var)) && all(fit$anc_var >= 0))
  expect_equal(unname(fit$anc_recon[!is.na(L)]), unname(L[!is.na(L)]))
  expect_equal(fit$pars$phylocov, t(fit$pars$phylocov), tolerance = 1e-8)
  expect_true(all(eigen(fit$pars$phylocov, symmetric = TRUE)$values > 0))
  # Optimizer convergence is platform-dependent. When both calls fall back,
  # compare the result with the real tree-based single-pass implementation.
  if (length(warnings) == 2L) {
    ref <- fit_mvn_bm_inhouse(L, tree = tree, sigma_method = "single_pass",
                              predict_method = "per_column", max_iter = 1L)
    expect_equal(fit$pars$phylocov, ref$pars$phylocov, tolerance = 1e-8)
    expect_equal(fit$anc_recon, ref$anc_recon, tolerance = 1e-8)
  }
})

test_that("large-tree memory warning precedes graph allocation", {
  tree <- ape::stree(10001L, type = "star")
  expect_s3_class(tree, "phylo")
  expect_length(tree$tip.label, 10001L)
  # If the warning disappears or moves past cache loading, this invalid
  # cache is a safety fuse: readRDS errors before dense graph construction.
  cache <- tempfile(fileext = ".rds")
  on.exit(unlink(cache), add = TRUE)
  writeLines("not an RDS file", cache)
  warning_text <- NULL
  stopped <- tryCatch(
    withCallingHandlers({
      build_phylo_graph(tree, cache_path = cache)
      FALSE
    }, warning = function(w) {
      warning_text <<- conditionMessage(w)
      stop(errorCondition("stop before graph allocation", class = "pigauto_warning_probe"))
    }),
    pigauto_warning_probe = function(e) TRUE
  )
  expect_true(stopped)
  expect_match(warning_text, "tree has 10001 tips", fixed = TRUE)
  expect_match(warning_text, "cophenetic() requires O(n^2) memory (~800 MB)", fixed = TRUE)
  expect_match(warning_text, "Consider subsetting your tree", fixed = TRUE)
})
