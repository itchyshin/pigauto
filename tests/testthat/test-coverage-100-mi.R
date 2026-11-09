# Focused MI contracts: these fixtures do not establish inferential calibration.
mi_cov_data <- function() {
  x <- seq(-2, 2, length.out = 24)
  d <- data.frame(y = 1 + 0.5 * x + sin(seq_along(x)) / 3,
                  x = x, z = cos(seq_along(x)))
  d$x[c(3, 8, 15)] <- NA_real_
  d
}

mi_cov_tree_data <- function() {
  set.seed(104)
  tree <- ape::rtree(8)
  d <- data.frame(x = c(1, 2, 4, 3, 7, 5, 8, 6),
                  y = c(2, 5, 3, 6, 4, 8, 7, 9),
                  row.names = tree$tip.label)
  list(tree = tree, data = d)
}

test_that("analysis-aware input guards reject unsuitable data", {
  d <- mi_cov_data()
  call <- function(data = d, formula = y ~ x + z, missing = "x", ...)
    multi_impute_analysis(data, formula, missing, model = "lm", m = 2, ...)
  expect_error(call(d[1, ]), "at least two rows")
  expect_error(call(formula = ~ x), "two-sided")
  expect_error(call(missing = "absent"), "exactly one column")
  complete <- d; complete$x[is.na(complete$x)] <- 0
  expect_error(call(complete), "at least one missing")
  nonfinite <- d; nonfinite$x[1] <- Inf
  expect_error(call(nonfinite), "Observed values.*finite")
  expect_error(call(auxiliary = c("z", "z")), "unique")
  expect_error(call(auxiliary = "absent"), "Unknown")
  expect_error(call(auxiliary = "x"), "cannot also")
  expect_error(call(formula = y ~ x + unknown), "additive")
  expect_error(call(formula = log(y) ~ x), "untransformed")
  expect_error(call(formula = y ~ x + (1 | z)), "Random-effect")
  expect_error(call(control = list(1)), "uniquely named")
  categorical <- d; categorical$y <- factor(rep(c("a", "b"), 12))
  expect_error(call(categorical), "Gaussian outcomes")
  expect_identical(.analysis_mi_binary_outcome(factor(c("b", "a"), levels = c("a", "b", "unused"))), c(1L, 0L))
  expect_error(.analysis_mi_binary_outcome(factor(c("a", "a"))), "exactly two")
  expect_error(.analysis_mi_binary_outcome(c("a", "b")), "coded 0/1")
  expect_error(.analysis_mi_binary_outcome(c(0, 2)), "coded 0/1")
})

test_that("completed datasets preserve observed cells, names and classes", {
  d <- mi_cov_data()
  mi <- multi_impute_analysis(d, y ~ x + z, "x", model = "lm", m = 3, seed = 4)
  expect_length(mi$datasets, 3)
  for (one in mi$datasets) {
    expect_identical(one$y, d$y)
    expect_identical(one$z, d$z)
    expect_identical(one$x[!is.na(d$x)], d$x[!is.na(d$x)])
    expect_true(all(is.finite(one$x)))
  }
  valid <- mi$datasets[[1]]
  check <- function(x) .analysis_mi_validate_datasets(x, d, "x", 1L)
  expect_error(check(list()), "exactly")
  bad <- valid; names(bad)[1] <- "other"
  expect_error(check(list(bad)), "shape and names")
  bad <- valid; bad$z <- as.integer(bad$z)
  expect_error(check(list(bad)), "column classes")
  bad <- valid; bad$x[3] <- Inf
  expect_error(check(list(bad)), "non-finite")
  bad <- valid; bad$x[1] <- bad$x[1] + 1
  expect_error(check(list(bad)), "modified observed")
  bad <- valid; bad$z[1] <- bad$z[1] + 1
  expect_error(check(list(bad)), "non-imputed")
})

test_that("real fits pool aligned coefficients with Rubin variance", {
  d <- mi_cov_data()
  mi <- multi_impute_analysis(d, y ~ x + z, "x", model = "lm", m = 3, seed = 11)
  fits <- with_imputations(mi, function(x) lm(y ~ x + z, x), .progress = FALSE)
  pooled <- pool_mi(fits)
  coefs <- vapply(fits, coef, numeric(3))
  vars <- vapply(fits, function(f) diag(vcov(f)), numeric(3))
  W <- rowMeans(vars); B <- apply(coefs, 1, var)
  expect_identical(pooled$term, rownames(coefs))
  expect_equal(pooled$estimate, unname(rowMeans(coefs)))
  expect_equal(pooled$std.error, unname(sqrt(W + (1 + 1 / 3) * B)))
  expect_equal(pooled$riv, unname((1 + 1 / 3) * B / W))
  expect_match(paste(capture.output(print(pooled)), collapse = " "), "Rubin")
  expect_match(paste(capture.output(print(fits)), collapse = " "), "3/3 successful")
  expect_error(pool_mi(fits, conf.level = 1), "conf.level")
  expect_error(pool_mi(fits, coef_fun = 1), "coef_fun")
  expect_error(pool_mi(fits, vcov_fun = 1), "vcov_fun")
  expect_error(pool_mi(fits, df_fun = 1), "df_fun")
  expect_error(pool_mi(fits, tidy_fun = 1), "tidy_fun")
  expect_error(pool_mi(fits, df_fun = function(f) -1), "finite positive")
  expect_error(pool_mi(fits, coef_fun = function(f) if (identical(f, fits[[1]])) coef(f) else coef(f)[-1]), "Coefficient")
  expect_error(.pool_mi_validate_vcov(matrix(1, 1, 1, dimnames = list("x", NULL)), "x", 1), "both row and column")
  expect_equal(.pool_mi_validate_vcov(diag(2), c("a", "b"), 1), `dimnames<-`(diag(2), list(c("a", "b"), c("a", "b"))))
  expect_error(.pool_mi_validate_vcov(diag(1), c("a", "b"), 1), "exactly one row")
})

test_that("failed fits and diagnostic provenance are not silently admitted", {
  d <- mi_cov_data()
  mi <- multi_impute_analysis(d, y ~ x + z, "x", model = "lm", m = 3, seed = 3)
  expect_error(with_imputations(mi, 1), "must be a function")
  one <- mi; one$datasets <- one$datasets[1]
  expect_error(with_imputations(one, identity), "fewer than 2")
  expect_error(with_imputations(mi, function(d) stop("failure"), .progress = FALSE, .on_error = "stop"), "imputation 1.*failure")
  i <- 0L
  expect_warning(fits <- with_imputations(mi, function(d) {
    i <<- i + 1L
    if (i == 1L) stop("fixture failure")
    lm(y ~ x + z, d)
  }, .progress = FALSE), "1 of 3 fits failed")
  expect_s3_class(fits[[1]], "pigauto_mi_error")
  expect_match(paste(capture.output(print(fits)), collapse = " "), "2/3 successful")
  expect_warning(pooled <- pool_mi(fits), "Dropping 1 fit")
  expect_equal(attr(pooled, "m"), 2L)
  for (workflow in c("pigauto_diagnostic_mi", "pigauto_tree_sensitivity_diagnostic", "pigauto_posterior_plugin_diagnostic")) {
    obj <- structure(list(datasets = mi$datasets, mi_workflow = workflow), class = "pigauto_mi")
    expect_error(with_imputations(obj, identity, .progress = FALSE))
    fit_obj <- structure(unclass(fits), mi_workflow = workflow)
    expect_error(pool_mi(fit_obj))
  }
  expect_error(pool_mi(1), "list of model fits")
  unknown <- fits; attr(unknown, "mi_workflow") <- "unknown"
  expect_error(pool_mi(unknown), "unknown multiple-imputation provenance")
})

test_that("automatic and explicit posterior admission respect trait and covariate scope", {
  f <- mi_cov_tree_data(); d <- f$data; tr <- f$tree
  auto <- function(data = d, species_col = NULL, types = NULL, groups = NULL, covs = NULL)
    .mi_resolve_draws_auto(data, tr, species_col, types, groups, FALSE, covs, verbose = FALSE)
  expect_identical(auto(), "posterior")
  expect_message(expect_identical(auto(covs = data.frame(v = 1:8)), "conformal"), "covariates")
  expect_message(expect_identical(auto(species_col = "species"), "conformal"), "species_col")
  expect_message(expect_identical(auto(groups = list(col = c("x", "y"))), "conformal"), "multi_proportion_groups")
  mixed <- d; mixed$y <- factor(rep(c("a", "b"), 4))
  expect_message(expect_identical(auto(mixed), "conformal"), "not continuous")
  expect_error(multi_impute(d, tr, m = 1), "integer >= 2")
  expect_error(multi_impute(d, tr, m = 2, draws_method = "posterior", covariates = data.frame(v = 1:8)), "does not support `covariates`")
  expect_error(multi_impute(mixed, tr, m = 2, draws_method = "posterior", log_transform = FALSE), "continuous traits only")
  species <- cbind(species = rownames(d), d)
  expect_error(multi_impute(species, tr, m = 2, draws_method = "posterior", species_col = "species", log_transform = FALSE), "does not support `species_col`")
  species <- rbind(species, species[1, ])
  expect_error(multi_impute(species, tr, m = 2, draws_method = "posterior", species_col = "species", log_transform = FALSE), "multiple rows")
})

test_that("posterior controls and diagnostic statistics reject degenerate cases", {
  expect_error(.mip_resolve_control(1, 2), "must be a list")
  expect_error(.mip_resolve_control(list(n_chains = 0), 2), "n_chains")
  expect_error(.mip_resolve_control(list(burnin = -1), 2), "burnin")
  expect_error(.mip_resolve_control(list(n_chains = 1, n_iter = 2, thin = 2, keep_draws = 2), 2), "fewer than m")
  expect_warning(ctl <- .mip_resolve_control(list(n_chains = 1, n_iter = 2, thin = 1, keep_draws = 3), 2), "only 2 kept")
  expect_equal(ctl$n_iter, 2L)
  x <- matrix(1, 4, 2)
  expect_true(is.na(.mip_rhat(x)))
  expect_true(is.na(.mip_ess_bulk(x)))
  expect_true(is.na(.mip_rhat_basic(x)))
  expect_true(is.na(.mip_ess_raw(x)))
  short <- matrix(1:6, 3, 2)
  expect_identical(.mip_split_chains(short), short)
  expect_true(is.na(.mip_ess_raw(short)))
  f <- mi_cov_tree_data(); y <- as.matrix(f$data)
  expect_error(.mip_problem(y[8:1, ], f$tree), "tree tip order")
  expect_error(.mip_problem(unname(y[-1, ]), f$tree), "number of tips")
  y[-1, 1] <- NA_real_
  expect_error(.mip_problem(y, f$tree), "at least 2 observed")
})

test_that("conditional BM means match a dense Gaussian calculation", {
  f <- mi_cov_tree_data(); d <- f$data; d$x[c(2, 5)] <- NA_real_
  data <- suppressMessages(preprocess_traits(d, f$tree, log_transform = FALSE))
  result <- draw_conditional_bm(list(data = data, tree = f$tree), m = 2500, seed = 9)
  X <- data$X_scaled; R <- cov2cor(ape::vcv(f$tree))[rownames(X), rownames(X)]
  C <- kronecker(result$sigma_hat + diag(1e-8, ncol(X)), R)
  v <- as.vector(X); miss <- which(is.na(v)); obs <- which(!is.na(v))
  expected <- C[miss, obs, drop = FALSE] %*% solve(C[obs, obs, drop = FALSE], v[obs])
  expect_equal(as.vector(result$mu_cond)[miss], as.vector(expected), tolerance = 1e-6)
  expected_cov <- C[miss, miss, drop = FALSE] -
    C[miss, obs, drop = FALSE] %*% solve(C[obs, obs, drop = FALSE], C[obs, miss, drop = FALSE])
  simulated <- vapply(result$latent_draws, function(z) as.vector(z)[miss], numeric(length(miss)))
  expect_equal(unname(cov(t(simulated))), unname(expected_cov), tolerance = 0.10)
  expect_identical(result$datasets[[1]]$x[!is.na(d$x)], d$x[!is.na(d$x)])
  expect_identical(result$datasets[[1]]$y, d$y)
  expect_identical(result$imputed_mask, is.na(X))
  expect_identical(result$latent_draws[1:20], draw_conditional_bm(list(data = data, tree = f$tree), m = 20, seed = 9)$latent_draws)
  expect_error(draw_conditional_bm(list(data = data, tree = f$tree), m = 0), "positive integer")
  expect_error(draw_conditional_bm(list(data = data, tree = f$tree), max_cells = 1), "problem too large")
  complete <- suppressMessages(preprocess_traits(f$data, f$tree, log_transform = FALSE))
  out <- draw_conditional_bm(list(data = complete, tree = f$tree), m = 2)
  expect_false(any(out$imputed_mask))
  expect_identical(out$datasets[[1]], complete$X_original)
  expect_equal(out$latent_draws[[1]], complete$X_scaled)
  expect_s3_class(draw_conditional_bm(list(data = data, tree = f$tree), m = 2, sigma_method = "inhouse"), "pigauto_draws_conditional_bm")
})

test_that("tree inputs and pooling reject unsupported cases", {
  f <- mi_cov_tree_data(); tr <- f$tree
  expect_error(resolve_reference_tree(list(tr, tr), reference_tree = 1), "phylo")
  expect_identical(resolve_reference_tree(list(tr, tr), reference_tree = tr), tr)
  expect_error(multi_impute_trees(f$data, list(tr, tr), m_per_tree = 0), "positive integer")
  expect_error(multi_impute_trees(f$data, 1), "list of phylo")
  expect_error(multi_impute_trees(f$data, list(tr)), "at least 2")
  expect_error(multi_impute_trees(f$data, list(tr, 1)), "Non-phylo")
})

test_that("conformal transformed draws honor scale and preserve other cells", {
  pred <- list(imputed = data.frame(p = c(0.25, 0.8), v = c(4, 10)),
               se = cbind(p = c(0, 0), v = c(0, 0)), probabilities = list(),
               conformal_scores = NULL)
  mask <- matrix(c(TRUE, FALSE, TRUE, FALSE), 2, dimnames = list(NULL, c("p", "v")))
  tm <- list(list(name = "p", type = "proportion", mean = 0.3, sd = 2),
             list(name = "v", type = "continuous", mean = 1, sd = 2, log_transform = TRUE))
  out <- .sample_conformal_draw(pred, mask, tm, seed_i = 1)
  expect_equal(out$p, pred$imputed$p)
  expect_equal(out$v, pred$imputed$v)
  pred$se[1, ] <- c(0.2, 0.4)
  out <- .sample_conformal_draw(pred, mask, tm, seed_i = 7)
  expect_true(out$p[1] > 0 && out$p[1] < 1)
  expect_true(out$v[1] > 0)
  expect_identical(out[2, ], pred$imputed[2, ])
})

test_that("low-level guards protect numerical and structural assumptions", {
  expect_error(.pool_mi_validate_vcov(list(1), "a", 1), "base matrix")
  expect_error(.pool_mi_validate_tidy(data.frame(term = "x"), 1), "std.error")
  expect_error(.pool_mi_validate_backend_status(list(), "unknown"), "Unknown optional backend")
  expect_error(.analysis_mi_require("pigauto_package_that_does_not_exist", "test"), "requires optional package")
  d <- data.frame(y = c(0, 1, 2), x = c(1, 2, NA_real_))
  expect_error(.analysis_mi_bayes_norm(d, y ~ x, "x", character(), 2, 1), "positive residual degrees")
  f <- mi_cov_tree_data(); data <- suppressMessages(preprocess_traits(f$data, f$tree, log_transform = FALSE))
  bad <- data; bad$X_scaled <- bad$X_scaled[-1, , drop = FALSE]
  expect_error(draw_conditional_bm(list(data = bad, tree = f$tree)), "tip count")
  expect_error(.mip_resolve_control(NULL, m = 100000), "fewer than m|at least m")
})

test_that("progress and printed tree sensitivity state describe their limits", {
  d <- mi_cov_data()
  mi <- multi_impute_analysis(d, y ~ x + z, "x", model = "lm", m = 2, seed = 4)
  expect_message(fits <- with_imputations(mi, function(d) lm(y ~ x + z, d), .progress = TRUE), "fitting imputation")
  expect_equal(length(fits), 2)
  expect_message(expect_error(with_imputations(mi, function(d) stop("bad fit"), .progress = TRUE, .on_error = "stop"), "bad fit"), "fitting imputation")
  obj <- structure(list(n_trees = 2L, m_per_tree = 2L, m = 4L,
                        datasets = list(d), species_col = NULL,
                        imputed_mask = is.na(d), draws_method = "mc_dropout"),
                   class = "pigauto_mi_trees")
  label <- paste(capture.output(print(obj)), collapse = " ")
  expect_match(label, "baseline perturbations")
  expect_match(label, "unsupported")
})

test_that("tree completions remain descriptive when fitted without a GNN", {
  f <- mi_cov_tree_data(); d <- f$data; d$x[c(2, 5)] <- NA_real_
  warnings <- character()
  withCallingHandlers(suppressMessages(capture.output(
    tree_mi <- multi_impute_trees(d, list(f$tree, f$tree), m_per_tree = 2L,
                                  draws_method = "conformal", share_gnn = FALSE,
                                  gnn = FALSE, log_transform = FALSE, seed = 7,
                                  verbose = TRUE)
  )), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("Few tree-sensitivity", warnings)))
  expect_true(all(grepl("Few tree-sensitivity|Small validation set", warnings)))
  expect_equal(length(tree_mi$datasets), 4L)
  expect_identical(tree_mi$tree_index, c(1L, 1L, 2L, 2L))
  for (one in tree_mi$datasets) {
    expect_identical(one$x[!is.na(d$x)], d$x[!is.na(d$x)])
    expect_identical(one$y, d$y)
  }
  expect_error(with_imputations(tree_mi, function(x) lm(y ~ x, x)), "tree-sensitivity")
})

test_that("joint conditional draws preserve cross-trait covariance", {
  f <- mi_cov_tree_data()
  d <- f$data
  d$y <- 1.5 * d$x + c(0.1, -0.2, 0.3, -0.1, 0.2, -0.3, 0.1, -0.1)
  d$x[2] <- NA_real_
  d$y[2] <- NA_real_
  data <- suppressMessages(preprocess_traits(d, f$tree, log_transform = FALSE))
  n_draws <- 3000L
  draws <- draw_conditional_bm(list(data = data, tree = f$tree), m = n_draws, seed = 413)
  X <- data$X_scaled
  R <- cov2cor(ape::vcv(f$tree))[rownames(X), rownames(X)]
  C <- kronecker(draws$sigma_hat + diag(1e-8, ncol(X)), R)
  missing <- which(is.na(as.vector(X)))
  observed <- which(!is.na(as.vector(X)))
  conditional_cov <- C[missing, missing, drop = FALSE] -
    C[missing, observed, drop = FALSE] %*%
    solve(C[observed, observed, drop = FALSE], C[observed, missing, drop = FALSE])
  sampled <- vapply(draws$latent_draws, function(x) as.vector(x)[missing], numeric(2))
  sampled_cov <- cov(t(sampled))
  # For independent bivariate Normal draws, the sampling variance of the
  # unbiased sample covariance is (V11 * V22 + V12^2) / (n - 1).
  mcse <- sqrt((conditional_cov[1, 1] * conditional_cov[2, 2] +
                  conditional_cov[1, 2]^2) / (n_draws - 1))
  expect_gt(abs(conditional_cov[1, 2]), 0.01)
  expect_lt(abs(sampled_cov[1, 2] - conditional_cov[1, 2]), 4 * mcse)
  for (one in draws$datasets[c(1, n_draws)]) {
    expect_identical(one$x[!is.na(d$x)], d$x[!is.na(d$x)])
    expect_identical(one$y[!is.na(d$y)], d$y[!is.na(d$y)])
  }
})

test_that("single-draw posterior blocks agree with a dense conditional mean", {
  f <- mi_cov_tree_data(); Y <- as.matrix(f$data)
  Y[2, 1] <- NA_real_; Y[5, 2] <- NA_real_
  prob <- .mip_problem(Y, f$tree)
  SP <- matrix(c(1, 0.3, 0.3, 0.7), 2)
  SE <- matrix(c(0.5, 0.1, 0.1, 0.4), 2)
  mu <- c(2, 3)
  R <- cov2cor(ape::vcv(f$tree))[rownames(Y), rownames(Y)]
  C <- kronecker(SP, R) + kronecker(SE, diag(nrow(Y)))
  v <- as.vector(Y); missing <- which(is.na(v)); observed <- which(!is.na(v))
  mu_vec <- rep(mu, each = nrow(Y))
  reference <- mu_vec[missing] + C[missing, observed, drop = FALSE] %*%
    solve(C[observed, observed, drop = FALSE], v[observed] - mu_vec[observed])
  expect_equal(.mip_cond_mean(prob, SP, SE, mu), as.vector(reference), tolerance = 1e-8)
  set.seed(125)
  one <- .mip_fixed_draws(prob, SP, SE, n_draws = 1L, mu = mu, block = 1L)
  set.seed(125)
  again <- .mip_fixed_draws(prob, SP, SE, n_draws = 1L, mu = mu, block = 1L)
  expect_identical(one, again)
  expect_equal(dim(one$ymis), c(2L, 1L))
  expect_true(all(is.finite(one$ymis)))
  expect_equal(as.numeric(one$mu), mu)
  expect_null(.mip_reml_lambda(prob, f$tree, reml_max_n = 7L))
})

test_that("analysis validation detects one group and ill-conditioned inputs", {
  d <- mi_cov_data(); d$group <- factor(rep("only", nrow(d)))
  expect_error(multi_impute_analysis(d, y ~ x + z + (1 | group), "x",
                                     model = "lmer", m = 2), "at least two groups")
  expect_length(.analysis_mi_findbars_safe(y ~ x + (1 | group)), 1)
  # Rescaling a non-collinear auxiliary keeps QR rank while making the
  # design condition number exceed the numerical stability cutoff.
  d$tiny <- d$z * 1e-10
  expect_error(.analysis_mi_bayes_norm(d, y ~ x, "x", "tiny", 2, 1), "ill-conditioned")
})

test_that("pool printing highlights large missing-information fractions", {
  base <- data.frame(x = 1:8, noise = c(.02, -.01, .03, -.02, .01, -.03, .02, -.01))
  fits <- lapply(c(-2, 0, 2), function(shift) {
    lm(y ~ x, transform(base, y = 1 + 0.3 * x + shift + noise))
  })
  expect_warning(pooled <- pool_mi(fits), "unverified provenance")
  expect_gt(pooled$fmi[pooled$term == "(Intercept)"], 0.5)
  expect_match(paste(capture.output(print(pooled)), collapse = " "), "Consider M >=")
})

test_that("diagnostic sampling skips absent traits and decodes zero-inflated counts", {
  pred <- list(imputed = data.frame(count = c(99L, 7L), a = c(0.4, 0.2), b = c(0.6, 0.8)),
               se = matrix(0, 2, 3, dimnames = list(NULL, c("count", "a", "b"))),
               probabilities = list(count = c(1, 1)), conformal_scores = NULL,
               imputed_latent = cbind(gate = c(0, 0), magnitude = c(log1p(3), 0)))
  mask <- matrix(FALSE, 2, 3, dimnames = list(NULL, c("count", "a", "b")))
  mask[1, "count"] <- TRUE
  traits <- list(
    list(name = "not_present", type = "continuous"),
    list(name = "count", type = "zi_count", latent_cols = 1:2, mean = 0, sd = 1),
    list(name = "bad_group", type = "multi_proportion", input_cols = c("a", "missing")),
    list(name = "unmasked_group", type = "multi_proportion", input_cols = c("a", "b"))
  )
  out <- .sample_conformal_draw(pred, mask, traits, seed_i = 30)
  expect_identical(out$count, c(3L, 7L))
  expect_identical(out$a, pred$imputed$a)
  expect_identical(out$b, pred$imputed$b)
  # Present in predictions but absent from the input mask: no draw is made.
  mask2 <- mask[, "count", drop = FALSE]
  traits2 <- list(list(name = "a", type = "continuous"),
                  list(name = "partial", type = "multi_proportion", input_cols = c("a", "b")))
  expect_identical(.sample_conformal_draw(pred, mask2, traits2, seed_i = 30), pred$imputed)
})

test_that("conditional draws restore missing latent row names from data", {
  f <- mi_cov_tree_data(); d <- f$data; d$x[c(1, 2)] <- NA_real_
  data <- suppressMessages(preprocess_traits(d, f$tree, log_transform = FALSE))
  reference <- draw_conditional_bm(list(data = data, tree = f$tree), m = 2, seed = 21)
  rownames(data$X_scaled) <- NULL
  restored <- draw_conditional_bm(list(data = data, tree = f$tree), m = 2, seed = 21)
  expect_equal(restored$mu_cond, reference$mu_cond)
  expect_identical(restored$datasets, reference$datasets)
  expect_identical(rownames(restored$mu_cond), data$species_names)
})

test_that("posterior output assembles a single missing cell", {
  f <- mi_cov_tree_data(); d <- f$data; d$x[2] <- NA_real_
  # Deliberately short: this tests shape and preservation, not convergence.
  ctl <- list(n_chains = 1L, burnin = 2L, n_iter = 8L, thin = 1L,
              keep_draws = 8L, auto_extend = FALSE)
  warnings <- character()
  mi <- withCallingHandlers(
    multi_impute(d, f$tree, m = 2L, draws_method = "posterior",
                 log_transform = FALSE, posterior_control = ctl, seed = 3,
                 verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("did not meet the convergence rule", warnings)))
  expect_false(mi$posterior$converged)
  expect_equal(nrow(mi$posterior$cell_interval), 1L)
  expect_equal(mi$posterior$cell_interval$row, 2L)
  expect_identical(mi$posterior$cell_interval$trait, "x")
  expect_true(all(is.finite(unlist(mi$posterior$cell_interval[c("lower", "median", "upper")]))))
  expect_lte(mi$posterior$cell_interval$lower, mi$posterior$cell_interval$median)
  expect_lte(mi$posterior$cell_interval$median, mi$posterior$cell_interval$upper)
  expect_true(all(vapply(mi$datasets, function(z) identical(z$x[-2], d$x[-2]) &&
                           identical(z$y, d$y), logical(1))))
})
test_that("analysis imputation rejects missing auxiliaries and other incomplete inputs", {
  d <- mi_cov_data()
  expect_error(multi_impute_analysis(d, y ~ x + z, "x", auxiliary = "absent",
                                     model = "lm", m = 2), "Unknown.*absent")
  d$z[1] <- NA_real_
  expect_error(multi_impute_analysis(d, y ~ x + z, "x", model = "lm", m = 2),
               "incomplete columns: z")
  expect_error(multi_impute_analysis(d, absent ~ x, "x", model = "lm", m = 2),
               "incomplete columns: z")
  d$z[1] <- 0
  expect_error(multi_impute_analysis(d, absent ~ x, "x", model = "lm", m = 2),
               "Formula variables not found.*absent")
})

test_that("optional SMCFCS draws preserve observations and support real pooling", {
  if (requireNamespace("smcfcs", quietly = TRUE)) {
    set.seed(41)
    n <- 100L
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    d$y <- rbinom(n, 1, plogis(.2 + .4 * d$x - .2 * d$z))
    d$x[seq(3, n, by = 9)] <- NA_real_
    mi <- multi_impute_analysis(d, y ~ x + z, "x", model = "glm", m = 2L,
                                seed = 12L, control = list(numit = 3L))
    expect_identical(mi$engine, "smcfcs")
    expect_length(mi$datasets, 2L)
    for (one in mi$datasets) {
      expect_identical(one$x[!is.na(d$x)], d$x[!is.na(d$x)])
      expect_identical(one$y, d$y)
      expect_identical(one$z, d$z)
      expect_true(all(is.finite(one$x)))
    }
    fits <- with_imputations(mi, function(z) glm(y ~ x + z, z, family = binomial()),
                             .progress = FALSE)
    pooled <- pool_mi(fits)
    coefs <- vapply(fits, coef, numeric(3))
    vars <- vapply(fits, function(f) diag(as.matrix(vcov(f))), numeric(3))
    expect_equal(pooled$estimate, unname(rowMeans(coefs)))
    expect_equal(pooled$std.error,
                 unname(sqrt(rowMeans(vars) + 1.5 * apply(coefs, 1, var))))
  }
})

test_that("optional jomo draws preserve the rows of a real mixed model", {
  if (requireNamespace("jomo", quietly = TRUE) && requireNamespace("lme4", quietly = TRUE)) {
    set.seed(52)
    n <- 120L
    d <- data.frame(x = rnorm(n), z = rnorm(n), group = factor(rep(1:12, each = 10)))
    intercept <- rnorm(12, sd = .5)
    d$y <- 1 + .4 * d$x + .2 * d$z + intercept[d$group] + rnorm(n)
    d$x[seq(3, n, by = 9)] <- NA_real_
    mi <- multi_impute_analysis(d, y ~ x + z + (1 | group), "x", model = "lmer",
                                m = 2L, seed = 12L,
                                control = list(nburn = 5L, nbetween = 2L))
    expect_identical(mi$engine, "jomo")
    expect_length(mi$datasets, 2L)
    for (one in mi$datasets) {
      expect_identical(one$x[!is.na(d$x)], d$x[!is.na(d$x)])
      expect_identical(one[c("y", "z", "group")], d[c("y", "z", "group")])
      expect_true(all(is.finite(one$x)))
    }
    fits <- with_imputations(mi, function(z) lme4::lmer(y ~ x + z + (1 | group), z),
                             .progress = FALSE)
    pooled <- pool_mi(fits)
    coefs <- vapply(fits, lme4::fixef, numeric(3))
    vars <- vapply(fits, function(f) diag(as.matrix(vcov(f))), numeric(3))
    expect_equal(pooled$estimate, unname(rowMeans(coefs)))
    expect_equal(pooled$std.error,
                 unname(sqrt(rowMeans(vars) + 1.5 * apply(coefs, 1, var))))
  }
})

test_that("pooling keeps aligned metadata after a failed bare-list fit", {
  # Artificial attributes exercise alignment only; these are not tree draws.
  d <- data.frame(x = 1:8, y = c(1.2, 2.1, 2.7, 4.1, 4.8, 6.2, 6.8, 8.1))
  fits <- list(lm(y ~ x, d), structure("failed", class = "try-error"),
               lm(y ~ x, transform(d, y = y + .5)))
  attr(fits, "tree_index") <- c(1L, 1L, 2L)
  attr(fits, "n_trees") <- 2L
  attr(fits, "m_per_tree") <- 2L
  warnings <- character()
  pooled <- withCallingHandlers(pool_mi(fits), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("Dropping 1 fit", warnings)))
  expect_true(any(grepl("unverified provenance", warnings)))
  expect_length(warnings, 2L)
  expect_identical(attr(pooled, "tree_index"), c(1L, 2L))
  expect_identical(attr(pooled, "n_trees"), 2L)
  expect_identical(attr(pooled, "m_per_tree"), 2L)
  expect_equal(attr(pooled, "m"), 2L)
})

test_that("analysis imputation reports a missing optional engine", {
  real_require <- base::requireNamespace
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (package == "smcfcs") return(FALSE)
      real_require(package, ...)
    }, .package = "base")
  d <- mi_cov_data(); d$y <- as.integer(seq_len(nrow(d)) %% 2L)
  expect_error(multi_impute_analysis(d, y ~ x + z, "x", model = "glm", m = 2),
               "requires optional package 'smcfcs'")
})

test_that("SMCFCS adapter refuses malformed external completions", {
  if (requireNamespace("smcfcs", quietly = TRUE)) {
    d <- mi_cov_data(); d$y <- as.integer(seq_len(nrow(d)) %% 2L)
    control <- list(numit = 2L, rjlimit = 100L)
    testthat::local_mocked_bindings(smcfcs = function(...) list(impDatasets = list(d)),
                                    .package = "smcfcs")
    expect_error(.analysis_mi_smcfcs(d, y ~ x + z, "x", character(), 2L, 1L, control),
                 "requested number")
  }
})

test_that("SMCFCS adapter refuses completion with wrong row count", {
  if (requireNamespace("smcfcs", quietly = TRUE)) {
    d <- mi_cov_data(); d$y <- as.integer(seq_len(nrow(d)) %% 2L)
    control <- list(numit = 2L, rjlimit = 100L)
    testthat::local_mocked_bindings(smcfcs = function(...) list(impDatasets = list(d[-1, ], d)),
                                    .package = "smcfcs")
    expect_error(.analysis_mi_smcfcs(d, y ~ x + z, "x", character(), 2L, 1L, control),
                 "incompatible rows")
  }
})

test_that("jomo adapter refuses malformed long-format output", {
  if (requireNamespace("jomo", quietly = TRUE)) {
    d <- mi_cov_data(); d$group <- rep(1:4, each = 6)
    control <- list(nburn = 2L, nbetween = 1L)
    testthat::local_mocked_bindings(jomo.smc = function(...) data.frame(x = 1),
                                    .package = "jomo")
    expect_error(.analysis_mi_jomo(d, y ~ x + z + (1 | group), "x", character(),
                                   "group", 2L, 1L, control), "incompatible long-format")
  }
})

test_that("jomo adapter rejects duplicate external row identifiers", {
  if (requireNamespace("jomo", quietly = TRUE)) {
    d <- mi_cov_data(); d$group <- rep(1:4, each = 6)
    control <- list(nburn = 2L, nbetween = 1L)
    bad <- data.frame(Imputation = rep(1L, nrow(d)), id = rep(1L, nrow(d)), x = 1)
    testthat::local_mocked_bindings(jomo.smc = function(...) bad, .package = "jomo")
    expect_error(.analysis_mi_jomo(d, y ~ x + z + (1 | group), "x", character(),
                                   "group", 2L, 1L, control), "incompatible row identifiers")
  }
})

test_that("conditional BM draws report a missing Matrix dependency", {
  real_require <- base::requireNamespace
  testthat::local_mocked_bindings(requireNamespace = function(package, ...) {
    if (package == "Matrix") return(FALSE)
    real_require(package, ...)
  }, .package = "base")
  expect_error(draw_conditional_bm(list()), "requires the 'Matrix' package")
})
