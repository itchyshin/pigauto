test_that("species folds keep every observation and trait together", {
  tree <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  d <- data.frame(species = rep(tree$tip.label, each = 2),
                  x = seq_len(16) + 0.5, y = seq_len(16) + 20.5)
  pd <- preprocess_traits(d, tree, species_col = "species", log_transform = FALSE)
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 13L, "species")
  expect_length(folds, 4L)
  seen <- integer(nrow(pd$X_scaled))
  for (i in seq_along(folds)) {
    idx <- folds[[i]]
    rows <- unique(((idx - 1L) %% nrow(pd$X_scaled)) + 1L)
    expect_setequal(idx, as.vector(outer(rows, seq_len(ncol(pd$X_scaled)) - 1L,
                                         function(r, j) r + j * nrow(pd$X_scaled))))
    for (sp in unique(pd$obs_to_species[rows]))
      expect_true(all(which(pd$obs_to_species == sp) %in% rows))
    seen[rows] <- seen[rows] + 1L
  }
  expect_true(all(seen == 1L))
  expect_identical(folds, pigauto:::make_grouped_cv_folds(pd, tree, 4L, 13L, "species"))
})

test_that("clade folds partition the tree into rooted clades", {
  tree <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  d <- data.frame(x = seq_len(8) + 0.5, y = seq_len(8) + 10.5, row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 7L, "clade")
  species <- lapply(folds, function(idx) tree$tip.label[unique(((idx - 1L) %% 8L) + 1L)])
  expect_setequal(unlist(species), tree$tip.label)
  expect_true(all(vapply(species, length, integer(1)) == 2L))
  expect_true(all(vapply(species, function(s) {
    node <- ape::getMRCA(tree, s)
    setequal(ape::extract.clade(tree, node)$tip.label, s)
  }, logical(1))))
})

test_that("grouped cross-validation exposes holdout type and runs", {
  set.seed(21)
  tree <- ape::stree(128, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(x = rnorm(128), y = rnorm(128), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  for (mode in c("species", "clade")) {
    cv <- cross_validate(pd, tree, k = 4L, seeds = 2L, holdout = mode,
                         gnn = FALSE, verbose = FALSE)
    expect_s3_class(cv, "pigauto_cv")
    expect_identical(cv$holdout, mode)
    expect_gt(nrow(cv$results), 0L)
  }
})

test_that("grouped fold scaling uses only training values", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1,e:1,f:1);")
  d <- data.frame(x = c(1, 2, 3, 4, 5, 100), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  splits <- list(val_idx = 5L, test_idx = 6L, n = 6L, p = 1L)
  fold <- pigauto:::prepare_grouped_cv_data(pd, splits)
  expect_equal(fold$trait_map$x$mean, mean(1:4))
  expect_equal(fold$trait_map$x$sd, sd(1:4))
  expect_equal(fold$X_scaled[6, 1], (100 - mean(1:4)) / sd(1:4))
})

test_that("grouped folds reject a trait without enough training observations", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1,e:1,f:1);")
  d <- data.frame(x = 1:6 + 0.5, row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  splits <- list(val_idx = 2:3, test_idx = 4:6, n = 6L, p = 1L)
  expect_error(pigauto:::prepare_grouped_cv_data(pd, splits), "training observations")
})

test_that("grouped folds preserve explicit no-centering or scaling", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1,e:1,f:1);")
  d <- data.frame(x = 1:6 + 0.5, row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE, center = FALSE, scale = FALSE)
  fold <- pigauto:::prepare_grouped_cv_data(pd, list(val_idx = 5L, test_idx = 6L))
  expect_equal(fold$trait_map$x$mean, 0)
  expect_equal(fold$trait_map$x$sd, 1)
})

test_that("grouped folds expose the held-out species", {
  tree <- ape::stree(32, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  set.seed(4)
  d <- data.frame(x = rnorm(32), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  cv <- suppressWarnings(cross_validate(pd, tree, k = 4L, seeds = 1L,
                                        holdout = "clade", gnn = FALSE,
                                        verbose = FALSE))
  expect_named(cv$fold_assignments, c("rep", "fold", "species"))
  expect_setequal(cv$fold_assignments$species, tree$tip.label)
  expect_equal(nrow(cv$fold_assignments), length(tree$tip.label))
})

test_that("grouped test values cannot influence held-out predictions", {
  tree <- ape::stree(32, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  set.seed(25)
  d <- data.frame(x = rnorm(32), y = rnorm(32), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 7L, "species")
  splits <- list(test_idx = folds[[1]], val_idx = folds[[2]], n = 32L, p = 2L)
  altered <- pd
  test_rows <- unique(((folds[[1]] - 1L) %% 32L) + 1L)
  altered$X_original[test_rows, ] <- altered$X_original[test_rows, ] + 1000
  run <- function(x) {
    x <- pigauto:::prepare_grouped_cv_data(x, splits)
    fit <- suppressWarnings(fit_pigauto(x, tree, splits = splits,
                                         gnn = FALSE, seed = 3L, verbose = FALSE))
    pred <- predict(fit, .mask_observed_idx = c(splits$val_idx, splits$test_idx))
    pred$imputed_latent[splits$test_idx]
  }
  expect_equal(run(pd), run(altered), tolerance = 1e-10)
})

test_that("grouped preprocessing preserves mixed and compositional schemas", {
  tree <- ape::stree(16, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(x = seq_len(16) + 0.5,
                  b = factor(rep(c("a", "b"), 8)),
                  c1 = seq(0.1, 0.8, length.out = 16), row.names = tree$tip.label)
  d$c2 <- 1 - d$c1
  pd <- preprocess_traits(d, tree, log_transform = FALSE,
                          multi_proportion_groups = list(comp = c("c1", "c2")))
  idx <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 2L, "species")
  fold <- pigauto:::prepare_grouped_cv_data(pd, list(val_idx = idx[[1]], test_idx = idx[[2]]))
  expect_identical(colnames(fold$X_scaled), colnames(pd$X_scaled))
  expect_identical(fold$trait_map$b$levels, levels(d$b))
  expect_identical(fold$trait_map$comp$input_cols, c("c1", "c2"))
})

test_that("grouped CV works with repeated observations and mixed types", {
  tree <- ape::stree(32, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  set.seed(44)
  d <- data.frame(species = rep(tree$tip.label, each = 2),
                  x = rnorm(64), b = factor(rep(c("a", "b"), 32)))
  pd <- preprocess_traits(d, tree, species_col = "species", log_transform = FALSE)
  cv <- suppressWarnings(cross_validate(pd, tree, k = 4L, seeds = 4L,
                                       holdout = "species", gnn = FALSE,
                                       verbose = FALSE))
  expect_s3_class(cv, "pigauto_cv")
  expect_setequal(cv$fold_assignments$species, tree$tip.label)
  expect_true(all(c("x", "b") %in% cv$results$trait))
})

test_that("grouped CV supports the optional GNN", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("libtorch not installed")
  tree <- ape::stree(32, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  set.seed(45)
  d <- data.frame(x = rnorm(32), y = rnorm(32), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  cv <- suppressWarnings(cross_validate(pd, tree, k = 4L, seeds = 4L,
                                       holdout = "clade", gnn = TRUE,
                                       epochs = 1L, hidden_dim = 8L,
                                       k_eigen = 4L, verbose = FALSE))
  expect_gt(nrow(cv$results), 0L)
})

test_that("grouped CV can recover a known Brownian trait on a balanced tree", {
  tree <- ape::stree(64, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  set.seed(502)
  trait <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  d <- data.frame(x = unname(trait[tree$tip.label]), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  for (mode in c("species", "clade")) {
    cv <- suppressWarnings(cross_validate(pd, tree, k = 4L, seeds = 8L,
                                         holdout = mode, gnn = FALSE,
                                         lambda_mode = "fixed_1", verbose = FALSE))
    rmse <- cv$results$value[cv$results$trait == "x" & cv$results$metric == "rmse"]
    expect_length(rmse, 4L)
    expect_true(all(is.finite(rmse)))
    expect_lt(mean(rmse), 1.5)
  }
})

test_that("repeated rows from one species do not count as multiple training species", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  d <- data.frame(species = rep(tree$tip.label, each = 3), x = seq_len(12) + 0.5)
  pd <- preprocess_traits(d, tree, species_col = "species", log_transform = FALSE)
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 1L, "species")
  splits <- list(val_idx = folds[[1]], test_idx = c(folds[[2]], folds[[3]]))
  expect_error(pigauto:::prepare_grouped_cv_data(pd, splits), "training species")
})

test_that("cross-validation keeps compositional and zero-count metrics", {
  metrics <- data.frame(trait = "comp", type = "multi_proportion",
                        rmse = NA_real_, pearson_r = NA_real_, coverage_95 = NA_real_,
                        mae = NA_real_, spearman_rho = NA_real_, accuracy = NA_real_,
                        brier = NA_real_, zero_accuracy = 0.75, aitchison = 0.4,
                        rmse_clr = 0.3, simplex_mae = 0.2)
  result <- pigauto:::pivot_eval_to_long(metrics, 1L, 1L)
  expect_setequal(result$metric, c("zero_accuracy", "aitchison", "rmse_clr", "simplex_mae"))
})

test_that("grouped folds do not silently score empty holdouts", {
  tree <- ape::stree(8, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(x = c(1, 2, 3, rep(NA_real_, 5)), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  expect_error(pigauto:::make_grouped_cv_folds(pd, tree, 4L, 1L, "species"),
               "no observed trait values")
})

test_that("clade holdouts reject unrooted trees", {
  set.seed(32)
  tree <- ape::unroot(ape::rtree(12))
  d <- data.frame(x = seq_len(12) + 0.5, row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  expect_false(ape::is.rooted(tree))
  expect_error(pigauto:::make_grouped_cv_folds(pd, tree, 3L, 1L, "clade"),
               "rooted tree")
})

test_that("clade folds support deep rooted trees without recursion limits", {
  tree <- ape::stree(1500, type = "left")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(x = seq_len(1500) + 0.5, row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 1L, "clade")
  expect_length(folds, 4L)
  expect_setequal(unlist(folds), seq_len(1500))
})

test_that("grouped fold scoring indices exclude naturally missing values", {
  tree <- ape::stree(16, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(x = seq_len(16) + 0.5, y = seq_len(16) + 1.5,
                  row.names = tree$tip.label)
  d$x[1:5] <- NA
  d$y[6:10] <- NA
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 2L, "species")
  expect_true(all(vapply(folds, function(i) all(!is.na(pd$X_scaled[i])), logical(1))))
  expect_setequal(unlist(folds), which(!is.na(pd$X_scaled)))
})

test_that("fold membership includes species with no observed traits", {
  tree <- ape::stree(16, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(x = seq_len(16) + 0.5, row.names = tree$tip.label)
  d$x[1] <- NA_real_
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  cv <- suppressWarnings(cross_validate(pd, tree, k = 4L, seeds = 2L,
                                       holdout = "species", gnn = FALSE,
                                       verbose = FALSE))
  expect_setequal(cv$fold_assignments$species, tree$tip.label)
  expect_equal(nrow(cv$fold_assignments), 16L)
})

test_that("sparse mixed folds exclude missing latent entries but preserve observed units", {
  tree <- ape::stree(16, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  d <- data.frame(z = rep(c(0L, 1L, 3L, NA_integer_), 4),
                  cat = factor(rep(c("a", "b", "c", NA_character_), 4)),
                  c1 = rep(c(0.2, 0.4, 0.6, NA_real_), 4),
                  c2 = rep(c(0.8, 0.6, 0.4, NA_real_), 4),
                  row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, trait_types = c(z = "zi_count"),
                          multi_proportion_groups = list(comp = c("c1", "c2")))
  folds <- pigauto:::make_grouped_cv_folds(pd, tree, 4L, 2L, "species")
  expect_setequal(unlist(folds), which(!is.na(pd$X_scaled)))
  expect_true(all(vapply(folds, function(i) all(!is.na(pd$X_scaled[i])), logical(1))))
  expect_setequal(unlist(attr(folds, "species_folds")), seq_len(16))
})

test_that("greedy clade partition reports its limitation honestly", {
  tree <- ape::read.tree(text = "((a:1,(b:1,c:1,d:1):1):1,(e:1,f:1,g:1):1);")
  d <- data.frame(x = seq_along(tree$tip.label) + 0.5, row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  expect_error(pigauto:::make_grouped_cv_folds(pd, tree, 4L, 1L, "clade"),
               "greedy clade partition did not find")
})
