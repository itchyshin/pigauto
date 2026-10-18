# Cross-cutting public contract checks for the 2026-10 validation campaign.
# Keep this fixture small and use the default GNN-off baseline so the test
# remains useful on machines without torch/libtorch.

validation_contract_fixture <- function(n = 48L, seed = 20261018L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  species <- tree$tip.label
  comp <- matrix(stats::rgamma(n * 3L, shape = 2), nrow = n,
                 dimnames = list(species, c("comp_a", "comp_b", "comp_c")))
  comp <- comp / rowSums(comp)
  dat <- data.frame(
    continuous = as.numeric(ape::rTraitCont(tree)) + 3,
    count = as.integer(stats::rpois(n, lambda = 4)),
    binary = factor(rep(c("absent", "present"), length.out = n),
                    levels = c("absent", "present")),
    categorical = factor(rep(c("low", "mid", "high"), length.out = n),
                         levels = c("low", "mid", "high")),
    ordinal = ordered(rep(c("small", "medium", "large"), length.out = n),
                      levels = c("small", "medium", "large")),
    proportion = stats::runif(n, min = 0.08, max = 0.92),
    zi_count = as.integer(ifelse(seq_len(n) %% 3L == 0L, 0L,
                                 stats::rpois(n, lambda = 3) + 1L)),
    comp_a = comp[, "comp_a"], comp_b = comp[, "comp_b"],
    comp_c = comp[, "comp_c"],
    row.names = species,
    check.names = FALSE
  )
  # Distinct missing cells exercise each decoder, plus whole-row missingness
  # for the explicit composition group.
  for (j in seq_len(7L)) dat[(j + 1L) + 4L * (0:3), j] <- NA
  dat[c(5L, 17L, 29L, 41L), c("comp_a", "comp_b", "comp_c")] <- NA
  list(tree = tree, data = dat,
       multi_proportion_groups = list(composition = c("comp_a", "comp_b", "comp_c")))
}

validation_contract_impute <- function(dat, tree, groups, seed = 20261018L) {
  impute(
    dat, tree,
    trait_types = c(proportion = "proportion", zi_count = "zi_count"),
    multi_proportion_groups = groups,
    missing_frac = 0, lambda_mode = "fixed_1",
    verbose = FALSE, seed = seed
  )
}

test_that("each supported trait route imputes on its own through the public API", {
  fx <- validation_contract_fixture()
  cases <- list(
    continuous = list(columns = "continuous", types = NULL, groups = NULL,
                      type = "continuous"),
    count = list(columns = "count", types = NULL, groups = NULL, type = "count"),
    binary = list(columns = "binary", types = NULL, groups = NULL, type = "binary"),
    categorical = list(columns = "categorical", types = NULL, groups = NULL,
                       type = "categorical"),
    ordinal = list(columns = "ordinal", types = NULL, groups = NULL,
                   type = "ordinal"),
    proportion = list(columns = "proportion",
                      types = c(proportion = "proportion"), groups = NULL,
                      type = "proportion"),
    zi_count = list(columns = "zi_count", types = c(zi_count = "zi_count"),
                    groups = NULL, type = "zi_count"),
    multi_proportion = list(columns = c("comp_a", "comp_b", "comp_c"),
                            types = NULL,
                            groups = fx$multi_proportion_groups,
                            type = "multi_proportion")
  )

  for (nm in names(cases)) {
    spec <- cases[[nm]]
    dat <- fx$data[, spec$columns, drop = FALSE]
    res <- impute(
      dat, fx$tree,
      trait_types = spec$types,
      multi_proportion_groups = spec$groups,
      missing_frac = 0, lambda_mode = "fixed_1",
      verbose = FALSE, seed = 20261018L
    )
    map_name <- if (identical(spec$type, "multi_proportion")) "composition" else nm
    expect_identical(res$data$trait_map[[map_name]]$type, spec$type,
                     info = paste("dispatch type:", nm))
    expect_true(all(stats::complete.cases(res$completed)), info = nm)
    for (column in spec$columns) {
      observed <- !is.na(dat[[column]])
      expect_equal(res$completed[[column]][observed], dat[[column]][observed],
                   info = paste("observed values preserved:", nm, column))
    }
    if (identical(spec$type, "multi_proportion")) {
      expect_equal(unname(rowSums(res$completed[spec$columns])), rep(1, nrow(dat)),
                   tolerance = 1e-8)
    }
  }
})

test_that("mixed public imputation preserves all eight trait contracts", {
  fx <- validation_contract_fixture()
  original <- fx$data
  result <- validation_contract_impute(original, fx$tree,
                                       fx$multi_proportion_groups)
  completed <- result$completed

  expect_s3_class(result, "pigauto_result")
  expect_named(result$data$trait_map,
                c("continuous", "count", "binary", "categorical", "ordinal",
                  "proportion", "zi_count", "composition"), ignore.order = TRUE)
  expect_identical(vapply(result$data$trait_map, `[[`, character(1), "type")
                   [c("continuous", "count", "binary", "categorical", "ordinal",
                      "proportion", "zi_count", "composition")],
                   c(continuous = "continuous", count = "count", binary = "binary",
                     categorical = "categorical", ordinal = "ordinal",
                     proportion = "proportion", zi_count = "zi_count",
                     composition = "multi_proportion"))

  expect_identical(class(completed$count), class(original$count))
  expect_identical(class(completed$binary), class(original$binary))
  expect_identical(class(completed$categorical), class(original$categorical))
  expect_identical(class(completed$ordinal), class(original$ordinal))
  expect_true(all(stats::complete.cases(completed)))
  for (nm in names(original)) {
    observed <- !is.na(original[[nm]])
    expect_equal(completed[[nm]][observed], original[[nm]][observed],
                 info = paste("observed values preserved:", nm))
  }

  expect_true(all(completed$proportion >= 0 & completed$proportion <= 1))
  expect_true(all(completed$zi_count >= 0 & completed$zi_count == floor(completed$zi_count)))
  expect_true(all(as.matrix(completed[c("comp_a", "comp_b", "comp_c")]) >= 0))
  expect_equal(unname(rowSums(completed[c("comp_a", "comp_b", "comp_c")])),
               rep(1, nrow(completed)), tolerance = 1e-8)

  pred <- result$prediction
  expect_true(all(is.finite(pred$conformal_lower[, "continuous"])))
  expect_true(all(is.finite(pred$conformal_upper[, "continuous"])))
  expect_true(all(pred$conformal_lower[, "continuous"] <=
                  pred$conformal_upper[, "continuous"]))
  expect_true(all(pred$conformal_lower[, "count"] >= 0))
  expect_true(all(pred$conformal_lower[, "proportion"] >= 0 &
                  pred$conformal_upper[, "proportion"] <= 1))
})

test_that("mixed public imputation honors default lambda estimation and repeats by seed", {
  fx <- validation_contract_fixture(n = 80L)
  original <- fx$data
  run_default <- function() {
    impute(
      original, fx$tree,
      trait_types = c(proportion = "proportion", zi_count = "zi_count"),
      multi_proportion_groups = fx$multi_proportion_groups,
      missing_frac = 0.2, verbose = FALSE, seed = 20261018L
    )
  }
  first <- run_default()
  second <- run_default()

  expect_identical(first$fit$model_config$lambda_mode, "estimate")
  expect_false(isTRUE(first$fit$model_config$gnn))
  expect_gt(length(first$fit$splits$val_idx), 0L)
  expect_gt(length(first$fit$splits$test_idx), 0L)
  expect_true(all(stats::complete.cases(first$completed)))
  for (nm in names(original)) {
    observed <- !is.na(original[[nm]])
    expect_equal(first$completed[[nm]][observed], original[[nm]][observed],
                 info = paste("observed values preserved:", nm))
  }
  numeric_traits <- c("continuous", "count", "proportion", "zi_count",
                      "comp_a", "comp_b", "comp_c")
  for (nm in numeric_traits) {
    missing <- is.na(original[[nm]])
    expect_true(all(is.finite(first$completed[[nm]][missing])),
                info = paste("finite missing predictions:", nm))
  }
  expect_identical(class(first$completed$count), class(original$count))
  expect_identical(class(first$completed$ordinal), class(original$ordinal))
  expect_equal(unname(rowSums(first$completed[c("comp_a", "comp_b", "comp_c")])),
               rep(1, nrow(original)), tolerance = 1e-8)
  expect_equal(second$completed, first$completed, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(second$prediction$imputed, first$prediction$imputed,
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("row order and RDS round-trip preserve named baseline predictions", {
  fx <- validation_contract_fixture()
  first <- validation_contract_impute(fx$data, fx$tree,
                                      fx$multi_proportion_groups)
  permutation <- rev(seq_len(nrow(fx$data)))
  reordered <- validation_contract_impute(
    fx$data[permutation, , drop = FALSE], fx$tree,
    fx$multi_proportion_groups
  )
  expect_identical(rownames(reordered$completed), rownames(fx$data)[permutation])
  expect_equal(reordered$completed[order(rownames(reordered$completed)), ],
               first$completed[order(rownames(first$completed)), ],
               tolerance = 1e-10, ignore_attr = TRUE)

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(first$fit, path)
  restored <- readRDS(path)
  before <- predict(first$fit, return_se = FALSE)
  after <- predict(restored, return_se = FALSE)
  expect_equal(after$imputed, before$imputed, tolerance = 1e-12,
               ignore_attr = TRUE)
})
