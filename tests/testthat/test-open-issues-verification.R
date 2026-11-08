test_that("issue 259 public impute path is invariant to factor coding", {
  # Reporter DGP and seed, checking coding invariance rather than claiming
  # lambda recovery from one simulated dataset.
  set.seed(1)
  n <- 200
  tree <- ape::rtree(n)
  tree$tip.label <- paste0("sp", seq_len(n))
  R <- ape::vcv(tree)
  R <- R / max(diag(R))
  g <- factor(sample(c("a", "b"), n, replace = TRUE))
  y <- as.numeric(t(chol(0.5 * R + 0.5 * diag(n))) %*% rnorm(n)) +
    1.5 * (g == "b")
  y[sample(n, 40)] <- NA
  traits <- data.frame(y = y, row.names = tree$tip.label)
  run <- function(covs) suppressWarnings(impute(
    traits, tree, covariates = covs, log_transform = FALSE,
    seed = 1, verbose = FALSE, gnn = FALSE))
  a <- run(data.frame(g = g, row.names = tree$tip.label))
  b <- run(data.frame(g_b = as.numeric(g == "b"), row.names = tree$tip.label))
  expect_true(all(is.finite(a$baseline$lambda_per_trait)))
  expect_equal(a$baseline$lambda_per_trait, b$baseline$lambda_per_trait,
               tolerance = 1e-6)
  expect_equal(a$completed$y, b$completed$y, tolerance = 1e-6)
  expect_equal(a$prediction$se, b$prediction$se, tolerance = 1e-6)
})

test_that("issue 258 data-only rows survive memory and CSV routes unchanged", {
  set.seed(258)
  tree <- ape::rtree(18)
  traits <- data.frame(y = rnorm(20), row.names = c(tree$tip.label, "extra1", "extra2"))
  traits$y[c(1, 19)] <- NA_real_
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv))
  write.csv(data.frame(species = rownames(traits), traits, row.names = NULL),
            csv, row.names = FALSE)
  for (d in list(traits, read_traits(csv, species_col = "species"))) {
    warnings <- character()
    result <- withCallingHandlers(
      impute(d, tree, gnn = FALSE, verbose = FALSE, seed = 1),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_true(any(grepl("not.*modeled|not.*modelled|excluded from", warnings)))
    expect_false(any(grepl("-- dropped", warnings, fixed = TRUE)))
    expect_equal(completed_data(result)[c("extra1", "extra2"), , drop = FALSE],
                 d[c("extra1", "extra2"), , drop = FALSE])
    expect_true(is.finite(completed_data(result)$y[1]))
  }
})
