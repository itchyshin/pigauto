test_that("read_traits returns data.frame with row names", {
  df <- data.frame(species = c("Sp_a", "Sp_b"), mass = c(10.0, 20.0))
  out <- read_traits(df, species_col = "species")
  expect_s3_class(out, "data.frame")
  expect_equal(rownames(out), c("Sp_a", "Sp_b"))
  expect_named(out, "mass")
})

test_that("read_traits errors on missing species column", {
  df <- data.frame(name = c("A", "B"), x = c(1, 2))
  expect_error(read_traits(df, species_col = "species"), "not found")
})

test_that("preprocess_traits output has correct dimensions", {
  set.seed(1)
  tree <- ape::rtree(10)
  sp   <- tree$tip.label
  df   <- data.frame(
    row.names = sp,
    tr1 = abs(rnorm(10)) + 0.1,
    tr2 = abs(rnorm(10)) + 0.1
  )
  pd <- preprocess_traits(df, tree)
  expect_equal(nrow(pd$X_scaled), 10)
  expect_equal(ncol(pd$X_scaled), 2)
  expect_equal(pd$species_names, tree$tip.label)
  expect_s3_class(pd, "pigauto_data")
})

test_that("preprocess_traits rows match tree tip order", {
  set.seed(2)
  tree <- ape::rtree(8)
  sp   <- sample(tree$tip.label)
  df   <- data.frame(row.names = sp, tr = abs(rnorm(8)) + 0.1)
  pd <- preprocess_traits(df, tree)
  expect_equal(pd$species_names, tree$tip.label)
})

test_that("preprocess_traits with negative values skips log-transform", {
  set.seed(3)
  tree <- ape::rtree(5)
  df   <- data.frame(row.names = tree$tip.label, tr = c(-1, 1, 2, 3, 4))
  pd   <- preprocess_traits(df, tree, log_transform = TRUE)
  # Negative values prevent auto-log; no error thrown
  expect_false(pd$trait_map$tr$log_transform)
})

test_that("preprocess_traits warns about species not in tree", {
  set.seed(4)
  tree <- ape::rtree(5)
  sp   <- c(tree$tip.label, "Ghost_species")
  df   <- data.frame(row.names = sp, tr = abs(rnorm(6)) + 0.1)
  expect_warning(preprocess_traits(df, tree), "not found in tree")
})

test_that("preprocess_traits creates trait_map", {
  set.seed(5)
  tree <- ape::rtree(10)
  df   <- data.frame(
    row.names = tree$tip.label,
    mass = abs(rnorm(10)) + 0.1,
    count_trait = as.integer(rpois(10, 3))
  )
  pd <- preprocess_traits(df, tree)
  expect_true(!is.null(pd$trait_map))
  expect_equal(length(pd$trait_map), 2L)
  expect_equal(pd$trait_map$mass$type, "continuous")
  expect_equal(pd$trait_map$count_trait$type, "count")
})

test_that("preprocess_traits handles factor columns", {
  set.seed(6)
  tree <- ape::rtree(10)
  df <- data.frame(
    row.names = tree$tip.label,
    diet = factor(sample(c("herbivore", "carnivore", "omnivore"), 10,
                         replace = TRUE)),
    migr = factor(sample(c("yes", "no"), 10, replace = TRUE))
  )
  pd <- preprocess_traits(df, tree)
  expect_equal(pd$trait_map$diet$type, "categorical")
  expect_equal(pd$trait_map$migr$type, "binary")
  expect_equal(pd$trait_map$diet$n_latent, 3L)
  expect_equal(pd$trait_map$migr$n_latent, 1L)
  expect_equal(pd$p_latent, 4L)  # 3 (diet one-hot) + 1 (migr binary)
})
