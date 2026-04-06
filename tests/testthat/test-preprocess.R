test_that("read_traits returns data.frame with row names", {
  df <- data.frame(species = c("Sp_a", "Sp_b"), mass = c(10.0, 20.0))
  out <- read_traits(df, species_col = "species")
  expect_s3_class(out, "data.frame")
  expect_equal(rownames(out), c("Sp_a", "Sp_b"))
  expect_named(out, "mass")
})

test_that("read_traits errors on missing species column", {
  df <- data.frame(name = c("A", "B"), x = c(1, 2))
  expect_snapshot(read_traits(df, species_col = "species"), error = TRUE)
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
})

test_that("preprocess_traits rows match tree tip order", {
  set.seed(2)
  tree <- ape::rtree(8)
  sp   <- sample(tree$tip.label)
  df   <- data.frame(row.names = sp, tr = abs(rnorm(8)) + 0.1)
  pd <- preprocess_traits(df, tree)
  expect_equal(pd$species_names, tree$tip.label)
})

test_that("preprocess_traits stops on non-positive values with log_transform", {
  set.seed(3)
  tree <- ape::rtree(5)
  df   <- data.frame(row.names = tree$tip.label, tr = c(-1, 1, 2, 3, 4))
  expect_snapshot(
    preprocess_traits(df, tree, log_transform = TRUE),
    error = TRUE
  )
})

test_that("preprocess_traits warns about species not in tree", {
  set.seed(4)
  tree <- ape::rtree(5)
  sp   <- c(tree$tip.label, "Ghost_species")
  df   <- data.frame(row.names = sp, tr = abs(rnorm(6)) + 0.1)
  expect_snapshot(preprocess_traits(df, tree))
})
