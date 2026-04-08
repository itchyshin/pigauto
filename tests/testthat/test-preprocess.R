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


# ---- Multi-observation per species tests ------------------------------------

test_that("preprocess_traits with species_col single-obs", {
  set.seed(10)
  tree <- ape::rtree(10)
  df <- data.frame(
    sp = tree$tip.label,
    mass = abs(rnorm(10)) + 0.1
  )
  pd <- preprocess_traits(df, tree, species_col = "sp")
  expect_true(pd$multi_obs)
  expect_equal(pd$n_obs, 10L)
  expect_equal(pd$n_species, 10L)
  expect_equal(nrow(pd$X_scaled), 10)
  expect_equal(length(pd$obs_to_species), 10)
  expect_equal(range(pd$obs_to_species), c(1L, 10L))
})

test_that("preprocess_traits with species_col multi-obs", {
  set.seed(11)
  tree <- ape::rtree(5)
  # 3 obs per species = 15 rows
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    mass = abs(rnorm(15)) + 0.1
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  expect_true(pd$multi_obs)
  expect_equal(pd$n_obs, 15L)
  expect_equal(pd$n_species, 5L)
  expect_equal(nrow(pd$X_scaled), 15)
  expect_equal(length(pd$obs_to_species), 15)
  expect_equal(length(pd$obs_species), 15)
  expect_equal(length(pd$species_names), 5)
  # Each species should map to unique index
  sp_idx <- tapply(pd$obs_to_species, pd$obs_species, unique)
  expect_true(all(lengths(sp_idx) == 1))
})

test_that("preprocess_traits backward compatible (no species_col, no multi_obs)", {
  set.seed(12)
  tree <- ape::rtree(8)
  df <- data.frame(row.names = tree$tip.label, tr = rnorm(8))
  pd <- preprocess_traits(df, tree)
  expect_false(pd$multi_obs)
  expect_null(pd$obs_species)
  expect_null(pd$obs_to_species)
  expect_equal(pd$n_obs, 8L)
  expect_equal(pd$n_species, 8L)
})


# ---- avonet300 categorical variables ----------------------------------------

test_that("avonet300 contains categorical and ordinal traits", {
  data(avonet300, package = "pigauto")
  expect_true("Trophic.Level" %in% names(avonet300))
  expect_true("Primary.Lifestyle" %in% names(avonet300))
  expect_true("Migration" %in% names(avonet300))
  expect_true(is.factor(avonet300$Trophic.Level))
  expect_true(is.factor(avonet300$Primary.Lifestyle))
  expect_true(is.ordered(avonet300$Migration))
  expect_equal(nlevels(avonet300$Trophic.Level), 4)
  expect_equal(nlevels(avonet300$Primary.Lifestyle), 5)
  expect_equal(nlevels(avonet300$Migration), 3)
})

test_that("avonet300 preprocesses with mixed types", {
  data(avonet300, tree300, package = "pigauto")
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  pd <- preprocess_traits(df, tree300)
  expect_s3_class(pd, "pigauto_data")
  types <- vapply(pd$trait_map, "[[", character(1), "type")
  expect_true("categorical" %in% types)
  expect_true("ordinal" %in% types)
  expect_true("continuous" %in% types)
  # 4 continuous + 4 (Trophic.Level one-hot) + 5 (Primary.Lifestyle one-hot)
  # + 1 (Migration ordinal) = 14
  expect_equal(pd$p_latent, 14L)
})
