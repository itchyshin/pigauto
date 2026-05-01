# Phase G (2026-05-01): tests for the optional clamp_outliers argument.
#
# Targets the Mass tail-extrapolation mode documented in
# useful/MEMO_2026-05-01_avonet_mass_diag.md: a +3-4 SD latent overshoot
# becomes a 50x-100x value error after expm1().  The clamp caps
# back-transformed predictions at tm$obs_max * clamp_factor (default 5).

test_that("[Phase G] preprocess_traits records obs_max / obs_min for log-cont, count, zi_count", {
  set.seed(2050L)
  n <- 30L
  tree <- ape::rtree(n)
  # Three at-risk types in one fixture
  df <- data.frame(
    mass     = exp(stats::rnorm(n, mean = 5, sd = 1)),       # log-cont
    counts   = as.integer(stats::rpois(n, lambda = 10)),     # count
    zicount  = {
      cnt <- as.integer(stats::rpois(n, lambda = 6))
      cnt[sample.int(n, round(n * 0.3))] <- 0L
      cnt
    },
    plain    = stats::rnorm(n),                              # un-log cont
    row.names = tree$tip.label
  )
  pd <- pigauto::preprocess_traits(
    df, tree,
    trait_types = c(zicount = "zi_count"),
    log_transform = TRUE
  )

  # Look up the trait_map entries by name
  tmap <- pd$trait_map
  by_name <- function(nm) Filter(function(tm) identical(tm$name, nm), tmap)[[1]]

  mass_tm <- by_name("mass")
  expect_true(isTRUE(mass_tm$log_transform))
  expect_true(is.finite(mass_tm$obs_max),
              info = "[Phase G] log-cont must record obs_max")
  expect_equal(mass_tm$obs_max, max(df$mass), tolerance = 1e-9)
  expect_equal(mass_tm$obs_min, min(df$mass), tolerance = 1e-9)

  count_tm <- by_name("counts")
  expect_true(is.finite(count_tm$obs_max),
              info = "[Phase G] count must record obs_max")
  expect_equal(count_tm$obs_max, max(df$counts), tolerance = 1e-9)

  zicount_tm <- by_name("zicount")
  # zi_count clamp uses non-zero observations
  nz_max <- max(df$zicount[df$zicount > 0L])
  nz_min <- min(df$zicount[df$zicount > 0L])
  expect_equal(zicount_tm$obs_max, nz_max, tolerance = 1e-9)
  expect_equal(zicount_tm$obs_min, nz_min, tolerance = 1e-9)

  # Untransformed continuous: NA (not at risk of expm1 amplification)
  plain_tm <- by_name("plain")
  expect_false(isTRUE(plain_tm$log_transform))
  expect_true(is.na(plain_tm$obs_max),
              info = "[Phase G] un-log continuous should NOT record obs_max")
})

test_that("[Phase G] clamp_outliers = FALSE is the default and is a no-op", {
  skip_if_not_installed("torch")
  set.seed(2051L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, mean = 5, sd = 1)),
    row.names = tree$tip.label
  )
  df$mass[c(2L, 5L, 9L)] <- NA

  res_default <- pigauto::impute(df, tree,
                                  epochs = 10L, n_imputations = 1L,
                                  verbose = FALSE, seed = 2051L)
  res_explicit_off <- pigauto::impute(df, tree,
                                       clamp_outliers = FALSE,
                                       epochs = 10L, n_imputations = 1L,
                                       verbose = FALSE, seed = 2051L)
  # Identical predictions on the imputed cells (no-op confirmation)
  miss <- which(is.na(df$mass))
  expect_equal(res_default$completed$mass[miss],
               res_explicit_off$completed$mass[miss],
               tolerance = 1e-9,
               info = "[Phase G] default and clamp_outliers=FALSE must be identical")
})

test_that("[Phase G] clamp_outliers = TRUE caps a synthetic tail outlier", {
  skip_if_not_installed("torch")
  # Reach into the decoder directly with a contrived latent vector that
  # back-transforms to a 50x-the-max blow-up.  Verifies the clamp catches
  # it when ON, and lets it through when OFF.
  set.seed(2052L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, mean = 4, sd = 0.5)),  # range ~ [10, 200]
    row.names = tree$tip.label
  )
  pd <- pigauto::preprocess_traits(df, tree, log_transform = TRUE)
  tm <- pd$trait_map[[1L]]
  # Build a latent matrix where row 1 is z = +6 SD (massive overshoot
  # after exp); other rows are at z=0.
  latent_mat <- matrix(0, nrow = n, ncol = 1L)
  latent_mat[1L, 1L] <- 6
  rownames(latent_mat) <- tree$tip.label

  decode_fn <- getFromNamespace("decode_from_latent", "pigauto")

  # OFF: no clamp; row 1 should be exp(6 * sd + mean), enormous
  out_off <- decode_fn(latent_mat, pd$trait_map, tree$tip.label,
                       clamp_outliers = FALSE)
  # ON: clamp at obs_max * 5
  out_on <- decode_fn(latent_mat, pd$trait_map, tree$tip.label,
                      clamp_outliers = TRUE, clamp_factor = 5)

  expect_gt(out_off$imputed$mass[1L], tm$obs_max * 5,
            label = "[Phase G] OFF: tail prediction must exceed obs_max * 5 (uncapped)")
  expect_equal(out_on$imputed$mass[1L], tm$obs_max * 5,
               tolerance = 1e-6,
               info = "[Phase G] ON: tail prediction must be exactly obs_max * 5")

  # Inside-range rows must be unchanged by the clamp
  inside_idx <- 2:n
  expect_equal(out_off$imputed$mass[inside_idx],
               out_on$imputed$mass[inside_idx],
               tolerance = 1e-9,
               info = "[Phase G] inside-range predictions must be unchanged")
})

test_that("[Phase G] clamp_factor controls the cap multiplier", {
  skip_if_not_installed("torch")
  set.seed(2053L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, mean = 4, sd = 0.5)),
    row.names = tree$tip.label
  )
  pd <- pigauto::preprocess_traits(df, tree, log_transform = TRUE)
  tm <- pd$trait_map[[1L]]
  # Force a guaranteed-to-overshoot prediction by using a large latent
  # value: pre-clamp prediction is exp(20 * sd + mean) which dwarfs any
  # reasonable obs_max * factor.
  latent_mat <- matrix(0, nrow = n, ncol = 1L)
  latent_mat[1L, 1L] <- 20
  rownames(latent_mat) <- tree$tip.label

  decode_fn <- getFromNamespace("decode_from_latent", "pigauto")
  out_2 <- decode_fn(latent_mat, pd$trait_map, tree$tip.label,
                     clamp_outliers = TRUE, clamp_factor = 2)
  out_10 <- decode_fn(latent_mat, pd$trait_map, tree$tip.label,
                      clamp_outliers = TRUE, clamp_factor = 10)

  expect_equal(out_2$imputed$mass[1L], tm$obs_max * 2, tolerance = 1e-6)
  expect_equal(out_10$imputed$mass[1L], tm$obs_max * 10, tolerance = 1e-6)
})

test_that("[Phase G] clamp_factor < 1 errors", {
  skip_if_not_installed("torch")
  set.seed(2054L)
  n <- 20L
  tree <- ape::rtree(n)
  df <- data.frame(mass = exp(stats::rnorm(n)),
                    row.names = tree$tip.label)
  expect_error(
    pigauto::impute(df, tree,
                     clamp_outliers = TRUE, clamp_factor = 0.5,
                     epochs = 5L, n_imputations = 1L,
                     verbose = FALSE, seed = 2054L),
    "clamp_factor"
  )
})
