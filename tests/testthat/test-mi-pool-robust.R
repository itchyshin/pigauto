# Regression tests for the MI-pool robustness fix:
# log-transformed continuous + count + zi_count + proportion now use
# median pooling across decoded values, so a single dropout-noisy draw
# with a high-tail latent value (which decodes via expm1 / plogis to
# an absurd magnitude) cannot blow up the pooled pigauto prediction.
#
# Repro for the bug that motivated this: FishBase n=10,654 bench
# (2026-04-21) with n_imputations=20 produced Length RMSE 3718 vs
# mean baseline 41 -- a decode artefact, not a modelling failure.

test_that("pool_imputations uses median for log-transformed continuous", {
  # Construct three "imputation draws" for one trait with one outlier
  # cell in draw 2. n=4 species, p=1 trait.
  sp <- paste0("s", 1:4)
  draws <- list(
    list(imputed = data.frame(row.names = sp, mass = c(10, 20, 30, 40))),
    list(imputed = data.frame(row.names = sp, mass = c(12, 22, 30, 1e6))),
    list(imputed = data.frame(row.names = sp, mass = c(11, 21, 29, 38)))
  )
  latents <- list(matrix(0, nrow = 4, ncol = 1),
                   matrix(0, nrow = 4, ncol = 1),
                   matrix(0, nrow = 4, ncol = 1))
  # log-transform = TRUE triggers median pooling for continuous
  trait_map <- list(list(
    name = "mass", type = "continuous",
    log_transform = TRUE,
    latent_cols = 1L
  ))

  # Call the internal pool_imputations directly
  pool_fn <- getFromNamespace("pool_imputations", "pigauto")
  out <- pool_fn(draws, latents, trait_map)

  # Median of {40, 1e6, 38} = 40. Mean would be ~3.3e5.
  expect_equal(as.numeric(out$imputed$mass[4]), 40)
  expect_lt(as.numeric(out$imputed$mass[4]), 100)   # sanity: not in outlier range
})

test_that("pool_imputations uses mean for untransformed continuous", {
  # Without log_transform, mean pooling is the backward-compat default.
  sp <- paste0("s", 1:3)
  draws <- list(
    list(imputed = data.frame(row.names = sp, x = c(1, 2, 3))),
    list(imputed = data.frame(row.names = sp, x = c(3, 4, 5))),
    list(imputed = data.frame(row.names = sp, x = c(5, 6, 7)))
  )
  latents <- list(matrix(0, nrow = 3, ncol = 1),
                   matrix(0, nrow = 3, ncol = 1),
                   matrix(0, nrow = 3, ncol = 1))
  trait_map <- list(list(
    name = "x", type = "continuous",
    log_transform = FALSE,
    latent_cols = 1L
  ))
  pool_fn <- getFromNamespace("pool_imputations", "pigauto")
  out <- pool_fn(draws, latents, trait_map)
  # Mean of {1,3,5} = 3, {2,4,6} = 4, {3,5,7} = 5
  expect_equal(as.numeric(out$imputed$x), c(3, 4, 5))
})

test_that("pool_imputations uses median for count traits", {
  # Matches Issue #40 / PR #41 rationale extended to all code paths.
  sp <- paste0("s", 1:3)
  draws <- list(
    list(imputed = data.frame(row.names = sp, n = c(5, 10, 15))),
    list(imputed = data.frame(row.names = sp, n = c(6, 11, 2000))),
    list(imputed = data.frame(row.names = sp, n = c(4, 9, 14)))
  )
  latents <- list(matrix(0, nrow = 3, ncol = 1),
                   matrix(0, nrow = 3, ncol = 1),
                   matrix(0, nrow = 3, ncol = 1))
  trait_map <- list(list(
    name = "n", type = "count",
    log_transform = TRUE,
    latent_cols = 1L
  ))
  pool_fn <- getFromNamespace("pool_imputations", "pigauto")
  out <- pool_fn(draws, latents, trait_map)
  # Row 3: median of {15, 2000, 14} = 15. Mean would be ~676.
  expect_equal(as.integer(out$imputed$n[3]), 15L)
})

test_that("pool_imputations uses median for proportion traits", {
  # Proportion decoded via plogis can amplify near-extreme latents.
  sp <- paste0("s", 1:3)
  draws <- list(
    list(imputed = data.frame(row.names = sp, p = c(0.3, 0.5, 0.7))),
    list(imputed = data.frame(row.names = sp, p = c(0.31, 0.51, 0.999))),
    list(imputed = data.frame(row.names = sp, p = c(0.29, 0.49, 0.69)))
  )
  latents <- list(matrix(0, nrow = 3, ncol = 1),
                   matrix(0, nrow = 3, ncol = 1),
                   matrix(0, nrow = 3, ncol = 1))
  trait_map <- list(list(
    name = "p", type = "proportion",
    log_transform = FALSE,
    latent_cols = 1L
  ))
  pool_fn <- getFromNamespace("pool_imputations", "pigauto")
  out <- pool_fn(draws, latents, trait_map)
  # Row 3: median of {0.7, 0.999, 0.69} = 0.7. Mean would be ~0.80.
  expect_equal(as.numeric(out$imputed$p[3]), 0.7)
})
