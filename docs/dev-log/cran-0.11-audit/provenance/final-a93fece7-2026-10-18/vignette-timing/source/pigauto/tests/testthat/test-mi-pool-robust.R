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

# ---- Phase H (2026-05-01): mode pooling for ordinal -----------------------

test_that("[Phase H] pool_imputations 'mode' picks majority class for ordinal", {
  # Construct M=5 imputation draws for K=3 ordinal trait with 4 species.
  # Truth pattern: each species' "true" class is the most frequent across
  # draws; default mean+round can flip species 1 from class 1 to class 2
  # because mean(1,1,1,2,3)=1.6 rounds to 2; mode picks 1.
  sp <- paste0("s", 1:4)
  L  <- c("low", "med", "high")
  mk <- function(v) ordered(L[v], levels = L)
  draws <- list(
    list(imputed = data.frame(row.names = sp, m = mk(c(1, 2, 3, 1)))),
    list(imputed = data.frame(row.names = sp, m = mk(c(1, 2, 3, 1)))),
    list(imputed = data.frame(row.names = sp, m = mk(c(1, 2, 3, 2)))),
    list(imputed = data.frame(row.names = sp, m = mk(c(2, 2, 3, 2)))),
    list(imputed = data.frame(row.names = sp, m = mk(c(3, 2, 3, 2))))
  )
  latents <- replicate(5L, matrix(0, nrow = 4, ncol = 1), simplify = FALSE)
  trait_map <- list(list(
    name = "m", type = "ordinal",
    levels = L,
    latent_cols = 1L
  ))

  pool_fn <- getFromNamespace("pool_imputations", "pigauto")

  # MODE: per-cell majority vote
  out_mode <- pool_fn(draws, latents, trait_map, pool_method = "mode")
  expect_equal(as.character(out_mode$imputed$m[1]), "low",
               label = "[Phase H] species 1: 3 votes for 'low' should win")
  expect_equal(as.character(out_mode$imputed$m[2]), "med",
               label = "[Phase H] species 2: all 5 'med', mode 'med'")
  expect_equal(as.character(out_mode$imputed$m[3]), "high",
               label = "[Phase H] species 3: all 5 'high', mode 'high'")
  expect_equal(as.character(out_mode$imputed$m[4]), "med",
               label = "[Phase H] species 4: 3 votes for 'med' should win")

  # MEDIAN (default): mean+round.  Verify the documented "wrong" cases
  # to lock in the regression-test pattern.
  out_median <- pool_fn(draws, latents, trait_map, pool_method = "median")
  # Species 1 votes (1,1,1,2,3) on integer scale: mean=1.6, round=2 -> "med"
  # vs mode "low".  This documents the failure mode that mode fixes.
  expect_equal(as.character(out_median$imputed$m[1]), "med",
               label = "[Phase H] mean+round biases species 1 to 'med' (mode would pick 'low')")
})

test_that("[Phase H] mode falls back to median for log-transformed continuous", {
  # When pool_method = "mode" is passed for a continuous-family trait,
  # the function should treat it like "median" (since mode does not
  # apply to continuous decoders).
  sp <- paste0("s", 1:4)
  draws <- list(
    list(imputed = data.frame(row.names = sp, mass = c(10, 20, 30, 40))),
    list(imputed = data.frame(row.names = sp, mass = c(12, 22, 30, 1e6))),
    list(imputed = data.frame(row.names = sp, mass = c(11, 21, 29, 38)))
  )
  latents <- replicate(3L, matrix(0, nrow = 4, ncol = 1), simplify = FALSE)
  trait_map <- list(list(
    name = "mass", type = "continuous",
    log_transform = TRUE, latent_cols = 1L
  ))
  pool_fn <- getFromNamespace("pool_imputations", "pigauto")
  out <- pool_fn(draws, latents, trait_map, pool_method = "mode")
  # Median of {40, 1e6, 38} = 40 -- mode falls back to median, NOT mean.
  expect_equal(as.numeric(out$imputed$mass[4]), 40,
               label = "[Phase H] mode pool method must fall back to median for log-cont")
})
