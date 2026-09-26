# The ordinal label-propagation candidate must decode the 0..K-1 class coding
# used by preprocess_traits() (a Phase F bug decoded 1..K).

test_that("[ordinal-lp] all observed species in the lowest class predict the lowest class", {
  tm <- list(levels = c("lo", "mid", "hi"), mean = 1, sd = 1)
  z <- c(-1, -1, -1, NA)                    # classes 0, 0, 0, missing
  sim <- matrix(1, 4, 4)                    # everyone equally similar
  out <- pigauto:::.ordinal_lp_candidate(z, sim, tm)
  expect_false(is.null(out))
  cls <- round(out$pred_z * tm$sd + tm$mean)
  expect_true(all(cls == 0))
})

test_that("[ordinal-lp] the candidate uses every class and reproduces a known mixture", {
  tm <- list(levels = c("a", "b", "c"), mean = 1, sd = 0.5)
  cls <- c(0, 1, 2, 2, NA)
  z <- (cls - tm$mean) / tm$sd
  sim <- matrix(1, 5, 5)
  out <- pigauto:::.ordinal_lp_candidate(z, sim, tm)
  # Equal weights over observed classes 0, 1, 2, 2: E[class] = 1.25.
  expect_equal(out$pred_z[5] * tm$sd + tm$mean, 1.25, tolerance = 1e-4)
  expect_equal(unname(out$se_z[5] * tm$sd), sqrt(mean((c(0, 1, 2, 2) - 1.25)^2)), tolerance = 1e-4)
})

test_that("[ordinal-lp] returns NULL with fewer observed cells than classes", {
  tm <- list(levels = c("a", "b", "c"), mean = 1, sd = 1)
  expect_null(pigauto:::.ordinal_lp_candidate(c(-1, 0, NA, NA), matrix(1, 4, 4), tm))
})
