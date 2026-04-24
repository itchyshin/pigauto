# tests/testthat/test-worldclim-covariates.R
# Smoke canary for WorldClim per-species bioclim covariate helper.
# See specs/2026-04-24-worldclim-covariates-design.md.

test_that(".wc_cache_key sanitises species names for filesystem", {
  expect_equal(pigauto:::.wc_cache_key("Quercus alba"),
               "Quercus_alba")
  expect_equal(pigauto:::.wc_cache_key("Quercus alba L."),
               "Quercus_alba_L_")
  expect_equal(pigauto:::.wc_cache_key("Pinus x contorta"),
               "Pinus_x_contorta")
  # Empty input: defensive default
  expect_error(pigauto:::.wc_cache_key(""),
               "must be a non-empty character scalar")
})

test_that(".wc_aggregate_one computes median + IQR per column", {
  vals <- data.frame(
    bio1 = c(10, 12, 11, 15, NA),
    bio2 = c(5,  6,  7,  8,  9),
    stringsAsFactors = FALSE)
  out <- pigauto:::.wc_aggregate_one(vals)
  expect_type(out, "list")
  expect_named(out,
               c("bio_median", "bio_iqr", "n_extracted"),
               ignore.order = TRUE)
  # bio1: median of {10, 12, 11, 15} = 11.5 (NA dropped);
  # IQR = quantile(., 0.75) - quantile(., 0.25)
  expect_equal(out$bio_median[["bio1"]], 11.5)
  expect_equal(out$bio_median[["bio2"]], 6.5)  # row 5 (bio1=NA) dropped: {5,6,7,8}
  expect_equal(out$n_extracted, 4L)  # only rows with all-non-NA counted
  expect_gte(out$bio_iqr[["bio1"]], 0)
  expect_gte(out$bio_iqr[["bio2"]], 0)
})

test_that(".wc_aggregate_one returns NA bio values when all rows NA", {
  vals <- data.frame(bio1 = rep(NA_real_, 3),
                      bio2 = rep(NA_real_, 3),
                      stringsAsFactors = FALSE)
  out <- pigauto:::.wc_aggregate_one(vals)
  expect_equal(out$n_extracted, 0L)
  expect_true(all(is.na(out$bio_median)))
  expect_true(all(is.na(out$bio_iqr)))
})

test_that(".wc_aggregate_one handles single-point species (IQR = 0)", {
  vals <- data.frame(bio1 = 15, bio2 = 20, stringsAsFactors = FALSE)
  out <- pigauto:::.wc_aggregate_one(vals)
  expect_equal(out$n_extracted, 1L)
  expect_equal(out$bio_median[["bio1"]], 15)
  expect_equal(out$bio_iqr[["bio1"]], 0)   # IQR of a singleton is 0
})
