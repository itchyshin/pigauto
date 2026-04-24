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

test_that(".wc_download_rasters is a no-op when sentinel + rasters present", {
  tmp <- tempfile("wc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  wc_dir <- file.path(tmp, "wc2.1_10m")
  dir.create(wc_dir)
  # Create fake rasters (empty files, only presence is checked)
  for (i in 1:19) {
    file.create(file.path(wc_dir, sprintf("wc2.1_10m_bio_%d.tif", i)))
  }
  file.create(file.path(wc_dir, ".wc_complete"))
  # Should return wc_dir without attempting download
  out <- pigauto:::.wc_download_rasters(tmp, resolution = "10m",
                                          verbose = FALSE)
  expect_equal(normalizePath(out), normalizePath(wc_dir))
})

test_that(".wc_download_rasters detects missing sentinel and would re-download", {
  tmp <- tempfile("wc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  # Sentinel absent -> would download; we only test the intent, not
  # the actual HTTP call.
  expect_error(
    pigauto:::.wc_download_rasters(tmp, resolution = "10m",
                                      verbose = FALSE,
                                      .download_fn = function(...)
                                        stop("MOCK: would download")),
    "MOCK: would download")
})

test_that(".wc_download_rasters errors cleanly for unsupported resolution", {
  tmp <- tempfile("wc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  expect_error(
    pigauto:::.wc_download_rasters(tmp, resolution = "99m", verbose = FALSE),
    "resolution must be one of")
})

test_that(".wc_extract_one reads from cache when present", {
  tmp <- tempfile("wc_cache_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  extracts_dir <- file.path(tmp, "extracts")
  dir.create(extracts_dir)
  key <- pigauto:::.wc_cache_key("Quercus alba")
  med <- setNames(runif(19, 0, 20), paste0("bio", 1:19))
  iqr <- setNames(runif(19, 0, 5),  paste0("bio", 1:19))
  saveRDS(list(species = "Quercus alba",
                bio_median = med, bio_iqr = iqr, n_extracted = 42L,
                extracted_at = Sys.time()),
           file.path(extracts_dir, paste0(key, ".rds")))
  out <- pigauto:::.wc_extract_one(
    sp = "Quercus alba",
    gbif_cache_dir = tempdir(),  # won't be read on cache hit
    wc_extracts_dir = extracts_dir,
    rast_stack = NULL,
    refresh_cache = FALSE)
  expect_equal(out$n_extracted, 42L)
  expect_equal(out$bio_median, med)
})

test_that(".wc_extract_one returns NA row when species has no GBIF cache", {
  tmp <- tempfile("wc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gbif_dir <- file.path(tmp, "gbif"); dir.create(gbif_dir)
  extracts_dir <- file.path(tmp, "extracts"); dir.create(extracts_dir)
  # No GBIF cache for "Ghost species"
  out <- pigauto:::.wc_extract_one(
    sp = "Ghost species",
    gbif_cache_dir = gbif_dir,
    wc_extracts_dir = extracts_dir,
    rast_stack = NULL,
    refresh_cache = FALSE)
  expect_equal(out$n_extracted, 0L)
  expect_true(all(is.na(out$bio_median)))
})
