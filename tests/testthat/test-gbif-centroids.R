# tests/testthat/test-gbif-centroids.R
# Smoke canary for GBIF centroid helper.
# See specs/2026-04-23-gbif-centroids-design.md.

test_that(".gbif_cache_key sanitises species names for filesystem", {
  expect_equal(pigauto:::.gbif_cache_key("Quercus alba"),
               "Quercus_alba")
  expect_equal(pigauto:::.gbif_cache_key("Quercus alba L."),
               "Quercus_alba_L_")
  expect_equal(pigauto:::.gbif_cache_key("Pinus x contorta"),
               "Pinus_x_contorta")
  # Unicode handling: converted to underscore-safe ASCII chars
  expect_true(nzchar(pigauto:::.gbif_cache_key("Abies Sp\u00e9cies")))
})

test_that(".gbif_centroid_one aggregates valid records to median lat/lon", {
  records <- data.frame(
    decimalLatitude  = c(10, 12, 11, 15, NA),
    decimalLongitude = c(-80, -82, -81, -85, -83),
    hasGeospatialIssues = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    basisOfRecord    = rep("PRESERVED_SPECIMEN", 5),
    stringsAsFactors = FALSE)
  out <- pigauto:::.gbif_centroid_one(records)
  expect_type(out, "list")
  expect_named(out, c("centroid_lat", "centroid_lon", "n_occurrences"),
               ignore.order = TRUE)
  expect_equal(out$centroid_lat, 11.5)    # median of 10, 11, 12, 15
  expect_equal(out$centroid_lon, -81.5)   # median of -85, -82, -81, -80 = (-82 + -81)/2
  expect_equal(out$n_occurrences, 4L)     # one row with NA lat excluded
})

test_that(".gbif_centroid_one filters geospatial-issue + fossil records", {
  records <- data.frame(
    decimalLatitude  = c(10, 12, 90, 11),
    decimalLongitude = c(-80, -82, -81, -85),
    hasGeospatialIssues = c(FALSE, FALSE, TRUE,  FALSE),
    basisOfRecord    = c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION",
                           "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN"),
    stringsAsFactors = FALSE)
  out <- pigauto:::.gbif_centroid_one(records)
  # geospatial-issue row (90, -81) and fossil row (11, -85) dropped
  expect_equal(out$n_occurrences, 2L)
  expect_equal(out$centroid_lat, 11.0)    # median of 10, 12
  expect_equal(out$centroid_lon, -81.0)   # median of -80, -82
})

test_that(".gbif_centroid_one returns NA when no valid records", {
  out <- pigauto:::.gbif_centroid_one(
    data.frame(decimalLatitude = numeric(0), decimalLongitude = numeric(0),
               hasGeospatialIssues = logical(0),
               basisOfRecord = character(0),
               stringsAsFactors = FALSE))
  expect_true(is.na(out$centroid_lat))
  expect_true(is.na(out$centroid_lon))
  expect_equal(out$n_occurrences, 0L)
})

test_that(".gbif_centroid_one rejects coords out of range", {
  records <- data.frame(
    decimalLatitude  = c(10, 91, 11),
    decimalLongitude = c(-80, -82, -200),
    hasGeospatialIssues = rep(FALSE, 3),
    basisOfRecord    = rep("PRESERVED_SPECIMEN", 3),
    stringsAsFactors = FALSE)
  out <- pigauto:::.gbif_centroid_one(records)
  # (91, -82) and (11, -200) both out of range; only (10, -80) valid
  expect_equal(out$n_occurrences, 1L)
  expect_equal(out$centroid_lat, 10)
  expect_equal(out$centroid_lon, -80)
})

test_that(".gbif_fetch_one uses cache when present", {
  tmp <- tempfile("gbif_cache_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  key <- pigauto:::.gbif_cache_key("Quercus alba")
  saveRDS(list(species = "Quercus alba",
                centroid_lat = 41.5, centroid_lon = -80.0,
                n_occurrences = 120L,
                fetched_at = as.POSIXct("2026-01-01 00:00:00", tz = "UTC")),
           file.path(tmp, paste0(key, ".rds")))
  out <- pigauto:::.gbif_fetch_one("Quercus alba",
                                      cache_dir = tmp,
                                      occurrence_limit = 500L,
                                      sleep_ms = 0L,
                                      refresh_cache = FALSE)
  expect_equal(out$centroid_lat, 41.5)
  expect_equal(out$centroid_lon, -80.0)
  expect_equal(out$n_occurrences, 120L)
})

test_that(".gbif_fetch_one reports helpful error when rgbif missing", {
  testthat::skip_if(requireNamespace("rgbif", quietly = TRUE),
                     "rgbif IS installed; cannot test missing-rgbif branch")
  tmp <- tempfile("gbif_cache_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  expect_error(
    pigauto:::.gbif_fetch_one("Quercus alba",
                                cache_dir = tmp,
                                occurrence_limit = 500L,
                                sleep_ms = 0L,
                                refresh_cache = FALSE),
    "rgbif")
})
