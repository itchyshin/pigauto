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

test_that("pull_gbif_centroids returns data.frame with expected columns from cache", {
  tmp <- tempfile("gbif_cache_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  sp_list <- c("Quercus alba", "Pinus taeda", "Acer saccharum")
  for (sp in sp_list) {
    key <- pigauto:::.gbif_cache_key(sp)
    saveRDS(list(species = sp, centroid_lat = runif(1, 30, 50),
                  centroid_lon = runif(1, -100, -70),
                  n_occurrences = sample(100:500, 1)),
             file.path(tmp, paste0(key, ".rds")))
  }
  out <- pigauto::pull_gbif_centroids(sp_list,
                                        cache_dir = tmp,
                                        verbose = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  expect_named(out,
               c("species", "centroid_lat", "centroid_lon", "n_occurrences"))
  expect_equal(rownames(out), sp_list)
  expect_true(all(is.finite(out$centroid_lat)))
  expect_true(all(out$n_occurrences > 0))
})

test_that("pull_gbif_centroids handles species with no GBIF hits gracefully", {
  tmp <- tempfile("gbif_cache_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  sp <- "Nonexistent species"
  key <- pigauto:::.gbif_cache_key(sp)
  saveRDS(list(species = sp, centroid_lat = NA_real_,
                centroid_lon = NA_real_, n_occurrences = 0L),
           file.path(tmp, paste0(key, ".rds")))
  out <- pigauto::pull_gbif_centroids(c(sp), cache_dir = tmp, verbose = FALSE)
  expect_equal(nrow(out), 1L)
  expect_true(is.na(out$centroid_lat[1]))
  expect_true(is.na(out$centroid_lon[1]))
  expect_equal(out$n_occurrences[1], 0L)
})

# ---- Task 1 (v1.1): store_points arg ----

test_that("pull_gbif_centroids(store_points = TRUE) writes points field", {
  tmp <- tempfile("gbif_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  # Mock rgbif: pretend 5 occurrences for one species
  fake_records <- data.frame(
    decimalLatitude  = c(40, 41, 39, 40.5, 41.5),
    decimalLongitude = c(-82, -81, -83, -82.5, -81.5),
    hasGeospatialIssues = rep(FALSE, 5),
    basisOfRecord = rep("PRESERVED_SPECIMEN", 5),
    stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(
    name_backbone = function(name, ...)
      list(usageKey = 1L, matchType = "EXACT"),
    occ_search = function(...) list(data = fake_records),
    .package = "rgbif")
  out <- pigauto::pull_gbif_centroids(
    species = "TestSpecies one",
    cache_dir = tmp,
    sleep_ms = 0L,
    verbose = FALSE,
    occurrence_limit = 5L,
    store_points = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$n_occurrences, 5L)
  # Confirm the cache RDS has the points field
  key <- pigauto:::.gbif_cache_key("TestSpecies one")
  cached <- readRDS(file.path(tmp, paste0(key, ".rds")))
  expect_true(!is.null(cached$points))
  expect_s3_class(cached$points, "data.frame")
  expect_equal(nrow(cached$points), 5L)
  expect_named(cached$points, c("lat", "lon"), ignore.order = TRUE)
  expect_true(all(cached$points$lat %in% fake_records$decimalLatitude))
})

test_that("pull_gbif_centroids(store_points = FALSE) omits points (default)", {
  tmp <- tempfile("gbif_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fake_records <- data.frame(
    decimalLatitude  = c(40, 41, 39),
    decimalLongitude = c(-82, -81, -83),
    hasGeospatialIssues = rep(FALSE, 3),
    basisOfRecord = rep("PRESERVED_SPECIMEN", 3),
    stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(
    name_backbone = function(name, ...)
      list(usageKey = 1L, matchType = "EXACT"),
    occ_search = function(...) list(data = fake_records),
    .package = "rgbif")
  out <- pigauto::pull_gbif_centroids(
    species = "TestSpecies two",
    cache_dir = tmp,
    sleep_ms = 0L,
    verbose = FALSE)   # store_points defaults to FALSE
  key <- pigauto:::.gbif_cache_key("TestSpecies two")
  cached <- readRDS(file.path(tmp, paste0(key, ".rds")))
  expect_true(is.null(cached$points))
})

test_that("legacy GBIF cache (no points field) still reads cleanly", {
  tmp <- tempfile("gbif_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  # Hand-write legacy cache (pre-v1.1 shape — no points field)
  key <- pigauto:::.gbif_cache_key("Legacy species")
  saveRDS(list(species = "Legacy species",
                centroid_lat = 40.0, centroid_lon = -82.0,
                n_occurrences = 25L,
                fetched_at = Sys.time(),
                match_type = "EXACT"),
           file.path(tmp, paste0(key, ".rds")))
  out <- pigauto::pull_gbif_centroids(
    species = "Legacy species",
    cache_dir = tmp,
    verbose = FALSE)
  expect_equal(out$centroid_lat, 40.0)
  expect_equal(out$centroid_lon, -82.0)
  expect_equal(out$n_occurrences, 25L)
})

test_that("end-to-end: GBIF centroids don't regress plants RMSE (>=1.10x guardrail on env traits)", {
  skip_if(Sys.getenv("NOT_CRAN") == "", "slow integration - set NOT_CRAN=TRUE to run")
  fx_path <- file.path(testthat::test_path("fixtures"), "gbif_plants_300.rds")
  if (!file.exists(fx_path))
    fx_path <- "tests/testthat/fixtures/gbif_plants_300.rds"
  skip_if_not(file.exists(fx_path), "GBIF fixture not found")

  cache_trait <- "script/data-cache/bien_trait_means.rds"
  cache_tree  <- "script/data-cache/bien_tree.rds"
  if (!file.exists(cache_trait)) {
    cache_trait <- file.path("/Users/z3437171/Dropbox/Github Local/pigauto",
                              "script/data-cache", "bien_trait_means.rds")
    cache_tree  <- file.path("/Users/z3437171/Dropbox/Github Local/pigauto",
                              "script/data-cache", "bien_tree.rds")
  }
  skip_if_not(file.exists(cache_trait) && file.exists(cache_tree),
              "BIEN cache not found")

  cov_all <- readRDS(fx_path)
  trait_means <- readRDS(cache_trait)
  tree_raw    <- readRDS(cache_tree)
  tree_all <- if (is.list(tree_raw) && !inherits(tree_raw, "phylo")) {
    t <- tree_raw$scenario.3; t$tip.label <- gsub("_", " ", t$tip.label); t
  } else tree_raw

  all_species <- Reduce(union, lapply(trait_means,
    function(d) if (!is.null(d)) d$species else character(0)))
  wide <- data.frame(species = all_species, stringsAsFactors = FALSE)
  for (nm in names(trait_means)) {
    d <- trait_means[[nm]]
    if (is.null(d)) { wide[[nm]] <- NA_real_; next }
    m <- match(wide$species, d$species)
    wide[[nm]] <- suppressWarnings(as.numeric(d$mean_value[m]))
  }

  # Intersect: has cov + has tree tip + has >= 1 trait
  has_cov_sp <- rownames(cov_all)[!is.na(cov_all$centroid_lat)]
  matched <- Reduce(intersect,
                     list(wide$species, tree_all$tip.label, has_cov_sp))
  skip_if_not(length(matched) >= 100L, "insufficient matched species")

  set.seed(2026L)
  sp_s <- sample(matched, min(200L, length(matched)))
  wide_s <- wide[wide$species %in% sp_s, , drop = FALSE]
  rownames(wide_s) <- wide_s$species; wide_s$species <- NULL
  cov_s <- cov_all[sp_s, c("centroid_lat", "centroid_lon"), drop = FALSE]
  tree_s <- ape::keep.tip(tree_all, sp_s)

  df <- wide_s
  cont_cols <- colnames(df)
  mask <- matrix(FALSE, nrow = nrow(df), ncol = length(cont_cols),
                  dimnames = list(NULL, cont_cols))
  for (v in cont_cols) {
    ok <- which(!is.na(wide_s[[v]]))
    if (length(ok) < 20L) next
    idx <- sample(ok, round(0.30 * length(ok)))
    mask[idx, v] <- TRUE
    df[[v]][idx] <- NA_real_
  }

  # Fit without cov
  res_nocov <- pigauto::impute(df, tree_s,
                                  epochs = 80L, n_imputations = 20L,
                                  verbose = FALSE, seed = 2026L)
  # Fit with GBIF centroids as cov
  res_withcov <- pigauto::impute(df, tree_s, covariates = cov_s,
                                    epochs = 80L, n_imputations = 20L,
                                    verbose = FALSE, seed = 2026L)

  env_traits <- intersect(c("sla", "leaf_area"), cont_cols)
  skip_if(length(env_traits) == 0L, "no environment-driven trait present")

  any_within_guardrail <- FALSE
  for (v in env_traits) {
    if (!any(mask[, v])) next
    truth <- wide_s[[v]][mask[, v]]
    ok <- is.finite(truth)
    if (sum(ok) < 5L) next
    r_nocov   <- sqrt(mean((res_nocov$completed[[v]][mask[, v]][ok] - truth[ok])^2))
    r_withcov <- sqrt(mean((res_withcov$completed[[v]][mask[, v]][ok] - truth[ok])^2))
    ratio <- r_withcov / r_nocov
    cat(sprintf("[end-to-end] %s: no-cov RMSE = %.3g, with-cov = %.3g, ratio = %.4f\n",
                 v, r_nocov, r_withcov, ratio))
    if (ratio <= 1.10) any_within_guardrail <- TRUE
  }

  expect_true(any_within_guardrail,
              info = "at least one of sla/leaf_area should be within 10% of no-cov RMSE (centroids not catastrophically hurting)")
})
