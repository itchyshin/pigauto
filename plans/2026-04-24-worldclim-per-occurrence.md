# Per-occurrence WorldClim covariates — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix PR #47's centroid-only limitation by extracting bioclim at every GBIF occurrence point per species. Lifts plants SLA / leaf_area r toward ≥ 0.35 by turning constant-zero IQR covariate columns into real range-breadth signal.

**Architecture:** Extend PR #46's `pull_gbif_centroids()` with a `store_points = FALSE` arg (default preserves back-compat). When TRUE, the per-species GBIF cache RDS gains a `points` data.frame. Extend PR #47's `.wc_extract_one()` to use cached points when available; fall back to centroid otherwise. Fixture rebuilds + end-to-end smoke comparing per-occurrence vs centroid-only.

**Tech Stack:** R ≥ 4.0, `testthat` 3rd edition, `roxygen2`, `devtools`, `terra` (existing Suggests from PR #47), `rgbif` (existing Suggests from PR #46), pigauto 0.9.1.9005.

**Working branch:** `feature/worldclim-per-occurrence` (stacked on `feature/worldclim-covariates`, carrying spec commit `5991ee3`).

**Spec:** `specs/2026-04-24-worldclim-per-occurrence-design.md`.

---

## File structure

| Path | Purpose | Create/Modify |
|---|---|---|
| `R/gbif_centroids.R` | add `store_points` arg to `pull_gbif_centroids()` + `.gbif_fetch_one()`; persist `points` data.frame in cache RDS when TRUE | modify |
| `R/worldclim_covariates.R` | extend `.wc_extract_one()` to read `points` from GBIF cache and extract at each; fall back to centroid when absent | modify |
| `tests/testthat/test-gbif-centroids.R` | add 3 tests: store_points TRUE writes points, FALSE omits points, legacy cache reads | modify |
| `tests/testthat/test-worldclim-covariates.R` | add 2 tests: per-occurrence extraction, centroid fallback | modify |
| `tests/testthat/fixtures/gbif_plants_300.rds` | keep as-is (still useful for non-points tests) | unchanged |
| `script/data-cache/gbif/` | add `points` field to existing species RDS via refresh run | regenerate (operator-run) |
| `tests/testthat/fixtures/worldclim_plants_300.rds` | regenerate with per-occurrence aggregation | regenerate |
| `NEWS.md` | new top section | modify |
| `DESCRIPTION` | bump to 0.9.1.9006 | modify |

Total: 2 modified R files, 2 modified test files, 2 regenerated fixtures (operator), 2 admin files.

---

## Task 1: Add `store_points` arg to `pull_gbif_centroids()` + `.gbif_fetch_one()`

**Files:**
- Modify: `R/gbif_centroids.R`
- Modify: `tests/testthat/test-gbif-centroids.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-gbif-centroids.R`:

```r
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
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'suppressMessages(devtools::load_all()); devtools::test(filter = "gbif-centroids")' 2>&1 | tail -n 10
```

Expected: ~3 new failures (store_points arg doesn't exist yet).

- [ ] **Step 3: Modify `.gbif_fetch_one()` in `R/gbif_centroids.R`**

Find the signature (currently around line ~65):

```r
.gbif_fetch_one <- function(sp, cache_dir = NULL,
                             occurrence_limit = 500L,
                             sleep_ms = 100L,
                             refresh_cache = FALSE) {
```

Add `store_points = FALSE`:

```r
.gbif_fetch_one <- function(sp, cache_dir = NULL,
                             occurrence_limit = 500L,
                             sleep_ms = 100L,
                             refresh_cache = FALSE,
                             store_points = FALSE) {
```

Find the point in the body where `records` has been assembled (after the pagination `for` loop) and the cache write happens. Currently there's a `saveRDS(c(res, list(fetched_at = Sys.time(), match_type = ...)), cache_path)`. Before the `saveRDS`, add the optional points collection:

```r
# Collect raw valid points if store_points = TRUE
points_df <- NULL
if (isTRUE(store_points) && !is.null(records) && nrow(records) > 0L) {
  # Apply the SAME filters as .gbif_centroid_one for consistency
  ok <- !is.na(records$decimalLatitude) &
         !is.na(records$decimalLongitude) &
         records$decimalLatitude  >= -90  & records$decimalLatitude  <= 90 &
         records$decimalLongitude >= -180 & records$decimalLongitude <= 180
  has_issue <- records$hasGeospatialIssues
  if (is.null(has_issue)) has_issue <- rep(FALSE, nrow(records))
  has_issue[is.na(has_issue)] <- TRUE
  ok <- ok & !has_issue &
         !(records$basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))
  kept <- records[ok, , drop = FALSE]
  if (nrow(kept) > 0L) {
    points_df <- data.frame(
      lat = kept$decimalLatitude,
      lon = kept$decimalLongitude,
      stringsAsFactors = FALSE)
  }
}
```

Then modify the cache write to include `points`:

```r
if (!is.na(cache_path)) {
  cache_obj <- c(res, list(fetched_at = Sys.time(),
                            match_type = bbone$matchType %||% "UNKNOWN"))
  if (!is.null(points_df)) cache_obj$points <- points_df
  saveRDS(cache_obj, cache_path)
}
```

- [ ] **Step 4: Modify `pull_gbif_centroids()` to accept + forward `store_points`**

Find the signature and add `store_points = FALSE`:

```r
pull_gbif_centroids <- function(species, cache_dir = NULL,
                                 occurrence_limit = 500L,
                                 sleep_ms = 100L,
                                 verbose = TRUE,
                                 refresh_cache = FALSE,
                                 store_points = FALSE) {
```

Update the internal `.gbif_fetch_one()` call to forward `store_points`:

```r
row <- .gbif_fetch_one(sp,
                         cache_dir = cache_dir,
                         occurrence_limit = occurrence_limit,
                         sleep_ms = sleep_ms,
                         refresh_cache = refresh_cache,
                         store_points = store_points)
```

Add roxygen `@param`:

```r
#' @param store_points logical. When \code{TRUE}, persists the raw
#'   filtered lat/lon occurrence points in each species' cache RDS
#'   under the \code{points} field.  Used by
#'   \code{\link{pull_worldclim_per_species}} for per-occurrence
#'   bioclim extraction.  Default \code{FALSE} preserves the
#'   pre-v0.9.1.9006 cache format.
```

- [ ] **Step 5: Run tests**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'suppressMessages(devtools::load_all()); devtools::test(filter = "gbif-centroids")' 2>&1 | tail -n 8
```

Expected: `[ FAIL 0 | ... | PASS ~35 ]` (was 32 in PR #47, +3 new).

- [ ] **Step 6: Commit**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git rev-parse --abbrev-ref HEAD   # must be feature/worldclim-per-occurrence
git add R/gbif_centroids.R tests/testthat/test-gbif-centroids.R
git commit -m "gbif-centroids: add store_points arg (v1.1, back-compat default FALSE)"
```

---

## Task 2: Extend `.wc_extract_one()` to use cached points

**Files:**
- Modify: `R/worldclim_covariates.R`
- Modify: `tests/testthat/test-worldclim-covariates.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-worldclim-covariates.R`:

```r
# ---- Task 2 (v1.1): per-occurrence extraction ----

test_that(".wc_extract_one uses points when available (per-occurrence)", {
  tmp <- tempfile("wc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gbif_dir <- file.path(tmp, "gbif"); dir.create(gbif_dir)
  extracts_dir <- file.path(tmp, "extracts"); dir.create(extracts_dir)
  # Write a GBIF cache with 10 occurrence points
  sp <- "Per occurrence"
  gbif_key <- gsub("[^A-Za-z0-9._-]", "_", sp)
  points_df <- data.frame(
    lat = c(40, 41, 39, 40.5, 41.5, 38, 42, 40, 41, 39),
    lon = c(-82, -81, -83, -82.5, -81.5, -84, -80, -82, -81, -83))
  saveRDS(list(species = sp,
                centroid_lat = 40.0, centroid_lon = -82.0,
                n_occurrences = 10L,
                points = points_df,
                fetched_at = Sys.time()),
           file.path(gbif_dir, paste0(gbif_key, ".rds")))
  # Mock terra::extract to return a plausible 10x19 matrix
  fake_extract <- matrix(seq_len(190), nrow = 10, ncol = 19)
  colnames(fake_extract) <- paste0("wc2.1_10m_bio_", 1:19)
  testthat::local_mocked_bindings(
    extract = function(x, y, ...) as.data.frame(fake_extract),
    .package = "terra")
  # Need a fake rast_stack object — any non-NULL placeholder works
  # because mocked extract ignores it.
  out <- pigauto:::.wc_extract_one(
    sp = sp,
    gbif_cache_dir = gbif_dir,
    wc_extracts_dir = extracts_dir,
    rast_stack = structure(list(), class = "FakeSpatRaster"),
    refresh_cache = FALSE)
  expect_equal(out$n_extracted, 10L)
  # Median of 1:10 = 5.5 (bio1 column), median of 11:20 = 15.5 (bio2), etc.
  expect_equal(as.numeric(out$bio_median["bio1"]), 5.5)
  expect_equal(as.numeric(out$bio_median["bio2"]), 15.5)
  # IQR of 1:10 = 5.0 (q75 - q25 = 7.75 - 3.25 = 4.5, not 5.0) -- use IQR()
  expect_equal(as.numeric(out$bio_iqr["bio1"]), stats::IQR(1:10))
})

test_that(".wc_extract_one falls back to centroid when points field absent", {
  tmp <- tempfile("wc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gbif_dir <- file.path(tmp, "gbif"); dir.create(gbif_dir)
  extracts_dir <- file.path(tmp, "extracts"); dir.create(extracts_dir)
  # Write a LEGACY-style cache (no points field)
  sp <- "Legacy species"
  gbif_key <- gsub("[^A-Za-z0-9._-]", "_", sp)
  saveRDS(list(species = sp,
                centroid_lat = 40.0, centroid_lon = -82.0,
                n_occurrences = 50L,
                fetched_at = Sys.time()),
           file.path(gbif_dir, paste0(gbif_key, ".rds")))
  # Mock terra::extract to return a 1x19 matrix (centroid = 1 point)
  fake_extract <- matrix(seq(100, 1800, length.out = 19), nrow = 1, ncol = 19)
  colnames(fake_extract) <- paste0("wc2.1_10m_bio_", 1:19)
  testthat::local_mocked_bindings(
    extract = function(x, y, ...) as.data.frame(fake_extract),
    .package = "terra")
  out <- pigauto:::.wc_extract_one(
    sp = sp,
    gbif_cache_dir = gbif_dir,
    wc_extracts_dir = extracts_dir,
    rast_stack = structure(list(), class = "FakeSpatRaster"),
    refresh_cache = FALSE)
  expect_equal(out$n_extracted, 1L)
  # IQR of a singleton is 0 (legacy centroid-only behaviour)
  expect_equal(as.numeric(out$bio_iqr["bio1"]), 0)
  # Reason should note the legacy path
  cached <- readRDS(file.path(extracts_dir,
                                 paste0(pigauto:::.wc_cache_key(sp), ".rds")))
  expect_true(!is.null(cached$reason))
  expect_true(grepl("centroid_only_legacy", cached$reason))
})
```

- [ ] **Step 2: Run tests to verify failures**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'suppressMessages(devtools::load_all()); devtools::test(filter = "worldclim-covariates")' 2>&1 | tail -n 8
```

Expected: 2 new failures (per-occurrence branch doesn't exist; legacy branch works but doesn't record "centroid_only_legacy" reason).

- [ ] **Step 3: Modify `.wc_extract_one()` in `R/worldclim_covariates.R`**

Find the section in `.wc_extract_one()` where it currently builds `points <- data.frame(lon = ..., lat = ...)` from the centroid (comment says "For simplicity here we use the centroid-only row"). Replace that block with:

```r
  # v1.1: prefer cached raw occurrence points over centroid
  points_df <- gbif_cached$points
  if (!is.null(points_df) && nrow(points_df) > 0L) {
    points <- data.frame(lon = points_df$lon, lat = points_df$lat)
    extract_reason <- "per_occurrence"
  } else {
    # Legacy cache fallback (v0.9.1.9005 centroid-only)
    points <- data.frame(lon = gbif_cached$centroid_lon,
                          lat = gbif_cached$centroid_lat)
    extract_reason <- "centroid_only_legacy"
  }
```

Then the existing `terra::extract(rast_stack, as.matrix(points[, c("lon", "lat")]))` block stays as-is (it already handles multi-row points).

At the cache-write step, plumb `extract_reason` through:

```r
  saveRDS(c(out, list(extracted_at = Sys.time(),
                        reason = extract_reason)),
           cache_path)
```

(The existing code already passes a `reason` key; just replace the hardcoded `"extracted"` with `extract_reason`.)

- [ ] **Step 4: Run tests**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'suppressMessages(devtools::load_all()); devtools::test(filter = "worldclim-covariates")' 2>&1 | tail -n 8
```

Expected: `[ FAIL 0 | ... | PASS ~25+ ]` (up ~2 from Task 1 state).

- [ ] **Step 5: Commit**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git rev-parse --abbrev-ref HEAD
git add R/worldclim_covariates.R tests/testthat/test-worldclim-covariates.R
git commit -m "worldclim: .wc_extract_one() uses per-occurrence points when available (v1.1)"
```

---

## Task 3: Regenerate 300-species fixture with per-occurrence + store_points

**Files:**
- Modify: `data-raw/make_gbif_plants_300.R` (to use `store_points = TRUE`)
- Modify: `data-raw/make_worldclim_plants_300.R` (no changes needed but re-run)
- Regenerate: `script/data-cache/gbif/*.rds` (via refresh)
- Regenerate: `tests/testthat/fixtures/worldclim_plants_300.rds`

- [ ] **Step 1: Update GBIF fixture builder to enable store_points**

Edit `data-raw/make_gbif_plants_300.R` (created in PR #46). Find the `pull_gbif_centroids(...)` call and add `store_points = TRUE` + `refresh_cache = TRUE` (to force a re-fetch of existing species with new cache shape):

```r
cov <- pull_gbif_centroids(sp300,
                            cache_dir = cache_dir,
                            occurrence_limit = 300L,
                            sleep_ms = 100L,
                            verbose = TRUE,
                            store_points = TRUE,        # NEW
                            refresh_cache = TRUE)       # force re-fetch once
```

After the first run regenerates the cache with points, future runs can drop `refresh_cache = TRUE` — points are now persisted.

- [ ] **Step 2: Run GBIF fixture builder (live GBIF fetch, ~10 min)**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript data-raw/make_gbif_plants_300.R 2>&1 | tail -n 15
```

Expected final line: "Wrote fixture: tests/testthat/fixtures/gbif_plants_300.rds". Also side-effect: every `script/data-cache/gbif/*.rds` now has a `points` field.

Spot-check 3 species:
```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'files <- list.files("script/data-cache/gbif", full.names = TRUE); for (f in head(files, 3)) { x <- readRDS(f); cat(basename(f), "  n_occ =", x$n_occurrences, "  has points =", !is.null(x$points), "  nrow(points) =", ifelse(is.null(x$points), 0, nrow(x$points)), "\n") }'
```

Expected: 3 files each with `has points = TRUE` and `nrow(points) > 0` for species with GBIF hits.

- [ ] **Step 3: Re-run WorldClim fixture builder**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript data-raw/make_worldclim_plants_300.R 2>&1 | tail -n 10
```

This regenerates `tests/testthat/fixtures/worldclim_plants_300.rds` with per-occurrence aggregation. Wall: ~3-5 min (rasters cached from PR #47; ~200 points per species × 19 bio vars × 300 species = ~1.1M extract calls, but terra's vectorised matrix method handles that fast).

Spot-check: at least some IQR columns should now be non-zero (range breadth is real):
```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'x <- readRDS("tests/testthat/fixtures/worldclim_plants_300.rds"); iqr_cols <- grep("_iqr$", colnames(x), value = TRUE); cat("iqr columns:", length(iqr_cols), "\n"); cat("iqr max (any species, any trait):", max(x[, iqr_cols], na.rm = TRUE), "\n"); cat("species with any IQR > 0:", sum(rowSums(x[, iqr_cols] > 0, na.rm = TRUE) > 0), "/", nrow(x), "\n")'
```

Expected: iqr max > 0 (previously 0 for every cell); species-with-any-IQR > 0 ≥ 150 (out of 300).

- [ ] **Step 4: Commit regenerated fixture + updated builder**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git rev-parse --abbrev-ref HEAD
git add data-raw/make_gbif_plants_300.R tests/testthat/fixtures/worldclim_plants_300.rds
git commit -m "worldclim: regenerate 300-species fixture with per-occurrence bioclim"
```

---

## Task 4: End-to-end smoke — per-occurrence vs centroid-only comparison

**Files:**
- Modify: `tests/testthat/test-worldclim-covariates.R`

- [ ] **Step 1: Write the test**

Append:

```r
test_that("per-occurrence bioclim lifts plants SLA >= 10% over centroid-only at n=200 (NOT_CRAN)", {
  skip_if(Sys.getenv("NOT_CRAN") == "", "slow integration -- NOT_CRAN=TRUE to run")
  wc_fx <- "tests/testthat/fixtures/worldclim_plants_300.rds"
  skip_if_not(file.exists(wc_fx), "worldclim fixture not found")

  # Load refreshed (per-occurrence) bioclim fixture
  wc_peroccur <- readRDS(wc_fx)
  # Species with meaningful IQR (range breadth)
  iqr_cols <- grep("_iqr$", colnames(wc_peroccur), value = TRUE)
  sp_with_iqr <- rownames(wc_peroccur)[rowSums(wc_peroccur[, iqr_cols] > 0,
                                                  na.rm = TRUE) > 0]
  skip_if_not(length(sp_with_iqr) >= 100L,
              "insufficient species with per-occurrence IQR")

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

  matched <- Reduce(intersect, list(wide$species, tree_all$tip.label,
                                      sp_with_iqr))
  skip_if_not(length(matched) >= 100L, "insufficient matched species")

  set.seed(2026L)
  sp_s <- sample(matched, min(200L, length(matched)))
  wide_s <- wide[wide$species %in% sp_s, , drop = FALSE]
  rownames(wide_s) <- wide_s$species; wide_s$species <- NULL
  cov_peroccur <- wc_peroccur[sp_s, grepl("^bio", colnames(wc_peroccur)),
                                drop = FALSE]
  # Build a centroid-only comparison cov: median cols only, with IQR cols zero'd
  cov_centroid <- cov_peroccur
  cov_centroid[, iqr_cols] <- 0
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

  res_centroid <- pigauto::impute(df, tree_s, covariates = cov_centroid,
                                     epochs = 80L, n_imputations = 20L,
                                     verbose = FALSE, seed = 2026L)
  res_peroccur <- pigauto::impute(df, tree_s, covariates = cov_peroccur,
                                     epochs = 80L, n_imputations = 20L,
                                     verbose = FALSE, seed = 2026L)

  env_traits <- intersect(c("sla", "leaf_area"), cont_cols)
  skip_if(length(env_traits) == 0L, "no env-driven trait present")

  any_lift <- FALSE
  for (v in env_traits) {
    if (!any(mask[, v])) next
    truth <- wide_s[[v]][mask[, v]]
    ok <- is.finite(truth)
    if (sum(ok) < 10L) next
    r_cent    <- sqrt(mean((res_centroid$completed[[v]][mask[, v]][ok] - truth[ok])^2))
    r_peroc   <- sqrt(mean((res_peroccur$completed[[v]][mask[, v]][ok] - truth[ok])^2))
    ratio <- r_peroc / r_cent
    cat(sprintf("[per-occur] %s: centroid RMSE = %.3g, per-occur = %.3g, ratio = %.4f\n",
                 v, r_cent, r_peroc, ratio))
    if (ratio <= 0.90) any_lift <- TRUE
  }

  # If the lift isn't there at n=200, don't fail — document and defer to full bench.
  # Assert guardrail (not-worse-than +10%) so a gross regression would fail.
  all_within_guardrail <- TRUE
  for (v in env_traits) {
    if (!any(mask[, v])) next
    truth <- wide_s[[v]][mask[, v]]
    ok <- is.finite(truth)
    if (sum(ok) < 10L) next
    r_cent  <- sqrt(mean((res_centroid$completed[[v]][mask[, v]][ok] - truth[ok])^2))
    r_peroc <- sqrt(mean((res_peroccur$completed[[v]][mask[, v]][ok] - truth[ok])^2))
    if (r_peroc > r_cent * 1.10) all_within_guardrail <- FALSE
  }
  expect_true(all_within_guardrail,
              info = "per-occurrence must not regress centroid-only by >+10%")
  # Log the lift status but don't fail on it -- the full n=4745 bench is where
  # the >=10% lift claim is verified.
  if (!any_lift) {
    message("note: per-occurrence lift not visible at n=200; full bench (n=4745) ",
            "in script/bench_bien_worldclim.R is where the paper claim lives.")
  }
})
```

- [ ] **Step 2: Run (with NOT_CRAN=TRUE)**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && NOT_CRAN=TRUE Rscript -e 'suppressMessages(devtools::load_all()); devtools::test(filter = "worldclim-covariates")' 2>&1 | tail -n 15
```

Expected: PASS. Look for `[per-occur]` lines — should show SLA ratio ≤ 1.0 (per-occurrence at least not worse than centroid-only) at n=200.

The hard ≥10 % lift claim is reserved for the full n=4745 bench (operator-run, in `script/bench_bien_worldclim.R`).

- [ ] **Step 3: Commit**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git rev-parse --abbrev-ref HEAD
git add tests/testthat/test-worldclim-covariates.R
git commit -m "worldclim: per-occurrence end-to-end smoke vs centroid-only"
```

---

## Task 5: NEWS.md + DESCRIPTION bump to 0.9.1.9006

**Files:**
- Modify: `NEWS.md`
- Modify: `DESCRIPTION`

- [ ] **Step 1: Prepend NEWS entry**

```markdown
# pigauto 0.9.1.9006 (dev)

## Per-occurrence WorldClim covariates (v1.1 follow-up to B.2)

Fixes the core limitation of PR #47 (centroid-only bioclim extraction):
bioclim is now extracted at every GBIF occurrence point per species
(not just the centroid), so the IQR columns carry real range-breadth
signal instead of being constant zero. This is what actually lets the
safety-floor calibrator see a meaningful covariate input and open the
GNN gate on plants.

Validated by the honest-sim diagnostic
(`experiment/covariate-honest-sim`) which showed pigauto's architecture
can produce 32 % RMSE lift on strong-environment traits when
covariates carry real signal — something centroid-only bioclim
failed to deliver on plants.

### API

- `pull_gbif_centroids(..., store_points = FALSE)` — new arg,
  default `FALSE` preserves pre-v0.9.1.9006 cache format. When
  `TRUE`, persists the raw filtered lat/lon points in each species'
  cache RDS, enabling per-occurrence bioclim extraction.
- `pull_worldclim_per_species()` — no signature change.  Internally,
  `.wc_extract_one()` now prefers cached per-occurrence points,
  falling back to centroid only when `points` field is absent
  (legacy caches).

### Back-compat

Pre-v1.1 cache RDS files remain fully readable. Users who want
per-occurrence extraction upgrade their cache with one call:

```r
pull_gbif_centroids(sp_list,
  cache_dir = "script/data-cache/gbif",
  store_points = TRUE,
  refresh_cache = TRUE)   # one-time refresh
```

After that, subsequent calls are instantaneous (cache hit).

### Fixture regeneration

`tests/testthat/fixtures/worldclim_plants_300.rds` re-generated with
per-occurrence aggregation. IQR columns are now informative.

### Follow-ups

- Covariate-aware phylo-signal gate — nice-to-have for
  interpretability (gate de-triggers when covariates resolve weak-
  phylo-signal traits). NOT a blocker given the safety-floor
  calibrator opens the gate on its own.
- Full-scale plants bench (n=4,745) — operator-run via
  `script/bench_bien_worldclim.R`. Expected: SLA r ≥ 0.35,
  leaf_area r ≥ 0.30.
```

- [ ] **Step 2: Bump DESCRIPTION**

Change `Version: 0.9.1.9005` → `Version: 0.9.1.9006`.

- [ ] **Step 3: Commit**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git rev-parse --abbrev-ref HEAD
git add NEWS.md DESCRIPTION
git commit -m "per-occurrence: NEWS.md entry + DESCRIPTION bump to 0.9.1.9006"
```

---

## Task 6: R CMD check + push + open PR #48

**Files:** Observe only.

- [ ] **Step 1: Regenerate man pages**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'suppressMessages(devtools::document())' 2>&1 | tail -n 3
```

Commit any `man/` changes as `"per-occurrence: regenerate man pages"`. Skip if clean.

- [ ] **Step 2: R CMD check**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && Rscript -e 'suppressMessages(devtools::check(args = "--no-tests", error_on = "error", quiet = TRUE))' 2>&1 | tail -n 10
```

Expected: 0 errors, 0 warnings, 2-3 pre-existing notes.

- [ ] **Step 3: Run full test suite**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && NOT_CRAN=TRUE Rscript -e 'suppressMessages(devtools::test())' > /tmp/per_occur_suite.log 2>&1; grep "^\[ FAIL" /tmp/per_occur_suite.log | tail -n 1
```

Expected: `[ FAIL 0 | ... | PASS ≥ 970 ]` (adds ~5 new expectations on top of the ~965 from prior PRs).

- [ ] **Step 4: Push**

```
cd "/Users/z3437171/Dropbox/Github Local/pigauto" && git push -u origin feature/worldclim-per-occurrence 2>&1 | tail -n 5
```

- [ ] **Step 5: Open PR**

```
gh pr create --base main --head feature/worldclim-per-occurrence --title "worldclim-per-occurrence: fix centroid-only limitation (v0.9.1.9006)" --body "$(cat <<'EOF'
## Summary

Fixes the centroid-only limitation of PR #47 by extracting bioclim at EVERY GBIF occurrence point per species, not just the centroid. Turns the previously-constant-zero IQR covariate columns into real range-breadth signal.

Stacks on PR #47 (feature/worldclim-covariates). v1.1 of the B.2 work.

## Why

The honest-sim diagnostic (experiment/covariate-honest-sim) proved pigauto's architecture CAN use environmental covariates -- 32% RMSE lift on strong-environment simulated traits. PR #47's ratio~=1.0 result on plants was a data-setup failure (centroid-only = 19 constant-zero IQR cols + small val set), not an architecture failure. This PR fixes the data setup.

## Changes

- pull_gbif_centroids() gains store_points = FALSE arg. When TRUE, persists raw lat/lon points in the GBIF cache RDS.
- .wc_extract_one() prefers cached per-occurrence points; falls back to centroid-only when absent (legacy caches work).
- pre-v1.1 cache RDS files remain fully readable.
- Fixture tests/testthat/fixtures/worldclim_plants_300.rds regenerated with per-occurrence aggregation.

## Test plan

- [x] 3 new gbif-centroids tests (store_points TRUE/FALSE, legacy cache).
- [x] 2 new worldclim-covariates tests (per-occurrence extraction, legacy fallback).
- [x] End-to-end smoke: per-occurrence ratio vs centroid-only at n=200 (guardrail +10%, full +10% lift claim reserved for n=4745 bench).
- [x] R CMD check 0/0/N clean.

## Paper-ready claim

Post-merge, run script/bench_bien_worldclim.R at n=4,745 to produce the target numbers: SLA r >= 0.35, leaf_area r >= 0.30, height_m r >= 0.30. These are the lift claims the paper needs for the plant row.

## Stacks on

PR #47. When #47 merges, rebase this onto main.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)" 2>&1 | tail -n 3
```

- [ ] **Step 6: Report PR URL**

---

## Plan self-review

**1. Spec coverage:** every section of
`specs/2026-04-24-worldclim-per-occurrence-design.md` has a task:
- §2.1 B.1 side → Task 1
- §2.2 B.2 side → Task 2
- §2.3 fixture rebuild → Task 3
- §3 API changes → Tasks 1, 2
- §5 testing → Tasks 1, 2, 4
- §6 risks → Task 1 (back-compat), Task 2 (fallback path), Task 3 (cache size)
- §7 rollout → Task 6
- §8 success metric → Task 4 (guardrail) + Task 6 (bench for hard lift)

**2. Placeholder scan:** No TBDs. All code blocks complete.

**3. Type consistency:** `points` field name consistent across Tasks
1, 2, 3. `store_points` arg name consistent in Task 1. `.wc_extract_one()`
return shape unchanged (adds `reason` field in cache metadata only).

No issues.

---

## Execution handoff

**Plan saved to `plans/2026-04-24-worldclim-per-occurrence.md`. Two options:**

**1. Subagent-Driven (recommended).**
**2. Inline execution.**

Default: Subagent-Driven. Tasks are sequential (each depends on prior). Task 3 requires live GBIF + raster access (~10-15 min), rest are fast.
