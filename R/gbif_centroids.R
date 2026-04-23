# R/gbif_centroids.R
# Fetch and aggregate GBIF species-occurrence records to per-species
# centroid lat/lon covariates for pigauto.
# See specs/2026-04-23-gbif-centroids-design.md.

# Sanitise a species name for use as a filesystem-safe cache key.
# Any character that's not alphanumeric, `_`, `.`, or `-` is replaced
# with `_`. Guaranteed non-empty; guaranteed to not contain `/`.
#
# @param sp character scalar, species name
# @return character scalar
# @noRd
.gbif_cache_key <- function(sp) {
  stopifnot(is.character(sp), length(sp) == 1L, nzchar(sp))
  key <- gsub("[^A-Za-z0-9_-]", "_", sp)
  if (!nzchar(key)) key <- "unnamed"
  key
}

# Aggregate a GBIF occurrence data.frame to a single centroid.
#
# Input must have at least these columns:
#   decimalLatitude, decimalLongitude, hasGeospatialIssues, basisOfRecord
#
# Applies filters (dropping fossil specimens, living specimens, records
# flagged with geospatial issues, NA coords, out-of-range coords), then
# returns list(centroid_lat, centroid_lon, n_occurrences).
#
# @noRd
.gbif_centroid_one <- function(records) {
  if (nrow(records) == 0L) {
    return(list(centroid_lat = NA_real_,
                 centroid_lon = NA_real_,
                 n_occurrences = 0L))
  }
  ok <- !is.na(records$decimalLatitude) &
         !is.na(records$decimalLongitude) &
         records$decimalLatitude  >= -90  & records$decimalLatitude  <= 90 &
         records$decimalLongitude >= -180 & records$decimalLongitude <= 180 &
         !isTRUE(records$hasGeospatialIssues) &
         !(records$basisOfRecord %in% c("FOSSIL_SPECIMEN",
                                           "LIVING_SPECIMEN"))
  # hasGeospatialIssues may be a logical vector of NA's; treat NA as TRUE (conservative)
  has_issue <- records$hasGeospatialIssues
  has_issue[is.na(has_issue)] <- TRUE
  ok <- ok & !has_issue
  kept <- records[ok, , drop = FALSE]
  if (nrow(kept) == 0L) {
    return(list(centroid_lat = NA_real_,
                 centroid_lon = NA_real_,
                 n_occurrences = 0L))
  }
  list(centroid_lat = stats::median(kept$decimalLatitude),
        centroid_lon = stats::median(kept$decimalLongitude),
        n_occurrences = nrow(kept))
}

# Fetch (and cache) centroids for one species.
#
# Cache check first: if an RDS exists at {cache_dir}/{cache_key}.rds,
# load it unless refresh_cache = TRUE.
#
# Otherwise: rgbif::name_backbone() -> rgbif::occ_search() -> filter -> aggregate.
# Cache the result (including n_occ = 0 cases so we don't re-hit the API
# for a species with no records).
#
# On network error: retry 3x with exponential backoff; if still failing,
# cache as (NA, NA, 0, fetched_at = NA) so subsequent calls short-circuit
# until refresh_cache = TRUE.
#
# @noRd
.gbif_fetch_one <- function(sp, cache_dir = NULL,
                             occurrence_limit = 500L,
                             sleep_ms = 100L,
                             refresh_cache = FALSE) {
  stopifnot(is.character(sp), length(sp) == 1L, nzchar(sp))
  # Cache lookup
  cache_path <- if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    file.path(cache_dir, paste0(.gbif_cache_key(sp), ".rds"))
  } else {
    NA_character_
  }
  if (!refresh_cache && !is.na(cache_path) && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$species, sp)) {
      return(list(species = sp,
                   centroid_lat = cached$centroid_lat,
                   centroid_lon = cached$centroid_lon,
                   n_occurrences = cached$n_occurrences))
    }
  }
  # Require rgbif for live fetch
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("pull_gbif_centroids() requires the 'rgbif' package: ",
         "install.packages('rgbif')", call. = FALSE)
  }
  # Resolve taxon
  bbone <- tryCatch(rgbif::name_backbone(name = sp),
                     error = function(e) NULL)
  if (is.null(bbone) || is.null(bbone$usageKey)) {
    res <- list(species = sp, centroid_lat = NA_real_, centroid_lon = NA_real_,
                 n_occurrences = 0L)
    if (!is.na(cache_path)) {
      saveRDS(c(res, list(fetched_at = Sys.time(), match_type = "NO_MATCH")),
              cache_path)
    }
    return(res)
  }
  # Fetch occurrences (with pagination to support occurrence_limit > 300)
  per_call <- 300L
  pages <- ceiling(occurrence_limit / per_call)
  records <- NULL
  for (p in seq_len(pages)) {
    offset <- (p - 1L) * per_call
    limit  <- min(per_call, occurrence_limit - offset)
    if (limit <= 0L) break
    ret <- .gbif_fetch_page(bbone$usageKey, offset, limit)
    if (is.null(ret) || nrow(ret) == 0L) break
    records <- rbind(records, ret)
    if (sleep_ms > 0L) Sys.sleep(sleep_ms / 1000)
  }
  # Aggregate
  centroid <- .gbif_centroid_one(records %||% data.frame(
    decimalLatitude = numeric(0), decimalLongitude = numeric(0),
    hasGeospatialIssues = logical(0), basisOfRecord = character(0)))
  res <- list(species = sp,
               centroid_lat = centroid$centroid_lat,
               centroid_lon = centroid$centroid_lon,
               n_occurrences = centroid$n_occurrences)
  if (!is.na(cache_path)) {
    saveRDS(c(res, list(fetched_at = Sys.time(),
                          match_type = bbone$matchType %||% "UNKNOWN")),
            cache_path)
  }
  res
}

# Inner page-fetcher, retriable.  Separated so tests can mock.
# @noRd
.gbif_fetch_page <- function(taxon_key, offset, limit) {
  tries <- 0L
  while (tries < 3L) {
    ret <- tryCatch(
      rgbif::occ_search(taxonKey = taxon_key,
                         hasCoordinate = TRUE,
                         limit = limit, start = offset),
      error = function(e) {
        Sys.sleep(2 ^ tries)
        NULL
      })
    if (!is.null(ret) && !is.null(ret$data)) return(ret$data)
    tries <- tries + 1L
  }
  NULL
}

# In-file %||% shim (pigauto defines one elsewhere, but making this file
# self-contained helps testthat mock_bindings for rgbif).
`%||%` <- function(a, b) if (is.null(a)) b else a
