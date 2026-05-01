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
  # hasGeospatialIssues may be absent from the GBIF response (column NULL) —
  # the v3 API does not always return this field.  Treat absent as FALSE
  # (pass through; we already filter on coord range and hasCoordinate=TRUE).
  # For a present column, treat NA as TRUE (conservative).
  has_issue <- records$hasGeospatialIssues
  if (is.null(has_issue)) {
    has_issue <- rep(FALSE, nrow(records))
  } else {
    has_issue[is.na(has_issue)] <- TRUE
  }
  # basisOfRecord may also be absent; treat absent as not-excluded.
  bor <- records$basisOfRecord
  if (is.null(bor)) bor <- rep("UNKNOWN", nrow(records))
  ok <- !is.na(records$decimalLatitude) &
         !is.na(records$decimalLongitude) &
         records$decimalLatitude  >= -90  & records$decimalLatitude  <= 90 &
         records$decimalLongitude >= -180 & records$decimalLongitude <= 180 &
         !has_issue &
         !(bor %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))
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

#' Fetch species range-centroid covariates from GBIF
#'
#' Pulls occurrence records from the Global Biodiversity Information
#' Facility (GBIF) for each species in \code{species}, aggregates them
#' to a median lat / lon centroid, and returns a data.frame ready to
#' be passed to \code{\link{impute}} via its \code{covariates} argument.
#'
#' Caching is strongly recommended — GBIF rate-limits anonymous calls
#' and the per-species fetch is expensive. With \code{cache_dir} set,
#' each species gets one RDS file; subsequent calls skip the API.
#'
#' For each species: resolve taxon via
#' \code{\link[rgbif]{name_backbone}}, fetch up to
#' \code{occurrence_limit} records via \code{\link[rgbif]{occ_search}}
#' (paginated at 300 per GBIF call), filter out records with
#' \code{hasGeospatialIssues = TRUE} and \code{basisOfRecord} in
#' \code{c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")}, drop
#' out-of-range coordinates, then take the median latitude and
#' longitude as the species centroid.
#'
#' Species with zero post-filter records receive \code{NA} centroids;
#' their rows are still included in the returned data.frame so it
#' aligns with the input species list.
#'
#' @param species           character vector of species names.
#' @param cache_dir         directory to cache per-species RDS files.
#'   \code{NULL} (default) disables caching — not recommended for
#'   production use.
#' @param occurrence_limit  integer, maximum number of occurrences to
#'   fetch per species (default 500; paginated if > 300).
#' @param sleep_ms          integer, polite delay between API calls in
#'   milliseconds (default 100).
#' @param verbose           logical, print progress every 50 species.
#' @param refresh_cache     logical, force re-fetch even when cache
#'   exists.
#'
#' @return A data.frame with columns \code{species}, \code{centroid_lat},
#'   \code{centroid_lon}, \code{n_occurrences}. Rownames are set to
#'   \code{species}. NA centroids are returned for species with no
#'   GBIF hits; the caller should decide how to handle them (typical:
#'   drop or impute).
#'
#' @seealso \code{\link{impute}} (pass the return value as
#'   \code{covariates}).
#'
#' @examples
#' \dontrun{
#' # Plants ecology example: pull centroids for a species list.
#' sp <- c("Quercus alba", "Pinus taeda", "Acer saccharum")
#' cov <- pull_gbif_centroids(sp, cache_dir = "script/data-cache/gbif")
#' # Use as covariates (drop the bookkeeping cols)
#' cov_num <- cov[, c("centroid_lat", "centroid_lon"), drop = FALSE]
#' # Then: impute(traits, tree, covariates = cov_num)
#' }
#'
#' @export
pull_gbif_centroids <- function(species, cache_dir = NULL,
                                 occurrence_limit = 500L,
                                 sleep_ms = 100L,
                                 verbose = TRUE,
                                 refresh_cache = FALSE) {
  stopifnot(is.character(species), length(species) >= 1L)
  n_sp <- length(species)
  rows <- vector("list", n_sp)
  n_success <- 0L
  for (i in seq_along(species)) {
    sp <- species[[i]]
    row <- .gbif_fetch_one(sp,
                             cache_dir = cache_dir,
                             occurrence_limit = occurrence_limit,
                             sleep_ms = sleep_ms,
                             refresh_cache = refresh_cache)
    rows[[i]] <- row
    if (is.finite(row$centroid_lat)) n_success <- n_success + 1L
    if (verbose && (i %% 50L == 0L || i == n_sp)) {
      cat(sprintf("[gbif] %d/%d species fetched, %d with valid centroids\n",
                   i, n_sp, n_success))
    }
  }
  out <- data.frame(
    species       = vapply(rows, function(r) r$species, character(1L)),
    centroid_lat  = vapply(rows, function(r) r$centroid_lat, numeric(1L)),
    centroid_lon  = vapply(rows, function(r) r$centroid_lon, numeric(1L)),
    n_occurrences = vapply(rows, function(r) r$n_occurrences, integer(1L)),
    stringsAsFactors = FALSE,
    row.names = species)
  out
}

# In-file %||% shim (pigauto defines one elsewhere, but making this file
# self-contained helps testthat mock_bindings for rgbif).
`%||%` <- function(a, b) if (is.null(a)) b else a
