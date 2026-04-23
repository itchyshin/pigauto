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
