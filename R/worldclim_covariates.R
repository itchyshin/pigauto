# R/worldclim_covariates.R
# Fetch and aggregate WorldClim v2.1 bioclim values at each species'
# GBIF occurrence points, producing per-species covariate rows for
# impute(covariates = ...).  Stacks on pull_gbif_centroids() (B.1).
# See specs/2026-04-24-worldclim-covariates-design.md.

# Sanitise a species name for use as a filesystem-safe cache key.
# Uses the SAME pattern as .gbif_cache_key (B.1) so GBIF + WorldClim
# caches key off the same string.  Any char not in [A-Za-z0-9_-]
# becomes `_`.
#
# @param sp character scalar, species name
# @return character scalar
# @noRd
.wc_cache_key <- function(sp) {
  if (!is.character(sp) || length(sp) != 1L || !nzchar(sp)) {
    stop(".wc_cache_key: species must be a non-empty character scalar")
  }
  gsub("[^A-Za-z0-9_-]", "_", sp)
}

# Aggregate a raster-extract data.frame (n_points x 19 bioclim cols)
# to per-species median and IQR.  Rows where ANY bio value is NA are
# treated as outside raster extent and dropped.  Returns NA-filled
# bio vectors when no rows survive.
#
# @param vals data.frame with columns "bio1" .. "bio19" (or whichever
#             subset was extracted); each row is one occurrence point
# @return list(bio_median, bio_iqr, n_extracted)
# @noRd
.wc_aggregate_one <- function(vals) {
  stopifnot(is.data.frame(vals))
  ncols <- ncol(vals)
  col_names <- colnames(vals)
  if (nrow(vals) == 0L) {
    na_vec <- rep(NA_real_, ncols); names(na_vec) <- col_names
    return(list(bio_median = na_vec, bio_iqr = na_vec, n_extracted = 0L))
  }
  # Drop rows where ANY bio value is NA (occurrence outside raster extent)
  ok <- stats::complete.cases(vals)
  kept <- vals[ok, , drop = FALSE]
  if (nrow(kept) == 0L) {
    na_vec <- rep(NA_real_, ncols); names(na_vec) <- col_names
    return(list(bio_median = na_vec, bio_iqr = na_vec, n_extracted = 0L))
  }
  med <- vapply(kept, stats::median, numeric(1L), na.rm = FALSE)
  iqr <- vapply(kept, function(x) stats::IQR(x, na.rm = FALSE),
                 numeric(1L))
  names(med) <- col_names
  names(iqr) <- col_names
  list(bio_median = med, bio_iqr = iqr, n_extracted = nrow(kept))
}
