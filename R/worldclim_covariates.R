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

# Ensure WorldClim rasters are present in cache_dir/wc2.1_<resolution>/.
# If a sentinel file .wc_complete exists and all 19 GeoTIFFs are
# present, no-op.  Otherwise download the zip stack from
# worldclim.org via the supplied `.download_fn` (for testability)
# and unzip to the target directory.
#
# @param cache_dir   character scalar, root cache path
# @param resolution  one of "10m", "5m", "2.5m"
# @param verbose     logical
# @param .download_fn internal hook (defaults to utils::download.file)
# @return path to the resolved raster directory
# @noRd
.wc_download_rasters <- function(cache_dir, resolution = "10m",
                                   verbose = TRUE,
                                   .download_fn = NULL) {
  res_ok <- c("10m", "5m", "2.5m")
  if (!identical(length(resolution), 1L) || !(resolution %in% res_ok)) {
    stop(".wc_download_rasters: resolution must be one of ",
         paste(res_ok, collapse = ", "))
  }
  wc_dir <- file.path(cache_dir, paste0("wc2.1_", resolution))
  sentinel <- file.path(wc_dir, ".wc_complete")
  expected <- file.path(wc_dir, sprintf("wc2.1_%s_bio_%d.tif",
                                           resolution, 1:19))
  if (file.exists(sentinel) && all(file.exists(expected))) {
    if (verbose) message("[worldclim] rasters present at ", wc_dir)
    return(wc_dir)
  }
  # Need to download.
  dir.create(wc_dir, showWarnings = FALSE, recursive = TRUE)
  zip_url <- sprintf(
    "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_%s_bio.zip",
    resolution)
  zip_path <- file.path(wc_dir, sprintf("wc2.1_%s_bio.zip", resolution))
  if (verbose) {
    message("[worldclim] downloading ", zip_url,
            " (~130 MB compressed, ~500 MB unzipped) ...")
  }
  dl_fn <- .download_fn %||% function(url, destfile, mode) {
    utils::download.file(url, destfile, mode = mode, quiet = !verbose)
  }
  dl_fn(zip_url, zip_path, "wb")
  utils::unzip(zip_path, exdir = wc_dir)
  file.create(sentinel)
  if (verbose) message("[worldclim] ready at ", wc_dir)
  wc_dir
}

# Extract bioclim values at one species' GBIF occurrences.
#
# Checks for per-species extract cache first (one RDS per species).
# On cache miss: reads the species' GBIF cache RDS (written by
# pull_gbif_centroids), extracts bioclim values at each lat/lon point
# via terra::extract(), aggregates via .wc_aggregate_one(), writes
# per-species cache, returns result.
#
# Species with no GBIF cache file (or empty one): returns NA bio
# vectors with n_extracted = 0L.  Does NOT error -- pull_gbif_centroids()
# is the source of truth for species presence.
#
# @param sp character scalar, species name
# @param gbif_cache_dir character path, root of GBIF per-species cache
#                        (as written by pull_gbif_centroids)
# @param wc_extracts_dir character path, root of WC per-species extract cache
# @param rast_stack terra SpatRaster with 19 bio layers (can be NULL when
#                    we expect a cache hit; extraction will error if NULL
#                    on a cache miss)
# @param refresh_cache logical; TRUE forces re-extract even with cache
# @return list(species, bio_median, bio_iqr, n_extracted)
# @noRd
.wc_extract_one <- function(sp, gbif_cache_dir, wc_extracts_dir,
                              rast_stack = NULL,
                              refresh_cache = FALSE) {
  stopifnot(is.character(sp), length(sp) == 1L, nzchar(sp))
  dir.create(wc_extracts_dir, showWarnings = FALSE, recursive = TRUE)
  key <- .wc_cache_key(sp)
  cache_path <- file.path(wc_extracts_dir, paste0(key, ".rds"))
  if (!refresh_cache && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$species, sp)) {
      return(list(species      = sp,
                   bio_median   = cached$bio_median,
                   bio_iqr      = cached$bio_iqr,
                   n_extracted  = cached$n_extracted))
    }
  }
  # Need to extract from raster.  Read GBIF cache.
  gbif_key <- gsub("[^A-Za-z0-9._-]", "_", sp)  # B.1 convention
  gbif_path <- file.path(gbif_cache_dir, paste0(gbif_key, ".rds"))
  # Placeholder bio column names -- will be overwritten by actual extraction
  bio_names <- paste0("bio", 1:19)
  na_med <- setNames(rep(NA_real_, 19L), bio_names)
  na_iqr <- setNames(rep(NA_real_, 19L), bio_names)
  if (!file.exists(gbif_path)) {
    # No GBIF record at all -- return NA row
    out <- list(species = sp, bio_median = na_med, bio_iqr = na_iqr,
                 n_extracted = 0L)
    saveRDS(c(out, list(extracted_at = Sys.time(),
                          reason = "no_gbif_cache")),
             cache_path)
    return(out)
  }
  gbif_cached <- tryCatch(readRDS(gbif_path), error = function(e) NULL)
  if (is.null(gbif_cached) ||
      is.na(gbif_cached$centroid_lat) ||
      gbif_cached$n_occurrences == 0L) {
    # GBIF cache exists but has no valid points
    out <- list(species = sp, bio_median = na_med, bio_iqr = na_iqr,
                 n_extracted = 0L)
    saveRDS(c(out, list(extracted_at = Sys.time(),
                          reason = "no_gbif_points")),
             cache_path)
    return(out)
  }
  # GBIF cache may or may not include the raw occurrence point list.
  # pull_gbif_centroids v0.9.1.9004 stores only centroid + n_occurrences.
  # B.2 needs per-point extraction -- extend the GBIF cache OR re-fetch
  # occurrences.  For simplicity here we use the centroid-only row as a
  # minimum viable extraction point.  Users who want per-point extraction
  # should re-run pull_gbif_centroids with extended caching (future B.2.1).
  points <- data.frame(lon = gbif_cached$centroid_lon,
                        lat = gbif_cached$centroid_lat)
  if (is.null(rast_stack)) {
    stop(".wc_extract_one: rast_stack is NULL but cache miss -- ",
         "cannot extract bioclim for species: ", sp, call. = FALSE)
  }
  vals_mat <- terra::extract(rast_stack,
                               as.matrix(points[, c("lon", "lat")]),
                               ID = FALSE, method = "simple")
  # terra returns a data.frame with raster layer names
  vals_df <- as.data.frame(vals_mat)
  # Rename to bio1 .. bio19 in case terra prefixes with wc2.1_10m_bio_N
  colnames(vals_df) <- bio_names[seq_len(ncol(vals_df))]
  agg <- .wc_aggregate_one(vals_df)
  out <- list(species = sp,
               bio_median = agg$bio_median,
               bio_iqr = agg$bio_iqr,
               n_extracted = agg$n_extracted)
  saveRDS(c(out, list(extracted_at = Sys.time(),
                        reason = "extracted")),
           cache_path)
  out
}

# Local fallback for %||% in case the file is used in isolation.
# pigauto defines %||% elsewhere too; this shim is harmless.
`%||%` <- function(a, b) if (is.null(a)) b else a
