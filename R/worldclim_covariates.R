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

# Local fallback for %||% in case the file is used in isolation.
# pigauto defines %||% elsewhere too; this shim is harmless.
`%||%` <- function(a, b) if (is.null(a)) b else a
