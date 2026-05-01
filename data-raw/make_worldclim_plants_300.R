# data-raw/make_worldclim_plants_300.R
# One-shot fixture builder.  Requires internet + ~30 min wall.
# Re-run only when the fixture species set changes materially.

library(devtools)
devtools::load_all(".", quiet = TRUE)

gbif_fx <- "tests/testthat/fixtures/gbif_plants_300.rds"
if (!file.exists(gbif_fx))
  stop("Run data-raw/make_gbif_plants_300.R first (B.1 fixture).")

cov_gbif <- readRDS(gbif_fx)
sp300 <- rownames(cov_gbif)

gbif_cache_dir <- "script/data-cache/gbif"
wc_cache_dir   <- "script/data-cache/worldclim"
dir.create(wc_cache_dir, showWarnings = FALSE, recursive = TRUE)

wc_df <- pull_worldclim_per_species(
  species             = sp300,
  gbif_cache_dir      = gbif_cache_dir,
  worldclim_cache_dir = wc_cache_dir,
  resolution          = "10m",
  verbose             = TRUE)

out_path <- "tests/testthat/fixtures/worldclim_plants_300.rds"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(wc_df, out_path, compress = "xz")
cat("Wrote fixture:", out_path,
    " -- rows =", nrow(wc_df),
    " with valid bioclim =", sum(wc_df$n_extracted > 0), "\n")
