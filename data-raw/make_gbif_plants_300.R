# data-raw/make_gbif_plants_300.R
# One-shot fixture builder for the end-to-end plants centroid test.
# Re-run only when the plants species set or GBIF data changes materially.

library(devtools)
devtools::load_all(".", quiet = TRUE)

bien_path <- "script/data-cache/bien_trait_means.rds"
if (!file.exists(bien_path))
  stop("Run script/bench_bien.R first to build the BIEN cache.")

trait_means <- readRDS(bien_path)
all_species <- Reduce(union, lapply(trait_means,
  function(d) if (!is.null(d)) d$species else character(0)))
set.seed(2026L)
sp300 <- sample(all_species, 300L)

cache_dir <- "script/data-cache/gbif"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cov <- pull_gbif_centroids(sp300,
                            cache_dir = cache_dir,
                            occurrence_limit = 300L,
                            sleep_ms = 100L,
                            verbose = TRUE,
                            store_points = TRUE,        # NEW (v1.1)
                            refresh_cache = TRUE)       # one-time refetch to add points to existing cache

out_path <- "tests/testthat/fixtures/gbif_plants_300.rds"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(cov, out_path, compress = "xz")
cat("Wrote fixture:", out_path,
    " -- rows =", nrow(cov),
    " with valid centroid =", sum(is.finite(cov$centroid_lat)), "\n")
