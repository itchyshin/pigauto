# GBIF range-centroid covariates — design spec

**Status:** approved (brainstorming 2026-04-23).
**Branch:** `feature/gbif-centroids` (stacked on `feature/safety-floor-mean-gate`).
**Plan:** `plans/2026-04-23-gbif-centroids.md` (to be written next).
**Stacks on:** `specs/2026-04-23-safety-floor-mean-gate-design.md` (PR #43).

---

## 1. Goal

Ship `pull_gbif_centroids(species, cache_dir)` — a new exported helper
that fetches species occurrence records from GBIF and returns a
2-column numeric data.frame (`centroid_lat`, `centroid_lon`) indexed
by species name, ready to pass into `impute(..., covariates = ...)`.
Validates the covariate pipeline end-to-end on real biogeographic
data before we commit to the heavier B.2 bioclim/SoilGrids raster
extraction.

## 2. Background

The plants bench at commit `16ce036` showed r = −0.02 to 0.43 across
five BIEN traits. PR #43 (safety floor) keeps pigauto from being
worse than the grand mean on these, but does not *lift* them toward
meaningful r without a covariate signal the GNN can exploit. Plant
traits (SLA, leaf area, seed mass) are fundamentally
environment-driven — bioclim and soil variables predict them better
than phylogeny.

B.1 is the proof-of-pipe: even the coarsest geographic covariate
(range centroid: 2 numbers per species) should lift SLA and leaf
area on plants noticeably. If it does, B.2 (full bioclim + soil)
has a clear path. If it doesn't, we need a different architectural
approach (foundation model, trait-similarity attention) before
committing to raster infrastructure.

## 3. Non-goals

- Bioclim variable extraction (WorldClim v2, CHELSA) — that is B.2.
- SoilGrids extraction — also B.2.
- Facial-recognition / images / non-GBIF species occurrence sources
  (e.g. BIEN's own occurrence endpoints) — out of scope.
- Changing pigauto internals. B.1 adds ONE exported data-fetching
  helper plus caching infrastructure; pigauto's existing
  `covariates` input already handles numeric data.frames.
- Any phylogenetic / taxonomic tree manipulation. The helper is
  pigauto-agnostic; it accepts a character vector of species names
  and returns a data.frame.

## 4. Architecture

### 4.1 Public API

```r
pull_gbif_centroids(
  species,                                 # character vector of binomials
  cache_dir       = NULL,                  # optional cache directory
  occurrence_limit = 500L,                 # max occurrences per species
  sleep_ms        = 100L,                  # polite delay between calls
  verbose         = TRUE,                  # print progress
  refresh_cache   = FALSE                  # force re-fetch even if cached
)
# → data.frame with columns `species` (character),
#   `centroid_lat` (numeric), `centroid_lon` (numeric),
#   `n_occurrences` (integer). Rownames = species.
#   Species with no GBIF hits get NA for lat/lon and 0 for n_occurrences.
```

### 4.2 GBIF request semantics

Per species:

1. `rgbif::name_backbone(name = sp)` to resolve to a GBIF `taxonKey`.
   Handles synonyms, misspellings, and taxonomic rank ambiguity.
2. `rgbif::occ_search(taxonKey = key, hasCoordinate = TRUE,
   limit = min(300, occurrence_limit))` — GBIF caps each call at 300.
   For occurrence_limit > 300, paginate via `offset`.
3. Filter out:
   - Records with `hasGeospatialIssues = TRUE`
   - Records with `basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")`
     (keep PRESERVED_SPECIMEN, HUMAN_OBSERVATION, MACHINE_OBSERVATION,
     MATERIAL_SAMPLE).
   - Records outside `[-90, 90]` lat or `[-180, 180]` lon (defensive).
4. Aggregate: median latitude, median longitude (centroids more
   robust to outliers than means — common-garden specimens in wrong
   latitude).

### 4.3 Caching

If `cache_dir` is set: one RDS per species at
`{cache_dir}/{underscore_safe_binomial}.rds` containing a list
`list(species = sp, n_occ = n, centroid_lat = lat, centroid_lon = lon,
fetched_at = Sys.time())`. Subsequent calls read from cache.

`refresh_cache = TRUE` forces re-fetch (for updating old caches).

Species with 0 occurrences are cached as
`list(centroid_lat = NA, centroid_lon = NA, n_occ = 0L)` — still cached,
so we don't re-hit the API for a known-empty species.

Cache file naming: species names may contain spaces and unicode. Use
`gsub("[^A-Za-z0-9._-]", "_", sp)` to produce a filesystem-safe key.
Collisions avoided by including the full raw species name inside the
cached RDS — on cache hit we verify `cache$species == sp` before
trusting the centroid.

### 4.4 Error handling

- Network error / GBIF timeout: retry with exponential backoff (3
  attempts, base 1s); after 3 fails, emit a warning and cache as
  `n_occ = 0L` with `fetched_at = NA_real_` so subsequent calls
  don't re-try until `refresh_cache = TRUE`.
- Unresolved `name_backbone` match (no taxonKey): cache as
  `n_occ = 0L` with a `match_type = "NO_MATCH"` note.
- Ambiguous backbone match (multiple taxonKeys with same confidence):
  pick the first and emit a warning with both candidates.

### 4.5 Progress reporting

When `verbose = TRUE`, print one line per 50 species:

```
[gbif] 50/300 species fetched, 41 with valid centroids (9 no hits)
```

### 4.6 Downstream integration

No pigauto code change. User does:

```r
cov <- pull_gbif_centroids(rownames(df), cache_dir = "script/data-cache/gbif")
# drop the bookkeeping col, keep only numeric covariates for pigauto
cov <- cov[, c("centroid_lat", "centroid_lon")]
res <- impute(df, tree, covariates = cov)
```

pigauto's existing `covariates` path in `preprocess_traits()` /
`fit_pigauto()` handles this as it does any numeric covariate matrix.

### 4.7 File layout

```
R/gbif_centroids.R           # new, ~180 lines
  ├── pull_gbif_centroids()  # public, @export
  ├── .gbif_cache_key()      # internal helper (name sanitisation)
  ├── .gbif_fetch_one()      # internal per-species fetch
  └── .gbif_centroid_one()   # internal aggregation from occurrences
```

Internal helpers prefixed with `.` to keep them off the user-facing
surface area.

### 4.8 Dependencies

Add `rgbif` to **Suggests** (not Imports — covariate fetching is
optional). Fail fast if not installed:

```r
if (!requireNamespace("rgbif", quietly = TRUE)) {
  stop("pull_gbif_centroids() requires the 'rgbif' package: ",
       "install.packages('rgbif')")
}
```

## 5. Testing

### 5.1 Offline unit (`tests/testthat/test-gbif-centroids.R`)

1. **Cache key sanitisation.** `.gbif_cache_key("Quercus alba L.")`
   returns a filesystem-safe string.
2. **Centroid aggregation.** Given a synthetic occurrence data.frame
   with known lat/lon values, `.gbif_centroid_one()` returns the
   median.
3. **Record filtering.** Synthetic occurrences with
   `hasGeospatialIssues = TRUE` + `basisOfRecord = "FOSSIL_SPECIMEN"`
   + out-of-range coordinates are dropped; valid records are kept.
4. **Cache read/write round-trip.** Write a fake cached RDS, read it
   back via `pull_gbif_centroids()` with `cache_dir` pointing at it,
   verify no network call is made (mock `rgbif::occ_search`).
5. **Graceful degradation when rgbif absent.** `pull_gbif_centroids()`
   stops with a clear install hint.

Mocks: use `testthat::local_mocked_bindings(occ_search = function(...)
list(data = fake_df))` to avoid any real network call in CI.

### 5.2 Online integration (skipped by default, `skip_if_offline()`)

1. **10-species live fetch.** For a hand-picked list of 10
   taxonomically-diverse, occurrence-rich species (e.g. `Quercus_alba`,
   `Pinus_taeda`, `Zea_mays`, ...), fetch centroids; verify every
   species has `centroid_lat` in `[-90, 90]` and `centroid_lon` in
   `[-180, 180]`, and that `n_occurrences >= 10`.
2. **Cache persistence.** Fetch once; call again; verify zero
   network activity on second call (by asserting wall time < 1 s).

### 5.3 End-to-end smoke (`tests/testthat/test-gbif-centroids.R`)

**Skipped if** the BIEN cache is absent OR online fixture RDS is
absent.

**Fixture approach** (avoids hitting GBIF in CI): ship
`tests/testthat/fixtures/gbif_plants_300.rds` — a pre-fetched
centroid table for 300 BIEN species, built once by
`data-raw/make_gbif_plants_300.R`.

Test:

1. Load plants 300 subset (BIEN cache).
2. Load the 300-species GBIF centroid fixture.
3. Fit pigauto with + without covariates.
4. Assert: on at least ONE of `sla` / `leaf_area` /
   `seed_mass`, the with-cov RMSE is **≥ 5 % lower** than
   without-cov RMSE. This is the proof-of-pipe claim.

## 6. Risks and mitigations

### 6.1 GBIF rate limits

GBIF enforces ~1 rps soft limit for anonymous calls. At 300 species
with `sleep_ms = 100` the floor wall is 300 × 0.4s ≈ 2 min. For
larger species pools (1000s), add documentation recommending the
`GBIF_USER` + `GBIF_PWD` env-var path (rgbif uses it for download
endpoints) and letting users fetch overnight.

### 6.2 Cache staleness

GBIF data changes (new records, taxon reclassification).
Document that `refresh_cache = TRUE` should be used periodically
(quarterly is reasonable for comparative-biology analyses).
Cache RDS stores `fetched_at` for inspection.

### 6.3 Tree alignment

Species in the user's dataset may have GBIF synonyms that match a
different `taxonKey`. We accept this — `name_backbone()` resolves
synonyms. Emit a warning if `match_type != "EXACT"` so users know
a name substitution happened.

### 6.4 Centroids can be meaningless

A cosmopolitan species (e.g. *Linum usitatissimum* cultivated
worldwide) has a centroid in the middle of Europe that says
nothing about its native range. Document this as a known limitation
in the roxygen; recommend B.2 (bioclim variance + centroid) for
cosmopolitan taxa.

### 6.5 Non-plant benches

The helper is taxon-agnostic. It works fine for birds, mammals,
fish, amphibians. But for marine taxa (fish), the lat/lon centroid
may be less informative than a *depth* covariate. We accept this
— non-plant users can use it as a cheap baseline or ignore it.

## 7. Rollout

Single PR on `feature/gbif-centroids`. Stacks on #43's work.
After #43 merges:

1. Rebase `feature/gbif-centroids` onto `main`.
2. Regenerate man pages + run R CMD check.
3. Run the online integration test once locally to prime the
   fixture RDS.
4. Merge.

Version bump: `0.9.1.9002` → `0.9.1.9004` (reserves `.9003` for
the phylo-signal gate PR).

## 8. Success metric

After merge:

- `pull_gbif_centroids()` is callable from a clean R session on any
  character-vector species list; cached results persist across
  sessions; online test passes; offline tests run in CI.
- Integration test on plants 300-subset shows at least one
  environment-driven trait (sla / leaf_area) improving by ≥ 5 %
  RMSE when centroids are supplied — the paper claim "pigauto
  accepts arbitrary covariates; here's a proof-of-pipe with GBIF
  centroids".
- Package install size increases by < 500 KB (the fixture RDS +
  new R file; no new Imports).

## 9. Open follow-ups (out of scope)

- B.2 — full WorldClim + SoilGrids extraction at each occurrence,
  aggregated per species. Separate spec after B.1 lands clean.
- Trait-specific covariate selection (e.g. only use centroids for
  traits with λ < threshold, per the phylo-signal gate PR).
- Multi-observation path: use per-observation lat/lon instead of
  species median. Useful for ctmax_sim-style acclimation
  studies. Separate spec.
