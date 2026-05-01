# WorldClim bioclim covariates — design spec

**Status:** approved (brainstorming 2026-04-24).
**Branch:** `feature/worldclim-covariates` (stacked on `feature/gbif-centroids`).
**Plan:** `plans/2026-04-24-worldclim-covariates.md` (to be written next).
**Stacks on:** PR #46 (`feature/gbif-centroids`, v0.9.1.9004).
**Target version:** v0.9.1.9005.

---

## 1. Goal

Ship `pull_worldclim_per_species(species, gbif_cache_dir, worldclim_cache_dir)` —
a new exported helper that extends the B.1 GBIF centroid pipeline by
extracting **19 WorldClim v2.1 bioclim variables** at each species'
GBIF occurrence points, aggregating per species (median + IQR) to 38
numeric covariates, and returning a data.frame ready for
`impute(..., covariates = ...)`.

This is the B.2 follow-up to B.1 (GBIF centroids alone). Centroids
don't carry enough climate signal to resolve environment-driven plant
traits. Full bioclim variables are expected to lift plants SLA and
leaf_area RMSE from the grand-mean floor to **r ≥ 0.35** on a BIEN
1000-species subset — the proof that pigauto can crack weak-phylo-
signal data when given real climate covariates.

## 2. Background

The plants benchmark sequence:

1. **Pre-safety-floor (v0.9.1.9000)**: pigauto 15-101 % worse than
   grand mean on 4 of 5 BIEN traits. Honest boundary case.
2. **After PR #43 (safety floor)**: never worse than mean by
   construction. But the 3 weak-signal traits (height_m, sla,
   leaf_area) sit AT the grand-mean floor — no lift.
3. **After PR #45 (phylo-signal gate)**: those 3 traits are now
   *explicitly* flagged as λ ≈ 0 and routed to pure mean. The "why"
   becomes diagnostic-based.
4. **After PR #46 (GBIF centroids)**: covariate pipeline validated.
   Centroids themselves don't help at n=200 (ratio 1.01–1.09);
   proof-of-pipe only.
5. **This spec (B.2)**: full bioclim per-species. Target: lift
   SLA + leaf_area + height_m **to r ≥ 0.35 via covariates**.

The scientific premise is well-established: plant traits correlate
strongly with climate (Leaf Economics Spectrum, Wright et al. 2004).
Centroid lat/lon discards this information; bioclim preserves it.

## 3. Non-goals

This spec does **not** cover:

- SoilGrids variables (soil pH, texture, SOC). Deferred to a
  potential B.3 spec if B.2 lifts plants but still leaves headroom.
- CHELSA v2 as an alternative to WorldClim. WorldClim is the
  canonical ecology standard; switching raster sources can be a
  separate user choice later.
- Spatial resolution finer than 10-arc-minute. The 10-min raster
  (~18 km) is adequate for range-level aggregation; 5-min / 2.5-min
  would 4× or 16× the disk footprint without adding signal at this
  scale.
- Changes to pigauto's internals. B.2 only adds one exported helper
  plus a raster cache; pigauto's existing `covariates` path already
  handles numeric data.frames.
- Per-observation covariates (multi-obs mode). This helper is
  species-level aggregation only.
- Runtime WorldClim-raster distribution. The helper downloads from
  `worldclim.org` on first call; no raster ships with the package.

## 4. Architecture

### 4.1 Public API

```r
pull_worldclim_per_species(
  species,                                 # character vector of binomials
  gbif_cache_dir,                          # REQUIRED: path to B.1 GBIF cache
  worldclim_cache_dir,                     # REQUIRED: path to WorldClim raster cache
  resolution       = "10m",                # "10m" (default), "5m", or "2.5m"
  variables        = paste0("bio", 1:19),  # which bioclim vars (default all 19)
  aggregation      = c("median", "iqr"),   # median + IQR by default
  verbose          = TRUE,
  refresh_cache    = FALSE                 # force re-extraction even when cached
)
# → data.frame with 38 numeric columns (19 vars × 2 aggregations):
#   bio1_median, bio1_iqr, bio2_median, bio2_iqr, ..., bio19_iqr,
#   plus n_extracted (integer, how many GBIF points contributed).
#   Rownames = species.
#   Species with no GBIF hits or no raster-valid points get all-NA covariates.
```

### 4.2 Data flow

```
species ──► pull_gbif_centroids()                        [B.1, already shipped in PR #46]
               │
               └──> gbif_cache_dir: one RDS per species
                    with occurrence records (lat/lon points)
                           │
              ┌────────────┘
              ▼
   pull_worldclim_per_species()
              │
              ├─ [first call only] download WorldClim raster stack
              │    to worldclim_cache_dir/wc2.1_10m_bio/
              │    (~500 MB global stack, 19 GeoTIFFs)
              │
              ├─ for each species:
              │    ├─ read GBIF cache RDS
              │    ├─ extract raster values at each occurrence point
              │    │   via terra::extract()
              │    ├─ aggregate (median + IQR per variable) across points
              │    └─ write per-species extract cache (RDS)
              │
              └─ assemble wide data.frame → return
```

### 4.3 Per-species extract cache

One RDS per species at `worldclim_cache_dir/extracts/{safe_key}.rds`
storing:

```r
list(
  species       = character(1),
  bio_median    = named numeric (19),      # names = "bio1".."bio19"
  bio_iqr       = named numeric (19),
  n_extracted   = integer(1),
  extracted_at  = POSIXct
)
```

Same key-sanitisation as B.1 (`[^A-Za-z0-9._-]` → `_`).

Species with zero GBIF occurrences (cached as NA by B.1) get a
per-species extract cache with `n_extracted = 0L` and all-NA bio values.

### 4.4 Raster download + cache

On first call, `pull_worldclim_per_species()` checks for
`{worldclim_cache_dir}/wc2.1_10m/wc2.1_10m_bio_1.tif` through
`_bio_19.tif`. If any are missing, downloads the zipped stack from:

```
https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip
```

(~130 MB compressed, ~500 MB unzipped). Uses `utils::download.file()`
with `method = "libcurl"` or `"auto"`, then `utils::unzip()`.

On success, writes a sentinel file
`{worldclim_cache_dir}/wc2.1_10m/.wc_complete` so subsequent calls
skip the download check entirely.

**Failure modes:**

- Network error / URL change: emit a clear error pointing to manual-download instructions.
- Partial unzip: detect missing GeoTIFFs and re-download.
- Disk full: let the download error propagate; user sees `cannot write file`.

### 4.5 Raster I/O with terra

`terra` moves to **Suggests** (not Imports) — B.2 is optional
functionality. Fail fast with install hint if absent:

```r
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("pull_worldclim_per_species() requires the 'terra' package: ",
       "install.packages('terra')", call. = FALSE)
}
```

Load the 19 bio rasters once per invocation via `terra::rast(list_of_paths)`
→ a single `SpatRaster` with 19 layers. Extract via:

```r
vals <- terra::extract(
  x = rast_stack,
  y = occ_matrix,       # n_points × 2, column order (lon, lat)
  ID = FALSE,
  method = "simple"     # nearest-neighbour; no bilinear interpolation
)
```

`vals` is a data.frame (n_points × 19). Aggregate per species via:

```r
med <- vapply(vals, stats::median, numeric(1L), na.rm = TRUE)
iqr <- vapply(vals, function(x) stats::IQR(x, na.rm = TRUE),
              numeric(1L))
```

### 4.6 Integration with pigauto

No pigauto internals change. User pattern:

```r
# One-time (cached):
gbif_df <- pull_gbif_centroids(rownames(df),
                                 cache_dir = "script/data-cache/gbif")
wc_df   <- pull_worldclim_per_species(
  species             = rownames(df),
  gbif_cache_dir      = "script/data-cache/gbif",
  worldclim_cache_dir = "script/data-cache/worldclim")

# Combined covariates (38 bioclim + 2 centroid = 40 cols):
cov <- cbind(gbif_df[, c("centroid_lat", "centroid_lon")],
              wc_df[, grep("^bio", colnames(wc_df))])

# Fit:
res <- impute(df, tree, covariates = cov)
```

### 4.7 File layout

```
R/worldclim_covariates.R             # new, ~220 lines
  ├── pull_worldclim_per_species()   # exported @export
  ├── .wc_download_rasters()         # internal: one-time raster download
  ├── .wc_extract_one()              # internal: per-species extract + cache
  ├── .wc_aggregate_one()            # internal: median + IQR from raster values
  └── .wc_cache_key()                # internal: filesystem-safe species key
```

### 4.8 Dependencies

Add to `DESCRIPTION` **Suggests**:
- `terra` — raster I/O

No new **Imports**. `utils::download.file()`, `utils::unzip()`, and
`stats::median` / `stats::IQR` are base-only.

## 5. Testing

### 5.1 Offline unit tests (`tests/testthat/test-worldclim-covariates.R`)

1. **Cache key sanitisation.** Reuse the same pattern as B.1 —
   confirm `.wc_cache_key("Quercus alba L.")` returns
   `"Quercus_alba_L_"` (drop `.` to match GBIF key convention).
2. **Extract aggregation.** Given a synthetic raster-value matrix
   (hand-crafted), `.wc_aggregate_one()` returns the correct
   median and IQR per variable.
3. **NA handling in aggregation.** Rows with any NA raster value
   (occurrence point outside raster extent) are dropped before
   aggregation; remaining rows produce correct median/IQR.
4. **Per-species cache read/write round-trip.** Hand-write a fake
   extract-cache RDS for 3 species; call
   `pull_worldclim_per_species()` with mocked raster stack;
   confirm no raster extraction is triggered.
5. **Fails cleanly when terra absent.**
   `skip_if(requireNamespace("terra", quietly = TRUE))`. Otherwise
   expect `stop()` with install hint.
6. **Fails cleanly when gbif_cache_dir missing.** No GBIF cache →
   per-species extract caches are NA, returned data.frame has
   all-NA bio columns for those species (graceful, not error).

All mocking via `testthat::local_mocked_bindings(rast = ..., extract = ...)`
to avoid touching `worldclim.org`.

### 5.2 Fixture build (`data-raw/make_worldclim_plants_300.R`)

One-shot builder that:

1. Reads the BIEN 300-species subset used by B.1 fixture.
2. Requires `script/data-cache/gbif/` to exist (B.1 fixture).
3. Downloads WorldClim 10m to `script/data-cache/worldclim/` (~500 MB,
   ~2–5 min).
4. Extracts per-species bioclim for all 300 species (~15–30 min wall;
   most time is raster I/O, not compute).
5. Writes `tests/testthat/fixtures/worldclim_plants_300.rds` (~50 KB
   compressed).

### 5.3 End-to-end smoke (`tests/testthat/test-worldclim-covariates.R`)

Gated by `NOT_CRAN=TRUE`. Uses the fixture above + BIEN cache + GBIF
fixture:

1. Load all three.
2. Subset to 200 species intersecting all fixtures.
3. 30 % MCAR mask on continuous traits.
4. Fit three times:
   - `res_none`: no covariates
   - `res_centroid`: GBIF centroids only (2 cols)
   - `res_bioclim`: GBIF centroids + 38 bioclim cols
5. Assert on **SLA and leaf_area** (the two most environment-driven):
   - `res_bioclim` RMSE ≤ `res_none` RMSE × 0.90 (**≥ 10 % lift**) on
     at least one of them.
   - `res_bioclim` r ≥ 0.30 on at least one of them.
6. Print the per-trait RMSE trio for manual review (centroid vs
   bioclim contribution).

The assertion threshold here (10 % lift, r ≥ 0.30 at n=200) is
**looser than the spec-level target** (r ≥ 0.35) because n=200 noise
is substantial. The full n=4,745 canary rerun enforces r ≥ 0.35 in
task 7.

### 5.4 Full plants bench rerun (manual, post-implementation)

`script/bench_bien_worldclim.R` (new): extends existing
`script/bench_bien.R` to include WorldClim covariates. Run at
n=4,745 full scale. Expected:

| trait | pre-B.2 RMSE | post-B.2 target | r target |
|---|---:|---:|---:|
| sla | 19.99 (at mean) | **≤ 15** | **≥ 0.35** |
| leaf_area | 12152 (at mean) | **≤ 9000** | **≥ 0.30** |
| height_m | 11.68 (at mean) | ≤ 9 | ≥ 0.30 |
| wood_density | 0.167 | ≤ 0.167 (non-regression) | ≥ 0.43 |
| seed_mass | 1750 (at BM) | ≤ 1500 | ≥ 0.20 |

This rerun happens AFTER the PR lands and WorldClim rasters are
downloaded. Committed as a benchmark run, not a hard test.

Success metric: `phylo_gate_triggered` should DE-trigger on
`sla`, `leaf_area`, `height_m` — the phylo-signal gate currently
fires on these three, but once the GNN has real covariate signal,
λ-only gating is no longer the right call. The gate-DE-triggering
is the architectural proof that the gate + covariate story is
coherent.

**Important nuance:** the phylo-signal gate in PR #45 is computed on
the trait's raw values, not on residuals-after-covariates. So λ will
still be ≈ 0 on those three traits even after adding bioclim — the
gate will still fire and force `r_cal_mean = 1`. For B.2 to actually
lift those traits, we need either:

**Option A (simpler):** Users disable the phylo-signal gate
(`phylo_signal_gate = FALSE`) when they have strong covariates.
Document this interaction in the B.2 NEWS entry and roxygen.

**Option B (more clever):** Compute λ on residuals-after-covariates.
Requires re-visiting phylo-signal-gate's internals.

**B.2 ships Option A** — simpler, transparent, low-risk. Option B
becomes a follow-up spec ("covariate-aware phylo-signal gating") if
the paper wants the two mechanisms to compose cleanly.

## 6. Risks and mitigations

### 6.1 WorldClim download fragility

`biogeo.ucdavis.edu` can be slow or down. Mitigation:

- Sentinel file `.wc_complete` means we download once per cache_dir.
- Document manual download path (user can drop the unzipped
  GeoTIFFs into `worldclim_cache_dir/wc2.1_10m/` themselves).
- Integration tests that require WorldClim are `skip_if_not()`.

### 6.2 Memory at full scale

19 rasters × global 10m resolution = ~40 MB each in memory when
loaded. Plus per-species extract cache (300 species × 19 vals × 2
aggregations = small).

At n=4,745 plants × up to 300 occurrences per species × 19 vars =
~27 M numeric cells at peak extract time = ~200 MB. Not a concern
on a 16 GB laptop.

### 6.3 Occurrence points falling outside raster extent

Ocean, island-not-mapped, extreme latitudes → raster returns NA.
`.wc_aggregate_one()` strips NA before median/IQR. If all points
are NA, per-species bio values are NA (graceful).

### 6.4 IQR = 0 for single-point species

Species with only 1 GBIF occurrence → IQR is 0 (not NA). That's
correct — range breadth is zero for a singleton. But 38 covariates
× several hundred IQR=0 columns is redundant information.
Acceptable for v1; users can drop `_iqr` columns if preferred.

### 6.5 Phylo-signal gate interaction (see §5.4)

Covered: ship Option A (user opt-out); Option B deferred.

### 6.6 Dependency creep

`terra` has a nontrivial install footprint (GDAL, PROJ). Only pulls
in as Suggests; users who don't need covariates aren't affected.

### 6.7 Licence compliance

WorldClim v2.1 is free for non-commercial research (CC BY-NC). We
don't redistribute rasters — users download directly. Add a note in
the roxygen: `@references` pointing to Fick & Hijmans 2017.

## 7. Rollout

Single PR on `feature/worldclim-covariates`. Stacks on PR #46.

1. Spec committed (this doc).
2. Plan committed (`plans/2026-04-24-worldclim-covariates.md`, next).
3. Implementation commits per plan (TDD, one task per commit).
4. Fixture build (manual, ~30 min live raster download + extract).
5. R CMD check clean; offline tests green.
6. End-to-end smoke on n=200 BIEN subset; assertions pass.
7. **Post-PR operator step:** full plants bench rerun at n=4,745 with
   bioclim, numbers added to the paper.

Version bump: `0.9.1.9004` → `0.9.1.9005`.

## 8. Success metrics

### 8.1 Hard (verified by PR tests)

- Offline unit tests: ≥ 10 expectations pass; R CMD check clean.
- End-to-end smoke at n=200: bioclim lifts at least one of
  SLA / leaf_area by ≥ 10 % RMSE over no-covariate baseline;
  r ≥ 0.30 on at least one.

### 8.2 Soft (verified by post-merge manual bench)

- Full n=4,745 plants bench: SLA RMSE drops from 19.99 (grand
  mean) to ≤ 15, r ≥ 0.35. leaf_area drops from 12152 to ≤ 9000,
  r ≥ 0.30.
- phylo-signal gate opt-out pattern documented with a worked
  example in the vignette or roxygen.

### 8.3 Paper claim

After B.2: "pigauto lifts plant trait prediction to r ≥ 0.35 on
SLA and leaf area when provided WorldClim covariates, compared
to r ≤ 0.21 on phylogeny alone. The safety floor + phylogenetic-
signal gate remain active: the covariate lift is *additional*
to the guaranteed-never-worse-than-mean floor."

## 9. Open follow-ups (separate specs)

- **B.3 SoilGrids** — add ~5-10 soil variables (pH, clay, SOC, bulk
  density). Harder pipeline (COG-based, not GeoTIFF), but
  complementary to bioclim.
- **Covariate-aware phylo-signal gating (Option B)** — recompute λ
  on residuals-after-covariates so the gate de-triggers
  automatically when covariates resolve weak-phylo-signal traits.
- **Per-observation bioclim (multi-obs mode)** — for benches like
  `ctmax_sim` where acclimation temperature varies within species.
- **CHELSA v2 as an alternative raster source** — user-selectable via
  a `raster_source = c("worldclim", "chelsa")` argument.
