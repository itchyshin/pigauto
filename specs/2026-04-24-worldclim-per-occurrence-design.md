# Per-occurrence WorldClim covariates — design spec (v1.1)

**Status:** approved (per user request 2026-04-24 "go with Track 2").
**Branch:** `feature/worldclim-per-occurrence` (stacked on `feature/worldclim-covariates`, PR #47).
**Plan:** `plans/2026-04-24-worldclim-per-occurrence.md` (to be written next).
**Target version:** v0.9.1.9006.

---

## 1. Goal

Fix the core limitation of B.2 (PR #47): bioclim was extracted at each
species' GBIF **centroid** (1 point per species), so the IQR columns
were all zero and the median was a single-point value. **Lift SLA /
leaf_area on plants to r ≥ 0.35** by extracting bioclim at *every*
GBIF occurrence point (up to 300 per species) and computing genuine
median + IQR across them.

Motivated by the honest-sim diagnostic (commit `04d11f3`, same
branch family) which showed pigauto's architecture CAN lift 32 % RMSE
on strong-environment traits — but only when covariates carry real
signal. Per-occurrence bioclim is the data-side fix.

## 2. What changes

### 2.1 B.1 side (GBIF cache)

The B.1 GBIF cache from PR #46 stores `(species, centroid_lat,
centroid_lon, n_occurrences)` per species. It does NOT store the
raw `(lat, lon)` points. To enable per-occurrence extraction we add
a **new optional cache field** `points`:

```r
# Per-species GBIF cache RDS layout, v1.1:
list(
  species          = character(1),
  centroid_lat     = numeric(1),      # kept for back-compat
  centroid_lon     = numeric(1),      # kept
  n_occurrences    = integer(1),      # kept
  points           = data.frame(lat = numeric, lon = numeric),  # NEW
  fetched_at       = POSIXct
)
```

Introduce a new optional argument `store_points = FALSE` on
`pull_gbif_centroids()`. When TRUE, `.gbif_fetch_one()` saves the raw
lat/lon data.frame alongside the centroid.

**Default stays FALSE** for back-compat: pre-v1.1 cache RDS files stay
readable unchanged; PR #46's end-to-end test keeps passing.

### 2.2 B.2 side (raster extract)

`.wc_extract_one()` currently reads `centroid_lat` / `centroid_lon` and
extracts at one point. v1.1:

- If GBIF cache has a `points` data.frame: extract bioclim at **every
  point** via `terra::extract(rast_stack, as.matrix(points[, c("lon",
  "lat")]))`, aggregate via `.wc_aggregate_one()`. IQR column is now
  meaningful.
- Else (legacy cache or `store_points` was FALSE): fall back to
  centroid-only extraction, log a warning note in cache metadata
  (`reason = "centroid_only_legacy"`). This preserves PR #47's smoke
  test behaviour.

### 2.3 Fixture rebuild

`data-raw/make_gbif_plants_300.R` (from PR #46) gets a rebuild with
`store_points = TRUE`. 300 species × up to 300 points each = up to 90k
lat/lon pairs stored. Cache-dir size grows from ~10 MB (centroid-only)
to ~40-80 MB (with points). Acceptable — the cache is
script/data-cache, not shipped with the package.

`data-raw/make_worldclim_plants_300.R` (from PR #47) re-runs to
regenerate `tests/testthat/fixtures/worldclim_plants_300.rds` with real
per-occurrence aggregation. The fixture RDS itself stays tiny (~20 KB
for 300 species × 38 numeric cols).

## 3. API changes

### 3.1 `pull_gbif_centroids()` (PR #46)

Add one argument (back-compatible):

```r
pull_gbif_centroids(
  species,
  cache_dir,
  occurrence_limit = 500L,
  sleep_ms         = 100L,
  verbose          = TRUE,
  refresh_cache    = FALSE,
  store_points     = FALSE   # NEW — default FALSE preserves v0.9.1.9004
)
```

Returned data.frame shape is unchanged (still just centroids +
n_occurrences). The `store_points = TRUE` flag ONLY changes what's
persisted to the cache RDS; downstream consumers of the return value
see the same columns.

### 3.2 `pull_worldclim_per_species()` (PR #47)

No signature change. Internal behaviour of `.wc_extract_one()` changes
to prefer per-occurrence when available. Return shape unchanged.

## 4. Non-goals

- Covariate-aware phylo-signal gate (Track 2 Option B) — demoted to
  nice-to-have after the honest-sim showed the safety-floor calibrator
  opens the GNN gate on its own when covariates help.
- SoilGrids extraction (B.3) — separate spec after v1.1 lands.
- Multi-obs per-observation bioclim — separate spec.
- CHELSA as alternative raster source — separate spec.

## 5. Testing

### 5.1 Offline unit tests

Add to `tests/testthat/test-gbif-centroids.R`:

- "pull_gbif_centroids(store_points = TRUE) writes points field" —
  mock `rgbif::occ_search`, confirm the per-species cache RDS has
  a `points` data.frame.
- "pull_gbif_centroids(store_points = FALSE) omits points" — default
  behaviour; cache RDS has no `points` field.
- "legacy cache (no points field) still reads cleanly" —
  hand-crafted legacy RDS without `points`; `pull_gbif_centroids()`
  reads it and returns correct centroids.

Add to `tests/testthat/test-worldclim-covariates.R`:

- ".wc_extract_one uses points when available" — fake GBIF cache with
  `points = data.frame(lat, lon)` of length 10; mock
  `terra::extract` to return a 10×19 matrix; aggregation uses all 10
  rows; `n_extracted = 10`.
- ".wc_extract_one falls back to centroid when points missing" —
  legacy-style GBIF cache; confirm it still works (single-point
  extraction), with `reason = "centroid_only_legacy"` recorded in
  the extract cache.

### 5.2 Fixture refresh

Two fixture rebuilds (manual, operator-run):

- `data-raw/make_gbif_plants_300.R` with `store_points = TRUE` to
  extend existing GBIF cache. Should be fast — we re-fetch the raw
  points but the cache is already partially populated (centroids),
  so only `points` arrays get added to existing RDS files.
- `data-raw/make_worldclim_plants_300.R` re-runs to re-extract at
  per-occurrence. Uses existing raster stack from PR #47. ~5 min
  for 300 species × ~200 points each.

### 5.3 End-to-end smoke

Extends the existing PR #47 smoke test with a new comparison:

```r
test_that("per-occurrence bioclim lifts plants SLA >= 10% over centroid-only", {
  # ... setup as before (NOT_CRAN, fixtures, n=200) ...
  # Fit once with centroid-only bioclim (old fixture contents) vs
  # per-occurrence bioclim (new fixture contents).
  # Assert: with-per-occurrence RMSE on SLA <= 0.90 * centroid-only RMSE.
})
```

If the fixture was regenerated in Task 6, both halves of this
comparison use the same data source but different aggregation. If lift
is real, the per-occurrence half wins by ≥ 10 % RMSE on SLA.

### 5.4 Full plants bench rerun (operator task, post-merge)

`script/bench_bien_worldclim.R` (already exists from PR #47) runs at
n=4,745. With per-occurrence bioclim:

| trait | target RMSE | target r |
|---|---:|---:|
| sla | ≤ 15 (baseline 20) | ≥ 0.35 |
| leaf_area | ≤ 9,000 (baseline 12,000) | ≥ 0.30 |
| height_m | ≤ 9 (baseline 12) | ≥ 0.30 |

These targets are inherited from PR #47's spec but were not achievable
with centroid-only. Per-occurrence is the unlock.

## 6. Risks and mitigations

### 6.1 Cache RDS size growth

Centroid-only cache: ~1 KB per species. With 300 points per species at
~16 bytes per (lat, lon) pair: ~5 KB per species. At scale (19,000
BIEN species), ~100 MB total in `script/data-cache/gbif/`. Cache-dir
only, not shipped. Acceptable.

### 6.2 Live fetch time

B.1's `pull_gbif_centroids()` already fetches raw occurrences per
species (needed to compute centroid). Storing them is a ~free side
effect of the existing fetch — no extra GBIF API calls.

### 6.3 Edge cases

- Species with < 10 valid occurrences: use all available points; don't
  gate. IQR will be noisier but still informative.
- Species with 0 occurrences: cache has `points = data.frame()`
  (zero-row), `.wc_extract_one()` handles it as current NA-row path.

### 6.4 Back-compat

Pre-v1.1 cache RDS files (no `points` field) keep working via the
fallback path. Users who fetched with `store_points = FALSE` can
re-run `pull_gbif_centroids(..., store_points = TRUE,
refresh_cache = TRUE)` to upgrade their cache.

## 7. Rollout

Single PR on `feature/worldclim-per-occurrence`. Stacks on PR #47.

1. Spec (this doc) committed.
2. Plan committed next.
3. Implementation commits per plan.
4. Fixture rebuilds (operator).
5. R CMD check clean.
6. Smoke test passes.

Version bump: `0.9.1.9005` → `0.9.1.9006`.

## 8. Success metric

After merge + fixture rebuild:

- Unit tests pass: `store_points = TRUE` writes `points` field;
  `.wc_extract_one()` uses them when available; legacy cache still
  works.
- End-to-end smoke: per-occurrence bioclim on SLA shows RMSE
  improvement ≥ 10 % over centroid-only at n=200 with
  `safety_floor = TRUE` (gate opens on the richer covariate input).
- Post-merge operator bench at n=4,745 reports **SLA r ≥ 0.35** on
  BIEN plants — the paper claim pigauto can lift weak-phylo-signal
  plant traits with environmental covariates.

## 9. Open follow-ups

- Covariate-aware phylo-signal gate (Option B) — now nice-to-have for
  interpretability, not a blocker. Can land after v1.1 if desired.
- SoilGrids (B.3) — separate spec.
- Foundation-model pretraining — long-term.
- Multi-obs per-observation bioclim — for ctmax_sim-style datasets.
