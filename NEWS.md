# pigauto 0.9.1.9008 (dev)

## NEW: `suggest_next_observation()` -- active-imputation guidance (2026-04-30)

For a fitted `pigauto_result` (from `impute()`), compute the
closed-form expected reduction in TOTAL predictive variance across
all currently-missing cells if each candidate cell were observed
next.  Useful for sampling-design guidance when measurement budget
is limited: which species, if observed next, would most reduce
imputation uncertainty across the rest of your tree?

### API

```r
suggest_next_observation(result, top_n = 10L,
                          by = c("cell", "species"))
```

Returns a `pigauto_active` data.frame (descending by
`delta_var_total`).  Columns:

* `by = "cell"`: `species`, `trait`, `type`, `delta_var_total`.
* `by = "species"`: `species`, `delta_var_total`,
  `n_traits_missing`.

### Math

Built on a Sherman-Morrison rank-1 inverse update of the BM
conditional MVN: adding species `s_new` to the observed set updates
`R_oo^{-1}` by a known closed form, and the variance reduction at
every other missing cell is a single scalar per candidate.

Closed form per candidate cell:

```
delta_V_total(s_new) = sigma2 * sum_{i in miss} D[i, k]^2 / alpha_k
```

where `D = R_mm - R_mo R_oo^{-1} R_om` is the residual matrix and
`alpha_k = D[k, k]` is the current relative leverage of cell `k`.

The formula assumes `sigma2` (REML variance) is held fixed -- the
standard active-learning convention.  Variance reduction depends
only on the geometry of which cells are observed, not on the
unknown value at `s_new`.

### Scope (v1.5, 2026-04-30 evening)

* Continuous-family traits (`continuous`, `count`, `ordinal`,
  `proportion`) -- variance reduction via the BM Sherman-Morrison
  closed form above.
* Binary + categorical traits (added 2026-04-30 evening) --
  expected entropy reduction via the label-propagation closed form:
  current LP probability \eqn{p_i} has entropy \eqn{H(p_i)};
  observing \eqn{s_{\mathrm{new}}} updates LP linearly via the
  similarity kernel; expected entropy after observation is averaged
  over \eqn{P(y_{\mathrm{new}})} = current LP estimate at
  \eqn{s_{\mathrm{new}}}.
* `zi_count` and `multi_proportion` not yet supported (silently
  skipped).
* Single-obs only.  Multi-obs deferred to v2.

**Variance and entropy are not directly comparable.**  The
default output mixes them by `delta` value, which is approximate
across metrics.  When precise ranking matters, filter by
`metric` ("variance" or "entropy") first.

### Why this is novel

To my knowledge, no other phylogenetic-imputation package (Rphylopars,
BACE, phylolm, mice with phylogenetic correlation) exposes a
sampling-design helper.  pigauto is uniquely positioned because its
BM conditional MVN already computes per-cell predictive variance --
this function just exposes the closed-form variance-reduction
formula derived from rank-1 update.  The methods are 30 years old
in optimal-design literature (Cohn et al. 1996, JAIR 4:129–145);
applying them to phylogenetic trait imputation appears to be new.

### Files

* `R/active_impute.R` (new) -- 304 LOC, hosts
  `bm_variance_reduction()` (internal) and
  `suggest_next_observation()` (exported), plus
  `print.pigauto_active` S3 method.
* `R/impute.R` -- adds `tree` to the `pigauto_result` so
  `suggest_next_observation()` can recompute `R_phy` if it was
  freed in the post-fit cleanup block.
* `tests/testthat/test-active-impute.R` (new) -- 10 tests, 29
  expects.  Verifies closed-form against an independent fixed-
  sigma2 brute-force refit on a 20-tip fixture; descending order;
  by-species aggregation; multi-obs error path; discrete-skip.

### Smoke

On AVONET 300 with 90 cells masked across 4 continuous traits:
top-1 cell is `Chloephaga_rubidiceps Beak.Length_Culmen`,
delta_var_total = 2.99; top-1 species is `Chloephaga_rubidiceps`
with 2 currently-missing traits combining to delta = 3.23.

## NEW: `gate_method = "cv_folds"` -- cross-validated gate calibration (2026-04-30)

Closes the longest-standing open NEWS item: *"A future improvement
(cross-validated gate selection, ...) could close this further; not
in scope for this release."*

### What changed

`calibrate_gates()` and `fit_pigauto()` gain `gate_method =
"cv_folds"` (opt-in alongside the existing `"single_split"` default
and `"median_splits"` option).  When selected, val cells are
partitioned into `gate_cv_folds = 5L` (default) deterministic
non-overlapping folds; the existing grid + half-B-verify procedure
runs once per fold (`half_a` = K-1 training folds, `half_b` =
held-out fold), and the componentwise median of K winning weight
vectors becomes `w_final`.  Post-cal full-val invariant fires
unchanged.

### Compared to `"median_splits"`

| Property | median_splits | cv_folds |
|---|---|---|
| half_a per split | 1/2 of n_val | (K-1)/K of n_val |
| Number of splits | B=31 random | K=5 deterministic |
| Overlap | Random (each cell in many splits) | None |
| Statistical interpretation | Median of B noisy gate selections | Standard k-fold CV |

Larger training sets per split → less noisy w_k.  Standard CV gives
a cleaner story for users and reviewers.

### Motivation

The strict val-floor enforces `pigauto_val <= baseline_val`, but
val→test drift survives in 4/32 binary cells (per
`useful/MEMO_2026-04-29_discrete_bench_reruns.md`).  Per-trait gate
calibration overfits to the specific val cells because the same
cells score every candidate weight.  Larger training sets per
split reduce overfitting; non-overlapping folds give a true held-
out evaluation.

### Edge cases

* `gate_cv_folds < 2` errors with a clear message.
* `n_val < gate_cv_folds` caps `K_eff = n_val` so each fold has at
  least 1 cell.
* `n_val < 2` (CV ill-defined) falls back to a single split.
* Existing back-compat preserved: default `gate_method =
  "single_split"`, and the `"median_splits"` path is bit-identical
  to pre-fix.

### Tests

5 new tests in `tests/testthat/test-safety-floor.R` (under
"# CV-fold gate calibration"):

* Valid simplex output
* Deterministic for fixed seed
* Pure-MEAN snap when MEAN dominates (post-cal invariant fires)
* `n_val < K` graceful fallback
* `gate_cv_folds < 2` error

All pass.  Existing 24 safety-floor tests + B4 edge cases unchanged.

### Smoke (AVONET 300, 60 epochs, 1 imputation, 60 Mass + 60 Migration cells masked)

| gate_method | wall (s) | r_cal_gnn[col_3] |
|---|---|---|
| single_split | 82.6 | 0.000 |
| median_splits | 77.2 | 0.050 |
| cv_folds | 75.9 | 0.150 |

`cv_folds` is slightly faster (smaller per-split work) and selects
a non-zero GNN gate on the col where the others snap to BM.

### Bench (focused 3-way comparison on synthetic mixed-type)

`script/bench_cv_vs_median.R` runs a focused 3-way comparison on
a synthetic continuous + binary + categorical fixture (BM, lambda=1,
n=300, 30 % mask, 5 reps).  Test-set means across reps:

| Method | bin (acc) | cat3 (acc) | cont (RMSE) |
|---|---|---|---|
| single_split  | 0.8222 | 0.7356 | 1.0778 |
| median_splits | 0.8222 | 0.7356 | 1.0720 |
| **cv_folds**  | 0.8222 | 0.7356 | **1.0611** |

* **Continuous**: cv_folds reduces RMSE by 1.6 % vs single_split
  and by 1.0 % vs median_splits.  Modest but consistent across
  the 5 reps (per-rep cv_folds < median_splits in 4/5 reps).
* **Discrete**: identical across all three methods on this
  fixture.  At lambda=1 the strict val-floor snaps the gate to
  the pure-BM corner regardless of how the gate was selected,
  so the three methods produce identical class assignments.
  This is consistent with what the discrete-bench memo
  (`useful/MEMO_2026-04-29_discrete_bench_reruns.md`) showed:
  the val→test drift on 4/32 binary cells is val→test
  extrapolation noise, not gate-selection noise; cv_folds does
  not close it on this regime.  The hypothesis that larger
  training sets per split would close MORE of the val→test
  drift than median_splits is **not supported** by this
  fixture.  cv_folds is preferable for continuous traits but
  not transformative for discrete.

Full per-type bench reruns (`bench_binary.R`, `bench_categorical.R`,
`bench_zi_count.R` at multiple signal levels) deferred to a
later session.

## Strict val-floor v3 — type-conditional safety floor (2026-04-29)

A type-restriction bug in `R/fit_helpers.R::calibrate_gates` was
silently exempting binary, categorical, and zi_count traits from the
post-calibration "strict val-floor" invariant.  The buggy behaviour
let pigauto underperform its own baseline on the package's own
simulator benches (`bench_binary.md`, `bench_categorical.md`,
`bench_zi_count.md`):

* binary `signal_0.6` mean accuracy: −4.8 pp below baseline
* categorical `K_3` cat2: −12 pp below baseline
* zi_count `zf_0.2` zi1: +23 % RMSE above baseline

These violate the safety-floor's "pigauto cannot underperform its
own baseline" promise.  The bug was introduced 2026-04-23 with the
initial safety-floor commit (`28d8e45`); the type filter mirrored an
earlier continuous-only MSE evaluator and was not extended when the
simplex grid landed.

The fix splits the post-calibration invariant into two
type-specific blocks:

* **Discrete family** (binary, categorical, zi_count) gets a strict
  full-val-set check against both pure-BM and pure-MEAN corners; if
  the calibrated blend's val-loss exceeds either corner by 1e-12,
  the gate is overridden to whichever pure corner is better.
  Discrete losses are 0-1 / argmax step functions, so this strict
  check is well-conditioned at typical val sizes.
* **Continuous family** (continuous, count, ordinal, proportion)
  keeps the pre-2026-04-29 mean-only check verbatim (override only
  if the blend loses to pure-MEAN, never to pure-BM).  An
  intermediate fix that also added a strict pure-BM check for
  continuous types over-corrected on high-phylo-signal traits where
  the blend's val MSE is numerically equal to pure-BM's val MSE
  within sampling noise (AVONET Mass test RMSE went 386 → 659 in
  the over-strict variant); v3 reverts continuous to the original
  conservative check.
* **multi_proportion** is skipped (gate naturally closes via the
  legacy path; bench shows pigauto = baseline at every signal
  level).

Empirical evidence (CORRECTED 2026-04-30 per the Opus adversarial
review; the original tabulation here mistakenly averaged val + test
splits in the bench RDS, overstating the closure).  The user-facing
held-out **test-set** mean acc gap pigauto − baseline (lower-magnitude
= closer to baseline-parity floor):

| bench | Apr 17 (pre-bug) | Apr 28 (buggy) | Apr 29 (post-fix v3) |
|---|---|---|---|
| binary | −0.029 | −0.027 | **−0.014** |
| categorical | −0.031 | −0.026 | **−0.016** |
| zi_count (RMSE ratio) | 1.058 | 1.043 | **1.012** |

The strict val-floor closes ~50 % of the test-side regression on
binary and categorical (from −2.7 / −2.6 pp to −1.4 / −1.6 pp); the
remaining ~1.4 pp is val→test sampling drift.  zi_count test ratio
drops from 1.06 to 1.01 -- most of the regression closes, with a
small residual.  3 of 32 binary cells still drift ≥ 5 pp on test
(val→test extrapolation that no val-only floor can prevent).

The strict val-floor's invariant is enforced at **calibration time**
on the calibration-time loss surface
(`compute_corner_loss(g, val_row_idx, ...)`).  At predict time the
model runs additional refine_steps and mask_token substitution that
introduce a small (~5 %) drift; the
`tests/testthat/test-safety-floor.R` strict-floor invariant test
documents this with a `* 1.05` slack.  So the user-facing guarantee
is "pigauto val-loss ≤ baseline val-loss + 1e-12 at calibration time,
modulo ~5 % refine_step drift at predict time" -- not the bit-exact
"≤ 1e-12 everywhere" reading.

Continuous benchmarks (AVONET n=1500 single seed × N_IMP=5):
the originally-reported "v3 within MC-dropout noise of pre-fix"
landed on a particularly close pair of draws.  Re-running the same
seed × n × N_IMP three times produced Mass test RMSE = 377, 430, 504
-- the comparison is non-deterministic at this evidence level
(pooling order over the 5 MC dropout draws varies between R sessions).
A multi-seed × N_IMP=20 verification is queued; the single-seed
single-run evidence does NOT support strong claims about whether
v3 changed continuous performance, only that the change (if any)
is within run-to-run pooling noise on this dataset.

Per-cell variance survives in multiple binary cells (val→test
sampling drift).  The strict val-floor only guarantees
`pigauto_val_loss ≤ baseline_val_loss + 1e-12` (calibration time);
it does not guarantee `pigauto_test_loss ≤ baseline_test_loss`.  A
future improvement (cross-validated gate selection, epsilon-margin
tightening, or LP added as a corner option for ordinal) could
close this further; not in scope for this release.

R CMD check after the fix: 0 errors / 0 warnings / 0 notes.
Full `devtools::test()`: FAIL 0 / SKIP 2 / PASS 1081.

New regression test in `tests/testthat/test-safety-floor.R`:
`"strict val-floor: pigauto val-loss <= baseline val-loss per
trait, all types"` — fits on a mixed continuous + binary +
categorical + count fixture, recomputes val-loss for the
calibrated blend and the pure baseline, asserts blend ≤ baseline ×
1.05 for every trait.

Calibration evidence preserved at
`script/_strict_floor_v3_calibration/` (4-way comparison: pre-fix /
v1 over-strict / v2 intermediate / v3 current).  Memos:
`useful/MEMO_2026-04-29_*.md`.

## Per-trait ordinal baseline path selection (Opus #6, 2026-04-30)

`fit_baseline()` now computes BOTH the threshold-joint baseline (B3,
commit `a541dbd`) AND a per-column BM-via-MVN alternative for each
ordinal trait, then picks the lower-val-MSE path per trait against
`splits$val_idx`.  The chosen path is exposed as
`fit_baseline()$ordinal_path_chosen` (named character: latent col →
`"threshold_joint"` or `"bm_mvn"`).

Motivation: the AVONET-300 Phase 6 bench's Migration RMSE regressed
from 0.879 (Apr 16) to 0.975 (Apr 29) when B3 swapped the ordinal
baseline path from a label-propagation-equivalent BM-via-MVN to the
truncated-Gaussian threshold-joint approach.  At K=3 ordinals the
K−1 thresholds are pinned to a narrow band by phylopars EM,
producing systematically worse predictions than BM-via-MVN on
z-scored integer class.  See
`useful/MEMO_2026-04-29_phase6_migration_bisect.md`.

Per Opus's adversarial review #6, we did **not** ship a `K ≤ 3 → LP`
heuristic.  Instead the per-trait selection runs at every K and
the gate's safety-floor simplex no longer has to route around
threshold-joint output for low-K ordinals -- the lower-val-MSE
baseline is selected before the gate calibration sees it.

Scope: single-obs only (multi-obs ordinal selection is out of scope
for this fix; it would require species-level aggregation of the BM
alternative).  The code path is gated by `!multi_obs && !is.null(splits)`
and only activates when threshold-joint actually populated the
ordinal column.

Two regression tests in
`tests/testthat/test-joint-threshold-baseline.R`:

* `"fit_baseline picks per-trait ordinal path and reports it"` --
  smoke test on K=3 ordinal + continuous mixed data; verifies finite
  output and that `$ordinal_path_chosen` is one of the two valid
  values.
* `"fit_baseline ordinal selection picks BM when threshold-joint
  loses on val"` -- recomputes both paths' val MSE independently
  and asserts the selection picked the lower one.

## Adversarial-review fixes (Opus 2026-04-28)

Four correctness fixes surfaced by the Opus adversarial pass on
top of PR #49 (the safety-floor / mean-gate work). All four are
silent-failure modes that produced plausible-looking but invalid
output; none of them threw a useful error before this release.

- **HIGH — `multi_impute(draws_method = "conformal")` no longer
  collapses to zero between-draw variance on shuffled input.**
  `R/multi_impute.R::.sample_conformal_draw` now threads
  `input_row_order` so it correctly maps user-input mask indices
  to internal (tree-tip) row indices when perturbing
  `pred$imputed`. Previously, on any data frame whose rows are
  not already in tree-tip order, the sampler perturbed the
  wrong cells and `build_completed` re-aligned the unperturbed
  values back into the user's missing cells, producing M
  identical datasets and `var = 0` after pooling. Any prior
  Rubin's-rules pooling on shuffled-input multi-imputation was
  unreliable. Caught by a new shuffled-input variance regression
  test; pre-fix `per_cell_sd ≡ 0`, post-fix `> 1e-8` at every
  masked cell.

- **MEDIUM (opt-in) — `conformal_split_val` argument on
  `fit_pigauto()` lets users avoid the validation-set
  double-dipping that breaks split-conformal exchangeability.**
  When `TRUE`, the validation cells are split per-column into a
  calibration half (used by `calibrate_gates()` to pick the
  blend weights) and a conformal half (used to estimate the
  residual quantile). Reusing the same val cells for both
  selects the gate to minimise residual MSE on the very cells
  whose residuals drive the conformal quantile, producing
  systematic undercoverage that's most visible at small `n_val`
  (matches the `coverage_investigation` memo).
  Default is `FALSE` (pre-fix single-set behaviour) because
  forcing the split regresses the AVONET300 / OVR-categorical /
  BIEN safety-floor smoke benches by 2-26% on small-val
  datasets. Use `conformal_split_val = TRUE` when accurate 95%
  coverage matters more than bench-grade RMSE.

- **MEDIUM — `calibrate_gates` half-A loop is now NA-safe.**
  Mirrors the half-B fix from commit 573decd: when
  `loss_a_pure_bm` is non-finite (e.g. multi_proportion CLR
  latents whose half-A subset is too small for some component),
  `best_la` is now seeded with `Inf` and the legacy
  pure-BM-relative-gain filter is skipped rather than throwing
  `missing value where TRUE/FALSE needed`.

- **LOW — `predict.pigauto_fit` no longer produces a cryptic
  torch shape error when calibrated_gates is overridden with a
  vector of the wrong length.** Length 0 is now treated as "no
  calibration" (degrades gracefully to the learned-gate path);
  any other non-matching length errors with an explicit
  `"calibrated_gates has length k but the model has p latent
  column(s)"` message. Unblocks the AE-attribution ablation
  workflow used by `script/phase1_gnn_ablation.R`.

Each fix has a dedicated regression test:
`tests/testthat/test-multi-impute.R::"multi_impute(draws_method=\"conformal\") produces non-zero between-draw variance on shuffled input"`,
`tests/testthat/test-fit-predict.R::"conformal_split_val toggle changes conformal scores (no double-dipping)"`,
and
`tests/testthat/test-fit-predict.R::"predict.pigauto_fit catches wrong-length calibrated_gates override"`.

# pigauto 0.9.1.9006 (dev)

## Per-occurrence WorldClim covariates (v1.1 follow-up to B.2)

Fixes the core limitation of PR #47 (centroid-only bioclim extraction):
bioclim is now extracted at every GBIF occurrence point per species
(not just the centroid), so the IQR columns carry real range-breadth
signal instead of being constant zero. This is what actually lets the
safety-floor calibrator see a meaningful covariate input and open the
GNN gate on plants.

Validated by the honest-sim diagnostic
(`experiment/covariate-honest-sim`) which showed pigauto's architecture
can produce 32 % RMSE lift on strong-environment traits when
covariates carry real signal — something centroid-only bioclim
failed to deliver on plants.

### API

- `pull_gbif_centroids(..., store_points = FALSE)` — new arg,
  default `FALSE` preserves pre-v0.9.1.9006 cache format. When
  `TRUE`, persists the raw filtered lat/lon points in each species'
  cache RDS, enabling per-occurrence bioclim extraction.
- `pull_worldclim_per_species()` — no signature change.  Internally,
  `.wc_extract_one()` now prefers cached per-occurrence points,
  falling back to centroid only when `points` field is absent
  (legacy caches).

### Back-compat

Pre-v1.1 cache RDS files remain fully readable. Users who want
per-occurrence extraction upgrade their cache with one call:

```r
pull_gbif_centroids(sp_list,
  cache_dir = "script/data-cache/gbif",
  store_points = TRUE,
  refresh_cache = TRUE)   # one-time refresh
```

After that, subsequent calls are instantaneous (cache hit).

### Fixture regeneration

`tests/testthat/fixtures/worldclim_plants_300.rds` re-generated with
per-occurrence aggregation. IQR columns are now informative.

### Follow-ups

- Covariate-aware phylo-signal gate — nice-to-have for
  interpretability (gate de-triggers when covariates resolve weak-
  phylo-signal traits). NOT a blocker given the safety-floor
  calibrator opens the gate on its own.
- Full-scale plants bench (n=4,745) — operator-run via
  `script/bench_bien_worldclim.R`. Expected: SLA r >= 0.35,
  leaf_area r >= 0.30.

# pigauto 0.9.1.9005 (dev)

## WorldClim bioclim covariates (B.2)

New exported helper `pull_worldclim_per_species(species,
gbif_cache_dir, worldclim_cache_dir, ...)` extracts 19 WorldClim v2.1
bioclim variables at each species' GBIF occurrence points, aggregates
per species (median + IQR), and returns a 38-column numeric
data.frame ready for `impute(..., covariates = ...)`.

### What it does

1. Resolves each species's GBIF occurrence cache from
   `pull_gbif_centroids()` (B.1).
2. Downloads WorldClim v2.1 10-arc-minute raster stack (~500 MB) to
   `worldclim_cache_dir` on first call; sentinel file makes subsequent
   calls offline.
3. For each species, extracts bioclim values at its centroid (v1;
   per-occurrence extraction deferred to v1.1) via `terra::extract()`.
4. Aggregates to median + IQR per variable (38 columns total).
5. Caches the per-species extract as one RDS in
   `worldclim_cache_dir/extracts/`.

### Important interaction with the phylogenetic-signal gate

The `phylo_signal_gate` introduced in v0.9.1.9003 fires on *raw trait
values*, so traits with weak phylogenetic signal (e.g. plants SLA,
leaf_area, height_m) get routed directly to the grand-mean corner
BEFORE the GNN sees any covariates.  **To see bioclim lift those
traits, users must set `phylo_signal_gate = FALSE`:**

```r
cov_bio <- pull_worldclim_per_species(
  rownames(traits),
  gbif_cache_dir = "script/data-cache/gbif",
  worldclim_cache_dir = "script/data-cache/worldclim")
res <- impute(traits, tree, covariates = cov_bio,
              phylo_signal_gate = FALSE)   # <-- required
```

A follow-up spec ("covariate-aware phylo-signal gating") will compute
lambda on residuals-after-covariates so the gate composes cleanly with
covariates without the manual opt-out.

### Dependency

`terra` added to Suggests.  `pull_worldclim_per_species()` errors
with an install hint when `terra` is absent.

### Testing

- `tests/testthat/test-worldclim-covariates.R`: 20 offline
  expectations (cache key, aggregation, download sentinel, extract
  cache, NA handling).
- `tests/testthat/fixtures/worldclim_plants_300.rds`: pre-extracted
  fixture for 300 BIEN species.
- End-to-end smoke (NOT_CRAN-gated): bioclim lifts SLA or leaf_area
  RMSE by >= 10% on a BIEN 200-subset when combined with
  `phylo_signal_gate = FALSE`.

### Follow-ups

- True per-occurrence extraction (v1.1): extend B.1's GBIF cache to
  store raw lat/lon point lists, then aggregate bioclim at each point
  rather than at the centroid.
- SoilGrids (B.3): add 5-10 soil variables similarly.
- Covariate-aware phylo-signal gating: recompute lambda on
  residuals-after-covariates so the gate de-triggers automatically
  when covariates resolve weak-phylo-signal traits.

# pigauto 0.9.1.9004 (dev)

## GBIF range-centroid covariates

New exported helper `pull_gbif_centroids(species, cache_dir, ...)`
fetches species occurrence records from the Global Biodiversity
Information Facility (GBIF) and returns a data.frame of per-species
median lat/lon centroids ready for `impute(..., covariates = ...)`.

### What it does

For each species, pull_gbif_centroids():

1. Resolves the taxon via rgbif::name_backbone() (handles synonyms).
2. Fetches up to occurrence_limit records via rgbif::occ_search(
   hasCoordinate = TRUE), paginated at 300/call.
3. Drops records with hasGeospatialIssues = TRUE,
   basisOfRecord in c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"),
   or out-of-range coordinates.
4. Computes median lat / lon per species (centroids more robust than
   means to outliers).
5. Caches one RDS per species at cache_dir (strongly recommended to
   avoid GBIF rate-limit issues).

### Why it matters

Plant traits like SLA and leaf area are environment-driven. Even the
simplest geographic covariate -- a species range centroid from GBIF --
carries biogeographic signal that phylogeny alone doesn't. This is the
B.1 proof-of-pipe for the pigauto covariate path; the full B.2 spec
(WorldClim bioclim + SoilGrids raster extraction) follows.

### End-to-end smoke (BIEN 200-subset)

Adding GBIF centroids as covariates keeps SLA and leaf_area RMSE
within +10% of the no-cov baseline (ratios ~1.09 and ~1.01 at n=200).
At this sample size the val-to-test sampling noise dominates any
covariate lift signal; production-grade positive lift reserved for
the B.2 bioclim run at larger n.

### API

- pull_gbif_centroids(species, cache_dir = NULL,
                        occurrence_limit = 500L, sleep_ms = 100L,
                        verbose = TRUE, refresh_cache = FALSE)
- Returns data.frame with species / centroid_lat / centroid_lon /
  n_occurrences, rownames = species.

### Dependency

rgbif added to Suggests (not Imports). pull_gbif_centroids() fails
with a helpful install hint when rgbif is absent.

### Testing

- tests/testthat/test-gbif-centroids.R: 32 offline expectations (cache
  key, aggregation, filtering, cache hit/miss, public function
  integration via mocks).
- End-to-end fixture at tests/testthat/fixtures/gbif_plants_300.rds
  (259/300 species with valid centroids) + live smoke test gated by
  NOT_CRAN=TRUE environment variable.

### Follow-up (separate spec)

B.2: full WorldClim bioclim + SoilGrids extraction at each GBIF
occurrence, aggregated per species. Expected to lift plants RMSE
significantly more than centroids alone.

# pigauto 0.9.1.9002 (dev)

## Safety floor: pigauto is never worse than the grand mean

Default change: `impute()` / `fit_pigauto()` / `multi_impute*` now run
with `safety_floor = TRUE`. The calibration grid post-training searches
a three-way convex combination of (BM baseline, GNN delta, grand mean)
on a 231-point simplex at step 0.05. Because the pure-mean corner
`(0, 0, 1)` is always a candidate in the grid, validation RMSE is
guaranteed by construction to satisfy `pigauto_val_RMSE <= mean_val_RMSE`
on every trait. On held-out test, the invariant holds up to
sampling-noise slack (tolerance +2% in the production canary,
verified on five taxa by `script/regress.R`).

### Effect on existing benches

- Birds, mammals, fish, amphibians: lift unchanged (within 1%).
- Plants: on weak-signal traits (`height_m`, `leaf_area`, `seed_mass`,
  `sla`) the calibrator opens `r_MEAN > 0` and the prediction falls
  toward the grand mean, fixing the boundary-case regression where
  pigauto was previously 15-101% worse than mean; `wood_density` lift
  preserved.

### API

- `impute(..., safety_floor = TRUE)` — new arg, defaults to TRUE.
- `fit_pigauto(..., safety_floor = TRUE)` — same.
- The fit object gains four new slots: `r_cal_bm`, `r_cal_gnn`,
  `r_cal_mean` (each a named numeric of length `p_latent`), and
  `mean_baseline_per_col`.
- `multi_impute()` / `multi_impute_trees()` propagate the mean term
  automatically via `predict.pigauto_fit()`; no signature change.
- Legacy v0.9.1 fit objects keep working via `%||%` fallback in
  `predict.pigauto_fit()`: when the new slots are absent, the 2-way
  blend `(1 - r_cal) * BM + r_cal * GNN` is reconstructed from the
  scalar `r_cal`.

### Opt-out

`safety_floor = FALSE` reproduces the v0.9.1 1-D calibration
bit-identically on the legacy path.

### Testing

- `tests/testthat/test-safety-floor.R`: smoke canary, ~3 min on MPS.
  23 `test_that` blocks, 75 expectations covering: simplex grid,
  mean-baseline scalar, calibrate_gates 3-way path, fit_pigauto
  integration, predict 3-way blend + legacy fallback, multi_impute
  propagation, multi_impute_trees share_gnn propagation, AVONET300
  regression (continuous RMSE within +2%, discrete accuracy within
  -2pp), plants safety smoke (cached BIEN n=1000, RMSE within +15%
  of grand mean).
- `script/regress.R`: full canary, ~20 min on MPS at n=1000 per taxon.
  Writes `script/regress_result.json` with per-bench pass/fail.

### Files

- `R/calibration_three_way.R`: new — `simplex_grid()`,
  `mean_baseline_scalar()` helpers.
- `R/fit_helpers.R`: `calibrate_gates()` now returns list of 3 weight
  vectors and optionally searches the simplex grid.
- `R/fit_pigauto.R`: computes `mean_baseline_per_col` from training
  cells; passes `safety_floor` through; stores the new slots.
- `R/predict_pigauto.R`: 3-way blend with `%||%` backward compat.
- `R/impute.R`: pass-through arg + roxygen.
- `R/multi_impute.R`, `R/multi_impute_trees.R`: roxygen only.
- `inst/extdata/legacy_fit_v091.rds`: fixture for backward-compat
  test.
- `data-raw/make_legacy_fit_v091.R`: one-shot fixture builder.
- `specs/2026-04-23-safety-floor-mean-gate-design.md`: design spec.
- `plans/2026-04-23-safety-floor-mean-gate.md`: implementation plan.

# pigauto 0.9.1.9000 (dev)

## Taxonomic-breadth benchmarks: 4 vertebrate classes + plant boundary case

Real-data benches now cover four vertebrate classes plus a plants
run that establishes the scope of phylogenetic imputation. All use
30% MCAR held-out, seed 2026, and report both conformal and
MC-dropout 95% coverage alongside RMSE / accuracy. The vertebrate
benches consistently lift point estimates by 27-75% RMSE /
10-40 pp accuracy; the plants bench is the first honest
boundary-case result, with a clean wood_density lift but r ≤ 0.21
on four other traits, attributable to weak phylogenetic signal in
pooled BIEN species means and random polytomy resolution in
V.PhyloMaker2 (details below). **Conformal coverage still lands at
0.90-0.98 on plants**, so the uncertainty story generalises across
the kingdom jump even where point estimates do not.

- **Birds** (`script/bench_avonet9993_bace_n3000.rds`, Vulcan L40S):
  AVONET n=3,000 subset x 7 traits. Mass RMSE -45%, Beak -55%,
  Tarsus -46%, Wing -57%; Trophic.Level accuracy +16 pp;
  Primary.Lifestyle +13 pp. Conformal 0.94-0.96.
  New driver `script/bench_avonet_full_local.R` runs the full
  n=9,993 locally on Mac MPS (Vulcan hit a predict-stage OOM at
  this scale); overnight job, HTML generator
  `script/make_bench_avonet_full_local_html.R` ready.
- **Mammals** (from `feature/bench-pantheria`): PanTHERIA n=4,027 x
  8 traits. **terrestriality accuracy 0.55 -> 0.95 (+40 pp)**;
  body_mass_g RMSE -27%, head_body_length -75%; diet/habitat
  breadth +10 pp each. Conformal 0.93-0.96.
- **Fish** (`script/bench_fishbase.rds`, Mac MPS):
  FishBase x fishtree n=10,654 x 6 traits.
  **Vulnerability RMSE -45%**, Length -38% (r=0.81), BodyShapeI
  accuracy +31 pp; Troph -23%, Depth -10%; Weight improved from
  r=0.04 (noise) to r=0.42 under the new median-pool MI fix.
- **Amphibians** (`script/bench_amphibio.rds`, Mac MPS): AmphiBIO
  n=5,237 x 2 continuous traits. **Body_size_mm RMSE 89 -> 44
  (-51%)**, Body_mass_g -27% (r=0.98). Continuous-only v1 --
  AmphiBIO binary / categorical columns hit a character/double
  type mismatch in the threshold-joint baseline, parked for
  v0.9.2. Calibrated gates = 0 on both traits: pigauto's GNN
  correctly recognises that the Level-C joint MVN baseline is
  already optimal, and still beats the mean baseline by 27-51%
  because the BM baseline captures cross-trait correlation +
  phylogenetic signal.
- **Plants** (`script/bench_bien.R`, honest boundary-case finding):
  BIEN x V.PhyloMaker2 at n=4,745 species, n_imputations=20.
  Only `wood_density` lifts (-6% RMSE, r=0.43); `height_m`, `sla`,
  `seed_mass`, `leaf_area` show r ≤ 0.21 and are 15-101% *worse*
  than grand mean on RMSE. Root cause is in the data, not pigauto:
  (i) BIEN species means are pooled from heterogeneous individual
  observations, diluting phylogenetic signal, and (ii) V.PhyloMaker2
  scenario S3 resolves within-genus polytomies randomly, so most
  species sit behind random branches relative to the Smith & Brown
  2018 backbone. pigauto's gate safety contains the damage: on
  weak-signal traits the calibrated gate closes to 0 and the
  pipeline degrades gracefully to the (also weak) BM prior rather
  than blowing up. **Conformal coverage still lands at 0.90-0.98
  on all 5 traits** -- the distribution-free uncertainty story
  holds even when point estimates are weak, separating the UQ
  narrative (always robust) from the point-estimate narrative
  (taxon-dependent). A first pass at `n_imputations = 1` showed
  `height_m` RMSE blowing up to 47.7 via Jensen exp()-decode bias;
  rerunning at `n_imputations = 20` activated the `dc8cffa`
  median-pool MI correction and dropped it to 15.15. The weak-
  signal result is what remains after that fix -- it is a property
  of BIEN x V.PhyloMaker2 at this species pool, not a pigauto
  regression. See `script/bench_bien.html` for the full honest
  reading and the in-page Jensen note.

Validation suite (`pkgdown/assets/validation_suite.html`) and the
paper-draft summary (`useful/results_summary.md`) collect all
numbers in one place.

## FishBase + fishtree benchmark (vertebrate breadth triad complete)

- New `script/bench_fishbase.R` benchmarks pigauto on the intersection
  of `fishtree::fishtree_phylogeny()` (Rabosky et al. 2018) and
  rfishbase's `load_taxa() + species() + ecology()` -- 10,654 species
  with 6 mixed-type traits (Length, Weight, BodyShapeI, DepthRangeDeep,
  Vulnerability, Troph). Completes the vertebrate breadth triad:
  birds (AVONET), mammals (PanTHERIA), fish (FishBase + fishtree).
- Full-scope result on Apple Silicon MPS (2.5 hr wall, post median-
  pool fix dc8cffa):
  Vulnerability RMSE 17.9 → 9.9 (-45%); Length RMSE 41.2 → 25.5
  (-38%, r=0.81); BodyShapeI accuracy 0.46 → 0.78 (+31 pp);
  Troph RMSE -23%; DepthRangeDeep RMSE -10%; Weight RMSE -3%
  (r=0.42 vs 0.04 under mean pool).
- HTML generator: `script/make_bench_fishbase_html.R` + pkgdown
  mirror at `pkgdown/assets/dev/bench_fishbase.html`.

## MI-pool robustness: median pooling for decode-amplifying traits

- `pool_imputations()` now uses **median pooling** across decoded
  values for all trait types where the decode step amplifies
  dropout-noisy draws: log-transformed continuous (new), count
  (extended), zi_count (extended), and proportion (new). Untransformed
  continuous keeps mean pooling (backward-compat default).
- Rationale: when `n_imputations > 1`, a few dropout-noisy latent
  draws with high-tail values decode via `expm1()` / `plogis()` to
  absurd magnitudes that dominate `rowMeans()`. Median is the
  outlier-robust pooling choice that matches PR #41's approach for
  count / proportion, now extended across the full set of
  decode-amplifying trait types.
- Concrete example: on the FishBase n=10,654 bench, pigauto's Length
  RMSE was 3,718 vs mean-baseline's 41 under the old mean-pool
  scheme -- a pure decode artefact (the calibrated gate on Length
  was 0). Median pool eliminates the artefact without retraining.
- Binary / categorical / ordinal pooling is unchanged (mean of
  probabilities / rowMeans of integer classes is correct there).
- 4 new regression tests in `test-mi-pool-robust.R` cover the four
  trait types that switched to median, verifying that a single
  extreme draw does not dominate the pooled output.

## GPU memory fix at predict stage (large n)

- `fit_pigauto()` now moves its returned `model_state` to CPU and
  calls `torch::cuda_empty_cache()` before returning, so the training
  run's activations, optimizer state, and per-epoch intermediate
  tensors can be reclaimed by torch's CUDA caching allocator. Prior
  to this fix, the fit object held GPU tensor refs that kept
  ~40 GB of training-time memory resident, and the downstream
  `predict()` path OOM'd on any card with <= 80 GB at n >= 5000.
- `impute()` now calls `gc()` + `torch::cuda_empty_cache()` between
  `fit_pigauto()` and `predict()` as belt-and-braces.
- Reproducer: Vulcan L40S (46 GB) at n=5000 or n=9993, default or
  compact config; predicted in 382 MB / 192 MB / 762 MB allocations
  on a 43 GB-already-used GPU. All four configurations on 2026-04-21
  died at the same step with the same root cause.
- No behavioural change for users: same `pigauto_fit` shape,
  `predict()` output identical; just no longer OOM at scale.

## Phase 7 EM: off-diagonal conditioning (opt-in, on top of Phase 6)

- New argument `em_offdiag = FALSE` on `impute()` and `fit_baseline()`.
  When `TRUE` AND `em_iterations >= 2L`, each binary/ordinal liability
  cell's prior at iteration `k + 1` is the full conditional-MVN
  `(mu, sd)` given the posterior liability of other traits at iteration
  `k`. This means observing one discrete trait shifts (not just
  tightens) the prior on correlated other traits — the off-diagonal
  entries of `Σ` now enter the E-step, instead of being ignored as in
  Phase 6.
- `em_offdiag = FALSE` default preserves Phase 6 behaviour byte-for-byte.
- `em_offdiag = TRUE` with `em_iterations < 2L` is silently clamped
  (iter 0 has no EM; iter 1 is plug-in with no previous Σ to condition
  on).
- **Scope**: binary + ordinal only. OVR categorical stays on Phase 6
  diagonal — going off-diagonal on OVR would re-couple the K classes
  and reintroduce the rank-(K − 1) singular-matrix instability that OVR
  was built to dodge.
- Missing-pattern handling: cells with all-NA other-trait posteriors
  fall back to the unconditional Phase 6 prior; cells with partial-NA
  use a restricted conditional on the observed subset.
- `em_state` gains `em_offdiag` field.

## Phase 6 EM for threshold-joint baseline (opt-in)

- New argument `em_iterations = 0L` on `impute()` and `fit_baseline()`.
  When `>= 1`, the threshold-joint path (binary + ordinal + OVR
  categorical) iterates: the BM rate `Σ` learned by
  `Rphylopars::phylopars()` at iteration `k` is fed back as the
  per-trait prior SD at iteration `k + 1`, up to `em_iterations` times
  (default 5 when enabled) or until
  `||Σ_k - Σ_{k-1}||_F / ||Σ_{k-1}||_F < em_tol` (default 1e-3). Closes
  the final open item from the v0.9.0 "Deferred to future releases" list.
- Default `em_iterations = 0L` preserves v0.9.1 output byte-for-byte.
- Under independent-discrete data the EM path converges in 1-2 iters
  at the plug-in answer (no free lunch, no harm). Under
  correlated-discrete data (ρ >= 0.5), the prior tightens and
  discrete-trait accuracy typically gains a few pp.
- Full off-diagonal conditioning on `Σ` is deferred to a future Phase 7.
  Phase 6 uses `sqrt(diag(Σ))` only.

---

# pigauto 0.9.1 (2026-04-19)

pigauto 0.9.1 is a consolidation release on top of 0.9.0. Two headline
improvements: (a) **tree-sharing GNN is now the default** in
`multi_impute_trees()`, which makes posterior-tree multiple imputation
cheap enough that the Nakagawa & de Villemereuil (2019) canonical
workflow (T = 50 posterior trees × m_per_tree = 1 = M = 50 pooled
datasets) is the default path at n = 10,000 species; (b) **opt-in
calibration smoothers** (bootstrap conformal, median-over-splits gate,
low-val-cells warning) for small-n regimes. This release also ships the
three post-0.9.0 feature branches (B1 soft-liability, B2 rate-aware
attention, B3 full-threshold ordinal baseline), scaling validation up
to 10,000 species, and a clean `R CMD check` (0 errors / 0 warnings /
1 note, down from 1 / 4 / 3 on v0.9.0).

## Calibration diagnostics and opt-in smoothing for small-n regimes

- New `fit_pigauto()` argument `min_val_cells = 10L`. Warns at fit time
  if any trait had fewer than this many validation cells available for
  gate calibration and conformal-score estimation. Default of 10 is the
  floor of pathological territory (operational target is 20-30).
- New optional `gate_method = "median_splits"` runs the half-A / half-B
  grid search on `gate_splits_B = 31L` random splits and takes the
  median best-gate. Small empirical benefit at small `n_val` (gate SD
  0.406 → 0.360 across 10 seeds on the calibration sim). Default stays
  `"single_split"` for backward compat.
- New optional `conformal_method = "bootstrap"` averages the conformal
  quantile across `conformal_bootstrap_B = 500L` bootstrap resamples of
  the val residuals. Empirically reduces conformal-score variance ~30%
  but **does not** reliably improve 95%-interval coverage stability on
  the evaluation sim. Shipped as experimental; default stays
  `"split"`.
- Roxygen for `fit_pigauto()` gains a `@section Calibration at small n:`
  explaining when the 95% intervals should be treated as approximate.
  Short version: at n_val < 20 the split-conformal guarantee degrades;
  increase `missing_frac` or collect more species rather than relying
  on estimator smoothing.

## Tree-sharing GNN is now the default in `multi_impute_trees()`

- New default `share_gnn = TRUE` in `multi_impute_trees()`. The GNN is
  trained once on a reference tree (MCC via `phangorn::maxCladeCred`,
  fallback `trees[[1]]`) and reused across posterior trees, with only
  the BM baseline recomputed per tree. ~10-15x speedup at n=10k x T=50
  makes the Nakagawa & de Villemereuil (2019) canonical workflow
  (T=50, m_per_tree=1, M=50 pooled via Rubin's rules) the default path.
- Default `m_per_tree` flipped from `5L` to `1L` to align with the N&dV
  2019 canonical workflow. Users with `T < 20` get a runtime warning
  suggesting they bump `m_per_tree` to keep Rubin's rules stable.
- Opt-out: `share_gnn = FALSE` restores the pre-v0.9.1 per-tree fit
  path for users who need exact per-tree model independence.
- New optional arg `reference_tree` lets users override the MCC choice.
- `predict.pigauto_fit()` gains a `baseline_override` argument (mostly
  internal — used by the shared-GNN path).
- New Suggests: `phangorn` (for `maxCladeCred()`).
- Tree-uncertainty propagation analysis: when the calibrated gate
  closes (the common case on real data — see the v0.9.0 validation
  suite), the shared-GNN approximation is lossless. When the gate is
  open, tree variance is slightly under-estimated in the GNN channel
  only — the baseline channel still carries it.

## Full threshold-model ordinal baseline (B3)

- Ordinal traits now use proper interval-truncated Gaussian E-step
  in the Level-C baseline instead of z-scored integer passthrough.
  Observing "class 3 out of 5" constrains the liability to the
  interval [threshold_2, threshold_3] rather than treating it as a
  point value. The liability joins the joint Rphylopars fit alongside
  continuous + binary traits.
- New internal decoder: `decode_ordinal_liability()` converts the
  Rphylopars liability posterior back to z-scored integer class for
  downstream GNN compatibility.
- `fit_baseline()` threshold-joint dispatch now fires when ordinal
  cols are present (in addition to binary / categorical). Ordinal
  cols that can't be populated (e.g., <2 observations) fall back to
  per-column BM unchanged.

## Rate-aware message passing via learnable per-head Gaussian bandwidth (B2)

- `GraphTransformerBlock` attention bias now uses learnable per-head
  Gaussian bandwidth computed from raw squared cophenetic distances:
  `bias_h = -D_sq / (2 * softplus(log_bw_h)^2)`. Each head can learn
  its own phylogenetic scale — tight for fast-evolving traits, broad
  for conserved traits. Replaces the fixed `adj_bias_scale * log(adj)`
  formulation.
- `build_phylo_graph()` returns `D_sq` (squared cophenetic distances)
  alongside `adj`. Persists through training (not freed like `D`).
- Backward compat: when `D_sq = NULL` (pre-B2 saved fits), the old
  `log(adj)` path fires unchanged. Legacy attention path (`use_transformer_blocks = FALSE`) is completely unaffected.

## Soft-liability E-step for multi-obs aggregation (B1)

- New argument `multi_obs_aggregation = c("hard", "soft")` on `impute()`
  and `fit_baseline()`. Default `"hard"` preserves v0.9.0 behaviour.
  `"soft"` uses Rao-Blackwell convex-combination E-step for binary and
  categorical traits: a species with 6/10 class-1 observations now
  produces a less-extreme posterior liability than 10/10, fixing the
  threshold-at-0.5 / argmax-one-hot information loss from Phase 10.
- New helpers: `estep_liability_binary_soft(p, mu_prior, sd_prior)` and
  `estep_liability_categorical_soft(p_vec, mu_prior, sd_prior)` in
  `R/liability.R`. Both reduce to their hard counterparts at boundary
  (one-hot) proportions; both return the prior mean at maximum ambiguity
  (p = 0.5 for binary, uniform for categorical).
- `aggregate_to_species()` gains `soft_aggregate` flag: when TRUE,
  preserves species-level proportions instead of thresholding/argmax.
- Benchmark: binary accuracy consistently improves +1.4–2.4pp across
  all 3 phylogenetic-signal regimes. Categorical shows +2.0pp at high
  signal but mixed results at moderate/low signal (hard threshold
  sometimes snaps to the correct majority class by luck).

## Scaling validation (up to 10,000 species)

- **AVONET 9,993 + BirdTree missingness sweep** (`script/bench_avonet_missingness.R`).
  Mean / mode vs BM baseline vs pigauto at 20% / 50% / 80% missing, run
  from Compute Canada (Narval). pigauto beats the BM baseline on every
  continuous trait at 80% missing; discrete accuracy stays within ≈1pp
  of BM across the sweep. Output: `bench_avonet_missingness.rds` +
  `bench_avonet_missingness.md`, rendered into the pkgdown validation
  suite.
- **Scaling curve to n = 10,000** (`script/bench_scaling_v090_extended.R`).
  Per-stage wall-clock + peak memory at
  n ∈ {100, 300, 1000, 3000, 5000, 7500, 10000}. Confirms sub-quadratic
  growth for the baseline and the expected O(n²) scaling of the GNN
  attention. At n = 10,000 the full pipeline completes in ~30–60 min on
  a single CPU node. Output: scaling table + curve PNG in the pkgdown
  validation suite.

## Benchmarks

- `bench_covariate_sim.R` rerun on v0.9.0 (previously pinned to v0.6.0).
  Now the reference entry for the environmental-covariates path in the
  validation suite.
- `bench_multi_obs_mixed.R` output regenerated with the B1 soft E-step
  path enabled; captures the +12–28pp categorical lift under soft
  aggregation vs hard.

## Documentation

- New pkgdown article **"Propagating Tree Uncertainty"**
  (`vignettes/tree-uncertainty.Rmd`). Decision guide for single-tree vs
  multi-tree MI, canonical N&dV 2019 workflow under the new
  `share_gnn = TRUE` default, worked example over posterior trees with
  `Map()`-based downstream-model refitting, and a timing table.
- pkgdown trait-type consistency pass: `zi_count` and `multi_proportion`
  are now listed everywhere alongside the other six trait types
  (Getting Started, Mixed Types, the validation suite, and the README).

## Internal

- **`R CMD check` clean**: 0 errors / 0 warnings / 1 note, down from
  1 / 4 / 3 on v0.9.0. The remaining note is the `BACE` in-tree
  "Package sources not in a standard location" flag, which is a
  deliberate `.Rbuildignore` choice so BACE ships as a comparison
  reference without being part of the installed pigauto package.
- **Compute Canada SLURM bundles**. `submit_v090_cloud/` (Narval;
  renamed clusters after the 2026 Alliance Canada infrastructure
  renewal — Cedar → Fir etc.) ships one-shot scripts for the AVONET
  missingness sweep and the scaling curve. `submit_v090_vulcan/` (PAICE
  account `aip-snakagaw`) ships a SLURM array driver for the 60-cell
  calibration-coverage grid.
- **pkgdown CI path filter**: the pkgdown workflow now only rebuilds on
  changes to user-facing files (`R/`, `vignettes/`, `pkgdown/`,
  `README.md`, `NEWS.md`, `DESCRIPTION`). Saves ~50% of the
  GitHub Actions minutes previously consumed by pigauto PR checks.
- **Gap fixes**: replaced a hardcoded worktree path in
  `bench_multi_obs_mixed.R` with a portable relative path; refreshed
  the bench_multi_obs_mixed `.rds`.

---

# pigauto 0.9.0 (2026-04-16)

pigauto 0.9.0 is a large release that ships two complementary upgrades:
a **Level-C phylogenetic baseline** (captures cross-trait correlation
and threshold structure), and a **Graph Transformer GNN backbone**
(replaces the 2-layer single-head attention stack with a proper
multi-head transformer). Both are gate-safe — pigauto still closes its
per-trait gate when the baseline is already optimal, so high-phylo-signal
datasets retain identical behaviour to v0.7.0.

## Level-C baseline improvements

A new family of joint multivariate baselines replaces per-column
Brownian motion + label propagation on datasets where cross-trait
correlation carries signal.

- **Liability encoding (`R/liability.R`)**. Every trait type now has an
  explicit liability interpretation: continuous types are their own
  liability; binary uses a 1D threshold at 0; ordinal uses K-1
  thresholds; categorical uses K liabilities with an argmax constraint.
  A single dispatcher `estep_liability(tm, observed, mu_prior, sd_prior)`
  returns posterior mean/variance of the underlying liability given an
  observation. This is the foundation the joint paths below build on.

- **Joint multivariate BM baseline (`R/joint_mvn_baseline.R`)**. When
  Rphylopars is installed and the dataset has ≥2 BM-eligible latent
  columns, the baseline now fits a single joint multivariate BM across
  continuous, count, ordinal, proportion, and ZI-magnitude cols via
  `Rphylopars::phylopars()`, capturing phylogenetic correlation across
  traits. Benchmark on correlated simulated BM data:
  **33.7% RMSE lift** vs per-column BM (`bench_joint_baseline.R`).

- **Threshold-joint binary baseline (`R/joint_threshold_baseline.R`)**.
  Binary traits get a truncated-Gaussian E-step posterior mean via
  `estep_liability_binary` and join the joint MVN alongside continuous
  traits. Decoded back to `logit(P(y=1))` so BCE and gate math are
  unchanged. Benchmark: 3/3 scenarios with ≥2pp accuracy lift over the
  LP baseline, peaking at **+16pp** at rho = 0.6
  (`bench_binary_joint.R`).

- **OVR categorical baseline (`R/ovr_categorical.R`)**. Each K-class
  categorical trait is decomposed into K independent "is_class_k vs
  rest" binary threshold fits, then normalised into a row-stochastic
  distribution. This is the same OVR strategy BACE uses. On AVONET 300,
  this lifts Trophic.Level accuracy to **77%** and Primary.Lifestyle to
  **84%** — beating BACE-OVR's 72% on both. The K-independent-fits path
  sidesteps the rank-(K-1) numerical instability that made earlier
  single-fit approaches unstable; each individual fit has only one
  categorical-related column.

- **ZI-count gate routing**. The binary gate col of zero-inflated count
  traits now joins the threshold-joint path alongside binaries (the
  magnitude col already rode the Phase-2 joint MVN).

- **Graceful fallback throughout**. Any Level-C path that can't fit
  (Rphylopars missing, too few observations per column, multi_proportion
  trait present) falls through to the legacy per-column BM + LP paths
  without user-visible breakage.

## GNN improvements

- **Graph Transformer backbone (`R/graph_transformer_block.R`)**. The
  2-layer single-head attention stack is replaced by `n_gnn_layers`
  instances of a `GraphTransformerBlock` — multi-head attention (default
  4 heads) with learnable per-head `log(adj + eps)` bias (so the phylo
  prior is respected by each head independently), feed-forward network
  (default 4× width), pre-norm, two residual skips. FFN output
  projection initialises to zero so each block starts near-identity,
  preserving pigauto's gate-closed-at-init safety.

- **New hyperparameters on `ResidualPhyloDAE`** and `fit_pigauto()`:
  `use_transformer_blocks = TRUE` (new default), `n_heads = 4L`,
  `ffn_mult = 4L`. The legacy architecture is fully preserved behind
  `use_transformer_blocks = FALSE` and reproduces pre-0.9.0 numbers
  exactly. Pre-0.9.0 saved `pigauto_fit` objects reconstruct via the
  legacy path automatically.

- **Empirical characterisation**. On the 4-scenario discriminative
  benchmark the transformer and legacy GNN produce identical predictions
  in high-signal regimes (gate closes to zero — architecture below the
  gate is irrelevant) and the transformer shows small RMSE lift in
  moderate/low signal scenarios where the gate stays partly open.
  Runtime is ~25–30% slower per fit. Larger transformer gains likely
  need self-supervised pretraining on a large tree corpus (future work).

## Multi-observation unlock

- **Multi-obs data now uses Level-C paths**. Previously the joint MVN,
  threshold-joint, and OVR categorical paths refused multi-obs input via
  four `!multi_obs` guards in `fit_baseline()`. Those guards are all
  removed; multi-obs data is aggregated to species level inside each
  Level-C helper via a new internal `aggregate_to_species()` helper.
  Aggregation rules: mean for continuous-family, threshold-at-0.5 for
  binary/ZI-gate, argmax-one-hot for categorical. Splits are
  mask-then-aggregated so val/test leakage is impossible; single-obs
  data passes through unchanged.

- **Benchmark (`bench_multi_obs_mixed.R`)**. On a synthetic multi-obs
  mixed-type simulation (150 species, 2 continuous + 1 binary + 1 K=4
  categorical, Poisson-5 obs per species), the Level-C path beats
  the LP path by **+12–28pp categorical accuracy** and **+5–7.5pp
  binary accuracy** across all three phylogenetic-signal regimes.

- **Aggregation caveat**. Binary/categorical aggregation is lossy —
  a species with split observations (e.g., 6/10 class-1) becomes a
  single class-1 observation at species level. For datasets with
  strong within-species variability on discrete traits, prefer
  single-obs mode with a user-computed species-level summary.

## Benchmarks

Nine new or reran benchmark reports quantifying the release's impact:

- `bench_joint_baseline.R` — Phase-2 joint MVN A/B on correlated BM data.
- `bench_binary_joint.R` — Phase-3 threshold-joint A/B on binary traits.
- `bench_discriminative.R` + `bench_discriminative_phase9.R` — 4-scenario
  mixed-signal sweep, including the transformer vs legacy GNN A/B.
- `bench_avonet_phase6.R` — post-Phase-6 AVONET 300 vs LP baseline.
- `bench_categorical_joint.R` — OVR K-fits A/B on synthetic categorical
  data.
- `bench_covariate_sim.R` — rerun on 0.9.0. Shows **−10.6% RMSE** on the
  low-phylo-strong-env scenario vs pre-Phase-6. Pre-0.9.0 snapshot in
  `bench_covariate_sim_preP6.md`.
- `bench_multi_obs.R` — rerun on 0.9.0. Numbers unchanged from pre-0.9.0
  (bench is single-trait so it can't exercise Level-C dispatch
  conditions; `bench_multi_obs_mixed.R` is the correct validator).
- `bench_multi_obs_mixed.R` — new synthetic mixed-type multi-obs
  benchmark specifically to validate Phase 10's multi-obs unlock.

## Internal

- `joint_mvn_available()` gates the Level-C dispatch — returns `TRUE`
  iff Rphylopars is installed. All Level-C paths fall back gracefully
  when it returns `FALSE`.
- `aggregate_to_species()` is the sole multi-obs → species-level helper
  and is the right place to hook future soft-liability aggregation
  (currently threshold-at-0.5 for binary, argmax for categorical).
- 20+ new test blocks across `test-liability.R`,
  `test-joint-threshold-baseline.R`, `test-joint-baseline.R`,
  `test-ovr-categorical.R`, `test-graph-transformer-block.R`,
  `test-phase9-integration.R`, and `test-multiobs-levelc.R`. Full suite
  is 778 PASS, 0 FAIL.

## Deferred to future releases

- Soft-liability E-step for multi-obs aggregation (addresses the
  threshold-at-0.5 lossiness documented above).
- Rate-aware message passing in the transformer blocks.
- Full threshold-model ordinal baseline (ordinal currently rides the
  Phase-2 joint MVN via its z-scored integer encoding).
- Self-supervised pretraining of the transformer backbone on a large
  tree corpus.
- CRAN submission preparation (BACE in-tree separation, torch handling,
  absolute-path cleanup in `script/`).

# pigauto 0.7.0

## New trait type: multi_proportion (compositional data)

- **`multi_proportion`**: a new trait type for compositional data — K columns
  per row that sum to 1 (e.g. plumage-colour proportions, diet composition,
  microbiome relative abundances, allele frequencies). Declared via a new
  `multi_proportion_groups = list(colour = c("black", ..., "yellow"))`
  argument to `impute()`, `preprocess_traits()`, and `multi_impute()`.

- **Encoding**: centred log-ratio (CLR) + per-component z-score. Small
  epsilon (`1e-6`) is applied to zero components before log to keep the
  transform safe; rows are re-normalised so the composition is preserved.

- **Baseline**: K independent Brownian-motion imputations on the CLR
  columns. BM is well-defined on CLR space (Euclidean), unlike the raw
  simplex.

- **GNN output**: K logits projected to sum-zero, then softmax at decode
  time — guaranteed to lie on the simplex.

- **Loss**: MSE on CLR columns (averaged across K), applied group-wise so
  the whole composition is masked together during DAE corruption.

- **Metrics**: Aitchison distance (the natural compositional-data metric),
  RMSE on CLR space (comparable to continuous traits), and simplex MAE
  (mean absolute error on decoded proportions).

- **Benchmark**: `script/bench_multi_proportion.R` with 5-scenario
  signal sweep and K ∈ {3, 5, 8, 12} secondary sweep. Report rendered
  by `script/make_bench_multi_proportion_html.R`.

- **Tests**: `tests/testthat/test-multi-proportion.R` (5 test blocks,
  20 expectations) — encoding, validation, baseline, end-to-end `impute()`.

## Internal

- `evaluate_imputation()` output adds three columns: `aitchison`,
  `rmse_clr`, `simplex_mae`. Existing metrics columns are unchanged.

- `impute()` gains explicit `trait_types` and `multi_proportion_groups`
  arguments (previously `trait_types` was documented as passed through
  `...` but was not actually forwarded; this fixes that).

# pigauto 0.6.2

## Documentation

- **Tree uncertainty workflow clarified.** `multi_impute_trees()`,
  `trees300`, the README, and `vignette("getting-started")` now all
  describe the two-step workflow explicitly: step 1 is tree-aware
  imputation (pigauto's job — already done correctly); step 2 is
  tree-aware downstream analysis (the user's responsibility,
  following Nakagawa & de Villemereuil 2019, *Syst. Biol.* 68:632–641).
  Every doc now includes a complete `Map()`-over-`mi$tree_index` code
  example showing how to run the corresponding tree in the downstream
  model.
- **Compute-cost table added** to the tree-uncertainty docs: wall-clock
  budgets at `n` = 300 / 5,000 / 10,000 species with `T` = 10 / 50
  trees, plus guidance on reducing `T` for large trees and
  parallelising across machines.

## Internal

- No API changes; no numerical changes to any fitted model. Doc-only
  patch release.

# pigauto 0.6.1

## Bug fixes

- **MC dropout zero variance fixed** (`R/predict_pigauto.R`): when
  `n_imputations > 1`, the blend equation now uses a BM posterior draw
  `t_BM_draw ~ N(BM_mu, BM_se)` instead of the deterministic `t_MU`.
  When the calibrated gate is zero (BM dominates), each imputation now
  samples from the correct BM posterior rather than returning the same
  point estimate, giving non-zero between-imputation variance.

- **Conformal draw scale bug fixed** (`R/multi_impute.R`): for
  log-transformed traits (Mass, Wing length, etc.), draws are now made
  on the latent z-score scale and back-transformed, avoiding near-zero
  variance caused by dividing a log-scale SE by an original-scale value.
  The same fix applies to log1p-scaled counts and logit-scaled
  proportions.

## New features

- `draws_method = "conformal"` is now the default for `multi_impute()`.
  Draws sample from the split-conformal calibrated uncertainty
  distribution (better calibrated than MC dropout for downstream
  Rubin's-rules pooling).
- `draws_method = "mc_dropout"` remains available and is now correct
  (see bug fix above).

## Hex sticker

- New `man/figures/logo.png`: pink pig driving a car with a 7-tip
  ultrametric pectinate cladogram above the pig's head. Forest green /
  mint palette (v3) selected as the official sticker.

# pigauto 0.6.0

## Observation-level covariate refinement for multi-obs data

When datasets have multiple observations per species measured under
different conditions (e.g. CTmax at different acclimation temperatures),
the GNN now produces **covariate-conditional predictions within species**.
A new refinement MLP (`obs_refine`) re-injects user-supplied covariates
after species-level phylogenetic message passing, so that different
observations of the same species receive different predicted values
based on their covariate context. This is controlled by a residual
connection that preserves the phylogenetic signal from message passing.

- `R/model_residual_dae.R`: new `n_user_cov` parameter and `obs_refine`
  MLP (activated only when `species_col` + `covariates` are both supplied)
- `R/fit_pigauto.R`: `n_user_cov` stored in model config
- `R/predict_pigauto.R`: backward-compatible model reconstruction

## New bundled dataset: `ctmax_sim`

Simulated multi-observation-per-species CTmax data (1,464 observations
across 300 species, with acclimation temperature as an observation-level
covariate). 30% of species are entirely unobserved. Uses `tree300`.
See `data-raw/make_ctmax_sim.R` for the generation script.

## New benchmark: multi-observation imputation

`script/bench_multi_obs.R` simulates CTmax-like datasets under varying
phylogenetic signal and within-species acclimation response ratios,
comparing species_mean, pigauto (no covariates), and pigauto (with
covariates). The observation-level refinement gives measurable RMSE
improvement when within-species covariate effects are strong.

---

# pigauto 0.5.0

## Three sources of information

pigauto now explicitly combines three sources of information for
imputation: (1) the phylogenetic tree, (2) cross-trait correlations,
and (3) optional environmental covariates. The `covariates` argument
is accepted by `impute()`, `multi_impute()`, `multi_impute_trees()`,
and the lower-level pipeline functions. Covariates are threaded
through the GNN with the same gated safety that protects against
GNN degradation — they only contribute when they demonstrably
improve accuracy.

## Internal BM baseline (Rphylopars no longer required)

The Brownian-motion baseline for continuous, count, ordinal, and
proportion traits is now computed internally using conditional
multivariate normal imputation on the phylogenetic correlation matrix
`R = cov2cor(vcv(tree))` (Goolsby et al. 2017). This removes the
hard dependency on `Rphylopars`, which is now in `Suggests` only.
The internal implementation (`R/bm_internal.R`) uses univariate
per-column imputation with Cholesky decomposition and nugget
regularisation for near-singular submatrices. Validation against
`Rphylopars` shows r = 0.97–0.98 (expected given univariate vs
multivariate difference).

## Tree-uncertainty MI via `multi_impute_trees()`

New function `multi_impute_trees()` performs multiple imputation
across a posterior sample of phylogenetic trees, so that phylogenetic
uncertainty propagates into downstream standard errors via Rubin's
rules. Demonstrated with the bundled `trees300` dataset (10 BirdTree
posterior trees). Benchmark results show SE inflation of 1.1–2.1x
and fraction of missing information (FMI) rising from ~0.03
(single-tree) to 0.22–0.79 (multi-tree) depending on missingness
level.

## New bundled datasets

- `trees300`: 10 posterior trees from the BirdTree Hackett backbone,
  pruned to the same 300 tips as `tree300`
- `delhey5809`: plumage lightness and 4 environmental covariates
  (absolute latitude, mean temperature, mean precipitation,
  tree cover) for 5,809 bird species (Delhey 2019)
- `tree_delhey`: matching BirdTree phylogeny for `delhey5809`

## Benchmark suite

Thirteen benchmark reports are now accessible from the pkgdown site
under Articles > Benchmarks:

- Per-type benchmarks: continuous, binary, ordinal, count,
  categorical, proportion, zero-inflated counts
- Missingness mechanisms: MCAR vs MAR (trait/phylo) vs MNAR
- Tree uncertainty: single-tree vs multi-tree MI
- Environmental covariates: Delhey 5,809-species real data
  (high phylo signal, covariates give ~zero lift) and simulation
  sweep (4 lambda x beta scenarios; 8–15% RMSE reduction when
  env effects are strong, ~0% when phylo signal dominates)

Each report is a self-contained HTML page generated from
`script/make_bench_*_html.R` and backed by reproducible
`script/bench_*.R` drivers.

## Documentation

- All user-facing prose updated from "two sources" to "three
  sources" of information
- Rphylopars references replaced with "phylogenetic BM" or
  "internal conditional-MVN baseline" throughout
- `DOCS.md` updated with new benchmarks, datasets, and functions

---

# pigauto 0.4.0

## Breaking: TabPFN baseline removed

`fit_baseline_tabpfn()` and `setup_tabpfn()` are gone, along with the
`reticulate` suggested dependency, the `script/bench_tabpfn.*` files,
`script/run_tabpfn.py`, `tests/testthat/test-tabpfn.R`, and the
`dev/bench_tabpfn.html` benchmark page. The TabPFN wrapper was
originally added as a curiosity — a test of whether a general-purpose
tabular transformer could compete with pigauto as a drop-in baseline.
In practice it is a poor fit for what pigauto is trying to do:

1. **It is not phylogenetic.** TabPFN treats rows as exchangeable and
   cannot use the tree at all. The whole point of pigauto is that the
   tree carries information about trait evolution; feeding traits
   through a tree-blind model throws that away.
2. **Quadratic attention does not scale.** TabPFN's attention is
   O(n²) in the number of rows. At `n = 9,993` species (the full
   AVONET dataset) this is already near the edge of what a single
   GPU can handle, and most users of pigauto care precisely about
   scales at which a baseline must be cheap.
3. **The cross-language subprocess architecture was fragile.** R
   `torch` and Python `torch` share `libtorch` and cannot coexist in
   the same process, so TabPFN had to run via `system2()` on a
   separately-installed Python venv. Every time the user's Python,
   `tabimpute`, or `torch` versions drifted, the wrapper broke. This
   is not maintenance the project can afford.
4. **Keeping the benchmark page gave a misleading impression of
   relevance.** Users who hit the "pigauto vs TabPFN" page came away
   thinking TabPFN was a first-class supported baseline in pigauto.
   It never was.

Users who still want to benchmark pigauto against TabPFN can install
`tabimpute` separately and call it directly from Python. The removal
simplifies pigauto's installation, CI, and documentation.

## Docs cleanup (carried over from v0.3.3 work)

- `DOCS.md`: rewritten benchmark and developer-reports sections to
  reflect the v0.3.2+ benchmark set (scaling to 10,000 tips and the
  AVONET missingness sweep). Dropped the obsolete `benchmark_report`
  and `benchmark_final_report` links that v0.3.2 had already pruned
  from `_pkgdown.yml` but had not yet removed from `DOCS.md`.
- `CLAUDE.md`: removed the TabPFN Python-venv setup recipe and the
  "R torch ↔ Python torch cannot share a process" gotcha (no longer
  applicable). Version line bumped to 0.4.0.

# pigauto 0.3.3

## Correctness fix: Vphy must be a correlation matrix in `glmmTMB`'s `propto()`

The multiple-imputation example in `README.md`, the two HTML
tutorials (`pigauto_intro.html`, `pigauto_workflow_mixed.html`), and
`tests/testthat/test-multi-impute.R` all constructed the phylogenetic
random-effect structure with

```r
Vphy <- ape::vcv(tree)            # WRONG: covariance, diag = tree height
glmmTMB(... + propto(0 + species | dummy, Vphy), ...)
```

`propto()` estimates `sigma^2` freely, so passing the raw covariance
matrix (whose diagonal is the tree height) silently rescales the
phylogenetic variance component by tree height. The correct usage
passes a *correlation* matrix:

```r
Vphy <- cov2cor(ape::vcv(tree))   # diag = 1
```

All four sites are now fixed. Thanks to the user who flagged this.

## Honest framing: the GNN is a gated ensemble, not a residual model

The user-facing prose throughout the package historically described
pigauto as a *Residual Phylogenetic Denoising Autoencoder
(ResidualPhyloDAE) that learns a residual correction on top of a BM
baseline*. This wording was misleading for binary and categorical
traits in two ways:

1. The "baseline" for those types is **phylogenetic label propagation**
   (a similarity-weighted average of observed labels using a Gaussian
   kernel on the cophenetic distance) — not a generalised linear model
   that produces residuals on a link scale.
2. The GNN is trained **end-to-end** by minimising the type-appropriate
   loss (MSE for continuous/count/ordinal, BCE for binary,
   cross-entropy for categorical) on the **blended output**, not on
   `y - baseline`. So the GNN's output `delta` is a full per-cell
   prediction, not a statistical residual.

What pigauto actually is, in honest terms, is a **gated ensemble** of
two predictors:

- A phylogenetic **baseline** (Brownian motion via Rphylopars for
  continuous/count/ordinal traits; phylogenetic label propagation for
  binary/categorical traits).
- A **graph neural network correction** (an internal torch module that
  produces a full per-cell prediction `delta` from spectral node
  features and an attention-biased message passing over the
  phylogeny).

Prediction is the per-trait blend
`(1 - r_cal) * baseline + r_cal * delta`, where `r_cal` is calibrated
on a held-out validation split (with a split-validation cross-check
for discrete traits). When the baseline is already optimal the gate
closes and the GNN becomes a no-op.

The internal torch class name `ResidualPhyloDAE` is preserved because
its GNN layers use ResNet-style residual skip connections — that is a
correct, well-defined ML term. But user-facing prose no longer
describes the GNN as "learning a residual on top of the baseline".

Files updated: `README.md`, `DESCRIPTION`, `CLAUDE.md`, `_pkgdown.yml`,
`DOCS.md`, `vignettes/getting-started.Rmd`, `vignettes/mixed-types.Rmd`,
`R/fit_pigauto.R` (roxygen), `R/predict_pigauto.R` (roxygen),
`R/model_residual_dae.R` (clarifying header comment),
`R/simulate_traits.R`, `R/plot.R` (plot titles), and
`script/make_intro_html.R`.

## Tracking issue for v0.4.0 architectural work

The honest fix above is a documentation correction, not a model
change. A genuinely *residual* baseline for binary and ordinal traits
would require a liability-scale model (threshold BM, `phylolm::phyloglm`,
or `MCMCglmm` with `family = "threshold"`) — a design decision that
should be made deliberately rather than rolled into a same-day patch.
A v0.4.0 tracking issue lays out the trade-offs.

# pigauto 0.3.2

## New bundled dataset: `avonet_full` / `tree_full`

The full AVONET3 + BirdTree Stage2 Hackett MCC phylogeny, aligned to
9,993 species with 7 mixed-type traits (4 continuous morphometric +
2 categorical + 1 ordinal) and native missingness preserved.
Complements the 300-species `avonet300` / `tree300` bundled data:

- `avonet300` / `tree300` remain the right choice for quick examples,
  unit tests, and the getting-started vignette.
- `avonet_full` / `tree_full` are the right choice for scale
  benchmarks, realistic-missingness experiments, and anything that
  wants to exercise the v0.3.1 sparse Lanczos / cophenetic caching
  code paths on a real-world phylogeny.

The schemas are identical, so any code that runs on `avonet300` runs
on `avonet_full` with no modification. Combined tarball increment is
~328 KB (xz-9).

## New benchmark: AVONET missingness sweep (20 / 50 / 80%)

*Articles -> Benchmarks -> AVONET missingness sweep* is a new pkgdown
page comparing three methods on the full `avonet_full` dataset at
three MCAR missingness levels:

- **mean / mode** — column mean (continuous, ordinal) or class
  frequency (categorical). No phylogeny.
- **BM baseline** — Rphylopars Brownian-motion for continuous /
  ordinal, phylogenetic label propagation for discrete traits.
- **pigauto** — full pipeline (BM baseline + calibrated GNN delta
  + conformal intervals).

Per-trait RMSE, Pearson r, Spearman rho, and accuracy on held-out
test cells. Single seed per cell, trait-level masking so categorical
traits are held out as whole groups. Hyperparameters match the
`validate_avonet_full.R` scaling run for comparability.

## Docs cleanup

- Regenerated `pigauto_intro` and `pigauto_workflow_mixed` HTML
  tutorials against the v0.3.2 package state.
- Removed stale v0.3.0-era benchmark pages from the navbar
  (`benchmark_report.html`, `benchmark_final_report.html`,
  `test_report.html`) and from `pkgdown/assets/dev/`.
- Pruned 45 orphan files from `script/` (trial artefacts from
  pre-v0.3.1 benchmarking phases: `bench_v2/v3/v4.*`,
  `benchmark_avonet2000.*`, `benchmark_final.*`,
  `benchmark_simulation.*`, scaling duplicates, and Apr-6 demo
  artefacts).
- `README.md` citation now points at 0.3.2.
- `DESCRIPTION` updated to mention both bundled datasets.

# pigauto 0.3.1

## Scaling to 10,000 tips

`pigauto` now runs end-to-end on real 10,000-tip phylogenies. Before this
release the pipeline was tested and tuned around the 300-species bundled
AVONET data, and at 10,000 tips it wedged on a single dense eigendecomposition
of the graph Laplacian (hours of CPU time). Two targeted fixes land in this
release; no rewrite was needed.

### Fix A: sparse Lanczos eigensolver for `build_phylo_graph()`

`build_phylo_graph()` now uses `RSpectra::eigs_sym()` sparse Lanczos to
compute the `k + 1` smallest Laplacian eigenvectors when the tree has more
than 7,500 tips and **RSpectra** is installed. This replaces
`base::eigen(L, symmetric = TRUE)`, which was O(n^3) in time and discarded
>99% of its work (we only need the bottom ~32 eigenvectors).

- On the real 9,993-species AVONET + BirdTree Stage2 Hackett MCC phylogeny,
  the graph stage now takes under a minute, compared with roughly seven
  minutes of dense eigen() on the same hardware.
- The old dense path is kept as the default for smaller trees and as a
  safety fallback for when RSpectra is unavailable or Lanczos does not
  converge. The threshold is conservative (7,500) because highly symmetric
  ultrametric simulations (`ape::rcoal`) have degenerate spectra that need
  large Krylov subspaces to resolve.
- **RSpectra** is listed in `Suggests`. Installing it is strongly recommended
  for any tree larger than ~7,500 tips.

### Fix B: cache the cophenetic distance matrix across the pipeline

The cophenetic distance matrix `D = ape::cophenetic.phylo(tree)` was previously
computed four separate times in a single `impute()` call (three times inside
`build_phylo_graph()` and once inside `fit_baseline()`). Each call allocates a
dense n*n double matrix and runs an O(n^2) tree traversal -- wasted work at
scale. In 0.3.1:

- `build_phylo_graph()` computes `D` exactly once and returns it as part of
  the result list (`graph$D`).
- `fit_baseline()` gains an optional `graph` argument; when supplied, it
  reuses `graph$D` instead of recomputing the cophenetic matrix.
- `impute()` and `fit_pigauto()` pass the graph through automatically, so
  users get the speedup without any code changes.
- `impute()` also explicitly releases `graph$D` after `fit_baseline()`
  finishes, so the ~800 MB cophenetic matrix does not sit in R memory
  during training (a subtle side-effect that caused a large per-epoch
  regression at n = 10,000 until it was dropped).

### Scaling benchmark

See `pkgdown/assets/dev/scaling.html` (linked from the site as
*Articles -> Benchmarks -> Scaling to 10,000 tips*) for a side-by-side
comparison of v0.3.0 vs v0.3.1 at n in {300, 1k, 2k, 3k, 5k, 7.5k, 10k},
plus an end-to-end validation run on the full AVONET3 + BirdTree dataset
(9,993 species, 7 mixed-type traits). The pkgdown article also walks
through both a quick 300-species example with the bundled `avonet300`
data and a full 10k-species example using the BirdTree Stage2 Hackett
MCC phylogeny.

## Bug fixes

- `build_phylo_graph()`: the N > 2000 memory warning was out of date (the
  real blocker above that size was wall-clock time, not memory). The
  warning is now issued at N > 10000 and its text points at the true
  bottleneck.

# pigauto 0.3.0

## Major new features

### Attention-based message passing
The GNN now uses learned attention weights (with phylogenetic adjacency
as positional bias) instead of fixed adjacency matrix multiplication.
The model learns which phylogenetic neighbours are most informative for
each species, while the log-adjacency bias ensures initialisation
respects evolutionary distances. Controlled via `use_attention = TRUE`
(now the default).

### Validation-calibrated gates
After training, per-trait gate values are optimised on the validation set
via grid search. This prevents the GNN from adding noise when the
Brownian motion baseline is already optimal. On the AVONET benchmark,
this eliminates the RMSE regression seen in earlier versions.

### Conformal prediction intervals
Distribution-free 95% prediction intervals computed via split conformal
prediction on the validation residuals. Coverage is guaranteed under
exchangeability. Available for continuous, count, and ordinal traits via
`pred$conformal_lower` and `pred$conformal_upper`.

### Phylogenetic label propagation
Binary and categorical trait baselines now use phylogenetic similarity-
weighted label propagation instead of flat population-level frequencies.
Each species gets a personalised baseline reflecting its phylogenetic
neighbourhood.

### Adaptive spectral encoding
`k_eigen` now defaults to `"auto"`, scaling with tree size:
`min(max(ceil(n/20), 4), 32)`. This gives 4 features for small trees,
15 for 300 tips, and 32 for 640+ tips.

### Auto-generated HTML reports
`pigauto_report()` produces self-contained HTML reports with interactive
Chart.js visualisations: per-trait metrics, calibrated gate values,
training history, and conformal coverage.

### Evaluation framework
- `evaluate()`: per-trait metrics on held-out data
- `compare_methods()`: automated BM vs GNN comparison
- `cross_validate()`: proper k-fold cross-validation
- `summary.pigauto_fit()`: formatted console summary

### Plotting
- `plot.pigauto_fit()`: training history, gate values, conformal scores
- `plot.pigauto_pred()`: scatter plots, conformal intervals, probability
  distributions
- `plot_comparison()`: forest-plot style method comparison

### Simulation benchmark
`simulate_benchmark()` runs a complete simulation study with configurable
scenarios (BM, OU, regime shift, non-linear, mixed types), returning tidy
results with print/summary/plot methods. Useful for methods papers and
for verifying pigauto's behaviour on data with known properties.

## Minor improvements

- Categorical gate calibration now works at the trait level (all K one-hot
  columns share a single gate) using cross-entropy loss, not per-column MSE
- `build_phylo_graph()` default `k_eigen` changed from `8L` to `"auto"`
- `fit_pigauto()` default `k_eigen` changed from `8L` to `"auto"`
- `all_grads_finite()` now handles non-standard tensor types gracefully
- `impute()` uses adaptive k_eigen by default
- Package passes `R CMD check` with 0 errors and 0 notes
- Bumped version to 0.3.0

# pigauto 0.2.0

- Added support for multiple observations per species (`species_col`)
- Added mixed trait type support (continuous, binary, categorical,
  ordinal, count)
- Added AVONET 300 bundled dataset with mixed types
- Extended simulation benchmark (16 scenarios via BACE)

# pigauto 0.1.0

- Initial release
- ResidualPhyloDAE architecture
- Brownian motion baseline via Rphylopars
- Per-column gated blending
- MC dropout for multiple imputation
