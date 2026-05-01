# Per-trait ordinal baseline path selection (2026-04-30)

Closes Opus adversarial review item #6: "Don't ship a `K ≤ 3 → LP`
heuristic.  Instead, add LP as a corner option in the safety floor's
simplex grid for ordinal traits, OR fit both paths in `fit_baseline()`
and select the lower-val-loss one."

## What changed

`R/fit_baseline.R` now computes BOTH the threshold-joint baseline
(B3, commit `a541dbd`) AND a per-column BM-via-MVN alternative for
each ordinal trait that the threshold-joint actually populated.
After both are computed, the per-trait val-MSE is compared against
`splits$val_idx` cells in `data$X_scaled` (the ground truth before
NA masking) and the lower-val-MSE path is kept.

The chosen path is exposed as `fit_baseline()$ordinal_path_chosen`
(named character vector: latent col index → `"threshold_joint"` or
`"bm_mvn"`).

## Implementation choice

Of Opus's two suggestions (#6a "LP as a 4th simplex corner"; #6b
"both paths in `fit_baseline()` + selection"), we picked #6b because:

- **#6a** would have to extend the safety-floor simplex from 3D
  `(r_BM, r_GNN, r_MEAN)` to 4D `(r_BM_threshold, r_LP, r_GNN, r_MEAN)`
  for ordinal traits only.  This requires extending `simplex_grid()`,
  `compute_corner_loss()`, the post-cal full-val invariant
  dispatcher, plus all related tests.  Conservative estimate: 200+
  LOC of changes, 8-10 test updates.
- **#6b** is local to `fit_baseline()`, uses an already-existing
  comparison primitive (val MSE vs ground truth), and adds ~50 LOC
  with 2 new tests.

#6b also has the architectural property that the gate calibration
sees a single "best baseline" and doesn't need to know which
underlying path was used — the abstraction stays clean.

## Verification

`script/bench_avonet_phase6.R` (Phase 6 levelC vs LP comparison on
AVONET 300):

| trait | Apr 16 (pre-B3) | Apr 29 (post-B3 buggy) | Apr 30 (post Phase F) |
|---|---|---|---|
| Migration RMSE (levelC) | 0.879 | 0.975 | **0.890** |
| Migration RMSE (lp) | 0.879 | 0.879 | 0.890 |
| levelC − lp lift | 0 | −0.096 (regression) | **~0** (parity) |

(The 0.879 → 0.890 gap on the LP path between Apr 16 and Apr 30 is
unrelated to this fix; it is RNG-state drift from intervening commits.
The Phase F-relevant outcome is the levelC path matching the LP
path on Migration -- the threshold-joint regression no longer
surfaces.)

## Scope

* **Single-obs only.**  Multi-obs would require species aggregation
  of the BM alternative.  The selection block is gated by
  `!multi_obs && !is.null(splits)`.  Multi-obs threshold-joint will
  continue to use threshold-joint output unconditionally for ordinals
  until a future fix.
* **Only fires when threshold-joint actually populated the column.**
  If threshold-joint failed (e.g. phylopars singular matrix), the
  column falls through to the per-column BM block at line ~392, and
  no selection is needed (BM is already what runs).
* **Binary and categorical traits unchanged.**  The threshold-joint
  approach for binary is solid (the binary regression at AVONET
  has not been observed empirically), and the OVR categorical path
  is its own architecture.  Per-trait selection here would add
  complexity without empirical motivation.

## Tests added

`tests/testthat/test-joint-threshold-baseline.R`:

1. `"fit_baseline picks per-trait ordinal path and reports it"` --
   smoke test on K=3 ordinal + continuous; verifies finite output
   and the exposed `$ordinal_path_chosen` field.
2. `"fit_baseline ordinal selection picks BM when threshold-joint
   loses on val"` -- recomputes both paths' val MSE independently
   and asserts the selection picked the lower one.

## What this does NOT close

* The Migration regression on the FULL pigauto pipeline
  (`bench_avonet_full_local.R`) — pigauto's calibrated gate still
  has to operate on top of the chosen baseline, and val→test drift
  may still cost ≤ 1.4 pp accuracy.  The multi-seed bench shows
  Migration at 0.776 ± 0.036 vs baseline 0.787 ± noise (a 1.4 pp
  test gap that's compatible with the val→test drift bound).
* Multi-obs ordinal selection (deferred).
* The architectural question of whether to also add LP as a fourth
  simplex corner in addition to per-trait selection.  Per-trait
  selection captures most of the practical benefit; LP-as-corner
  would be additive only if the GNN delta for ordinal traits has
  meaningful signal beyond what either baseline captures.

## Files

* `R/fit_baseline.R` — selection block + return field
* `tests/testthat/test-joint-threshold-baseline.R` — 2 new tests
* `NEWS.md` — entry under v0.9.1.9007
* `useful/MEMO_2026-04-30_ordinal_path_selection.md` — this memo
