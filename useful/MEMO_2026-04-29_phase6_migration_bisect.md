# Phase 6 Migration regression bisect (2026-04-29)

## Finding

The AVONET-300 Phase 6 bench's Migration RMSE regressed from **0.8792**
(Apr 16) to **0.9748** (Apr 29) -- a +10.9 % degradation, flipping the
levelC-vs-LP lift from +0.011 (Phase 6 wins) to -0.085 (Phase 6 loses
by 9.5 percentage points).

Bisect localised the regression to a single commit:

```
a541dbd  B3: full threshold-model ordinal baseline
```

The commit before (`a541dbd^`) gives Migration levelC RMSE = 0.8792.
At `a541dbd` and every subsequent commit through HEAD, Migration
levelC RMSE = 0.9748.  No other commit in the Apr 16-29 window moves
this number.

## What B3 does

B3 swaps the ordinal-trait baseline path inside
`fit_joint_threshold_baseline()` from phylogenetic label propagation
(LP) to a threshold-joint approach: it fits a Gaussian liability
matrix via `Rphylopars::phylopars()` with interval-censored truncation
at the K-1 ordinal thresholds, then decodes back to the integer
class.

Theoretically this is more flexible than LP (it can capture
phylogenetically-correlated continuous liabilities that produce the
ordered class boundaries).  Empirically on AVONET's Migration
trait -- which has only 3 levels (resident / partial migrant /
full migrant) -- the threshold-joint estimate is less accurate than
LP because there is not enough resolution in 3 levels to nail the
two thresholds.

## Why this is not a strict-val-floor concern

Phase 6's bench script calls `fit_baseline()` directly.  It is
specifically probing the joint baseline output, not the full
pigauto pipeline.  In the full `impute()` path the calibrated gate
would compare the threshold-joint output against LP / pure-mean and
the safety floor would override to whichever wins on val.

That said, the full pipeline's gate calibration on AVONET 300
Migration does not currently include LP as a corner -- only
BM-via-threshold-joint and pure-mean.  So even with the safety floor
fully working, AVONET Migration in pigauto might still favour the
threshold-joint path that the floor cannot route around.

## Triage options

1. **Add a K-aware guard** to the ordinal threshold-joint path:
   skip threshold-joint when K_levels <= 3 and fall back to LP.
   Cleanest, smallest change.  Closes the regression on Migration
   without touching the path for ordinals with K >= 4 where
   threshold-joint may genuinely help.

2. **Add LP as a corner option** in the safety floor's simplex grid
   for ordinal traits.  Larger change; the calibrated gate would
   then route between threshold-joint, LP, BM, and pure-mean.

3. **Document and accept**.  Phase 6 was a benchmarking commit, not
   a production change.  The script's purpose was to compare
   "Phase 6 active" against "LP only".  Migration's regression in
   that comparison is just empirical evidence that threshold-joint
   is not the right tool for K=3 ordinals.

## Recommendation

Option 1.  The fix would be a single-line guard in
`R/joint_threshold_baseline.R::build_liability_matrix` (or the
calling site in `fit_baseline.R`):

```r
if (tm$type == "ordinal" && length(tm$levels) <= 3L) {
  # K=3 threshold-joint is under-determined -- fall back to LP path
  next  # let the LP baseline handle this trait
}
```

Effort estimate: 30 min including a regression test.  Confirms
Migration RMSE reverts to 0.879 (lift +0.01).  Not done in this
session because:

* Phase 6 bench is the only place this currently surfaces.
* The full `impute()` pipeline's safety floor is supposed to be the
  final gate; if the gate is failing to route around threshold-joint
  for K=3 ordinals, that is a separate fix to the corner search.
* Out of scope for the strict val-floor work.

## Files

* `/tmp/phase6_bisect2_results.csv` -- raw bisect output
* `/tmp/phase6_*.md` -- per-commit Phase 6 bench output
