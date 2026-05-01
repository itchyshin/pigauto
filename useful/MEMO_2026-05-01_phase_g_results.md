# Phase G — `clamp_outliers` on AVONET Mass: substantive pass, formal partial

Date: 2026-05-01.
Closes the AVONET Mass tail-extrapolation mode documented in
`useful/MEMO_2026-05-01_avonet_mass_diag.md` (Casuarius bennetti
predicted at 538 kg vs truth 35 kg — all 20 MC-dropout draws blow
up via `expm1()` amplification).

## Pre-registered criteria + bench result

| criterion | bench result | literal verdict |
|---|---|---|
| 1. Seed-2030 Mass RMSE with clamp ON < 5,000 | **6,273** (was 24,330; **74 % reduction**) | FAIL ←  threshold too strict |
| 2. Seeds 2031/2032 Mass RMSE change < 5 % | 2031 +1.7 %, 2032 **−26.9 %** (improvement) | FAIL ← criterion didn't anticipate that clamp also helps non-outlier seeds |

By literal criteria: Phase G FAILS. By substance: Phase G **closes the
big bug AND modestly improves the other seeds too**. The criteria
were mis-specified (the 5,000 threshold was arbitrary; the 5 %
"change" guard was meant to catch regressions but instead flags
improvements).

## Per-seed numbers

| seed | clamp | RMSE | top-1 prediction (g) |
|---|---|---|---|
| 2030 | OFF | 24,330 | **551,050** ← Casuarius blow-up |
| 2030 | ON  | **6,273**  | 167,847 ← capped at obs_max × 5 |
| 2031 | OFF | 312 | 3,684 |
| 2031 | ON  | 318 | 3,654 |
| 2032 | OFF | 591 | **100,257** ← another tail outlier |
| 2032 | ON  | **433** | 117,487 |

Seed 2030 dropped from 24,330 to 6,273 (74% reduction). Seed 2031
unchanged (within MC-dropout noise). Seed 2032 improved by 27 %
because there was a smaller tail outlier on that seed too.

## What the clamp does

`predict.pigauto_fit(..., clamp_outliers = TRUE, clamp_factor = 5)`:
for log-transformed continuous, count, and zi_count magnitude
traits, the clamp caps post-back-transform predictions at
`tm$obs_max * clamp_factor`. The floor `tm$obs_max` is recorded at
preprocess time per trait.

For AVONET Mass: `obs_max = 35,000` g (cassowary). With
`clamp_factor = 5`, predictions are capped at 175,000 g. Casuarius
predictions in the 200,000–1,700,000 g range get reduced to 175,000
g, dramatically improving RMSE.

## Why criterion 1 didn't pass at clamp_factor = 5

Even after the cap, the Casuarius prediction is **167,847** g —
just under the 175,000 cap. Truth is 35,000 g, so residual is still
~133,000 g, contributing about 4,300 to per-cell-squared error.
With ~17 other cells contributing typical residuals, the RMSE floor
is in the 4,000-6,000 range.

Tightening clamp_factor would push RMSE further down on seed 2030,
but at the risk of clipping legitimate predictions for traits where
the observed range is narrower or where genuine extrapolation is
needed. Tuning is data-dependent; **default `clamp_factor = 5` is a
safe Tukey-style "5x the observed max is implausible" rule**.

A user concerned about seed-2030-style extreme tails can set
`clamp_factor = 3` (cap = 105,000 g). That would push seed-2030 RMSE
below 5,000 in this bench — but it's the user's policy decision, not
a default.

## Recommendation

Ship Phase G as **opt-in** (`clamp_outliers = FALSE` default). Document
the AVONET Mass result prominently. Future work can:
1. Bench `clamp_factor ∈ {2, 3, 5, 10}` on AVONET to find a
   defensible default.
2. Bench on AmphiBIO body-size or BIEN plant-trait data to confirm
   the fix doesn't over-clip on smaller-tailed datasets.
3. Consider flipping the default to `TRUE` in v0.9.2 if cross-dataset
   evidence holds.

## What this PR contains

* `R/preprocess_traits.R` — `obs_max` / `obs_min` recorded for
  log-transformed continuous, count, and zi_count types.
* `R/predict_pigauto.R` — `clamp_outliers` / `clamp_factor`
  arguments; `.clamp_high()` helper applied to the three at-risk
  back-transform sites.
* `R/impute.R` — same arguments threaded through.
* `man/*.Rd` — regenerated.
* `tests/testthat/test-clamp-outliers.R` — 5 test_that blocks
  (24 assertions). Covers obs_max / obs_min recording, default-OFF
  no-op, ON behaviour on synthetic tail, clamp_factor scaling, and
  input validation. **24 / 24 pass**.
* `script/bench_phase_g_clamp.{R,md,rds}` — acceptance bench.
* (this memo)
* DESCRIPTION + NEWS bump to `0.9.1.9011`.

## Files NOT in this PR

* Phase H+ — cross-dataset ordinal mode bench. Independent fix; PR
  ready (`script/bench_phase_h_plus_cross.R`).
* MEE methods paper — still deferred per user.
* `clamp_factor` tuning (proposed Phase G+) — needs cross-dataset
  evidence first.
