# AVONET Mass instability diagnosis — seed 2030, N_IMP = 20

Date: 2026-05-01.
Triggers `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md`
finding "Mass on AVONET is currently unstable: at N_IMP=20 the BM joint
baseline emits singular-Σ warnings on Mass, and one of three seeds
produces a tail outlier with RMSE ≈ 11× the column-mean baseline".

This memo documents the diagnostic run that localised the failure mode.

## What was run

`script/diag_avonet_mass.R` reproduced the seed-2030 AVONET n=1500
N_IMP=20 setup exactly (same subset, same mask, same impute() call)
and dumped per-imputation Mass predictions on the 450 held-out cells
to `script/diag_avonet_mass.rds`.

Wall: ~7 min.

## Headline numbers

| metric | value |
|---|---|
| Per-imputation Mass RMSE — min   | 10,156 |
| Per-imputation Mass RMSE — q25   | 15,250 |
| Per-imputation Mass RMSE — median | 23,728 |
| Per-imputation Mass RMSE — q75   | 29,743 |
| Per-imputation Mass RMSE — max   | 79,844 |
| **Median-pooled Mass RMSE**      | **23,723** |
| Mean-pooled Mass RMSE            | 26,517 |
| Per-imputation max / median ratio | 3.36 |
| Imputations with RMSE > 5× median | 0 / 20 |

## Top-5 worst pooled-prediction residuals

| species | truth (g) | pooled median (g) | residual | per-draw range |
|---|---|---|---|---|
| **Casuarius bennetti**       | 35,000.0 | **538,166** | −503,166 | 250,202 to 1,728,734 |
| Plectropterus gambensis      | 3,869.0  | 7,860       | −3,991   | 3,711 to 18,877 |
| Ardea insignis               | 2,024.0  | 5,488       | −3,464   | 2,317 to 17,381 |
| Gyps himalayensis            | 9,798.0  | 6,625       | +3,173   | 4,100 to 9,391 |
| Platalea flavipes            | 1,741.3  | 3,869       | −2,128   | 1,469 to 17,515 |

## Verdict

**It is NOT a single rogue MC-dropout draw.** All 20 draws for
Casuarius bennetti are in the 250 kg–1.7 *megagram* (1,700 kg) range;
the per-imputation RMSE max is only 3.36× the median, and **0 of 20**
imputations have RMSE > 5× the median. So median pooling will not
rescue this — the model is consistently extrapolating Casuarius into
implausible territory, every single draw.

This is **tail-extrapolation amplified by the `expm1()` back-transform**.
On the latent log-z scale, the prediction sits ~3–4 SDs above any
observed Casuariiformes value. After `expm1()` on the un-z'd log
scale, that 3–4 SD overshoot becomes a 7×-50× value error in original
grams.

Casuarius bennetti is a cassowary — flightless, plate-2 ratite,
isolated phylogenetic position (only ~3 cassowary species in the
tree). With ~30 % cell mask, Casuarius's Mass is a held-out cell and
the BM joint baseline has nothing close on the tree to constrain
it. The GNN delta on top adds further noise. Result: every draw
hallucinates a Mass that's 7×–50× too large.

The other top-5 species (Plectropterus, Ardea, Gyps, Platalea) are
moderately phylogenetically isolated large-bodied birds; their
residuals are in the thousands-of-grams range — meaningful but not
catastrophic. The pooled RMSE of 23,723 is mostly driven by
Casuarius alone.

## Implications for users

This is **not** a bug introduced by Fix A-H, the safety-floor work,
or the multi-seed re-verify. It's a long-standing tail-extrapolation
behaviour that only became visible at N_IMP = 20 because the noise
envelope tightened enough to expose the consistent over-prediction.

The README "Caveats from multi-seed evidence" added in PR #51 already
warns users that "If you are imputing Mass, run multiple seeds and
inspect the distribution of imputed values for outliers". This
diagnosis confirms the warning is the right one — multi-seed alone
won't fix it for Casuarius, but it will at least flag the species as
an outlier.

## Proposed fix (Phase G — separate from Phase F)

A sensible value clamp would catch the worst offender:

```r
# In predict_pigauto.R, after expm1() back-transform for log-transformed
# continuous traits:
obs_max <- max(observed_values, na.rm = TRUE)
obs_min <- min(observed_values, na.rm = TRUE)
clamp_factor <- 5  # tunable; 5× the max plausible sticks Casuarius
                   # at ~175 kg cap rather than 538 kg
pred[pred > obs_max * clamp_factor] <- obs_max * clamp_factor
pred[pred < obs_min / clamp_factor] <- obs_min / clamp_factor
```

This is opt-in (`clamp_outliers = TRUE` argument on `impute()` and
`predict()`), backward-compatible (default `FALSE`), and trivially
reversible if it harms a benchmark.

**Pre-registered acceptance criterion for Phase G**:
- AVONET seed-2030 Mass RMSE drops below 5,000 (vs current 23,723)
  with `clamp_factor = 5`.
- No regression on AVONET seeds 2031, 2032 (which already beat
  baseline) — clamped predictions on those seeds should change by
  < 5 % of pooled RMSE.
- Tested on 3 simulated continuous traits with known tails to verify
  the clamp doesn't bias predictions inside the observed range.

Tracked separately; not part of Phase F.

## Files

* `script/diag_avonet_mass.R` — the diagnostic driver
* `script/diag_avonet_mass.rds` — per-imputation predictions for 450 Mass cells × 20 draws
* (this memo)
* upstream: `script/bench_avonet_full_local_n20_seed2030.{rds,md}` — the original bench
* upstream: `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md` — the multi-seed evidence that prompted this diagnosis
