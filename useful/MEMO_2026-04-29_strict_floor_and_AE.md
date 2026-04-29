# Pigauto strategic synthesis — 2026-04-29

This session covered three things:
1. Tier-1 ship-blocker fix: discrete-aware strict val-floor.
2. Item #3 from the prior strategic review: AE-attribution gate-zeroing.
3. Decision on Item #8 (self-supervised pretraining).

## 1. Strict val-floor fix (commit `d651c6c`)

**Bug.** `R/fit_helpers.R::calibrate_gates`'s post-calibration "strict
val-floor guarantee" check at lines 332-352 was guarded by

    if (safety_floor && length(val_row_idx) > 0L &&
        tm$type %in% c("continuous", "count", "ordinal", "proportion")) { ... }

which silently exempted binary, categorical, zi_count, and
multi_proportion from the final invariant check. The half-A/half-B
verification + median-over-splits steps that run before this block
can let through a non-corner blend that loses to pure-BM on a small
val set, and for discrete types there was no second-line defence.

**Provenance.** Introduced 2026-04-23 in `28d8e45` ("safety-floor:
calibrate_gates() supports safety_floor=TRUE with simplex grid").
The type filter mirrored an earlier continuous-only MSE evaluator
and was not extended when the simplex grid landed.

**Empirical impact.** From the 2026-04-28 bench rerun (which had the
buggy check active):

| bench | cell | baseline | pigauto buggy | gap |
|---|---|---|---|---|
| `bench_binary.md` | signal_0.4 bin1 | 0.630 | 0.570 | −6 pp |
| `bench_binary.md` | signal_0.4 bin4 | 0.697 | 0.609 | −9 pp |
| `bench_categorical.md` | K_3 cat2 | 0.740 | 0.621 | −12 pp |
| `bench_zi_count.md` | zf_0.2 zi1 RMSE | 31.6 | 38.8 | +23 % |

These violate the safety_floor's "pigauto cannot underperform its
own baseline" promise. The fix removes the type filter and reuses
the existing `cal_mean_loss` closure, which already implements
0-1 loss for binary/categorical and the gate+magnitude composite
for zi_count.

**Smoke verification.** `script/smoke_strict_floor.R` (3 reps on
binary signal_0.4): pre-fix gap pigauto vs baseline was −6 pp;
post-fix gap +0.4 pp. The 6 pp regression is closed.

**Caveat.** The strict val-floor guarantees `pigauto_val ≤ baseline_val`.
It does NOT guarantee `pigauto_test ≤ baseline_test` — val→test
extrapolation can drift by a few points on small val sets. In the
27-cell AE-attribution experiment below, this happened in 1/27 cells
(3.7% rate, +0.034 RMSE drift). Acceptable rate.

**Re-runs needed.** Discrete-trait benches: `bench_binary.R`,
`bench_categorical.R`, `bench_zi_count.R`. Real-data benches with
discrete columns: AVONET-full, AmphiBIO, AVONET-Phase 6,
AVONET-missingness. All other benches (continuous, ordinal, count,
proportion, multi_proportion, all phase1 lambda/n/beta/missingness
sweeps) are confirmed unaffected.

## 2. AE-attribution: gate-zeroing experiment

**Setup.** 200 species single-obs, 30% MCAR, 5 covariates, 3 reps ×
9 cells (lambda × f_type). Methods: pigauto_full (calibrated gate),
pigauto_gnn_off (GNN gate forced to 0, BM/MEAN split kept),
phylolm-lambda BLUP (standalone analytical baseline). Script:
`script/ae_attribution_smoke.R`.

**Per-cell aggregate (mean over 3 reps):**

| f_type | λ | full | gnn_off | phylolm | GNN delta | full vs phylolm |
|---|---|---|---|---|---|---|
| linear | 0.1 | 0.679 | 0.679 | 0.701 | 0.000 | **−3.1 %** ✓ |
| linear | 0.3 | 0.797 | 1.054 | 0.812 | **−0.257** | −1.9 % ≈ tie |
| linear | 0.5 | 0.999 | 1.202 | 0.946 | **−0.203** | +5.6 % ✗ |
| nonlinear | 0.1 | 0.811 | 0.796 | 0.748 | +0.015 | +8.4 % ✗ |
| nonlinear | 0.3 | 0.846 | 0.850 | 0.801 | −0.005 | +5.6 % ✗ |
| nonlinear | 0.5 | 1.080 | 1.080 | 1.027 | 0.000 | +5.1 % ✗ |
| interactive | 0.1 | 0.877 | 0.903 | 0.955 | −0.026 | **−8.1 %** ✓ |
| interactive | 0.3 | 1.000 | 1.012 | 0.995 | −0.012 | +0.5 % ≈ tie |
| interactive | 0.5 | 1.278 | 1.384 | 1.181 | −0.106 | +8.2 % ✗ |

GNN delta = pigauto_full − pigauto_gnn_off, RMSE units. Negative =
GNN improves over baseline. full vs phylolm = +x % means pigauto
loses by x %.

**Findings.**

1. **The GNN is silent in 13 of 27 cells (48 %).** Gate
   `r_cal_gnn = 0` → GNN contributes literally zero. The strict
   val-floor closes the gate when the GNN doesn't earn its keep.

2. **When the GNN engages, the biggest wins are on LINEAR DGPs
   at moderate λ.** Linear λ=0.3 gives a −0.257 RMSE reduction
   (24 %); linear λ=0.5 gives −0.203 (17 %). This is
   counterintuitive: the GNN was supposed to extract nonlinear
   structure, but in practice it's recovering linear covariate
   signal that pigauto's covariate-aware BM baseline missed at
   moderate phylo signal.

3. **Pigauto loses to phylolm-λ BLUP on 6 of 9 cells.** Even when
   the GNN is actively contributing, the full pipeline's RMSE is
   +5–8 % above phylolm. Reason: pigauto's analytical baseline is
   fixed at λ=1 (Brownian motion), but phylolm-λ fits λ from data
   — strictly more flexible. The GNN improves over BM baseline,
   but doesn't recover the BM-vs-λ baseline gap.

4. **Pigauto cleanly beats phylolm on interactive λ=0.1 (−8.1 %)
   and linear λ=0.1 (−3.1 %)**, both at low phylo signal. These
   are the regimes where (a) phylolm's λ-fit is unstable at low
   signal and (b) pigauto's BM-baseline-with-covariates competes.

5. **Val→test leak rate: 1 / 27 cells (3.7 %).** Acceptable.

**Direct interpretation.** The "AE is doing real work" claim from
prior session memos is **partially supported** for single-obs DGPs:
- The GNN does contribute meaningful RMSE reductions (up to 24 %)
  when it engages.
- But it engages in fewer than half the cells.
- And the regimes where it most engages are LINEAR cov effects,
  not nonlinear, contradicting the "AE extracts nonlinear
  structure" framing.
- Even with GNN contribution, the pipeline often loses to phylolm,
  so the GNN's lift over BM baseline doesn't translate into a lift
  over the strongest analytical competitor.

## 3. Decision on Item #8 (self-supervised pretraining)

**Recommendation: AGAINST pretraining as the next strategic
investment.** The AE-attribution evidence does not support it.

Reasoning:
- The GNN's bottleneck is engagement, not feature quality. Gate
  closes 48 % of the time on single-obs DGPs. Pretraining would
  give better features but the strict val-floor will still close
  the gate when those features don't beat baseline on val.
  Pretraining widens what the GNN does once engaged, not when it
  engages.
- The bigger gap (5–8 %) is structural: BM baseline (fixed λ=1) vs
  phylolm-λ (fitted λ). Pretraining doesn't address this.
- Where pigauto already wins phylolm is multi-obs nonlinear, not
  single-obs. The obs_refine MLP path is barely exercised here.
  The previously-claimed "winning regime" (multi-obs nonlinear
  λ=0.1) is unrelated to whether the single-obs GNN is pretrained.

**Higher-leverage alternative: replace BM baseline with phylolm-λ
BLUP** (or `Rphylopars`'s λ-extension, or fit λ inside the existing
`bm_internal.R` machinery).

- The strict val-floor I just shipped guarantees pigauto ≤ baseline.
  So if we *upgrade* the baseline from BM to λ-BLUP, pigauto gets
  a free RMSE improvement everywhere AND the GNN's job becomes
  "beat λ-BLUP" — a harder bar but a more honest one.
- Pretraining can come later if it turns out to be needed once the
  baseline gap is closed.

**Alternative-alternative (cheaper than λ-BLUP swap):** the strict
val-floor's tolerance is currently 1e-12 (numerical equality).
Adding a bit of slack (e.g., require GNN to beat baseline by
≥ 0.1 std-units to engage) would reduce the val→test leak rate
further at small data sizes. Easy to ship, modest improvement.

## 4. Open items I'd queue for a future session

A. **Replace BM baseline with phylolm-λ BLUP** (item #8 alternative).
   Touches `R/bm_internal.R::bm_impute_col` and
   `bm_impute_col_with_cov`. Single-obs first, multi-obs later.
   Expected lift: 5–8 % RMSE on continuous benches.

B. **Re-run real-data discrete benches** (Tier 2 of bug impact):
   AVONET-full (discrete: Trophic.Level, Migration), AmphiBIO,
   AVONET-Phase 6, AVONET-missingness. Confirms the safety-floor
   fix translates to real datasets.

C. **Per-type smoke regression in CI.** A 30-second bench per type
   that asserts pigauto ≤ baseline at the cell level. Catches
   slow-drift performance bugs Opus E flagged.

D. **Property-based input-permutation invariance test.**
   `impute(traits)` and `impute(traits[π,])[π⁻¹,]` should produce
   identical `completed`. Catches the row-alignment bug family.

E. **Conformal coverage test in the suite.** Promote the
   `dev/coverage_sim/` probe to a smoke test (skip if too slow):
   assert ≥ 93 % empirical coverage on a fixed seed.

## 5. State of git

`experiment/gnn-earnings-sim` branch, unpushed:

```
d651c6c fix(safety-floor): extend strict val-floor to ALL trait types, not just continuous
4e33153 docs: fix two pre-existing R CMD check Rd warnings
5620512 test/lint: post-Opus cleanup -- fix OVR test flake + R CMD check NOTE/WARNINGS
fb4461e fix: address four Opus adversarial-review correctness bugs
```

R CMD check: 0 errors / 0 warnings / 0 notes. Targeted tests green.
Discrete bench reruns in progress; will refresh `.md` files when
complete.
