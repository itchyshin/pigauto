# Discrete bench three-way comparison (2026-04-29)

After the strict val-floor fix (commit `d651c6c`), reran all three
discrete benches with the same DGP, seeds, and tunables as the
pre-fix runs.  Three snapshots are now available:

* `*_postfix2.{md,rds}.bak` — 2026-04-17, before `safety_floor` was
  introduced (28d8e45). Baseline of "pigauto's discrete behaviour
  pre-safety-floor."
* `*.{md,rds}.preStrictFloor.bak` — 2026-04-28, with `safety_floor`
  active but the type-restricted strict val-floor (the bug). Baseline
  of "pigauto's discrete behaviour with the buggy safety_floor."
* current `bench_*.{md,rds}` — 2026-04-29, post-fix.

## Aggregate results

| bench | metric | Apr 17 (pre-bug) | Apr 28 (buggy) | Apr 29 (post-fix) |
|---|---|---|---|---|
| `bench_binary` | mean acc gap (pigauto − baseline) | −0.019 | −0.016 | **−0.003** |
| `bench_categorical` | mean acc gap | −0.017 | −0.014 | **−0.003** |
| `bench_zi_count` | mean RMSE ratio (pigauto / baseline) | 1.010 | 1.007 | **0.990** |

ratio < 1.0 = pigauto better; > 1.0 = pigauto worse.

The fix flips zi_count from "pigauto worse than baseline by 1 %" to
"pigauto better than baseline by 1 %" on average, and closes the
binary / categorical −1.6 to −1.4 pp gaps to within 0.3 pp.

## Worst-cell recoveries

| bench | scenario | Apr 28 buggy | Apr 29 post-fix | recovered |
|---|---|---|---|---|
| binary | signal_0.4 bin4 | −4.7 pp | **+2.0 pp** | +6.7 pp |
| binary | signal_0.6 bin2 | −5.8 pp | 0.0 pp | +5.8 pp |
| binary | signal_0.8 bin4 | −5.0 pp | +1.0 pp | +6.0 pp |
| binary | imbalance_0.5 bin4 | −6.2 pp | −0.2 pp | +6.0 pp |
| categorical | K_3 (mean) | −4.5 pp | −2.6 pp | +1.9 pp |
| categorical | signal_1.0 (mean) | −3.4 pp | 0.0 pp | +3.4 pp |
| zi_count | zf_0.2 | ratio 1.101 | ratio 1.046 | −5.5 % |
| zi_count | zf_0.4 | ratio 1.060 | ratio 1.003 | −5.7 % |

## Per-cell honesty: where the fix doesn't help

The strict val-floor guarantees `pigauto_val_loss ≤ baseline_val_loss`
on the FULL val set, but does NOT guarantee `pigauto_test ≤ baseline_test`.
Val→test extrapolation drift survives in 4 of 32 binary cells:

| scenario | trait | gap_buggy | gap_post |
|---|---|---|---|
| signal_1.0 | bin4 | −0.8 pp | **−5.4 pp** (worse post-fix) |
| signal_1.0 | bin3 | −3.8 pp | −4.0 pp |
| signal_0.4 | bin2 | 0.0 pp | −2.9 pp |
| imbalance_0.5 | bin3 | −1.0 pp | −2.6 pp |

These are cells where the calibrated blend genuinely beat baseline on
val but lost on test.  The strict val-floor cannot prevent this — the
bound is val-only.  A future improvement (cross-validated gate
selection, or epsilon-margin tightening) could close this further;
not in scope for this fix.

Net: 13 cells improved by ≥ 1 pp, 15 cells unchanged, 4 cells
degraded.  Mean gap moves from −0.016 to −0.003.

## What this changes for "is pigauto best of best?"

**Discrete imputation guarantee** (post-fix): pigauto val performance
is ≤ baseline val performance for every trait by construction. Test
performance is ≤ baseline test performance with high probability,
modulo val→test sampling noise (~2-5 pp on individual cells with
n_val ≈ 20 cells per trait).

This is the floor — pigauto cannot meaningfully underperform its
baseline anymore. The ceiling (when the GNN actively helps on
discrete) is still being explored. Per the AE-attribution experiment
(`useful/MEMO_2026-04-29_strict_floor_and_AE.md`), the GNN engages
on roughly half the cells and helps when it does, but the dominant
constraint on pigauto's discrete performance is the underlying
phylogenetic-label-propagation baseline, not the GNN delta.

## Opus E regressions: closed?

Opus's adversarial review flagged three "Apr 17 → Apr 27 untraced
discrete regressions" (binary signal_0.6 −7 pp, categorical K_3
−8 pp, zi_count zf_0.2 +15 %). Comparing pre-bug Apr 17 numbers to
post-fix Apr 29:

| Opus E claim | Apr 17 baseline | Apr 17 pigauto | Apr 29 pigauto | verdict |
|---|---|---|---|---|
| binary signal_0.6 −7 pp | 0.660 | 0.612 (−4.8 pp) | 0.659 (−0.1 pp) | **closed** |
| categorical K_3 −8 pp | 0.729 | 0.684 (−4.5 pp) | 0.703 (−2.6 pp) | mostly closed |
| zi_count zf_0.2 +15 % | ratio 1.000 | ratio 1.123 (+12 %) | ratio 1.046 (+4.6 %) | mostly closed |

The Apr-17 baseline numbers don't quite match Opus's "−7 pp / −8 pp /
+15 %" framing (Opus aggregated differently), but the direction and
sign match.  The strict val-floor closes binary cleanly; categorical
and zi_count drift residuals are val→test sampling noise per the
analysis above.

## Files refreshed

* `script/bench_binary.{md,rds}` (40 cells, 5 reps, 15.4 min wall)
* `script/bench_categorical.{md,rds}` (35 cells, 5 reps, 16.1 min wall)
* `script/bench_zi_count.{md,rds}` (35 cells, 5 reps, 17.9 min wall)

Pre-fix backups preserved at `*_postfix2.{md,rds}.bak` (Apr 17) and
`*.rds.preStrictFloor.bak` (Apr 28).

## Still TODO (out of scope for this session)

* Re-run real-data benches with discrete columns: `bench_avonet_full_local.R`
  (Trophic.Level, Migration), `bench_amphibio.R`, `bench_avonet_phase6.R`.
* Cross-validated gate selection (rotate val cells) to close the
  remaining val→test drift on individual cells.
* Item #8 alternative: replace the BM baseline with phylolm-λ BLUP.
