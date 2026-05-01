# AVONET n=1500 multi-seed re-verification (2026-04-30)

Closes Opus adversarial review item E2 / #4: "Re-verify the AVONET
continuous evidence either at N_IMP=20 OR run 3-5 seeds at N_IMP=5
and report mean ± SD per trait."

## Method

`script/bench_avonet_full_local.R`, AVONET n=1500 random subset,
miss_frac = 0.30, N_IMP = 5, three seeds: 2026, 2027, 2028.  Same
binary (script unchanged except for `PIGAUTO_SEED` env-var support
in commit `9bfce30`).  Each run: ~8 min wall on Apple MPS.

The script was run sequentially (not parallel) to avoid GPU memory
contention.  Each seed produced one `script/bench_avonet_full_local_seed${seed}.rds`
+ `.md`.

## Results (pigauto vs mean baseline)

| trait | metric | mean_baseline | pigauto mean ± SD | lift % mean ± SD |
|---|---|---|---|---|
| Beak.Length_Culmen | RMSE | 20.4 | **8.89 ± 1.95** | 55.2 ± 14.6 |
| Mass               | RMSE | 1041 | **702 ± 403** | 35.2 ± 23.3 |
| Tarsus.Length      | RMSE | 28.5 | **16.0 ± 5.2** | 44.2 ± 17.4 |
| Wing.Length        | RMSE | 100  | **34.5 ± 8.4** | 65.8 ± 5.5  |
| Migration          | acc  | 0.787 | **0.776 ± 0.036** | -1.4 ± 1.1 |
| Primary.Lifestyle  | acc  | 0.573 | **0.809 ± 0.029** | 41.4 ± 7.1  |
| Trophic.Level      | acc  | 0.557 | **0.826 ± 0.007** | 48.5 ± 5.9  |

(Lift = 100 × (baseline − pigauto)/baseline for RMSE; 100 × (pigauto −
baseline)/baseline for accuracy.  Positive = pigauto better.)

## Per-seed pigauto values (continuous, RMSE)

| trait | seed 2026 | seed 2027 | seed 2028 |
|---|---|---|---|
| Beak.Length_Culmen | 6.97 | 10.86 | 8.84 |
| Mass               | 642.3 | 332.0 | **1131.0** |
| Tarsus.Length      | 11.01 | 21.46 | 15.42 |
| Wing.Length        | 35.09 | 25.92 | 42.60 |

## Per-seed pigauto values (discrete, accuracy)

| trait | seed 2026 | seed 2027 | seed 2028 |
|---|---|---|---|
| Migration          | 0.815 | 0.744 | 0.769 |
| Primary.Lifestyle  | 0.838 | 0.809 | 0.780 |
| Trophic.Level      | 0.829 | 0.831 | 0.818 |

## Conformal coverage95 (target = 0.95) and MC-dropout coverage95

| trait | conformal mean ± SD | MC-dropout mean ± SD |
|---|---|---|
| Beak.Length_Culmen | 0.955 ± 0.007 | 0.679 ± 0.038 |
| Mass               | 0.970 ± 0.024 | 0.702 ± 0.012 |
| Tarsus.Length      | 0.931 ± 0.018 | 0.684 ± 0.018 |
| Wing.Length        | 0.919 ± 0.031 | 0.686 ± 0.021 |

Conformal is on target; MC-dropout under-covers at ~0.69, consistent
with the known compactness of MC-dropout intervals at this
N_IMP (5) and pooling structure.

## Interpretation

1. **Continuous-trait pigauto values are noisy across seeds at
   N_IMP=5.**  Mass shows the strongest swing (332→642→1131 over the
   three seeds, CV ≈ 57 %).  This is consistent with the Opus E2
   prediction that the v3 memo's single-seed evidence (Mass=377)
   was within the noise envelope.

2. **The strict val-floor v3 split (continuous mean-only check;
   discrete strict both-corners check) does not appear to hurt
   continuous performance**: pigauto beats baseline on every
   continuous trait at every seed (lift bands all positive, never
   negative).  But the per-seed RMSE numbers themselves are unstable
   enough that we cannot use single-seed continuous numbers as
   evidence to distinguish v1 vs v3 vs pre-fix in either direction.

3. **Discrete traits are essentially deterministic across seeds.**
   Trophic.Level CV = 0.86 %, Primary.Lifestyle CV = 3.6 %, Migration
   CV = 4.7 % — orders of magnitude smaller than continuous.  This
   matches the trait-type bisection: model state and gate values are
   deterministic; what varies is the MC-dropout pooling order on
   continuous predict-time refines.  Discrete predictions go
   through `argmax` and `pmax(p, 0.01)` clamping that absorbs
   sub-bit-level noise.

4. **Migration regresses by 1.4 pp** (vs baseline 0.787 → pigauto
   0.776).  This is the known B3 (commit a541dbd) threshold-joint
   ordinal regression (`useful/MEMO_2026-04-29_phase6_migration_bisect.md`).
   K=3 ordinal threshold-joint is under-determined; the calibrated
   gate cannot route around it because LP is not currently a corner
   in the safety-floor simplex grid.  Phase F (LP corner for ordinal
   traits) addresses this.

5. **Lift is robustly positive on the bird traits where the GNN
   has signal**: Wing 65.8 ± 5.5 % (n=3 seeds bracket).  Even Mass,
   the noisiest, has lower bound 35.2 − 23.3 = +11.9 % over baseline.

## What this changes vs the v3 single-seed memo

The numerical headline in `MEMO_2026-04-29_strict_floor_v3.md` of
"Mass −2.75 RMSE mean change" is now confirmed to be within
pooling noise.  v3 should be defended on its discrete-bench
closures and architectural reasoning, not on the AVONET continuous
single-seed numbers.  The multi-seed evidence here neither rejects
nor confirms v3 vs pre on continuous; it confirms the noise envelope
is too wide to distinguish them at N_IMP=5.

## Files

* `script/bench_avonet_full_local_seed{2026,2027,2028}.{rds,md}`
* `/tmp/bench_avonet_seed{2026,2027,2028}.log` — per-seed training logs

## Open

* True N_IMP=20 multi-seed re-verify deferred (each run is ~32 min;
  3 seeds × 32 min = 1.6 h; not in this session's budget).
* Phase F (LP corner for ordinal) would close the Migration
  regression and is queued next.
