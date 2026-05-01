# Phase F smoke bench: AVONET Migration ordinal regression

Run: 2026-05-01 07:38:56
n=1500 species random subset, miss_frac=0.30, N_IMP=1 (single pass)
Pre-registered target: Migration accuracy >= 0.78 on >= 2 of 3 seeds

## Results

```
 seed pigauto_acc baseline_acc   lift_pp path_chosen n_test
 2030   0.7555556    0.7888889 -3.333333      bm_mvn    450
 2031   0.7755556    0.7977778 -2.222222      bm_mvn    450
 2032   0.7711111    0.8133333 -4.222222      bm_mvn    450
```

## Verdict

**PHASE F INCONCLUSIVE**: only 0 / 3 seeds reach Migration acc >= 0.78

- Mean pigauto Migration acc: 0.767 (SD 0.011)
- Mean baseline Migration acc: 0.800 (SD 0.012)
- Mean lift over baseline: -3.3 pp (SD 1.0)

## Pre-Phase-F baseline (from `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md`)

Migration acc with pigauto N_IMP=20: 0.713 +/- 0.032
Migration acc with mean baseline:    0.800 +/- 0.013
Pre-Phase-F regression: -10.9 +/- 3.1 pp

## Path selection

Phase F adds 'lp' as a third option in the per-trait ordinal path
selection (alongside 'threshold_joint' and 'bm_mvn'). The path_chosen
column above shows which option was selected on each seed.
