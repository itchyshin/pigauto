# CV-folds vs single_split: default-flip evidence

Run: 2026-05-01 06:03:52
n=300 species, 10 reps x 2 gate methods x 3 lambdas = 60 cells

## Per-trait per-lambda mean +/- SD across reps

```
  gate_method trait   metric lambda      mean         sd  n
     cv_folds   bin accuracy    0.4 0.7911111 0.07253966 10
 single_split   bin accuracy    0.4 0.7933333 0.07429596 10
     cv_folds   bin accuracy    0.7 0.7944444 0.06589054 10
 single_split   bin accuracy    0.7 0.7911111 0.06453378 10
     cv_folds   bin accuracy    1.0 0.8300000 0.06749486 10
 single_split   bin accuracy    1.0 0.8233333 0.07036062 10
     cv_folds  cat3 accuracy    0.4 0.7144444 0.07389794 10
 single_split  cat3 accuracy    0.4 0.7022222 0.09569614 10
     cv_folds  cat3 accuracy    0.7 0.7855556 0.06286485 10
 single_split  cat3 accuracy    0.7 0.7833333 0.06225088 10
     cv_folds  cat3 accuracy    1.0 0.6955556 0.05930554 10
 single_split  cat3 accuracy    1.0 0.6933333 0.05766371 10
     cv_folds  cont     rmse    0.4 0.4352492 0.05170721 10
 single_split  cont     rmse    0.4 0.4357323 0.05193312 10
     cv_folds  cont     rmse    0.7 0.7206674 0.04974633 10
 single_split  cont     rmse    0.7 0.7344590 0.06534026 10
     cv_folds  cont     rmse    1.0 1.0858668 0.10126730 10
 single_split  cont     rmse    1.0 1.0893233 0.09782274 10
```

## Paired delta (cv_folds - single_split) per cell, summarised

Lower delta is better for RMSE; higher is better for accuracy.
n_pos / n_neg = number of reps where cv_folds was bigger / smaller
than single_split.

```
 trait lambda    delta_mean    delta_sd n_pos n_neg
   bin    0.4 -0.0022222222 0.007027284     0     1
   bin    0.7  0.0033333333 0.010540926     1     0
   bin    1.0  0.0066666667 0.021081851     1     0
  cat3    0.4  0.0122222222 0.038650060     1     0
  cat3    0.7  0.0022222222 0.007027284     1     0
  cat3    1.0  0.0022222222 0.007027284     1     0
  cont    0.4 -0.0004831863 0.003636568     1     3
  cont    0.7 -0.0137915442 0.027280146     0     4
  cont    1.0 -0.0034564557 0.020187381     1     2
```

## Decision

- Mean continuous RMSE lift (cv - single, negated): +0.00591
- Mean discrete accuracy gap (cv - single): +0.004074
- Continuous lift >= 1% of baseline RMSE? NO
- Discrete gap >= -0.5 pp? YES

**RECOMMENDATION: keep single_split as default; cv_folds remains opt-in.**
