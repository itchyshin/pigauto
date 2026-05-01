# CV-folds vs median_splits vs single_split bench

Run on: 2026-04-30 18:44:24
n=300 species, 5 reps, 3 gate methods, mask_frac=0.30

## Per-trait mean +/- SD across reps

```
   gate_method trait   metric      mean         sd
      cv_folds   bin accuracy 0.8222222 0.04303315
 median_splits   bin accuracy 0.8222222 0.04303315
  single_split   bin accuracy 0.8222222 0.04303315
      cv_folds  cat3 accuracy 0.7355556 0.07511308
 median_splits  cat3 accuracy 0.7355556 0.07511308
  single_split  cat3 accuracy 0.7355556 0.07511308
      cv_folds  cont     rmse 1.0610641 0.09601671
 median_splits  cont     rmse 1.0720463 0.10430196
  single_split  cont     rmse 1.0778249 0.11061847
```

## Comparison: cv_folds vs median_splits per trait

- **bin** (accuracy): single=0.8222, median=0.8222, cv=0.8222
- **cat3** (accuracy): single=0.7356, median=0.7356, cv=0.7356
- **cont** (rmse): single=1.078, median=1.072, cv=1.061
