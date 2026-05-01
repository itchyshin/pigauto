# Phase H+ smoke bench: cross-K ordinal pooling (median vs mode)

Run: 2026-05-01 12:07:47
n_species=200, miss_frac=0.30, N_IMP=20, K=3,5,7, lambda=0.6, reps=3

## Per-cell results

```
 K rep pool_method pigauto_acc baseline_acc   lift_pp
 3   1      median   0.6666667   0.28333333 38.333333
 3   1        mode   0.7166667   0.28333333 43.333333
 3   2      median   0.3666667   0.30000000  6.666667
 3   2        mode   0.3666667   0.30000000  6.666667
 3   3      median   0.5166667   0.30000000 21.666667
 3   3        mode   0.5166667   0.30000000 21.666667
 5   1      median   0.4000000   0.13333333 26.666667
 5   1        mode   0.4333333   0.13333333 30.000000
 5   2      median   0.5500000   0.13333333 41.666667
 5   2        mode   0.3166667   0.13333333 18.333333
 5   3      median   0.4666667   0.13333333 33.333333
 5   3        mode   0.4666667   0.13333333 33.333333
 7   1      median   0.3500000   0.08333333 26.666667
 7   1        mode   0.3500000   0.08333333 26.666667
 7   2      median   0.2500000   0.10000000 15.000000
 7   2        mode   0.2333333   0.10000000 13.333333
 7   3      median   0.1833333   0.11666667  6.666667
 7   3        mode   0.3000000   0.11666667 18.333333
```

## Aggregated by K and pool_method

```
 K pool_method pigauto_acc.mean pigauto_acc.sd lift_pp.mean lift_pp.sd
 3      median       0.51666667     0.15000000    22.222222  15.840642
 5      median       0.47222222     0.07515416    33.888889   7.515416
 7      median       0.26111111     0.08388705    16.111111  10.046190
 3        mode       0.53333333     0.17559423    23.888889  18.434067
 5        mode       0.40555556     0.07876359    27.222222   7.876359
 7        mode       0.29444444     0.05853141    19.444444   6.735753
```

Mode wins (>= median acc) on 7 / 9 cells (78%).
Worst per-cell mode-vs-median delta: -23.3 pp.

## Verdict

**PHASE H+ INCONCLUSIVE**: mode wins on only 7 / 9 cells (78%, target >= 80%).  Default should remain 'median'.

