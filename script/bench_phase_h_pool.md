# Phase H smoke bench: ordinal pooling (mean+round vs mode)

Run: 2026-05-01 09:21:09
AVONET n=1500 random subset, miss_frac=0.30, seeds=2030,2031,2032, N_IMP=1,20, pools=mode

## Per-seed Migration accuracy

```
 pool_method seed n_imp pigauto_acc baseline_acc    lift_pp
        mode 2030     1   0.7555556    0.7888889 -3.3333333
        mode 2030    20   0.7666667    0.7888889 -2.2222222
        mode 2031     1   0.7755556    0.7977778 -2.2222222
        mode 2031    20   0.7888889    0.7977778 -0.8888889
        mode 2032     1   0.7711111    0.8133333 -4.2222222
        mode 2032    20   0.7822222    0.8133333 -3.1111111
```

## Aggregated (mean +/- SD over 3 seeds)

```
 pool_method n_imp pigauto_acc.mean pigauto_acc.sd lift_pp.mean lift_pp.sd
        mode     1       0.76740741     0.01050181    -3.259259   1.002056
        mode    20       0.77925926     0.01140356    -2.074074   1.118494
```

## Verdict

PHASE H PASSES: with pool_method = 'mode', N_IMP=20 acc >= N_IMP=1 acc on 3 / 3 seeds

