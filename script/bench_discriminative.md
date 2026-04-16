# Phase 8 discriminative benchmark

Run on: 2026-04-16 05:38:04
Species per scenario: 200, reps per scenario: 3

Level-C baseline (Phase 6 active) vs LP (force_lp).

Continuous RMSE (lower better; lift = LP - LevelC):

```
                   scenario   rmse_lc   rmse_lp     rmse_lift
     high_signal_correlated 0.4158045 0.5189941  0.1031895859
    high_signal_independent 0.4441700 0.4435881 -0.0005818633
      low_signal_correlated 0.7598528 0.8401790  0.0803261731
 moderate_signal_correlated 0.5241917 0.6212241  0.0970324258
```

Binary accuracy (higher better; lift = LevelC - LP):

```
                   scenario  acc_y_lc  acc_y_lp acc_y_lift
     high_signal_correlated 0.8388619 0.6717775 0.16708438
    high_signal_independent 0.7811686 0.7494226 0.03174603
      low_signal_correlated 0.7382181 0.6734483 0.06476977
 moderate_signal_correlated 0.7582191 0.7477517 0.01046734
```

Categorical (K=4) accuracy (higher better; lift = LevelC - LP):

```
                   scenario  acc_z_lc  acc_z_lp  acc_z_lift
     high_signal_correlated 0.7079365 0.6134921  0.09444444
    high_signal_independent 0.6285714 0.5968254  0.03174603
      low_signal_correlated 0.4095238 0.4674603 -0.05793651
 moderate_signal_correlated 0.7261905 0.5833333  0.14285714
```

