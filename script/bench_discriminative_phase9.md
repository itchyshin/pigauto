# Phase 9 A/B: Level-C baseline + transformer GNN vs legacy GNN

Run on: 2026-04-16 07:16:18
Species per scenario: 200, reps: 2, epochs: 100, total wall: 826 s

Three methods:
- **lp**: LP baseline + legacy GNN (use_transformer_blocks = FALSE)
- **levelC_legacy**: Phase 6 joint-threshold baseline + legacy GNN
- **levelC_transformer**: Phase 6 joint-threshold baseline + Phase 9 transformer GNN

Metrics: rmse = continuous RMSE (lower=better), accb = binary accuracy,
accc = categorical accuracy (both higher=better), wall = training seconds.

Wide summary (mean over reps):

```
               scenario lp_rmse lp_accb lp_accc lp_wall levelC_legacy_rmse
 high_signal_correlated  0.5581  0.6253  0.5631   34.30             0.4415
      high_signal_indep  0.5579  0.7757  0.5143   32.05             0.5682
  low_signal_correlated  0.9137  0.6278  0.4869   29.20             0.8539
  mod_signal_correlated  0.7842  0.6516  0.4417   30.50             0.7509
 levelC_legacy_accb levelC_legacy_accc levelC_legacy_wall
             0.8759             0.7048              30.15
             0.9236             0.5619              31.85
             0.7544             0.3643              29.25
             0.7807             0.4857              29.45
 levelC_transformer_rmse levelC_transformer_accb levelC_transformer_accc
                  0.4415                  0.8759                  0.7048
                  0.5682                  0.9236                  0.5619
                  0.8349                  0.7544                  0.3643
                  0.7174                  0.7807                  0.4857
 levelC_transformer_wall
                   38.50
                   41.05
                   36.95
                   37.55
```

Raw per-rep results:

```
               scenario rep             method rmse_cont   acc_bin   acc_cat
 high_signal_correlated   1                 lp 0.4694347 0.6190476 0.6500000
 high_signal_correlated   1      levelC_legacy 0.3960822 0.8571429 0.6000000
 high_signal_correlated   1 levelC_transformer 0.3960822 0.8571429 0.6000000
 high_signal_correlated   2                 lp 0.6466765 0.6315789 0.4761905
 high_signal_correlated   2      levelC_legacy 0.4869542 0.8947368 0.8095238
 high_signal_correlated   2 levelC_transformer 0.4869542 0.8947368 0.8095238
  mod_signal_correlated   1                 lp 0.7430108 0.6190476 0.5500000
  mod_signal_correlated   1      levelC_legacy 0.6222303 0.6666667 0.4000000
  mod_signal_correlated   1 levelC_transformer 0.5553552 0.6666667 0.4000000
  mod_signal_correlated   2                 lp 0.8254424 0.6842105 0.3333333
  mod_signal_correlated   2      levelC_legacy 0.8794802 0.8947368 0.5714286
  mod_signal_correlated   2 levelC_transformer 0.8794802 0.8947368 0.5714286
  low_signal_correlated   1                 lp 0.8959158 0.5714286 0.4500000
  low_signal_correlated   1      levelC_legacy 0.7498539 0.6666667 0.3000000
  low_signal_correlated   1 levelC_transformer 0.7118622 0.6666667 0.3000000
  low_signal_correlated   2                 lp 0.9314040 0.6842105 0.5238095
  low_signal_correlated   2      levelC_legacy 0.9579329 0.8421053 0.4285714
  low_signal_correlated   2 levelC_transformer 0.9579329 0.8421053 0.4285714
      high_signal_indep   1                 lp 0.5966872 0.7619048 0.6000000
      high_signal_indep   1      levelC_legacy 0.5831859 0.9523810 0.6000000
      high_signal_indep   1 levelC_transformer 0.5831859 0.9523810 0.6000000
      high_signal_indep   2                 lp 0.5191563 0.7894737 0.4285714
      high_signal_indep   2      levelC_legacy 0.5531421 0.8947368 0.5238095
      high_signal_indep   2 levelC_transformer 0.5531421 0.8947368 0.5238095
 wall
 36.3
 29.6
 38.4
 32.3
 30.7
 38.6
 30.2
 29.8
 37.1
 30.8
 29.1
 38.0
 30.0
 29.4
 36.8
 28.4
 29.1
 37.1
 31.3
 31.0
 37.3
 32.8
 32.7
 44.8
```
