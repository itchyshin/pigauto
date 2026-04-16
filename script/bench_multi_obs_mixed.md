# Multi-obs + mixed-type benchmark (Phase 10 validation)

Run on: 2026-04-16 11:56:08
Species per scenario: 150, reps: 2, epochs: 100
Total wall time: 7.7 min

Three methods, obs-level metrics on held-out cells:
- **species_mean**: observed-species mean (continuous) or mode (discrete); reference.
- **pigauto_levelC**: `impute()` with Phase 10 Level-C default (joint MVN / threshold-joint / OVR via species-level aggregation).
- **pigauto_LP**: same call with `joint_mvn_available()` forced FALSE (legacy LP path).

Lift column: `pigauto_LP - pigauto_levelC` for RMSE (positive = Level-C better); `pigauto_levelC - pigauto_LP` for accuracy (positive = Level-C better).

```
              scenario trait   metric species_mean pigauto_levelC pigauto_LP
 high_phylo_correlated     y accuracy       0.8838         0.7794     0.7046
 high_phylo_correlated     z accuracy       0.8556         0.6745     0.4861
 high_phylo_correlated    x1     RMSE       0.3638         0.9234     0.9649
 high_phylo_correlated    x2     RMSE       0.3476         0.8840     0.8987
  low_phylo_correlated     y accuracy       0.8907         0.6882     0.6351
  low_phylo_correlated     z accuracy       0.8770         0.5193     0.3921
  low_phylo_correlated    x1     RMSE       0.3486         0.7791     0.7838
  low_phylo_correlated    x2     RMSE       0.3166         0.6941     0.7220
  mod_phylo_correlated     y accuracy       0.8861         0.7147     0.6588
  mod_phylo_correlated     z accuracy       0.8770         0.6494     0.3651
  mod_phylo_correlated    x1     RMSE       0.3486         0.8508     0.8386
  mod_phylo_correlated    x2     RMSE       0.3166         0.7652     0.7699
 lift_vs_LP
   0.074844
   0.188408
   0.041448
   0.014622
   0.053069
   0.127212
   0.004747
   0.027845
   0.055873
   0.284242
  -0.012263
   0.004687
```

Exit criterion (>=2pp lift on >=1 trait in >=2 scenarios): **MET**
