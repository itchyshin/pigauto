# Multi-obs + mixed-type benchmark (Phase 10 + B1 soft E-step)

Run on: 2026-04-18 06:12:47
Species per scenario: 150, reps: 2, epochs: 100
Total wall time: 11.2 min

Four methods, obs-level metrics on held-out cells:
- **species_mean**: observed-species mean (continuous) or mode (discrete); reference.
- **pigauto_hard**: `impute(..., multi_obs_aggregation = 'hard')` — Phase 10 default (threshold-at-0.5 / argmax).
- **pigauto_soft**: `impute(..., multi_obs_aggregation = 'soft')` — B1 Rao-Blackwell soft E-step.
- **pigauto_LP**: same call with `joint_mvn_available()` forced FALSE (legacy LP path).

hard_vs_LP: `pigauto_LP - pigauto_hard` for RMSE, `pigauto_hard - pigauto_LP` for accuracy (positive = hard better).
soft_vs_hard: `pigauto_hard - pigauto_soft` for RMSE, `pigauto_soft - pigauto_hard` for accuracy (positive = soft better).

```
              scenario trait   metric species_mean pigauto_hard pigauto_soft
 high_phylo_correlated     y accuracy       0.8838       0.7525       0.7688
 high_phylo_correlated     z accuracy       0.8556       0.7300       0.7328
 high_phylo_correlated    x1     RMSE       0.3638       0.9073       0.9073
 high_phylo_correlated    x2     RMSE       0.3476       0.8740       0.8739
  low_phylo_correlated     y accuracy       0.8907       0.6882       0.7021
  low_phylo_correlated     z accuracy       0.8770       0.5193       0.5163
  low_phylo_correlated    x1     RMSE       0.3486       0.7791       0.7755
  low_phylo_correlated    x2     RMSE       0.3166       0.6876       0.6910
  mod_phylo_correlated     y accuracy       0.8861       0.7147       0.7390
  mod_phylo_correlated     z accuracy       0.8770       0.6494       0.6331
  mod_phylo_correlated    x1     RMSE       0.3486       0.8508       0.8508
  mod_phylo_correlated    x2     RMSE       0.3166       0.7652       0.7643
 pigauto_LP hard_vs_LP soft_vs_hard
     0.6759   0.076527    1.635e-02
     0.5854   0.144637    2.790e-03
     0.9371   0.029732    6.746e-05
     0.9232   0.049284    9.229e-05
     0.6168   0.071361    1.392e-02
     0.3714   0.147983   -3.049e-03
     0.7838   0.004747    3.530e-03
     0.7220   0.034405   -3.427e-03
     0.6588   0.055873    2.427e-02
     0.3651   0.284242   -1.632e-02
     0.8386  -0.012263    7.954e-05
     0.7699   0.004687    8.540e-04
```

Phase 10 exit criterion (>=2pp lift on >=1 trait in >=2 scenarios): **MET**

B1 exit criterion (soft >= hard on >=1 trait in >=2 scenarios): **MET**
