# Multi-obs + mixed-type benchmark (Phase 10 + B1 soft E-step)

Run on: 2026-04-16 17:42:10
Species per scenario: 150, reps: 2, epochs: 100
Total wall time: 12.0 min

Four methods, obs-level metrics on held-out cells:
- **species_mean**: observed-species mean (continuous) or mode (discrete); reference.
- **pigauto_hard**: `impute(..., multi_obs_aggregation = 'hard')` — Phase 10 default (threshold-at-0.5 / argmax).
- **pigauto_soft**: `impute(..., multi_obs_aggregation = 'soft')` — B1 Rao-Blackwell soft E-step.
- **pigauto_LP**: same call with `joint_mvn_available()` forced FALSE (legacy LP path).

hard_vs_LP: `pigauto_LP - pigauto_hard` for RMSE, `pigauto_hard - pigauto_LP` for accuracy (positive = hard better).
soft_vs_hard: `pigauto_hard - pigauto_soft` for RMSE, `pigauto_soft - pigauto_hard` for accuracy (positive = soft better).

```
              scenario trait   metric species_mean pigauto_hard pigauto_soft
 high_phylo_correlated     y accuracy       0.8838       0.7643       0.7823
 high_phylo_correlated     z accuracy       0.8651       0.7418       0.7614
 high_phylo_correlated    x1     RMSE       0.3638       0.9425       0.9461
 high_phylo_correlated    x2     RMSE       0.3476       0.9219       0.9252
  low_phylo_correlated     y accuracy       0.8907       0.6882       0.7021
  low_phylo_correlated     z accuracy       0.8770       0.5193       0.5163
  low_phylo_correlated    x1     RMSE       0.3486       0.7791       0.7744
  low_phylo_correlated    x2     RMSE       0.3166       0.6910       0.6876
  mod_phylo_correlated     y accuracy       0.8861       0.7147       0.7390
  mod_phylo_correlated     z accuracy       0.8770       0.6494       0.6331
  mod_phylo_correlated    x1     RMSE       0.3486       0.8508       0.8508
  mod_phylo_correlated    x2     RMSE       0.3166       0.7652       0.7643
 pigauto_LP hard_vs_LP soft_vs_hard
     0.7214   0.042857    1.804e-02
     0.6392   0.102549    1.962e-02
     0.9670   0.024508   -3.664e-03
     0.9409   0.018971   -3.340e-03
     0.6199   0.068313    1.392e-02
     0.3995   0.119794   -3.049e-03
     0.7838   0.004747    4.672e-03
     0.7220   0.030931    3.430e-03
     0.6588   0.055873    2.427e-02
     0.3651   0.284242   -1.632e-02
     0.8386  -0.012263    7.954e-05
     0.7699   0.004687    8.540e-04
```

Phase 10 exit criterion (>=2pp lift on >=1 trait in >=2 scenarios): **MET**

B1 exit criterion (soft >= hard on >=1 trait in >=2 scenarios): **MET**
