# Phase C v2 cross-dataset bench: PanTHERIA + AmphiBIO (multi-seed)

Run: 2026-05-01 18:39:36
Seeds: 2026, 2027, 2028.  miss_frac = 0.30, N_IMP = 20.  AmphiBIO Diu / Noc dropped (presence-only encoding).

## Aggregated results (mean +/- SD across seeds)

```
   dataset               trait                type  config   metric
  amphibio         body_mass_g continuous_or_count    both     rmse
  amphibio         body_mass_g continuous_or_count   clamp     rmse
  amphibio         body_mass_g continuous_or_count default     rmse
  amphibio         body_mass_g continuous_or_count    mode     rmse
  amphibio        body_size_mm continuous_or_count    both     rmse
  amphibio        body_size_mm continuous_or_count   clamp     rmse
  amphibio        body_size_mm continuous_or_count default     rmse
  amphibio        body_size_mm continuous_or_count    mode     rmse
  amphibio        diet_breadth             ordinal    both accuracy
  amphibio        diet_breadth             ordinal   clamp accuracy
  amphibio        diet_breadth             ordinal default accuracy
  amphibio        diet_breadth             ordinal    mode accuracy
  amphibio             habitat         categorical    both accuracy
  amphibio             habitat         categorical   clamp accuracy
  amphibio             habitat         categorical default accuracy
  amphibio             habitat         categorical    mode accuracy
 pantheria         body_mass_g continuous_or_count    both     rmse
 pantheria         body_mass_g continuous_or_count   clamp     rmse
 pantheria         body_mass_g continuous_or_count default     rmse
 pantheria         body_mass_g continuous_or_count    mode     rmse
 pantheria        diet_breadth             ordinal    both accuracy
 pantheria        diet_breadth             ordinal   clamp accuracy
 pantheria        diet_breadth             ordinal default accuracy
 pantheria        diet_breadth             ordinal    mode accuracy
 pantheria         gestation_d continuous_or_count    both     rmse
 pantheria         gestation_d continuous_or_count   clamp     rmse
 pantheria         gestation_d continuous_or_count default     rmse
 pantheria         gestation_d continuous_or_count    mode     rmse
 pantheria     habitat_breadth             ordinal    both accuracy
 pantheria     habitat_breadth             ordinal   clamp accuracy
 pantheria     habitat_breadth             ordinal default accuracy
 pantheria     habitat_breadth             ordinal    mode accuracy
 pantheria head_body_length_mm continuous_or_count    both     rmse
 pantheria head_body_length_mm continuous_or_count   clamp     rmse
 pantheria head_body_length_mm continuous_or_count default     rmse
 pantheria head_body_length_mm continuous_or_count    mode     rmse
 pantheria         litter_size continuous_or_count    both     rmse
 pantheria         litter_size continuous_or_count   clamp     rmse
 pantheria         litter_size continuous_or_count default     rmse
 pantheria         litter_size continuous_or_count    mode     rmse
 pantheria     max_longevity_m continuous_or_count    both     rmse
 pantheria     max_longevity_m continuous_or_count   clamp     rmse
 pantheria     max_longevity_m continuous_or_count default     rmse
 pantheria     max_longevity_m continuous_or_count    mode     rmse
 pantheria      terrestriality              binary    both accuracy
 pantheria      terrestriality              binary   clamp accuracy
 pantheria      terrestriality              binary default accuracy
 pantheria      terrestriality              binary    mode accuracy
         mean           sd lift_pct_mean
 6.208817e+01 3.606728e+01    -8.5426052
 6.026134e+01 3.919227e+01    -5.3489239
 5.720166e+01 3.767265e+01     0.0000000
 6.437201e+01 3.574598e+01   -12.5352135
 5.894673e+01 1.614004e+01     1.0009623
 5.803515e+01 1.424760e+01     2.5319348
 5.954273e+01 1.890864e+01     0.0000000
 6.306131e+01 1.935083e+01    -5.9093388
 8.720930e-01 5.068487e-02     0.0000000
 8.720930e-01 5.068487e-02     0.0000000
 8.720930e-01 5.068487e-02     0.0000000
 8.720930e-01 5.068487e-02     0.0000000
 8.468254e-01 1.221810e-02     0.0000000
 8.468254e-01 1.221810e-02     0.0000000
 8.468254e-01 1.221810e-02     0.0000000
 8.468254e-01 1.221810e-02     0.0000000
 3.360240e+06 2.971154e+05     9.1568890
 3.430203e+06 3.807456e+05     7.2654728
 3.698949e+06 8.549679e+04     0.0000000
 3.365899e+06 3.372474e+05     9.0039100
 4.052632e-01 1.051315e-01    -0.1754386
 4.315789e-01 4.823764e-02     2.4561404
 4.070175e-01 3.083929e-02     0.0000000
 3.859649e-01 1.111926e-01    -2.1052632
 4.135812e+01 9.537802e+00     2.9379251
 4.017183e+01 1.088486e+01     5.7219874
 4.260997e+01 1.087209e+01     0.0000000
 3.685799e+01 7.173452e+00    13.4991369
 7.833333e-01 2.253640e-02     1.6666667
 7.807692e-01 3.003942e-02     1.4102564
 7.666667e-01 3.108809e-02     0.0000000
 7.820513e-01 1.554405e-02     1.5384615
 4.860769e+02 2.179838e+02     3.5946925
 4.872813e+02 1.774512e+02     3.3558135
 5.042013e+02 1.822635e+02     0.0000000
 4.606650e+02 1.737740e+02     8.6347048
 1.059290e+00 1.314449e-01     0.7825336
 1.095959e+00 7.955094e-02    -2.6519799
 1.067645e+00 7.342275e-02     0.0000000
 1.074286e+00 7.177261e-02    -0.6220172
 1.769672e+02 1.101648e+02   -28.8582317
 1.453222e+02 5.646448e+01    -5.8160113
 1.373348e+02 5.488878e+01     0.0000000
 1.308022e+02 4.029450e+01     4.7567573
 9.146667e-01 8.326664e-03    -0.2666667
 9.133333e-01 6.110101e-03    -0.4000000
 9.173333e-01 1.285820e-02     0.0000000
 9.173333e-01 1.285820e-02     0.0000000
```

## Raw per-seed results

```
   dataset seed               trait  config   metric        value     baseline
  amphibio 2026         body_mass_g    both     rmse 3.369999e+01 7.421740e+01
  amphibio 2027         body_mass_g    both     rmse 4.989227e+01 8.511810e+01
  amphibio 2028         body_mass_g    both     rmse 1.026723e+02 1.042163e+02
  amphibio 2026         body_mass_g   clamp     rmse 3.153014e+01 7.421740e+01
  amphibio 2027         body_mass_g   clamp     rmse 4.434620e+01 8.511810e+01
  amphibio 2028         body_mass_g   clamp     rmse 1.049077e+02 1.042163e+02
  amphibio 2026         body_mass_g default     rmse 3.531651e+01 7.421740e+01
  amphibio 2027         body_mass_g default     rmse 3.558646e+01 8.511810e+01
  amphibio 2028         body_mass_g default     rmse 1.007020e+02 1.042163e+02
  amphibio 2026         body_mass_g    mode     rmse 2.998521e+01 7.421740e+01
  amphibio 2027         body_mass_g    mode     rmse 6.179330e+01 8.511810e+01
  amphibio 2028         body_mass_g    mode     rmse 1.013375e+02 1.042163e+02
  amphibio 2026        body_size_mm    both     rmse 7.522400e+01 1.239157e+02
  amphibio 2027        body_size_mm    both     rmse 4.294751e+01 8.021672e+01
  amphibio 2028        body_size_mm    both     rmse 5.866868e+01 1.072698e+02
  amphibio 2026        body_size_mm   clamp     rmse 7.278771e+01 1.239157e+02
  amphibio 2027        body_size_mm   clamp     rmse 4.435288e+01 8.021672e+01
  amphibio 2028        body_size_mm   clamp     rmse 5.696485e+01 1.072698e+02
  amphibio 2026        body_size_mm default     rmse 7.195060e+01 1.239157e+02
  amphibio 2027        body_size_mm default     rmse 3.778019e+01 8.021672e+01
  amphibio 2028        body_size_mm default     rmse 6.889741e+01 1.072698e+02
  amphibio 2026        body_size_mm    mode     rmse 7.772506e+01 1.239157e+02
  amphibio 2027        body_size_mm    mode     rmse 4.112856e+01 8.021672e+01
  amphibio 2028        body_size_mm    mode     rmse 7.033032e+01 1.072698e+02
  amphibio 2026        diet_breadth    both accuracy 9.302326e-01 9.302326e-01
  amphibio 2027        diet_breadth    both accuracy 8.488372e-01 8.488372e-01
  amphibio 2028        diet_breadth    both accuracy 8.372093e-01 8.372093e-01
  amphibio 2026        diet_breadth   clamp accuracy 9.302326e-01 9.302326e-01
  amphibio 2027        diet_breadth   clamp accuracy 8.488372e-01 8.488372e-01
  amphibio 2028        diet_breadth   clamp accuracy 8.372093e-01 8.372093e-01
  amphibio 2026        diet_breadth default accuracy 9.302326e-01 9.302326e-01
  amphibio 2027        diet_breadth default accuracy 8.488372e-01 8.488372e-01
  amphibio 2028        diet_breadth default accuracy 8.372093e-01 8.372093e-01
  amphibio 2026        diet_breadth    mode accuracy 9.302326e-01 9.302326e-01
  amphibio 2027        diet_breadth    mode accuracy 8.488372e-01 8.488372e-01
  amphibio 2028        diet_breadth    mode accuracy 8.372093e-01 8.372093e-01
  amphibio 2026             habitat    both accuracy 8.500000e-01 8.071429e-01
  amphibio 2027             habitat    both accuracy 8.333333e-01 7.833333e-01
  amphibio 2028             habitat    both accuracy 8.571429e-01 7.785714e-01
  amphibio 2026             habitat   clamp accuracy 8.500000e-01 8.071429e-01
  amphibio 2027             habitat   clamp accuracy 8.333333e-01 7.833333e-01
  amphibio 2028             habitat   clamp accuracy 8.571429e-01 7.785714e-01
  amphibio 2026             habitat default accuracy 8.500000e-01 8.071429e-01
  amphibio 2027             habitat default accuracy 8.333333e-01 7.833333e-01
  amphibio 2028             habitat default accuracy 8.571429e-01 7.785714e-01
  amphibio 2026             habitat    mode accuracy 8.500000e-01 8.071429e-01
  amphibio 2027             habitat    mode accuracy 8.333333e-01 7.833333e-01
  amphibio 2028             habitat    mode accuracy 8.571429e-01 7.785714e-01
 pantheria 2026         body_mass_g    both     rmse 3.311722e+06 4.420580e+06
 pantheria 2027         body_mass_g    both     rmse 3.090370e+06 4.405038e+06
 pantheria 2028         body_mass_g    both     rmse 3.678628e+06 4.419351e+06
 pantheria 2026         body_mass_g   clamp     rmse 3.664786e+06 4.420580e+06
 pantheria 2027         body_mass_g   clamp     rmse 2.990893e+06 4.405038e+06
 pantheria 2028         body_mass_g   clamp     rmse 3.634928e+06 4.419351e+06
 pantheria 2026         body_mass_g default     rmse 3.791068e+06 4.420580e+06
 pantheria 2027         body_mass_g default     rmse 3.683634e+06 4.405038e+06
 pantheria 2028         body_mass_g default     rmse 3.622144e+06 4.419351e+06
 pantheria 2026         body_mass_g    mode     rmse 3.247876e+06 4.420580e+06
 pantheria 2027         body_mass_g    mode     rmse 3.103524e+06 4.405038e+06
 pantheria 2028         body_mass_g    mode     rmse 3.746296e+06 4.419351e+06
 pantheria 2026        diet_breadth    both accuracy 4.947368e-01 3.684211e-01
 pantheria 2027        diet_breadth    both accuracy 2.894737e-01 4.263158e-01
 pantheria 2028        diet_breadth    both accuracy 4.315789e-01 3.631579e-01
 pantheria 2026        diet_breadth   clamp accuracy 3.894737e-01 3.684211e-01
 pantheria 2027        diet_breadth   clamp accuracy 4.842105e-01 4.263158e-01
 pantheria 2028        diet_breadth   clamp accuracy 4.210526e-01 3.631579e-01
 pantheria 2026        diet_breadth default accuracy 3.842105e-01 3.684211e-01
 pantheria 2027        diet_breadth default accuracy 4.421053e-01 4.263158e-01
 pantheria 2028        diet_breadth default accuracy 3.947368e-01 3.631579e-01
 pantheria 2026        diet_breadth    mode accuracy 4.578947e-01 3.684211e-01
 pantheria 2027        diet_breadth    mode accuracy 2.578947e-01 4.263158e-01
 pantheria 2028        diet_breadth    mode accuracy 4.421053e-01 3.631579e-01
 pantheria 2026         gestation_d    both     rmse 3.059723e+01 9.283943e+01
 pantheria 2027         gestation_d    both     rmse 4.876886e+01 9.818021e+01
 pantheria 2028         gestation_d    both     rmse 4.470828e+01 9.816027e+01
 pantheria 2026         gestation_d   clamp     rmse 2.772429e+01 9.283943e+01
 pantheria 2027         gestation_d   clamp     rmse 4.488754e+01 9.818021e+01
 pantheria 2028         gestation_d   clamp     rmse 4.790366e+01 9.816027e+01
 pantheria 2026         gestation_d default     rmse 3.116544e+01 9.283943e+01
 pantheria 2027         gestation_d default     rmse 5.280093e+01 9.818021e+01
 pantheria 2028         gestation_d default     rmse 4.386354e+01 9.816027e+01
 pantheria 2026         gestation_d    mode     rmse 2.861497e+01 9.283943e+01
 pantheria 2027         gestation_d    mode     rmse 4.168512e+01 9.818021e+01
 pantheria 2028         gestation_d    mode     rmse 4.027388e+01 9.816027e+01
 pantheria 2026     habitat_breadth    both accuracy 7.576923e-01 7.192308e-01
 pantheria 2027     habitat_breadth    both accuracy 7.923077e-01 7.461538e-01
 pantheria 2028     habitat_breadth    both accuracy 8.000000e-01 7.269231e-01
 pantheria 2026     habitat_breadth   clamp accuracy 7.461538e-01 7.192308e-01
 pantheria 2027     habitat_breadth   clamp accuracy 7.961538e-01 7.461538e-01
 pantheria 2028     habitat_breadth   clamp accuracy 8.000000e-01 7.269231e-01
 pantheria 2026     habitat_breadth default accuracy 7.384615e-01 7.192308e-01
 pantheria 2027     habitat_breadth default accuracy 7.615385e-01 7.461538e-01
 pantheria 2028     habitat_breadth default accuracy 8.000000e-01 7.269231e-01
 pantheria 2026     habitat_breadth    mode accuracy 7.730769e-01 7.192308e-01
 pantheria 2027     habitat_breadth    mode accuracy 7.730769e-01 7.461538e-01
 pantheria 2028     habitat_breadth    mode accuracy 8.000000e-01 7.269231e-01
 pantheria 2026 head_body_length_mm    both     rmse 3.043600e+02 5.948335e+02
 pantheria 2027 head_body_length_mm    both     rmse 7.277702e+02 7.047560e+02
 pantheria 2028 head_body_length_mm    both     rmse 4.261003e+02 1.207376e+03
 pantheria 2026 head_body_length_mm   clamp     rmse 3.529642e+02 5.948335e+02
 pantheria 2027 head_body_length_mm   clamp     rmse 6.884477e+02 7.047560e+02
 pantheria 2028 head_body_length_mm   clamp     rmse 4.204319e+02 1.207376e+03
 pantheria 2026 head_body_length_mm default     rmse 3.758493e+02 5.948335e+02
 pantheria 2027 head_body_length_mm default     rmse 7.128221e+02 7.047560e+02
 pantheria 2028 head_body_length_mm default     rmse 4.239326e+02 1.207376e+03
 pantheria 2026 head_body_length_mm    mode     rmse 2.971440e+02 5.948335e+02
 pantheria 2027 head_body_length_mm    mode     rmse 6.431382e+02 7.047560e+02
 pantheria 2028 head_body_length_mm    mode     rmse 4.417129e+02 1.207376e+03
 pantheria 2026         litter_size    both     rmse 9.117147e-01 1.738224e+00
 pantheria 2027         litter_size    both     rmse 1.163800e+00 1.693335e+00
 pantheria 2028         litter_size    both     rmse 1.102357e+00 1.775852e+00
 pantheria 2026         litter_size   clamp     rmse 1.012579e+00 1.738224e+00
 pantheria 2027         litter_size   clamp     rmse 1.171029e+00 1.693335e+00
 pantheria 2028         litter_size   clamp     rmse 1.104269e+00 1.775852e+00
 pantheria 2026         litter_size default     rmse 9.829775e-01 1.738224e+00
 pantheria 2027         litter_size default     rmse 1.106178e+00 1.693335e+00
 pantheria 2028         litter_size default     rmse 1.113780e+00 1.775852e+00
 pantheria 2026         litter_size    mode     rmse 9.915253e-01 1.738224e+00
 pantheria 2027         litter_size    mode     rmse 1.119448e+00 1.693335e+00
 pantheria 2028         litter_size    mode     rmse 1.111884e+00 1.775852e+00
 pantheria 2026     max_longevity_m    both     rmse 2.987612e+02 1.581080e+02
 pantheria 2027     max_longevity_m    both     rmse 1.478656e+02 1.985683e+02
 pantheria 2028     max_longevity_m    both     rmse 8.427490e+01 1.364567e+02
 pantheria 2026     max_longevity_m   clamp     rmse 1.888145e+02 1.581080e+02
 pantheria 2027     max_longevity_m   clamp     rmse 1.656422e+02 1.985683e+02
 pantheria 2028     max_longevity_m   clamp     rmse 8.151003e+01 1.364567e+02
 pantheria 2026     max_longevity_m default     rmse 1.842362e+02 1.581080e+02
 pantheria 2027     max_longevity_m default     rmse 1.508026e+02 1.985683e+02
 pantheria 2028     max_longevity_m default     rmse 7.696569e+01 1.364567e+02
 pantheria 2026     max_longevity_m    mode     rmse 1.586905e+02 1.581080e+02
 pantheria 2027     max_longevity_m    mode     rmse 1.491121e+02 1.985683e+02
 pantheria 2028     max_longevity_m    mode     rmse 8.460388e+01 1.364567e+02
 pantheria 2026      terrestriality    both accuracy 9.080000e-01 5.960000e-01
 pantheria 2027      terrestriality    both accuracy 9.120000e-01 5.520000e-01
 pantheria 2028      terrestriality    both accuracy 9.240000e-01 5.960000e-01
 pantheria 2026      terrestriality   clamp accuracy 9.080000e-01 5.960000e-01
 pantheria 2027      terrestriality   clamp accuracy 9.120000e-01 5.520000e-01
 pantheria 2028      terrestriality   clamp accuracy 9.200000e-01 5.960000e-01
 pantheria 2026      terrestriality default accuracy 9.080000e-01 5.960000e-01
 pantheria 2027      terrestriality default accuracy 9.120000e-01 5.520000e-01
 pantheria 2028      terrestriality default accuracy 9.320000e-01 5.960000e-01
 pantheria 2026      terrestriality    mode accuracy 9.080000e-01 5.960000e-01
 pantheria 2027      terrestriality    mode accuracy 9.120000e-01 5.520000e-01
 pantheria 2028      terrestriality    mode accuracy 9.320000e-01 5.960000e-01
```

