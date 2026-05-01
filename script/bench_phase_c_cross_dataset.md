# Phase C cross-dataset bench: PanTHERIA + AmphiBIO

Run: 2026-05-01 16:23:38
Seed = 2026, miss_frac = 0.30, N_IMP = 20

## Per-cell results

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
  amphibio                 diu              binary    both accuracy
  amphibio                 diu              binary   clamp accuracy
  amphibio                 diu              binary default accuracy
  amphibio                 diu              binary    mode accuracy
  amphibio             habitat         categorical    both accuracy
  amphibio             habitat         categorical   clamp accuracy
  amphibio             habitat         categorical default accuracy
  amphibio             habitat         categorical    mode accuracy
  amphibio                 noc              binary    both accuracy
  amphibio                 noc              binary   clamp accuracy
  amphibio                 noc              binary default accuracy
  amphibio                 noc              binary    mode accuracy
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
        value     baseline    lift_pct
 5.289177e+01 7.421740e+01   9.5541094
 5.049154e+01 7.421740e+01  13.6585375
 5.847891e+01 7.421740e+01   0.0000000
 5.467942e+01 7.421740e+01   6.4971839
 7.014187e+01 1.239157e+02   4.9029195
 7.583089e+01 1.239157e+02  -2.8101488
 7.375818e+01 1.239157e+02   0.0000000
 6.830356e+01 1.239157e+02   7.3952735
 9.069767e-01 9.069767e-01   0.0000000
 9.069767e-01 9.069767e-01   0.0000000
 9.069767e-01 9.069767e-01   0.0000000
 9.069767e-01 9.069767e-01   0.0000000
 8.488889e-01 8.488889e-01   0.0000000
 8.444444e-01 8.488889e-01  -0.4444444
 8.488889e-01 8.488889e-01   0.0000000
 4.511111e-01 8.488889e-01 -39.7777778
 8.666667e-01 8.142857e-01   0.0000000
 8.619048e-01 8.142857e-01  -0.4761905
 8.666667e-01 8.142857e-01   0.0000000
 8.547619e-01 8.142857e-01  -1.1904762
 7.666667e-01 7.666667e-01   0.0000000
 7.666667e-01 7.666667e-01   0.0000000
 7.666667e-01 7.666667e-01   0.0000000
 7.688889e-01 7.666667e-01   0.2222222
 3.221207e+06 4.420580e+06   0.2720725
 3.093925e+06 4.420580e+06   4.2127066
 3.229995e+06 4.420580e+06   0.0000000
 3.881528e+06 4.420580e+06 -20.1713155
 4.000000e-01 3.684211e-01   0.5263158
 3.789474e-01 3.684211e-01  -1.5789474
 3.947368e-01 3.684211e-01   0.0000000
 4.368421e-01 3.684211e-01   4.2105263
 2.879498e+01 9.283943e+01  18.4154576
 2.880225e+01 9.283943e+01  18.3948765
 3.529465e+01 9.283943e+01   0.0000000
 3.192939e+01 9.283943e+01   9.5347597
 7.730769e-01 7.192308e-01   5.7692308
 7.000000e-01 7.192308e-01  -1.5384615
 7.153846e-01 7.192308e-01   0.0000000
 7.884615e-01 7.192308e-01   7.3076923
 3.501442e+02 5.948335e+02 -14.5277994
 3.361448e+02 5.948335e+02  -9.9487657
 3.057285e+02 5.948335e+02   0.0000000
 3.582857e+02 5.948335e+02 -17.1907949
 9.700146e-01 1.738224e+00  -4.5532044
 1.002107e+00 1.738224e+00  -8.0123450
 9.277713e-01 1.738224e+00   0.0000000
 9.936507e-01 1.738224e+00  -7.1008321
 1.670913e+02 1.581080e+02  25.6899399
 1.697955e+02 1.581080e+02  24.4872988
 2.248568e+02 1.581080e+02   0.0000000
 2.495511e+02 1.581080e+02 -10.9822293
 9.080000e-01 5.960000e-01   0.0000000
 9.080000e-01 5.960000e-01   0.0000000
 9.080000e-01 5.960000e-01   0.0000000
 9.080000e-01 5.960000e-01   0.0000000
```

