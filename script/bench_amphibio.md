# AmphiBIO x pigauto

n = 5237 amphibian species x 2 traits (AmphiBIO + taxonomic tree).
Seed = 2026, miss_frac = 0.30, n_imputations = 1.

## Per-trait metrics

```
          method        trait    metric        value n_cells  wall_s
   mean_baseline Body_size_mm      rmse   89.0839479    1569   0.086
   mean_baseline Body_size_mm pearson_r           NA    1569   0.086
   mean_baseline  Body_mass_g      rmse 1947.7361511     178   0.086
   mean_baseline  Body_mass_g pearson_r           NA     178   0.086
 pigauto_default Body_size_mm      rmse   43.5734904    1569 237.081
 pigauto_default Body_size_mm pearson_r    0.8723709    1569 237.081
 pigauto_default  Body_mass_g      rmse 1422.3058875     178 237.081
 pigauto_default  Body_mass_g pearson_r    0.9779864     178 237.081
```
