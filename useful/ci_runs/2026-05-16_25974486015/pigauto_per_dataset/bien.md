# bien — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: height_m, leaf_area, sla, seed_mass, wood_density
Log-transformed: height_m, leaf_area, sla, seed_mass, wood_density
Masked cells: 913
Wall time: 1065.1 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE, fit scale; log for log_traits)

|trait        |type       |   rmse|
|:------------|:----------|------:|
|height_m     |continuous | 1.9284|
|leaf_area    |continuous | 3.1353|
|seed_mass    |continuous | 3.5475|
|sla          |continuous | 7.4561|
|wood_density |continuous | 0.4559|

