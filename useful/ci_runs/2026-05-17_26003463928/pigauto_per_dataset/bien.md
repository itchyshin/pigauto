# bien — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=20
N species used: 2000
Traits: height_m, leaf_area, sla, seed_mass, wood_density
Log-transformed: height_m, leaf_area, sla, seed_mass, wood_density
Masked cells: 913
Wall time: 998.7 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE + 95% conformal coverage / width)

|trait        |type       |   rmse| coverage_95| interval_width|
|:------------|:----------|------:|-----------:|--------------:|
|height_m     |continuous | 1.6371|      0.8276|         4.3609|
|leaf_area    |continuous | 1.7384|      0.9667|         8.0431|
|seed_mass    |continuous | 2.1045|      0.8871|         6.0451|
|sla          |continuous | 2.5571|      1.0000|        26.9881|
|wood_density |continuous | 0.2453|      0.9958|         1.5045|

