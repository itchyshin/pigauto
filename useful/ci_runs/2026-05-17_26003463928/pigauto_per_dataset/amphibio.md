# amphibio — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=20
N species used: 2000
Traits: body_size_mm, body_mass_g
Log-transformed: body_size_mm, body_mass_g
Masked cells: 518
Wall time: 784.3 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE + 95% conformal coverage / width)

|trait        |type       |   rmse| coverage_95| interval_width|
|:------------|:----------|------:|-----------:|--------------:|
|body_mass_g  |continuous | 1.3541|      1.0000|          6.216|
|body_size_mm |continuous | 0.4097|      0.9461|          1.660|

