# globtherm — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=20
N species used: 1969
Traits: Tmax, Tmin, lat_max, long_max, elevation_max
Log-transformed: (none)
Masked cells: 2074
Wall time: 840.2 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE + 95% conformal coverage / width)

|trait         |type       |     rmse| coverage_95| interval_width|
|:-------------|:----------|--------:|-----------:|--------------:|
|elevation_max |continuous | 625.7727|      0.9255|      2259.6573|
|lat_max       |continuous |  28.3093|      0.9439|       123.3226|
|long_max      |continuous |  73.8305|      0.9283|       341.4861|
|Tmax          |continuous |   4.6588|      0.9666|        23.6839|
|Tmin          |continuous |   6.3303|      0.9609|        27.3291|

