# avonet — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=20
N species used: 2000
Traits: mass_g, wing_length_mm, beak_length_culmen_mm, tarsus_length_mm, trophic_level, primary_lifestyle, migration
Log-transformed: mass_g, wing_length_mm, beak_length_culmen_mm, tarsus_length_mm
Masked cells: 4197
Wall time: 1391.3 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE + 95% conformal coverage / width)

|trait                 |type       |   rmse| coverage_95| interval_width|
|:---------------------|:----------|------:|-----------:|--------------:|
|beak_length_culmen_mm |continuous | 0.2848|      0.9667|         1.4121|
|mass_g                |continuous | 0.5963|      0.9517|         2.3591|
|tarsus_length_mm      |continuous | 0.2757|      0.9567|         1.1542|
|wing_length_mm        |continuous | 0.4676|      0.9600|         0.9918|

## Discrete-family (accuracy + Brier)

|trait             |type        | accuracy|  brier|
|:-----------------|:-----------|--------:|------:|
|primary_lifestyle |categorical |   0.8283| 0.5135|
|trophic_level     |categorical |   0.8247| 0.4608|
|migration         |ordinal     |   0.7475|     NA|

