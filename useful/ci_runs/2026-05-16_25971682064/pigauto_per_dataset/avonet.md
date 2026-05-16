# avonet — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: mass_g, wing_length_mm, beak_length_culmen_mm, tarsus_length_mm, trophic_level, primary_lifestyle, migration
Log-transformed: mass_g, wing_length_mm, beak_length_culmen_mm, tarsus_length_mm
Masked cells: 4197
Wall time: 1344.5 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE, fit scale; log for log_traits)

|trait                 |type       |   rmse|
|:---------------------|:----------|------:|
|beak_length_culmen_mm |continuous | 0.4208|
|mass_g                |continuous | 0.7768|
|tarsus_length_mm      |continuous | 0.4152|
|wing_length_mm        |continuous | 0.5547|

## Discrete-family (accuracy)

|trait             |type        | accuracy|
|:-----------------|:-----------|--------:|
|primary_lifestyle |categorical |   0.3500|
|trophic_level     |categorical |   0.3573|
|migration         |ordinal     |   0.6446|

