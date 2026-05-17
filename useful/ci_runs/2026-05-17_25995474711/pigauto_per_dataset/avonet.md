# avonet — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: mass_g, wing_length_mm, beak_length_culmen_mm, tarsus_length_mm, trophic_level, primary_lifestyle, migration
Log-transformed: mass_g, wing_length_mm, beak_length_culmen_mm, tarsus_length_mm
Masked cells: 4197
Wall time: 1520.2 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE, fit scale; log for log_traits)

|trait                 |type       |   rmse|
|:---------------------|:----------|------:|
|beak_length_culmen_mm |continuous | 0.2848|
|mass_g                |continuous | 0.5965|
|tarsus_length_mm      |continuous | 0.2757|
|wing_length_mm        |continuous | 0.4676|

## Discrete-family (accuracy)

|trait             |type        | accuracy|
|:-----------------|:-----------|--------:|
|primary_lifestyle |categorical |   0.8283|
|trophic_level     |categorical |   0.8247|
|migration         |ordinal     |   0.7458|

