# pantheria — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: body_mass_g, head_body_length_mm, gestation_d, max_longevity_m, litter_size, diet_breadth, habitat_breadth, terrestriality
Log-transformed: body_mass_g, head_body_length_mm, gestation_d, max_longevity_m
Masked cells: 2258
Wall time: 1497.3 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE, fit scale; log for log_traits)

|trait               |type       |   rmse|
|:-------------------|:----------|------:|
|body_mass_g         |continuous | 1.2545|
|gestation_d         |continuous | 0.6940|
|head_body_length_mm |continuous | 0.5257|
|max_longevity_m     |continuous | 0.8563|
|litter_size         |count      | 1.8544|

## Discrete-family (accuracy)

|trait           |type    | accuracy|
|:---------------|:-------|--------:|
|terrestriality  |binary  |   0.7193|
|diet_breadth    |ordinal |   0.3256|
|habitat_breadth |ordinal |   0.6610|

