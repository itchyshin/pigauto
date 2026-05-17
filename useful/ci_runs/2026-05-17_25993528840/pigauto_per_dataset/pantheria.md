# pantheria — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: body_mass_g, head_body_length_mm, gestation_d, max_longevity_m, litter_size, diet_breadth, habitat_breadth, terrestriality
Log-transformed: body_mass_g, head_body_length_mm, gestation_d, max_longevity_m
Masked cells: 2258
Wall time: 1423.5 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE, fit scale; log for log_traits)

|trait               |type       |   rmse|
|:-------------------|:----------|------:|
|body_mass_g         |continuous | 1.4305|
|gestation_d         |continuous | 0.6698|
|head_body_length_mm |continuous | 0.5770|
|max_longevity_m     |continuous | 0.8288|
|litter_size         |count      | 1.8796|

## Discrete-family (accuracy)

|trait           |type    | accuracy|
|:---------------|:-------|--------:|
|terrestriality  |binary  |   0.7018|
|diet_breadth    |ordinal |   0.3256|
|habitat_breadth |ordinal |   0.6610|

