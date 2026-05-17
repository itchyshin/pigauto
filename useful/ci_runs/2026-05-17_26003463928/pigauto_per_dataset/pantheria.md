# pantheria — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=20
N species used: 2000
Traits: body_mass_g, head_body_length_mm, gestation_d, max_longevity_m, litter_size, diet_breadth, habitat_breadth, terrestriality
Log-transformed: body_mass_g, head_body_length_mm, gestation_d, max_longevity_m
Masked cells: 2258
Wall time: 1432.9 s (fit) + 0.3 s (eval)

## Continuous-family (RMSE + 95% conformal coverage / width)

|trait               |type       |   rmse| coverage_95| interval_width|
|:-------------------|:----------|------:|-----------:|--------------:|
|body_mass_g         |continuous | 0.9067|      0.9682|         4.3868|
|gestation_d         |continuous | 0.3886|      0.9659|         2.0831|
|head_body_length_mm |continuous | 0.3747|      0.9585|         1.7690|
|max_longevity_m     |continuous | 0.5289|      0.9612|         2.4308|
|litter_size         |count      | 1.1285|      0.9843|         5.3736|

## Discrete-family (accuracy + Brier)

|trait           |type    | accuracy| brier|
|:---------------|:-------|--------:|-----:|
|terrestriality  |binary  |   0.9327|    NA|
|diet_breadth    |ordinal |   0.5039|    NA|
|habitat_breadth |ordinal |   0.7684|    NA|

