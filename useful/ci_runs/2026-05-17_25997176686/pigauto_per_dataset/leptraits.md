# leptraits — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: wingspan_lower, forewing_length_lower, flight_duration, n_hostplant_families
Log-transformed: wingspan_lower, forewing_length_lower, flight_duration, n_hostplant_families
Masked cells: 1080
Wall time: 784.2 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE + 95% conformal coverage / width)

|trait                 |type       |   rmse| coverage_95| interval_width|
|:---------------------|:----------|------:|-----------:|--------------:|
|flight_duration       |continuous | 0.8433|      0.9558|         3.0747|
|forewing_length_lower |continuous | 0.2772|      0.9259|         0.8870|
|n_hostplant_families  |continuous | 0.5297|      0.9799|         3.0663|
|wingspan_lower        |continuous | 0.3047|      0.9012|         1.0433|

