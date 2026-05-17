# leptraits — pigauto CI bench (BACE-compat)
Seed=2026, subset_n=2000, missing_frac=0.30, n_imputations=10
N species used: 2000
Traits: wingspan_lower, forewing_length_lower, flight_duration, n_hostplant_families
Log-transformed: wingspan_lower, forewing_length_lower, flight_duration, n_hostplant_families
Masked cells: 1080
Wall time: 768.4 s (fit) + 0.1 s (eval)

## Continuous-family (RMSE, fit scale; log for log_traits)

|trait                 |type       |   rmse|
|:---------------------|:----------|------:|
|flight_duration       |continuous | 0.8433|
|forewing_length_lower |continuous | 0.2772|
|n_hostplant_families  |continuous | 0.5297|
|wingspan_lower        |continuous | 0.3047|

