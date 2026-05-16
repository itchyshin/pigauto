# pigauto vs BACE — cross-dataset head-to-head

Run date: 2026-05-16 20:35:13 (UTC).
Datasets: avonet, pantheria, amphibio, bien, globtherm, leptraits.

## Wins / ties / losses per dataset

|dataset   |winner  |  n|
|:---------|:-------|--:|
|amphibio  |pigauto |  1|
|avonet    |pigauto |  4|
|bien      |pigauto |  1|
|globtherm |pigauto |  1|
|leptraits |pigauto |  2|
|pantheria |pigauto |  4|
|amphibio  |tie     |  0|
|avonet    |tie     |  0|
|bien      |tie     |  0|
|globtherm |tie     |  0|
|leptraits |tie     |  0|
|pantheria |tie     |  0|
|amphibio  |bace    |  1|
|avonet    |bace    |  2|
|bien      |bace    |  4|
|globtherm |bace    |  4|
|leptraits |bace    |  2|
|pantheria |bace    |  2|

## Per-trait detail (medians across imputations)

|dataset   |trait                 |type        | pigauto|    bace|winner  |
|:---------|:---------------------|:-----------|-------:|-------:|:-------|
|amphibio  |body_mass_g           |continuous  |   1.962|   1.731|bace    |
|amphibio  |body_size_mm          |continuous  |   0.610|   2.353|pigauto |
|avonet    |beak_length_culmen_mm |continuous  |   0.421|   1.772|pigauto |
|avonet    |mass_g                |continuous  |   0.777|   2.726|pigauto |
|avonet    |migration             |ordinal     |      NA|      NA|NA      |
|avonet    |primary_lifestyle     |categorical |   0.350|   0.828|bace    |
|avonet    |tarsus_length_mm      |continuous  |   0.415|   1.864|pigauto |
|avonet    |trophic_level         |categorical |   0.357|   0.816|bace    |
|avonet    |wing_length_mm        |continuous  |   0.555|   2.987|pigauto |
|bien      |height_m              |continuous  |   2.111|   0.967|bace    |
|bien      |leaf_area             |continuous  |   2.669|   5.260|pigauto |
|bien      |seed_mass             |continuous  |   2.841|   1.978|bace    |
|bien      |sla                   |continuous  |   7.456|   1.554|bace    |
|bien      |wood_density          |continuous  |   0.408|   0.124|bace    |
|globtherm |elevation_max         |continuous  | 916.678| 599.094|bace    |
|globtherm |lat_max               |continuous  |  42.405|  29.234|bace    |
|globtherm |long_max              |continuous  | 113.592|  74.975|bace    |
|globtherm |Tmax                  |continuous  |   7.558|   5.200|bace    |
|globtherm |Tmin                  |continuous  |   9.281|  14.154|pigauto |
|leptraits |flight_duration       |continuous  |   1.109|   3.522|pigauto |
|leptraits |forewing_length_lower |continuous  |   0.357|   0.306|bace    |
|leptraits |n_hostplant_families  |continuous  |   0.929|   0.348|bace    |
|leptraits |wingspan_lower        |continuous  |   0.402|   0.602|pigauto |
|pantheria |body_mass_g           |continuous  |   1.254|   4.509|pigauto |
|pantheria |diet_breadth          |ordinal     |      NA|      NA|NA      |
|pantheria |gestation_d           |continuous  |   0.694|   2.766|pigauto |
|pantheria |habitat_breadth       |ordinal     |      NA|      NA|NA      |
|pantheria |head_body_length_mm   |continuous  |   0.526|   3.825|pigauto |
|pantheria |litter_size           |count       |   1.854|   1.521|bace    |
|pantheria |max_longevity_m       |continuous  |   0.856|   3.207|pigauto |
|pantheria |terrestriality        |binary      |   0.719|   0.901|bace    |
