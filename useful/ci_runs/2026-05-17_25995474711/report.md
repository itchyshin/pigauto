# pigauto vs BACE — cross-dataset head-to-head

Run date: 2026-05-17 16:19:18 (UTC).
Datasets: avonet, pantheria, amphibio, bien, globtherm, leptraits.

## Wins / ties / losses per dataset

|dataset   |winner  |  n|
|:---------|:-------|--:|
|amphibio  |pigauto |  2|
|avonet    |pigauto |  6|
|bien      |pigauto |  1|
|globtherm |pigauto |  4|
|leptraits |pigauto |  3|
|pantheria |pigauto |  8|
|amphibio  |tie     |  0|
|avonet    |tie     |  0|
|bien      |tie     |  0|
|globtherm |tie     |  0|
|leptraits |tie     |  0|
|pantheria |tie     |  0|
|amphibio  |bace    |  0|
|avonet    |bace    |  1|
|bien      |bace    |  4|
|globtherm |bace    |  1|
|leptraits |bace    |  1|
|pantheria |bace    |  0|

## Per-trait detail (medians across imputations)

|dataset   |trait                 |type        |metric |better | pigauto|    bace|winner  |
|:---------|:---------------------|:-----------|:------|:------|-------:|-------:|:-------|
|amphibio  |body_mass_g           |continuous  |rmse   |lower  |   1.354|   1.731|pigauto |
|amphibio  |body_size_mm          |continuous  |rmse   |lower  |   0.410|   2.353|pigauto |
|avonet    |beak_length_culmen_mm |continuous  |rmse   |lower  |   0.285|   1.772|pigauto |
|avonet    |mass_g                |continuous  |rmse   |lower  |   0.597|   2.726|pigauto |
|avonet    |migration             |ordinal     |acc    |higher |   0.746|   0.819|bace    |
|avonet    |primary_lifestyle     |categorical |acc    |higher |   0.828|   0.828|pigauto |
|avonet    |tarsus_length_mm      |continuous  |rmse   |lower  |   0.276|   1.864|pigauto |
|avonet    |trophic_level         |categorical |acc    |higher |   0.825|   0.816|pigauto |
|avonet    |wing_length_mm        |continuous  |rmse   |lower  |   0.468|   2.987|pigauto |
|bien      |height_m              |continuous  |rmse   |lower  |   1.637|   0.967|bace    |
|bien      |leaf_area             |continuous  |rmse   |lower  |   1.738|   5.260|pigauto |
|bien      |seed_mass             |continuous  |rmse   |lower  |   2.104|   1.978|bace    |
|bien      |sla                   |continuous  |rmse   |lower  |   2.557|   1.554|bace    |
|bien      |wood_density          |continuous  |rmse   |lower  |   0.245|   0.124|bace    |
|globtherm |elevation_max         |continuous  |rmse   |lower  | 625.773| 599.094|bace    |
|globtherm |lat_max               |continuous  |rmse   |lower  |  28.309|  29.234|pigauto |
|globtherm |long_max              |continuous  |rmse   |lower  |  73.830|  74.975|pigauto |
|globtherm |Tmax                  |continuous  |rmse   |lower  |   4.659|   5.200|pigauto |
|globtherm |Tmin                  |continuous  |rmse   |lower  |   6.330|  14.154|pigauto |
|leptraits |flight_duration       |continuous  |rmse   |lower  |   0.843|   3.522|pigauto |
|leptraits |forewing_length_lower |continuous  |rmse   |lower  |   0.277|   0.306|pigauto |
|leptraits |n_hostplant_families  |continuous  |rmse   |lower  |   0.530|   0.348|bace    |
|leptraits |wingspan_lower        |continuous  |rmse   |lower  |   0.305|   0.602|pigauto |
|pantheria |body_mass_g           |continuous  |rmse   |lower  |   0.907|   4.509|pigauto |
|pantheria |diet_breadth          |ordinal     |acc    |higher |   0.504|   0.430|pigauto |
|pantheria |gestation_d           |continuous  |rmse   |lower  |   0.389|   2.766|pigauto |
|pantheria |habitat_breadth       |ordinal     |acc    |higher |   0.768|   0.766|pigauto |
|pantheria |head_body_length_mm   |continuous  |rmse   |lower  |   0.375|   3.825|pigauto |
|pantheria |litter_size           |count       |rmse   |lower  |   1.129|   1.521|pigauto |
|pantheria |max_longevity_m       |continuous  |rmse   |lower  |   0.529|   3.207|pigauto |
|pantheria |terrestriality        |binary      |acc    |higher |   0.933|   0.901|pigauto |
