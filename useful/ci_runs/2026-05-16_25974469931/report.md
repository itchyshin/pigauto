# pigauto vs BACE — cross-dataset head-to-head

Run date: 2026-05-16 22:22:26 (UTC).
Datasets: avonet, pantheria, amphibio, bien, globtherm, leptraits.

## Wins / ties / losses per dataset

|dataset   |winner  |  n|
|:---------|:-------|--:|
|amphibio  |pigauto |  0|
|avonet    |pigauto |  0|
|bien      |pigauto |  0|
|globtherm |pigauto |  0|
|leptraits |pigauto |  0|
|pantheria |pigauto |  0|
|amphibio  |tie     |  0|
|avonet    |tie     |  0|
|bien      |tie     |  0|
|globtherm |tie     |  0|
|leptraits |tie     |  0|
|pantheria |tie     |  0|
|amphibio  |bace    |  0|
|avonet    |bace    |  0|
|bien      |bace    |  0|
|globtherm |bace    |  0|
|leptraits |bace    |  0|
|pantheria |bace    |  0|

## Per-trait detail (medians across imputations)

|dataset   |trait                 |type        |metric |better | pigauto|    bace|winner |
|:---------|:---------------------|:-----------|:------|:------|-------:|-------:|:------|
|amphibio  |body_mass_g           |continuous  |rmse   |lower  |      NA|   1.731|NA     |
|amphibio  |body_size_mm          |continuous  |rmse   |lower  |      NA|   2.353|NA     |
|avonet    |beak_length_culmen_mm |continuous  |rmse   |lower  |      NA|   1.772|NA     |
|avonet    |mass_g                |continuous  |rmse   |lower  |      NA|   2.726|NA     |
|avonet    |migration             |ordinal     |rmse   |lower  |      NA|      NA|NA     |
|avonet    |primary_lifestyle     |categorical |acc    |higher |      NA|   0.828|NA     |
|avonet    |tarsus_length_mm      |continuous  |rmse   |lower  |      NA|   1.864|NA     |
|avonet    |trophic_level         |categorical |acc    |higher |      NA|   0.816|NA     |
|avonet    |wing_length_mm        |continuous  |rmse   |lower  |      NA|   2.987|NA     |
|bien      |height_m              |continuous  |rmse   |lower  |      NA|   0.967|NA     |
|bien      |leaf_area             |continuous  |rmse   |lower  |      NA|   5.260|NA     |
|bien      |seed_mass             |continuous  |rmse   |lower  |      NA|   1.978|NA     |
|bien      |sla                   |continuous  |rmse   |lower  |      NA|   1.554|NA     |
|bien      |wood_density          |continuous  |rmse   |lower  |      NA|   0.124|NA     |
|globtherm |elevation_max         |continuous  |rmse   |lower  |      NA| 599.094|NA     |
|globtherm |lat_max               |continuous  |rmse   |lower  |      NA|  29.234|NA     |
|globtherm |long_max              |continuous  |rmse   |lower  |      NA|  74.975|NA     |
|globtherm |Tmax                  |continuous  |rmse   |lower  |      NA|   5.200|NA     |
|globtherm |Tmin                  |continuous  |rmse   |lower  |      NA|  14.154|NA     |
|leptraits |flight_duration       |continuous  |rmse   |lower  |      NA|   3.522|NA     |
|leptraits |forewing_length_lower |continuous  |rmse   |lower  |      NA|   0.306|NA     |
|leptraits |n_hostplant_families  |continuous  |rmse   |lower  |      NA|   0.348|NA     |
|leptraits |wingspan_lower        |continuous  |rmse   |lower  |      NA|   0.602|NA     |
|pantheria |body_mass_g           |continuous  |rmse   |lower  |      NA|   4.509|NA     |
|pantheria |diet_breadth          |ordinal     |rmse   |lower  |      NA|      NA|NA     |
|pantheria |gestation_d           |continuous  |rmse   |lower  |      NA|   2.766|NA     |
|pantheria |habitat_breadth       |ordinal     |rmse   |lower  |      NA|      NA|NA     |
|pantheria |head_body_length_mm   |continuous  |rmse   |lower  |      NA|   3.825|NA     |
|pantheria |litter_size           |count       |rmse   |lower  |      NA|   1.521|NA     |
|pantheria |max_longevity_m       |continuous  |rmse   |lower  |      NA|   3.207|NA     |
|pantheria |terrestriality        |binary      |acc    |higher |      NA|   0.901|NA     |
