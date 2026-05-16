# pigauto vs BACE — cross-dataset head-to-head

Run date: 2026-05-16 23:51:23 (UTC).
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

|dataset   |trait                 |type        |metric |better | pigauto|    bace|winner  |
|:---------|:---------------------|:-----------|:------|:------|-------:|-------:|:-------|
|amphibio  |body_mass_g           |continuous  |rmse   |lower  |   1.962|   1.731|bace    |
|amphibio  |body_size_mm          |continuous  |rmse   |lower  |   0.664|   2.353|pigauto |
|avonet    |beak_length_culmen_mm |continuous  |rmse   |lower  |   0.507|   1.772|pigauto |
|avonet    |mass_g                |continuous  |rmse   |lower  |   0.654|   2.726|pigauto |
|avonet    |migration             |ordinal     |rmse   |lower  |      NA|      NA|NA      |
|avonet    |primary_lifestyle     |categorical |acc    |higher |   0.322|   0.828|bace    |
|avonet    |tarsus_length_mm      |continuous  |rmse   |lower  |   0.424|   1.864|pigauto |
|avonet    |trophic_level         |categorical |acc    |higher |   0.377|   0.816|bace    |
|avonet    |wing_length_mm        |continuous  |rmse   |lower  |   0.559|   2.987|pigauto |
|bien      |height_m              |continuous  |rmse   |lower  |   1.928|   0.967|bace    |
|bien      |leaf_area             |continuous  |rmse   |lower  |   3.135|   5.260|pigauto |
|bien      |seed_mass             |continuous  |rmse   |lower  |   3.548|   1.978|bace    |
|bien      |sla                   |continuous  |rmse   |lower  |   7.456|   1.554|bace    |
|bien      |wood_density          |continuous  |rmse   |lower  |   0.456|   0.124|bace    |
|globtherm |elevation_max         |continuous  |rmse   |lower  | 984.383| 599.094|bace    |
|globtherm |lat_max               |continuous  |rmse   |lower  |  45.820|  29.234|bace    |
|globtherm |long_max              |continuous  |rmse   |lower  | 111.124|  74.975|bace    |
|globtherm |Tmax                  |continuous  |rmse   |lower  |   9.455|   5.200|bace    |
|globtherm |Tmin                  |continuous  |rmse   |lower  |  13.238|  14.154|pigauto |
|leptraits |flight_duration       |continuous  |rmse   |lower  |   1.091|   3.522|pigauto |
|leptraits |forewing_length_lower |continuous  |rmse   |lower  |   0.404|   0.306|bace    |
|leptraits |n_hostplant_families  |continuous  |rmse   |lower  |   0.929|   0.348|bace    |
|leptraits |wingspan_lower        |continuous  |rmse   |lower  |   0.525|   0.602|pigauto |
|pantheria |body_mass_g           |continuous  |rmse   |lower  |   1.703|   4.509|pigauto |
|pantheria |diet_breadth          |ordinal     |rmse   |lower  |      NA|      NA|NA      |
|pantheria |gestation_d           |continuous  |rmse   |lower  |   1.459|   2.766|pigauto |
|pantheria |habitat_breadth       |ordinal     |rmse   |lower  |      NA|      NA|NA      |
|pantheria |head_body_length_mm   |continuous  |rmse   |lower  |   0.791|   3.825|pigauto |
|pantheria |litter_size           |count       |rmse   |lower  |   2.629|   1.521|bace    |
|pantheria |max_longevity_m       |continuous  |rmse   |lower  |   1.309|   3.207|pigauto |
|pantheria |terrestriality        |binary      |acc    |higher |   0.683|   0.901|bace    |
