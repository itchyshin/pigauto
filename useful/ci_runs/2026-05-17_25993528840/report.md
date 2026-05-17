# pigauto vs BACE — cross-dataset head-to-head

Run date: 2026-05-17 14:53:49 (UTC).
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
|avonet    |bace    |  3|
|bien      |bace    |  4|
|globtherm |bace    |  4|
|leptraits |bace    |  2|
|pantheria |bace    |  4|

## Per-trait detail (medians across imputations)

|dataset   |trait                 |type        |metric |better | pigauto|    bace|winner  |
|:---------|:---------------------|:-----------|:------|:------|-------:|-------:|:-------|
|amphibio  |body_mass_g           |continuous  |rmse   |lower  |   1.996|   1.731|bace    |
|amphibio  |body_size_mm          |continuous  |rmse   |lower  |   0.594|   2.353|pigauto |
|avonet    |beak_length_culmen_mm |continuous  |rmse   |lower  |   0.461|   1.772|pigauto |
|avonet    |mass_g                |continuous  |rmse   |lower  |   0.834|   2.726|pigauto |
|avonet    |migration             |ordinal     |acc    |higher |   0.645|   0.819|bace    |
|avonet    |primary_lifestyle     |categorical |acc    |higher |   0.372|   0.828|bace    |
|avonet    |tarsus_length_mm      |continuous  |rmse   |lower  |   0.403|   1.864|pigauto |
|avonet    |trophic_level         |categorical |acc    |higher |   0.428|   0.816|bace    |
|avonet    |wing_length_mm        |continuous  |rmse   |lower  |   0.531|   2.987|pigauto |
|bien      |height_m              |continuous  |rmse   |lower  |   1.991|   0.967|bace    |
|bien      |leaf_area             |continuous  |rmse   |lower  |   2.703|   5.260|pigauto |
|bien      |seed_mass             |continuous  |rmse   |lower  |   2.596|   1.978|bace    |
|bien      |sla                   |continuous  |rmse   |lower  |   7.456|   1.554|bace    |
|bien      |wood_density          |continuous  |rmse   |lower  |   0.455|   0.124|bace    |
|globtherm |elevation_max         |continuous  |rmse   |lower  | 856.945| 599.094|bace    |
|globtherm |lat_max               |continuous  |rmse   |lower  |  42.140|  29.234|bace    |
|globtherm |long_max              |continuous  |rmse   |lower  | 114.672|  74.975|bace    |
|globtherm |Tmax                  |continuous  |rmse   |lower  |   7.580|   5.200|bace    |
|globtherm |Tmin                  |continuous  |rmse   |lower  |   9.355|  14.154|pigauto |
|leptraits |flight_duration       |continuous  |rmse   |lower  |   1.156|   3.522|pigauto |
|leptraits |forewing_length_lower |continuous  |rmse   |lower  |   0.364|   0.306|bace    |
|leptraits |n_hostplant_families  |continuous  |rmse   |lower  |   0.929|   0.348|bace    |
|leptraits |wingspan_lower        |continuous  |rmse   |lower  |   0.405|   0.602|pigauto |
|pantheria |body_mass_g           |continuous  |rmse   |lower  |   1.430|   4.509|pigauto |
|pantheria |diet_breadth          |ordinal     |acc    |higher |   0.326|   0.430|bace    |
|pantheria |gestation_d           |continuous  |rmse   |lower  |   0.670|   2.766|pigauto |
|pantheria |habitat_breadth       |ordinal     |acc    |higher |   0.661|   0.766|bace    |
|pantheria |head_body_length_mm   |continuous  |rmse   |lower  |   0.577|   3.825|pigauto |
|pantheria |litter_size           |count       |rmse   |lower  |   1.880|   1.521|bace    |
|pantheria |max_longevity_m       |continuous  |rmse   |lower  |   0.829|   3.207|pigauto |
|pantheria |terrestriality        |binary      |acc    |higher |   0.702|   0.901|bace    |
