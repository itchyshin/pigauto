# pigauto vs BACE — cross-dataset head-to-head

Run date: 2026-05-17 17:30:29 (UTC).
Datasets: avonet, pantheria, amphibio, bien, globtherm, leptraits.

## Wins / ties / losses per dataset

|dataset   |winner  |  n|
|:---------|:-------|--:|
|amphibio  |pigauto |  2|
|avonet    |pigauto |  4|
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
|avonet    |bace    |  3|
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
|avonet    |mass_g                |continuous  |rmse   |lower  |   0.582|   2.726|pigauto |
|avonet    |migration             |ordinal     |acc    |higher |   0.746|   0.819|bace    |
|avonet    |primary_lifestyle     |categorical |acc    |higher |   0.657|   0.828|bace    |
|avonet    |tarsus_length_mm      |continuous  |rmse   |lower  |   0.276|   1.864|pigauto |
|avonet    |trophic_level         |categorical |acc    |higher |   0.731|   0.816|bace    |
|avonet    |wing_length_mm        |continuous  |rmse   |lower  |   0.467|   2.987|pigauto |
|bien      |height_m              |continuous  |rmse   |lower  |   1.642|   0.967|bace    |
|bien      |leaf_area             |continuous  |rmse   |lower  |   1.739|   5.260|pigauto |
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

## Brier score (lower is better)

|   |dataset   |trait             |type        | brier_pigauto| brier_bace|brier_winner |
|:--|:---------|:-----------------|:-----------|-------------:|----------:|:------------|
|5  |avonet    |migration         |ordinal     |            NA|      0.299|NA           |
|6  |avonet    |primary_lifestyle |categorical |         0.430|      0.318|bace         |
|8  |avonet    |trophic_level     |categorical |         0.405|      0.338|bace         |
|25 |pantheria |diet_breadth      |ordinal     |            NA|      0.727|NA           |
|27 |pantheria |habitat_breadth   |ordinal     |            NA|      0.366|NA           |
|31 |pantheria |terrestriality    |binary      |            NA|      0.169|NA           |

## Pigauto 95% conformal coverage + interval width (BACE has none)

Target coverage = 0.95. Coverage close to target with smaller
interval width indicates well-calibrated, tight predictions.

|   |dataset   |trait                 |type       | coverage_95| interval_width|
|:--|:---------|:---------------------|:----------|-----------:|--------------:|
|1  |amphibio  |body_mass_g           |continuous |       1.000|          6.216|
|2  |amphibio  |body_size_mm          |continuous |       0.946|          1.660|
|3  |avonet    |beak_length_culmen_mm |continuous |       0.967|          1.412|
|4  |avonet    |mass_g                |continuous |       0.955|          2.376|
|7  |avonet    |tarsus_length_mm      |continuous |       0.957|          1.151|
|9  |avonet    |wing_length_mm        |continuous |       0.960|          0.991|
|10 |bien      |height_m              |continuous |       0.828|          4.384|
|11 |bien      |leaf_area             |continuous |       0.950|          7.945|
|12 |bien      |seed_mass             |continuous |       0.887|          6.045|
|13 |bien      |sla                   |continuous |       1.000|         26.988|
|14 |bien      |wood_density          |continuous |       0.996|          1.505|
|15 |globtherm |elevation_max         |continuous |       0.926|       2259.657|
|16 |globtherm |lat_max               |continuous |       0.944|        123.323|
|17 |globtherm |long_max              |continuous |       0.928|        341.486|
|18 |globtherm |Tmax                  |continuous |       0.967|         23.684|
|19 |globtherm |Tmin                  |continuous |       0.961|         27.329|
|20 |leptraits |flight_duration       |continuous |       0.956|          3.075|
|21 |leptraits |forewing_length_lower |continuous |       0.926|          0.887|
|22 |leptraits |n_hostplant_families  |continuous |       0.980|          3.066|
|23 |leptraits |wingspan_lower        |continuous |       0.901|          1.043|
|24 |pantheria |body_mass_g           |continuous |       0.968|          4.387|
|26 |pantheria |gestation_d           |continuous |       0.966|          2.083|
|28 |pantheria |head_body_length_mm   |continuous |       0.959|          1.769|
|29 |pantheria |litter_size           |count      |       0.984|          5.374|
|30 |pantheria |max_longevity_m       |continuous |       0.961|          2.431|
