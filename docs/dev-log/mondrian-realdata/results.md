# Mondrian real-data confirmation: results

Pre-registration: [00-preregistration.md](00-preregistration.md), Amendments 1 and 2.
Source SHA: b4d34e7
Generated: 2026-09-23 17:49:17 UTC from RESULTS_ROOT=script/mondrian_confirmation/returned

Method notes. Table 1 pools across masks (seeds) within each dataset x arm x
trait x stratum. Coverage and Winkler score pool exactly, as n_test-weighted
sums. Median half-width, and so width_ratio, pools as an n_test-weighted average
of the per-mask medians, not a literal pooled median. mcse is the paired-
difference MCSE, sqrt(mcse_mondrian^2 + mcse_split^2), each term the analytic
formula from 02_summarise_masked_confirmation.R applied to the pooled
coverage and n. Amendment 2: cond1_eligible marks traits with real_missing_frac
at least 0.05; Table 2's structured-arm far-stratum aggregates use eligible
traits only and name the traits excluded. Between-mask SD uses up to 3
masks (2 df); NA when fewer than 2 masks contributed.

## Missing or incomplete cells

- fishbase-mcar-m20260818: not run (pre-registered/amended)

## Table 1: per dataset x arm x trait x stratum

| dataset | arm | trait | stratum | n_masks | split_cov | mondrian_cov | paired_diff | mcse | width_ratio | winkler_split | winkler_mondrian | n_test_split | n_test_mondrian | n_val | n_near | n_far | fallback | real_missing_frac | cond1_eligible |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| avonet | mcar | Beak.Length_Culmen | far | 3 | 0.9385 | 0.9716 | 0.0331 | 0.0248 | 1.1918 | 59.3305 | 60.1373 | 423 | 423 | 228 | 229 | 228 | FALSE | 0.000 | FALSE |
| avonet | mcar | Beak.Length_Culmen | near | 3 | 0.9874 | 0.9644 | -0.0231 | 0.0226 | 0.7909 | 34.3101 | 33.5385 | 477 | 477 | 229 | 229 | 228 | FALSE | 0.000 | FALSE |
| avonet | mcar | Mass | far | 3 | 0.9234 | 0.9582 | 0.0348 | 0.0256 | 1.4956 | 3620.2827 | 3258.4830 | 431 | 431 | 236 | 236 | 236 | FALSE | 0.000 | FALSE |
| avonet | mcar | Mass | near | 3 | 0.9808 | 0.9467 | -0.0341 | 0.0234 | 0.6325 | 511.7989 | 530.2727 | 469 | 469 | 236 | 236 | 236 | FALSE | 0.000 | FALSE |
| avonet | mcar | Tarsus.Length | far | 3 | 0.9301 | 0.9510 | 0.0210 | 0.0262 | 2.0482 | 156.7398 | 179.1296 | 429 | 429 | 221 | 223 | 221 | FALSE | 0.000 | FALSE |
| avonet | mcar | Tarsus.Length | near | 3 | 0.9851 | 0.9766 | -0.0085 | 0.0224 | 0.8612 | 33.2919 | 29.7879 | 471 | 471 | 223 | 223 | 221 | FALSE | 0.000 | FALSE |
| avonet | mcar | Wing.Length | far | 3 | 0.9242 | 0.9747 | 0.0505 | 0.0258 | 1.8359 | 266.9360 | 305.5028 | 396 | 396 | 222 | 224 | 222 | FALSE | 0.000 | FALSE |
| avonet | mcar | Wing.Length | near | 3 | 0.9841 | 0.9563 | -0.0278 | 0.0231 | 0.8247 | 126.4058 | 122.3304 | 504 | 504 | 224 | 224 | 222 | FALSE | 0.000 | FALSE |
| fishbase | structured | DepthRangeDeep | far | 1 | 0.9466 | 0.9630 | 0.0164 | 0.0240 | 1.4218 | 4799.2094 | 4899.9397 | 487 | 487 | 235 | 236 | 235 | FALSE | 0.553 | TRUE |
| fishbase | structured | DepthRangeDeep | near | 1 | 0.9867 | 0.9601 | -0.0266 | 0.0227 | 0.6642 | 3056.2543 | 2633.8760 | 451 | 451 | 236 | 236 | 235 | FALSE | 0.553 | TRUE |
| fishbase | structured | Length | far | 1 | 0.9090 | 0.9235 | 0.0145 | 0.0186 | 1.1216 | 110.1559 | 108.5412 | 967 | 967 | 505 | 506 | 505 | FALSE | 0.024 | FALSE |
| fishbase | structured | Length | near | 1 | 0.9667 | 0.9546 | -0.0120 | 0.0160 | 0.8733 | 86.8359 | 84.6619 | 1080 | 1080 | 506 | 506 | 505 | FALSE | 0.024 | FALSE |
| fishbase | structured | Troph | far | 1 | 0.9543 | 0.9688 | 0.0146 | 0.0234 | 1.0764 | 2.2375 | 2.2845 | 481 | 481 | 240 | 241 | 240 | FALSE | 0.539 | TRUE |
| fishbase | structured | Troph | near | 1 | 0.9237 | 0.9237 | 0.0000 | 0.0261 | 1.0000 | 2.5720 | 2.5720 | 485 | 485 | 241 | 241 | 240 | FALSE | 0.539 | TRUE |
| fishbase | structured | Vulnerability | far | 1 | 0.9310 | 0.9512 | 0.0202 | 0.0172 | 1.1385 | 57.6100 | 57.4163 | 942 | 942 | 530 | 531 | 530 | FALSE | 0.000 | FALSE |
| fishbase | structured | Vulnerability | near | 1 | 0.9732 | 0.9481 | -0.0251 | 0.0156 | 0.7766 | 55.6406 | 49.7721 | 1155 | 1155 | 531 | 531 | 530 | FALSE | 0.000 | FALSE |
| fishbase | structured | Weight | far | 1 | 0.9216 | 0.9755 | 0.0539 | 0.0392 | 2.2726 | 163050.4854 | 295725.5963 | 204 | 204 | 87 | 87 | 87 | FALSE | 0.820 | TRUE |
| fishbase | structured | Weight | near | 1 | 0.9827 | 0.9827 | 0.0000 | 0.0356 | 0.6111 | 588577.8756 | 360435.2821 | 173 | 173 | 87 | 87 | 87 | FALSE | 0.820 | TRUE |
| pantheria | mcar | body_mass_g | far | 3 | 0.9278 | 0.9532 | 0.0254 | 0.0172 | 1.1954 | 4.9963 | 5.2579 | 983 | 983 | 523 | 535 | 523 | FALSE | 0.138 | TRUE |
| pantheria | mcar | body_mass_g | near | 3 | 0.9773 | 0.9655 | -0.0118 | 0.0151 | 0.7871 | 5.4553 | 4.6473 | 1102 | 1102 | 535 | 535 | 523 | FALSE | 0.138 | TRUE |
| pantheria | mcar | gestation_d | far | 3 | 0.8943 | 0.9295 | 0.0352 | 0.0303 | 1.2697 | 2.2062 | 2.0901 | 369 | 369 | 195 | 198 | 195 | FALSE | 0.664 | TRUE |
| pantheria | mcar | gestation_d | near | 3 | 0.9505 | 0.9482 | -0.0023 | 0.0263 | 0.9155 | 1.6558 | 1.6439 | 444 | 444 | 198 | 198 | 195 | FALSE | 0.664 | TRUE |
| pantheria | mcar | head_body_length_mm | far | 3 | 0.9374 | 0.9652 | 0.0278 | 0.0222 | 1.2318 | 1.9623 | 2.0425 | 575 | 575 | 282 | 287 | 282 | FALSE | 0.526 | TRUE |
| pantheria | mcar | head_body_length_mm | near | 3 | 0.9930 | 0.9720 | -0.0210 | 0.0197 | 0.7198 | 1.6341 | 1.3403 | 571 | 571 | 287 | 287 | 282 | FALSE | 0.526 | TRUE |
| pantheria | mcar | litter_size | far | 3 | 0.9494 | 0.9747 | 0.0253 | 0.0187 | 1.2350 | 5.0018 | 5.2805 | 751 | 751 | 376 | 379 | 376 | FALSE | 0.386 | TRUE |
| pantheria | mcar | litter_size | near | 3 | 0.9823 | 0.9837 | 0.0014 | 0.0172 | 0.9862 | 4.6998 | 4.6308 | 734 | 734 | 379 | 379 | 376 | FALSE | 0.386 | TRUE |
| pantheria | mcar | max_longevity_m | far | 3 | 0.9263 | 0.9391 | 0.0128 | 0.0324 | 1.0615 | 3.3359 | 3.3255 | 312 | 312 | 145 | 148 | 145 | FALSE | 0.749 | TRUE |
| pantheria | mcar | max_longevity_m | near | 3 | 0.9796 | 0.9592 | -0.0204 | 0.0289 | 0.8298 | 2.7998 | 2.4904 | 294 | 294 | 148 | 148 | 145 | FALSE | 0.749 | TRUE |
| pantheria | structured | body_mass_g | far | 3 | 0.8918 | 0.9286 | 0.0368 | 0.0174 | 1.1791 | 4.5368 | 4.3243 | 1331 | 1331 | 528 | 530 | 528 | FALSE | 0.138 | TRUE |
| pantheria | structured | body_mass_g | near | 3 | 0.9682 | 0.9390 | -0.0292 | 0.0172 | 0.8382 | 4.0256 | 3.7483 | 754 | 754 | 530 | 530 | 528 | FALSE | 0.138 | TRUE |
| pantheria | structured | gestation_d | far | 3 | 0.9591 | 0.9795 | 0.0205 | 0.0248 | 1.2741 | 2.0594 | 2.1507 | 440 | 440 | 196 | 197 | 196 | FALSE | 0.664 | TRUE |
| pantheria | structured | gestation_d | near | 3 | 0.9786 | 0.9491 | -0.0295 | 0.0258 | 0.6480 | 1.8911 | 1.6640 | 373 | 373 | 197 | 197 | 196 | FALSE | 0.664 | TRUE |
| pantheria | structured | head_body_length_mm | far | 3 | 0.9107 | 0.9217 | 0.0110 | 0.0234 | 1.1079 | 1.9127 | 1.8966 | 728 | 728 | 282 | 287 | 282 | FALSE | 0.526 | TRUE |
| pantheria | structured | head_body_length_mm | near | 3 | 0.9856 | 0.9833 | -0.0024 | 0.0200 | 1.1154 | 1.3083 | 1.3976 | 418 | 418 | 287 | 287 | 282 | FALSE | 0.526 | TRUE |
| pantheria | structured | litter_size | far | 3 | 0.9350 | 0.9458 | 0.0108 | 0.0194 | 1.0531 | 5.3791 | 5.3438 | 923 | 923 | 372 | 383 | 372 | FALSE | 0.386 | TRUE |
| pantheria | structured | litter_size | near | 3 | 0.9662 | 0.9662 | 0.0000 | 0.0191 | 0.9327 | 4.5764 | 4.5883 | 562 | 562 | 383 | 383 | 372 | FALSE | 0.386 | TRUE |
| pantheria | structured | max_longevity_m | far | 3 | 0.9452 | 0.9562 | 0.0110 | 0.0300 | 1.1775 | 3.3531 | 3.6946 | 365 | 365 | 145 | 148 | 145 | FALSE | 0.749 | TRUE |
| pantheria | structured | max_longevity_m | near | 3 | 0.9793 | 0.9710 | -0.0083 | 0.0289 | 0.8175 | 3.3202 | 3.0076 | 241 | 241 | 148 | 148 | 145 | FALSE | 0.749 | TRUE |

## Table 2: dataset-level

| dataset | median_far_gain_structured | far_gain_between_mask_sd | far_gain_n_masks | min_mondrian_far_cov | far_mincov_between_mask_sd | near_noninferiority_p | near_p_between_mask_sd | near_width_ratio | near_width_between_mask_sd | ineligible_far_traits |
|---|---|---|---|---|---|---|---|---|---|---|
| avonet | NA | NA | 0 | NA | NA | 0.7431 | 0.5381 | 0.8078 | 0.2205 | none |
| fishbase | 0.0164 | NA | 1 | 0.9630 | NA | 0.2156 | NA | 0.7766 | NA | Length; Vulnerability |
| pantheria | 0.0110 | 0.0107 | 3 | 0.9217 | 0.0183 | 0.0083 | 0.0288 | 0.8340 | 0.1074 | none |

