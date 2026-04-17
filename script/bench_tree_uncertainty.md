# Tree uncertainty benchmark

- Species: 300 (avonet300)
- Trees: 10 posterior (BirdTree Hackett)
- Methods: single_tree (m=50) vs multi_tree (10 trees x 5 = 50)
- Downstream model: log(Mass) ~ log(Wing.Length) + log(Beak.Length_Culmen)
- Missingness fracs: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 279.4 min

## Pooled SE comparison

| method | miss_frac | term | mean_estimate | mean_SE | mean_FMI | mean_df |
|--------|-----------|------|---------------|---------|----------|---------|
| multi_tree | 0.20 | (Intercept) | -7.7447 | 0.3106 | 0.311 | 544.9 |
| single_tree | 0.20 | (Intercept) | -7.7734 | 0.4264 | 0.293 | 680.1 |
| multi_tree | 0.20 | log(Beak.Length_Culmen) | 0.3504 | 0.1103 | 0.389 | 348.2 |
| single_tree | 0.20 | log(Beak.Length_Culmen) | 0.3313 | 0.1301 | 0.322 | 512.9 |
| multi_tree | 0.20 | log(Wing.Length) | 2.2720 | 0.1115 | 0.389 | 339.2 |
| single_tree | 0.20 | log(Wing.Length) | 2.2892 | 0.1356 | 0.301 | 602.9 |
| multi_tree | 0.40 | (Intercept) | -7.2205 | 0.4985 | 0.628 | 136.7 |
| single_tree | 0.40 | (Intercept) | -7.3682 | 0.4327 | 0.368 | 379.4 |
| multi_tree | 0.40 | log(Beak.Length_Culmen) | 0.4260 | 0.1541 | 0.567 | 162.8 |
| single_tree | 0.40 | log(Beak.Length_Culmen) | 0.3292 | 0.1492 | 0.476 | 224.3 |
| multi_tree | 0.40 | log(Wing.Length) | 2.1008 | 0.1642 | 0.610 | 147.8 |
| single_tree | 0.40 | log(Wing.Length) | 2.1998 | 0.1470 | 0.451 | 247.3 |
| multi_tree | 0.60 | (Intercept) | -6.8702 | 0.7352 | 0.760 | 88.0 |
| single_tree | 0.60 | (Intercept) | -7.0156 | 0.4859 | 0.410 | 327.5 |
| multi_tree | 0.60 | log(Beak.Length_Culmen) | 0.5650 | 0.2022 | 0.687 | 112.0 |
| single_tree | 0.60 | log(Beak.Length_Culmen) | 0.4547 | 0.1382 | 0.405 | 308.1 |
| multi_tree | 0.60 | log(Wing.Length) | 1.9386 | 0.2334 | 0.753 | 91.7 |
| single_tree | 0.60 | log(Wing.Length) | 2.0410 | 0.1474 | 0.392 | 341.5 |

## SE ratio (multi_tree / single_tree)

| miss_frac | term | SE_single | SE_multi | ratio |
|-----------|------|-----------|----------|-------|
| 0.20 | (Intercept) | 0.4264 | 0.3106 | 0.73 |
| 0.20 | log(Beak.Length_Culmen) | 0.1301 | 0.1103 | 0.85 |
| 0.20 | log(Wing.Length) | 0.1356 | 0.1115 | 0.82 |
| 0.40 | (Intercept) | 0.4327 | 0.4985 | 1.15 |
| 0.40 | log(Beak.Length_Culmen) | 0.1492 | 0.1541 | 1.03 |
| 0.40 | log(Wing.Length) | 0.1470 | 0.1642 | 1.12 |
| 0.60 | (Intercept) | 0.4859 | 0.7352 | 1.51 |
| 0.60 | log(Beak.Length_Culmen) | 0.1382 | 0.2022 | 1.46 |
| 0.60 | log(Wing.Length) | 0.1474 | 0.2334 | 1.58 |

---
Generated:
2026-04-17 16:01
