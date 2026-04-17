# Tree uncertainty benchmark

- Species: 300 (avonet300)
- Trees: 10 posterior (BirdTree Hackett)
- Methods: single_tree (m=50) vs multi_tree (10 trees x 5 = 50)
- Downstream model: log(Mass) ~ log(Wing.Length) + log(Beak.Length_Culmen)
- Missingness fracs: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 149.5 min

## Pooled SE comparison

| method | miss_frac | term | mean_estimate | mean_SE | mean_FMI | mean_df |
|--------|-----------|------|---------------|---------|----------|---------|
| multi_tree | 0.20 | (Intercept) | -6.9852 | 0.6131 | 0.669 | 114.7 |
| single_tree | 0.20 | (Intercept) | -6.2273 | 0.6011 | 0.420 | 306.8 |
| multi_tree | 0.20 | log(Beak.Length_Culmen) | 0.5196 | 0.1758 | 0.648 | 124.9 |
| single_tree | 0.20 | log(Beak.Length_Culmen) | 0.5527 | 0.1655 | 0.465 | 256.5 |
| multi_tree | 0.20 | log(Wing.Length) | 1.9912 | 0.2143 | 0.746 | 92.9 |
| single_tree | 0.20 | log(Wing.Length) | 1.8237 | 0.1724 | 0.468 | 249.7 |
| multi_tree | 0.40 | (Intercept) | -5.8061 | 0.8071 | 0.698 | 103.4 |
| single_tree | 0.40 | (Intercept) | -5.2676 | 0.6114 | 0.371 | 372.0 |
| multi_tree | 0.40 | log(Beak.Length_Culmen) | 0.6206 | 0.1968 | 0.650 | 122.3 |
| single_tree | 0.40 | log(Beak.Length_Culmen) | 0.5729 | 0.1612 | 0.468 | 254.8 |
| multi_tree | 0.40 | log(Wing.Length) | 1.6633 | 0.2537 | 0.788 | 81.0 |
| single_tree | 0.40 | log(Wing.Length) | 1.5861 | 0.1672 | 0.471 | 255.8 |
| multi_tree | 0.60 | (Intercept) | -4.7812 | 0.8371 | 0.577 | 154.5 |
| single_tree | 0.60 | (Intercept) | -4.3014 | 0.7062 | 0.386 | 404.3 |
| multi_tree | 0.60 | log(Beak.Length_Culmen) | 0.6728 | 0.2403 | 0.719 | 98.6 |
| single_tree | 0.60 | log(Beak.Length_Culmen) | 0.8197 | 0.1675 | 0.412 | 304.9 |
| multi_tree | 0.60 | log(Wing.Length) | 1.4090 | 0.2458 | 0.721 | 95.9 |
| single_tree | 0.60 | log(Wing.Length) | 1.2229 | 0.1675 | 0.409 | 312.6 |

## SE ratio (multi_tree / single_tree)

| miss_frac | term | SE_single | SE_multi | ratio |
|-----------|------|-----------|----------|-------|
| 0.20 | (Intercept) | 0.6011 | 0.6131 | 1.02 |
| 0.20 | log(Beak.Length_Culmen) | 0.1655 | 0.1758 | 1.06 |
| 0.20 | log(Wing.Length) | 0.1724 | 0.2143 | 1.24 |
| 0.40 | (Intercept) | 0.6114 | 0.8071 | 1.32 |
| 0.40 | log(Beak.Length_Culmen) | 0.1612 | 0.1968 | 1.22 |
| 0.40 | log(Wing.Length) | 0.1672 | 0.2537 | 1.52 |
| 0.60 | (Intercept) | 0.7062 | 0.8371 | 1.19 |
| 0.60 | log(Beak.Length_Culmen) | 0.1675 | 0.2403 | 1.43 |
| 0.60 | log(Wing.Length) | 0.1675 | 0.2458 | 1.47 |

---
Generated:
2026-04-13 17:37
