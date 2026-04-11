# Tree uncertainty benchmark

- Species: 300 (avonet300)
- Trees: 10 posterior (BirdTree Hackett)
- Methods: single_tree (m=50) vs multi_tree (10 trees x 5 = 50)
- Downstream model: log(Mass) ~ log(Wing.Length) + log(Beak.Length_Culmen)
- Missingness fracs: 0.2, 0.4, 0.6
- Replicates: 3
- Total wall time: 149.2 min

## Pooled SE comparison

| method | miss_frac | term | mean_estimate | mean_SE | mean_FMI | mean_df |
|--------|-----------|------|---------------|---------|----------|---------|
| multi_tree | 0.20 | (Intercept) | -7.9818 | 0.2654 | 0.221 | 6275.3 |
| single_tree | 0.20 | (Intercept) | -8.0608 | 0.2366 | 0.038 | 49.0 |
| multi_tree | 0.20 | log(Beak.Length_Culmen) | 0.2938 | 0.0969 | 0.294 | 947.9 |
| single_tree | 0.20 | log(Beak.Length_Culmen) | 0.2854 | 0.0813 | 0.038 | 49.0 |
| multi_tree | 0.20 | log(Wing.Length) | 2.3602 | 0.0960 | 0.260 | 2832.3 |
| single_tree | 0.20 | log(Wing.Length) | 2.3830 | 0.0826 | 0.038 | 49.0 |
| multi_tree | 0.40 | (Intercept) | -7.9115 | 0.3038 | 0.431 | 270.1 |
| single_tree | 0.40 | (Intercept) | -8.1749 | 0.2203 | 0.038 | 49.0 |
| multi_tree | 0.40 | log(Beak.Length_Culmen) | 0.2402 | 0.1029 | 0.259 | 1752.5 |
| single_tree | 0.40 | log(Beak.Length_Culmen) | 0.0808 | 0.0855 | 0.038 | 49.0 |
| multi_tree | 0.40 | log(Wing.Length) | 2.3742 | 0.1040 | 0.310 | 952.3 |
| single_tree | 0.40 | log(Wing.Length) | 2.5378 | 0.0826 | 0.038 | 49.0 |
| multi_tree | 0.60 | (Intercept) | -8.0709 | 0.4339 | 0.710 | 110.3 |
| single_tree | 0.60 | (Intercept) | -7.8253 | 0.2313 | 0.027 | 1017765.4 |
| multi_tree | 0.60 | log(Beak.Length_Culmen) | 0.2282 | 0.1900 | 0.786 | 82.0 |
| single_tree | 0.60 | log(Beak.Length_Culmen) | 0.2110 | 0.0891 | 0.027 | 1170037.9 |
| multi_tree | 0.60 | log(Wing.Length) | 2.4146 | 0.1776 | 0.742 | 97.2 |
| single_tree | 0.60 | log(Wing.Length) | 2.3951 | 0.0894 | 0.027 | 1506952.5 |

## SE ratio (multi_tree / single_tree)

| miss_frac | term | SE_single | SE_multi | ratio |
|-----------|------|-----------|----------|-------|
| 0.20 | (Intercept) | 0.2366 | 0.2654 | 1.12 |
| 0.20 | log(Beak.Length_Culmen) | 0.0813 | 0.0969 | 1.19 |
| 0.20 | log(Wing.Length) | 0.0826 | 0.0960 | 1.16 |
| 0.40 | (Intercept) | 0.2203 | 0.3038 | 1.38 |
| 0.40 | log(Beak.Length_Culmen) | 0.0855 | 0.1029 | 1.20 |
| 0.40 | log(Wing.Length) | 0.0826 | 0.1040 | 1.26 |
| 0.60 | (Intercept) | 0.2313 | 0.4339 | 1.88 |
| 0.60 | log(Beak.Length_Culmen) | 0.0891 | 0.1900 | 2.13 |
| 0.60 | log(Wing.Length) | 0.0894 | 0.1776 | 1.99 |

---
Generated:
2026-04-10 16:36
