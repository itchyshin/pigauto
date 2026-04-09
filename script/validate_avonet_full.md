# AVONET full-dataset validation (post Fix A + Fix B)

Run on: 2026-04-09 13:20:20
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Commit:  435bddb0b39d03e22dfa8f955eff620ca2b0f35f

Species aligned to phylogeny: **9993**
Trait columns: Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length, Trophic.Level, Primary.Lifestyle, Migration
Held-out test cells: 10489

## Per-stage wall time

| Stage      | Wall (s) | Wall (min) | R heap max (MB) |
|------------|----------|------------|-----------------|
| preprocess |      0.0 |       0.00 |              62 |
| graph      |     39.9 |       0.66 |            5394 |
| splits     |      0.3 |       0.00 |            2432 |
| baseline   |     80.0 |       1.33 |            4465 |
| train      |    398.3 |       6.64 |             894 |
| predict    |      5.6 |       0.09 |             850 |
| **total** |    547.4 |       9.12 |                 |

## Held-out test-cell metrics (latent / z-score space)

| Trait | Type | n | BM RMSE | GNN RMSE | BM r | GNN r |
|---|---|---|---|---|---|---|
| Mass | continuous | 1058 | 0.151 | 0.151 | 0.989 | 0.989 |
| Beak.Length_Culmen | continuous | 1066 | 0.236 | 0.236 | 0.970 | 0.970 |
| Tarsus.Length | continuous | 1048 | 0.198 | 0.198 | 0.981 | 0.981 |
| Wing.Length | continuous | 1089 | 0.169 | 0.169 | 0.986 | 0.986 |
| Migration | ordinal | 1085 | 0.743 | 0.742 | 0.673 | 0.673 |
