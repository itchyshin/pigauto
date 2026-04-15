# Multi-proportion trait benchmark (compositional data)

Run on: 2026-04-15 08:03:59
Species: 300, reps: 5
Total wall: 25.3 min

## Methods

- **mean**: column-mean imputation on CLR-z latent scale.
- **baseline**: per-component Brownian-motion imputation on CLR-z latent scale.
- **pigauto**: full pipeline (BM baseline + calibrated GNN) on CLR-z latent scale.

Metrics:

- `aitchison`: average Aitchison distance (Euclidean in CLR space) across held-out rows.
- `rmse_clr`: RMSE on z-scored CLR values.
- `simplex_mae`: mean absolute error on the simplex after softmax decode.

### signal_0.2

(no metrics)

### signal_0.4

(no metrics)

### signal_0.6

(no metrics)

### signal_0.8

(no metrics)

### signal_1.0

(no metrics)

### K_3

(no metrics)

### K_5

(no metrics)

### K_8

(no metrics)

### K_12

(no metrics)

