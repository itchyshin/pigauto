# Multi-proportion trait benchmark (compositional data)

Run on: 2026-04-28 07:20:44
Species: 300, reps: 3
Total wall: 35.2 min

## Methods

- **mean**: column-mean imputation on CLR-z latent scale.
- **baseline**: per-component Brownian-motion imputation on CLR-z latent scale.
- **pigauto**: full pipeline (BM baseline + calibrated GNN) on CLR-z latent scale.

Metrics:

- `aitchison`: average Aitchison distance (Euclidean in CLR space) across held-out rows.
- `rmse_clr`: RMSE on z-scored CLR values.
- `simplex_mae`: mean absolute error on the simplex after softmax decode.

### signal_0.2

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  8.819142 0.9887473   0.2737606 0.1929825
 baseline  comp  9.298765 1.0630718   0.2721827 0.2573099
  pigauto  comp  9.298765 1.0630718   0.2721827 0.2573099
```

### signal_0.4

```
   method trait aitchison rmse_clr simplex_mae  accuracy
     mean  comp  6.404743 1.018898   0.2289592 0.3333333
 baseline  comp  6.424845 1.029778   0.2074157 0.4678363
  pigauto  comp  6.424845 1.029778   0.2074157 0.4678363
```

### signal_0.6

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  4.760185 0.9897147   0.2128985 0.3157895
 baseline  comp  4.282956 0.9036986   0.1674832 0.4970760
  pigauto  comp  4.282956 0.9036986   0.1674832 0.4970760
```

### signal_0.8

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  4.125447 1.0072362   0.1922229 0.4210526
 baseline  comp  3.004510 0.7386089   0.1307477 0.5438596
  pigauto  comp  3.004510 0.7386089   0.1307477 0.5438596
```

### signal_1.0

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  3.856234 1.0036552  0.19460312 0.3216374
 baseline  comp  1.858306 0.5122988  0.08997025 0.7076023
  pigauto  comp  1.858306 0.5122988  0.08997025 0.7076023
```

### K_3

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  2.885605 0.9706635   0.2971837 0.3450292
 baseline  comp  2.622369 0.8745383   0.2573747 0.5263158
  pigauto  comp  2.622369 0.8745383   0.2573747 0.5263158
```

### K_8

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  6.878183 1.0050894   0.1522016 0.1929825
 baseline  comp  5.758896 0.8558837   0.1226024 0.3918129
  pigauto  comp  5.758896 0.8558837   0.1226024 0.3918129
```

### K_12

```
   method trait aitchison  rmse_clr simplex_mae  accuracy
     mean  comp  8.961886 1.0043540  0.10811662 0.2923977
 baseline  comp  7.771377 0.8747852  0.09378061 0.2865497
  pigauto  comp  7.771377 0.8747852  0.09378061 0.2865497
```

