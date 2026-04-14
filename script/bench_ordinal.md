# Ordinal-trait benchmark

Run on: 2026-04-13 14:18:39
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 3, reps: 5, missing_frac: 0.25
Total wall: 41.8 min

## Methods

- **median**: column median imputation in latent space (rounded to nearest integer-z level).
- **baseline**: Brownian-motion baseline (Rphylopars on integer-z scale).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (varying number of ordinal levels, signal = 0.8)

### levels_3

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7761613    0.6614682
  pigauto  ord1 0.7707338    0.6619510
 baseline  ord2 0.7276893    0.6841231
  pigauto  ord2 0.7223312    0.6848615
 baseline  ord3 0.7443488    0.7096795
  pigauto  ord3 0.7431725    0.7084340
```

### levels_5

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7316746    0.6672903
  pigauto  ord1 0.7316746    0.6672903
 baseline  ord2 0.8167753    0.5815670
  pigauto  ord2 0.8116026    0.5844760
 baseline  ord3 0.7100621    0.6913265
  pigauto  ord3 0.7105481    0.6912657
```

### levels_7

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7024348    0.7337649
  pigauto  ord1 0.7008857    0.7337294
 baseline  ord2 0.7574172    0.6185209
  pigauto  ord2 0.7501433    0.6188688
 baseline  ord3 0.6995976    0.7135153
  pigauto  ord3 0.6946443    0.7133001
```

### levels_10

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7074618    0.7131729
  pigauto  ord1 0.7094305    0.7143017
 baseline  ord2 0.6876655    0.7395981
  pigauto  ord2 0.6842579    0.7366464
 baseline  ord3 0.7784372    0.6221665
  pigauto  ord3 0.7784372    0.6221665
```

## Secondary sweep (varying phylogenetic signal, n_levels = 5)

### signal_0.3

```
   method trait      rmse spearman_rho
 baseline  ord1 1.0718404    0.2277417
  pigauto  ord1 1.0513570    0.2318770
 baseline  ord2 0.9773393    0.3730459
  pigauto  ord2 0.9580944    0.3726906
 baseline  ord3 1.0184295    0.3036906
  pigauto  ord3 0.9903381    0.3022611
```

### signal_0.6

```
   method trait      rmse spearman_rho
 baseline  ord1 0.8802754    0.4843306
  pigauto  ord1 0.8570156    0.4837278
 baseline  ord2 0.8499001    0.5660927
  pigauto  ord2 0.8425191    0.5621443
 baseline  ord3 0.8795817    0.4965449
  pigauto  ord3 0.8755407    0.4968302
```

### signal_1.0

```
   method trait      rmse spearman_rho
 baseline  ord1 0.5500535    0.8339649
  pigauto  ord1 0.5500621    0.8346428
 baseline  ord2 0.5692030    0.8063849
  pigauto  ord2 0.5692030    0.8063849
 baseline  ord3 0.4939499    0.8672514
  pigauto  ord3 0.4939499    0.8672514
```

