# Ordinal-trait benchmark

Run on: 2026-04-12 06:51:53
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 3, reps: 5, missing_frac: 0.25
Total wall: 0.0 min

## Methods

- **median**: column median imputation in latent space (rounded to nearest integer-z level).
- **baseline**: Brownian-motion baseline (Rphylopars on integer-z scale).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (varying number of ordinal levels, signal = 0.8)

### levels_3

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7757250    0.6604434
  pigauto  ord1 0.7719374    0.6589602
 baseline  ord2 0.7468882    0.6588958
  pigauto  ord2 0.7388652    0.6568790
 baseline  ord3 0.7647998    0.6856596
  pigauto  ord3 0.7628889    0.6856596
```

### levels_5

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7490844    0.6568662
  pigauto  ord1 0.7490844    0.6568662
 baseline  ord2 0.8225636    0.5863699
  pigauto  ord2 0.8182460    0.5848799
 baseline  ord3 0.7228254    0.6866254
  pigauto  ord3 0.7220778    0.6854764
```

### levels_7

```
   method trait      rmse spearman_rho
 baseline  ord1 0.6982856    0.7359043
  pigauto  ord1 0.6976319    0.7365143
 baseline  ord2 0.7785534    0.6065098
  pigauto  ord2 0.7662446    0.6026244
 baseline  ord3 0.6991106    0.7072399
  pigauto  ord3 0.6943658    0.7066054
```

### levels_10

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7126037    0.7133052
  pigauto  ord1 0.7127960    0.7122916
 baseline  ord2 0.6870394    0.7370112
  pigauto  ord2 0.6846541    0.7364465
 baseline  ord3 0.7987482    0.6146538
  pigauto  ord3 0.7953441    0.6102829
```

## Secondary sweep (varying phylogenetic signal, n_levels = 5)

### signal_0.3

```
   method trait      rmse spearman_rho
 baseline  ord1 1.0777780    0.2005091
  pigauto  ord1 1.0633719    0.2026602
 baseline  ord2 0.9817494    0.3672892
  pigauto  ord2 0.9627043    0.3709890
 baseline  ord3 1.0272530    0.2951343
  pigauto  ord3 0.9871769    0.2932744
```

### signal_0.6

```
   method trait      rmse spearman_rho
 baseline  ord1 0.8795332    0.4869093
  pigauto  ord1 0.8539091    0.4895277
 baseline  ord2 0.8630739    0.5625234
  pigauto  ord2 0.8470137    0.5557060
 baseline  ord3 0.8899837    0.4907516
  pigauto  ord3 0.8853193    0.4901929
```

### signal_1.0

```
   method trait      rmse spearman_rho
 baseline  ord1 0.5551472    0.8320642
  pigauto  ord1 0.5547095    0.8337663
 baseline  ord2 0.5716622    0.8093641
  pigauto  ord2 0.5682102    0.8073447
 baseline  ord3 0.4970674    0.8688909
  pigauto  ord3 0.4970674    0.8688909
```

