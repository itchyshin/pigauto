# Ordinal-trait benchmark

Run on: 2026-04-17 07:38:20
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 3, reps: 5, missing_frac: 0.25
Total wall: 24.1 min

## Methods

- **median**: column median imputation in latent space (rounded to nearest integer-z level).
- **baseline**: Brownian-motion baseline (Rphylopars on integer-z scale).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (varying number of ordinal levels, signal = 0.8)

### levels_3

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7757250    0.6604434
  pigauto  ord1 0.7733120    0.6584425
 baseline  ord2 0.7468882    0.6588958
  pigauto  ord2 0.7468882    0.6588958
 baseline  ord3 0.7647998    0.6856596
  pigauto  ord3 0.7620027    0.6834185
```

### levels_5

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7490844    0.6568662
  pigauto  ord1 0.7490844    0.6568662
 baseline  ord2 0.8225636    0.5863699
  pigauto  ord2 0.8240919    0.5785261
 baseline  ord3 0.7228254    0.6866254
  pigauto  ord3 0.7228254    0.6866254
```

### levels_7

```
   method trait      rmse spearman_rho
 baseline  ord1 0.6982856    0.7359043
  pigauto  ord1 0.6982856    0.7359043
 baseline  ord2 0.7785534    0.6065098
  pigauto  ord2 0.7769562    0.6071776
 baseline  ord3 0.6991106    0.7072399
  pigauto  ord3 0.6991106    0.7072399
```

### levels_10

```
   method trait      rmse spearman_rho
 baseline  ord1 0.7126037    0.7133052
  pigauto  ord1 0.7119527    0.7128598
 baseline  ord2 0.6870394    0.7370112
  pigauto  ord2 0.6870394    0.7370112
 baseline  ord3 0.7987482    0.6146538
  pigauto  ord3 0.7987482    0.6146538
```

## Secondary sweep (varying phylogenetic signal, n_levels = 5)

### signal_0.3

```
   method trait      rmse spearman_rho
 baseline  ord1 1.0777780    0.2005091
  pigauto  ord1 1.0722228    0.2050722
 baseline  ord2 0.9817494    0.3672892
  pigauto  ord2 0.9713463    0.3600663
 baseline  ord3 1.0272530    0.2951343
  pigauto  ord3 1.0033333    0.2920804
```

### signal_0.6

```
   method trait      rmse spearman_rho
 baseline  ord1 0.8795332    0.4869093
  pigauto  ord1 0.8672322    0.4877125
 baseline  ord2 0.8630739    0.5625234
  pigauto  ord2 0.8581154    0.5580501
 baseline  ord3 0.8899837    0.4907516
  pigauto  ord3 0.8819857    0.4934346
```

### signal_1.0

```
   method trait      rmse spearman_rho
 baseline  ord1 0.5551472    0.8320642
  pigauto  ord1 0.5535415    0.8326979
 baseline  ord2 0.5716622    0.8093641
  pigauto  ord2 0.5716622    0.8094366
 baseline  ord3 0.4970674    0.8688909
  pigauto  ord3 0.4970674    0.8688909
```

