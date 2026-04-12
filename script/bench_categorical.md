# Categorical-trait benchmark

Run on: 2026-04-12 06:51:51
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 2, reps: 5, missing_frac: 0.25
Total wall: 0.0 min

## Methods

- **mode**: column-mean imputation in latent (one-hot) space; argmax determines predicted category.
- **baseline**: phylogenetic label propagation.
- **pigauto**: full pipeline (label-propagation baseline + calibrated GNN).

## Primary sweep (varying K, signal = 0.8)

### K_3

```
   method trait  accuracy
     mode  cat1 0.5344526
 baseline  cat1 0.6223499
  pigauto  cat1 0.6223499
     mode  cat2 0.5015637
 baseline  cat2 0.6678702
  pigauto  cat2 0.6678702
```

### K_5

```
   method trait  accuracy
     mode  cat1 0.4089613
 baseline  cat1 0.5043528
  pigauto  cat1 0.5043528
     mode  cat2 0.3097127
 baseline  cat2 0.4626480
  pigauto  cat2 0.4013577
```

### K_8

```
   method trait  accuracy
     mode  cat1 0.2885468
 baseline  cat1 0.4420211
  pigauto  cat1 0.4420211
     mode  cat2 0.2800575
 baseline  cat2 0.4055132
  pigauto  cat2 0.4055132
```

### K_12

```
   method trait  accuracy
     mode  cat1 0.2270906
 baseline  cat1 0.3794704
  pigauto  cat1 0.3794704
     mode  cat2 0.2377294
 baseline  cat2 0.3810834
  pigauto  cat2 0.3810834
```

## Secondary sweep (K = 5, varying signal)

### signal_0.3

```
   method trait  accuracy
     mode  cat1 0.3313959
 baseline  cat1 0.4182690
  pigauto  cat1 0.4182690
     mode  cat2 0.4491975
 baseline  cat2 0.5048619
  pigauto  cat2 0.5264003
```

### signal_0.6

```
   method trait  accuracy
     mode  cat1 0.2364399
 baseline  cat1 0.4247620
  pigauto  cat1 0.4247620
     mode  cat2 0.3907111
 baseline  cat2 0.5465322
  pigauto  cat2 0.5465322
```

### signal_1.0

```
   method trait  accuracy
     mode  cat1 0.3144947
 baseline  cat1 0.5584643
  pigauto  cat1 0.5584643
     mode  cat2 0.3794363
 baseline  cat2 0.5344037
  pigauto  cat2 0.5344037
```

## Grand summary (accuracy averaged across traits and reps)

```
   method   scenario  accuracy
     mode       K_12 0.2324100
 baseline       K_12 0.3802769
  pigauto       K_12 0.3802769
     mode        K_3 0.5180082
 baseline        K_3 0.6451101
  pigauto        K_3 0.6451101
     mode        K_5 0.3593370
 baseline        K_5 0.4835004
  pigauto        K_5 0.4528552
     mode        K_8 0.2843022
 baseline        K_8 0.4237672
  pigauto        K_8 0.4237672
     mode signal_0.3 0.3902967
 baseline signal_0.3 0.4615655
  pigauto signal_0.3 0.4723347
     mode signal_0.6 0.3135755
 baseline signal_0.6 0.4856471
  pigauto signal_0.6 0.4856471
     mode signal_1.0 0.3469655
 baseline signal_1.0 0.5464340
  pigauto signal_1.0 0.5464340
```

