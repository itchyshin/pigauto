# Categorical-trait benchmark

Run on: 2026-04-29 17:15:32
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 2, reps: 5, missing_frac: 0.25
Total wall: 15.2 min

## Methods

- **mode**: column-mean imputation in latent (one-hot) space; argmax determines predicted category.
- **baseline**: phylogenetic label propagation.
- **pigauto**: full pipeline (label-propagation baseline + calibrated GNN).

## Primary sweep (varying K, signal = 0.8)

### K_3

```
   method trait  accuracy
     mode  cat1 0.5344526
 baseline  cat1 0.7526704
  pigauto  cat1 0.7026704
     mode  cat2 0.5015637
 baseline  cat2 0.7396138
  pigauto  cat2 0.6777956
```

### K_5

```
   method trait  accuracy
     mode  cat1 0.4089613
 baseline  cat1 0.6440265
  pigauto  cat1 0.6440265
     mode  cat2 0.3097127
 baseline  cat2 0.6272433
  pigauto  cat2 0.6272433
```

### K_8

```
   method trait  accuracy
     mode  cat1 0.2885468
 baseline  cat1 0.5600532
  pigauto  cat1 0.5600532
     mode  cat2 0.2800575
 baseline  cat2 0.5885368
  pigauto  cat2 0.5569578
```

### K_12

```
   method trait  accuracy
     mode  cat1 0.2270906
 baseline  cat1 0.5272010
  pigauto  cat1 0.5272010
     mode  cat2 0.2377294
 baseline  cat2 0.5535055
  pigauto  cat2 0.5535055
```

## Secondary sweep (K = 5, varying signal)

### signal_0.3

```
   method trait  accuracy
     mode  cat1 0.3313959
 baseline  cat1 0.4673075
  pigauto  cat1 0.4066974
     mode  cat2 0.4491975
 baseline  cat2 0.5680947
  pigauto  cat2 0.5746521
```

### signal_0.6

```
   method trait  accuracy
     mode  cat1 0.2364399
 baseline  cat1 0.6149180
  pigauto  cat1 0.6149180
     mode  cat2 0.3907111
 baseline  cat2 0.6624808
  pigauto  cat2 0.6624808
```

### signal_1.0

```
   method trait  accuracy
     mode  cat1 0.3144947
 baseline  cat1 0.6924967
  pigauto  cat1 0.6653781
     mode  cat2 0.3794363
 baseline  cat2 0.6971368
  pigauto  cat2 0.6971368
```

## Grand summary (accuracy averaged across traits and reps)

```
   method   scenario  accuracy
     mode       K_12 0.2324100
 baseline       K_12 0.5403533
  pigauto       K_12 0.5403533
     mode        K_3 0.5180082
 baseline        K_3 0.7461421
  pigauto        K_3 0.6902330
     mode        K_5 0.3593370
 baseline        K_5 0.6356349
  pigauto        K_5 0.6356349
     mode        K_8 0.2843022
 baseline        K_8 0.5742950
  pigauto        K_8 0.5585055
     mode signal_0.3 0.3902967
 baseline signal_0.3 0.5177011
  pigauto signal_0.3 0.4906748
     mode signal_0.6 0.3135755
 baseline signal_0.6 0.6386994
  pigauto signal_0.6 0.6386994
     mode signal_1.0 0.3469655
 baseline signal_1.0 0.6948167
  pigauto signal_1.0 0.6812574
```

