# Continuous-trait benchmark

Run on: 2026-04-13 13:15:25
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 4, reps: 5
Total wall: 58.8 min

## Methods

- **mean**: column mean imputation (no phylogeny).
- **BM**: Brownian-motion baseline (Rphylopars).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (missing_frac = 0.25)

### BM

```
  method  trait      rmse pearson_r
      BM trait1 0.4967363 0.8663988
 pigauto trait1 0.4967363 0.8663988
      BM trait2 0.4544072 0.8864250
 pigauto trait2 0.4536278 0.8865820
      BM trait3 0.4495054 0.8731011
 pigauto trait3 0.4495054 0.8731011
      BM trait4 0.4754557 0.8779050
 pigauto trait4 0.4754557 0.8779050
```

### OU

```
  method  trait     rmse  pearson_r
      BM trait1 1.119439 0.14857626
 pigauto trait1 1.087567 0.14485961
      BM trait2 1.065164 0.27734385
 pigauto trait2 1.047763 0.27981444
      BM trait3 1.167358 0.06945645
 pigauto trait3 1.122307 0.07190917
      BM trait4 1.142055 0.13698177
 pigauto trait4 1.106968 0.13824618
```

### regime_shift

```
  method  trait      rmse pearson_r
      BM trait1 0.3850487 0.9142773
 pigauto trait1 0.3850487 0.9142773
      BM trait2 0.4171097 0.9097522
 pigauto trait2 0.4173008 0.9103828
      BM trait3 0.3314232 0.9284165
 pigauto trait3 0.3314232 0.9284165
      BM trait4 0.3033656 0.9507367
 pigauto trait4 0.3033656 0.9507367
```

### nonlinear

```
  method  trait      rmse pearson_r
      BM trait1 0.7147945 0.7189273
 pigauto trait1 0.7157085 0.7235380
      BM trait2 0.5717182 0.8187038
 pigauto trait2 0.5628034 0.8206180
      BM trait3 0.6967794 0.7731351
 pigauto trait3 0.6994407 0.7771120
      BM trait4 0.4488400 0.8748362
 pigauto trait4 0.4450811 0.8752356
```

## Secondary sweep (BM + OU, varying missingness)

### BM, missing_frac = 0.15

```
  method  trait      rmse pearson_r
      BM trait1 0.4419746 0.8919239
 pigauto trait1 0.4558371 0.8904821
      BM trait2 0.4592459 0.8989809
 pigauto trait2 0.4592459 0.8989809
      BM trait3 0.4315676 0.8861479
 pigauto trait3 0.4408927 0.8837973
      BM trait4 0.4541139 0.9021030
 pigauto trait4 0.4816927 0.8985276
```

### BM, missing_frac = 0.30

```
  method  trait      rmse pearson_r
      BM trait1 0.4904403 0.8652277
 pigauto trait1 0.4980983 0.8643317
      BM trait2 0.4787526 0.8727168
 pigauto trait2 0.4787526 0.8727168
      BM trait3 0.4707709 0.8702795
 pigauto trait3 0.4707709 0.8702795
      BM trait4 0.4836186 0.8766321
 pigauto trait4 0.4836186 0.8766321
```

### BM, missing_frac = 0.50

```
  method  trait      rmse pearson_r
      BM trait1 0.5412067 0.8427619
 pigauto trait1 0.5412067 0.8427619
      BM trait2 0.5056230 0.8462154
 pigauto trait2 0.5056230 0.8462154
      BM trait3 0.5232074 0.8428584
 pigauto trait3 0.5232074 0.8428584
      BM trait4 0.5168270 0.8454739
 pigauto trait4 0.5168270 0.8454739
```

### OU, missing_frac = 0.15

```
  method  trait      rmse  pearson_r
      BM trait1 1.0636690 0.25983298
 pigauto trait1 1.0309551 0.26062807
      BM trait2 1.0068474 0.32252923
 pigauto trait2 0.9879805 0.32762094
      BM trait3 1.2029089 0.07861814
 pigauto trait3 1.1467637 0.07936132
      BM trait4 1.1336748 0.10173831
 pigauto trait4 1.1103333 0.09856296
```

### OU, missing_frac = 0.30

```
  method  trait     rmse  pearson_r
      BM trait1 1.074144 0.15950978
 pigauto trait1 1.039044 0.15606319
      BM trait2 1.119834 0.17898792
 pigauto trait2 1.102825 0.18073749
      BM trait3 1.156460 0.09500669
 pigauto trait3 1.132107 0.09762647
      BM trait4 1.157110 0.09640980
 pigauto trait4 1.108821 0.09413957
```

### OU, missing_frac = 0.50

```
  method  trait     rmse  pearson_r
      BM trait1 1.123196 0.12727929
 pigauto trait1 1.101016 0.12348904
      BM trait2 1.112662 0.15549385
 pigauto trait2 1.091417 0.15564637
      BM trait3 1.145531 0.04548523
 pigauto trait3 1.122582 0.04909450
      BM trait4 1.151708 0.08848987
 pigauto trait4 1.114978 0.08888818
```

