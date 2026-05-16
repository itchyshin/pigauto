# Continuous-trait benchmark

Run on: 2026-05-15 13:40:42
Machine: Darwin 25.5.0 (arm64), R 4.5.2
Species: 300, traits: 4, reps: 5
Total wall: 16.2 min

## Methods

- **mean**: column mean imputation (no phylogeny).
- **BM**: Brownian-motion baseline (Rphylopars).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (missing_frac = 0.25)

### BM

```
  method  trait      rmse pearson_r
      BM trait1 0.4886809 0.8705976
 pigauto trait1 0.5066360 0.8589168
      BM trait2 0.4558267 0.8862792
 pigauto trait2 0.4621881 0.8808734
      BM trait3 0.4533913 0.8701876
 pigauto trait3 0.4533913 0.8701876
      BM trait4 0.4745146 0.8783615
 pigauto trait4 0.4509547 0.8865646
```

### OU

```
  method  trait     rmse   pearson_r
      BM trait1 1.134902  0.13963660
 pigauto trait1 1.274152 -0.15036545
      BM trait2 1.111751  0.23476098
 pigauto trait2 1.000936  0.38545108
      BM trait3 1.185897  0.05873885
      BM trait4 1.159863  0.12775233
 pigauto trait4 1.303441  0.03421722
```

### regime_shift

```
  method  trait      rmse pearson_r
      BM trait1 0.3842905 0.9152567
 pigauto trait1 0.3906924 0.9094763
      BM trait2 0.4204941 0.9082347
 pigauto trait2 0.4204941 0.9082347
      BM trait3 0.3414650 0.9255868
 pigauto trait3 0.3414650 0.9255868
      BM trait4 0.3066973 0.9492055
 pigauto trait4 0.3306409 0.9429195
```

### nonlinear

```
  method  trait      rmse pearson_r
      BM trait1 0.7208384 0.7191276
 pigauto trait1 0.7282694 0.7055843
      BM trait2 0.5528077 0.8333589
 pigauto trait2 0.6163727 0.8081554
      BM trait3 0.6565790 0.7959599
 pigauto trait3 0.6925087 0.7760900
      BM trait4 0.4504407 0.8735342
 pigauto trait4 0.4705268 0.8555554
```

## Secondary sweep (BM + OU, varying missingness)

### BM, missing_frac = 0.15

```
  method  trait      rmse pearson_r
      BM trait1 0.4439541 0.8908556
 pigauto trait1 0.4439541 0.8908556
      BM trait2 0.4558421 0.8991569
 pigauto trait2 0.4558421 0.8991569
      BM trait3 0.4318808 0.8854750
 pigauto trait3 0.4318808 0.8854750
      BM trait4 0.4478773 0.9047288
 pigauto trait4 0.5191595 0.8616040
```

### BM, missing_frac = 0.30

```
  method  trait      rmse pearson_r
      BM trait1 0.4877209 0.8669367
 pigauto trait1 0.5106114 0.8520637
      BM trait2 0.4864311 0.8683894
 pigauto trait2 0.4970997 0.8651051
      BM trait3 0.4754692 0.8672855
 pigauto trait3 0.4754692 0.8672855
      BM trait4 0.4821998 0.8782160
 pigauto trait4 0.4868573 0.8762852
```

### BM, missing_frac = 0.50

```
  method  trait      rmse pearson_r
      BM trait1 0.5461826 0.8388566
 pigauto trait1 0.5485897 0.8388566
      BM trait2 0.5051078 0.8476373
 pigauto trait2 0.5051078 0.8476373
      BM trait3 0.5265010 0.8414441
 pigauto trait3 0.5265010 0.8414441
      BM trait4 0.5207890 0.8435174
 pigauto trait4 0.5316908 0.8370655
```

### OU, missing_frac = 0.15

```
  method  trait     rmse  pearson_r
      BM trait1 1.068858 0.25996846
 pigauto trait1 1.068553 0.23821339
      BM trait2 1.033483 0.29242276
 pigauto trait2 1.199582 0.16121700
      BM trait3 1.204560 0.07539891
      BM trait4 1.147486 0.09273750
 pigauto trait4 1.071099 0.21997331
```

### OU, missing_frac = 0.30

```
  method  trait     rmse   pearson_r
      BM trait1 1.107821  0.12963658
 pigauto trait1 1.114486 -0.18121881
      BM trait2 1.155393  0.15163548
 pigauto trait2 1.021820  0.23044067
      BM trait3 1.169799  0.08961219
      BM trait4 1.156355  0.10404563
```

### OU, missing_frac = 0.50

```
  method  trait     rmse    pearson_r
      BM trait1 1.138535  0.141254888
 pigauto trait1 1.094618  0.003166585
      BM trait2 1.164806  0.119845142
      BM trait3 1.164028  0.037689113
      BM trait4 1.178730  0.085835027
 pigauto trait4 1.174580 -0.014931206
```

