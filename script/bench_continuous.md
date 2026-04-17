# Continuous-trait benchmark

Run on: 2026-04-17 08:46:48
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 4, reps: 5
Total wall: 0.0 min

## Methods

- **mean**: column mean imputation (no phylogeny).
- **BM**: Brownian-motion baseline (Rphylopars).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (missing_frac = 0.25)

### BM

```
  method  trait      rmse pearson_r
      BM trait1 0.4886809 0.8705976
 pigauto trait1 0.4949629 0.8700480
      BM trait2 0.4558267 0.8862792
 pigauto trait2 0.4534462 0.8866550
      BM trait3 0.4533913 0.8701876
 pigauto trait3 0.4533913 0.8701876
      BM trait4 0.4745146 0.8783615
 pigauto trait4 0.4763183 0.8784153
```

### OU

```
  method  trait     rmse  pearson_r
      BM trait1 1.134902 0.13963660
 pigauto trait1 1.116960 0.13025472
      BM trait2 1.111751 0.23476098
 pigauto trait2 1.103693 0.23879219
      BM trait3 1.185897 0.05873885
 pigauto trait3 1.147709 0.06064279
      BM trait4 1.159863 0.12775233
 pigauto trait4 1.120844 0.13638028
```

### regime_shift

```
  method  trait      rmse pearson_r
      BM trait1 0.3842905 0.9152567
 pigauto trait1 0.3842905 0.9152567
      BM trait2 0.4204941 0.9082347
 pigauto trait2 0.4197538 0.9086790
      BM trait3 0.3414650 0.9255868
 pigauto trait3 0.3414650 0.9255868
      BM trait4 0.3066973 0.9492055
 pigauto trait4 0.3236757 0.9487930
```

### nonlinear

```
  method  trait      rmse pearson_r
      BM trait1 0.7208384 0.7191276
 pigauto trait1 0.7208384 0.7191276
      BM trait2 0.5528077 0.8333589
 pigauto trait2 0.5457424 0.8342428
      BM trait3 0.6565790 0.7959599
 pigauto trait3 0.6474180 0.8009217
      BM trait4 0.4504407 0.8735342
 pigauto trait4 0.4468838 0.8725947
```

## Secondary sweep (BM + OU, varying missingness)

### BM, missing_frac = 0.15

```
  method  trait      rmse pearson_r
      BM trait1 0.4439541 0.8908556
 pigauto trait1 0.4588447 0.8907375
      BM trait2 0.4558421 0.8991569
 pigauto trait2 0.4558421 0.8991569
      BM trait3 0.4318808 0.8854750
 pigauto trait3 0.4318808 0.8854750
      BM trait4 0.4478773 0.9047288
 pigauto trait4 0.4627135 0.9036003
```

### BM, missing_frac = 0.30

```
  method  trait      rmse pearson_r
      BM trait1 0.4877209 0.8669367
 pigauto trait1 0.4877209 0.8669367
      BM trait2 0.4864311 0.8683894
 pigauto trait2 0.4864311 0.8683894
      BM trait3 0.4754692 0.8672855
 pigauto trait3 0.4754692 0.8672855
      BM trait4 0.4821998 0.8782160
 pigauto trait4 0.4964750 0.8781354
```

### BM, missing_frac = 0.50

```
  method  trait      rmse pearson_r
      BM trait1 0.5461826 0.8388566
 pigauto trait1 0.5461826 0.8388566
      BM trait2 0.5051078 0.8476373
 pigauto trait2 0.5040534 0.8479441
      BM trait3 0.5265010 0.8414441
 pigauto trait3 0.5265010 0.8414441
      BM trait4 0.5207890 0.8435174
 pigauto trait4 0.5207890 0.8435174
```

### OU, missing_frac = 0.15

```
  method  trait     rmse  pearson_r
      BM trait1 1.068858 0.25996846
 pigauto trait1 1.037773 0.25215555
      BM trait2 1.033483 0.29242276
 pigauto trait2 1.021008 0.29800038
      BM trait3 1.204560 0.07539891
 pigauto trait3 1.164650 0.07898685
      BM trait4 1.147486 0.09273750
 pigauto trait4 1.124111 0.09277594
```

### OU, missing_frac = 0.30

```
  method  trait     rmse  pearson_r
      BM trait1 1.107821 0.12963658
 pigauto trait1 1.091716 0.12375281
      BM trait2 1.155393 0.15163548
 pigauto trait2 1.133356 0.15191263
      BM trait3 1.169799 0.08961219
 pigauto trait3 1.147481 0.08541191
      BM trait4 1.156355 0.10404563
 pigauto trait4 1.139342 0.10394357
```

### OU, missing_frac = 0.50

```
  method  trait     rmse  pearson_r
      BM trait1 1.138535 0.14125489
 pigauto trait1 1.119471 0.14324702
      BM trait2 1.164806 0.11984514
 pigauto trait2 1.149506 0.11947169
      BM trait3 1.164028 0.03768911
 pigauto trait3 1.149039 0.04156518
      BM trait4 1.178730 0.08583503
 pigauto trait4 1.154306 0.08034537
```

