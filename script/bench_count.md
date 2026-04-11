# Count-trait benchmark

Run on: 2026-04-10 06:54:34
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 3, reps: 5, missing_frac: 0.25
Total wall: 6.8 min

## Methods

- **mean**: column-mean imputation on log1p-z scale (no phylogeny).
- **baseline**: Brownian-motion baseline on log1p-z scale (Rphylopars).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (x-axis = mean count)

### mean_5

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.7000593 0.5620885 0.7131791
  pigauto count1 0.6913485 0.5525000 0.7146812
 baseline count2 0.7405025 0.5939283 0.6853644
  pigauto count2 0.7405025 0.5939283 0.6853644
 baseline count3 0.6561535 0.5183689 0.7694328
  pigauto count3 0.6561535 0.5183689 0.7694328
```

### mean_20

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5091533 0.4146581 0.8460041
  pigauto count1 0.5091533 0.4146581 0.8460041
 baseline count2 0.5623869 0.4517998 0.8109011
  pigauto count2 0.5590139 0.4488119 0.8120454
 baseline count3 0.6046432 0.4751993 0.7961838
  pigauto count3 0.6072417 0.4779426 0.7948055
```

### mean_100

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5100635 0.3991978 0.8779281
  pigauto count1 0.5100635 0.3991978 0.8779281
 baseline count2 0.5329371 0.4185469 0.8398738
  pigauto count2 0.5364345 0.4218608 0.8394172
 baseline count3 0.5407721 0.4223225 0.8609755
  pigauto count3 0.5456272 0.4284844 0.8601965
```

### mean_500

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5111837 0.4084474 0.8763128
  pigauto count1 0.5223400 0.4164658 0.8766758
 baseline count2 0.5078690 0.3939186 0.8806644
  pigauto count2 0.5078690 0.3939186 0.8806644
 baseline count3 0.4316448 0.3390723 0.8998590
  pigauto count3 0.4316448 0.3390723 0.8998590
```

## Secondary sweep (Poisson vs NegBin, mean_count = 20)

### poisson

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.6076440 0.4718866 0.7811223
  pigauto count1 0.6076440 0.4718866 0.7811223
 baseline count2 0.5978294 0.4604726 0.8161929
  pigauto count2 0.6031275 0.4631175 0.8143913
 baseline count3 0.5503606 0.4199198 0.8404936
  pigauto count3 0.5503606 0.4199198 0.8404936
```

### negbin

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.8182266 0.6561378 0.5886514
  pigauto count1 0.8012019 0.6437680 0.5902301
 baseline count2 0.7500370 0.5874533 0.6893585
  pigauto count2 0.7441285 0.5854482 0.6880876
 baseline count3 0.7626462 0.6076534 0.6564664
  pigauto count3 0.7485059 0.6011484 0.6575390
```

