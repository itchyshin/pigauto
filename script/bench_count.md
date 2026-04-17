# Count-trait benchmark

Run on: 2026-04-17 07:59:49
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 3, reps: 5, missing_frac: 0.25
Total wall: 21.4 min

## Methods

- **mean**: column-mean imputation on log1p-z scale (no phylogeny).
- **baseline**: Brownian-motion baseline on log1p-z scale (Rphylopars).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (x-axis = mean count)

### mean_5

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.7000593 0.5620885 0.7131791
  pigauto count1 0.6934007 0.5549006 0.7146712
 baseline count2 0.7405025 0.5939283 0.6853644
  pigauto count2 0.7361405 0.5937839 0.6848626
 baseline count3 0.6561535 0.5183689 0.7694328
  pigauto count3 0.6561535 0.5183689 0.7694328
```

### mean_20

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5091533 0.4146581 0.8460041
  pigauto count1 0.5091533 0.4146581 0.8460041
 baseline count2 0.5623869 0.4517998 0.8109011
  pigauto count2 0.5553936 0.4474525 0.8120074
 baseline count3 0.6046432 0.4751993 0.7961838
  pigauto count3 0.6059701 0.4760237 0.7959815
```

### mean_100

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5100635 0.3991978 0.8779281
  pigauto count1 0.5100635 0.3991978 0.8779281
 baseline count2 0.5329371 0.4185469 0.8398738
  pigauto count2 0.5386653 0.4246798 0.8373602
 baseline count3 0.5407721 0.4223225 0.8609755
  pigauto count3 0.5399331 0.4230222 0.8611074
```

### mean_500

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5111837 0.4084474 0.8763128
  pigauto count1 0.5295014 0.4227677 0.8725269
 baseline count2 0.5078690 0.3939186 0.8806644
  pigauto count2 0.5160178 0.4034294 0.8795183
 baseline count3 0.4316448 0.3390723 0.8998590
  pigauto count3 0.4316448 0.3390723 0.8998590
```

## Secondary sweep (Poisson vs NegBin, mean_count = 20)

### poisson

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.6076440 0.4718866 0.7811223
  pigauto count1 0.6033339 0.4705716 0.7800903
 baseline count2 0.5978294 0.4604726 0.8161929
  pigauto count2 0.6015694 0.4619106 0.8143331
 baseline count3 0.5503606 0.4199198 0.8404936
  pigauto count3 0.5503606 0.4199198 0.8404936
```

### negbin

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.8182266 0.6561378 0.5886514
  pigauto count1 0.8059891 0.6483029 0.5892285
 baseline count2 0.7500370 0.5874533 0.6893585
  pigauto count2 0.7442200 0.5859880 0.6883807
 baseline count3 0.7626462 0.6076534 0.6564664
  pigauto count3 0.7563455 0.6032514 0.6552011
```

