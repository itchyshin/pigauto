# Count-trait benchmark

Run on: 2026-04-13 12:51:10
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 3, reps: 5, missing_frac: 0.25
Total wall: 34.5 min

## Methods

- **mean**: column-mean imputation on log1p-z scale (no phylogeny).
- **baseline**: Brownian-motion baseline on log1p-z scale (Rphylopars).
- **pigauto**: full pipeline (BM baseline + calibrated GNN).

## Primary sweep (x-axis = mean count)

### mean_5

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.6984047 0.5605624 0.7152117
  pigauto count1 0.6941426 0.5544905 0.7166518
 baseline count2 0.7379015 0.5993583 0.6870909
  pigauto count2 0.7310624 0.5961655 0.6884801
 baseline count3 0.6545848 0.5128453 0.7689963
  pigauto count3 0.6545848 0.5128453 0.7689963
```

### mean_20

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5076317 0.4141577 0.8479026
  pigauto count1 0.5076317 0.4141577 0.8479026
 baseline count2 0.5702204 0.4584912 0.8046237
  pigauto count2 0.5652111 0.4546780 0.8052225
 baseline count3 0.5946668 0.4689131 0.8033113
  pigauto count3 0.5972539 0.4710728 0.8021403
```

### mean_100

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5170611 0.4031619 0.8744981
  pigauto count1 0.5170611 0.4031619 0.8744981
 baseline count2 0.5287768 0.4169084 0.8441026
  pigauto count2 0.5309283 0.4188929 0.8438903
 baseline count3 0.5434353 0.4214173 0.8601037
  pigauto count3 0.5480791 0.4260032 0.8596304
```

### mean_500

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.5093199 0.4107953 0.8764482
  pigauto count1 0.5215935 0.4184907 0.8766837
 baseline count2 0.5105160 0.3966344 0.8797576
  pigauto count2 0.5105160 0.3966344 0.8797576
 baseline count3 0.4358218 0.3372652 0.8972600
  pigauto count3 0.4358218 0.3372652 0.8972600
```

## Secondary sweep (Poisson vs NegBin, mean_count = 20)

### poisson

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.6046351 0.4697995 0.7836874
  pigauto count1 0.5968869 0.4673437 0.7837450
 baseline count2 0.5998570 0.4621319 0.8146307
  pigauto count2 0.6035923 0.4630591 0.8131625
 baseline count3 0.5519187 0.4210210 0.8383651
  pigauto count3 0.5519187 0.4210210 0.8383651
```

### negbin

```
   method  trait      rmse       mae pearson_r
 baseline count1 0.8122090 0.6496536 0.5891892
  pigauto count1 0.7989813 0.6403532 0.5905081
 baseline count2 0.7445550 0.5844316 0.6871347
  pigauto count2 0.7394903 0.5814768 0.6862420
 baseline count3 0.7598400 0.6022959 0.6578320
  pigauto count3 0.7503484 0.5993211 0.6578217
```

