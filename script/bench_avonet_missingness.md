# AVONET full-scale missingness sweep

Run on: 2026-04-09 15:32:02
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 9993, traits: 7
Total wall: 59.4 min

## Methods

- **mean**: column mean / mode imputation (no phylogeny).
- **BM**: Brownian-motion baseline (Rphylopars + phylogenetic label
  propagation for discrete traits).
- **pigauto**: full pipeline (BM baseline + calibrated GNN delta).

## Test-set metrics (per trait)

### missing_frac = 0.20

```
  method              trait        type    n      rmse pearson_r spearman_rho
      BM Beak.Length_Culmen  continuous 1041 0.2072963 0.9791111           NA
    mean Beak.Length_Culmen  continuous 1041 1.0157984        NA           NA
 pigauto Beak.Length_Culmen  continuous 1041 0.2072963 0.9791111           NA
      BM               Mass  continuous  999 0.1470619 0.9891742           NA
    mean               Mass  continuous  999 1.0009498        NA           NA
 pigauto               Mass  continuous  999 0.1470619 0.9891742           NA
      BM          Migration     ordinal  918 0.7660908        NA    0.5508744
    mean          Migration     ordinal  918 0.9974487        NA           NA
 pigauto          Migration     ordinal  918 0.7660908        NA    0.5508744
      BM  Primary.Lifestyle categorical  977        NA        NA           NA
    mean  Primary.Lifestyle categorical  977        NA        NA           NA
 pigauto  Primary.Lifestyle categorical  977        NA        NA           NA
      BM      Tarsus.Length  continuous 1011 0.1972703 0.9794662           NA
    mean      Tarsus.Length  continuous 1011 0.9786237        NA           NA
 pigauto      Tarsus.Length  continuous 1011 0.1972703 0.9794662           NA
      BM      Trophic.Level categorical 1038        NA        NA           NA
    mean      Trophic.Level categorical 1038        NA        NA           NA
 pigauto      Trophic.Level categorical 1038        NA        NA           NA
      BM        Wing.Length  continuous 1009 0.1558959 0.9905010           NA
    mean        Wing.Length  continuous 1009 1.1017743        NA           NA
 pigauto        Wing.Length  continuous 1009 0.1558959 0.9905010           NA
  accuracy
        NA
        NA
        NA
        NA
        NA
        NA
        NA
        NA
        NA
 0.5957011
 0.5762538
 0.4298874
        NA
        NA
        NA
 0.5722543
 0.5558767
 0.5915222
        NA
        NA
        NA
```

### missing_frac = 0.50

```
  method              trait        type    n      rmse pearson_r spearman_rho
      BM Beak.Length_Culmen  continuous 2518 0.2637200 0.9654796           NA
    mean Beak.Length_Culmen  continuous 2518 1.0121054        NA           NA
 pigauto Beak.Length_Culmen  continuous 2518 0.2637200 0.9654796           NA
      BM               Mass  continuous 2498 0.1796257 0.9841747           NA
    mean               Mass  continuous 2498 1.0134188        NA           NA
 pigauto               Mass  continuous 2498 0.1796257 0.9841747           NA
      BM          Migration     ordinal 2470 0.7825293        NA    0.5268293
    mean          Migration     ordinal 2470 0.9849580        NA           NA
 pigauto          Migration     ordinal 2470 0.7825293        NA    0.5268293
      BM  Primary.Lifestyle categorical 2459        NA        NA           NA
    mean  Primary.Lifestyle categorical 2459        NA        NA           NA
 pigauto  Primary.Lifestyle categorical 2459        NA        NA           NA
      BM      Tarsus.Length  continuous 2462 0.2123486 0.9772153           NA
    mean      Tarsus.Length  continuous 2462 1.0003332        NA           NA
 pigauto      Tarsus.Length  continuous 2462 0.2123486 0.9772153           NA
      BM      Trophic.Level categorical 2522        NA        NA           NA
    mean      Trophic.Level categorical 2522        NA        NA           NA
 pigauto      Trophic.Level categorical 2522        NA        NA           NA
      BM        Wing.Length  continuous 2553 0.3001447 0.9509328           NA
    mean        Wing.Length  continuous 2553 0.9697099        NA           NA
 pigauto        Wing.Length  continuous 2553 0.3001447 0.9509328           NA
  accuracy
        NA
        NA
        NA
        NA
        NA
        NA
        NA
        NA
        NA
 0.6059374
 0.5884506
 0.6498577
        NA
        NA
        NA
 0.5800952
 0.5626487
 0.6974623
        NA
        NA
        NA
```

### missing_frac = 0.80

```
  method              trait        type    n      rmse pearson_r spearman_rho
      BM Beak.Length_Culmen  continuous 3900 0.3436194 0.9382612           NA
    mean Beak.Length_Culmen  continuous 3900 0.9930025        NA           NA
 pigauto Beak.Length_Culmen  continuous 3900 0.3436194 0.9382612           NA
      BM               Mass  continuous 3980 0.2495982 0.9691992           NA
    mean               Mass  continuous 3980 1.0127319        NA           NA
 pigauto               Mass  continuous 3980 0.2495982 0.9691992           NA
      BM          Migration     ordinal 4041 0.8450629        NA    0.4953706
    mean          Migration     ordinal 4041 1.0168020        NA           NA
 pigauto          Migration     ordinal 4041 0.8450629        NA    0.4953710
      BM  Primary.Lifestyle categorical 3986        NA        NA           NA
    mean  Primary.Lifestyle categorical 3986        NA        NA           NA
 pigauto  Primary.Lifestyle categorical 3986        NA        NA           NA
      BM      Tarsus.Length  continuous 4019 0.2904315 0.9572934           NA
    mean      Tarsus.Length  continuous 4019 1.0035918        NA           NA
 pigauto      Tarsus.Length  continuous 4019 0.2904315 0.9572934           NA
      BM      Trophic.Level categorical 4003        NA        NA           NA
    mean      Trophic.Level categorical 4003        NA        NA           NA
 pigauto      Trophic.Level categorical 4003        NA        NA           NA
      BM        Wing.Length  continuous 4042 0.2451478 0.9680665           NA
    mean        Wing.Length  continuous 4042 0.9778393        NA           NA
 pigauto        Wing.Length  continuous 4042 0.2451478 0.9680665           NA
  accuracy
        NA
        NA
        NA
        NA
        NA
        NA
        NA
        NA
        NA
 0.6043653
 0.5893126
 0.6555444
        NA
        NA
        NA
 0.5720709
 0.5548339
 0.5715713
        NA
        NA
        NA
```

