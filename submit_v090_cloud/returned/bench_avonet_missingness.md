# AVONET full-scale missingness sweep

Run on: 2026-04-18 15:41:35
Machine: Linux 5.14.0-611.42.1.el9_7.x86_64 (x86_64), R 4.4.0
Species: 9993, traits: 7
Total wall: 273.9 min

## Methods

- **mean**: column mean / mode imputation (no phylogeny).
- **BM**: Brownian-motion baseline (Rphylopars + phylogenetic label
  propagation for discrete traits).
- **pigauto**: full pipeline (BM baseline + calibrated GNN delta).

## Test-set metrics (per trait)

### missing_frac = 0.20

```
  method              trait        type    n      rmse pearson_r spearman_rho
      BM Beak.Length_Culmen  continuous 1041 0.2963931 0.9565484           NA
    mean Beak.Length_Culmen  continuous 1041 1.0157984        NA           NA
 pigauto Beak.Length_Culmen  continuous 1041 0.2963931 0.9565484           NA
      BM               Mass  continuous  999 0.2677125 0.9636107           NA
    mean               Mass  continuous  999 1.0009498        NA           NA
 pigauto               Mass  continuous  999 0.4056547 0.9646187           NA
      BM          Migration     ordinal  918 0.7747608        NA    0.5639290
    mean          Migration     ordinal  918 0.9974487        NA           NA
 pigauto          Migration     ordinal  918 0.7747608        NA    0.5639291
      BM  Primary.Lifestyle categorical  977        NA        NA           NA
    mean  Primary.Lifestyle categorical  977        NA        NA           NA
 pigauto  Primary.Lifestyle categorical  977        NA        NA           NA
      BM      Tarsus.Length  continuous 1011 0.2759708 0.9594579           NA
    mean      Tarsus.Length  continuous 1011 0.9786237        NA           NA
 pigauto      Tarsus.Length  continuous 1011 0.2759708 0.9594579           NA
      BM      Trophic.Level categorical 1038        NA        NA           NA
    mean      Trophic.Level categorical 1038        NA        NA           NA
 pigauto      Trophic.Level categorical 1038        NA        NA           NA
      BM        Wing.Length  continuous 1009 0.2724725 0.9698869           NA
    mean        Wing.Length  continuous 1009 1.1017743        NA           NA
 pigauto        Wing.Length  continuous 1009 0.3051831 0.9720296           NA
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
 0.6120778
        NA
        NA
        NA
 0.5722543
 0.5558767
 0.5173410
        NA
        NA
        NA
```

### missing_frac = 0.50

```
  method              trait        type    n      rmse pearson_r spearman_rho
      BM Beak.Length_Culmen  continuous 2518 0.3538667 0.9372087           NA
    mean Beak.Length_Culmen  continuous 2518 1.0121054        NA           NA
 pigauto Beak.Length_Culmen  continuous 2518 0.3538667 0.9372087           NA
      BM               Mass  continuous 2498 0.2821525 0.9605139           NA
    mean               Mass  continuous 2498 1.0134188        NA           NA
 pigauto               Mass  continuous 2498 0.2821525 0.9605139           NA
      BM          Migration     ordinal 2470 0.7848863        NA    0.5476193
    mean          Migration     ordinal 2470 0.9849580        NA           NA
 pigauto          Migration     ordinal 2470 0.7848863        NA    0.5476187
      BM  Primary.Lifestyle categorical 2459        NA        NA           NA
    mean  Primary.Lifestyle categorical 2459        NA        NA           NA
 pigauto  Primary.Lifestyle categorical 2459        NA        NA           NA
      BM      Tarsus.Length  continuous 2462 0.2893810 0.9572937           NA
    mean      Tarsus.Length  continuous 2462 1.0003332        NA           NA
 pigauto      Tarsus.Length  continuous 2462 0.2893810 0.9572937           NA
      BM      Trophic.Level categorical 2522        NA        NA           NA
    mean      Trophic.Level categorical 2522        NA        NA           NA
 pigauto      Trophic.Level categorical 2522        NA        NA           NA
      BM        Wing.Length  continuous 2553 0.3716355 0.9239390           NA
    mean        Wing.Length  continuous 2553 0.9697099        NA           NA
 pigauto        Wing.Length  continuous 2553 0.3716355 0.9239390           NA
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
 0.6217975
 0.5884506
 0.5811305
        NA
        NA
        NA
 0.5927835
 0.5626487
 0.6046788
        NA
        NA
        NA
```

### missing_frac = 0.80

```
  method              trait        type    n      rmse pearson_r spearman_rho
      BM Beak.Length_Culmen  continuous 3900 0.4272373 0.9029234           NA
    mean Beak.Length_Culmen  continuous 3900 0.9930025        NA           NA
 pigauto Beak.Length_Culmen  continuous 3900 0.4272373 0.9029234           NA
      BM               Mass  continuous 3980 0.3231053 0.9479004           NA
    mean               Mass  continuous 3980 1.0127319        NA           NA
 pigauto               Mass  continuous 3980 0.3231053 0.9479004           NA
      BM          Migration     ordinal 4041 0.8505674        NA    0.4960723
    mean          Migration     ordinal 4041 1.0168020        NA           NA
 pigauto          Migration     ordinal 4041 0.8505674        NA    0.4960762
      BM  Primary.Lifestyle categorical 3986        NA        NA           NA
    mean  Primary.Lifestyle categorical 3986        NA        NA           NA
 pigauto  Primary.Lifestyle categorical 3986        NA        NA           NA
      BM      Tarsus.Length  continuous 4019 0.3615476 0.9331197           NA
    mean      Tarsus.Length  continuous 4019 1.0035918        NA           NA
 pigauto      Tarsus.Length  continuous 4019 0.3615476 0.9331197           NA
      BM      Trophic.Level categorical 4003        NA        NA           NA
    mean      Trophic.Level categorical 4003        NA        NA           NA
 pigauto      Trophic.Level categorical 4003        NA        NA           NA
      BM        Wing.Length  continuous 4042 0.3185733 0.9454351           NA
    mean        Wing.Length  continuous 4042 0.9778393        NA           NA
 pigauto        Wing.Length  continuous 4042 0.3185733 0.9454351           NA
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
 0.6194180
 0.5893126
 0.6204215
        NA
        NA
        NA
 0.5910567
 0.5548339
 0.5910567
        NA
        NA
        NA
```

