# Binary-trait benchmark

Run on: 2026-04-10 06:36:16
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Species: 300, traits: 4, reps: 5, missing_frac: 0.25
Total wall: 10.0 min

## Methods

- **mode**: majority-class imputation (no phylogeny).
- **baseline**: phylogenetic label propagation.
- **pigauto**: full pipeline (baseline + calibrated GNN).

## Primary sweep (phylogenetic signal, threshold_quantile = 0.5)

### signal_0.2

```
   method trait  accuracy     brier
     mode  bin1 0.5106869 0.2931827
 baseline  bin1 0.5746902 0.2424626
  pigauto  bin1 0.5746902 0.2424626
     mode  bin2 0.4892830 0.2933657
 baseline  bin2 0.5801194 0.2425266
  pigauto  bin2 0.5801194 0.2425266
     mode  bin3 0.5155702 0.2763520
 baseline  bin3 0.6122610 0.2403516
  pigauto  bin3 0.6159647 0.2395552
     mode  bin4 0.4746424 0.2959516
 baseline  bin4 0.5537311 0.2456909
  pigauto  bin4 0.5537311 0.2456909
```

### signal_0.4

```
   method trait  accuracy     brier
     mode  bin1 0.4758091 0.3051731
 baseline  bin1 0.6301377 0.2315307
  pigauto  bin1 0.6268044 0.2305341
     mode  bin2 0.4943013 0.2962158
 baseline  bin2 0.6138919 0.2332056
  pigauto  bin2 0.6138919 0.2332056
     mode  bin3 0.5226393 0.2805976
 baseline  bin3 0.6080198 0.2371133
  pigauto  bin3 0.6080198 0.2371133
     mode  bin4 0.5078400 0.2589661
 baseline  bin4 0.6973949 0.2315774
  pigauto  bin4 0.6973949 0.2315774
```

### signal_0.6

```
   method trait  accuracy     brier
     mode  bin1 0.5282502 0.2656543
 baseline  bin1 0.6521893 0.2279319
  pigauto  bin1 0.6521893 0.2279319
     mode  bin2 0.4850633 0.2896968
 baseline  bin2 0.6068895 0.2345977
  pigauto  bin2 0.6097861 0.2319104
     mode  bin3 0.5075091 0.2836308
 baseline  bin3 0.6973778 0.2155216
  pigauto  bin3 0.7067528 0.2182580
     mode  bin4 0.5048352 0.2872203
 baseline  bin4 0.6492156 0.2292792
  pigauto  bin4 0.6492156 0.2292792
```

### signal_0.8

```
   method trait  accuracy     brier
     mode  bin1 0.4869412 0.2988403
 baseline  bin1 0.6629786 0.2241337
  pigauto  bin1 0.6629786 0.2241337
     mode  bin2 0.5664590 0.2655111
 baseline  bin2 0.7139051 0.2075373
  pigauto  bin2 0.7139051 0.2075373
     mode  bin3 0.5111140 0.2784045
 baseline  bin3 0.6322208 0.2278702
  pigauto  bin3 0.6322208 0.2278702
     mode  bin4 0.4780418 0.2814069
 baseline  bin4 0.6913322 0.2212722
  pigauto  bin4 0.6913322 0.2212722
```

### signal_1.0

```
   method trait  accuracy     brier
     mode  bin1 0.4899804 0.2932487
 baseline  bin1 0.6775067 0.2073715
  pigauto  bin1 0.7034326 0.2067986
     mode  bin2 0.4597599 0.3164407
 baseline  bin2 0.7201143 0.2109569
  pigauto  bin2 0.7201143 0.2109569
     mode  bin3 0.5246706 0.2759206
 baseline  bin3 0.7245681 0.2025963
  pigauto  bin3 0.7245681 0.2025963
     mode  bin4 0.5358681 0.2701922
 baseline  bin4 0.6616975 0.2177338
  pigauto  bin4 0.6616975 0.2177338
```

## Secondary sweep (class imbalance, signal = 0.6)

### imbalance_0.5

```
   method trait  accuracy     brier
     mode  bin1 0.5357446 0.2633183
 baseline  bin1 0.6586796 0.2255664
  pigauto  bin1 0.6586796 0.2231181
     mode  bin2 0.4892567 0.2880275
 baseline  bin2 0.5763727 0.2337068
  pigauto  bin2 0.5800090 0.2340663
     mode  bin3 0.4580955 0.3120753
 baseline  bin3 0.6649303 0.2194681
  pigauto  bin3 0.6649303 0.2194681
     mode  bin4 0.5351997 0.2883621
 baseline  bin4 0.6938771 0.2159868
  pigauto  bin4 0.6938771 0.2159868
```

### imbalance_0.7

```
   method trait  accuracy     brier
     mode  bin1 0.3172401 0.2500000
 baseline  bin1 0.6863963 0.1925513
  pigauto  bin1 0.6863963 0.1925513
     mode  bin2 0.2526189 0.2500000
 baseline  bin2 0.7624755 0.1738400
  pigauto  bin2 0.7624755 0.1738400
     mode  bin3 0.2598246 0.2500000
 baseline  bin3 0.7613300 0.1699391
  pigauto  bin3 0.7613300 0.1699391
     mode  bin4 0.2777733 0.2500000
 baseline  bin4 0.7601802 0.1650778
  pigauto  bin4 0.7601802 0.1650778
```

### imbalance_0.9

```
   method trait   accuracy      brier
     mode  bin1 0.08934216 0.25000000
 baseline  bin1 0.91065784 0.07348413
  pigauto  bin1 0.91065784 0.07348413
     mode  bin2 0.10829861 0.25000000
 baseline  bin2 0.89170139 0.09063631
  pigauto  bin2 0.89170139 0.09063631
     mode  bin3 0.10133759 0.25000000
 baseline  bin3 0.89866241 0.08442597
  pigauto  bin3 0.89866241 0.08442597
     mode  bin4 0.10124734 0.25000000
 baseline  bin4 0.89875266 0.09010593
  pigauto  bin4 0.89875266 0.09010593
```

