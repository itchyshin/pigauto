# Missingness-mechanism benchmark (MCAR vs MAR vs MNAR)

Run on: 2026-04-17 11:22:16
Species: 300, reps: 5
Total wall: 35.0 min

### MCAR

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4670105 0.2742888           NA
 baseline  bin1        NA        NA 0.7582888 0.1866384           NA
  pigauto  bin1        NA        NA 0.7582888 0.1866384           NA
     mean  cat1        NA        NA 0.4857650        NA           NA
 baseline  cat1        NA        NA 0.6873722        NA           NA
  pigauto  cat1        NA        NA 0.6873722        NA           NA
     mean  cnt1 0.9960974        NA        NA        NA           NA
 baseline  cnt1 0.5444791 0.8393451        NA        NA           NA
  pigauto  cnt1 0.5467374 0.8387235        NA        NA           NA
     mean cont1 0.9989989        NA        NA        NA           NA
 baseline cont1 0.4518238 0.8880068        NA        NA           NA
  pigauto cont1 0.4619473 0.8860719        NA        NA           NA
     mean cont2 1.0002461        NA        NA        NA           NA
 baseline cont2 0.4254635 0.9033542        NA        NA           NA
  pigauto cont2 0.4413946 0.8986965        NA        NA           NA
     mean  ord1 1.0135115        NA        NA        NA           NA
 baseline  ord1 0.7895720        NA        NA        NA    0.6302362
  pigauto  ord1 0.7828981        NA        NA        NA    0.6184182
```

### MAR_trait

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.5148810 0.2614484           NA
 baseline  bin1        NA        NA 0.7615858 0.1782132           NA
  pigauto  bin1        NA        NA 0.7615858 0.1782132           NA
     mean  cat1        NA        NA 0.5937631        NA           NA
 baseline  cat1        NA        NA 0.7673080        NA           NA
  pigauto  cat1        NA        NA 0.7673080        NA           NA
     mean  cnt1 0.9646366        NA        NA        NA           NA
 baseline  cnt1 0.5373236 0.8251132        NA        NA           NA
  pigauto  cnt1 0.5405253 0.8244310        NA        NA           NA
     mean cont1 1.1449939        NA        NA        NA           NA
 baseline cont1 0.5416414 0.8069876        NA        NA           NA
  pigauto cont1 0.5591713 0.8050055        NA        NA           NA
     mean cont2 1.0379911        NA        NA        NA           NA
 baseline cont2 0.5064056 0.8404554        NA        NA           NA
  pigauto cont2 0.5064056 0.8404554        NA        NA           NA
     mean  ord1 1.0029493        NA        NA        NA           NA
 baseline  ord1 0.7989974        NA        NA        NA    0.6164172
  pigauto  ord1 0.7826407        NA        NA        NA    0.6247962
```

### MAR_phylo

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4638517 0.2758331           NA
 baseline  bin1        NA        NA 0.7237387 0.1894481           NA
  pigauto  bin1        NA        NA 0.7237387 0.1894481           NA
     mean  cat1        NA        NA 0.6495571        NA           NA
 baseline  cat1        NA        NA 0.7824873        NA           NA
  pigauto  cat1        NA        NA 0.7824873        NA           NA
     mean  cnt1 1.0228798        NA        NA        NA           NA
 baseline  cnt1 0.5264138 0.8477127        NA        NA           NA
  pigauto  cnt1 0.5264138 0.8477127        NA        NA           NA
     mean cont1 1.0173087        NA        NA        NA           NA
 baseline cont1 0.4997146 0.8552217        NA        NA           NA
  pigauto cont1 0.4997146 0.8552217        NA        NA           NA
     mean cont2 1.1073704        NA        NA        NA           NA
 baseline cont2 0.5709941 0.8609475        NA        NA           NA
  pigauto cont2 0.5692082 0.8563648        NA        NA           NA
     mean  ord1 1.0076724        NA        NA        NA           NA
 baseline  ord1 0.7549669        NA        NA        NA     0.652993
  pigauto  ord1 0.7421424        NA        NA        NA     0.635041
```

### MNAR

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4790342 0.2705791           NA
 baseline  bin1        NA        NA 0.7563373 0.1807308           NA
  pigauto  bin1        NA        NA 0.7563373 0.1807308           NA
     mean  cat1        NA        NA 0.5367513        NA           NA
 baseline  cat1        NA        NA 0.7564655        NA           NA
  pigauto  cat1        NA        NA 0.7564655        NA           NA
     mean  cnt1 1.0265510        NA        NA        NA           NA
 baseline  cnt1 0.5600870 0.8474685        NA        NA           NA
  pigauto  cnt1 0.5662364 0.8461275        NA        NA           NA
     mean cont1 1.0377395        NA        NA        NA           NA
 baseline cont1 0.5489515 0.8372231        NA        NA           NA
  pigauto cont1 0.5766494 0.8324814        NA        NA           NA
     mean cont2 1.0778159        NA        NA        NA           NA
 baseline cont2 0.4993881 0.8849610        NA        NA           NA
  pigauto cont2 0.4993881 0.8849610        NA        NA           NA
     mean  ord1 1.0603380        NA        NA        NA           NA
 baseline  ord1 0.7664669        NA        NA        NA    0.6953626
  pigauto  ord1 0.7742290        NA        NA        NA    0.6964684
```

### MAR_beta_0.5

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.5325553 0.2569189           NA
 baseline  bin1        NA        NA 0.7333945 0.1800562           NA
  pigauto  bin1        NA        NA 0.7333945 0.1800562           NA
     mean  cat1        NA        NA 0.4301990        NA           NA
 baseline  cat1        NA        NA 0.6449442        NA           NA
  pigauto  cat1        NA        NA 0.6449442        NA           NA
     mean  cnt1 0.9820696        NA        NA        NA           NA
 baseline  cnt1 0.5682759 0.7989151        NA        NA           NA
  pigauto  cnt1 0.5682759 0.7989151        NA        NA           NA
     mean cont1 1.0054982        NA        NA        NA           NA
 baseline cont1 0.5239485 0.8430429        NA        NA           NA
  pigauto cont1 0.5271012 0.8415708        NA        NA           NA
     mean cont2 0.9948743        NA        NA        NA           NA
 baseline cont2 0.4731560 0.8742216        NA        NA           NA
  pigauto cont2 0.4758382 0.8746558        NA        NA           NA
     mean  ord1 0.9880507        NA        NA        NA           NA
 baseline  ord1 0.8604410        NA        NA        NA    0.5985687
  pigauto  ord1 0.8160567        NA        NA        NA    0.6013976
```

### MAR_beta_1.0

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.5323267 0.2575844           NA
 baseline  bin1        NA        NA 0.7583173 0.1822769           NA
  pigauto  bin1        NA        NA 0.7583173 0.1822769           NA
     mean  cat1        NA        NA 0.4834674        NA           NA
 baseline  cat1        NA        NA 0.7243077        NA           NA
  pigauto  cat1        NA        NA 0.7243077        NA           NA
     mean  cnt1 0.9874385        NA        NA        NA           NA
 baseline  cnt1 0.5644846 0.8199311        NA        NA           NA
  pigauto  cnt1 0.5660544 0.8199236        NA        NA           NA
     mean cont1 1.0685412        NA        NA        NA           NA
 baseline cont1 0.4792881 0.8509554        NA        NA           NA
  pigauto cont1 0.4792881 0.8509554        NA        NA           NA
     mean cont2 1.0203342        NA        NA        NA           NA
 baseline cont2 0.5597690 0.8296060        NA        NA           NA
  pigauto cont2 0.5597690 0.8296060        NA        NA           NA
     mean  ord1 0.9859620        NA        NA        NA           NA
 baseline  ord1 0.7687415        NA        NA        NA    0.6774013
  pigauto  ord1 0.7573879        NA        NA        NA    0.6653537
```

### MAR_beta_2.0

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4579036 0.2762312           NA
 baseline  bin1        NA        NA 0.7621593 0.1871643           NA
  pigauto  bin1        NA        NA 0.7621593 0.1871643           NA
     mean  cat1        NA        NA 0.4104719        NA           NA
 baseline  cat1        NA        NA 0.7024214        NA           NA
  pigauto  cat1        NA        NA 0.7024214        NA           NA
     mean  cnt1 0.9544252        NA        NA        NA           NA
 baseline  cnt1 0.6046590 0.7672103        NA        NA           NA
  pigauto  cnt1 0.6046590 0.7672103        NA        NA           NA
     mean cont1 1.0998915        NA        NA        NA           NA
 baseline cont1 0.5183713 0.8175250        NA        NA           NA
  pigauto cont1 0.5329744 0.8063057        NA        NA           NA
     mean cont2 1.0326010        NA        NA        NA           NA
 baseline cont2 0.5193352 0.8507810        NA        NA           NA
  pigauto cont2 0.5193352 0.8507810        NA        NA           NA
     mean  ord1 0.9832041        NA        NA        NA           NA
 baseline  ord1 0.7625936        NA        NA        NA    0.6746508
  pigauto  ord1 0.7580836        NA        NA        NA    0.6751926
```

### MAR_beta_4.0

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4911564 0.2686721           NA
 baseline  bin1        NA        NA 0.7712476 0.1708035           NA
  pigauto  bin1        NA        NA 0.7712476 0.1708035           NA
     mean  cat1        NA        NA 0.5977455        NA           NA
 baseline  cat1        NA        NA 0.8042719        NA           NA
  pigauto  cat1        NA        NA 0.8042719        NA           NA
     mean  cnt1 1.0463307        NA        NA        NA           NA
 baseline  cnt1 0.6421121 0.7673093        NA        NA           NA
  pigauto  cnt1 0.6421121 0.7673093        NA        NA           NA
     mean cont1 1.0505281        NA        NA        NA           NA
 baseline cont1 0.4908059 0.7674803        NA        NA           NA
  pigauto cont1 0.4908059 0.7674803        NA        NA           NA
     mean cont2 0.9813665        NA        NA        NA           NA
 baseline cont2 0.5060758 0.8449920        NA        NA           NA
  pigauto cont2 0.5060758 0.8449920        NA        NA           NA
     mean  ord1 0.9966297        NA        NA        NA           NA
 baseline  ord1 0.8977659        NA        NA        NA    0.5612406
  pigauto  ord1 0.8910799        NA        NA        NA    0.5587893
```

