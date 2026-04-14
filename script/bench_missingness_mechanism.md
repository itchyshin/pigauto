# Missingness-mechanism benchmark (MCAR vs MAR vs MNAR)

Run on: 2026-04-14 05:57:51
Species: 300, reps: 5
Total wall: 55.8 min

### MCAR

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4670105 0.2742888           NA
 baseline  bin1        NA        NA 0.6318672 0.2182153           NA
  pigauto  bin1        NA        NA 0.6511734 0.2155939           NA
     mean  cat1        NA        NA 0.4857650        NA           NA
 baseline  cat1        NA        NA 0.5799578        NA           NA
  pigauto  cat1        NA        NA 0.5799578        NA           NA
     mean  cnt1 0.9960974        NA        NA        NA           NA
 baseline  cnt1 0.5326569 0.8469958        NA        NA           NA
  pigauto  cnt1 0.5326569 0.8469958        NA        NA           NA
     mean cont1 0.9989989        NA        NA        NA           NA
 baseline cont1 0.4512863 0.8886455        NA        NA           NA
  pigauto cont1 0.4512863 0.8886455        NA        NA           NA
     mean cont2 1.0002461        NA        NA        NA           NA
 baseline cont2 0.4287694 0.9033301        NA        NA           NA
  pigauto cont2 0.4287016 0.9025320        NA        NA           NA
     mean  ord1 1.0135115        NA        NA        NA           NA
 baseline  ord1 0.7405531        NA        NA        NA    0.6675156
  pigauto  ord1 0.7345484        NA        NA        NA    0.6687516
```

### MAR_trait

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.5148810 0.2614484           NA
 baseline  bin1        NA        NA 0.6606975 0.2187222           NA
  pigauto  bin1        NA        NA 0.6606975 0.2187222           NA
     mean  cat1        NA        NA 0.5937631        NA           NA
 baseline  cat1        NA        NA 0.6832604        NA           NA
  pigauto  cat1        NA        NA 0.6930165        NA           NA
     mean  cnt1 0.9646366        NA        NA        NA           NA
 baseline  cnt1 0.5382185 0.8238082        NA        NA           NA
  pigauto  cnt1 0.5378895 0.8249164        NA        NA           NA
     mean cont1 1.1449939        NA        NA        NA           NA
 baseline cont1 0.5465430 0.7997457        NA        NA           NA
  pigauto cont1 0.5465430 0.7997457        NA        NA           NA
     mean cont2 1.0379911        NA        NA        NA           NA
 baseline cont2 0.5175196 0.8329426        NA        NA           NA
  pigauto cont2 0.5175196 0.8329426        NA        NA           NA
     mean  ord1 1.0029493        NA        NA        NA           NA
 baseline  ord1 0.7507718        NA        NA        NA    0.6660258
  pigauto  ord1 0.7489197        NA        NA        NA    0.6666945
```

### MAR_phylo

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4638517 0.2758331           NA
 baseline  bin1        NA        NA 0.7290806 0.2058607           NA
  pigauto  bin1        NA        NA 0.7290806 0.2058607           NA
     mean  cat1        NA        NA 0.6495571        NA           NA
 baseline  cat1        NA        NA 0.6972334        NA           NA
  pigauto  cat1        NA        NA 0.6972334        NA           NA
     mean  cnt1 1.0228798        NA        NA        NA           NA
 baseline  cnt1 0.5189051 0.8518085        NA        NA           NA
  pigauto  cnt1 0.5189051 0.8518085        NA        NA           NA
     mean cont1 1.0173087        NA        NA        NA           NA
 baseline cont1 0.4949984 0.8579254        NA        NA           NA
  pigauto cont1 0.4949984 0.8579254        NA        NA           NA
     mean cont2 1.1073704        NA        NA        NA           NA
 baseline cont2 0.5671176 0.8621851        NA        NA           NA
  pigauto cont2 0.5759583 0.8605572        NA        NA           NA
     mean  ord1 1.0076724        NA        NA        NA           NA
 baseline  ord1 0.6991375        NA        NA        NA    0.6866553
  pigauto  ord1 0.7109139        NA        NA        NA    0.6890149
```

### MNAR

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4790342 0.2705791           NA
 baseline  bin1        NA        NA 0.7044951 0.2172401           NA
  pigauto  bin1        NA        NA 0.7044951 0.2172401           NA
     mean  cat1        NA        NA 0.5367513        NA           NA
 baseline  cat1        NA        NA 0.6380085        NA           NA
  pigauto  cat1        NA        NA 0.6380085        NA           NA
     mean  cnt1 1.0265510        NA        NA        NA           NA
 baseline  cnt1 0.5520198 0.8514336        NA        NA           NA
  pigauto  cnt1 0.5520198 0.8514336        NA        NA           NA
     mean cont1 1.0377395        NA        NA        NA           NA
 baseline cont1 0.5479025 0.8391676        NA        NA           NA
  pigauto cont1 0.5479025 0.8391676        NA        NA           NA
     mean cont2 1.0778159        NA        NA        NA           NA
 baseline cont2 0.4980096 0.8854420        NA        NA           NA
  pigauto cont2 0.4980096 0.8854420        NA        NA           NA
     mean  ord1 1.0603380        NA        NA        NA           NA
 baseline  ord1 0.6966587        NA        NA        NA    0.7473917
  pigauto  ord1 0.7072554        NA        NA        NA    0.7470029
```

### MAR_beta_0.5

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.5325553 0.2569189           NA
 baseline  bin1        NA        NA 0.7240201 0.2040398           NA
  pigauto  bin1        NA        NA 0.7240201 0.2022681           NA
     mean  cat1        NA        NA 0.4301990        NA           NA
 baseline  cat1        NA        NA 0.5522915        NA           NA
  pigauto  cat1        NA        NA 0.5522915        NA           NA
     mean  cnt1 0.9820696        NA        NA        NA           NA
 baseline  cnt1 0.5605025 0.8052623        NA        NA           NA
  pigauto  cnt1 0.5605025 0.8052623        NA        NA           NA
     mean cont1 1.0054982        NA        NA        NA           NA
 baseline cont1 0.5257813 0.8412011        NA        NA           NA
  pigauto cont1 0.5305662 0.8363595        NA        NA           NA
     mean cont2 0.9948743        NA        NA        NA           NA
 baseline cont2 0.4823222 0.8691083        NA        NA           NA
  pigauto cont2 0.4823222 0.8691083        NA        NA           NA
     mean  ord1 0.9880507        NA        NA        NA           NA
 baseline  ord1 0.7944434        NA        NA        NA    0.6096374
  pigauto  ord1 0.7839771        NA        NA        NA    0.6093227
```

### MAR_beta_1.0

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.5323267 0.2575844           NA
 baseline  bin1        NA        NA 0.6833176 0.2219383           NA
  pigauto  bin1        NA        NA 0.6833176 0.2219383           NA
     mean  cat1        NA        NA 0.4834674        NA           NA
 baseline  cat1        NA        NA 0.6077158        NA           NA
  pigauto  cat1        NA        NA 0.6077158        NA           NA
     mean  cnt1 0.9874385        NA        NA        NA           NA
 baseline  cnt1 0.5637805 0.8184992        NA        NA           NA
  pigauto  cnt1 0.5670659 0.8173945        NA        NA           NA
     mean cont1 1.0685412        NA        NA        NA           NA
 baseline cont1 0.4841246 0.8465062        NA        NA           NA
  pigauto cont1 0.4906331 0.8463730        NA        NA           NA
     mean cont2 1.0203342        NA        NA        NA           NA
 baseline cont2 0.5550262 0.8330738        NA        NA           NA
  pigauto cont2 0.5550262 0.8330738        NA        NA           NA
     mean  ord1 0.9859620        NA        NA        NA           NA
 baseline  ord1 0.7060026        NA        NA        NA    0.7081839
  pigauto  ord1 0.7074096        NA        NA        NA    0.7061178
```

### MAR_beta_2.0

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4579036 0.2762312           NA
 baseline  bin1        NA        NA 0.6830079 0.2164515           NA
  pigauto  bin1        NA        NA 0.6830079 0.2164515           NA
     mean  cat1        NA        NA 0.4104719        NA           NA
 baseline  cat1        NA        NA 0.6020283        NA           NA
  pigauto  cat1        NA        NA 0.6050134        NA           NA
     mean  cnt1 0.9544252        NA        NA        NA           NA
 baseline  cnt1 0.6012930 0.7697882        NA        NA           NA
  pigauto  cnt1 0.5977731 0.7708690        NA        NA           NA
     mean cont1 1.0998915        NA        NA        NA           NA
 baseline cont1 0.5252857 0.8156549        NA        NA           NA
  pigauto cont1 0.5268149 0.8144149        NA        NA           NA
     mean cont2 1.0326010        NA        NA        NA           NA
 baseline cont2 0.5199821 0.8480105        NA        NA           NA
  pigauto cont2 0.5199821 0.8480105        NA        NA           NA
     mean  ord1 0.9832041        NA        NA        NA           NA
 baseline  ord1 0.6874322        NA        NA        NA    0.7142993
  pigauto  ord1 0.6874322        NA        NA        NA    0.7142993
```

### MAR_beta_4.0

```
   method trait      rmse pearson_r  accuracy     brier spearman_rho
     mean  bin1        NA        NA 0.4911564 0.2686721           NA
 baseline  bin1        NA        NA 0.6718707 0.2204616           NA
  pigauto  bin1        NA        NA 0.6825850 0.2163136           NA
     mean  cat1        NA        NA 0.5977455        NA           NA
 baseline  cat1        NA        NA 0.6560232        NA           NA
  pigauto  cat1        NA        NA 0.6637155        NA           NA
     mean  cnt1 1.0463307        NA        NA        NA           NA
 baseline  cnt1 0.6333241 0.7755526        NA        NA           NA
  pigauto  cnt1 0.6333241 0.7755526        NA        NA           NA
     mean cont1 1.0505281        NA        NA        NA           NA
 baseline cont1 0.4850512 0.7688647        NA        NA           NA
  pigauto cont1 0.4850512 0.7688647        NA        NA           NA
     mean cont2 0.9813665        NA        NA        NA           NA
 baseline cont2 0.5117789 0.8403427        NA        NA           NA
  pigauto cont2 0.5172698 0.8392131        NA        NA           NA
     mean  ord1 0.9966297        NA        NA        NA           NA
 baseline  ord1 0.8428918        NA        NA        NA    0.5826163
  pigauto  ord1 0.8346950        NA        NA        NA    0.5857121
```

