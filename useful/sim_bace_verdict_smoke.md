# sim_bace verdict (smoke tier)

Generated: 2026-04-26 06:27:56
Source: /Users/z3437171/Dropbox/Github Local/pigauto/script/bench_sim_bace_pigauto_smoke.rds

## Per-cell mean scores (RMSE for gaussian, accuracy for binary/threshold)

```
 response_type n_predictors n_species multi_obs phylo_signal beta_strength
      gaussian            3       100         1          0.2           0.0
      gaussian            3       100         1          0.2           0.0
      gaussian            3       100         1          0.2           0.0
      gaussian            3       100         1          0.2           0.0
      gaussian            3       100         1          0.2           0.5
      gaussian            3       100         1          0.2           0.5
      gaussian            3       100         1          0.2           0.5
      gaussian            3       100         1          0.2           0.5
      gaussian            3       100         1          0.6           0.0
      gaussian            3       100         1          0.6           0.0
      gaussian            3       100         1          0.6           0.0
      gaussian            3       100         1          0.6           0.0
      gaussian            3       100         1          0.6           0.5
      gaussian            3       100         1          0.6           0.5
      gaussian            3       100         1          0.6           0.5
      gaussian            3       100         1          0.6           0.5
      gaussian            3       100         4          0.2           0.0
      gaussian            3       100         4          0.2           0.0
      gaussian            3       100         4          0.2           0.0
      gaussian            3       100         4          0.2           0.0
      gaussian            3       100         4          0.2           0.5
      gaussian            3       100         4          0.2           0.5
      gaussian            3       100         4          0.2           0.5
      gaussian            3       100         4          0.2           0.5
      gaussian            3       100         4          0.6           0.0
      gaussian            3       100         4          0.6           0.0
      gaussian            3       100         4          0.6           0.0
      gaussian            3       100         4          0.6           0.0
      gaussian            3       100         4          0.6           0.5
      gaussian            3       100         4          0.6           0.5
      gaussian            3       100         4          0.6           0.5
      gaussian            3       100         4          0.6           0.5
              method mean_score   se_score score_metric
         column_mean   1.356902 0.14576178         rmse
                  lm   1.407760 0.19611502         rmse
 phylolm_lambda_blup   1.426515 0.20918082         rmse
     pigauto_cov_sfT   1.369745 0.13793199         rmse
         column_mean   1.884667 0.13326268         rmse
                  lm   1.809913 0.08639777         rmse
 phylolm_lambda_blup   1.772816 0.12021992         rmse
     pigauto_cov_sfT   1.917263 0.16498681         rmse
         column_mean   1.694410 0.11819431         rmse
                  lm   1.709470 0.15766987         rmse
 phylolm_lambda_blup   1.445290 0.08709554         rmse
     pigauto_cov_sfT   1.588220 0.24746088         rmse
         column_mean   2.266930 0.32143251         rmse
                  lm   2.093124 0.42813656         rmse
 phylolm_lambda_blup   1.713862 0.07649819         rmse
     pigauto_cov_sfT   2.202255 0.36274188         rmse
         column_mean   1.424931 0.03571113         rmse
                  lm   1.416345 0.03215567         rmse
 phylolm_lambda_blup   1.411333 0.03469191         rmse
     pigauto_cov_sfT   1.776951 0.07676583         rmse
         column_mean   1.588139 0.11094216         rmse
                  lm   1.446423 0.04559382         rmse
 phylolm_lambda_blup   1.481448 0.04759631         rmse
     pigauto_cov_sfT   2.069104 0.13220818         rmse
         column_mean   1.861553 0.01427752         rmse
                  lm   1.866139 0.01933336         rmse
 phylolm_lambda_blup   1.856114 0.04199258         rmse
     pigauto_cov_sfT   2.398941 0.15345285         rmse
         column_mean   1.905920 0.06680392         rmse
                  lm   1.809401 0.06732138         rmse
 phylolm_lambda_blup   1.804942 0.07726722         rmse
     pigauto_cov_sfT   2.631500 0.04357952         rmse
```

## Decisive comparison: pigauto / phylolm-lambda BLUP (gaussian only)

Below 1.0 = pigauto wins; above 1.0 = phylolm wins.

```
 response_type n_predictors n_species multi_obs phylo_signal beta_strength
      gaussian            3       100         1          0.2           0.0
      gaussian            3       100         1          0.2           0.5
      gaussian            3       100         1          0.6           0.0
      gaussian            3       100         1          0.6           0.5
      gaussian            3       100         4          0.2           0.0
      gaussian            3       100         4          0.2           0.5
      gaussian            3       100         4          0.6           0.0
      gaussian            3       100         4          0.6           0.5
 mean_score.column_mean mean_score.lm mean_score.phylolm_lambda_blup
               1.356902      1.407760                       1.426515
               1.884667      1.809913                       1.772816
               1.694410      1.709470                       1.445290
               2.266930      2.093124                       1.713862
               1.424931      1.416345                       1.411333
               1.588139      1.446423                       1.481448
               1.861553      1.866139                       1.856114
               1.905920      1.809401                       1.804942
 mean_score.pigauto_cov_sfT ratio_pig_phy ratio_pig_mean
                   1.369745     0.9602034      1.0094646
                   1.917263     1.0814788      1.0172955
                   1.588220     1.0988934      0.9373288
                   2.202255     1.2849662      0.9714705
                   1.776951     1.2590584      1.2470436
                   2.069104     1.3966766      1.3028485
                   2.398941     1.2924540      1.2886776
                   2.631500     1.4579412      1.3806983
```

### Headline summary
- pigauto / phylolm-lambda median ratio: **1.272**
- pigauto / column-mean median ratio:   **1.132**
- cells where pigauto beats phylolm:    1 / 8
- cells where pigauto beats column-mean: 2 / 8

