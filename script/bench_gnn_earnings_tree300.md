# GNN earnings sim (tree300, n=300)

DGP: y = sqrt(0.4)*phylo + sqrt(beta)*f(cov) + sqrt(residual)*eps
Sweep: beta ∈ {0, 0.2, 0.4, 0.6}, f ∈ {linear, nonlinear, interactive}, 3 reps, 30% MCAR.

## Mean across reps (RMSE, Pearson r)

```
          method beta      f_type      rmse   pearson_r
              lm  0.0 interactive 0.9759252  0.08650951
    phylolm_blup  0.0 interactive 0.9464743  0.37464671
 pigauto_cov_sfF  0.0 interactive 0.9878795  0.30629084
 pigauto_cov_sfT  0.0 interactive 1.0164604  0.24361829
  pigauto_no_cov  0.0 interactive 0.9937223  0.27421782
              lm  0.2 interactive 0.9886237  0.04921770
    phylolm_blup  0.2 interactive 1.0074168  0.34350985
 pigauto_cov_sfF  0.2 interactive 1.0270562  0.30795157
 pigauto_cov_sfT  0.2 interactive 0.9210266  0.38934056
  pigauto_no_cov  0.2 interactive 0.9212864  0.36493840
              lm  0.4 interactive 1.0729041 -0.01933349
    phylolm_blup  0.4 interactive 1.0672111  0.25423040
 pigauto_cov_sfF  0.4 interactive 1.0689108  0.28242306
 pigauto_cov_sfT  0.4 interactive 1.0104012  0.31299138
  pigauto_no_cov  0.4 interactive 1.0492092  0.22352846
              lm  0.6 interactive 0.8923359  0.00119317
    phylolm_blup  0.6 interactive 0.9100739  0.41343980
 pigauto_cov_sfF  0.6 interactive 0.9605193  0.40754393
 pigauto_cov_sfT  0.6 interactive 0.8861807  0.44523482
  pigauto_no_cov  0.6 interactive 0.9161745  0.40392452
              lm  0.0      linear 1.0545936  0.02707092
    phylolm_blup  0.0      linear 0.9946198  0.42531215
 pigauto_cov_sfF  0.0      linear 1.0321920  0.38462592
 pigauto_cov_sfT  0.0      linear 0.9974541  0.39626121
  pigauto_no_cov  0.0      linear 1.0159255  0.37308954
              lm  0.2      linear 0.9127237  0.36648669
    phylolm_blup  0.2      linear 0.8034577  0.58996190
 pigauto_cov_sfF  0.2      linear 0.9429547  0.42804845
 pigauto_cov_sfT  0.2      linear 0.8596836  0.41609946
  pigauto_no_cov  0.2      linear 0.8851483  0.43516278
              lm  0.4      linear 0.7970004  0.63170116
    phylolm_blup  0.4      linear 0.6245580  0.79566096
 pigauto_cov_sfF  0.4      linear 0.9520590  0.43918539
 pigauto_cov_sfT  0.4      linear 0.9776971  0.38713334
  pigauto_no_cov  0.4      linear 0.9843333  0.42771777
              lm  0.6      linear 0.6332963  0.78468107
    phylolm_blup  0.6      linear 0.3724965  0.93006229
 pigauto_cov_sfF  0.6      linear 1.0227186  0.31367304
 pigauto_cov_sfT  0.6      linear 0.9676375  0.37526593
  pigauto_no_cov  0.6      linear 1.0415951  0.27087086
              lm  0.0   nonlinear 0.9551300 -0.08536557
    phylolm_blup  0.0   nonlinear 1.0156896  0.21445631
 pigauto_cov_sfF  0.0   nonlinear 0.9970190  0.27435589
 pigauto_cov_sfT  0.0   nonlinear 1.0678895  0.23383022
  pigauto_no_cov  0.0   nonlinear 1.0280504  0.28692469
              lm  0.2   nonlinear 1.0016340  0.06793781
    phylolm_blup  0.2   nonlinear 0.9840022  0.37118241
 pigauto_cov_sfF  0.2   nonlinear 1.0438721  0.29817875
 pigauto_cov_sfT  0.2   nonlinear 1.0389406  0.28830315
  pigauto_no_cov  0.2   nonlinear 0.9924113  0.34309580
              lm  0.4   nonlinear 0.9798212  0.29037825
    phylolm_blup  0.4   nonlinear 0.9633263  0.38666366
 pigauto_cov_sfF  0.4   nonlinear 1.1078210  0.19110017
 pigauto_cov_sfT  0.4   nonlinear 1.0737534  0.25568611
  pigauto_no_cov  0.4   nonlinear 1.0591157  0.25381368
              lm  0.6   nonlinear 0.9269983  0.38794020
    phylolm_blup  0.6   nonlinear 0.9407666  0.40659238
 pigauto_cov_sfF  0.6   nonlinear 0.9905414  0.28240977
 pigauto_cov_sfT  0.6   nonlinear 0.9545983  0.25716090
  pigauto_no_cov  0.6   nonlinear 1.0320425  0.16594015
```
