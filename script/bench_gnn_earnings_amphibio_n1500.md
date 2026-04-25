# GNN earnings sim (amphibio_n1500, n=1500)

DGP: y = sqrt(0.4)*phylo + sqrt(beta)*f(cov) + sqrt(residual)*eps
Sweep: beta ∈ {0, 0.2, 0.4, 0.6}, f ∈ {linear, nonlinear, interactive}, 3 reps, 30% MCAR.

## Mean across reps (RMSE, Pearson r)

```
          method beta      f_type       rmse  pearson_r
              lm  0.0 interactive 0.97882587 0.04525832
    phylolm_blup  0.0 interactive 0.94360727 0.46594330
 pigauto_cov_sfF  0.0 interactive 0.98310585 0.43984187
 pigauto_cov_sfT  0.0 interactive 0.93770359 0.43461549
  pigauto_no_cov  0.0 interactive 0.94632752 0.41871210
              lm  0.2 interactive 0.98926521 0.03229927
    phylolm_blup  0.2 interactive 0.95394144 0.45564864
 pigauto_cov_sfF  0.2 interactive 1.08855217 0.40197922
 pigauto_cov_sfT  0.2 interactive 0.92571271 0.42411007
  pigauto_no_cov  0.2 interactive 0.91586757 0.43280121
              lm  0.4 interactive 0.98604060 0.01733384
    phylolm_blup  0.4 interactive 0.94379599 0.47245068
 pigauto_cov_sfF  0.4 interactive 1.00460435 0.43193898
 pigauto_cov_sfT  0.4 interactive 0.97525894 0.42860607
  pigauto_no_cov  0.4 interactive 0.97254596 0.40270344
              lm  0.6 interactive 0.97810708 0.02010989
    phylolm_blup  0.6 interactive 0.94375839 0.47660907
 pigauto_cov_sfF  0.6 interactive 0.86194826 0.54756238
 pigauto_cov_sfT  0.6 interactive 0.81861608 0.56346886
  pigauto_no_cov  0.6 interactive 1.01701836 0.39785928
              lm  0.0      linear 0.99884787 0.01866519
    phylolm_blup  0.0      linear 0.93473561 0.48812897
 pigauto_cov_sfF  0.0      linear 1.03670871 0.43264165
 pigauto_cov_sfT  0.0      linear 0.89039973 0.45382311
  pigauto_no_cov  0.0      linear 0.90476853 0.43912434
              lm  0.2      linear 0.87046153 0.47296149
    phylolm_blup  0.2      linear 0.75160263 0.69019559
 pigauto_cov_sfF  0.2      linear 0.96411771 0.47864305
 pigauto_cov_sfT  0.2      linear 0.87092762 0.52710817
  pigauto_no_cov  0.2      linear 0.90362467 0.47804453
              lm  0.4      linear 0.77897673 0.61954782
    phylolm_blup  0.4      linear 0.54897118 0.83930068
 pigauto_cov_sfF  0.4      linear 0.89757914 0.50338087
 pigauto_cov_sfT  0.4      linear 0.84999534 0.52666152
  pigauto_no_cov  0.4      linear 0.97055830 0.37129495
              lm  0.6      linear 0.63016412 0.77958900
    phylolm_blup  0.6      linear 0.04801041 0.99885734
 pigauto_cov_sfF  0.6      linear 0.75548268 0.66238356
 pigauto_cov_sfT  0.6      linear 0.75323890 0.69765296
  pigauto_no_cov  0.6      linear 0.97404649 0.43247184
              lm  0.0   nonlinear 0.98425783 0.02989994
    phylolm_blup  0.0   nonlinear 0.96630733 0.43723389
 pigauto_cov_sfF  0.0   nonlinear 1.03581863 0.40980143
 pigauto_cov_sfT  0.0   nonlinear 0.97616669 0.40247616
  pigauto_no_cov  0.0   nonlinear 0.93170245 0.44592052
              lm  0.2   nonlinear 0.99235501 0.14685574
    phylolm_blup  0.2   nonlinear 0.95984398 0.46335488
 pigauto_cov_sfF  0.2   nonlinear 1.03407306 0.41008288
 pigauto_cov_sfT  0.2   nonlinear 0.91224212 0.44351460
  pigauto_no_cov  0.2   nonlinear 0.93132809 0.42103660
              lm  0.4   nonlinear 0.96870422 0.20255574
    phylolm_blup  0.4   nonlinear 0.89836172 0.53040091
 pigauto_cov_sfF  0.4   nonlinear 1.07716172 0.39076869
 pigauto_cov_sfT  0.4   nonlinear 0.88339531 0.47100504
  pigauto_no_cov  0.4   nonlinear 0.90178522 0.43622591
              lm  0.6   nonlinear 0.95138184 0.29653154
    phylolm_blup  0.6   nonlinear 0.87617231 0.56372938
 pigauto_cov_sfF  0.6   nonlinear 0.96810991 0.46539316
 pigauto_cov_sfT  0.6   nonlinear 0.86893934 0.50060707
  pigauto_no_cov  0.6   nonlinear 0.91375064 0.44304372
```
