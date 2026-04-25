# GNN earnings sim (mammal_n1000, n=1000)

DGP: y = sqrt(0.4)*phylo + sqrt(beta)*f(cov) + sqrt(residual)*eps
Sweep: beta ∈ {0, 0.2, 0.4, 0.6}, f ∈ {linear, nonlinear, interactive}, 3 reps, 30% MCAR.

## Mean across reps (RMSE, Pearson r)

```
          method beta      f_type      rmse    pearson_r
              lm  0.0 interactive 0.9922820 -0.033475929
    phylolm_blup  0.0 interactive 0.9518406  0.429033293
 pigauto_cov_sfF  0.0 interactive 1.1672115  0.313921262
 pigauto_cov_sfT  0.0 interactive 1.0000090  0.386771516
  pigauto_no_cov  0.0 interactive 1.0072450  0.419636805
              lm  0.2 interactive 1.0328693 -0.031177488
    phylolm_blup  0.2 interactive 0.9255833  0.478414063
 pigauto_cov_sfF  0.2 interactive 1.1930644  0.346424531
 pigauto_cov_sfT  0.2 interactive 1.1441101  0.315632690
  pigauto_no_cov  0.2 interactive 1.1036103  0.346668341
              lm  0.4 interactive 1.0188363  0.039370157
    phylolm_blup  0.4 interactive 0.9746592  0.448535475
 pigauto_cov_sfF  0.4 interactive 1.1012105  0.381215553
 pigauto_cov_sfT  0.4 interactive 1.0021716  0.437485332
  pigauto_no_cov  0.4 interactive 1.0737619  0.350510944
              lm  0.6 interactive 0.9858532 -0.035792200
    phylolm_blup  0.6 interactive 0.9689626  0.423931877
 pigauto_cov_sfF  0.6 interactive 0.8203825  0.571849790
 pigauto_cov_sfT  0.6 interactive 0.7947290  0.594234321
  pigauto_no_cov  0.6 interactive 1.0711278  0.377760444
              lm  0.0      linear 0.9962081  0.019762822
    phylolm_blup  0.0      linear 0.9131628  0.501165947
 pigauto_cov_sfF  0.0      linear 1.3107488  0.299562216
 pigauto_cov_sfT  0.0      linear 1.0160568  0.380286125
  pigauto_no_cov  0.0      linear 1.0431609  0.377130539
              lm  0.2      linear 0.9501994  0.401248226
    phylolm_blup  0.2      linear 0.7405788  0.705823291
 pigauto_cov_sfF  0.2      linear 0.9874655  0.513485779
 pigauto_cov_sfT  0.2      linear 0.9337074  0.519268545
  pigauto_no_cov  0.2      linear 1.0685213  0.413374478
              lm  0.4      linear 0.7917576  0.621808624
    phylolm_blup  0.4      linear 0.5425140  0.843408178
 pigauto_cov_sfF  0.4      linear 0.8507913  0.566621221
 pigauto_cov_sfT  0.4      linear 0.8346473  0.583568611
  pigauto_no_cov  0.4      linear 1.1034660  0.406554884
              lm  0.6      linear 0.6392428  0.760462569
    phylolm_blup  0.6      linear 0.1437725  0.988645197
 pigauto_cov_sfF  0.6      linear 0.7675013  0.635728771
 pigauto_cov_sfT  0.6      linear 0.7052805  0.692036212
  pigauto_no_cov  0.6      linear 1.1714692  0.338695288
              lm  0.0   nonlinear 0.9960856 -0.002868077
    phylolm_blup  0.0   nonlinear 0.8833915  0.531086334
 pigauto_cov_sfF  0.0   nonlinear 1.0874114  0.436385349
 pigauto_cov_sfT  0.0   nonlinear 1.1791161  0.375464763
  pigauto_no_cov  0.0   nonlinear 1.0077169  0.466371725
              lm  0.2   nonlinear 0.9909862  0.167064035
    phylolm_blup  0.2   nonlinear 0.9078788  0.506054737
 pigauto_cov_sfF  0.2   nonlinear 1.1035026  0.385446010
 pigauto_cov_sfT  0.2   nonlinear 0.8827558  0.483596495
  pigauto_no_cov  0.2   nonlinear 0.9178028  0.409348259
              lm  0.4   nonlinear 0.9680426  0.239695465
    phylolm_blup  0.4   nonlinear 0.8794103  0.563019583
 pigauto_cov_sfF  0.4   nonlinear 1.0093825  0.472334475
 pigauto_cov_sfT  0.4   nonlinear 1.1667705  0.397215774
  pigauto_no_cov  0.4   nonlinear 1.2539000  0.316480169
              lm  0.6   nonlinear 0.9451869  0.285165757
    phylolm_blup  0.6   nonlinear 0.8496003  0.556834226
 pigauto_cov_sfF  0.6   nonlinear 1.2907086  0.307935195
 pigauto_cov_sfT  0.6   nonlinear 1.1076552  0.373790788
  pigauto_no_cov  0.6   nonlinear 1.0496177  0.359141387
```
