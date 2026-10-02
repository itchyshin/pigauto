# mi-posterior real-data summary (real data)

9/10 planned cells have an "ok" receipt (2 flagged non-converged; see below).

## Per-trait coverage of masked cells (model-based vs split vs Mondrian)

   dataset        arm     seed               trait n_masked n_matched
 pantheria       mcar 20260818         body_mass_g      695       695
 pantheria       mcar 20260818 head_body_length_mm      382       382
 pantheria       mcar 20260818         gestation_d      271       271
 pantheria       mcar 20260818     max_longevity_m      202       202
 pantheria       mcar 20260819         body_mass_g      695       695
 pantheria       mcar 20260819 head_body_length_mm      382       382
 pantheria       mcar 20260819         gestation_d      271       271
 pantheria       mcar 20260819     max_longevity_m      202       202
 pantheria       mcar 20260820         body_mass_g      695       695
 pantheria       mcar 20260820 head_body_length_mm      382       382
 pantheria       mcar 20260820         gestation_d      271       271
 pantheria       mcar 20260820     max_longevity_m      202       202
 pantheria structured 20260818         body_mass_g      695       695
 pantheria structured 20260818 head_body_length_mm      382       382
 pantheria structured 20260818         gestation_d      271       271
 pantheria structured 20260818     max_longevity_m      202       202
 pantheria structured 20260819         body_mass_g      695       695
 pantheria structured 20260819 head_body_length_mm      382       382
 pantheria structured 20260819         gestation_d      271       271
 pantheria structured 20260819     max_longevity_m      202       202
 pantheria structured 20260820         body_mass_g      695       695
 pantheria structured 20260820 head_body_length_mm      382       382
 pantheria structured 20260820         gestation_d      271       271
 pantheria structured 20260820     max_longevity_m      202       202
    avonet       mcar 20260818                Mass      300       300
    avonet       mcar 20260818  Beak.Length_Culmen      300       300
    avonet       mcar 20260818       Tarsus.Length      300       300
    avonet       mcar 20260818         Wing.Length      300       300
    avonet       mcar 20260819                Mass      300       300
    avonet       mcar 20260819  Beak.Length_Culmen      300       300
    avonet       mcar 20260819       Tarsus.Length      300       300
    avonet       mcar 20260819         Wing.Length      300       300
    avonet       mcar 20260820                Mass      300       300
    avonet       mcar 20260820  Beak.Length_Culmen      300       300
    avonet       mcar 20260820       Tarsus.Length      300       300
    avonet       mcar 20260820         Wing.Length      300       300
 model_coverage split_coverage mondrian_coverage
         0.9424         0.9597            0.9597
         0.9555         0.9817            0.9712
         0.9262         0.9151            0.9520
         0.9505         0.9455            0.9554
         0.9324         0.9482            0.9496
         0.9188         0.9503            0.9686
         0.9151         0.9299            0.9373
         0.9554         0.9406            0.9505
         0.9410         0.9540            0.9698
         0.9529         0.9634            0.9660
         0.9520         0.9299            0.9299
         0.9307         0.9703            0.9406
         0.9180         0.9324            0.9424
         0.9372         0.9346            0.9372
         0.9410         0.9705            0.9631
         0.9158         0.9752            0.9901
         0.8691         0.9252            0.9424
         0.9398         0.9476            0.9607
         0.9188         0.9520            0.9631
         0.9257         0.9703            0.9604
         0.8719         0.9007            0.9122
         0.9476         0.9319            0.9346
         0.9446         0.9815            0.9705
         0.9455         0.9307            0.9356
         0.9133         0.9633            0.9433
         0.9467         0.9733            0.9633
         0.9300         0.9867            0.9800
         0.9833         0.9533            0.9233
         0.9667         0.9900            0.9867
         0.9633         0.9567            0.9600
         0.9467         0.9533            0.9600
         0.9833         0.9400            0.9700
         0.9633         0.9067            0.9267
         0.9467         0.9633            0.9800
         0.9533         0.9367            0.9533
         0.9867         0.9800            1.0000

## Interval widths, like for like (mean and median of upper - lower, original scale)

   dataset        arm     seed               trait model_mean_width
 pantheria       mcar 20260818         body_mass_g           3.0260
 pantheria       mcar 20260818 head_body_length_mm           0.7539
 pantheria       mcar 20260818         gestation_d           1.0020
 pantheria       mcar 20260818     max_longevity_m           1.9490
 pantheria       mcar 20260819         body_mass_g           2.6430
 pantheria       mcar 20260819 head_body_length_mm           0.7477
 pantheria       mcar 20260819         gestation_d           0.8903
 pantheria       mcar 20260819     max_longevity_m           1.9350
 pantheria       mcar 20260820         body_mass_g           2.9670
 pantheria       mcar 20260820 head_body_length_mm           0.7674
 pantheria       mcar 20260820         gestation_d           0.9385
 pantheria       mcar 20260820     max_longevity_m           1.8520
 pantheria structured 20260818         body_mass_g           2.3700
 pantheria structured 20260818 head_body_length_mm           0.7440
 pantheria structured 20260818         gestation_d           0.9491
 pantheria structured 20260818     max_longevity_m           1.9010
 pantheria structured 20260819         body_mass_g           2.2890
 pantheria structured 20260819 head_body_length_mm           0.7575
 pantheria structured 20260819         gestation_d           0.9144
 pantheria structured 20260819     max_longevity_m           1.8330
 pantheria structured 20260820         body_mass_g           2.2120
 pantheria structured 20260820 head_body_length_mm           0.7696
 pantheria structured 20260820         gestation_d           0.9603
 pantheria structured 20260820     max_longevity_m           1.8960
    avonet       mcar 20260818                Mass         341.5000
    avonet       mcar 20260818  Beak.Length_Culmen          17.2600
    avonet       mcar 20260818       Tarsus.Length          17.8100
    avonet       mcar 20260818         Wing.Length          88.2800
    avonet       mcar 20260819                Mass         267.7000
    avonet       mcar 20260819  Beak.Length_Culmen          16.8000
    avonet       mcar 20260819       Tarsus.Length          18.4700
    avonet       mcar 20260819         Wing.Length          81.4800
    avonet       mcar 20260820                Mass         436.4000
    avonet       mcar 20260820  Beak.Length_Culmen          19.2800
    avonet       mcar 20260820       Tarsus.Length          16.2500
    avonet       mcar 20260820         Wing.Length          93.9000
 split_mean_width mondrian_mean_width model_median_width split_median_width
           4.6680              4.4900             2.3580             3.8270
           1.8560              1.6710             0.6574             1.7510
           1.0310              1.2050             0.9430             1.0130
           2.6220              2.7540             1.9200             2.7030
           4.3210              3.8650             2.2530             3.5870
           1.4390              1.4460             0.6473             1.3620
           1.1770              1.2250             0.8184             1.1500
           1.9350              1.9810             1.9550             1.9930
           4.6790              4.6680             2.4200             3.8070
           1.2500              1.2840             0.6898             1.1900
           0.8973              0.9042             0.9188             0.8813
           2.9020              2.2960             1.8480             2.9900
           3.2690              3.2490             1.8720             3.0300
           1.1980              1.4730             0.6149             1.1510
           1.3560              1.3140             0.8318             1.3440
           3.5600              3.8750             1.8610             3.6720
           3.1840              3.5410             1.7370             2.8860
           1.1650              1.2330             0.6261             1.0940
           1.5740              1.7800             0.7320             1.4870
           2.9730              2.7860             1.8250             3.0160
           2.8120              2.7720             1.7290             2.5990
           1.0280              1.0640             0.6588             0.9873
           1.5200              1.2450             0.8723             1.3970
           2.1280              2.2290             1.8540             2.1600
         733.6000            874.1000            40.8000           115.7000
          35.7000             32.2200            12.1900            28.5800
          43.4300             85.4300            12.0300            33.5700
         129.0000            129.8000            59.8500           102.8000
         611.7000            642.6000            39.5800           108.9000
          29.1400             28.9000            11.6200            23.9900
          33.8100             38.5500            11.9400            26.2300
          95.0100            132.1000            54.1800            71.1600
         531.6000            705.2000            39.8100            80.7300
          34.0600             38.7100            12.8400            26.5400
          22.1000             27.3700            12.1000            17.8600
         152.9000            245.6000            56.1600           110.2000
 mondrian_median_width
                4.1090
                1.4080
                1.1410
                2.8450
                3.4470
                1.4600
                1.2290
                2.0290
                4.0070
                1.2240
                0.8904
                2.2090
                2.7720
                1.3420
                1.1910
                3.8200
                3.2920
                1.3110
                1.9630
                2.6400
                2.5600
                1.0120
                1.2430
                2.1860
              100.7000
               23.6800
               31.7500
               95.3200
              114.8000
               22.8000
               29.7700
               93.3800
               86.8400
               27.3100
               22.1600
              152.3000

## Downstream slope check, per cell (reference PGLS vs MI-pooled, paired)

diff = MI - reference; rel_diff = diff / reference; diff_ref_se = diff / reference SE.
m_used of m_total completions entered the Rubin pool; n_nonfinite were non-finite on the
analysis scale. A row with m_used < m_total is a selective pool and fails G8 (03_acceptance.R).

   dataset        arm     seed           response           predictor    n
 pantheria       mcar 20260818        body_mass_g head_body_length_mm 1773
 pantheria       mcar 20260818        gestation_d         body_mass_g 1324
 pantheria       mcar 20260818    max_longevity_m         body_mass_g  997
 pantheria       mcar 20260819        body_mass_g head_body_length_mm 1773
 pantheria       mcar 20260819        gestation_d         body_mass_g 1324
 pantheria       mcar 20260819    max_longevity_m         body_mass_g  997
 pantheria       mcar 20260820        body_mass_g head_body_length_mm 1773
 pantheria       mcar 20260820        gestation_d         body_mass_g 1324
 pantheria       mcar 20260820    max_longevity_m         body_mass_g  997
 pantheria structured 20260818        body_mass_g head_body_length_mm 1773
 pantheria structured 20260818        gestation_d         body_mass_g 1324
 pantheria structured 20260818    max_longevity_m         body_mass_g  997
 pantheria structured 20260819        body_mass_g head_body_length_mm 1773
 pantheria structured 20260819        gestation_d         body_mass_g 1324
 pantheria structured 20260819    max_longevity_m         body_mass_g  997
 pantheria structured 20260820        body_mass_g head_body_length_mm 1773
 pantheria structured 20260820        gestation_d         body_mass_g 1324
 pantheria structured 20260820    max_longevity_m         body_mass_g  997
    avonet       mcar 20260818        Wing.Length                Mass 1500
    avonet       mcar 20260818      Tarsus.Length                Mass 1500
    avonet       mcar 20260818 Beak.Length_Culmen                Mass 1500
    avonet       mcar 20260819        Wing.Length                Mass 1500
    avonet       mcar 20260819      Tarsus.Length                Mass 1500
    avonet       mcar 20260819 Beak.Length_Culmen                Mass 1500
    avonet       mcar 20260820        Wing.Length                Mass 1500
    avonet       mcar 20260820      Tarsus.Length                Mass 1500
    avonet       mcar 20260820 Beak.Length_Culmen                Mass 1500
 ref_slope   ref_se mi_slope    mi_se m_used m_total n_nonfinite       diff
   2.85800 0.023180  2.92200 0.042380     20      20           0  6.428e-02
   0.06592 0.006423  0.05888 0.007291     20      20           0 -7.046e-03
   0.17980 0.009542  0.14840 0.011300     20      20           0 -3.138e-02
   2.85800 0.023180  2.91600 0.036260     20      20           0  5.806e-02
   0.06592 0.006423  0.05761 0.007821     20      20           0 -8.314e-03
   0.17980 0.009542  0.17440 0.011560     20      20           0 -5.401e-03
   2.85800 0.023180  2.95200 0.045600     20      20           0  9.377e-02
   0.06592 0.006423  0.05131 0.007359     20      20           0 -1.461e-02
   0.17980 0.009542  0.15500 0.012770     20      20           0 -2.483e-02
   2.85800 0.023180  2.84600 0.034980     20      20           0 -1.168e-02
   0.06592 0.006423  0.07285 0.007645     20      20           0  6.930e-03
   0.17980 0.009542  0.16730 0.011300     20      20           0 -1.244e-02
   2.85800 0.023180  2.78700 0.035350     20      20           0 -7.093e-02
   0.06592 0.006423  0.06599 0.006856     20      20           0  7.036e-05
   0.17980 0.009542  0.17560 0.010980     20      20           0 -4.153e-03
   2.85800 0.023180  2.79600 0.046640     20      20           0 -6.200e-02
   0.06592 0.006423  0.06809 0.008555     20      20           0  2.165e-03
   0.17980 0.009542  0.17610 0.012300     20      20           0 -3.685e-03
   0.33090 0.006988  0.33940 0.008964     20      20           0  8.523e-03
   0.32660 0.006610  0.33770 0.007632     20      20           0  1.110e-02
   0.35960 0.007465  0.36070 0.008570     20      20           0  1.104e-03
   0.33090 0.006988  0.32660 0.009065     20      20           0 -4.356e-03
   0.32660 0.006610  0.32950 0.007574     20      20           0  2.840e-03
   0.35960 0.007465  0.35690 0.008933     20      20           0 -2.701e-03
   0.33090 0.006988  0.32520 0.008317     20      20           0 -5.683e-03
   0.32660 0.006610  0.32810 0.007772     20      20           0  1.476e-03
   0.35960 0.007465  0.36070 0.008929     20      20           0  1.052e-03
  rel_diff diff_ref_se se_ratio within_5pct ref_status mi_status
  0.022490     2.77400    1.829        TRUE         ok        ok
 -0.106900    -1.09700    1.135       FALSE         ok        ok
 -0.174500    -3.28800    1.184       FALSE         ok        ok
  0.020320     2.50500    1.564        TRUE         ok        ok
 -0.126100    -1.29400    1.218       FALSE         ok        ok
 -0.030040    -0.56600    1.212        TRUE         ok        ok
  0.032810     4.04600    1.967        TRUE         ok        ok
 -0.221600    -2.27500    1.146       FALSE         ok        ok
 -0.138100    -2.60200    1.339       FALSE         ok        ok
 -0.004087    -0.50400    1.509        TRUE         ok        ok
  0.105100     1.07900    1.190       FALSE         ok        ok
 -0.069220    -1.30400    1.185       FALSE         ok        ok
 -0.024820    -3.06000    1.525        TRUE         ok        ok
  0.001067     0.01095    1.067        TRUE         ok        ok
 -0.023100    -0.43520    1.150        TRUE         ok        ok
 -0.021700    -2.67500    2.012        TRUE         ok        ok
  0.032830     0.33700    1.332        TRUE         ok        ok
 -0.020500    -0.38620    1.289        TRUE         ok        ok
  0.025750     1.22000    1.283        TRUE         ok        ok
  0.033990     1.68000    1.155        TRUE         ok        ok
  0.003069     0.14790    1.148        TRUE         ok        ok
 -0.013160    -0.62330    1.297        TRUE         ok        ok
  0.008695     0.42970    1.146        TRUE         ok        ok
 -0.007511    -0.36180    1.197        TRUE         ok        ok
 -0.017170    -0.81320    1.190        TRUE         ok        ok
  0.004519     0.22330    1.176        TRUE         ok        ok
  0.002925     0.14090    1.196        TRUE         ok        ok

## Downstream slope check, per pair over all cells (5% criterion reported, not gated)

Overall: 20 of 27 pair-cells with a finite relative difference are within 5%.

   dataset           response           predictor y_log x_log n_cells_attempted
 pantheria        body_mass_g head_body_length_mm FALSE FALSE                 6
 pantheria        gestation_d         body_mass_g FALSE FALSE                 6
 pantheria    max_longevity_m         body_mass_g FALSE FALSE                 6
    avonet        Wing.Length                Mass  TRUE  TRUE                 3
    avonet      Tarsus.Length                Mass  TRUE  TRUE                 3
    avonet Beak.Length_Culmen                Mass  TRUE  TRUE                 3
  fishbase             Weight              Length  TRUE  TRUE                 0
  fishbase     DepthRangeDeep              Length FALSE  TRUE                 0
  fishbase              Troph              Length FALSE  TRUE                 0
 n_cells_ok_slope n_cells_degraded has_ok_slope n_within_5pct n_rel_finite
                6                0         TRUE             6            6
                6                0         TRUE             2            6
                6                0         TRUE             3            6
                3                0         TRUE             3            3
                3                0         TRUE             3            3
                3                0         TRUE             3            3
                0                0        FALSE             0            0
                0                0        FALSE             0            0
                0                0        FALSE             0            0
 all_within_5pct max_abs_rel_diff max_abs_diff max_abs_diff_ref_se
            TRUE         0.032810     0.093770              4.0460
           FALSE         0.221600     0.014610              2.2750
           FALSE         0.174500     0.031380              3.2880
            TRUE         0.025750     0.008523              1.2200
            TRUE         0.033990     0.011100              1.6800
            TRUE         0.007511     0.002701              0.3618
              NA               NA           NA                  NA
              NA               NA           NA                  NA
              NA               NA           NA                  NA

## Convergence (max split R-hat / min bulk ESS over Sigma_P, Sigma_E, lambda)

Descriptive only; does not gate coverage or REALDATA_COMPLETE (see 03_acceptance.R).

   dataset        arm     seed                           name max_rhat min_ess
 pantheria       mcar 20260818       pantheria-mcar-m20260818    1.013   560.8
 pantheria       mcar 20260819       pantheria-mcar-m20260819    1.008   368.2
 pantheria       mcar 20260820       pantheria-mcar-m20260820    1.012   429.4
 pantheria structured 20260818 pantheria-structured-m20260818    1.013   436.4
 pantheria structured 20260819 pantheria-structured-m20260819    1.010   365.6
 pantheria structured 20260820 pantheria-structured-m20260820    1.011   581.0
    avonet       mcar 20260818          avonet-mcar-m20260818    1.003  1009.0
    avonet       mcar 20260819          avonet-mcar-m20260819    1.005  1001.0
    avonet       mcar 20260820          avonet-mcar-m20260820    1.002   854.7
 converged
      TRUE
     FALSE
      TRUE
      TRUE
     FALSE
      TRUE
      TRUE
      TRUE
      TRUE
