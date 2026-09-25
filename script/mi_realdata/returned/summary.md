# mi-posterior real-data summary (real data)

10/10 planned cells have an "ok" receipt (2 flagged non-converged; see below).

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
  fishbase structured 20260818              Length     2047      2047
  fishbase structured 20260818              Weight      377       377
  fishbase structured 20260818      DepthRangeDeep      938       938
  fishbase structured 20260818       Vulnerability     2097      2097
  fishbase structured 20260818               Troph      966       966
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
         0.9419         0.9394            0.9399
         0.9390         0.9496            0.9788
         0.9574         0.9659            0.9616
         0.9499         0.9542            0.9495
         0.9327         0.9389            0.9462

## Interval widths, like for like (mean and median of upper - lower, original scale)

   dataset        arm     seed               trait model_mean_width
 pantheria       mcar 20260818         body_mass_g        3.026e+00
 pantheria       mcar 20260818 head_body_length_mm        7.539e-01
 pantheria       mcar 20260818         gestation_d        1.002e+00
 pantheria       mcar 20260818     max_longevity_m        1.949e+00
 pantheria       mcar 20260819         body_mass_g        2.643e+00
 pantheria       mcar 20260819 head_body_length_mm        7.477e-01
 pantheria       mcar 20260819         gestation_d        8.903e-01
 pantheria       mcar 20260819     max_longevity_m        1.935e+00
 pantheria       mcar 20260820         body_mass_g        2.967e+00
 pantheria       mcar 20260820 head_body_length_mm        7.674e-01
 pantheria       mcar 20260820         gestation_d        9.385e-01
 pantheria       mcar 20260820     max_longevity_m        1.852e+00
 pantheria structured 20260818         body_mass_g        2.370e+00
 pantheria structured 20260818 head_body_length_mm        7.440e-01
 pantheria structured 20260818         gestation_d        9.491e-01
 pantheria structured 20260818     max_longevity_m        1.901e+00
 pantheria structured 20260819         body_mass_g        2.289e+00
 pantheria structured 20260819 head_body_length_mm        7.575e-01
 pantheria structured 20260819         gestation_d        9.144e-01
 pantheria structured 20260819     max_longevity_m        1.833e+00
 pantheria structured 20260820         body_mass_g        2.212e+00
 pantheria structured 20260820 head_body_length_mm        7.696e-01
 pantheria structured 20260820         gestation_d        9.603e-01
 pantheria structured 20260820     max_longevity_m        1.896e+00
    avonet       mcar 20260818                Mass        3.415e+02
    avonet       mcar 20260818  Beak.Length_Culmen        1.726e+01
    avonet       mcar 20260818       Tarsus.Length        1.781e+01
    avonet       mcar 20260818         Wing.Length        8.828e+01
    avonet       mcar 20260819                Mass        2.677e+02
    avonet       mcar 20260819  Beak.Length_Culmen        1.680e+01
    avonet       mcar 20260819       Tarsus.Length        1.847e+01
    avonet       mcar 20260819         Wing.Length        8.148e+01
    avonet       mcar 20260820                Mass        4.364e+02
    avonet       mcar 20260820  Beak.Length_Culmen        1.928e+01
    avonet       mcar 20260820       Tarsus.Length        1.625e+01
    avonet       mcar 20260820         Wing.Length        9.390e+01
  fishbase structured 20260818              Length        4.019e+01
  fishbase structured 20260818              Weight        1.368e+05
  fishbase structured 20260818      DepthRangeDeep        1.897e+03
  fishbase structured 20260818       Vulnerability        2.718e+01
  fishbase structured 20260818               Troph        1.746e+00
 split_mean_width mondrian_mean_width model_median_width split_median_width
        4.668e+00           4.490e+00             2.3580          3.827e+00
        1.856e+00           1.671e+00             0.6574          1.751e+00
        1.031e+00           1.205e+00             0.9430          1.013e+00
        2.622e+00           2.754e+00             1.9200          2.703e+00
        4.321e+00           3.865e+00             2.2530          3.587e+00
        1.439e+00           1.446e+00             0.6473          1.362e+00
        1.177e+00           1.225e+00             0.8184          1.150e+00
        1.935e+00           1.981e+00             1.9550          1.993e+00
        4.679e+00           4.668e+00             2.4200          3.807e+00
        1.250e+00           1.284e+00             0.6898          1.190e+00
        8.973e-01           9.042e-01             0.9188          8.813e-01
        2.902e+00           2.296e+00             1.8480          2.990e+00
        3.269e+00           3.249e+00             1.8720          3.030e+00
        1.198e+00           1.473e+00             0.6149          1.151e+00
        1.356e+00           1.314e+00             0.8318          1.344e+00
        3.560e+00           3.875e+00             1.8610          3.672e+00
        3.184e+00           3.541e+00             1.7370          2.886e+00
        1.165e+00           1.233e+00             0.6261          1.094e+00
        1.574e+00           1.780e+00             0.7320          1.487e+00
        2.973e+00           2.786e+00             1.8250          3.016e+00
        2.812e+00           2.772e+00             1.7290          2.599e+00
        1.028e+00           1.064e+00             0.6588          9.873e-01
        1.520e+00           1.245e+00             0.8723          1.397e+00
        2.128e+00           2.229e+00             1.8540          2.160e+00
        7.336e+02           8.741e+02            40.8000          1.157e+02
        3.570e+01           3.222e+01            12.1900          2.858e+01
        4.343e+01           8.543e+01            12.0300          3.357e+01
        1.290e+02           1.298e+02            59.8500          1.028e+02
        6.117e+02           6.426e+02            39.5800          1.089e+02
        2.914e+01           2.890e+01            11.6200          2.399e+01
        3.381e+01           3.855e+01            11.9400          2.623e+01
        9.501e+01           1.321e+02            54.1800          7.116e+01
        5.316e+02           7.052e+02            39.8100          8.073e+01
        3.406e+01           3.871e+01            12.8400          2.654e+01
        2.210e+01           2.737e+01            12.1000          1.786e+01
        1.529e+02           2.456e+02            56.1600          1.102e+02
        6.636e+01           6.526e+01            25.1200          4.206e+01
        3.400e+05           3.243e+05          6722.0000          6.422e+04
        2.237e+03           2.365e+03          1822.0000          2.237e+03
        4.645e+01           4.300e+01            23.6900          4.244e+01
        1.807e+00           1.877e+00             1.7420          1.833e+00
 mondrian_median_width
             4.109e+00
             1.408e+00
             1.141e+00
             2.845e+00
             3.447e+00
             1.460e+00
             1.229e+00
             2.029e+00
             4.007e+00
             1.224e+00
             8.904e-01
             2.209e+00
             2.772e+00
             1.342e+00
             1.191e+00
             3.820e+00
             3.292e+00
             1.311e+00
             1.963e+00
             2.640e+00
             2.560e+00
             1.012e+00
             1.243e+00
             2.186e+00
             1.007e+02
             2.368e+01
             3.175e+01
             9.532e+01
             1.148e+02
             2.280e+01
             2.977e+01
             9.338e+01
             8.684e+01
             2.731e+01
             2.216e+01
             1.523e+02
             4.187e+01
             6.504e+04
             3.180e+03
             3.962e+01
             1.908e+00

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
  fishbase structured 20260818             Weight              Length 1879
  fishbase structured 20260818     DepthRangeDeep              Length 4638
  fishbase structured 20260818              Troph              Length 4818
 ref_slope    ref_se  mi_slope     mi_se m_used m_total n_nonfinite       diff
   2.85800  0.023180   2.92200  0.042380     20      20           0  6.428e-02
   0.06592  0.006423   0.05888  0.007291     20      20           0 -7.046e-03
   0.17980  0.009542   0.14840  0.011300     20      20           0 -3.138e-02
   2.85800  0.023180   2.91600  0.036260     20      20           0  5.806e-02
   0.06592  0.006423   0.05761  0.007821     20      20           0 -8.314e-03
   0.17980  0.009542   0.17440  0.011560     20      20           0 -5.401e-03
   2.85800  0.023180   2.95200  0.045600     20      20           0  9.377e-02
   0.06592  0.006423   0.05131  0.007359     20      20           0 -1.461e-02
   0.17980  0.009542   0.15500  0.012770     20      20           0 -2.483e-02
   2.85800  0.023180   2.84600  0.034980     20      20           0 -1.168e-02
   0.06592  0.006423   0.07285  0.007645     20      20           0  6.930e-03
   0.17980  0.009542   0.16730  0.011300     20      20           0 -1.244e-02
   2.85800  0.023180   2.78700  0.035350     20      20           0 -7.093e-02
   0.06592  0.006423   0.06599  0.006856     20      20           0  7.036e-05
   0.17980  0.009542   0.17560  0.010980     20      20           0 -4.153e-03
   2.85800  0.023180   2.79600  0.046640     20      20           0 -6.200e-02
   0.06592  0.006423   0.06809  0.008555     20      20           0  2.165e-03
   0.17980  0.009542   0.17610  0.012300     20      20           0 -3.685e-03
   0.33090  0.006988   0.33940  0.008964     20      20           0  8.523e-03
   0.32660  0.006610   0.33770  0.007632     20      20           0  1.110e-02
   0.35960  0.007465   0.36070  0.008570     20      20           0  1.104e-03
   0.33090  0.006988   0.32660  0.009065     20      20           0 -4.356e-03
   0.32660  0.006610   0.32950  0.007574     20      20           0  2.840e-03
   0.35960  0.007465   0.35690  0.008933     20      20           0 -2.701e-03
   0.33090  0.006988   0.32520  0.008317     20      20           0 -5.683e-03
   0.32660  0.006610   0.32810  0.007772     20      20           0  1.476e-03
   0.35960  0.007465   0.36070  0.008929     20      20           0  1.052e-03
   2.75600  0.031060   2.73500  0.036070     20      20           0 -2.058e-02
 100.20000 14.060000 101.30000 18.530000     20      20           0  1.094e+00
   0.16600  0.011310   0.17090  0.014050     20      20           0  4.890e-03
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
 -0.007466    -0.66260    1.162        TRUE         ok        ok
  0.010920     0.07782    1.318        TRUE         ok        ok
  0.029450     0.43220    1.242        TRUE         ok        ok

## Downstream slope check, per pair over all cells (5% criterion reported, not gated)

Overall: 23 of 30 pair-cells with a finite relative difference are within 5%.

   dataset           response           predictor y_log x_log n_cells_attempted
 pantheria        body_mass_g head_body_length_mm FALSE FALSE                 6
 pantheria        gestation_d         body_mass_g FALSE FALSE                 6
 pantheria    max_longevity_m         body_mass_g FALSE FALSE                 6
    avonet        Wing.Length                Mass  TRUE  TRUE                 3
    avonet      Tarsus.Length                Mass  TRUE  TRUE                 3
    avonet Beak.Length_Culmen                Mass  TRUE  TRUE                 3
  fishbase             Weight              Length  TRUE  TRUE                 1
  fishbase     DepthRangeDeep              Length FALSE  TRUE                 1
  fishbase              Troph              Length FALSE  TRUE                 1
 n_cells_ok_slope n_cells_degraded has_ok_slope n_within_5pct n_rel_finite
                6                0         TRUE             6            6
                6                0         TRUE             2            6
                6                0         TRUE             3            6
                3                0         TRUE             3            3
                3                0         TRUE             3            3
                3                0         TRUE             3            3
                1                0         TRUE             1            1
                1                0         TRUE             1            1
                1                0         TRUE             1            1
 all_within_5pct max_abs_rel_diff max_abs_diff max_abs_diff_ref_se
            TRUE         0.032810     0.093770             4.04600
           FALSE         0.221600     0.014610             2.27500
           FALSE         0.174500     0.031380             3.28800
            TRUE         0.025750     0.008523             1.22000
            TRUE         0.033990     0.011100             1.68000
            TRUE         0.007511     0.002701             0.36180
            TRUE         0.007466     0.020580             0.66260
            TRUE         0.010920     1.094000             0.07782
            TRUE         0.029450     0.004890             0.43220

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
  fishbase structured 20260818  fishbase-structured-m20260818    1.007   830.2
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
      TRUE
