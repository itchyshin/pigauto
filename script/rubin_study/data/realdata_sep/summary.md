# mi-posterior real-data summary (real data)

10/10 planned cells have an "ok" receipt (0 flagged non-converged; see below).

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
         0.9453         0.9597            0.9597
         0.9555         0.9817            0.9712
         0.9336         0.9151            0.9520
         0.9505         0.9455            0.9554
         0.9381         0.9482            0.9496
         0.9215         0.9503            0.9686
         0.9151         0.9299            0.9373
         0.9703         0.9406            0.9505
         0.9453         0.9540            0.9698
         0.9555         0.9634            0.9660
         0.9446         0.9299            0.9299
         0.9208         0.9703            0.9406
         0.9165         0.9324            0.9424
         0.9319         0.9346            0.9372
         0.9373         0.9705            0.9631
         0.9158         0.9752            0.9901
         0.8763         0.9252            0.9424
         0.9372         0.9476            0.9607
         0.9225         0.9520            0.9631
         0.9257         0.9703            0.9604
         0.8734         0.9007            0.9122
         0.9424         0.9319            0.9346
         0.9446         0.9815            0.9705
         0.9455         0.9307            0.9356
         0.9167         0.9633            0.9433
         0.9467         0.9733            0.9633
         0.9300         0.9867            0.9800
         0.9800         0.9533            0.9233
         0.9567         0.9900            0.9867
         0.9700         0.9567            0.9600
         0.9467         0.9533            0.9600
         0.9767         0.9400            0.9700
         0.9700         0.9067            0.9267
         0.9467         0.9633            0.9800
         0.9467         0.9367            0.9533
         0.9767         0.9800            1.0000
         0.9424         0.9394            0.9399
         0.9416         0.9496            0.9788
         0.9574         0.9659            0.9616
         0.9504         0.9542            0.9495
         0.9369         0.9389            0.9462

## Interval widths, like for like (mean and median of upper - lower, original scale)

   dataset        arm     seed               trait model_mean_width
 pantheria       mcar 20260818         body_mass_g        3.046e+00
 pantheria       mcar 20260818 head_body_length_mm        7.563e-01
 pantheria       mcar 20260818         gestation_d        1.010e+00
 pantheria       mcar 20260818     max_longevity_m        1.964e+00
 pantheria       mcar 20260819         body_mass_g        2.642e+00
 pantheria       mcar 20260819 head_body_length_mm        7.479e-01
 pantheria       mcar 20260819         gestation_d        8.904e-01
 pantheria       mcar 20260819     max_longevity_m        1.935e+00
 pantheria       mcar 20260820         body_mass_g        2.974e+00
 pantheria       mcar 20260820 head_body_length_mm        7.677e-01
 pantheria       mcar 20260820         gestation_d        9.408e-01
 pantheria       mcar 20260820     max_longevity_m        1.868e+00
 pantheria structured 20260818         body_mass_g        2.376e+00
 pantheria structured 20260818 head_body_length_mm        7.437e-01
 pantheria structured 20260818         gestation_d        9.506e-01
 pantheria structured 20260818     max_longevity_m        1.905e+00
 pantheria structured 20260819         body_mass_g        2.306e+00
 pantheria structured 20260819 head_body_length_mm        7.581e-01
 pantheria structured 20260819         gestation_d        9.161e-01
 pantheria structured 20260819     max_longevity_m        1.849e+00
 pantheria structured 20260820         body_mass_g        2.223e+00
 pantheria structured 20260820 head_body_length_mm        7.706e-01
 pantheria structured 20260820         gestation_d        9.635e-01
 pantheria structured 20260820     max_longevity_m        1.894e+00
    avonet       mcar 20260818                Mass        3.395e+02
    avonet       mcar 20260818  Beak.Length_Culmen        1.724e+01
    avonet       mcar 20260818       Tarsus.Length        1.790e+01
    avonet       mcar 20260818         Wing.Length        8.826e+01
    avonet       mcar 20260819                Mass        2.688e+02
    avonet       mcar 20260819  Beak.Length_Culmen        1.684e+01
    avonet       mcar 20260819       Tarsus.Length        1.848e+01
    avonet       mcar 20260819         Wing.Length        8.088e+01
    avonet       mcar 20260820                Mass        4.412e+02
    avonet       mcar 20260820  Beak.Length_Culmen        1.921e+01
    avonet       mcar 20260820       Tarsus.Length        1.627e+01
    avonet       mcar 20260820         Wing.Length        9.368e+01
  fishbase structured 20260818              Length        4.015e+01
  fishbase structured 20260818              Weight        1.289e+05
  fishbase structured 20260818      DepthRangeDeep        1.901e+03
  fishbase structured 20260818       Vulnerability        2.719e+01
  fishbase structured 20260818               Troph        1.744e+00
 split_mean_width mondrian_mean_width model_median_width split_median_width
        4.668e+00           4.490e+00             2.3210          3.827e+00
        1.856e+00           1.671e+00             0.6666          1.751e+00
        1.031e+00           1.205e+00             0.9402          1.013e+00
        2.622e+00           2.754e+00             1.9190          2.703e+00
        4.321e+00           3.865e+00             2.2760          3.587e+00
        1.439e+00           1.446e+00             0.6521          1.362e+00
        1.177e+00           1.225e+00             0.8220          1.150e+00
        1.935e+00           1.981e+00             1.9370          1.993e+00
        4.679e+00           4.668e+00             2.4130          3.807e+00
        1.250e+00           1.284e+00             0.6831          1.190e+00
        8.973e-01           9.042e-01             0.9243          8.813e-01
        2.902e+00           2.296e+00             1.8680          2.990e+00
        3.269e+00           3.249e+00             1.8570          3.030e+00
        1.198e+00           1.473e+00             0.6209          1.151e+00
        1.356e+00           1.314e+00             0.8453          1.344e+00
        3.560e+00           3.875e+00             1.8720          3.672e+00
        3.184e+00           3.541e+00             1.7600          2.886e+00
        1.165e+00           1.233e+00             0.6242          1.094e+00
        1.574e+00           1.780e+00             0.7445          1.487e+00
        2.973e+00           2.786e+00             1.8620          3.016e+00
        2.812e+00           2.772e+00             1.7610          2.599e+00
        1.028e+00           1.064e+00             0.6519          9.873e-01
        1.520e+00           1.245e+00             0.8752          1.397e+00
        2.128e+00           2.229e+00             1.8480          2.160e+00
        7.336e+02           8.741e+02            43.0900          1.157e+02
        3.570e+01           3.222e+01            12.2400          2.858e+01
        4.343e+01           8.543e+01            12.0700          3.357e+01
        1.290e+02           1.298e+02            59.1800          1.028e+02
        6.117e+02           6.426e+02            40.4500          1.089e+02
        2.914e+01           2.890e+01            11.7800          2.399e+01
        3.381e+01           3.855e+01            12.2500          2.623e+01
        9.501e+01           1.321e+02            54.6900          7.116e+01
        5.316e+02           7.052e+02            39.1500          8.073e+01
        3.406e+01           3.871e+01            12.9200          2.654e+01
        2.210e+01           2.737e+01            12.0800          1.786e+01
        1.529e+02           2.456e+02            54.5100          1.102e+02
        6.636e+01           6.526e+01            25.4900          4.206e+01
        3.400e+05           3.243e+05          6341.0000          6.422e+04
        2.237e+03           2.365e+03          1815.0000          2.237e+03
        4.645e+01           4.300e+01            23.6600          4.244e+01
        1.807e+00           1.877e+00             1.7370          1.833e+00
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
   2.85800  0.023180   2.91800  0.039550     20      20           0  0.0605800
   0.06592  0.006423   0.06036  0.007641     20      20           0 -0.0055660
   0.17980  0.009542   0.15070  0.011890     20      20           0 -0.0291000
   2.85800  0.023180   2.91700  0.042960     20      20           0  0.0592500
   0.06592  0.006423   0.05695  0.006922     20      20           0 -0.0089720
   0.17980  0.009542   0.17270  0.011810     20      20           0 -0.0070670
   2.85800  0.023180   2.94900  0.043250     20      20           0  0.0907800
   0.06592  0.006423   0.05259  0.007659     20      20           0 -0.0133300
   0.17980  0.009542   0.15300  0.010400     20      20           0 -0.0267400
   2.85800  0.023180   2.84800  0.034430     20      20           0 -0.0097580
   0.06592  0.006423   0.07303  0.007096     20      20           0  0.0071090
   0.17980  0.009542   0.16890  0.011620     20      20           0 -0.0108400
   2.85800  0.023180   2.78500  0.038620     20      20           0 -0.0733100
   0.06592  0.006423   0.06574  0.007116     20      20           0 -0.0001888
   0.17980  0.009542   0.17770  0.010320     20      20           0 -0.0020580
   2.85800  0.023180   2.78900  0.040930     20      20           0 -0.0684400
   0.06592  0.006423   0.06721  0.007406     20      20           0  0.0012870
   0.17980  0.009542   0.17920  0.012560     20      20           0 -0.0006074
   0.33090  0.006988   0.33540  0.009044     20      20           0  0.0045000
   0.32660  0.006610   0.33700  0.007714     20      20           0  0.0103600
   0.35960  0.007465   0.35900  0.008219     20      20           0 -0.0006243
   0.33090  0.006988   0.32510  0.008629     20      20           0 -0.0058260
   0.32660  0.006610   0.32900  0.008120     20      20           0  0.0024030
   0.35960  0.007465   0.35650  0.009526     20      20           0 -0.0031270
   0.33090  0.006988   0.32350  0.008278     20      20           0 -0.0074490
   0.32660  0.006610   0.32680  0.007443     20      20           0  0.0001784
   0.35960  0.007465   0.35880  0.008922     20      20           0 -0.0008272
   2.75600  0.031060   2.74100  0.038020     20      20           0 -0.0155100
 100.20000 14.060000 101.00000 16.760000     20      20           0  0.8339000
   0.16600  0.011310   0.17040  0.013060     20      20           0  0.0043800
   rel_diff diff_ref_se se_ratio within_5pct ref_status mi_status
  0.0212000     2.61400    1.706        TRUE         ok        ok
 -0.0844300    -0.86660    1.190       FALSE         ok        ok
 -0.1619000    -3.05000    1.246       FALSE         ok        ok
  0.0207300     2.55600    1.854        TRUE         ok        ok
 -0.1361000    -1.39700    1.078       FALSE         ok        ok
 -0.0393100    -0.74060    1.238        TRUE         ok        ok
  0.0317700     3.91700    1.866        TRUE         ok        ok
 -0.2023000    -2.07600    1.193       FALSE         ok        ok
 -0.1487000    -2.80200    1.090       FALSE         ok        ok
 -0.0034150    -0.42100    1.486        TRUE         ok        ok
  0.1078000     1.10700    1.105       FALSE         ok        ok
 -0.0603000    -1.13600    1.217       FALSE         ok        ok
 -0.0256500    -3.16300    1.666        TRUE         ok        ok
 -0.0028640    -0.02939    1.108        TRUE         ok        ok
 -0.0114500    -0.21570    1.082        TRUE         ok        ok
 -0.0239500    -2.95300    1.766        TRUE         ok        ok
  0.0195300     0.20050    1.153        TRUE         ok        ok
 -0.0033780    -0.06365    1.316        TRUE         ok        ok
  0.0136000     0.64390    1.294        TRUE         ok        ok
  0.0317300     1.56800    1.167        TRUE         ok        ok
 -0.0017360    -0.08363    1.101        TRUE         ok        ok
 -0.0176000    -0.83360    1.235        TRUE         ok        ok
  0.0073560     0.36350    1.228        TRUE         ok        ok
 -0.0086960    -0.41890    1.276        TRUE         ok        ok
 -0.0225100    -1.06600    1.185        TRUE         ok        ok
  0.0005461     0.02698    1.126        TRUE         ok        ok
 -0.0023000    -0.11080    1.195        TRUE         ok        ok
 -0.0056280    -0.49950    1.224        TRUE         ok        ok
  0.0083260     0.05932    1.192        TRUE         ok        ok
  0.0263800     0.38720    1.154        TRUE         ok        ok

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
            TRUE         0.031770     0.090780             3.91700
           FALSE         0.202300     0.013330             2.07600
           FALSE         0.161900     0.029100             3.05000
            TRUE         0.022510     0.007449             1.06600
            TRUE         0.031730     0.010360             1.56800
            TRUE         0.008696     0.003127             0.41890
            TRUE         0.005628     0.015510             0.49950
            TRUE         0.008326     0.833900             0.05932
            TRUE         0.026380     0.004380             0.38720

## Convergence (max split R-hat / min bulk ESS over Sigma_P, Sigma_E, lambda)

Descriptive only; does not gate coverage or REALDATA_COMPLETE (see 03_acceptance.R).

   dataset        arm     seed                           name max_rhat min_ess
 pantheria       mcar 20260818       pantheria-mcar-m20260818    1.010   469.4
 pantheria       mcar 20260819       pantheria-mcar-m20260819    1.009   459.3
 pantheria       mcar 20260820       pantheria-mcar-m20260820    1.001   750.8
 pantheria structured 20260818 pantheria-structured-m20260818    1.003   404.1
 pantheria structured 20260819 pantheria-structured-m20260819    1.007   654.6
 pantheria structured 20260820 pantheria-structured-m20260820    1.008   627.1
    avonet       mcar 20260818          avonet-mcar-m20260818    1.008   918.0
    avonet       mcar 20260819          avonet-mcar-m20260819    1.012   487.9
    avonet       mcar 20260820          avonet-mcar-m20260820    1.007   625.4
  fishbase structured 20260818  fishbase-structured-m20260818    1.006   688.7
 converged
      TRUE
      TRUE
      TRUE
      TRUE
      TRUE
      TRUE
      TRUE
      TRUE
      TRUE
      TRUE
