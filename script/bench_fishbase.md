# FishBase + fishtree x pigauto + BACE

n = 500 species x 6 traits (fishtree x FishBase matched).
Seed = 2026, miss_frac = 0.30.
**BACE skipped** (not installed or opted out).

## Per-trait metrics

```
          method          trait               metric        value n_cells
   mean_baseline         Length                 rmse 3.833280e+01     146
   mean_baseline         Length            pearson_r           NA     146
   mean_baseline         Weight                 rmse 1.701923e+04      26
   mean_baseline         Weight            pearson_r           NA      26
   mean_baseline     BodyShapeI             accuracy 4.200000e-01     150
   mean_baseline DepthRangeDeep                 rmse 1.022985e+03      66
   mean_baseline DepthRangeDeep            pearson_r           NA      66
   mean_baseline  Vulnerability                 rmse 1.696239e+01     150
   mean_baseline  Vulnerability            pearson_r           NA     150
   mean_baseline          Troph                 rmse 6.609146e-01      66
   mean_baseline          Troph            pearson_r           NA      66
 pigauto_default         Length                 rmse 3.332091e+01     146
 pigauto_default         Length            pearson_r 6.608597e-01     146
 pigauto_default         Length coverage95_conformal 9.863014e-01     146
 pigauto_default         Length coverage95_mcdropout 8.424658e-01     146
 pigauto_default         Weight                 rmse 6.950333e+04      26
 pigauto_default         Weight            pearson_r 8.231800e-03      26
 pigauto_default         Weight coverage95_conformal 9.230769e-01      26
 pigauto_default         Weight coverage95_mcdropout 7.307692e-01      26
 pigauto_default     BodyShapeI             accuracy 6.000000e-01     150
 pigauto_default DepthRangeDeep                 rmse 8.693877e+02      66
 pigauto_default DepthRangeDeep            pearson_r 5.407489e-01      66
 pigauto_default DepthRangeDeep coverage95_conformal 1.000000e+00      66
 pigauto_default DepthRangeDeep coverage95_mcdropout 7.878788e-01      66
 pigauto_default  Vulnerability                 rmse 1.151914e+01     150
 pigauto_default  Vulnerability            pearson_r 7.373948e-01     150
 pigauto_default  Vulnerability coverage95_conformal 1.000000e+00     150
 pigauto_default  Vulnerability coverage95_mcdropout 8.333333e-01     150
 pigauto_default          Troph                 rmse 5.955235e-01      66
 pigauto_default          Troph            pearson_r 4.897656e-01      66
 pigauto_default          Troph coverage95_conformal 9.393939e-01      66
 pigauto_default          Troph coverage95_mcdropout 3.484848e-01      66
  wall_s
   0.004
   0.004
   0.004
   0.004
   0.004
   0.004
   0.004
   0.004
   0.004
   0.004
   0.004
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
 256.701
```
