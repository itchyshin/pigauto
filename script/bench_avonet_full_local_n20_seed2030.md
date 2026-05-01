# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 1500 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2030, miss_frac = 0.30, n_imputations = 20.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric        value n_cells
   mean_baseline               Mass                 rmse 1.819669e+03     450
   mean_baseline               Mass            pearson_r           NA     450
   mean_baseline Beak.Length_Culmen                 rmse 2.416229e+01     450
   mean_baseline Beak.Length_Culmen            pearson_r           NA     450
   mean_baseline      Tarsus.Length                 rmse 2.075521e+01     450
   mean_baseline      Tarsus.Length            pearson_r           NA     450
   mean_baseline        Wing.Length                 rmse 9.462514e+01     450
   mean_baseline        Wing.Length            pearson_r           NA     450
   mean_baseline      Trophic.Level             accuracy 5.777778e-01     450
   mean_baseline  Primary.Lifestyle             accuracy 5.822222e-01     450
   mean_baseline          Migration             accuracy 7.888889e-01     450
 pigauto_default               Mass                 rmse 1.995988e+04     450
 pigauto_default               Mass            pearson_r 9.148541e-01     450
 pigauto_default               Mass coverage95_conformal 9.866667e-01     450
 pigauto_default               Mass coverage95_mcdropout 8.822222e-01     450
 pigauto_default Beak.Length_Culmen                 rmse 1.262309e+01     450
 pigauto_default Beak.Length_Culmen            pearson_r 8.537535e-01     450
 pigauto_default Beak.Length_Culmen coverage95_conformal 9.022222e-01     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout 8.733333e-01     450
 pigauto_default      Tarsus.Length                 rmse 8.249088e+00     450
 pigauto_default      Tarsus.Length            pearson_r 9.285297e-01     450
 pigauto_default      Tarsus.Length coverage95_conformal 9.622222e-01     450
 pigauto_default      Tarsus.Length coverage95_mcdropout 8.911111e-01     450
 pigauto_default        Wing.Length                 rmse 3.092297e+01     450
 pigauto_default        Wing.Length            pearson_r 9.466659e-01     450
 pigauto_default        Wing.Length coverage95_conformal 9.711111e-01     450
 pigauto_default        Wing.Length coverage95_mcdropout 9.422222e-01     450
 pigauto_default      Trophic.Level             accuracy 8.155556e-01     450
 pigauto_default  Primary.Lifestyle             accuracy 8.044444e-01     450
 pigauto_default          Migration             accuracy 6.755556e-01     450
  wall_s
   0.012
   0.012
   0.012
   0.012
   0.012
   0.012
   0.012
   0.012
   0.012
   0.012
   0.012
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
 515.112
```
