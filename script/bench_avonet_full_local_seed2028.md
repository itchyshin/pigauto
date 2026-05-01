# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 1500 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2028, miss_frac = 0.30, n_imputations = 5.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric        value n_cells
   mean_baseline               Mass                 rmse 1396.5139641     450
   mean_baseline               Mass            pearson_r           NA     450
   mean_baseline Beak.Length_Culmen                 rmse   24.3437869     450
   mean_baseline Beak.Length_Culmen            pearson_r           NA     450
   mean_baseline      Tarsus.Length                 rmse   29.0702456     450
   mean_baseline      Tarsus.Length            pearson_r           NA     450
   mean_baseline        Wing.Length                 rmse  105.3888459     450
   mean_baseline        Wing.Length            pearson_r           NA     450
   mean_baseline      Trophic.Level             accuracy    0.5377778     450
   mean_baseline  Primary.Lifestyle             accuracy    0.5800000     450
   mean_baseline          Migration             accuracy    0.7777778     450
 pigauto_default               Mass                 rmse 1130.9987954     450
 pigauto_default               Mass            pearson_r    0.5873000     450
 pigauto_default               Mass coverage95_conformal    0.9844444     450
 pigauto_default               Mass coverage95_mcdropout    0.6933333     450
 pigauto_default Beak.Length_Culmen                 rmse    8.8429289     450
 pigauto_default Beak.Length_Culmen            pearson_r    0.9418040     450
 pigauto_default Beak.Length_Culmen coverage95_conformal    0.9600000     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout    0.6955556     450
 pigauto_default      Tarsus.Length                 rmse   15.4202563     450
 pigauto_default      Tarsus.Length            pearson_r    0.8481012     450
 pigauto_default      Tarsus.Length coverage95_conformal    0.9511111     450
 pigauto_default      Tarsus.Length coverage95_mcdropout    0.6822222     450
 pigauto_default        Wing.Length                 rmse   42.6012707     450
 pigauto_default        Wing.Length            pearson_r    0.9180308     450
 pigauto_default        Wing.Length coverage95_conformal    0.9111111     450
 pigauto_default        Wing.Length coverage95_mcdropout    0.6622222     450
 pigauto_default      Trophic.Level             accuracy    0.8177778     450
 pigauto_default  Primary.Lifestyle             accuracy    0.7800000     450
 pigauto_default          Migration             accuracy    0.7688889     450
  wall_s
   0.011
   0.011
   0.011
   0.011
   0.011
   0.011
   0.011
   0.011
   0.011
   0.011
   0.011
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
 486.085
```
