# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 1500 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2027, miss_frac = 0.30, n_imputations = 5.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric       value n_cells
   mean_baseline               Mass                 rmse 872.1299118     450
   mean_baseline               Mass            pearson_r          NA     450
   mean_baseline Beak.Length_Culmen                 rmse  17.5987042     450
   mean_baseline Beak.Length_Culmen            pearson_r          NA     450
   mean_baseline      Tarsus.Length                 rmse  28.8425825     450
   mean_baseline      Tarsus.Length            pearson_r          NA     450
   mean_baseline        Wing.Length                 rmse  85.1391383     450
   mean_baseline        Wing.Length            pearson_r          NA     450
   mean_baseline      Trophic.Level             accuracy   0.5866667     450
   mean_baseline  Primary.Lifestyle             accuracy   0.5444444     450
   mean_baseline          Migration             accuracy   0.7639198     449
 pigauto_default               Mass                 rmse 331.9215749     450
 pigauto_default               Mass            pearson_r   0.9340224     450
 pigauto_default               Mass coverage95_conformal   0.9822222     450
 pigauto_default               Mass coverage95_mcdropout   0.6977778     450
 pigauto_default Beak.Length_Culmen                 rmse  10.8582567     450
 pigauto_default Beak.Length_Culmen            pearson_r   0.8040792     450
 pigauto_default Beak.Length_Culmen coverage95_conformal   0.9466667     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout   0.7066667     450
 pigauto_default      Tarsus.Length                 rmse  21.4593883     450
 pigauto_default      Tarsus.Length            pearson_r   0.6703622     450
 pigauto_default      Tarsus.Length coverage95_conformal   0.9266667     450
 pigauto_default      Tarsus.Length coverage95_mcdropout   0.7022222     450
 pigauto_default        Wing.Length                 rmse  25.9162765     450
 pigauto_default        Wing.Length            pearson_r   0.9530096     450
 pigauto_default        Wing.Length coverage95_conformal   0.9533333     450
 pigauto_default        Wing.Length coverage95_mcdropout   0.7022222     450
 pigauto_default      Trophic.Level             accuracy   0.8311111     450
 pigauto_default  Primary.Lifestyle             accuracy   0.8088889     450
 pigauto_default          Migration             accuracy   0.7438753     449
  wall_s
   0.010
   0.010
   0.010
   0.010
   0.010
   0.010
   0.010
   0.010
   0.010
   0.010
   0.010
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
 496.467
```
