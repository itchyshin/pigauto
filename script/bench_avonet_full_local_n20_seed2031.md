# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 1500 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2031, miss_frac = 0.30, n_imputations = 20.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric       value n_cells
   mean_baseline               Mass                 rmse 676.3354368     450
   mean_baseline               Mass            pearson_r          NA     450
   mean_baseline Beak.Length_Culmen                 rmse  22.1370168     450
   mean_baseline Beak.Length_Culmen            pearson_r          NA     450
   mean_baseline      Tarsus.Length                 rmse  33.7868617     450
   mean_baseline      Tarsus.Length            pearson_r          NA     450
   mean_baseline        Wing.Length                 rmse  93.9752068     450
   mean_baseline        Wing.Length            pearson_r          NA     450
   mean_baseline      Trophic.Level             accuracy   0.5355556     450
   mean_baseline  Primary.Lifestyle             accuracy   0.5800000     450
   mean_baseline          Migration             accuracy   0.7977778     450
 pigauto_default               Mass                 rmse 292.7489552     450
 pigauto_default               Mass            pearson_r   0.9063355     450
 pigauto_default               Mass coverage95_conformal   0.9400000     450
 pigauto_default               Mass coverage95_mcdropout   0.8866667     450
 pigauto_default Beak.Length_Culmen                 rmse   9.0160230     450
 pigauto_default Beak.Length_Culmen            pearson_r   0.9225220     450
 pigauto_default Beak.Length_Culmen coverage95_conformal   0.9911111     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout   0.9222222     450
 pigauto_default      Tarsus.Length                 rmse  21.3358583     450
 pigauto_default      Tarsus.Length            pearson_r   0.7763619     450
 pigauto_default      Tarsus.Length coverage95_conformal   0.9155556     450
 pigauto_default      Tarsus.Length coverage95_mcdropout   0.8622222     450
 pigauto_default        Wing.Length                 rmse  27.3102326     450
 pigauto_default        Wing.Length            pearson_r   0.9566455     450
 pigauto_default        Wing.Length coverage95_conformal   0.9577778     450
 pigauto_default        Wing.Length coverage95_mcdropout   0.8777778     450
 pigauto_default      Trophic.Level             accuracy   0.8311111     450
 pigauto_default  Primary.Lifestyle             accuracy   0.8422222     450
 pigauto_default          Migration             accuracy   0.7311111     450
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
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
 523.913
```
