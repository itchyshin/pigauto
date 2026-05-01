# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 9993 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2026, miss_frac = 0.30, n_imputations = 20.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric        value n_cells
   mean_baseline               Mass                 rmse 1222.9605832    2998
   mean_baseline               Mass            pearson_r           NA    2998
   mean_baseline Beak.Length_Culmen                 rmse   24.3568713    2998
   mean_baseline Beak.Length_Culmen            pearson_r           NA    2998
   mean_baseline      Tarsus.Length                 rmse   22.1936918    2998
   mean_baseline      Tarsus.Length            pearson_r           NA    2998
   mean_baseline        Wing.Length                 rmse   95.8444425    2998
   mean_baseline        Wing.Length            pearson_r           NA    2998
   mean_baseline      Trophic.Level             accuracy    0.5448782    2997
   mean_baseline  Primary.Lifestyle             accuracy    0.5833889    2998
   mean_baseline          Migration             accuracy    0.8068182    2992
 pigauto_default               Mass                 rmse  678.9934636    2998
 pigauto_default               Mass            pearson_r    0.8554686    2998
 pigauto_default               Mass coverage95_conformal    0.9222815    2998
 pigauto_default               Mass coverage95_mcdropout    0.8662442    2998
 pigauto_default Beak.Length_Culmen                 rmse    7.1709845    2998
 pigauto_default Beak.Length_Culmen            pearson_r    0.9559389    2998
 pigauto_default Beak.Length_Culmen coverage95_conformal    0.9402935    2998
 pigauto_default Beak.Length_Culmen coverage95_mcdropout    0.8735824    2998
 pigauto_default      Tarsus.Length                 rmse    5.8572952    2998
 pigauto_default      Tarsus.Length            pearson_r    0.9649227    2998
 pigauto_default      Tarsus.Length coverage95_conformal    0.9586391    2998
 pigauto_default      Tarsus.Length coverage95_mcdropout    0.8942628    2998
 pigauto_default        Wing.Length                 rmse   21.5478495    2998
 pigauto_default        Wing.Length            pearson_r    0.9745037    2998
 pigauto_default        Wing.Length coverage95_conformal    0.9256171    2998
 pigauto_default        Wing.Length coverage95_mcdropout    0.8949300    2998
 pigauto_default      Trophic.Level             accuracy    0.8498498    2997
 pigauto_default  Primary.Lifestyle             accuracy    0.8689126    2998
 pigauto_default          Migration             accuracy    0.8368984    2992
   wall_s
    0.006
    0.006
    0.006
    0.006
    0.006
    0.006
    0.006
    0.006
    0.006
    0.006
    0.006
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
 9416.104
```
