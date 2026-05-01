# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 1500 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2026, miss_frac = 0.30, n_imputations = 5.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric       value n_cells
   mean_baseline               Mass                 rmse 853.5424927     450
   mean_baseline               Mass            pearson_r          NA     450
   mean_baseline Beak.Length_Culmen                 rmse  19.1112856     450
   mean_baseline Beak.Length_Culmen            pearson_r          NA     450
   mean_baseline      Tarsus.Length                 rmse  27.6360367     450
   mean_baseline      Tarsus.Length            pearson_r          NA     450
   mean_baseline        Wing.Length                 rmse 110.8497462     450
   mean_baseline        Wing.Length            pearson_r          NA     450
   mean_baseline      Trophic.Level             accuracy   0.5466667     450
   mean_baseline  Primary.Lifestyle             accuracy   0.5933333     450
   mean_baseline          Migration             accuracy   0.8195991     449
 pigauto_default               Mass                 rmse 659.0100823     450
 pigauto_default               Mass            pearson_r   0.7946685     450
 pigauto_default               Mass coverage95_conformal   0.9311111     450
 pigauto_default               Mass coverage95_mcdropout   0.7022222     450
 pigauto_default Beak.Length_Culmen                 rmse   6.3864151     450
 pigauto_default Beak.Length_Culmen            pearson_r   0.9435588     450
 pigauto_default Beak.Length_Culmen coverage95_conformal   0.9600000     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout   0.6666667     450
 pigauto_default      Tarsus.Length                 rmse  12.0439824     450
 pigauto_default      Tarsus.Length            pearson_r   0.9169316     450
 pigauto_default      Tarsus.Length coverage95_conformal   0.9200000     450
 pigauto_default      Tarsus.Length coverage95_mcdropout   0.7066667     450
 pigauto_default        Wing.Length                 rmse  35.1059552     450
 pigauto_default        Wing.Length            pearson_r   0.9509387     450
 pigauto_default        Wing.Length coverage95_conformal   0.9000000     450
 pigauto_default        Wing.Length coverage95_mcdropout   0.6955556     450
 pigauto_default      Trophic.Level             accuracy   0.8288889     450
 pigauto_default  Primary.Lifestyle             accuracy   0.8377778     450
 pigauto_default          Migration             accuracy   0.8151448     449
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
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
 501.728
```
