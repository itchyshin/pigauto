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
 pigauto_default               Mass                 rmse 377.2854278     450
 pigauto_default               Mass            pearson_r   0.8975241     450
 pigauto_default               Mass coverage95_conformal   0.9244444     450
 pigauto_default               Mass coverage95_mcdropout   0.7222222     450
 pigauto_default Beak.Length_Culmen                 rmse   6.5533521     450
 pigauto_default Beak.Length_Culmen            pearson_r   0.9390982     450
 pigauto_default Beak.Length_Culmen coverage95_conformal   0.9577778     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout   0.6666667     450
 pigauto_default      Tarsus.Length                 rmse  11.9181553     450
 pigauto_default      Tarsus.Length            pearson_r   0.9255312     450
 pigauto_default      Tarsus.Length coverage95_conformal   0.9133333     450
 pigauto_default      Tarsus.Length coverage95_mcdropout   0.6422222     450
 pigauto_default        Wing.Length                 rmse  34.1175818     450
 pigauto_default        Wing.Length            pearson_r   0.9515771     450
 pigauto_default        Wing.Length coverage95_conformal   0.9222222     450
 pigauto_default        Wing.Length coverage95_mcdropout   0.6977778     450
 pigauto_default      Trophic.Level             accuracy   0.8288889     450
 pigauto_default  Primary.Lifestyle             accuracy   0.8377778     450
 pigauto_default          Migration             accuracy   0.8151448     449
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
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
 474.474
```
