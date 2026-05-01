# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 1500 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2032, miss_frac = 0.30, n_imputations = 20.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric        value n_cells
   mean_baseline               Mass                 rmse 5268.5303553     450
   mean_baseline               Mass            pearson_r           NA     450
   mean_baseline Beak.Length_Culmen                 rmse   31.9209670     450
   mean_baseline Beak.Length_Culmen            pearson_r           NA     450
   mean_baseline      Tarsus.Length                 rmse   24.1038072     450
   mean_baseline      Tarsus.Length            pearson_r           NA     450
   mean_baseline        Wing.Length                 rmse   95.1762099     450
   mean_baseline        Wing.Length            pearson_r           NA     450
   mean_baseline      Trophic.Level             accuracy    0.5844444     450
   mean_baseline  Primary.Lifestyle             accuracy    0.5622222     450
   mean_baseline          Migration             accuracy    0.8133333     450
 pigauto_default               Mass                 rmse 2713.6276482     450
 pigauto_default               Mass            pearson_r    0.9920232     450
 pigauto_default               Mass coverage95_conformal    0.9777778     450
 pigauto_default               Mass coverage95_mcdropout    0.9288889     450
 pigauto_default Beak.Length_Culmen                 rmse   15.9178899     450
 pigauto_default Beak.Length_Culmen            pearson_r    0.8914881     450
 pigauto_default Beak.Length_Culmen coverage95_conformal    0.9622222     450
 pigauto_default Beak.Length_Culmen coverage95_mcdropout    0.8955556     450
 pigauto_default      Tarsus.Length                 rmse   10.1303878     450
 pigauto_default      Tarsus.Length            pearson_r    0.9077691     450
 pigauto_default      Tarsus.Length coverage95_conformal    0.9577778     450
 pigauto_default      Tarsus.Length coverage95_mcdropout    0.8911111     450
 pigauto_default        Wing.Length                 rmse   30.1493604     450
 pigauto_default        Wing.Length            pearson_r    0.9485843     450
 pigauto_default        Wing.Length coverage95_conformal    0.9577778     450
 pigauto_default        Wing.Length coverage95_mcdropout    0.9266667     450
 pigauto_default      Trophic.Level             accuracy    0.8288889     450
 pigauto_default  Primary.Lifestyle             accuracy    0.7844444     450
 pigauto_default          Migration             accuracy    0.7311111     450
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
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
 503.039
```
