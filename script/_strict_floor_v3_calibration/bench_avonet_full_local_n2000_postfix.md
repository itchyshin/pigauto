# AVONET full n=9,993 x pigauto (local Mac MPS)

n = 2000 species x 7 traits (bundled avonet_full + tree_full).
Seed = 2026, miss_frac = 0.30, n_imputations = 5.

## Per-trait metrics (pigauto_default vs mean_baseline)

```
          method              trait               metric        value n_cells
   mean_baseline               Mass                 rmse 2085.3933927     600
   mean_baseline               Mass            pearson_r           NA     600
   mean_baseline Beak.Length_Culmen                 rmse   20.7154799     600
   mean_baseline Beak.Length_Culmen            pearson_r           NA     600
   mean_baseline      Tarsus.Length                 rmse   36.0474140     600
   mean_baseline      Tarsus.Length            pearson_r           NA     600
   mean_baseline        Wing.Length                 rmse  106.9809056     600
   mean_baseline        Wing.Length            pearson_r           NA     600
   mean_baseline      Trophic.Level             accuracy    0.5333333     600
   mean_baseline  Primary.Lifestyle             accuracy    0.5783333     600
   mean_baseline          Migration             accuracy    0.8213689     599
 pigauto_default               Mass                 rmse 2460.8804342     600
 pigauto_default               Mass            pearson_r    0.7194757     600
 pigauto_default               Mass coverage95_conformal    0.9716667     600
 pigauto_default               Mass coverage95_mcdropout    0.6300000     600
 pigauto_default Beak.Length_Culmen                 rmse    8.1765194     600
 pigauto_default Beak.Length_Culmen            pearson_r    0.9193716     600
 pigauto_default Beak.Length_Culmen coverage95_conformal    0.9366667     600
 pigauto_default Beak.Length_Culmen coverage95_mcdropout    0.6916667     600
 pigauto_default      Tarsus.Length                 rmse   12.1712817     600
 pigauto_default      Tarsus.Length            pearson_r    0.9476972     600
 pigauto_default      Tarsus.Length coverage95_conformal    0.9583333     600
 pigauto_default      Tarsus.Length coverage95_mcdropout    0.6633333     600
 pigauto_default        Wing.Length                 rmse   33.6365000     600
 pigauto_default        Wing.Length            pearson_r    0.9516139     600
 pigauto_default        Wing.Length coverage95_conformal    0.9033333     600
 pigauto_default        Wing.Length coverage95_mcdropout    0.7700000     600
 pigauto_default      Trophic.Level             accuracy    0.8150000     600
 pigauto_default  Primary.Lifestyle             accuracy    0.8433333     600
 pigauto_default          Migration             accuracy    0.7946578     599
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
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
 608.755
```
