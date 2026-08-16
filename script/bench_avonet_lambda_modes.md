# AVONET300: pigauto lambda_mode fixed_1 vs bayes

n = 300 species x 7 traits, seed = 2026, MCAR 30%, m = 20.
Same data/seed/mask as `script/bench_avonet_bace.md`, so BACE's numbers there are
directly comparable. Single seed -- a direction, not an interval.

```
          method              trait    metric        value n_cells  wall_s
 pigauto_fixed_1               Mass      rmse 2625.3151262      90 913.359
 pigauto_fixed_1               Mass pearson_r    0.1477451      90 913.359
 pigauto_fixed_1 Beak.Length_Culmen      rmse   22.1924846      90 913.359
 pigauto_fixed_1 Beak.Length_Culmen pearson_r    0.8219531      90 913.359
 pigauto_fixed_1      Tarsus.Length      rmse   27.0578539      90 913.359
 pigauto_fixed_1      Tarsus.Length pearson_r    0.8827318      90 913.359
 pigauto_fixed_1        Wing.Length      rmse   67.7668863      90 913.359
 pigauto_fixed_1        Wing.Length pearson_r    0.8789303      90 913.359
 pigauto_fixed_1      Trophic.Level  accuracy    0.7888889      90 913.359
 pigauto_fixed_1  Primary.Lifestyle  accuracy    0.6666667      90 913.359
 pigauto_fixed_1          Migration  accuracy    0.6222222      90 913.359
   pigauto_bayes               Mass      rmse 2301.6948989      90 552.503
   pigauto_bayes               Mass pearson_r    0.8290303      90 552.503
   pigauto_bayes Beak.Length_Culmen      rmse   21.2678199      90 552.503
   pigauto_bayes Beak.Length_Culmen pearson_r    0.8523104      90 552.503
   pigauto_bayes      Tarsus.Length      rmse   24.4169358      90 552.503
   pigauto_bayes      Tarsus.Length pearson_r    0.8314513      90 552.503
   pigauto_bayes        Wing.Length      rmse   67.8870915      90 552.503
   pigauto_bayes        Wing.Length pearson_r    0.8755065      90 552.503
   pigauto_bayes      Trophic.Level  accuracy    0.6000000      90 552.503
   pigauto_bayes  Primary.Lifestyle  accuracy    0.5777778      90 552.503
   pigauto_bayes          Migration  accuracy    0.7555556      90 552.503
```
