# Phase 8 MVP: AVONET 300 head-to-head (pigauto vs BACE)

Seed = 2026, miss_frac = 0.30, identical splits across methods.
**BACE skipped** (not installed or failed).

## Per-trait metrics

```
          method              trait    metric       value n_cells  wall_s
 pigauto_default               Mass      rmse 297.2611782      49 236.836
 pigauto_default               Mass pearson_r   0.9369211      49 236.836
 pigauto_default Beak.Length_Culmen      rmse  14.0765465      45 236.836
 pigauto_default Beak.Length_Culmen pearson_r   0.9323617      45 236.836
 pigauto_default      Tarsus.Length      rmse  20.0862436      48 236.836
 pigauto_default      Tarsus.Length pearson_r   0.9838924      48 236.836
 pigauto_default        Wing.Length      rmse  23.7708165      46 236.836
 pigauto_default        Wing.Length pearson_r   0.9677502      46 236.836
 pigauto_default      Trophic.Level  accuracy   0.8048780      41 236.836
 pigauto_default  Primary.Lifestyle  accuracy   0.7317073      41 236.836
 pigauto_default          Migration  accuracy   0.7777778      45 236.836
     pigauto_em5               Mass      rmse 311.5829208      49 440.491
     pigauto_em5               Mass pearson_r   0.9353836      49 440.491
     pigauto_em5 Beak.Length_Culmen      rmse  14.3740549      45 440.491
     pigauto_em5 Beak.Length_Culmen pearson_r   0.9277050      45 440.491
     pigauto_em5      Tarsus.Length      rmse  16.9755768      48 440.491
     pigauto_em5      Tarsus.Length pearson_r   0.9857766      48 440.491
     pigauto_em5        Wing.Length      rmse  23.6461224      46 440.491
     pigauto_em5        Wing.Length pearson_r   0.9679577      46 440.491
     pigauto_em5      Trophic.Level  accuracy   0.6585366      41 440.491
     pigauto_em5  Primary.Lifestyle  accuracy   0.6585366      41 440.491
     pigauto_em5          Migration  accuracy   0.7777778      45 440.491
```
