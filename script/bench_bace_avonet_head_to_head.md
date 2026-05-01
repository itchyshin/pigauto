# Phase 8 MVP: AVONET 300 head-to-head (pigauto vs BACE)

Seed = 2026, miss_frac = 0.30, identical splits across methods.
**BACE skipped** (not installed or failed).

## Per-trait metrics

```
          method              trait               metric       value n_cells
 pigauto_default               Mass                 rmse 345.2655707      49
 pigauto_default               Mass            pearson_r   0.9325543      49
 pigauto_default               Mass coverage95_mcdropout   0.8979592      49
 pigauto_default Beak.Length_Culmen                 rmse  14.9974919      45
 pigauto_default Beak.Length_Culmen            pearson_r   0.9211402      45
 pigauto_default Beak.Length_Culmen coverage95_mcdropout   0.7777778      45
 pigauto_default      Tarsus.Length                 rmse  18.3932497      48
 pigauto_default      Tarsus.Length            pearson_r   0.9820217      48
 pigauto_default      Tarsus.Length coverage95_mcdropout   0.8125000      48
 pigauto_default        Wing.Length                 rmse  25.0218286      46
 pigauto_default        Wing.Length            pearson_r   0.9636297      46
 pigauto_default        Wing.Length coverage95_mcdropout   0.9130435      46
 pigauto_default      Trophic.Level             accuracy   0.8048780      41
 pigauto_default  Primary.Lifestyle             accuracy   0.7317073      41
 pigauto_default          Migration             accuracy   0.7777778      45
     pigauto_em5               Mass                 rmse 294.5262990      49
     pigauto_em5               Mass            pearson_r   0.9444636      49
     pigauto_em5               Mass coverage95_mcdropout   0.8367347      49
     pigauto_em5 Beak.Length_Culmen                 rmse  14.4413883      45
     pigauto_em5 Beak.Length_Culmen            pearson_r   0.9269645      45
     pigauto_em5 Beak.Length_Culmen coverage95_mcdropout   0.8444444      45
     pigauto_em5      Tarsus.Length                 rmse  18.2349860      48
     pigauto_em5      Tarsus.Length            pearson_r   0.9832165      48
     pigauto_em5      Tarsus.Length coverage95_mcdropout   0.7916667      48
     pigauto_em5        Wing.Length                 rmse  23.8502730      46
     pigauto_em5        Wing.Length            pearson_r   0.9679491      46
     pigauto_em5        Wing.Length coverage95_mcdropout   0.8695652      46
     pigauto_em5      Trophic.Level             accuracy   0.6585366      41
     pigauto_em5  Primary.Lifestyle             accuracy   0.7317073      41
     pigauto_em5          Migration             accuracy   0.7777778      45
  wall_s
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 238.567
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
 468.576
```
