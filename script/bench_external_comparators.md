# pigauto vs external comparators -- AVONET300

## Regime

- Dataset: bundled AVONET300 (n = 300 species, single-obs).
- Scored traits (4 continuous, all methods): Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length.
- pigauto is fit on all 7 AVONET traits (its unified mixed-type
  imputation across continuous/categorical/ordinal is the point of
  the package) but scored only on the 4 continuous traits above,
  for comparability with the continuous-only rivals.
- Missingness: 30% MCAR on observed cells, 5 reps, mask seeds 2027..2031.
- Metric: z-scored RMSE (scale = training-portion mean/sd per
  trait/rep, no leakage from held-out truth) + Pearson r on
  held-out cells. Reported as mean +/- MCSE (sd/sqrt(n_reps_ok)).

## Methods

1. `pigauto::impute()` -- production defaults (all 7 traits fit; scored on 4).
2. `Rphylopars::phylopars(model = "BM")` standalone -- continuous traits only.
3. `phylolm(model = "lambda")` per trait, no covariates (y ~ 1) --
   phylogenetic BLUP, mirrors the arm in `script/bench_lambda_sweep.R`
   (`m_phylolm_blup`, ~line 155) with the fixed-effect covariates dropped.
4. Column mean -- floor.

## IMPORTANT: what the pigauto-vs-phylopars comparison actually measures

`Rphylopars::phylopars()` is ALSO the solver pigauto calls internally
for its joint-baseline path (`R/joint_mvn_baseline.R`,
`R/joint_threshold_baseline.R`, via `R/joint_mvn_solver.R`, both using
`model = "BM"`). So `pigauto` vs `rphylopars` here is NOT a comparison
against an unrelated package -- it measures what pigauto's GNN + gating
+ conformal layer adds on top of the raw phylopars solver it can
already call as its own baseline.

## Per-trait summary (mean +/- MCSE across reps)

```
   dataset              trait         method n_reps_ok          rmse_z
 avonet300 Beak.Length_Culmen    column_mean         5 1.204 +/- 0.196
 avonet300 Beak.Length_Culmen phylolm_lambda         5 0.672 +/- 0.083
 avonet300 Beak.Length_Culmen        pigauto         5 0.912 +/- 0.162
 avonet300 Beak.Length_Culmen     rphylopars         5 0.445 +/- 0.086
 avonet300               Mass    column_mean         5 1.734 +/- 0.491
 avonet300               Mass phylolm_lambda         5 1.629 +/- 0.473
 avonet300               Mass        pigauto         5 1.594 +/- 0.508
 avonet300               Mass     rphylopars         5 1.360 +/- 0.341
 avonet300      Tarsus.Length    column_mean         5 1.516 +/- 0.337
 avonet300      Tarsus.Length phylolm_lambda         5 1.046 +/- 0.213
 avonet300      Tarsus.Length        pigauto         5 1.220 +/- 0.270
 avonet300      Tarsus.Length     rphylopars         5 0.639 +/- 0.120
 avonet300        Wing.Length    column_mean         5 1.120 +/- 0.113
 avonet300        Wing.Length phylolm_lambda         5 0.675 +/- 0.105
 avonet300        Wing.Length        pigauto         5 0.688 +/- 0.083
 avonet300        Wing.Length     rphylopars         5 0.409 +/- 0.066
       pearson_r wall_s_mean errors
       NA +/- NA      0.0004   <NA>
 0.840 +/- 0.036      0.1108   <NA>
 0.738 +/- 0.045    877.7496   <NA>
 0.920 +/- 0.027      1.0168   <NA>
       NA +/- NA      0.0004   <NA>
 0.424 +/- 0.064      0.1108   <NA>
 0.579 +/- 0.131    877.7496   <NA>
 0.778 +/- 0.028      1.0168   <NA>
       NA +/- NA      0.0004   <NA>
 0.744 +/- 0.053      0.1108   <NA>
 0.713 +/- 0.076    877.7496   <NA>
 0.917 +/- 0.020      1.0168   <NA>
       NA +/- NA      0.0004   <NA>
 0.804 +/- 0.031      0.1108   <NA>
 0.833 +/- 0.035    877.7496   <NA>
 0.935 +/- 0.011      1.0168   <NA>
```

## Rival failures

None -- all methods produced output for all traits/reps.

Total wall time: 73.2 min.
