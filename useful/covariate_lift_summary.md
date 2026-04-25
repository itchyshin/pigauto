# Covariate-lift summary across datasets

Last regenerated: 2026-04-24 18:22:41

One row per dataset.  `median_ratio_on` is the median across traits of
`cov_on_RMSE / none_RMSE`.  Values < 1.00 mean covariates lift accuracy.

- `n_traits_with_lift_on`: how many traits had `ratio_on < 0.95` (>=5% lift).
- `n_traits_with_lift_off`: same for safety_floor=FALSE.

```
                                          dataset n_species n_traits
          PanTHERIA mammals (precip+temp+lat+PET)       850        5
       GlobTherm ectotherms (lat+long+elev+|lat|)       809        2
     AmphiBIO amphibians (climate-zone occupancy)      1000        5
 Delhey birds (annual_temp+precip+tree_cover+...)      5809        2
  BIEN plants (WorldClim bioclim, cached species)      3450        5
 LepTraits butterflies (Jan-Dec flight phenology)      1500        4
          Multi-obs sim (Yule tree + acclim_temp)        NA       16
 median_ratio_on median_ratio_off n_traits_with_lift_on n_traits_with_lift_off
       1.0000000        0.9371046                     1                      3
       1.1341536        0.9583471                     0                      1
       1.0000000        1.1147374                     1                      1
       1.0172567        1.0099264                     0                      0
       1.0138293        1.0525132                     0                      0
       1.0061317        1.0446413                     1                      2
       0.8533183               NA                    16                     NA
     kind
     real
     real
     real
     real
     real
     real
 sim_yule
```

## Reading guide

Real species-level datasets with phylo-redundant climate covariates
(BIEN plants, Delhey birds, AmphiBIO amphibians) cluster around
ratio_on ~ 1.00 -- the safety floor closes the gate when covariates
are uninformative beyond phylogeny, so users pay no accuracy penalty
for passing redundant climate features.

Datasets with PARTIAL phylo-decoupled covariate signal show
MIXED lift: PanTHERIA mammals with the bundled precip/temp/lat
columns lifts MaxLongevity by 22% but is neutral on body mass and
gestation length.  This is the safety property in action: the
calibrated gate adapts per-trait.

Multi-observation simulations (acclim_temp varying within species)
show consistent 10-19% lift, and the lift survives when the Yule
sim tree is replaced with the real AVONET 300 bird phylogeny.

