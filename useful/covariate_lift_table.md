# Covariate-lift results: per-trait detail

Last regenerated: 2026-04-24 18:22:41

Each row = one (dataset, trait) pair.  RMSE columns are the absolute RMSE on held-out cells.

- `mean_RMSE`: imputing every held-out cell with the column grand mean
- `none_RMSE`: pigauto baseline (no covariates)
- `cov_on_RMSE`:  pigauto with covariates, safety_floor = TRUE
- `cov_off_RMSE`: pigauto with covariates, safety_floor = FALSE
- `ratio_on`:  cov_on_RMSE / none_RMSE.  <1 = covariates help.
- `ratio_off`: cov_off_RMSE / none_RMSE. <1 = covariates help.

Pearson r columns describe the same fits' correlation between predicted and held-out truth.

```
                                          dataset n_species
          PanTHERIA mammals (precip+temp+lat+PET)       850
          PanTHERIA mammals (precip+temp+lat+PET)       850
          PanTHERIA mammals (precip+temp+lat+PET)       850
          PanTHERIA mammals (precip+temp+lat+PET)       850
          PanTHERIA mammals (precip+temp+lat+PET)       850
       GlobTherm ectotherms (lat+long+elev+|lat|)       809
       GlobTherm ectotherms (lat+long+elev+|lat|)       809
     AmphiBIO amphibians (climate-zone occupancy)      1000
     AmphiBIO amphibians (climate-zone occupancy)      1000
     AmphiBIO amphibians (climate-zone occupancy)      1000
     AmphiBIO amphibians (climate-zone occupancy)      1000
     AmphiBIO amphibians (climate-zone occupancy)      1000
 Delhey birds (annual_temp+precip+tree_cover+...)      5809
 Delhey birds (annual_temp+precip+tree_cover+...)      5809
  BIEN plants (WorldClim bioclim, cached species)      3450
  BIEN plants (WorldClim bioclim, cached species)      3450
  BIEN plants (WorldClim bioclim, cached species)      3450
  BIEN plants (WorldClim bioclim, cached species)      3450
  BIEN plants (WorldClim bioclim, cached species)      3450
 LepTraits butterflies (Jan-Dec flight phenology)      1500
 LepTraits butterflies (Jan-Dec flight phenology)      1500
 LepTraits butterflies (Jan-Dec flight phenology)      1500
 LepTraits butterflies (Jan-Dec flight phenology)      1500
                         trait n_held    mean_RMSE    none_RMSE cov_on_RMSE
          X5.1_AdultBodyMass_g    200 3.013661e+00 1.358738e+00   1.4155550
           X9.1_GestationLen_d     80 8.696069e-01 3.982843e-01   0.4344477
          X17.1_MaxLongevity_m     62 9.809042e-01 8.376791e-01   0.6515648
              X15.1_LitterSize    142 1.845225e+00 1.869379e+00   1.8693788
 X21.1_PopulationDensity_n.km2     50 3.382826e+00 2.229402e+00   2.1450661
                          Tmax    243 6.627142e+00 5.302985e+00   6.8067738
                          tmin    109 7.786461e+00 4.738468e+00   4.6661264
                  Body_size_mm    297 5.339422e-01 4.877409e-01   0.4036511
                   Body_mass_g     68 1.583490e+00 1.411466e+00   1.4389471
               Longevity_max_y     36 8.090344e-01 8.061949e-01   0.8061949
         Age_at_maturity_min_y     34 6.356585e-01 6.435208e-01   0.6435208
             Litter_size_min_n    172 2.024907e+00 1.513213e+00   1.6322339
                lightness_male   1743 3.736833e+01 2.698113e+01  27.5747203
              lightness_female   1743 3.292126e+01 2.443325e+01  24.7389846
                      height_m    161 1.186425e+01 1.094375e+01          NA
                     leaf_area    215 1.360660e+04 1.416518e+04          NA
                           sla    287 2.096544e+01 2.162134e+01          NA
                     seed_mass    548 2.481197e+04 2.477287e+04          NA
                  wood_density    401 1.764024e-01 1.596984e-01          NA
                          WS_L    169 4.199106e-01 1.098606e+00   0.3151293
                          FW_L     19 3.640642e-01 2.216503e-01   0.2276826
                FlightDuration    450 3.673510e+00 4.139918e+00   4.0947548
     NumberOfHostplantFamilies    287 4.627133e-01 4.606028e-01   0.4712760
 cov_off_RMSE  ratio_on ratio_off     none_r   cov_on_r  cov_off_r
    2.0163097 1.0418163 1.4839581 0.89395480 0.88398144  0.7873380
    0.3732341 1.0907979 0.9371046 0.89077888 0.87152094  0.9035544
    0.6142631 0.7778215 0.7332916 0.80160123 0.77237446  0.7845749
    1.1190651 1.0000000 0.5986294         NA         NA  0.7986510
    2.3709944 0.9621711 1.0635114 0.76810459 0.78505374  0.7545352
    4.8882001 1.2835741 0.9217828 0.66955640         NA  0.7251741
    4.7143556 0.9847332 0.9949114 0.80109757 0.80792431  0.7995171
    0.5474663 0.8275934 1.1224532 0.62444270 0.70113755  0.5706778
    1.2466121 1.0194699 0.8832038 0.58837878 0.58932304  0.6675585
    0.8479460 1.0000000 1.0517878         NA         NA  0.4087068
    0.7520112 1.0000000 1.1685889         NA         NA -0.1253187
    1.6868349 1.0786546 1.1147374 0.66673533 0.59238467  0.6219735
   27.3986815 1.0220004 1.0154758 0.72455557 0.71289622  0.7167262
   24.5401953 1.0125129 1.0043769 0.69529008 0.68625137  0.6936129
           NA 1.1706322 5.5976854 0.46354321         NA         NA
           NA 1.0138293 1.0437579 0.14200526         NA         NA
           NA 1.0279007 1.7198546 0.33171931         NA         NA
           NA 1.0000131 0.9999957 0.05809927         NA         NA
           NA 0.9882790 1.0525132 0.44699975         NA         NA
    0.3290526 0.2868447 0.2995183 0.29302459 0.72689794  0.7126155
    0.2574245 1.0272154 1.1613993 0.80233925 0.78265440  0.8266274
    3.8413608 0.9890909 0.9278834 0.13319404 0.08177033  0.2603063
    0.6025781 1.0231724 1.3082382 0.22133989 0.14880777  0.1619579
```
