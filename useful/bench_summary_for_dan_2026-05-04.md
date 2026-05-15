# pigauto bench summary for Dan (cross-dataset, 2026-05-04)

Single-seed default-config runs.  Sources cited inline.

## Trait selection per dataset

These are the traits pigauto was benched on; each is in
`phyloTraitData::<dataset>_traits` so BACE can target the same set.

| Dataset           | n       | Continuous traits                               | Discrete traits                                    |
|-------------------|--------:|--------------------------------------------------|----------------------------------------------------|
| AVONET (birds)    | 9,993   | Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length | Trophic.Level (cat K=4), Primary.Lifestyle (cat K=7), Migration (ord K=3) |
| PanTHERIA (mam.)  | 4,027   | body_mass_g, head_body_length_mm, gestation_d, max_longevity_m | litter_size (count), diet_breadth (ord K=5), habitat_breadth (ord K=3), terrestriality (bin) |
| AmphiBIO (amph.)  | 5,237   | Body_size_mm, Body_mass_g                       | (continuous-only in v1 bench; discrete cols added in phyloTraitData) |
| BIEN (plants)     | 4,745   | height_m, leaf_area, sla, seed_mass, wood_density | -- |
| FishBase (fish)   | 10,654  | Length, Weight, DepthRangeDeep, Vulnerability, Troph | BodyShapeI (cat) |

n = species (post tree alignment).  miss_frac = 0.30 MCAR per trait, seed = 2026.

## Headline per-trait metrics

`pigauto_default` vs `mean_baseline`.  Continuous: RMSE on raw
(back-transformed) scale + Pearson r against truth.  Discrete: accuracy.
Coverage = empirical 95% conformal coverage.

### AVONET birds (n = 9,993; `script/validate_avonet_full.md`)

Numbers below are on the latent z-score scale (test-cell metrics; full-scale validation run).

| Trait              | Type        | n_test | BM RMSE | GNN RMSE | r     |
|--------------------|-------------|-------:|--------:|---------:|------:|
| Mass               | continuous  | 1058   | 0.151   | 0.151    | 0.989 |
| Beak.Length_Culmen | continuous  | 1066   | 0.236   | 0.236    | 0.970 |
| Tarsus.Length      | continuous  | 1048   | 0.198   | 0.198    | 0.981 |
| Wing.Length        | continuous  | 1089   | 0.169   | 0.169    | 0.986 |
| Trophic.Level      | categorical | 1069   | -- (acc) | -- (acc) | --    |
| Primary.Lifestyle  | categorical | 1059   | -- (acc) | -- (acc) | --    |
| Migration          | ordinal K=3 | 1100   | -- (acc) | -- (acc) | --    |

(Cell-level discrete accuracy in the md file; gate stays mostly closed on
high-phylo-signal traits so BM and GNN coincide.)

### PanTHERIA mammals (n = 4,027; `script/bench_pantheria_full.md`)

| Trait               | Type     | n_test | mean_baseline | pigauto_default | r      |
|---------------------|----------|-------:|---------------|-----------------|-------:|
| body_mass_g         | log-cont | 1042   | 3.077 RMSE    | **0.725 RMSE**  | 0.972  |
| head_body_length_mm | log-cont |  573   | 1.021 RMSE    | **0.249 RMSE**  | 0.970  |
| gestation_d         | log-cont |  407   | 0.961 RMSE    | **0.387 RMSE**  | 0.917  |
| max_longevity_m     | log-cont |  303   | 0.914 RMSE    | (see md)        | (md)   |
| litter_size         | count    |  742   | 1.766 RMSE    | (see md)        | (md)   |
| diet_breadth        | ord K=5  |  646   | 0.385 acc     | (see md)        | --     |
| habitat_breadth     | ord K=3  |  828   | 0.742 acc     | (see md)        | --     |
| terrestriality      | binary   |  802   | 0.551 acc     | (see md)        | --     |

### AmphiBIO amphibians (n = 5,237; `script/bench_amphibio.md`)

| Trait        | Type     | n_test | mean_baseline | pigauto_default | r     |
|--------------|----------|-------:|---------------|-----------------|------:|
| Body_size_mm | log-cont | 1569   |   89.08 RMSE  | **43.57 RMSE**  | 0.872 |
| Body_mass_g  | log-cont |  178   | 1947.74 RMSE  | **1422.31 RMSE**| 0.978 |

(Discrete columns -- diurnal/nocturnal binary, diet_breadth ordinal,
habitat categorical -- are in `phyloTraitData::amphibio_traits` but
were skipped in the v1 amphibio bench because the threshold-joint
baseline hits an Rphylopars singular-matrix error on the AmphiBIO
taxonomic tree.  Phase F LP corner / OVR categorical fix that.)

### BIEN plants (n = 4,745; `script/bench_bien.md`)

The kingdom jump.  Coverage holds; point estimates struggle on this
many-thousands-of-clades tree without environmental covariates.

| Trait        | Type     | n_test | mean_baseline | pigauto_default | r       |
|--------------|----------|-------:|---------------|-----------------|--------:|
| height_m     | log-cont |  218   |  11.68 RMSE   | 15.15 RMSE      |  0.197  |
| leaf_area    | log-cont |  297   | 12152 RMSE    | 24401 RMSE      | -0.024  |
| sla          | log-cont |  386   |  19.99 RMSE   | 23.04 RMSE      |  0.205  |
| seed_mass    | log-cont |  743   | 1745 RMSE     | (see md)        | (md)    |
| wood_density | continuous | 560 |   0.177 RMSE  | (see md)        | (md)    |

Plants need the WorldClim covariate path (`pull_worldclim_per_species`)
to lift point estimates above the grand-mean baseline.

### FishBase fish (n = 10,654; `script/bench_fishbase.md`)

| Trait          | Type        | n_test | mean_baseline | pigauto_default | r      |
|----------------|-------------|-------:|---------------|-----------------|-------:|
| Length         | log-cont    | 3103   |  41.20 RMSE   | **25.47 RMSE**  | 0.809  |
| Weight         | log-cont    |  565   | 87111 RMSE    | **84804 RMSE**  | 0.418  |
| DepthRangeDeep | log-cont    | 1410   | 883.96 RMSE   | (see md)        | (md)   |
| Vulnerability  | continuous  | 3197   |  17.87 RMSE   | (see md)        | (md)   |
| Troph          | continuous  | 1457   |   0.645 RMSE  | (see md)        | (md)   |
| BodyShapeI     | categorical | 3183   |   0.464 acc   | **0.775 acc**   | --     |

## Files Dan can pull verbatim

- `script/validate_avonet_full.md` -- AVONET headline + wall time
- `script/bench_pantheria_full.md` -- PanTHERIA headline
- `script/bench_amphibio.md` -- AmphiBIO headline
- `script/bench_bien.md` -- BIEN headline
- `script/bench_fishbase.md` -- FishBase headline
- `script/bench_pantheria_bace_head_to_head.md` -- direct PanTHERIA pigauto-vs-BACE comparison (preview of what Dan's CI report will compare to)
- `script/bench_phase_c_cross_dataset.md` -- Phase C v2 multi-seed sweep over (default / clamp_outliers / pool_method=mode / both) on PanTHERIA + AmphiBIO

## Runtime ranges (single seed, default config, on Mac M-series)

| Dataset              | Wall time      |
|----------------------|----------------|
| AVONET    (n=9,993)  | ~9 min         |
| PanTHERIA (n=4,027)  | ~5 min         |
| AmphiBIO  (n=5,237)  | ~4 min         |
| BIEN      (n=4,745)  | ~6 min         |
| FishBase  (n=10,654) | ~12 min        |

GitHub Actions Ubuntu runner is ~2x slower; budget accordingly.
