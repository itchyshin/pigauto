# CTmax follow-up results

The two fixed three-seed averages improved aggregate point prediction in this exercise, but did not improve predictions for the singleton species whose CTmax was entirely withheld. Seed sensitivity remains a material limitation. No package default or runtime code changed.

## Regime and provenance

Model source: `033649b93115c6d8191e791f06b1ddf52eb82d19`. The runtime source, metadata, bundled data and extdata matched a 61-file SHA-256 manifest. Runtime: Totoro CPU, R 4.6.1, torch 0.17.0, two torch threads and one inter-op thread per GNN process, at most six concurrent fits.

Public inputs are the CTmax files linked in [issue #260](https://github.com/itchyshin/pigauto/issues/260), pinned to BACE commit `b4b3cde02c4b171947042100e26b13b9f7ac425d`. Experimental rows: 2,661 observations, 524 species. Raw input MD5 values: data `316183777c6e2db74f98f4f3056d4cb7`, tree `ee09cab017022d3da54d88cecb886854`.

Five new row folds were stratified within species using split seed 2601109; fold sizes were 536, 508, 544, 528 and 545. CTmax was hidden before preprocessing. The tree and six reported covariates remained available. All observations were scored once per method. There were 237 singleton species; repeated-observation species retained target values in other outer folds. This is a different regime from the reporter's 375 standard-condition targets and is not an external replication or uniform whole-species validation.

The complete roster was five baseline fits with covariates supplied and 30 GNN fits, using seeds 51001 through 51006 on the same outer cells. Internal validation splits changed with fit seed, so the study measures end-to-end seed sensitivity. Two thousand epochs were requested: 25 GNN fits reached 2,000, four stopped at 1,900 and one at 1,800 through ordinary early stopping. All 35 attempts succeeded and yielded finite point predictions and interval endpoints. GNN fit times ranged from 225 to 260 seconds, median 250 seconds. The launcher finished within its 40-minute budget.

Five captured warnings were the same expected advisory that the GNN is off and baseline covariate use depends on the selected route. There were no other captured warnings. The earlier pilot was timing evidence only; it was not mixed into the final result set. Repeating its fold/seed in production gave exactly the same point predictions.

## Predeclared point summaries

RMSE is on the original degrees C scale. Species-balanced RMSE averages squared errors within each species, gives each species equal weight, then takes the square root. The individual-seed row is a range across six complete cross-validation runs, not an uncertainty interval.

| Method | Observation-weighted RMSE | Species-balanced RMSE |
|---|---:|---:|
| Baseline, covariates supplied | 1.993 | 2.204 |
| Individual GNN seeds | 1.632–1.706 | 2.153–2.364 |
| Mean of seeds 51001–51003 | 1.556 | 2.072 |
| Mean of seeds 51004–51006 | 1.526 | 2.032 |

All six individual GNN runs improved observation-weighted RMSE, but five were worse than baseline under species-balanced RMSE. Both predeclared trio means improved both summaries. This does not establish that every three-seed set will improve prediction or justify a package-wide recommendation.

After excluding the 110 eastern-newt observations, baseline observation-weighted RMSE was 1.985, individual GNN values were 1.648–1.721, and the trio means were 1.570 and 1.539. The newt therefore does not account for the entire observed improvement in this exercise. This sensitivity analysis removes those rows from scoring, not from fitting.

## Seed spread and individual intervals

Across six GNN seeds, the per-row prediction range had median 0.987 degrees C, 90th percentile 3.363 degrees C and maximum 14.007 degrees C. For 343 of 2,661 rows, the range exceeded 3 degrees C; excluding the newt left 342 of 2,551. These six-seed ranges should not be compared directly with the reporter's three-seed ranges under a different holdout design.

Individual-fit intervals contained 2,576/2,661 values for baseline and 2,518–2,551/2,661 across GNN seeds. Mean widths were 9.074 degrees C for baseline and 6.930–7.770 for the individual GNN runs. These are descriptive cross-validation summaries with dependent rows. No general calibration, conditional coverage or ensemble interval guarantee follows. No ensemble intervals or downstream-inference datasets were constructed.

## Post-hoc breakdown by remaining target observations

This analysis was added after the primary results were examined. No model was refitted or tuned in response.

| Method | No outer-training CTmax: 237 rows / 237 singleton species | Outer-training CTmax retained: 2,424 rows / 287 species |
|---|---:|---:|
| Baseline | 2.289 | 1.962 |
| Individual GNN seeds | 2.399–2.830 | 1.536–1.582 |
| Mean of seeds 51001–51003 | 2.406 | 1.447 |
| Mean of seeds 51004–51006 | 2.375 | 1.415 |

The table gives observation-weighted RMSE in degrees C. Every individual run and both trio means were worse than baseline for these singleton species, and better among species retaining CTmax. Singleton species still occurred in the supplied tree and had available covariates. The breakdown is descriptive and does not identify why performance differs. It argues against presenting seed averaging as a general remedy, especially where an entire species' target trait is missing.

## Separate factor-covariate check for #259

Using the reporter's original five folds and the GNN off, factor and numeric treatment coding agreed on lambda and predictions. Lambda ranged from 0.826858 to 0.918489. The maximum absolute lambda difference was below 1.1e-12; completed-value differences were at most 2.57e-10 degrees C, and uncertainty differences at most 1.64e-12. This supports coding invariance for that real-data workflow. It is not a lambda-recovery result or a check of every covariate design.

## Independent verification and retained artifacts

A separate reviewer recomputed the complete 35-fit roster, row/truth alignment, prepared/input/runner identities, all point summaries, seed spread and interval counts from raw records. Agreement with the summaries was within 1e-12. The post-hoc grouping and all its metric rows were also independently recomputed to that tolerance.

Raw inputs, per-row predictions, logs and results remain on Totoro at `~/sim-archive/pigauto/seed260-20261108/`. The summary directory contains `attempts.csv`, `point_metrics.csv`, `seed_spread.csv`, `individual_intervals.csv`, `posthoc_strata.csv` and `provenance.rds`. The coding check has a separate result and warning/provenance record. No raw rows were copied into this repository.

Primary summary SHA-256 values:

- `attempts.csv`: `01c5b786773114574a9f4cc2320b364d1b644815b6c2a89019f894e9c1cad18b`
- `point_metrics.csv`: `4bd3d8a35c1ca02b9ce311888cb61129e2762f21fa8f021e249bd059986bb3cc`
- `seed_spread.csv`: `04929255e8cb58ca41ce8d087665fb406ccccb12e8acde5685c716f3b725c859`
- `individual_intervals.csv`: `30f02f38e55a39a49d5569cbdff455fc0166748befee493364dc32aa4892bc57`

## Remaining work

The evidence supports a bounded result for this dataset and row-holdout design. A remedy for species with no observed target values needs separate development and validation. Do not select a subgroup fallback on these test results and then score it on the same rows as if it were independently validated. Broad mixed-type comparisons in #135 remain separate.

The no-refit helpers for #250 and #251 are ready and tested, but the actual mammal reproducer and stored draws have not arrived. The CRAN archive and its owner remain separate from this study; independent Win-builder review and the final release gate are still unproven here. No new issue replies, public results, tag changes or submission were made in this follow-up.

Review note: self-review of the full report on 8 November 2026, perceived AI-like style 3/10 with moderate confidence. The main contrast between observation and species weighting and the visible post-hoc label prevent an overbroad conclusion; no repair needed. Science and factual gates pass within independently checked artifacts; source references point to the originating issue and pinned data. Formal exact-draft assessment compiler: NOTASSESSED.
