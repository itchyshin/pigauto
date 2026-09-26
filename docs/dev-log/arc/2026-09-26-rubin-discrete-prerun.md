# Discrete traits in the Rubin study: pre-run plan (2026-09-26)

Status (updated 2026-09-26, 09:15): Shinichi approved the plan ("go ahead", recommendations accepted). An
adversarial review then changed the frequentist rate models and the downstream target, and fixed bugs (section
"After approval: review and changes" below). Campaign LAUNCHED 09:11 at commit 55a0998 (Shinichi: "let's go ahead
and do the simulation"): nibi job 22721311 (BACE n = 100), fir 61677143 and rorqual 21855230 (BACE n = 300), Totoro
freq + castor (process group 1061449, 40 cores). The Totoro files record git_hash "HEAD" (their folder sits in a
git repository without commits); they ran 55a0998 (script/DISC_COMMIT, checksums verified at launch).

**What will be built, in one sentence.** The same 18 cells, seeds and datasets as the continuous campaign, re-run
with BACE (as shipped, chained, and the residual control) and the frequentist arms, saving all 20 imputations,
and scored on `bin`, `ord` and `cat3` per value and on a Rubin-pooled PGLS slope of `c1` on `bin`; n = 100 and
300 first.

## After approval: review and changes

A review workflow (four reviewers, each finding challenged by a skeptic) returned 15 findings; 13 survived. None
touches the BACE runs. Changes, all in commits 4a4395a and 55a0998 (lane suite FAIL 0, PASS 459; discrete tests
165, aggregation tests 67):

1. **Blocker, fixed.** Boundary castor fits (a stepwise rate at exactly 0) give a defective rate matrix, and the
   eigenvector route to the transition probabilities crashed. It hit freq A in about 30% of λ = 1, n = 100
   datasets and silently dropped all three discrete traits for that arm, on a selected, harder subset. The
   matrix exponential (`Matrix::expm`) now computes them, and a castor failure is isolated to its trait and
   recorded.
2. **Rate models changed (Q1 revised).** Equal rates (v1's model) run to the rate bound at low λ and then impute
   uniform classes, below the observed-frequency floor. Measured at n = 100, ρ = 0.5, 15 datasets per λ, with
   improper draws:

   | λ | trait | equal rates (ER, SUEDE) | flexible rates (ARD, SRD) | draws from observed frequencies |
   |---|---|---|---|---|
   | 0.3 | `bin` | 0.50 | 0.60 | 0.60 |
   | 0.3 | `cat3` | 0.36 | 0.43 | 0.46 |
   | 0.7 | `bin` | 0.68 | 0.74 | 0.68 |
   | 0.7 | `ord` | 0.52 | 0.58 | 0.45 |
   | 1 | `bin` | 0.96 | 0.97 | 0.74 |

   (accuracy of the modal class; 0 fit failures in 396 fits). Flexible rates are the frequentist mirror of the
   continuous arm: their stationary distribution can match the observed class frequencies, as Pagel's λ can fall
   back to the sample mean. They are now the primary freq A and freq B models (ARD for `bin` and `cat3`, SRD for
   `ord`). v1's equal rates stay as two comparison arms (`freqA_er`, `freqB_er`, seeds 707 and 808) on the same
   continuous draws, scored on discrete traits only. Cost: nothing extra in BACE; about 30 s per n = 300 dataset.
3. **Bootstrap rates the data cannot have come from.** A bootstrap refit can put a rate at exactly 0, under which
   the observed tips are impossible; the draw then returned NA silently. Such a refit is now rejected and redrawn
   (counted), with a counted fallback to the fitted rates after 20 tries. Smoke: 3 rejections, 0 fallbacks.
4. **Two targets for `c1 ~ bin`.** At λ = 1, ρ = 0.5 the dataset-specific value varies more than the standard
   error, and even complete-data intervals cover the population target in only 31% (n = 300) to 56% (n = 100) of
   datasets. Each result now stores the per-dataset target, ρ times the GLS slope (at the true λ) of `bin`'s
   liability on `bin`; complete-data coverage of it is about 0.96 in every cell. Coverage is read against it,
   with the population target reported beside it. Checked on the campaign's own 2,339 datasets with a defined
   estimand (n 100 and 300, seeds 1 to 200): complete-data coverage of the per-dataset target 0.905 to 0.985 by
   cell (lowest at n = 100, λ = 0.3, ρ = 0, the small-n shortfall of the complete-data analysis itself); of the
   population target 0.557 (n = 100) and 0.303 (n = 300) at λ = 1, ρ = 0.5, and 0.905 to 0.973 elsewhere.
5. Smaller fixes: `c1 ~ bin` is "undefined" (not an error) when `bin` has one observed class; per-value scores and
   the estimand are recorded separately; traits with one realised class are skipped for every trait type; ties at
   the edge of the 95% set and in the ordinal class distance count fractionally; cells whose draws are all NA are
   dropped and counted; the rescore script is safe per file and clears stale errors; the campaign driver builds
   its task list in a temporary file; result files record the code commit on the rsynced hosts.

Two review points stand as reported limitations, not fixes: BACE returns some NA draws for `cat3` (21 to 28 of
600 in the smoke; scored on the remaining draws and counted), and mean ECE over 30 to 90 cells is biased at small
n, so the study tables do not show it.

## What the smoke measured

| check | result |
|---|---|
| datasets rebuild from (n, λ, ρ, seed) | 12 of 12 stored datasets (n 100, 300, 1000): mask and discrete traits identical, continuous columns within 2e-12 (BLAS rounding across hosts) |
| frequentist continuous rows reproduce | freq A and freq B, n = 100 λ 0.7 ρ 0.5 seed 1: all scores within 5e-8 of the stored rds |
| BACE continuous rows reproduce | not exactly on the Mac (R 4.6.0; stored fit R 4.5.0 on nibi): pooled slope 0.597 against 0.590 (as shipped), 0.554 against 0.495 (chained), differences within one SE. The MCMC path diverges across R builds, so a Mac re-run is a replicate. Exact reproduction on the stored builds is checked at launch (gate 4) |
| castor arms, cost per dataset | freq A analogue 1.0 s (n 100), 4.9 s (n 300), 5.3 s (n 1000); freq B analogue under 0.3 s |
| joint draw is exact | tip marginals of 6,000 draws match brute-force exact marginals within 0.011 (ER 2 and 3 states, SUEDE 4 states) |
| castor's own marginals | match for ER, but off by up to 0.47 for SUEDE (ordinal); castor's rerooting method assumes a reversible model with its own root handling. Simulation v1's ordinal probabilities were therefore approximate |
| degenerate traits | λ = 1, seed 11 (n 100 and 1000): `bin` had one observed class and was one class in the complete data; skipped and counted; `ord` and `cat3` with empty levels scored normally |
| BACE discrete scoring in the run | works (same dataset): accuracy `bin` 0.73, `ord` 0.87, `cat3` 0.70 for BACE as shipped against 0.80, 0.77, 0.78 for freq A; `c1 ~ bin` pooled 0.88 (BACE) and 0.72 (freq A) against complete data 0.76. One dataset: no conclusion |
| lane suite | FAIL 0, PASS 330 (234 before, 96 new) |
| rds size with imputations | about 95 KB (n 100), 190 KB (n 300), 1.5 MB (n 1000) per frequentist file |

One slip, recorded: the first BACE smoke run was lost at the final save because `rubin_cell.R` was edited while
Rscript was reading it (Rscript parses a script incrementally). Cluster jobs must run from a committed checkout
that is not edited during the run.

## Decisions for Shinichi (recommendation first)

**Q1. Frequentist discrete arm.** (Revised after approval: flexible rates are primary, equal rates a comparison;
see above.) Recommend castor Mk (equal rates for `bin` and `cat3`, stepwise SUEDE for `ord`, as in v1), in two
forms that mirror the continuous arms:
- freq B analogue (improper): one fit, then 20 joint draws of the missing tips at the fitted rates;
- freq A analogue (proper): for each of the 20 draws, simulate the trait at the fitted rates, keep the observed
  tips, refit, then one joint draw at the refitted rates.
The joint draw is forward filtering, backward sampling on the tree (exact; built and tested). The castor draws
fill the freq A and freq B datasets after their continuous draws, under their own seeds, so the continuous rows
do not change. castor ignores `c1`, `c2` and the driver `d1`, so for `c1 ~ bin` the frequentist imputation
treats the two traits as independent. That is the honest frequentist option: no off-the-shelf phylogenetic
frequentist method imputes discrete traits jointly with continuous ones. A phylogenetic probit for `bin` in
`glmmTMB` is not recommended now: it would be a third arm for one trait only and would break the mirror.

**Q2. Scores.** Recommend, from each arm's 20 draws (every arm on the same basis):
- per value: accuracy of the modal class (ties counted fractionally), multi-class Brier score, top-label
  calibration error (ECE), coverage and size of the smallest class set holding 95% of the draws, and mean class
  distance for `ord`;
- downstream: the Rubin-pooled PGLS slope of `c1` on `bin` (`bin` = "yes" coded 1), scored on coverage, bias and
  the paired difference from the complete-data estimate.
The slope's truth is not ρ. The population target is the mean complete-data estimate over 5,000 fresh datasets
per cell (seeds 1e6 + 1 to 5,000; Monte Carlo SE at most 0.007; committed as
`script/rubin_study/data/discrete/truth_slope_c1_bin.csv`; the first version of this table used 2,000):

| n | λ | truth, ρ = 0 | truth, ρ = 0.5 | complete-data SD, ρ = 0.5 | datasets with one-class `bin` |
|---|---|---|---|---|---|
| 100 | 0.3 | -0.002 | 0.720 | 0.18 | 0% |
| 300 | 0.3 | -0.001 | 0.706 | 0.10 | 0% |
| 100 | 0.7 | -0.005 | 0.526 | 0.19 | 0.3% |
| 300 | 0.7 | -0.002 | 0.506 | 0.11 | 0.0% |
| 100 | 1 | -0.005 | 0.253 | 0.49 | 7.1% |
| 300 | 1 | -0.003 | 0.179 | 0.44 | 5.3% |

(After the review, coverage is read against the per-dataset target; see point 4 above.)

The ρ = 0 truths are 0 within Monte Carlo error, as they must be. At ρ = 0.5 the target falls with λ, and at
λ = 1 it also moves with n and the estimate is noisy: `bin` is nearly constant within clades, so the PGLS
contrasts carry little information about it. The estimand is well defined but belongs to this estimator,
not to the population. Accepting that is the recommendation; the alternative (a liability-scale correlation
from a threshold model) has a clean truth of ρ but has no fast frequentist estimator.

Handling rules (apply to every arm alike):
- a trait with one observed class: every arm imputes that class; the filled cells are counted;
- a trait with one class in the complete data: nothing to predict; skipped and counted;
- `c1 ~ bin` is undefined when `bin` has one class in the complete data, or one observed class (every imputed
  dataset then has a constant `bin`); left out and counted.

**Q3. Scale.** Recommend n = 100 and 300 first: all 12 cells, 200 datasets each, seeds 1 to 200 (the continuous
campaign's datasets). Decide n = 1000 afterwards, from these results.

## Compute estimate

Measured per dataset (medians from the continuous campaign, confirmed by the Mac smoke):

| part | n = 100 | n = 300 | total for n 100 + 300 |
|---|---|---|---|
| BACE fit + chain (1,200 datasets per n) | 0.64 h | 1.65 h | about 2,750 core-hours |
| freq A + B + castor, flexible and equal rates (1,200 per n) | 31 s | 112 s | about 48 core-hours |
| Monte Carlo truth | seconds | seconds | under 1 core-hour |

Wall clock: BACE is the whole cost. Placement (revised): all BACE on the DRAC clusters, which carry the same BACE
build (0.1.0) and R (4.5.0) as every stored n = 100 and 300 fit, so the re-run can reproduce them: n = 100 on
nibi (600 tasks of 2 fits, 3 h limit, 4 GB; measured peak 2.4 GB), n = 300 seeds 1 to 100 on fir and 101 to 200 on
rorqual (600 tasks each, 4 h limit, 10 GB; measured peak 6 GB). Totoro's BACE is a different build (0.0.0.9000,
R 4.5.3), so it runs only the frequentist arms (40 cores, about 1.2 h). If the arrays start promptly: 3 to 5 h;
nibi queue waits of 2 h and more were seen before. For n = 1000 later: about 3,400 core-hours of BACE, 5 to 7 h
per dataset, 19 GB each.

A run that overshoots these numbers by more than half stops and re-reports.

## Before launch (gates)

1. Shinichi approves Q1 to Q3 and the scale. DONE (2026-09-26).
2. Cluster copies at the committed hash; no edits during the run. DONE: 55a0998 on nibi, fir, rorqual and Totoro
   (checksums match; `script/DISC_COMMIT` on each host).
3. Launch with `--save_imp --discrete`, results to new folders (`results_disc/`), never over the continuous ones.
   Dry runs checked: 600 tasks per cluster, the campaign's BACE settings, output in `results_disc/`.
4. Changed for speed (Shinichi: "time matters"): the full arrays go in at once instead of a 20-fit pilot. The
   first completions are checked for wall time, memory, errors, and whether BACE on nibi reproduces the stored
   continuous rows exactly; any problem cancels the arrays. This risks little: the BACE path already ran 3,000
   fits, and every result keeps its imputations, so any scoring fix is applied afterwards without re-running.

Launch commands (from the Mac; they attach to the existing connections):

```
ssh nibi 'R=~/projects/def-snakagaw/snakagaw/pigauto_rubin; cd $R && RUBIN_OUT=$R/results_disc CELL_FLAGS="--save_imp --discrete" TAG=disc ARMSET=bace RUNS=5 NITT=50000 N=100 SEEDS=1-200 BLOCK=2 TIME=03:00:00 MEM=4G THROTTLE=600 CONFIRM=yes bash script/rubin_campaign_nibi.sh'
ssh fir 'R=~/pigauto_rubin; cd $R && RUBIN_ROOT=$R RUBIN_ENV=~/pigauto_sim/env.sh RUBIN_OUT=$R/results_disc CELL_FLAGS="--save_imp --discrete" TAG=disc ARMSET=bace RUNS=5 NITT=50000 N=300 SEEDS=1-100 BLOCK=1 TIME=04:00:00 MEM=10G THROTTLE=600 CONFIRM=yes bash script/rubin_campaign_nibi.sh'
ssh rorqual 'R=~/projects/def-snakagaw/snakagaw/pigauto_rubin; cd $R && RUBIN_ENV=$R/env.sh RUBIN_OUT=$R/results_disc CELL_FLAGS="--save_imp --discrete" TAG=disc ARMSET=bace RUNS=5 NITT=50000 N=300 SEEDS=101-200 BLOCK=1 TIME=04:00:00 MEM=10G THROTTLE=600 CONFIRM=yes bash script/rubin_campaign_nibi.sh'
ssh totoro 'cd ~/pigauto_rubin && for n in 100 300; do for l in 0.3 0.7 1; do for r in 0 0.5; do for s in $(seq 1 200); do echo "--n $n --lambda $l --rho $r --M 20 --arms freqA,freqB --seed $s --save_imp --discrete --out $HOME/pigauto_rubin/results_disc/freq"; done; done; done; done > logs/disc_freq_tasks.txt && bash script/rubin_totoro.sh logs/disc_freq_tasks.txt 40 3'
```

After the runs: `bash script/rubin_disc_pool.sh`, then `Rscript script/rubin_discrete_aggregate.R`, then
`Rscript script/rubin_study/build_data.R`, then render `script/rubin_study/study.qmd`.

## What this does not cover

- n = 1000 (decided later).
- Castor sees one trait at a time; a frequentist method that uses the correlated traits is not tested.
- One downstream estimand involving a discrete trait (`bin`); none for `ord` or `cat3`.
- The same DGP limits as the continuous study (Brownian motion with Pagel's λ, MCAR 30%, random trees).
