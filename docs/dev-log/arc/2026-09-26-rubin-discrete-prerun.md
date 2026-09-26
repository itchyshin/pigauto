# Discrete traits in the Rubin study: pre-run plan (2026-09-26)

Status: smoke done on the Mac; campaign NOT launched. It needs Shinichi's approval (over 3 h, D-287).

**What will be built, in one sentence.** The same 18 cells, seeds and datasets as the continuous campaign, re-run
with BACE (as shipped, chained, and the residual control) and the frequentist arms, saving all 20 imputations,
and scored on `bin`, `ord` and `cat3` per value and on a Rubin-pooled PGLS slope of `c1` on `bin`; n = 100 and
300 first.

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

**Q1. Frequentist discrete arm.** Recommend castor Mk (equal rates for `bin` and `cat3`, stepwise SUEDE for
`ord`, as in v1), in two forms that mirror the continuous arms:
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
The slope's truth is not ρ. It is the mean complete-data estimate over 2,000 fresh datasets per cell (seeds
1e6 + 1 to 2,000; Monte Carlo SE at most 0.011):

| n | λ | truth, ρ = 0 | truth, ρ = 0.5 | complete-data SD, ρ = 0.5 | datasets with one-class `bin` |
|---|---|---|---|---|---|
| 100 | 0.3 | 0.005 | 0.717 | 0.18 | 0% |
| 300 | 0.3 | -0.001 | 0.707 | 0.11 | 0% |
| 100 | 0.7 | 0.003 | 0.525 | 0.19 | 0.2% |
| 300 | 0.7 | 0.001 | 0.505 | 0.11 | 0% |
| 100 | 1 | 0.006 | 0.250 | 0.49 | 6.6% |
| 300 | 1 | 0.001 | 0.172 | 0.43 | 5.1% |

The ρ = 0 truths are 0 within Monte Carlo error, as they must be. At ρ = 0.5 the target falls with λ, and at
λ = 1 it also moves with n and the estimate is noisy: `bin` is nearly constant within clades, so the PGLS
contrasts carry little information about it. The estimand is well defined but belongs to this estimator,
not to the population. Accepting that is the recommendation; the alternative (a liability-scale correlation
from a threshold model) has a clean truth of ρ but has no fast frequentist estimator.

Handling rules (apply to every arm alike):
- a trait with one observed class: every arm imputes that class; the filled cells are counted;
- a trait with one class in the complete data: nothing to predict; skipped and counted;
- `c1 ~ bin` is undefined when `bin` has one class in the complete data; left out and counted.

**Q3. Scale.** Recommend n = 100 and 300 first: all 12 cells, 200 datasets each, seeds 1 to 200 (the continuous
campaign's datasets). Decide n = 1000 afterwards, from these results.

## Compute estimate

Measured per dataset (medians from the continuous campaign, confirmed by the Mac smoke):

| part | n = 100 | n = 300 | total for n 100 + 300 |
|---|---|---|---|
| BACE fit + chain (1,200 datasets per n) | 0.64 h | 1.65 h | about 2,750 core-hours |
| freq A + B + castor (1,200 per n) | 25 s | 75 s | about 35 core-hours |
| Monte Carlo truth | seconds | seconds | under 1 core-hour |

Wall clock: BACE is the whole cost. Split as before: n = 100 BACE on Totoro (1,200 × 0.64 h over up to 150
cores: about 5 to 6 h, if no other lane is using the 150-core allowance) and n = 300 BACE as SLURM arrays on
fir, nibi and rorqual (1,200 × 1.65 h; queue waits of 2 h and more were seen on nibi). Expected finish: 12 to
24 h after launch. The frequentist runs and the truth go on Totoro in under an hour. For n = 1000 later: about
3,400 core-hours of BACE, 5 to 7 h per dataset, 19 GB each.

A run that overshoots these numbers by more than half stops and re-reports.

## Before launch (gates)

1. Shinichi approves Q1 to Q3 and the scale.
2. Cluster checkouts at the committed hash of this plan; no edits during the run.
3. Launch with `--save_imp --discrete`, results to new folders (`results_disc/`), never over the continuous ones.
4. After the first 20 fits: check wall time, memory and errors, and whether BACE on nibi reproduces the stored
   continuous rows exactly. If it does, the re-run confirms and extends the continuous results; if not, it is
   reported as an independent replicate of them. Then release the rest.

## What this does not cover

- n = 1000 (decided later).
- Castor sees one trait at a time; a frequentist method that uses the correlated traits is not tested.
- One downstream estimand involving a discrete trait (`bin`); none for `ord` or `cat3`.
- The same DGP limits as the continuous study (Brownian motion with Pagel's λ, MCAR 30%, random trees).
