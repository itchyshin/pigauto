# Rubin study: pigauto posterior MI arm (pig_post), pre-run results

Lane `claude:pigauto-mi-posterior` · branch `arc/rubin-freq-bace` · 2026-10-01 · Claude Code (Opus 5.5)

## What was run

- Arm `pig_post`: `pigauto::multi_impute(draws_method = "posterior", m = 20, log_transform = FALSE)` on the
  same continuous block as freqA (c1, c2, qlogis(prp), driver d1; counts excluded), written back with
  `apply_block_draw()` (`script/rubin_pigauto.R`). Opt-in in `script/rubin_cell.R`, seed offset 909.
- pigauto build: Main at b565cad (#189 merged), installed from GitHub into a private Totoro library, so each fit
  records `RemoteSha`.
- No conformal arm: Conformal draws already failed the downstream test (`arc/mi-gls-attenuation`, 16 regimes:
  slope bias -0.20 to -0.46, coverage 0 to 17 percent). Shinichi, 2026-10-01.
- Pre-run: n {100, 300, 1000} x lambda {0.3, 0.7, 1} x rho {0, 0.5} x seeds {1, 2} = 36 fits, MCAR 30 percent,
  Totoro, 36 cores, 1 h cap per cell. Started 19:22, finished 19:53. Same datasets as the stored campaign.
- Mac smoke: n = 100, lambda 0.7, rho 0.5, seed 1, all six arms.

## Gates (`.unlazy/pigauto-rubin-arm/GATES.md`)

| Gate | Result |
|---|---|
| GB1 smoke: all 6 arms finite slope and cor | `ARMS_OK 6` |
| GB2 smoke dataset equals the stored campaign dataset | `IDENTITY_OK` (numeric traits within 1e-12; mask and discrete traits identical) |
| GB4 freqA, freqB and complete data reproduce the stored rows | `REPRO_OK` (relative 1e-6) |
| GB3 provenance: sha b565cad, proper posterior MI, no logged trait | `PROV_OK`, 36 of 36 |
| GB5 pre-run fits that succeeded | `PRERUN 36` |
| GB6 fresh review | Rose: three blocking findings, all fixed (below) |

Two gate changes, both disclosed. First, GB2 and GB4 compared with `identical()`, which fails across platforms:
Mac against fir differs by at most 1.7e-15 in the data and 2.4e-9 (relative) in freqA's estimates. Tolerances
now sit at 1e-12 and 1e-6. A perturbation of 1e-6 in one data value, or 1e-4 in one estimate, still fails both.
Second, GB3 first required convergence as well. Convergence is an outcome of the run, not provenance, so it is now
reported beside the gate and flagged on every estimand row (`converged` column).

## Results

**Time per fit (single core, Totoro).**

| n | median | max |
|---|---|---|
| 100 | 327 s | 1,045 s |
| 300 | 616 s | 638 s |
| 1000 | 1,771 s | 1,812 s |

Convergence: 32 of 36 fits met pigauto's rule (split R-hat < 1.05, bulk ESS > 400). All 24 fits at n = 300
and 1000 converged with no extension. The 4 failures are all n = 100, seed 1 (lambda 0.3 and 0.7, both rho), and
each ran 3 extensions (R-hat 1.04 to 1.13, ESS 21 to 106). On the Mac smoke dataset (also seed 1), the poorly
mixing parameters are those of the fully observed driver d1 (Sigma_P[2,4], lambda[4], Sigma_P[4,4]). The same
dataset without d1 converged in 151 s with no extension. So the problem is one tree and the d1 trait, not n = 100
in general.

**Downstream estimands against the other arms, same 36 datasets** (mean over 12 datasets per n; d = estimate
minus the complete-data estimate; SE ratio = MI SE / complete-data SE; 2 seeds per cell, so read as a smoke, not
as evidence).

| n | arm | slope d | slope SE ratio | slope covered | cor d | cor covered |
|---|---|---|---|---|---|---|
| 100 | pig_post | -0.040 | 1.30 | 1.00 | -0.037 | 1.00 |
| 100 | freqA | -0.029 | 1.40 | 1.00 | -0.030 | 1.00 |
| 100 | bace_chain | -0.015 | 1.54 | 1.00 | -0.015 | 1.00 |
| 300 | pig_post | -0.053 | 1.45 | 1.00 | -0.047 | 0.92 |
| 300 | freqA | -0.063 | 1.45 | 0.92 | -0.054 | 0.75 |
| 300 | bace_chain | -0.085 | 1.46 | 0.92 | -0.068 | 0.92 |
| 1000 | pig_post | -0.028 | 1.44 | 0.92 | -0.031 | 0.92 |
| 1000 | freqA | -0.015 | 1.30 | 1.00 | -0.018 | 1.00 |
| 1000 | bace_chain | -0.021 | 1.37 | 0.92 | -0.023 | 1.00 |

No arm separates from the others at this size. Nothing here supports a claim in either direction.

## Campaign proposal (not launched; needs Shinichi's approval, D-287)

- Size: 18 cells x 200 seeds = 3,600 pig_post fits, paired with the stored freqA, freqB and BACE rows (BACE
  has 100 seeds at n = 1000).
- Cost from the measured times: About 183 core-hours at n = 100 (using the observed share of slow fits),
  205 at n = 300 and 590 at n = 1000: about 980 core-hours, or 1,225 with a 25 percent margin.
- Target: Totoro at up to 110 cores (the 150-core cap less the other lanes), about 11 to 12 h wall, no queue.
  DRAC (nibi or fir arrays) is the fallback if Totoro is busy. Read the first finished n = 1000 file within the
  first hour; stop and re-report if fits run more than 30 percent over these times.
- Open choice on convergence: (a) Keep the matched block with d1 and flag unconverged rows (expect a few
  percent at n = 100); (b) raise `max_extend` for n = 100 only, which adds cost only to the slow fits there (not yet measured); (c) drop d1 from pig_post's
  block, which converges faster but gives pigauto less information than freqA. Recommendation: (b), which keeps
  the comparison matched and costs little.

## Does NOT cover

Discrete traits (posterior MI is continuous-only); n = 1000 discrete (on hold, D-291); MAR or clade missingness;
real trees; the cause of the poor d1 mixing inside the sampler (a pigauto question for later).
