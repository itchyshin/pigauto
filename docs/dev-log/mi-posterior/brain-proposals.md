# Proposed brain-vault additions (staged; NOT written to the vault)

The vault needs Shinichi's explicit approval before any write (pigauto `AGENTS.md`, brain-write
boundary). These are the durable lessons from the posterior-MI arc (2026-09-24 to 25), each with its
evidence in this repository. Approve, edit or drop each one.

## 1. A gate script must fail closed: non-zero exit on failure, and the pass token only on the pass line

`gate-check.mjs` counts a gate as met when the process exits 0 and the EXPECT token appears anywhere
in the combined output (substring match). `script/mi_realdata/03_acceptance.R` exited 0 on failure
and printed its token inside a header line ("does not gate REALDATA_COMPLETE"). The ledger recorded a
PASS on a report with 19 failed conditions.
- Fix: commit 041c2e0.
- Evidence: `results.md` (gate section) and the after-task report, section 6.
- Suggested home: `CROSS-REPO-GUARDS` or `LESSONS`, with a mechanised check. For example, a lint that
  flags an EXPECT token appearing in any non-final `cat()` of a gate script.

## 2. Export BLAS/OpenMP thread caps in the launching shell, never only inside R

`Sys.setenv(OPENBLAS_NUM_THREADS = 1)` inside R runs after OpenBLAS has already started its threads.
A 75-process launch without shell-level caps pushed Totoro's load to about 18,000 for about 10
minutes (2026-09-24; `evidence/README.md`, "Incident"). The GATES CHECK lines were fixed the same
day to prefix `env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
- Suggested home: the compute-routing skill or `tools/totoro-setup.md`.

## 3. A Monte Carlo SE for a ratio of two SDs from the same replicates must account for their correlation

Treating the MI and complete-data empirical SDs as independent overstated the MCSE of their ratio by
about 1.8x (correlation 0.64 to 0.93). That turned three real G6 failures into apparent noise.
- Correct delta-method form: rel * sqrt((1 - rho^2) / (R - 1)), confirmed by a bootstrap.
- Evidence: `script/mi_gls/07_se_ratio_noise.R`, `se_ratio_noise.txt`, and the M2 review.
- Suggested home: the validation-harness skill (paired-comparison statistics).

## 4. In-model twins separate model-family mismatch from estimator error

When a simulation's DGP lies outside the fitted model family, project the same draws (same seeds,
trees and masks) onto the family and run both. In this arc the twin removed 54% to 86% of an apparent
estimator bias. The rest could then be attributed (`diagnosis.md`, `results.md`).
- Suggested home: `WHAT-WORKS`.

## 5. pigauto's `cov2cor(vcv(tree))` convention gives every tip equal variance

On non-ultrametric trees this is a real modelling choice for every imputation path, not a detail.
Simulations that draw from the raw `vcv(tree)` are misspecified for pigauto. For ultrametric (dated)
trees the two are proportional and the issue does not arise (`diagnosis.md`).
- Suggested home: `projects/pigauto.md` (a design invariant with a known consequence).

## 6. D-280 head-to-head verdict

Already recorded in the vault (DECISIONS.md D-280 and MODEL-ROUTING.md, checked 2026-09-24). No
action.
