# After-task report: posterior multiple imputation (`multi_impute(draws_method = "posterior")`)

2026-09-24 to 25. Branch `arc/mi-posterior`, draft PR itchyshin/pigauto#189 (base `main`). Platform:
Claude Code (Opus 5.5 orchestrator). Plan:
`~/.claude/plans/read-agents-md-and-docs-dev-log-handover-curious-eagle.md`.

## 1. Goal

Deliver `multi_impute(draws_method = "posterior")` for continuous traits: joint draws from
Sigma_P %x% R + Sigma_E %x% I with full Sigma_E, where Sigma, lambda and the means are redrawn by
Gibbs data augmentation (proper imputation) and the GNN is not used. It ships with tests, docs and
NEWS, a simulation report, a real-data report, and a draft PR. The headline has two parts, both
judged against complete data under the same analysis model:
- per-cell predictive intervals that cover the truth about 95% of the time;
- a pooled downstream slope that is unbiased with an honest SE.

## 2. Implemented

Sampler: `R/mi_posterior.R` implements the model above, following `design.md` sections 1 to 2.5:
- the Hadfield and Nakagawa sparse node precision, rescaled so the tip block is R (Qc = D Q D);
- a collapsed (a, mu) | y_obs block, then y_mis | a, mu;
- MCMCglmm-style parameter expansion, inverse-Wishart updates, and collapsed Metropolis moves for
  mixing near lambda = 0 or 1;
- internal split R-hat and bulk ESS;
- automatic chain extension: exact continuation from saved state and RNG, up to 4x;
- plug-in modes `"none"` and `"both"`, for validation.

Wiring: `R/multi_impute.R` adds `draws_method = "posterior"` and `posterior_control`, clear errors
for unsupported inputs, a print method, and roxygen and man pages.

Provenance: `R/with_imputations.R` and `R/pool_mi.R` accept the provenance marker
`pigauto_posterior_mi_v1`. Plug-in draws carry a diagnostic marker, and both functions refuse it.

Tip depths: `R/henderson_s_inv.R` gains an optional `tip_depths` argument, so the posterior path
builds no dense n x n matrix; the default path is unchanged.

Tests: `tests/testthat/test-mi-posterior.R` has 118 expectations and
`tests/testthat/test-multi-impute-posterior.R` has 196. They include:
- a kernel-agreement test, Metropolis-only vs Gibbs-only;
- exact continuation;
- a folded R-hat test;
- sampler row slicing;
- a first-run regression test pinned to commit 69670d4.

Harnesses: - simulation, in `script/mi_gls/`: 40 regimes, a Totoro runner, fail-closed gates G4, G6, G7 and G3
  (calibration), results tables and the SE-ratio noise analysis;
- real data, in `script/mi_realdata/`: 10 Mondrian-mask cells, a Totoro runner, and the fail-closed
  gate G8.

**Records**, in `docs/dev-log/mi-posterior/`:
- `design.md` (every change of plan, with who decided and when);
- `review-design.md`, `review.md`, `diagnosis.md`, `results.md`, `results_tables.md`,
  `se_ratio_noise.txt`;
- `evidence/` and `real_preview/`;
- `note-for-lambda-lane.md`;
- the plan-vs-actual reconciliation, `docs/dev-log/plan-actual/2026-09-24-mi-posterior.md`.

Campaigns on Totoro: - first campaign: 24 regimes x 200 reps at 69670d4;
- G3 calibration: 200 fits;
- diagnosis;
- round 2: 3,200 twin cells plus 44 reruns at 9597e18, merged to 8,000 cells;
- real data: 10 cells at 69670d4.

## 3a. Decisions and Rejected Alternatives

Blind head-to-head (D-280): Shinichi picked Plan 1 (Opus 5.5) over Plan 2 (Fable 5.1).

Locked by Shinichi before the build: - full Sigma_E;
- the GNN ignored;
- weak inverse-Wishart priors with parameter expansion;
- the name `"posterior"`.

Orchestrator decisions: - Draw (a, mu) with y_mis integrated out, rather than alternating a and y_mis. Alternation mixes badly
  when lambda is near 1.
- Provenance option B: a new marker. Rejected: reusing `pigauto_analysis_mi_v1` (it misstates the
  provenance) and bare lists (no provenance at all).

At CP1 (before any result): - G3 changed from a per-fit tolerance of 0.1 to a calibration check;
- the SE ratio is judged relative to complete data;
- proper-vs-plug-in is reported, not gated;
- the real-data 5% slope criterion is reported, not gated;
- the campaign moved to Totoro (Shinichi's compute routing);
- the real-data overrun was allowed to continue.

At CP2 (after the round-1 diagnosis): - regimes 1 to 16 re-run as in-model twins, kept as a stress test;
- automatic chain extension added.

Rejected: (`diagnosis.md`):
- switching the sampler to the raw covariance, which would break pigauto's `cov2cor` convention;
- a 100x smaller Sigma_E prior, which failed numerically in 2 of 24 fits;
- doubling the default chain length, which costs every user twice the wall time.

The G6 relative SE-ratio rule was not changed after the results: The three failures are
presented to Shinichi as a decision, with options (`results.md`).

## 4. Files Touched

Diff of `origin/main...arc/mi-posterior`.

Added: - `R/mi_posterior.R`
- `R/draws_conditional.R` (the prototype conditional draw, used by G2)
- `tests/testthat/test-mi-posterior.R`
- `tests/testthat/test-multi-impute-posterior.R`
- `tests/testthat/test-draws-conditional.R`
- simulation scripts in `script/mi_gls/`: `01_cell.R`, `01_cell_v2.R`, `02_summarise.R`,
  `03_combine.R`, `03_summarise_v2.R`, `04_acceptance.R`, `05_cell_coverage.R`,
  `06_results_tables.R`, `07_se_ratio_noise.R`, `11_fir_array.sbatch`, `12_totoro_campaign.sh`,
  `dgp_v2.R`, `regimes.R`
- gate scripts in `script/mi_gls/`: `gate_calibration.R`, `gate_convergence.R`,
  `gate_exactness.R`, `gate_recovery.R`
- `script/mi_gls/diag/` (diagnosis scripts)
- real-data scripts in `script/mi_realdata/`: `00_fetch_masks.R`, `01_run.R`, `02_summarise.R`,
  `03_acceptance.R`, `10_fir.sbatch`, `12_totoro_run.sh`, `lib.R`, `pairs.R`
- records in `docs/dev-log/mi-posterior/`: `design.md`, `review-design.md`, `review.md`,
  `diagnosis.md`, `results.md`, `results_tables.md`, `se_ratio_noise.txt`, `sim_summary.csv`,
  `sim_summary.md`, `cell_coverage.csv`, `note-for-lambda-lane.md`, `evidence/`, `real_preview/`
- `docs/dev-log/plan-actual/2026-09-24-mi-posterior.md`
- this report

Modified: - R: `NEWS.md`, `R/henderson_s_inv.R`, `R/multi_impute.R`, `R/multi_impute_trees.R` (seealso only),
  `R/pool_mi.R`, `R/predict_pigauto.R` (roxygen text only), `R/with_imputations.R`
- `README.md` (one passage on the supported inference routes)
- man pages: `man/multi_impute.Rd`, `man/multi_impute_trees.Rd`, `man/pool_mi.Rd`,
  `man/predict.pigauto_fit.Rd`, `man/with_imputations.Rd`

Not committed (local): - `.unlazy/mi-posterior/GATES.md` (the acceptance ledger; `.unlazy/` is in `info/exclude`)
- Totoro run folders under `/home/snakagaw/pigauto_mi_posterior/`

Not touched: `R/joint_mvn_solver.R` (G5b).

## 5. Checks Run

Ledger: (`gate-check.mjs --reverify`, worktree root, 1 h timeout): FINAL COUNTS AT THE END OF THIS
REPORT (section "Final gate state").

Per gate: (details in `results.md`):
- G1: two posterior test files, 118 + 196 expectations, 0 failures.
- G2: `EXACTNESS_OK`.
- G3: `RECOVERY_OK` (200 fits, 69670d4).
- G4: `CONVERGENCE_OK`.
- G5a: full suite 2,840 pass, 0 fail (9059f87).
- G5b: `SOLVER_UNTOUCHED`.
- G5c: R CMD check --as-cran, 0 errors and 0 warnings, 1 NOTE for the dev version (9059f87).
- G6: fails 3 in-model rows, relative SE ratio.
- G7: `CELL_COVERAGE_PASS`.
- G8: pending the FishBase cell.

CI on PR #189: R CMD check passed at 3cd3139 on macOS (R release) and Ubuntu (R release, R devel). That commit includes the last change to `R/` and `tests/` (8d2f612); later pushes were docs-only and their runs were superseded. pkgdown is skipped on PRs by design.

Campaign integrity: - `SETTINGS 0 of 8000` non-campaign files;
- `MIXED_CODE_SHA`, as expected;
- 16 comparisons of 15 distinct converged cells were byte-identical between 69670d4 and 9597e18.

## 6. Tests of the Tests

- Kernel-agreement test: caught 8 of 8 sampler mutants (dropped Jacobians, bare Q, the wrong
  degrees of freedom, row shifts); max |z| 5 to 25 against 1.8 on clean code.
- Continuation test: after the checker's repair, it catches a mutant that skips restoring the
  adapted step sizes.
- Folded R-hat test: fails when `.mip_rhat` computes the bulk half only.
- Row-slicing and interval-provenance tests: each caught its own mutants in scratch copies.
- Twin DGP self-check: catches a mask-RNG shift and rescaling by the wrong matrix.
- Gate selftests: they include fixtures that must fail, and do. These cover a missing regime, a
  short rep count, a non-finite value, more than 2% non-converged, and exactly 2% (which must pass).
  They also cover a stress-regime failure (must not gate), a gated twin failure, and the pairing of
  proper and plug-in rows.
- A real false PASS was found and fixed: G8's acceptance script exited 0 on failure and printed
  its token inside a header line, so gate-check recorded a pass on 19 failed conditions. Fixed in
  041c2e0. Every other gate script was audited for the same two faults; none had them.

## 7a. Issue Ledger

Fixed: - Design review: B1 (Qc scaling), B2 (parameter-expansion equations), B3 (convergence gating),
  R1 to R6.
- Review of the simulation and real-data scripts: 28 findings: fail-closed completeness, the coverage truth for regimes
  17 to 24, the conformal arm, and runner bugs.
- S5 review: 28 findings, including plug-in provenance, the K = 1 crash, the congeniality scope, and
  the kernel test.
- Round-2 review: mixed-format summarise, settings guard, package versions, fir array, pinned test.
- Claims audit: 18 documentation findings.
- M2 review:
  - the SE-ratio noise MCSE, which wrongly assumed independence; the corrected analysis shows the G6
    failures are not noise;
  - the lambda = 1 residual-bias caveat;
  - the under-coverage caveat;
  - the real-data preview source and numbers;
  - the rule-2 float tolerance.
- The G8 false pass.
- The thread-cap hole in the GATES check lines.

Deferred (for Shinichi): - the G6 relative-SE-ratio decision;
- the default draws method;
- a `log_transform = FALSE` sensitivity run for PanTHERIA;
- the Sigma_E prior's pull on the residual correlation at lambda near 1 (possible remedy: an
  off-diagonal-aware scale);
- whether tip-variance heterogeneity on non-ultrametric trees matters for users. It is a pigauto-wide
  convention, not specific to this method.

## 8. Consistency Audit

- Gate scripts: the token and exit-code audit covered every gate script after the G8 false pass.
  Thread caps were checked on every CHECK line and in both Totoro runners.
- Byte-identity: the new code reproduces converged 69670d4 cells exactly, shown by 16 comparisons
  and the pinned test.
- Documentation: the claims audit covered NEWS, roxygen, man, print, messages, README and
  `design.md` against the code, and found 18 inaccuracies, all fixed. `multi_impute_trees()`'s
  seealso was updated as a neighbour.
- The joint baseline is unchanged: `build_henderson_S_inv()`'s default path is untouched and
  `R/joint_mvn_solver.R` is untouched.
- The plug-in covariance shrinkage: found in the MI-GLS lane was handed to the lambda lane as a
  written note, `note-for-lambda-lane.md`.

## 9. What Did Not Go Smoothly

- Thread-cap incident: The first G3 calibration launch exported no thread caps. Setting them
  inside R is too late for BLAS, so the Totoro load reached about 18,000 for about 10 minutes. The
  run was killed and relaunched with caps.
- Campaign overruns: Both campaigns took longer than estimated because Totoro was contended (load
  300 to 500 from other users), and other lanes of this account launched uncapped jobs (drmTMB,
  `exact_prerun_cell.R`), taking the account above 150 cores. Shinichi chose to keep ours running.
  - The FishBase real-data cell ran for more than 15 h.
- G8 false pass: Found only because the gate "passed" with no receipts present.
- Rep-spec bug: The round-2 driver passed single rep numbers ("118", "37"), which the runner
  reads as counts. It re-ran 14 converged cells serially before this was caught. They proved
  byte-identical, so they became extra evidence.
- The first SE-ratio noise argument was wrong: (it assumed independence). The M2 review caught it.
- Minor slips: - `closeout.py` wrote the report template into the brain vault's folder, because it resolves
    relative paths against the vault. The file was moved out, so no vault content changed.
  - Two `rcmdcheck` launches failed on an empty `check_dir` variable.
  - A `pkill -f` pattern matched its own ssh shell.

## 10. Known Residuals

- G6: is not met on 3 in-model rows (relative SE ratio 1.153 to 1.186). The absolute MI SE ratio
  in those rows is 0.98 to 1.06; the complete-data analysis there is over-confident (0.85 to 0.90).
  Decision pending.
- G8: is pending the FishBase cell. FINAL STATE AT THE END OF THIS REPORT.
- **At lambda = 1**, all 16 twin rows show a small negative paired bias (-0.004 to -0.014, 2.7 to 8.7
  MCSE). They pass only through the 0.02 floor.
- Per-cell coverage: is slightly below 0.95 in regimes 17 to 24 and in the lambda = 0.5 twins
  (0.933 to 0.949).
- Evidence that predates the final code: - G3 calibration evidence is from 69670d4 (the sampler is identical on converged fits);
  - real data ran at 69670d4, where two PanTHERIA cells are ESS-short (366 and 368).
- Cost of large trees: FishBase-sized trees (about 10,000 tips, 5 traits) need hours per fit.
- Minor: the runner's ok counter reported 0 for one completed rerun.
- K > 2: is exercised only on real data (4 to 5 traits), not in the simulation.

## 11. Team Learning

- Proposed brain additions: (staged, not written; the vault needs Shinichi's approval):
  - a gate script must exit non-zero on failure and print its token only on the pass line, because
    `gate-check` matches substrings;
  - export BLAS and OpenMP thread caps in the launching shell, never only inside R;
  - the runner's rep spec treats a single number as a count; use `a-a` for one rep;
  - an MCSE for a ratio of SDs computed on the same replicates must account for their correlation.
- What worked: - adversarial review workflows with two skeptics per finding caught every real defect in this arc;
  - in-model twins (same seeds, with the DGP projected onto the model family) cleanly separated
    model-family mismatch from estimator error.

Memory receipt: loaded the pigauto `AGENTS.md` and `CLAUDE.md`, the hub `AGENTS.md` (D-139 estimate
before running, D-143 core cap, thread caps, the brain-write boundary, the Rose principle, D-43),
`protocols/after-task.md` and `protocols/handoff.md`. They shaped the CP1 and CP2 stop points, the
overrun re-reports and the thread-cap fix.

Golden Set: not in scope (no memory-regression class touched).

## 12. Cross-Product Coverage

Covers:
- continuous traits on pigauto's latent scale (log-transformed when all positive);
- single-observation data;
- ultrametric and non-ultrametric trees, with R = cov2cor(vcv(tree));
- K >= 1 traits (K = 1 tested);
- `with_imputations()` and `pool_mi()` for the fit classes `pool_mi()` supports (gls and lm used here; phylolm was pooled by hand with Rubin's rules in the harnesses);
- `param_uncertainty` `"full"` (default), `"none"` (refused downstream) and `"both"` (validation);
- automatic chain extension;
- print, errors and docs.

Does NOT cover:
- binary, categorical, ordinal, count, proportion, zi_count or multi_proportion traits (clear errors);
- multi-observation data (error);
- covariates in the imputation model;
- GNN blending;
- `multi_impute_trees()` (no posterior method across trees);
- `multi_impute_analysis()`;
- analyses with nonlinear terms, interactions or external covariates (the congeniality scope);
- a default change (still `"conformal"`);
- CRAN release checks beyond `--as-cran` on one Mac.

## Final gate state

Ledger (`.unlazy/mi-posterior/GATES.md`, local; re-verified 2026-09-25):

| Gate | State |
|---|---|
| G1 to G5c | met at the final code |
| G6 | not met: 3 in-model rows of the relative SE-ratio rule; decision owed by Shinichi |
| G7 | met |
| G8 | see the handover's Final state (depends on the FishBase cell) |
| M1 | met (design review) |
| M2 | pending: findings fixed; the D-43 panel withheld nothing; the formal PROCEED waits for the final real-data section |

`check-after-task.R` correctly refuses to call the work finished while G6, G8 and M2 are unmet.
