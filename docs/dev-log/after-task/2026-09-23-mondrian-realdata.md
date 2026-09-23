# After-task: Mondrian conformal real-data confirmation (2026-09-23)

Branch `arc/mondrian-realdata`, worktree `pigauto-mondrian-realdata`. Platform: Claude
(Opus 5.5 session). Status: DRAFT until the FishBase cell, the decision rule and the
review panel are in; sections marked PENDING are filled at close.

## 1. Goal

Decide whether `conformal_method = "mondrian"` should become the default, by testing it
against split conformal on real trait databases under a pre-registered rule; write the
paper's uncertainty section; and answer whether the Mondrian half-width can serve as the
draw SD in multiple imputation.

## 2. Implemented

- Stage-C harness brought onto main's line and extended: two mask arms (random, and a
  structured arm masking by a phylo-eigenvector propensity fitted to the real
  missingness), per-cell stratum labels, paired per-stratum summary with Winkler score
  and MCSE, AVONET input branch, smoke gate, results-doc generator with a byte-compare
  check, and the pre-registered decision-rule script.
- `compute_conformal_scores()` records `n_val`, `n_near`, `n_far` and names the realised
  stratum sizes on fallback.
- `mondrian_cell_scores()` extracted from `predict()`; `multi_impute()` conformal draws
  use per-cell Mondrian scores when the fit is Mondrian.
- Paper section 8 (uncertainty quantification) with five references checked; one
  citation corrected (Boström and Johansson 2020, PMLR 128).
- Pre-registration plus two amendments, each committed before the data it governs were
  read.
- MI draw-scale memo and a 500-replicate simulation; a one-tree diagnostic localising
  the GLS attenuation.
- PENDING: results table with FishBase; decision; NEWS entry.

## 3a. Decisions and Rejected Alternatives

- Two mask arms (Gauss review): a random mask of observed cells cannot show the failure
  Mondrian repairs; kept only as a no-harm control.
- Width cap read on the near stratum only (Shinichi, before any result).
- Retired the "within 3 MCSE of 0.95" rule as powerless at n_far of 19 to 60.
- Amendment 1: AVONET has almost no real missingness, so it runs the random arm only.
- Amendment 2: the structured-arm condition uses traits with at least 5% real
  missingness.
- FishBase moved from Tamia to kohaku (Shinichi), after a user-space R and CUDA 12.8
  install proved the GPU route.
- MI: did not wire a fix into `multi_impute()`; the attenuation is pre-existing and was
  spun into its own lane (`arc/mi-gls-attenuation`).

## 4. Files Touched

PENDING (generated from `git diff --stat origin/main...arc/mondrian-realdata` at close).

## 5. Checks Run

- `devtools::test(filter = "mondrian")`: FAIL 0, PASS 43.
- Full suite on 9f8f2b8: FAIL 0, PASS 2509, SKIP 8 (23.4 min).
- Smoke gate: SMOKE_OK for both arms. Decision-rule self-test: SELFTEST_OK.
- Results CSV regenerated from receipts: ROWS_MATCH.
- PENDING: gate-check --reverify on the full ledger; review panel.

## 6. Tests of the Tests

- G2b negative control: the stratum-size fields are absent on origin/main
  (LACKS_FIELDS) and present on the branch (FIELDS_OK).
- Decision-rule self-test exercises both verdicts on synthetic fixtures.
- PENDING: one coverage and one MCSE re-derived by hand from a receipt (reviewer M2).

## 7a. Issue Ledger

- MI draws halve a phylogenetic GLS slope (pre-existing; default split path): own lane.
- Default `predict_method = "per_column"` point prediction ignores co-observed traits.
- fir: `/project` quota full; cu126 libtorch needs `cuda/12.6` plus the runtime path.
- Mask separation warnings in the structured propensity fit on near-complete traits.

## 8. Consistency Audit

- The Boström citation error was fixed in both the paper section and the source
  paragraph file.
- The oracle scale error in the MI memo was fixed in the memo and in the diagnostic
  script.
- Swept R/, man/, vignettes/, NEWS.md and README.md for Mondrian wording that implies a
  coverage guarantee: none found. The roxygen already says the conformal options
  "support nominal held-out diagnostics, not package-certified coverage" and calls the
  /1.96 draw conversion a heuristic. The NEWS entry for the Mondrian MI draw change must
  carry the measured phylogenetic-GLS caveat.

## 9. What Did Not Go Smoothly

- Memory requests on fir guessed from AVONET scale: PanTHERIA cells and the MI array ran
  out of memory once; resized from measured peaks.
- The first MI pre-run launched with `R --vanilla`, which dropped the torch library;
  caught at one minute by the smoke check.
- Totoro was at its 150-core cap because of another lane; the MI work moved to fir.
- FishBase ran single-threaded for 2 h in its dense phase before being restarted with
  8 BLAS threads.
- My own MI oracle used the wrong trait scale; corrected and recorded.

## 10. Known Residuals

- Coverage on masked observed cells is a proxy for coverage on cells users impute.
- FishBase has one mask and is descriptive only.
- PENDING: anything the review panel raises.

## 11. Team Learning

- Size Slurm `--mem` from a measured peak on the target dataset, never from a smaller one.
- On fir, cu126 libtorch needs `module load cuda/12.6` and
  `LD_LIBRARY_PATH=$EBROOTCUDA/targets/x86_64-linux/lib`, even on CPU nodes.
- A pool-and-GLS downstream check belongs in any MI validation; marginal draw
  calibration and OLS can both look fine while GLS is badly biased.

## 12. Cross-Product Coverage

Covers: single-observation continuous-family traits, gnn on, three databases, two mask
arms. Does NOT cover: multi-observation data (Mondrian stops there by design), discrete
traits (no conformal intervals), `gnn = FALSE` (Mondrian stops there), MNAR beyond what
the structured arm captures, trees other than the three used.
