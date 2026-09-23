# After-task: Mondrian conformal real-data confirmation (2026-09-23)

Branch `arc/mondrian-realdata`, worktree `pigauto-mondrian-realdata`. Platform: Claude
(Opus 5.5 session). Status: complete; results, decision, NEWS and review panel are in.

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
- Pre-registration plus two amendments, each committed before any outcome it governs
  existed (exact receipt and commit times are in the dated clarifications under each
  amendment header).
- MI draw-scale memo and a 500-replicate simulation; a one-tree diagnostic localising
  the GLS attenuation.
- Results table from all 20 method receipts (PanTHERIA 12, AVONET 6, FishBase 2); rule
  verdict KEEP_SPLIT (condition 2, near-stratum non-inferiority, failed; conditions 1
  and 3 passed on every dataset); NEWS records the verdict and the MI GLS caveat;
  paper section 8.3 reports Table S-UQ2.
- Fallback message now names the realised stratum sizes (claim-gate blocking item).

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
- The decision-rule script was rewritten at `5709a7f`, after all results were read, to
  evaluate conditions 1 and 3 per dataset as the pre-registration says, to fail closed on
  condition 2 with no evidence, and to gate on mask completeness. The verdict is
  KEEP_SPLIT under both versions; conditions 1 and 3 pass and condition 2 fails under
  both readings.
- Applied the rule as registered rather than re-reading it after the data: its near
  non-inferiority condition penalises removal of over-coverage, which is Shinichi's call
  to revisit, not this arc's.

## 4. Files Touched

From `git diff --stat origin/main...HEAD` (39 files, excluding receipts):

- R/: `fit_helpers.R`, `multi_impute.R`, `predict_pigauto.R`.
- tests/: `test-mondrian-conformal.R`, `test-mondrian-mi-draws.R` (new).
- script/mondrian_confirmation/: `00` to `13b` (harness, launchers, generators, rule,
  diagnostics) and `returned/` (20 method receipts plus mask receipts).
- docs/dev-log/: the 08-16 and 08-18/19 records; `mondrian-realdata/` (pre-registration,
  recon, run log, kohaku install, results, MI summary and diagnostic log, M2 audit);
  `review/` (traceability, claim gate); this report.
- useful/: `paper_section_draft.md`, `mondrian-methods-paragraph.md`,
  `mondrian-mi-se-justification.md`.
- NEWS.md.

## 5. Checks Run

- `devtools::test(filter = "mondrian")`: FAIL 0, PASS 43.
- Full suite on 9f8f2b8: FAIL 0, PASS 2509, SKIP 8 (23.4 min).
- Smoke gate: SMOKE_OK for both arms. Decision-rule self-test: SELFTEST_OK.
- Results CSV regenerated from receipts: ROWS_MATCH.
- Review panel: method audit (Sonnet; hand re-derivation matched to 4 decimals; three
  rule-code deviations fixed), traceability (Sonnet; 0 mismatches), claim gate (Fable;
  1 blocking and 14 required items, all addressed).
- gate-check --reverify: see the PR description.

## 6. Tests of the Tests

- G2b negative control: the stratum-size fields are absent on origin/main
  (LACKS_FIELDS) and present on the branch (FIELDS_OK).
- Decision-rule self-test exercises both verdicts on synthetic fixtures.
- Manual gate M2: the method auditor re-derived PanTHERIA gestation_d far-stratum
  coverage and paired MCSE from the raw receipts with independent code; they match the
  table to 4 decimals.
- The fallback test was first changed to expect the realised stratum sizes and seen to
  fail (FAIL 2) before the fix.

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
- The acceptance ledger was first run without `--cwd`, so its checks ran from `.unlazy/`
  and could not find repository files; approvals then stayed bound to that directory and
  `--approve` skipped gates the ledger already marked met. Fixed by resetting runnable
  gates to pending and approving from the worktree root.
- Fan-out exceeded the plan's cap of six new children without a recorded amendment
  (about ten: recon, instrumentation, MI draws, paper, kohaku install, results builder,
  MI-GLS builder, method audit, traceability, claim gate, reconciliation). Each was a
  bounded Sonnet or Haiku slice except the planned Fable claim gate. Recorded as drift in
  the plan-actual reconciliation.

## 10. Known Residuals

- Coverage on masked observed cells is a proxy for coverage on cells users impute.
- FishBase has one mask and is descriptive only.
- The rule's near non-inferiority condition penalises removal of over-coverage; the
  decision stands as registered.
- Six of nineteen near-stratum trait rows fall below 0.95 under Mondrian.
- Ordinal traits were masked but not scored.

## 11. Team Learning

- Size Slurm `--mem` from a measured peak on the target dataset, never from a smaller one.
- On fir, cu126 libtorch needs `module load cuda/12.6` and
  `LD_LIBRARY_PATH=$EBROOTCUDA/targets/x86_64-linux/lib`, even on CPU nodes.
- A pool-and-GLS downstream check belongs in any MI validation; marginal draw
  calibration and OLS can both look fine while GLS is badly biased.

## 12. Cross-Product Coverage

Covers: single-observation continuous traits and one count trait (litter_size), gnn on,
three databases, two mask arms. Ordinal traits were masked but not scored. Does NOT cover: multi-observation data (Mondrian stops there by design), discrete
traits (no conformal intervals), `gnn = FALSE` (Mondrian stops there), MNAR beyond what
the structured arm captures, trees other than the three used.
