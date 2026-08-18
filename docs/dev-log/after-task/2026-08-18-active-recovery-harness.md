## 1. Goal

Implement the pre-registered Stage-A active-recovery harness and integrity smoke for pigauto's next evidence programme.

## 2. Implemented

Added a deterministic Stage-A harness under `script/active_recovery/`: the
locked DGP, 60/20/20 split, common-candidate active/random/uncertainty policies,
one-cell restoration, per-replicate RDS receipts, and paired summary script.
Added tests for split protection, exact restoration, changed-data receipts,
decision-time uncertainty, continuous/binary active-ranking maximality, and
finite fixed-test predictions. Added pre-registrations for Stages A--C and one
experimental validation-ledger entry.

## 3a. Decisions and Rejected Alternatives

The harness uses cell-level, single-observation rankings and reports continuous
and binary families separately. It does not pool variance and entropy metrics,
rank species, make a static top-10 claim, or attribute an effect to the GNN.
Stage B remains a protocol because exact-solver PR #173 is not yet green and
merged; Stage C remains a protocol because real-data masking is a separate
calibration claim.

## 4. Files Touched

- `script/active_recovery/00_prepare_active_recovery.R`
- `script/active_recovery/01_run_active_recovery_cell.R`
- `script/active_recovery/02_summarise_active_recovery.R`
- `tests/testthat/test-active-recovery.R`
- `VALIDATION_LEDGER.md`
- `docs/dev-log/2026-08-18-active-recovery-{design,results}.md`
- `docs/dev-log/2026-08-18-{exact-comparator,mondrian-realdata}-protocol.md`
- `docs/dev-log/2026-08-18-active-prior-art-audit.md`
- `docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md`

## 5. Checks Run

`Rscript -e 'parse(...)'` passed for all Stage-A scripts and the test file.
The pure split/restoration smoke passed. Continuous and binary 30-tip,
five-epoch structural runs completed with finite protected-test predictions and
changed data receipts. `testthat::test_file("tests/testthat/test-active-recovery.R")`
passed 12 expectations. `devtools::test(filter = "active")` completed without
test failures (existing small-sample/trait-detection warnings remain).

## 6. Tests of the Tests

The restoration test would fail if a test cell were made observed or if zero or
multiple cells changed. The selection tests independently reconstruct the
candidate-filtered public ranking and would fail if the selected continuous or
binary candidate were not maximal. The finite-prediction assertions would fail
for all-NA or non-finite fixed-test output.

## 7a. Issue Ledger

Fixed: the first harness version did not retain a no-acquisition secondary row;
it now does. Deferred: the 100-replicate-per-cell Totoro pilot, full campaign,
post-merge comparator execution, and real-data Mondrian execution.

## 8. Consistency Audit

Checked the existing active-imputation formula tests, prediction output shapes,
the validation ledger, and FishBase/PanTHERIA benchmark conventions. The Stage-A
runner does not modify public package defaults or BACE. The repository's
existing active tests already verify the independent fixed-sigma BM and LP
entropy mathematics; this slice verifies treatment application to a common,
fixed test set.

## 9. What Did Not Go Smoothly

NotebookLM's token check could not resolve DNS inside the sandbox; the approved
outside-sandbox retry authenticated successfully and started the quarantined
third-party research task. The closeout helper initially resolved a relative
path into the Shinichi hub, so it was rerun with an explicit worktree path.

## 10. Known Residuals

No pilot or full evidence result exists. The local 30-tip smoke's two-cell
validation warning means it is not an effect-size or calibration result. The
third-party scan is now primary-source audited but establishes adjacent work,
not novelty. The source-content receipt is a deterministic two-part checksum,
not a cryptographic hash.

## 11. Team Learning

Memory receipt: recalled the repo handover and validation ledger; applied the
LOAD-FIRST guards for prediction-path correctness, `r_cal = 0` preservation,
and compute gating. The brain was not modified because no separate approval to
write it was given.

Golden Set: not in scope; this adds campaign machinery rather than a known
cross-repo regression class.

## 12. Cross-Product Coverage

Covers: continuous and binary single-observation, one-cell active acquisition;
candidate/test disjointness; matched random and uncertainty comparators;
per-replicate failures and source receipts.

Does NOT cover: mixed-type or species-level rankings, batch selection, real-data
benefit, GNN attribution, a full Monte-Carlo result, Stage-B external parity,
or Stage-C real-data interval calibration.
