## 1. Goal

Run the explicitly approved 7,200-replicate Stage-A extension on Totoro,
combine it with the audited pilot, and report the registered active-recovery
claim gate without broadening the result.

## 2. Implemented

Ran replicates 101--1,000 in every one of the eight registered cells at
extension SHA `ad3990c`, using 96 workers with all BLAS/OpenMP thread caps at
one.  Combined these receipts with the SHA `beee5df` pilot, producing 8,000
receipts and 32,000 policy rows.  Added an independent receipt-audit script
and updated the results, reconciliation, and validation ledger.

## 3a. Decisions and Rejected Alternatives

Kept the locked DGP, 60/20/20 split, one-cell budget, same fitting seed,
four policies, 100 epochs, and all cells.  Did not pool continuous and binary
outcomes, report an interim result, alter defaults, or claim cross-type
active-imputation benefit.  The registered broad gate fails because binary
lambda=1 paired intervals include zero.

## 4. Files Touched

- `script/active_recovery/06_audit_active_recovery_receipts.R`
- `docs/dev-log/2026-08-18-active-recovery-results.md`
- `docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md`
- `VALIDATION_LEDGER.md`
- `docs/dev-log/after-task/2026-08-18-active-recovery-full.md`

Raw pilot and extension receipts, config, session information, launcher log,
independent `receipt_audit.rds`, and `combined_summary.rds` are retained on
Totoro under `/home/snakagaw/pigauto-active-recovery-results/`.

## 5. Checks Run

The extension completed with 7,200 receipts and zero logged errors; no worker
remained.  The independent audit required 800 pilot plus 7,200 extension
receipts, 32,000 policy rows, exact registered replicate grid, four policies
per replicate, two expected source SHAs, zero error rows, finite primary
scores, correct split sizes, and correct changed/unchanged treatment hashes;
it passed.  The registered combined summary then produced 1,000 paired
replicates per cell and zero failure rate for every policy.

## 6. Tests of the Tests

The receipt audit fails closed on a missing receipt, incomplete policy set,
wrong source SHA, errored fit, invalid grid, non-finite score, mutation in
the no-acquisition branch, unchanged acquisition data, missing selected
species, or split-size mismatch.  It is separate from the campaign summary.

## 7a. Issue Ledger

No campaign execution failures occurred.  The initially planned broad
headline is not supported: binary lambda=1 is inconclusive.  Deferred:
Stage B requires the green PR #173 to be merged; Stage C requires executable
real-data confirmation and a separately approved compute run.

## 8. Consistency Audit

The receipts match the registered eight cells, 1,000 replicates per cell,
one-cell restoration, fixed-test protection, and source/treatment receipt
requirements.  The claim verdict matches the registered rule: continuous
lambda=1 passes against random at both sizes, but binary lambda=1 does not;
therefore no cross-family headline is made.

## 9. What Did Not Go Smoothly

The full run took about 56 minutes rather than the 1.04-hour estimate but
remained within it.  Expected n=100 small-validation warnings were retained
in the logs; they were not treated as errors or removed from evidence.

## 10. Known Residuals

This evidence supports only continuous BM, lambda=1, single-observation,
cell-level, one-step recovery versus random acquisition.  It does not show
binary efficacy, superiority over uncertainty-only selection, weak-signal
benefit, mixed-type or batch optimality, real-data benefit, GNN attribution,
or a default change.

## 11. Team Learning

Memory receipt: the pilot-derived wall-time estimate, frozen source bridge,
single-threaded 96-worker Totoro launch, and per-replicate RDS receipts made
the full campaign auditable without hidden failures.  Independent claim
review prevented the favourable continuous result from becoming an unsupported
cross-family statement.  The brain vault was not modified: no separate
approval was given for a vault write.  Golden Set: not in scope.

## 12. Cross-Product Coverage

Covers continuous and binary single-observation families, n=100/300,
lambda=1/.2, active/random/uncertainty/no-acquisition policies, 1,000
replicates per cell, and source/treatment/failure/process receipts.

Does NOT cover mixed types, multi-observation data, species or batch ranking,
real data, external comparators, exact-solver competitiveness, or Mondrian
interval calibration.
