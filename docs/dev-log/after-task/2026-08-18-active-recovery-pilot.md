## 1. Goal

Run the explicitly approved, pre-registered Stage-A 100-replicate-per-cell
Totoro pilot and produce a durable, non-claim-bearing pilot receipt.

## 2. Implemented

Ran eight cells at source SHA `beee5df`, using 96 workers and one BLAS/OpenMP
thread per worker.  The campaign retained 800 per-replicate RDS receipts and
generated its registered paired summary.  Updated the pilot results report,
plan-versus-actual reconciliation, and validation ledger.

## 3a. Decisions and Rejected Alternatives

Kept the locked 60/20/20 split, one-cell restoration, 100 epochs, and all
four policies.  Did not start the 1,000-replicate campaign, alter defaults,
or turn the pilot pattern into a public benefit claim: the binary lambda=1
pilot intervals do not yet reliably favour active selection.

## 4. Files Touched

- `script/active_recovery/03_launch_totoro_pilot.sh`
- `docs/dev-log/2026-08-18-active-recovery-results.md`
- `docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md`
- `VALIDATION_LEDGER.md`
- `docs/dev-log/after-task/2026-08-18-active-recovery-pilot.md`

Raw receipts, configuration, session information, launcher log, and summary
are retained on Totoro at
`/home/snakagaw/pigauto-active-recovery-results/active-recovery-pilot-20260818-beee5df`.

## 5. Checks Run

`bash -n script/active_recovery/03_launch_totoro_pilot.sh` and `git diff
--check` passed.  Totoro loaded the package at the recorded SHA before launch.
The independent post-run audit required 800 receipt files, 3,200 rows, four
policy rows per replicate, one source SHA, zero error rows, correct changed
vs unchanged data hashes, and an empty campaign process group; all passed.

## 6. Tests of the Tests

The receipt audit would fail on a missing replicate, a partial four-policy
set, a source mismatch, an errored fit, a no-acquisition treatment mutation,
or a restored-policy data hash that matched the initial data.  It tested raw
receipts independently of the campaign summariser.

## 7a. Issue Ledger

The first detached launch failed before compute because its output parent
directory did not exist.  It left no process or campaign directory; the
parent was created and the verified launch then completed.  Deferred: the
full Stage-A campaign requires post-pilot approval; Stage B awaits PR #173;
Stage C remains protocol-only.

## 8. Consistency Audit

Checked the pre-registration's eight cells, pilot count, source receipt,
policy treatment rule, failure retention, and full-campaign formula against
the generated outputs.  All cells retained 100 complete replicates, no
failure rate exceeded zero, and the 1,000-replicate minimum applies to every
cell.  No public package code or defaults changed.

## 9. What Did Not Go Smoothly

The generic Totoro campaign-status wrapper lacked its profile file, so the
canonical SSH connection and an explicit campaign-specific launcher were used.
One monitoring hash comparison initially recycled differently sized vectors;
it was corrected to aligned policy subsets before being treated as evidence.

## 10. Known Residuals

This is a pilot only.  It cannot establish the headline active-imputation
claim, particularly for binary lambda=1 conditions.  No full 1,000-replicate
campaign, Stage-B external comparison, or Stage-C real-data Mondrian study
has been run.  PR #173 remains open and its checks were still running when
this report was prepared.

## 11. Team Learning

Memory receipt: loaded the repository LOAD-FIRST compute, recovery-to-truth,
and `r_cal = 0` guards; they required the 96-worker ceiling, single-threaded
BLAS, retained failures, and no default changes.  The durable operational
result is that single-threaded BLAS plus a receipt per replicate allows a
96-worker Totoro pilot to be audited without hidden failures or orphaned
processes.  The brain vault was not modified: no separate approval was given
for a vault write.  Golden Set: not in scope; this work exercised a new
evidence harness rather than a known cross-repo regression class.

## 12. Cross-Product Coverage

Covers: continuous and binary single-observation families; n=100 and n=300;
lambda=1 and lambda=0.2; active/random/uncertainty/no-acquisition policies;
source, treatment, failure, and process-cleanup receipts.

does NOT cover: a full Monte-Carlo claim; mixed-type, multi-observation, or
batch acquisition; real-data benefit; GNN attribution; default changes;
Stage-B exact-solver competitiveness; or Stage-C interval calibration.
