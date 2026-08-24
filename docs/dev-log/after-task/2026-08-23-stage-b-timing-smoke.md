## 1. Goal

Run and retain a bounded one-mask Stage-B timing smoke after the exact
conditional solver merged, without making a comparator-performance claim.

## 2. Implemented

Merged current `main` into `codex/active-recovery-evidence`, refreshed the
isolated Totoro checkout, verified `impute()` exposes `predict_method`, and
completed one locked-mask continuous and mixed timing receipt. Recorded the
result and the Tamia four-GPU ladder correction in the evidence ledger.

## 3a. Decisions and Rejected Alternatives

Used a normal merge of `origin/main` rather than rewriting the published
evidence branch. Retained BACE as an unavailable error rather than silently
dropping it. Rejected a full five-mask run because BACE has no verified source
or timing receipt in the isolated library.

## 4. Files Touched

- `docs/dev-log/2026-08-19-mondrian-scalability-correction.md`
- `docs/dev-log/2026-08-23-stage-b-timing-smoke.md`
- `docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md`
- `docs/dev-log/after-task/2026-08-23-stage-b-timing-smoke.md`

## 5. Checks Run

- Local exact-conditional test file: passed.
- Totoro private-library check: `impute()` contained `predict_method`.
- Continuous receipt: default and exact pigauto arms completed with no error
  rows; aggregate and per-mask RDS files exist.
- Mixed receipt: default and exact pigauto arms completed with no error rows;
  aggregate and per-mask RDS files exist.
- `git diff --check`: passed before documentation commits.

## 6. Tests of the Tests

The initial smoke failed first when `PIGAUTO_OUT_DIR` was omitted, then when
source loading required an absent `devtools`; each error was retained and a
fresh output directory was used for the repaired attempt. The previous branch
lacked the merged exact route and produced `unused argument (predict_method)`;
the refreshed checkout's explicit formals check would have caught that same
parameter-drop class again.

## 7a. Issue Ledger

- Fixed: timing checkout did not contain merged PR #173; merged `main`, pushed,
  refreshed, and reinstalled in the isolated library.
- Fixed: source-loading path unnecessarily required `devtools`; used the
  already verified installed package for the receipt.
- Deferred: BACE is unavailable; CRAN and public `itchyshin` repository
  searches did not yield a source. Full Stage B cannot proceed unchanged.

## 8. Consistency Audit

Checked both continuous and mixed runners, their locked seed behavior,
separate receipt paths, process cap, primary pigauto method rows, and all
external method error rows. `Rphylopars`, `phylolm`, and `missForest` were
present and ran; BACE was consistently retained as unavailable.

## 9. What Did Not Go Smoothly

The first two smoke attempts were setup failures, not model results: missing
output configuration and absent source-loading dependency. The original
evidence branch also lagged `main`, so it lacked the exact route until merged.

## 10. Known Residuals

This is one-mask timing evidence only. It does NOT cover five-mask precision,
BACE performance/timing, competitiveness, parity, a default change, or any
scientific conclusion. Stage C's completed ladder likewise does NOT cover
calibration or a CPU-versus-GPU comparison.

## 11. Team Learning

For a post-merge benchmark, verify the exact argument in the remote installed
package before spending compute; a branch name alone is insufficient. A
private, minimal library makes missing optional comparators visible early, so
retain unavailable rows and resolve the source before costing a full campaign.
Memory receipt: the routed LOAD-FIRST manifest and lane preflight shaped this
work; the prediction-path correctness and retained-failure guards were applied.
Golden Set: not in scope because no package implementation was changed.

## 12. Cross-Product Coverage

`predict_method` covers ✓ the installed `impute()` formals check and both
continuous/mixed timing runners. It does NOT cover five masks, alternative
datasets, BACE availability, performance equivalence, interval calibration,
or user-facing defaults. GPU ladder coverage ✓ distinct device assignment and
bounded split inputs; it does NOT cover GPU acceleration, CPU comparison, or
a full Stage-C campaign.
