# After-task report: current artifact and site reconciliation

## 1. Goal

Record the October 15 exact-artifact checks and reconcile the release ledger with the newest deployment and the newly found Getting Started discrepancy.

## 2. Implemented

- Updated the ledger to the 0.11.0 tarball built from merged source `5e46c46` and retained its checksum-bound local receipts.
- Corrected the testthat output pointer and hash to the separate retained `testthat.Rout` file.
- Refreshed deployed-site evidence to the latest verified commit and recorded the Getting Started mismatch.
- Added the scoped source correction PR #234, local article-build result, independent review, and limitations.
- Marked the 4b32f669 tarball as historical for final release if PR #234 is merged.

## 3a. Decisions and Rejected Alternatives

- Kept PR #228 unmerged because its exact-artifact and release verdict remain NOT READY.
- Did not call the prior archive the final candidate after discovering the reader-facing defect.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-15-current-artifact-refresh.md`

## 5. Checks Run

- JSON parse with `python3 -m json.tool`: passed.
- `cran_release_gate.py --selftest`: passed all positive and negative controls.
- `cran_release_gate.py release-ledger.json`: ready for the claimed source-clean rung only.
- Exact offline archive check: passed for hash `4b32f669bb92474f6ec9acb45eb7c7806467ed104bb4338371e4b14369cf230f`.
- Separate `--as-cran` attempt: blocked before package tests by CRAN/Bioconductor DNS.
- Four-cell optional backend matrix: passed for that hash.
- Chrome verified deployment `5e46c46` and identified a reader inconsistency; PR #234 checks were queued at the last observation.

## 6. Tests of the Tests

The release validator's self-test passed negative controls for stale hashes, rights, forbidden archive entries, missing incoming checks, missing timing, unconfirmed uploads, retired public URLs, omitted compiled-code gates, fictional live URLs, reachable non-CRAN URLs, incomplete panel review, invalid release type, and omitted vignette budgets. These controls validate ledger rejection behavior, not the release artifact itself.

## 7a. Issue Ledger

- Fixed: current ledger pointed to an older tarball and deployment.
- Fixed: testthat output location was described as embedded in the check log.
- Open: source correction PR #234, checks, merge, deployment, and post-deployment review.
- Open: checksum-bound Windows checks and a network-capable `--as-cran` incoming check on the final artifact.
- Open: independent review of the final artifact, not its predecessor.
- Open: redistribution evidence for bundled tree data remains bounded by the maintainer's confirmation and cited public pages.

## 8. Consistency Audit

Compared the archive identity, inventory, check logs, testthat output, MI matrix manifest, site deployment receipt, and release ledger. Their hash and path references now agree for the predecessor artifact. The deployed page mismatch is separately recorded and its correction is in PR #234. The predecessor artifact is not promoted to final status.

## 9. What Did Not Go Smoothly

The fresh `--as-cran` attempt stopped at incoming feasibility because package-index DNS failed. PR checks were queued during the final browser observation. A later correction to the Getting Started source means the current frozen tarball cannot serve as the final release artifact if that correction is merged.

## 10. Known Residuals

- Highest proven rung is source-clean for the merged predecessor source; tarball-clean for the final corrected source is not established.
- Windows checks are not checksum-bound to the current archive.
- The network-dependent `--as-cran` gate did not complete.
- Final exact-hash review and deployed-site reader consistency remain open.
- The rights record is maintainer confirmation with attribution; this is not an independent legal determination.
- No merge of the evidence PR or CRAN submission occurred.

## 11. Team Learning

Treat an independently found page mismatch as a new source correction and rebuild boundary. Keep exact-artifact results attached to the exact hash that was actually tested, then rerun them after source changes.

### Memory receipt

Loaded the pigauto route manifest, repository AGENTS instructions, and after-task protocol. The route guard to diff main and preserve the predecessor artifact identity shaped this work. No cross-project scouting was done. The Golden Set was not run because no package code changed.

## 12. Cross-Product Coverage

- Covers: the predecessor tarball identity, its bounded local checks, optional MI backend matrix, current deployed route controls, and the Getting Started source correction status.
- Does NOT cover: the final merged source, a final tarball, Windows checks, a complete network-capable `--as-cran` run, final deployed-page review, or general statistical validity.
