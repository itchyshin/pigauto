# After-task report: release-ledger validator and PR refresh, 2026-10-10

## 1. Goal

Make the CRAN audit ledger executable under its fail-closed validator and refresh the evidence PR receipt after the latest push.

## 2. Implemented

The ledger now uses the validator's rung vocabulary and links the exact archive, inventory, clean source checkout, local check, installed test output, policy, rights, and rendered-site receipts. It records `tarball-clean` as the highest proven rung, `platform-clean` as the next unproven rung, and keeps the overall audit verdict NOT READY. The source checkout status receipt binds the current clean checkout to the recorded source commit and archive identity.

## 3a. Decisions and Rejected Alternatives

The earlier free-form `status_claim` was invalid under the executable validator. It has been replaced by the highest rung supported by the exact local archive checks. The ledger does not claim platform-clean because Windows diagnostics fail and remain unbound to the archive hash. The full release remains NOT READY, and PR #228 stays Draft.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/post-merge-source-checkout-status-2026-10-10.txt`
- This report.

No package source, tests, or frozen archive changed.

## 5. Checks Run

- Verified source worktree commit `d76804e768bf77f430f43dcf591633f9cd900dab` and empty `git status --porcelain=v1 --untracked-files=all`.
- Recomputed archive SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`, size 5,135,324 bytes, and inventory SHA-256 `10a31a50b2504f42514b0ec8b5368861543cf1f902d9ae075ebb5279cd5a777e` for 250 entries.
- `cran_release_gate.py --selftest`: passed all planted negative controls.
- `cran_release_gate.py docs/dev-log/cran-0.11-audit/release-ledger.json`: `READY FOR CLAIMED RUNG` at `tarball-clean`.
- JSON parse and `git diff --check`: passed.
- Chrome observed PR #228 before this report commit at 66 commits, latest `427ef1f`, Draft and unmerged. Workflow run `38009821247` was skipped under the pull-request pkgdown guard; no artifact or deployment was listed. This report and ledger reconciliation only add audit documentation.
- Attempted to select the exact archive in Chrome's Win-builder upload form. The extension blocked file selection because it lacks access to file URLs. No upload occurred and no browser permission was changed.

## 6. Tests of the Tests

The executable validator's self-test exercised its negative controls for unlicensed data, predecessor hashes, forbidden tarball paths, missing incoming checks or timing, unconfirmed upload, fictional public URLs, panel status, and omitted conditional gates. All controls failed closed as expected.

## 7a. Issue Ledger

- Resolved: the release ledger now has a valid highest-rung claim and evidence objects accepted by the executable validator.
- Open: R-release and R-devel results for the exact archive are still needed.
- Open: independent exact-artifact review remains NOT READY.
- Open: Chrome file selection requires the extension's file-URL access permission before the authorized Win-builder upload can proceed.

## 8. Consistency Audit

The ledger's exact archive hash, size, source commit, and 250-entry inventory match the frozen build receipt. The next unproven rung is platform-clean. `audit_verdict` remains NOT READY, G8 and G9 remain open, and the PR remains unmerged.

## 9. What Did Not Go Smoothly

The first validator run exposed a ledger-schema mismatch: the file used prose where the validator requires a rung and omitted several typed fields. The ledger was updated and rechecked. Chrome then refused the local file chooser because its extension lacks file-URL access; the permission was left unchanged.

## 10. Known Residuals

The source-checkout status was captured on 2026-10-10 at the recorded commit. The build receipt identifies the same source path and commit, but this receipt does not claim a status command was captured during the original build. Windows results and independent exact-artifact closure remain outstanding.

## 11. Team Learning

Keep the release status as a named evidence rung and the overall submission verdict as a separate field. A validator pass for `tarball-clean` does not imply platform passage or submission readiness.

## 12. Cross-Product Coverage

Covers release-ledger schema, frozen archive identity, current source-checkout status, executable-validator controls, and current evidence-PR state.

Does not cover Windows results, source or site changes, deployment, merge, or CRAN submission.
