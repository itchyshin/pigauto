# After-task report: exact MI run provenance, 2026-10-09

## 1. Goal

Make the exact-tarball optional MI adapter matrix independently replayable from the audit PR by retaining the original run manifest and logs with their recorded hashes.

## 2. Implemented

Added the original four-cell manifest and four original logs to the audit provenance directory. Updated the matrix receipt, release ledger, and gate record with the manifest hash, original-log hashes, test script hash, and bounded Grace and Rose verdicts. The frozen archive and package source were not changed.

## 3a. Decisions and Rejected Alternatives

Grace and Rose both passed this optional-adapter subcheck after checking the archive and run-level provenance. Retaining the original manifest and full logs is necessary because the previously committed summary logs omit install commands, backend visibility, and before/after archive hashes. The evidence remains limited to these four installed-library cells and does not address the separate Windows failure.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-manifest-2026-10-09.json`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-neither-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-drmTMB-only-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-gllvmTMB-only-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-both-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-backend-matrix-2026-10-09.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-exact-mi-manifest.md`

## 5. Checks Run

- Verified the manifest SHA-256 is `fb9af9cbd7337358497ac8c72ecfcb8256f9504d01f8f696fac8ccdd39fb062e`.
- Verified its archive path identifies `pigauto_0.11.0.tar.gz`, 5,135,324 bytes, SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`.
- Verified all four original run-log hashes match the manifest and each before/after archive hash matches the frozen tarball.
- Verified test script SHA-256 `d2b0b5f01e3470ce1f0f0d2116cadbd1a9a2beffcd864a5644e71c0b7d53fb7d` matches `script/cran-0.11-integration/check-adapters.R` at commit `2c212ab`.
- Grace and Rose independently reviewed the bounded matrix and run-level provenance; both returned PASS for the optional-adapter subcheck.
- `check-after-task.R` passed its required-section check. Its full closure check remains blocked: the CRAN scope ledger has G3/G4 open, and the scope wrapper reports approval required for G6. The explicit status check showed 5/7 gates met and G3/G4 unmet; no approval-gated command was run.
- The repository-wide after-task check also finds open leaves in the separate imputation simulation scope. Those files and campaign gates are outside this audit lane and were not changed.
- Release status remains NOT READY because Windows logs are not bound to this archive and exact-artifact release closure is incomplete.

## 6. Tests of the Tests

No package behavior or test script changed. The independent reviewers checked that the real backend cells use three successful fits per backend, compare drmTMB output to a Gaussian oracle and pooled estimates/standard errors to Rubin calculations, and retain saved-fit reload checks. The neither-backend cell checks fixture behavior, provenance refusal, and the missing-package guard.

## 7a. Issue Ledger

No new package issue was found. The original manifest and logs improve replayability of the bounded MI integration check. Windows and broader release blockers remain open.

## 8. Consistency Audit

The manifest hash and four original-log hashes are recorded in both the matrix receipt and release ledger. Archive identity and test script identity agree across the manifest, receipt, and evidence checkout. G8 and G9 remain unmet, and the overall status claim remains NOT READY.

## 9. What Did Not Go Smoothly

The managed worktree was outside the default writable roots, so the copy and edits required the approved worktree-write escalation. The first attempted default write was rejected without changing files. The initial lane-preflight invocation also used the repository name instead of its path; rerunning it with the current worktree path showed the audit lane's active lease and a foreign direct-to-main lane. Work remained limited to the audit lane's leased paths.

## 10. Known Residuals

The MI evidence does not certify all downstream models, pooled estimands, or imputation regimes. Windows results still report failures and lack archive-hash binding. G8 and G9 remain open. The evidence PR remains Draft and unmerged; no CRAN submission occurred.

## 11. Team Learning

Grace and Rose independently reviewed the same bounded integration claim and its archive provenance. Their separate review records strengthen confidence in this subcheck without converting the overall release verdict to ready.

## 12. Cross-Product Coverage

Covers: frozen archive identity, optional backend presence and absence, four installed-library configurations, original test script and log provenance, bounded reviewer closure.

Does not cover: Windows compatibility, full exact-archive release checks, every downstream model or pooled estimand, deployment changes, merge, or CRAN acceptance.
