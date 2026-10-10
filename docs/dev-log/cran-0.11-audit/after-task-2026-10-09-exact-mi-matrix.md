## 1. Goal

Bind the optional multiple-imputation backend matrix to the exact frozen pigauto 0.11.0 source archive and preserve reviewable evidence.

## 2. Implemented

Recorded the four installed-library cells (both optional model packages, each package alone, and neither), retained five raw logs with hashes, added an artifact-specific receipt, and linked that receipt from the release ledger and gate record. No package source or frozen tarball changed.

## 3a. Decisions and Rejected Alternatives

Reused the completed exact-archive matrix after verifying the archive hash, per-cell logs, installed package paths, package version, and optional-package visibility in separate R sessions. A further fit campaign was unnecessary for this bounded adapter check. Assumption: the hash-qualified archive in the current release ledger is the intended candidate; if that ledger points to the wrong archive, this evidence must be retargeted and rerun.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-backend-matrix-2026-10-09.md`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-both-backends-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-neither-fixtures-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-neither-reload-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-drmTMB-only-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/exact-mi-gllvmTMB-only-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-exact-mi-matrix.md`

## 5. Checks Run

- Recomputed SHA-256 of the frozen archive before recording: `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`; size 5,135,324 bytes. The archive is recorded as source commit `d76804e768bf77f430f43dcf591633f9cd900dab` with 250 entries.
- The installed matrix logs report: both backends, 50 expectations; drmTMB only, 28; gllvmTMB only, 22; neither, 50 adapter fixtures plus 55 provenance expectations and the saved-object missing-package guard. All real fits reported convergence 0 and positive-definite Hessians.
- Fresh R 4.6.0 sessions loaded pigauto 0.11.0 from each cell's own library and confirmed backend visibility: both, neither, drmTMB only, and gllvmTMB only.
- `python3 -m json.tool docs/dev-log/cran-0.11-audit/release-ledger.json`: passed.
- `slop_check.py` on the matrix receipt: 0 findings.
- `node .../gate-check.mjs --status .../GATES.md`: 11 gates, 9 met and 2 unmet (G8, G9). No release-completion claim is made.

## 6. Tests of the Tests

The fixture and provenance suites include adapter-availability, malformed-object, and provenance checks; the real-backend cells compare pooled output with independent Gaussian and Rubin-rule calculations. No deliberate code mutation was used to measure fault sensitivity during this evidence-only slice.

## 7a. Issue Ledger

No new package defect was found. The exact-artifact optional-adapter subcheck passes. Existing Windows logs report five failed expectations in each of R-release and R-devel but do not bind those results to this archive hash; they remain diagnostic, not exact-hash platform evidence. Independent exact-artifact review is also pending.

## 8. Consistency Audit

Confirmed the receipt's artifact identity matches `release-ledger.json`, the four library paths match their declared backend cells, all five copied logs match their recorded hashes, and the gate record leaves G8 and G9 open. The package remains independent of drmTMB and gllvmTMB at installation time while supporting their optional fixed-effect adapters when users install them.

## 9. What Did Not Go Smoothly

The first library visibility probe reused one R process, so namespace caching made its later labels invalid. I replaced it with four fresh R processes and confirmed each cell. The managed worktree initially rejected writes; the authorized worktree write was then approved. A relative closeout path also resolved under the brain directory; I removed only the new skeleton and recreated the report at the absolute worktree path.

## 10. Known Residuals

G8 remains open for checksum-bound Windows results and remaining exact-artifact release checks. G9 remains open for independent verdicts on the exact archive and site evidence. This matrix does not establish broad inferential validity for every model or estimand. The evidence PR remains unmerged and nothing was submitted to CRAN.

## 11. Team Learning

Memory receipt: `route.py pigauto` and the repository LOAD-FIRST manifest were read. Their emphasis on exact artifact identity and distinguishing measured results from claims shaped this receipt. No Golden Set lookup was needed because this slice changed evidence records, not prediction code.

Golden Set: not checked; no known-mistake class in package behavior was modified.

## 12. Cross-Product Coverage

Covers: exact-archive installation, optional backend presence/absence, automatic fixed-effect adapters, the multiple-imputation to downstream-fit to `pool_mi()` path, saved-object reload behavior, and independent arithmetic checks.

Does NOT cover: Windows compatibility, a full exact-archive platform matrix, every downstream model family or pooled estimand, all trait/imputation regimes, or CRAN acceptance.
