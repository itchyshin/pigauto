# Pat final release review

**Verdict: READY** for the stated 0.11.0 release gates, bound to archive SHA-256 `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb` from source commit `b2e6054685e6bdedb073ba46e058af2c0a619a18`.

The content-addressed archive read back at 5,137,837 bytes with the expected SHA-256. The exact archive's `impute` help defaults to `gnn = FALSE`, and its report help describes baseline-only output for that default and GNN comparison when `gnn = TRUE`. Installed help receipts agree. The installed novice workflow configured the GNN off, filled both missing cells, and preserved observed values and their classes. The report and impute descriptions therefore agree with the repaired source and the deployed-page receipt.

The records for four native checks report `Status: OK`, zero errors, warnings, and NOTEs. The separate optional-model matrix passed its 12 recorded steps. Six vignette diagnostics passed with disabled teaching chunks left disabled. The live-site receipt binds the reviewed deployment to this source commit and records successful route checks plus Chrome review of the changed references and Getting Started page.

There are no blockers in this review. The evidence describes bounded checks and does not claim broad statistical recovery, interval coverage, clean-room dependency installation, CRAN timing guarantees, or CRAN acceptance. The disclosed DOI resolver 403 remains an access denial; Crossref confirmed the article identity. The separate warning capture has two coherent timing records, but the identical RDS bytes do not identify which capture produced the audited objects, and overwritten history may be incomplete. These are accurately disclosed limits, not contradictions in the native check results.

## Evidence checked

- Frozen archive readback, size, source identity, inventory, and independent freeze review.
- Actual installed `impute` and `pigauto_report` help text, plus both corresponding Rd files in the frozen archive.
- Installed novice workflow receipt, output, and run log.
- Native platform summary and direct logs, optional-model matrix supervisor, and six-vignette timing receipt.
- Deployed-site verdict, live route receipts, and final CRAN comments.
- DOI and Crossref disclosure, warning-content audit and timing provenance amendment, current policy receipt, and rights/source summary.

No source, archive, installed package, browser state, or external service was changed during this review. This review did not rerun checks, fits, captures, or simulations.
