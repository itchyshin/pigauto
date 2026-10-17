# After-task report: current-hash Win-builder submissions and evidence reconciliation

## 1. Goal
Continue the approved pigauto 0.11 CRAN audit, preserve the frozen exact artifact, and advance the platform gate using authorized Win-builder checks.

## 2. Implemented
Submitted the frozen `pigauto_0.11.0.tar.gz` (SHA-256 `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`, 5,136,083 bytes) to R-release and R-devel without rebuilding. Both returned HTTP 200 and displayed the exact filename and size. The R-devel form says “Currently unavailable.” Recorded the acknowledgements in the source binding, release ledger, and GATES; retained the previous artifact incident report as a superseded historical record.

## 3a. Decisions and Rejected Alternatives
The prior plan authorizes exact-artifact platform checks and evidence updates. The archive, source, website, dependencies, and package behavior were not changed. PR #236 remains draft/unmerged; no site deployment or CRAN submission was attempted.

## 4. Files Touched
- `docs/dev-log/cran-0.11-audit/provenance/winbuilder-upload-fe1600d7-2026-10-17.md`
- `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-source-binding.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-17-artifact-replacement.md` (supersession notice)
- This report.

## 5. Checks Run
- Recomputed the frozen archive SHA-256 after upload: `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`.
- Official R-release form response: HTTP 200; exact filename and 5,136,083-byte size displayed.
- Official R-devel form response: HTTP 200; exact filename and 5,136,083-byte size displayed; form marked R-devel currently unavailable.
- Prior exact-artifact `R CMD check --as-cran --run-donttest --timings`: completed with `Status: 1 NOTE` for the 6.142-second `compare_methods` timing; 3,194 pass, 0 fail, 161 warnings, 86 skips. See source-binding receipt and hashed logs.
- The upload responses are acknowledgements only. No Windows check result or server-side digest was returned.

## 6. Tests of the Tests
No source or test code changed in this slice. Existing package-test and as-cran results remain bound to the same frozen archive and are cited above; this update does not rerun them.

## 7a. Issue Ledger
- Exact-hash local package gate: met at `tarball-clean`, with one timing NOTE disclosed.
- Windows/platform gate: unmet; upload acknowledgements are present, result logs are not.
- Exact-hash MI matrix: met within its bounded 12-step optional-backend workflow, as independently reviewed.
- Site gate: bounded receipts exist for the deployed site at source commit `780c084`; this slice did not replay the site audit.
- Independent closure and final release verdict: unmet. Overall release state remains `NOT_READY`.

## 8. Consistency Audit
The tarball hash and size agree across the local file, upload acknowledgement, source binding, and ledger. NOTICE plus all six shipped `.rda` streams were independently rehashed from the frozen archive and recorded in a dedicated receipt. The ledger and GATES distinguish acknowledgement from check result and keep `windows_results_bound_to_hash` false. Older replacement-incident prose is explicitly marked superseded. No source-tree modification was made.

## 9. What Did Not Go Smoothly
The first uploader call used version labels that the form did not recognize; it made no submission. The form’s actual labels are `R-release` and `R-devel`, and the corrected calls succeeded. The R-devel page reports that this service is currently unavailable. Gmail automation was blocked by an open Chrome extension UI, so no result email was read in this slice.

## 10. Known Residuals
No Windows result has arrived. The after-task structure validator also reported the co-located `.unlazy/imputation-sim` acceptance ledger as not runnable because its node or gate-check is missing; that ledger belongs to a separate scope and was left untouched. This is not a pass for that ledger. The service does not expose a server-side SHA-256, so a later result can only be tied by filename, exact byte size, timing, and any run-specific email/log metadata. R-devel availability needs follow-up. The local check has one timing NOTE. Live PR head/state has not been refreshed after the latest local evidence changes. No merge, deployment, or submission is claimed.

## 11. Team Learning
Grace independently verified the exact tarball and completed as-cran log, including the sole timing NOTE and all 236 package-file matches. The MI reviewer independently verified the bounded replacement-hash matrix. Rose correctly kept platform status open because no result was bound to this archive; upload acknowledgements do not change that verdict. Pat’s website review is bounded to its supplied route inventory and NEWS exception.

Memory receipt: used the pigauto release memory entry to preserve exact-artifact rung discipline and the no-merge/no-submission boundary; no brain files were written.

Golden Set: no new Golden Set item was created; this slice is release-evidence bookkeeping, not a model-behavior change.

## 12. Cross-Product Coverage
This slice covers exact-artifact upload acknowledgements and local evidence consistency. It does NOT cover Windows check completion, server-side archive identity, full deployed-site replay, a clean-room optional-dependency install, broader statistical recovery, platform-clean status, CRAN acceptance, merge, or publication.
