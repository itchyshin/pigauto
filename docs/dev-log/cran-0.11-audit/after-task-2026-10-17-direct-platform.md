# Direct frozen-archive platform checks, 2026-10-17

## 1. Goal

Close the bounded platform-clean evidence gap for PR #236 by checking the original frozen pigauto 0.11.0 gzip archive on Windows release, Windows devel and an independent Linux runner. Preserve separate merge, deployment and submission decisions.

## 2. Implemented

Added a dedicated platform workflow and build-excluded helpers. It downloads the original archive, verifies SHA-256, size and inventory, installs all declared dependencies and check tools, exercises libtorch, and checks the original absolute gzip path. Identity is checked before and after execution, including failures. Raw evidence is retained on every job.

The dependency parser now removes empty fields from the frozen DESCRIPTION's trailing comma and is tested against R's dependency parser. Manual-tool setup explicitly supplies missing Courier metrics, V8 and Linux HTML Tidy. A font compilation canary and JavaScript engine canary fail before the package check if those tools do not work. The later index probe records executable paths and command results, conditionally installs makeindex only when absent, and requires a nonempty index and an indexed PDF through R's own tools.

## 3a. Decisions and Rejected Alternatives

Kept the existing frozen archive unchanged. Did not replace its bytes with a rebuild, disable manual checks, waive HTML validation, or treat R-hub's rebuilt-source results as checks of the original gzip. All fixes are checker or evidence changes.

Assumptions made without asking: ordinary ephemeral GitHub runners are acceptable for these authorized public package checks. This is bounded package validation, not a simulation campaign. The first cold estimate was 45-75 minutes per job. Attempt 2 measured about 14-16 minutes; retry estimate was 20-40 minutes including setup and review.

## 4. Files Touched

- `.github/workflows/frozen-tarball-check.yaml`.
- `script/cran-0.11-platform/`: `artifact.py`, `check.R`, `test-dependencies.R`, `verify_results.py`, `manual-tools.R`, `README.md`.
- `docs/dev-log/cran-0.11-audit/release-ledger.json` and `GATES.md`.
- This report and `provenance/platform-c2e768ad-direct-2026-10-17/`. The retained manifests enumerate each evidence path and its digest. The current-artifact evidence bundle also receives a correction note above its verbatim historical snapshots.

## 5. Checks Run

Direct original-gzip run 39070922656 completed Status OK on Windows release, Windows devel and Linux release. The matrix verifier returned DIRECT_ARCHIVE_MATRIX_VERIFIED_3_OF_3. The independent reviewer verified all 96 retained raw/compressed files, actual manual indexes, all 29 dependencies, archive identity and source/run linkage. Tests report 3,195 passing expectations, zero failures, 161 unclassified warnings and 86 skips on each runner. Totals: 599.377, 656.047 and 324.048 seconds, with 379, 417 and 192 seconds for tests.

Source CI 359, run 39070922883, passed all three jobs at checker head a550071, including focused macOS prediction tests. Python/R/YAML parsing and diff whitespace checks passed. The portable release-validator's eight tests and selftest passed. Claimed-rung execution and prepublication unlazy reverification are retained in the platform receipt. The final publication/readback result is recorded in the PR description. The final evidence commit's applicable check/readback state is recorded separately in the PR description.

## 6. Tests of the Tests

Archive identity controls exercise a valid archive plus changed SHA, size, inventory, entry count and same-size byte mutation. The original archive is unchanged by those controls. Parser tests include the actual frozen DESCRIPTION and a trailing-comma fixture against R's parser. The console-streaming canary verifies stdout and stderr are both retained. The manual canary exercises the exact Courier metrics named by the failed checks. Every runner recorded missing makeindex (exit 127) before conditional installation and a working nonempty index afterward. Matrix-oracle controls passed a positive copy and rejected wrong hash, failed status, wrong CRAN mode and wrong command; original evidence stayed unchanged.

## 7a. Issue Ledger

Attempt 1 exposed an empty dependency token before R CMD check. Attempt 2 passed package tests, examples and vignette rebuilding but failed Windows PDF generation due missing Courier metrics; Linux reported missing Tidy/V8. Those failed attempts remain failure evidence. Attempt 3 passed the complete Linux check and Windows HTML/no-index PDF checks, while Windows indexed PDF still failed. Its devel result also recorded a 5.92-second elapsed example timing NOTE, with 2.82 CPU seconds. A short retained index probe precedes attempt 4; no warning is waived. Independently reviewed tool repairs precede each retry.

Reader review also found a stale nested rung, stale panel wording and a one-character duplicate trees300 hash typo. All six current data hashes plus NOTICE were recomputed from the exact archive. The authoritative map was already correct; the duplicate and current wording are corrected, with historical snapshots preserved and annotated.

## 8. Consistency Audit

Compared the current checkout with the original archive: 215 common regular files are byte-identical; DESCRIPTION differs by expected build whitespace normalization. No package, test, help, data or vignette source changed in this slice. Existing R-hub results retain their rebuilt-source qualification; older Win-builder results remain associated with their predecessor artifact.

## 9. What Did Not Go Smoothly

The first parser assumed dependency fields could not contain a trailing empty token. The first full check exposed incomplete manual tools despite TinyTeX setup. The next run isolated an additional index-processing gap; R overwrites the first manual error log during its fallback, so the checker now captures an independent index probe before the package check. These were checker defects and environmental failures, not passing platform results. The original tarball retained its pinned identity throughout.

## 10. Known Residuals

There are 161 testthat warnings whose messages the CRAN-mode reporter suppresses; the counts remain visible and messages are unclassified. The 86 skips comprise 54 explicit CRAN exclusions, 2 unavailable BIEN-cache tests, 3 slow/live integration tests, 24 source-only documentation/helper tests and 3 alternate missing-dependency branches in environments where the dependencies are installed. No claim of full feature coverage or warning-free tests follows.

Current candidate deployment, final post-merge freeze and final release closure remain open. This slice performs no merge, deployment or CRAN submission.

## 11. Team Learning

Check setup actions do not establish tool completeness. Test the actual font metrics and rendering engines, then retain the complete final check verdict. A fallback PDF is not proof that the manual check passed. Trace original compressed bytes into the check process and verify them afterward.

Memory receipt: recalled the pigauto release history before continuing; route.py and the project LOAD-FIRST manifest shaped the exact-artifact and optional-runtime boundaries. The portable validator's regression controls were rerun in the previous slice. Golden Set: no new package behavior or prediction-path defect is changed here; source equivalence and existing artifact controls protect that boundary. No brain-vault edits were made in this slice.

## 12. Cross-Product Coverage

Covers the unchanged frozen gzip's CRAN-mode installation, examples, runnable tests, vignette rebuilding and complete PDF/HTML manual checks on the named runner/R matrix, subject to the retained final verdicts. It does NOT cover skipped tests, broad scientific recovery, all optional backend families, loading or variance-component pooling, future R versions, deployed-site parity, a new post-merge archive, CRAN acceptance or publication.

Writing assessment: separate exact-draft self-assessment is retained beside this report. It distinguishes style judgment from scientific, factual and reference validation.
