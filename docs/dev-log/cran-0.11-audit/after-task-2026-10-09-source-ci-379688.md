## 1. Goal

Record the completed candidate-source CI run for pigauto PR #231 and update the release evidence without claiming exact-artifact or deployment readiness.

## 2. Implemented

Verified GitHub Actions run #37968809743 on PR #231 head 3702d5d533e101dd3583b30703200fe525d1db4f. All three platform jobs completed successfully. Each full R CMD check reported Status: OK and testthat reported 3,196 passes, 0 failures, 175 warnings, and 83 skips. The focused macOS MPS test reported 201 passes, 0 failures, 50 warnings, and 0 skips. Archived the three job logs in one deterministic gzip receipt and verified byte-for-byte decompression against the fetched job logs. Updated the release gate ledger with the completed result and exact source SHA. PR #231's description was verified against the completed run and exact head. PR #228's Draft description was refreshed in the Codex browser with the completed result and verified rendered text.

## 3a. Decisions and Rejected Alternatives

Recorded this as candidate-source CI only. The run used NOT_CRAN=true and _R_CHECK_FORCE_SUGGESTS_=false, so it does not satisfy the force-Suggests exact-artifact gate. Kept G7, G8, and G9 open. Preserved the full logs instead of only recording a green status. No package tests were rerun because the tested source delta since the preceding passing head changed only the audit ledger and its after-task record.

## 4. Files Touched

- docs/dev-log/cran-0.11-audit/GATES.md
- docs/dev-log/cran-0.11-audit/provenance/source-ci-37968809743.log.gz
- docs/dev-log/cran-0.11-audit/after-task-2026-10-09-source-ci-379688.md

The external PR description refreshes are recorded here as browser-verified metadata changes; no repository source or website input changed.

## 5. Checks Run

- Browser check of Actions run #37968809743: completed successfully in 18m02s; all three matrix jobs succeeded.
- GitHub workflow-job records: all three jobs completed successfully; the R CMD check action step succeeded on each platform.
- Downloaded and inspected all three job logs. Each reported Status: OK and [ FAIL 0 | WARN 175 | SKIP 83 | PASS 3196 ]. The macOS focused MPS test reported [ FAIL 0 | WARN 50 | SKIP 0 | PASS 201 ].
- Confirmed NOT_CRAN=true and _R_CHECK_FORCE_SUGGESTS_=false in all check environments.
- Confirmed source diff from cff6ea685d2788687d085e24c2a8ef83ce50f7bb to 3702d5d533e101dd3583b30703200fe525d1db4f changes only GATES.md and the prior after-task record.
- Original raw job log sizes and SHA-256 values: Ubuntu release, 1,048,748 bytes, f84c55b1cf04dfd7f635ab23c0e505433c6db459bbe90d2872327081f0bfed75; Ubuntu devel, 1,068,934 bytes, 0cd0784f594a9eb10cf81f7b9b28b3b598f62bc40c481cb42b23e386d0562f04; macOS arm64 release, 955,330 bytes, c22ed0b0e7a73ed90332d50f398271e5244ced4208bf5b2a3b1ed819ece1c293.
- Combined raw job logs: 3,073,186 bytes, SHA-256 ba14e68809ff5d206ca37f5bddb3aea0905fd44fc3f451183d4da40db04aab6f. Deterministic gzip: 409,760 bytes, SHA-256 291a3c2b3cadd07b2e854a528a54934a03ac14080c1a7808f990e4053d78a4e0. Decompression matched the combined source bytes exactly.
- git diff --check passed for the source-CI-only delta.
- The after-task structural check passed. The full closeout acceptance check remains red on five existing leaves in .unlazy/imputation-sim/gates; this slice owns none of them.
- slop_check.py reported 0 findings for this report.

## 6. Tests of the Tests

The existing full test suite ran on Ubuntu R-release, Ubuntu R-devel, and macOS arm64 R-release. All three results showed zero failures. This receipt did not add or mutate tests, so no mutation test was run.

## 7a. Issue Ledger

Resolved in the ledger: the run #37968809743 result and correct full source SHA are recorded. A reviewer caught two older paragraphs that still read as current while naming PR #231 head dee8f028; both were explicitly relabelled as point-in-time snapshots, preserving their historical run #379525 binding. The current PR #231 description was verified in the browser, and the stale PR #228 description was refreshed and checked for the correct full SHA, completed run, and open release gates. Deferred: exact tarball checks, deployment verification, and independent review of the final artifact.

## 8. Consistency Audit

The Actions run, three job records, and raw logs agree on the successful result and test counts. The tested commit is the current PR #231 head. Its delta from cff6ea consists only of audit files, so the source result remains a candidate-source check and does not imply a behavior change or deployed-site update. PR #231 and PR #228 remain Draft and unmerged. PR #231's body shows the correct head and run; PR #228's refreshed body was verified rendered in the Codex browser. Both descriptions preserve the candidate-source versus exact-artifact distinction.

## 9. What Did Not Go Smoothly

The first log-transfer attempt used a terminal input path that could not safely carry large logs. It was stopped, and the logs were transferred in bounded chunks with per-job hashes and a verified gzip round trip. The GitHub blob and contents API endpoints returned 403, so the PR description was refreshed through the Codex browser. The original ledger, report, and compressed logs were committed as `71038b6` and pushed to `release/cran-0.11-gate`. The post-push browser check confirmed PR #228 remains Draft, its commit list includes `71038b6`, the updated description is rendered, and the branch has no deployment.

## 10. Known Residuals

G7 remains open because PR #231 is not merged or deployed and the live retired route has not been verified after deployment. G8 remains open because no final post-merge tarball has been frozen and checked with force-Suggests enabled on that exact artifact. G9 remains open because reviewers have not approved that final artifact. This source run does NOT establish Windows results, final-artifact results, or CRAN readiness.

## 11. Team Learning

A source-CI receipt must bind the full run result to the full tested commit SHA and retain enough log provenance for later review. The distinction between candidate-source checks and exact-artifact checks shaped this update.

Memory receipt: the pigauto LOAD-FIRST manifest and release-audit record were consulted; they require exact-artifact evidence and prohibit upgrading candidate-source CI into CRAN readiness.

Golden Set: no package source or runtime behavior changed, so no source known-mistake class was in scope.

## 12. Cross-Product Coverage

Covers: Ubuntu R-release, Ubuntu R-devel, macOS arm64 R-release, full testthat output, and focused macOS MPS tests for the named PR head.

Does NOT cover: Windows checks, force-Suggests exact-artifact checks, deployed routes after merge, the final tarball, or CRAN submission.
