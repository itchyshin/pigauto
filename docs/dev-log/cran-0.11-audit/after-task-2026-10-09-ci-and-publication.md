# Candidate-source CI and publication refresh: after-task report

## 1. Goal

Record final candidate-source CI and refresh the public release-history check for PR #231 while keeping the remaining CRAN gates explicit.

## 2. Implemented

Added the completed CI matrix and current publication-history receipt to `docs/dev-log/cran-0.11-audit/GATES.md`. This is evidence documentation only; no package source, tests, website source, or generated help changed.

## 3a. Decisions and Rejected Alternatives

Kept 0.11.0 because the live CRAN record and the latest GitHub release both remain at 0.10.0. Treated the green GitHub matrix as candidate-source evidence only because it used `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`, and ran on GitHub's synthetic PR merge commit. Did not create or label a candidate tarball as the final artifact.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-ci-and-publication.md`

## 5. Checks Run

- GitHub Actions run [#37876149667](https://github.com/itchyshin/pigauto/actions/runs/37876149667) completed successfully for PR #231 head `3f45b0a9f74f3db616046d36844a7810e335f364`, tested through merge commit `101a6308f6ff4639a8faf9cd45911374b28475d1`.
- Ubuntu R release (R 4.6.1), Ubuntu R-devel (R Under development, 2026-10-06 r90643), and macOS arm64 R release (R 4.6.1) each reported `R CMD check` `Status: OK`.
- Each full test run reported 3,196 passes, 0 failures, 175 warnings, and 83 skips. The macOS focused MPS prediction tests reported 201 passes, 0 failures, 50 warnings, and 0 skips.
- The run used `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`. The pkgdown pull-request workflow was skipped by repository design.
- Chrome showed CRAN 0.10.0 and GitHub latest release v0.10.0. The local worktree was clean at the start of this evidence update; lane preflight reported one pigauto lane.
- The after-task structure check passed. The complete acceptance scan exited 1 because five pre-existing `.unlazy/imputation-sim/gates/leaf-*.md` gates remain unmet. They belong to a separate campaign record and were not changed or reclassified in this evidence slice.

## 6. Tests of the Tests

No test code changed in this slice, so I did not run a negative control. The three platform jobs exercised the existing full package test suite; the focused macOS MPS test also completed.

## 7a. Issue Ledger

Resolved: the candidate-source CI result and public release-version check are now recorded against exact source and merge SHAs.

Open: G0 rights and provenance, G7 deployed-site verification, G8 the exact post-merge tarball and its checks, and G9 independent review of final artifact evidence. The gate tally remains 7 met and 4 open.

## 8. Consistency Audit

Compared the live CRAN version and GitHub latest release with `DESCRIPTION` version 0.11.0. Confirmed the later commit after the exact-head local site build changed only audit files under `docs/dev-log/cran-0.11-audit/`, so it did not alter website inputs. The public record cannot establish a private or pending CRAN submission.

## 9. What Did Not Go Smoothly

The macOS job took about 25 minutes and remained active longer than prior matrix runs, but it eventually completed successfully. While it was active, the logs endpoint returned a temporary missing-blob response. After completion, the full job logs were available and their summaries were checked.

The brain repository's `closeout.py new` resolves relative paths against the repository containing the script. It created a template outside pigauto. I left that file untouched under the brain-write boundary and used the explicit pigauto worktree path for this report.

## 10. Known Residuals

These are candidate-source checks, not checks of one frozen tarball. The run did not force all Suggests packages to install and does not replace the planned exact-artifact checks. No Windows result, site deployment, live retired-route verification, final tarball, or exact-artifact independent verdict is available. The CRAN redistribution-warranty basis remains unresolved pending the maintainer's clarification.

## 11. Team Learning

Memory receipt: loaded the pigauto `route.py` manifest and checked D-318, which keeps the CRAN gate after the bug-fix work closes. That shaped the separation between source readiness and publication. Golden Set: not checked because this evidence-only slice changed no code and did not touch a known-mistake class.

## 12. Cross-Product Coverage

Covers candidate-source Linux and macOS release checks, the R-devel check, the focused MPS test, and current CRAN/GitHub public version records.

Does NOT cover the frozen source tarball, force-Suggests artifact checks, Windows, deployed website state, redistribution-rights evidence, final independent artifact review, or CRAN acceptance.
