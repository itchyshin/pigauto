## 1. Goal

Verify whether pigauto 0.11.0 is already published on CRAN and add direct evidence to G0.

## 2. Implemented

Added the direct 0.11.0 tarball URL result to the CRAN 0.11 gate ledger. The evidence supports 0.11.0 as unpublished at the time checked.

## 3a. Decisions and Rejected Alternatives

Used the official CRAN package page and direct candidate tarball URL in Chrome. A missing archive directory alone would not establish publication history, so the check included the exact `pigauto_0.11.0.tar.gz` URL. No version was changed.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-cran-version-check.md`

## 5. Checks Run

- Chrome opened `https://cran.r-project.org/web/packages/pigauto/index.html`. It lists version 0.10.0, published 2026-07-30.
- Chrome opened `https://cran.r-project.org/src/contrib/pigauto_0.11.0.tar.gz`. The page displayed CRAN's `Error 404` response.
- `rg -n '^Version:' DESCRIPTION` on the candidate evidence worktree reported `Version: 0.11.0`.
- `git diff --check` passed.
- Lane preflight confirmed one active pigauto lane and reported the GATES file is not being worked on by a missing ref.

## 6. Tests of the Tests

No package tests changed. The direct URL check tested the exact source archive name for the candidate version, alongside the current official package record.

## 7a. Issue Ledger

- Resolved: the available CRAN page lists 0.10.0, and the direct 0.11.0 source tarball URL returned an error page.
- Open: full archived release history, BirdTree redistribution rights, final merged source identity, deployment, and exact-artifact platform checks.

## 8. Consistency Audit

The live CRAN record, direct 0.11.0 archive URL, and candidate `DESCRIPTION` version agree with selecting 0.11.0 as the candidate version. The 404 is a point-in-time check and is not evidence about the complete CRAN archive history or the rights of bundled data.

## 9. What Did Not Go Smoothly

The GitHub CLI could not reach its API in this environment. GitHub state was not needed for this version check; PR state remains sourced from the Chrome check recorded in the active gate evidence.

## 10. Known Residuals

This check does not establish the complete publication archive, data redistribution rights, release readiness, or validity of any source tarball. The overall G0 gate remains open.

## 11. Team Learning

For a candidate package version, check the current CRAN record and the exact expected tarball URL. Keep a 404 result bounded to that URL and time.

## 12. Cross-Product Coverage

Covers the current CRAN package record, the exact 0.11.0 source tarball URL, and the candidate `DESCRIPTION` version. This check does NOT cover complete archive history, BirdTree data rights, deployment, or release-artifact validation.
