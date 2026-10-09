# Exact-head candidate CI: after-task report

## 1. Goal

Verify the candidate-source matrix after ledger commit `21576ac` and preserve its result in the CRAN 0.11 gate record.

## 2. Implemented

Recorded GitHub Actions run #37885180439 against PR #231 head `21576acb089eeb09b57e6d7891efeb872a4ba213`. All three configured jobs succeeded: Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release. The run completed in 16m35s. The ledger states that this is candidate-source evidence and does not close the exact-tarball gate.

## 3a. Decisions and Rejected Alternatives

Kept G7, G8, and G9 open. The green matrix does not establish the deployed site's state, an exact post-merge tarball, force-Suggests results, Windows results, or independent review of that artifact. Did not change PR #231's Draft state or merge it.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-exact-head-ci.md`

## 5. Checks Run

- Memory receipt: `python3 ~/shinichi-brain/tools/route.py pigauto` refreshed the project manifest. The candidate-source versus exact-artifact distinction shaped the decision to keep G8 open despite the green matrix.
- Chrome run #37885180439 showed overall `Success`, three completed matrix jobs, and successful Ubuntu R release, Ubuntu R-devel, and macOS R release jobs.
- The Ubuntu R release and R-devel job pages showed successful completion in 10m48s and 10m47s. The macOS R release job showed successful completion in 16m27s. Run summary duration was 16m35s.
- Git diff from prior tested head `044c1b6274543d39d67cccb73d10f28979c3883b` to `21576acb089eeb09b57e6d7891efeb872a4ba213` contains only the CRAN audit ledger and the preceding ledger-reconciliation report; package source and site inputs are unchanged.
- GitHub reported runner-image migration notices on the Ubuntu jobs and a macOS arm64 capacity notice. The jobs still succeeded.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md`: 11 gates, 8 met, 3 unmet (G7, G8, G9).
- No package tests or site build were rerun in this evidence-recording slice; the matrix itself performed candidate-source checks.

## 6. Tests of the Tests

The run summary explicitly reports overall success and three completed jobs; the individual job views show each platform succeeded. This verifies the matrix completion status, not the contents of the exact post-merge tarball.

Golden Set: No source-code known-mistake class was in scope. The relevant control was the per-platform GitHub Actions status.

## 7a. Issue Ledger

Resolved: the ledger's exact-head candidate-source CI receipt now covers the current PR #231 head `21576ac`.

Open: G7 deployed-site verification after maintainer-controlled merge and deployment; G8 exact post-merge tarball, its inventory, local force-Suggests checks, and platform results; G9 independent review of the final artifact and site evidence.

## 8. Consistency Audit

Confirmed that commit `21576ac` differs from previously checked source head `044c1b6` only in audit-record files. The updated ledger identifies the current run as candidate-source CI and explicitly leaves the exact-artifact gates open. The PR remains Draft and unmerged.

## 9. What Did Not Go Smoothly

The GitHub CLI could not connect to `api.github.com` in this environment. Chrome provided the authoritative run summary and individual job completion evidence. No failed or cancelled jobs were observed.

## 10. Known Residuals

The new ledger entry and report are not yet pushed. The source PR still requires Shinichi's review and merge before deployed-site verification. No final tarball, Windows result, force-Suggests check, exact-artifact verdict, deployment, or CRAN submission exists for this candidate.

## 11. Team Learning

For a source-only docs commit, bind the CI receipt to the exact head and state the unchanged package-input scope. Keep source-matrix success separate from frozen-artifact evidence.

## 12. Cross-Product Coverage

This slice covers candidate-source CI status and the release ledger receipt. It does NOT cover runtime behavior changes, the deployed website, a frozen post-merge tarball, CRAN-style force-Suggests checks, Windows, exact-artifact independent review, or CRAN acceptance.
