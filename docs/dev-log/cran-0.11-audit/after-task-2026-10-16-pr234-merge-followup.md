# After-task: merge and post-merge state check for PR #234

## 1. Goal

Record the user-authorized merge of the Getting Started correction and distinguish merged source from the deployed website and settled checks.

## 2. Implemented

- PR #234 merged to `main` as `8082ef8273177ff664d49ac4bd643f7def671f20`.
- Updated the release ledger and gate narrative with the merge SHA, main workflow run IDs, and observed live-site state.
- Created a clean detached checkout of the merge commit at `/private/tmp/pigauto-0.11-merged-8082ef` for post-deployment artifact work.

## 3a. Decisions and Rejected Alternatives

- Followed the user's explicit instruction to merge PR #234 despite queued non-required checks; did not bypass any repository protection.
- Kept G7 open because the live Getting Started page still showed the pre-correction `X_scaled` 300 × 4 claim.
- Did not freeze a tarball before the deployment state is verified.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-pr234-merge-followup.md`

## 5. Checks Run

- GitHub PR page confirmed merged status and merge commit `8082ef8`.
- The PR showed 1 of 7 checks passed at merge time; three R checks were queued and pkgdown PR check was skipped by its workflow guard.
- Commit checks confirmed main R-CMD-check run `38843648553` and pkgdown run `38843648852`, both queued at observation.
- GitHub Status showed Actions operational.
- Chrome rendered the deployed Getting Started page and still showed the old `X_scaled` 300 × 4 statement; this confirms G7 is unmet at that observation.
- `release-ledger.json` parses as valid JSON.

## 6. Tests of the Tests

No test code changed. Merge state, commit SHA, workflow IDs, and live page contents were checked directly in GitHub and Chrome.

## 7a. Issue Ledger

- Resolved: source PR #234 merged.
- Open: main checks and Pages deployment; current live site remains stale; G7, G8, and G9 remain unmet.

## 8. Consistency Audit

The merge event, release ledger, GitHub commit checks, and live-site content agree that source is merged but deployment has not yet reached the corrected page.

## 9. What Did Not Go Smoothly

GitHub's R-CMD-check and pkgdown jobs remained queued after the merge. The run queue is the current external dependency.

## 10. Known Residuals

No deployed-site verification, final post-merge tarball, exact-hash platform results, or final artifact review exists yet. CRAN submission remains unauthorized.

## 11. Team Learning

A successful merge and a main-branch check entry do not establish that the documentation deployment completed. Keep source, deployment, and artifact states separately recorded. No brain-vault notes were edited.

## 12. Cross-Product Coverage

Covers: merge state, post-merge run IDs, and one deployed Getting Started page observation.

Does NOT cover: settled main checks, deployed retired/search/sitemap routes, the exact post-merge tarball, platform checks, or CRAN submission.
