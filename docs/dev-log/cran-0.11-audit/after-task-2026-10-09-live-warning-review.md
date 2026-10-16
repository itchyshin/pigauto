# After-task: live homepage warning review

## 1. Goal

Recheck the known homepage warning in the public deployment, review the proposed source correction, and preserve the distinction between local rendering and a live deployment.

## 2. Implemented

Recorded the current Chrome observation and Pat's bounded review in the CRAN 0.11 gate ledger. The live homepage still displays the literal `[!WARNING]` token. Pat's local build of commit `45df106ba58f5bc97d7a274bc0a22a19d202fa17` renders the replacement warning as a blockquote with a bold label and removes the literal token.

## 3a. Decisions and Rejected Alternatives

Keep the deployment fix open until the correction reaches `main`, the Pages workflow succeeds, and the deployed homepage is checked in Chrome. Do not count the local build as live-site evidence. Keep the release evidence PR Draft and unmerged.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- This after-task report.
- No package source, tests, website source, deployment settings, or other lane files changed.

## 5. Checks Run

- Opened `https://itchyshin.github.io/pigauto/` in Chrome. The page shows version 0.11.0 and its warning text begins with the literal `[!WARNING]` marker.
- Checked local Git ancestry: `45df106ba58f5bc97d7a274bc0a22a19d202fa17` is not an ancestor of `origin/main`, whose current recorded commit is `d76804e768bf77f430f43dcf591633f9cd900dab`.
- Verified Pat's recorded exact-head build-log SHA-256: `ba7e2f361fa0c6980bbbbcca9fb342388602968edeee60963450e540f3c49db6`.
- No deployment, source merge, or artifact check was run.
- `check-after-task.R` passed the required report structure check, then withheld overall completion because the release-audit and separate simulation acceptance ledgers still have unmet gates. This review is recorded as a bounded completed slice; the CRAN audit goal remains active.

## 6. Tests of the Tests

No tests were added or changed. Pat's static candidate-site verifier and local pkgdown build are recorded as independent candidate evidence; they do not test the deployed page.

## 7a. Issue Ledger

- Open: the live homepage warning is rendered as raw Markdown syntax.
- Reviewed fix: local source/rendering pass on `codex/pkgdown-warning-callout`, commit `45df106`.
- Required closure evidence: merge to `main`, successful Pages run on that commit, and a fresh Chrome check of the deployed homepage.
- Unchanged: the frozen tarball's Win-builder diagnostics remain failed and are not checksum-bound; exact-artifact closure remains open.

## 8. Consistency Audit

The public page observed in Chrome matches the ledger's open-warning description. The candidate source build is described separately from the live deployment. The source correction is not present in `origin/main`, and the audit PR has no deployment.

## 9. What Did Not Go Smoothly

The warning source fix is complete in a separate branch, but the public site cannot reflect it until that source reaches `main` and Pages deploys it. This audit lane did not change that branch.

## 10. Known Residuals

This check does not establish the visual rendering of the proposed fix in the deployed site, future deployment permanence, Windows behavior, CRAN acceptance, or resolution of the frozen artifact's failed platform diagnostics.

## 11. Team Learning

A correct local HTML build can isolate a display defect, but only a successful deployment and a fresh browser check establish the public result.

## 12. Cross-Product Coverage

This review covers the live homepage display and the recorded local candidate build. It does NOT cover all live routes, the full CSS cascade, an integrated source change, a new deployment, the frozen artifact, or CRAN submission.
