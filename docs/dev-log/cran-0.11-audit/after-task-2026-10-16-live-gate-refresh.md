# After-task: live release-gate refresh

## 1. Goal
Refresh the authoritative external state for the CRAN 0.11 audit, confirm the audit branch is published, and preserve the live-site and CI status without claiming release readiness.

## 2. Implemented
Updated `GATES.md` and `release-ledger.json` with the successful branch push, current post-merge job states, current live-page discrepancy, and the still-unmet final evidence PR state. No package code, website source, or tarball was changed.

## 3a. Decisions and Rejected Alternatives
- Kept the final evidence PR unopened until the exact tarball, platform results, and final independent reviews are available.
- Kept the deployed-site and artifact gates unmet because the Pages workflow is queued and the live page still shows stale text.
- Treated lane preflight's named second branch as a merged PR #234 source branch after verifying its commit is an ancestor of `origin/main`; the branch is not an unmerged source change.
- Assumption: none.

## 4. Files Touched
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-live-gate-refresh.md`

## 5. Checks Run
- Bundled Git: current audit branch is clean except for pre-existing untracked `graft/` and `script/cran-0.11-site/__pycache__/`; remote ref resolves to `26048be2405960318b164a3f3912bc01ebd77c05`.
- Bundled Git: `codex/cran-011-getting-started-correction` at `6ebc676` is an ancestor of `origin/main` at `8082ef8`; PR #234 is merged.
- GitHub CLI: R-CMD-check run `38843648553` is queued overall; macOS R-release is running `check-r-package`, while Ubuntu R-release and R-devel remain queued. GitHub reports logs are unavailable until the live job finishes.
- GitHub CLI: Pages run `38843648852` remains queued; no open PR exists for `codex/cran-011-gate-refresh`.
- Chrome: refreshed live Getting Started page still shows `X_scaled (300 × 4)` and the obsolete unmatched-species warning.
- Chrome: GitHub Status reports Actions and Pages operational; the queue cause is unknown.
- `python3 -m json.tool docs/dev-log/cran-0.11-audit/release-ledger.json`: passed.
- Bundled `git diff --check`: passed.

## 6. Tests of the Tests
No package or verifier logic changed. JSON parsing and Git whitespace checks validate the ledger update. The live page was checked in Chrome, and job state came from the authoritative GitHub run records.

## 7a. Issue Ledger
- The final evidence branch is now pushed, resolving the earlier transport blocker. The final evidence PR is still not open.
- The deployed-site gate remains unmet: the Pages run is queued and live content remains stale.
- G3 remains unmet: no final frozen tarball or checksum-bound R-release/R-devel results exist.
- G4 remains unmet: final artifact and deployed-site review remain outstanding.
- G5 remains unmet: the unmerged final evidence PR awaits the artifact and review evidence.
- Scoped audit count remains 4/7 met and 3/7 unmet.

## 8. Consistency Audit
The exact-main local site receipt remains bound to `8082ef8273177ff664d49ac4bd643f7def671f20`. The current live page differs from that local output. The successful branch push does not change package source, artifact identity, or release readiness. Highest proven rung remains `source-clean`; next is `tarball-clean`.

## 9. What Did Not Go Smoothly
Ubuntu release/devel and Pages remain queued despite GitHub reporting the services operational. The macOS check is active but has not exposed logs yet. Lane preflight reported the already-merged PR #234 source branch as a second recent Codex lane.

## 10. Known Residuals
No final tarball or final evidence PR exists. The live site remains stale, Ubuntu checks and Pages deployment have not completed, and exact-artifact checks and final review are outstanding. The cause of the workflow queue is unknown.

## 11. Team Learning
The prior branch-push failure was transport-specific: the verified SSH-over-port-443 retry succeeded, and the pack for `origin/main..HEAD` was 119 KB. Treat a merged branch still checked out in another worktree as a preflight signal to verify, not as proof of an active conflicting edit. Golden Set: not checked; no package behavior changed. No brain files were modified.

## 12. Cross-Product Coverage
Covers ✓: audit branch publication, current Actions and Pages state, lane-preflight branch verification, and the live Getting Started discrepancy.

Does NOT cover: package runtime, final tarball, CRAN checks, final artifact platform evidence, all live site routes/search/sitemap, exact-artifact independent review, final evidence PR, or CRAN submission.
