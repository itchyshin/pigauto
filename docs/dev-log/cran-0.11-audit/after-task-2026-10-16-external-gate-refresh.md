# After-task: external gate refresh

## 1. Goal
Refresh the release gate state after PR #234 merged and test whether the previously blocked evidence branch can now be published.

## 2. Implemented
Recorded the latest Actions/deployment/live-site state and the failed evidence-branch push in `GATES.md` and `release-ledger.json`. No package source, website source, or tarball changed.

## 3a. Decisions and Rejected Alternatives
- Kept G7 unmet because GitHub still has not deployed the merged source; the live page demonstrably contains stale text.
- Did not freeze a tarball before the deployed-site gate and exact platform checks can be completed.
- Assumption: none.

## 4. Files Touched
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-external-gate-refresh.md`

## 5. Checks Run
- Chrome: R-CMD-check run `38843648553` still 0/3 jobs; Pages run `38843648852` queued.
- Chrome: Deployments marks `8082ef8` “Ready to deploy (idle)”. The live Getting Started page still says `X_scaled` is 300 × 4 and unmatched species are dropped with a warning.
- Bundled Git remote read: no `codex/cran-011-gate-refresh` ref exists.
- Authorized `git push -u origin codex/cran-011-gate-refresh`: rejected with `remote: fatal: early EOF` and `index-pack failed`; no ref was created.
- `python3 -m json.tool docs/dev-log/cran-0.11-audit/release-ledger.json`: passed.
- Bundled `git diff --check`: passed.
- `lane_preflight.sh` and `route.py` could not recognize this managed worktree as a Git repository/project manifest. The closeout generator also could not write its scaffold because it resolved the path outside the writable worktree. These were recorded; no lane takeover or brain write was attempted.

## 6. Tests of the Tests
No package test or verifier logic changed. The JSON parser and Git whitespace check validate the two ledger edits; the live page and deployment UI directly expose the stale/deployed state.

## 7a. Issue Ledger
- The deployed-site gate in the release ledger remains unmet: successful corrected Pages deployment and live route/discovery verification. The separate scoped Unlazy G7 (local candidate-site verifier) remains met.
- G3 remains unmet: exact frozen-artifact checks on R-release/R-devel and available independent platform.
- G4 remains unmet: independent verdict on the exact artifact and deployed site.
- G5 remains unmet: publish and leave unmerged the final evidence PR.
- Scoped audit count remains 4/7 met, 3/7 unmet.

## 8. Consistency Audit
The merged-source local site receipt remains bound to commit `8082ef8273177ff664d49ac4bd643f7def671f20`; this refresh does not alter it. Current live output still differs from that receipt. No artifact hash or release-readiness claim was advanced.

## 9. What Did Not Go Smoothly
The two GitHub workflows remain queued. A push retry reached GitHub but the remote rejected its pack with early EOF. Local lane and closeout helpers did not resolve this managed worktree correctly.

## 10. Known Residuals
No evidence PR was created. The live site remains stale, platform checks are absent, and no final tarball exists. The cause of GitHub's queue and pack rejection is unknown.

## 11. Team Learning
Treat “Ready to deploy (idle)” as pending deployment, not proof of a live update. Recheck the actual route before closing G7. Repository routing/preflight must handle managed worktrees before relying on those helpers here. Memory receipt: used the repo LOAD-FIRST contract supplied in this task; `route.py` had no manifest for this worktree. Golden Set: not checked; no package behavior changed. No brain files were modified.

## 12. Cross-Product Coverage
Covers ✓: remote branch existence, GitHub workflow/deployment state, the live Getting Started reader surface, and JSON/whitespace integrity for the release ledger.

Does NOT cover: package runtime, a built tarball, CRAN checks, Windows or macOS platform evidence, exact-artifact review, all live site routes/search/sitemap, or CRAN submission.
