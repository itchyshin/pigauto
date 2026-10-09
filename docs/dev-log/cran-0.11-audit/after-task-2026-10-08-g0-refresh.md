## 1. Goal

Refresh the CRAN 0.11 release ledger against the current source PR head, confirm the active lane count, and add a repeatable fresh-site build check.

## 2. Implemented

Updated G0 with current PR #231 head and head-specific CI state, current PR #228 status, the publication-history evidence path, and the one-active-lane preflight result. Added a source-pinned site-build script and recorded its fresh build, retirement, link, and pkgdown checks under G5/G5b. G0 remains open.

## 3a. Decisions and Rejected Alternatives

Kept 0.11.0 as the candidate because the current public CRAN record is 0.10.0. Did not treat the public record as proof that no private or pending submission exists. Did not treat checks attached to older source commit `7cbbcbc` as validation of the current `f0f6388` head.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-g0-refresh.md`
- `script/cran-0.11-site/build-release-site.sh`

## 5. Checks Run

- `bash ~/shinichi-brain/tools/lane_preflight.sh /Users/z3437171/.codex/worktrees/cran-011-pr231-verified/pigauto`: one active lane; 58 total worktrees.
- Chrome: inspected PR #231 head and its specific checks page; R release, R devel, and macOS jobs all passed in run #37864204151 (22m49s); pkgdown PR run was skipped by workflow design.
- Chrome: inspected PR #228; it remains Draft with 28 commits, no reviews, and no deployment. Its body refers to an older PR #231 commit.
- Chrome cache-busted recheck: PR #231 now has 20 commits and points to pushed commit `977cca41c29018ed88bac96bc3e444fd395cca66`; its checks were queued when inspected. The previous head's three platform checks passed.
- Unlazy status-only check: 11 gates, 6 met and 5 unmet (G0, G6, G7, G8, G9). This parsed checkboxes and ran no gate commands.
- `git diff --check`: passed.
- `CHECK_AFTER_TASK_ACTIVE=1 Rscript ~/shinichi-brain/tools/check-after-task.R "$PWD/docs/dev-log/cran-0.11-audit/after-task-2026-10-08-g0-refresh.md"`: after-task structure check passed.
- `python3 ~/shinichi-brain/tools/slop_check.py /Users/z3437171/.codex/worktrees/cran-011-pr231-verified/pigauto/docs/dev-log/cran-0.11-audit/after-task-2026-10-08-g0-refresh.md`: 0 hits / 776 words; no em dashes.
- `CHECK_AFTER_TASK_ACTIVE=1 Rscript ~/shinichi-brain/tools/check-after-task.R "$PWD/docs/dev-log/cran-0.11-audit/after-task-2026-10-08-g0-refresh.md"`: after-task structure check passed. `python3 ~/shinichi-brain/tools/slop_check.py "$PWD/docs/dev-log/cran-0.11-audit/after-task-2026-10-08-g0-refresh.md"`: 0 hits / 840 words, no em dashes.
- `bash -n script/cran-0.11-site/build-release-site.sh`: passed.
- `bash script/cran-0.11-site/build-release-site.sh`: passed on source `f0f6388dafa4d9f8840fed19970439c6f33dbe83`; 62 HTML pages, 3,443 local references, 34 retired routes, 613 search entries, zero crawler errors, and `pkgdown::check_pkgdown()` reported no problems. Build-log SHA-256 `e72daf8b93d5ab2a223079753bf905754ffdec15f46d204e792282758e107e0a`; crawl-log SHA-256 `2981f54a2056d62a71e37ea4f85eb8a8d6e22f242fc8e2d47d54daa1ab560ece`.
- Reran the committed site-build script on pushed `HEAD` `977cca41c29018ed88bac96bc3e444fd395cca66`: 62 HTML pages, 3,443 local references, 34 retired routes, 613 search entries, zero crawl errors, and no pkgdown problems. Build-log SHA-256 `7f35ddcc34e563995882d72a75fe9eabff6e1b9a762182a1a7a1c11516989b0f`; crawl-log SHA-256 `2981f54a2056d62a71e37ea4f85eb8a8d6e22f242fc8e2d47d54daa1ab560ece`.
- Local commit `977cca41c29018ed88bac96bc3e444fd395cca66` was pushed to the authorized draft PR #231 branch; Chrome confirmed the PR head and Draft state.
- Independent read-only review of the site-build chain found no blocker. It verified that all 34 manifest routes and four simulation-study variants are absent from generated files, search, and sitemap. It noted the crawler covers HTML `href`/`src`, not CSS, `srcset`, or JavaScript-generated routes, and the default SRI-checked asset cache path is host-specific.
- Normalized rendered-content assertions passed for BirdTree retrieval and citations, Newick/NEXUS guidance, inverse-Wishart historical-use limits, and NEWS missingness results. A first raw-HTML substring check missed text split by an anchor; the normalized text check passed.
- Official BirdTree Downloads and Subsets pages and current CRAN Repository Policy were checked. BirdTree requires attribution for research use; the reviewed pages do not state a redistribution grant. CRAN requires clear rights for all package components, including data.

## 6. Tests of the Tests

This change records current-state evidence and does not change package behavior. The Unlazy `--status` command was used only to verify the ledger parse and tally; it was not treated as a rerun of any gate.

## 7a. Issue Ledger

- G0 remains open: current public history supports the 0.11.0 candidate, but pending-submission status and the redistribution basis for bundled BirdTree data remain unresolved.
- G6, G7, G8, and G9 remain open as recorded in the release ledger.
- PR #231 checks passed on commit `f0f6388`; the exact-head run for `977cca4` is in progress.
- The 3,441-reference site count in an earlier update is superseded by the pinned 3,443-reference reruns.
- The PR #228 body has a stale source-head summary and needs reconciliation before it is used as the current release snapshot.

## 8. Consistency Audit

Compared the current branch commit with GitHub's PR commit list and the head-specific checks page. Compared the publication-history record with the CRAN and GitHub release/tag pages already recorded there. Compared the release ledger tally with its checkbox states. The source PR remains unmerged and undeployed, and the release artifact described by PR #228 predates the current source candidate.

## 9. What Did Not Go Smoothly

The first lane-preflight path did not exist in this checkout; the hub script produced the valid one-lane result. An initial PR page load showed stale 13-commit state, so I loaded a cache-busted page; it confirmed 19 commits at `f0f6388`. After pushing, a fresh page confirmed the PR at 20 commits and `977cca4`. Shell `git ls-remote` could not reach GitHub over SSH/DNS, so Chrome provided the current remote branch and PR state. The closeout compiler's structure check passed, but its integrated acceptance check found five open gates under `.unlazy/imputation-sim`; that separate goal card says the simulation arc is complete. I left its ledger untouched. The pigauto CRAN ledger is `docs/dev-log/cran-0.11-audit/GATES.md`; the `OWNS` list previously named a nonexistent `.unlazy/cran-0.11-audit/GATES.md`, which I removed. The repo-local Unlazy status check parses the actual CRAN ledger.

## 10. Known Residuals

No exact final tarball is frozen. Local visual review and deployed-site verification are incomplete. CRAN redistribution rights for the bundled BirdTree objects remain unverified. PR #231 checks for `977cca4` and a fresh PR #228 summary are pending. No merge, deployment, or CRAN submission occurred.

## 11. Team Learning

Bind every CI statement to the exact commit's checks page. A PR-level checks page can retain a previously selected commit and mislead a status refresh. Memory receipt: loaded the pigauto operating contract and used its one-lane preflight, exact-artifact, and unlazy gate rules. Brain search for a durable pigauto decision returned no result; the repository ledger was used as technical truth. Golden Set: not in scope for this evidence-ledger refresh.

## 12. Cross-Product Coverage

Covers the release gate ledger, current source-PR status, publication-history record, and lane census. Does NOT cover the final merged source, deployment, frozen artifact, platform checks, exact-artifact independent reviews, or CRAN submission.

Assessment: 2/10, moderate confidence. Current GitHub state and lane count were checked directly. The site-build chain received an independent review; exact-head CI is still underway. BirdTree and CRAN policy statements were checked against their official pages. No scientific claim is made.
