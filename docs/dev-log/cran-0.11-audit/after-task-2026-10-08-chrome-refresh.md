## 1. Goal

Refresh the pigauto 0.11 release and Unlazy ledgers with current Chrome observations of PRs #228 and #231 and the independent local-site review.

## 2. Implemented

Added dated receipts to the release gate and local Unlazy gate. PR #228 is recorded as Draft with seven commits, no review or deployment, and three merge conflicts. PR #231 is recorded as Draft with nine commits, three successful checks, one skipped check, no review, conflicts, or deployment. The local site reviewer’s crawler result is recorded with its limits. A fresh official BirdTree page check supports research use with attribution while leaving package redistribution unresolved.

## 3a. Decisions and Rejected Alternatives

Used Chrome as requested. Kept the candidate tarball and Windows results explicitly unbound where the visible receipts lack checksum evidence. Did not resolve PR conflicts, merge either PR, deploy the site, or submit to CRAN.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `.unlazy/cran-0.11-audit/GATES.md` (ignored local ledger)
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-chrome-refresh.md`

## 5. Checks Run

- Lane preflight and a scoped lease covered the two ledgers and this report.
- `git diff --check` passed.
- `python3 -m json.tool docs/dev-log/cran-0.11-audit/release-ledger.json` parsed successfully; it still records empty artifact and evidence objects and `audit_verdict: NOT_READY`.
- Unlazy status parsing reports 5 met and 6 unmet release gates, and 12 met with I7 unmet in the detailed implementation ledger.
- Chrome showed PR #228 conflicts in `pre-pr-check.md`, `rights-and-policy.md`, and `surface-inventory.md`. PR #231 shows three successful checks and one skipped.
- The independent site review reran the local crawler: 62 HTML pages, 3,548 local references, 34 retired routes, zero errors; a planted Markdown route failed as expected.
- Chrome checked BirdTree’s [downloads](https://birdtree.org/downloads/), [FAQ](https://birdtree.org/faq/), and [methods](https://birdtree.org/methods/) pages. The downloads page requires Jetz et al. (2012) for research use and BirdTree.org for use of its web tool. The pages checked do not state a licence or package redistribution grant.
- Chrome rechecked the official [CRAN pigauto record](https://cran.r-project.org/web/packages/pigauto/index.html) on 2026-10-08. It still lists version 0.10.0, published 2026-07-30; the candidate `DESCRIPTION` declares 0.11.0. This check confirms the current listing; the complete CRAN archive history remains unverified.
- Chrome opened the official [CRAN per-package archive URL](https://cran.r-project.org/src/contrib/Archive/pigauto/), which returned CRAN’s 404 page. The parent archive index loaded, but Chrome accessibility timed out before its contents could be checked; opening the source `PACKAGES` index returned `net::ERR_BLOCKED_BY_CLIENT`. These results leave the full archived version history unverified.
- A read-only Gmail search for `in:anywhere after:2026/10/07 "win-builder"` returned no messages, so no new Win-builder result email is available to bind to the current candidate.
- Chrome refreshed PRs [#228](https://github.com/itchyshin/pigauto/pull/228) and [#231](https://github.com/itchyshin/pigauto/pull/231). #228 remains Draft with seven commits, three conflicts, no reviews or deployment; its conversation records macOS checks for hashes `d34d4699…` and `10573d4f…`, plus a separate Windows upload hash `5fbbec8c…` whose results are pending. #231 remains Draft with nine commits, three successful checks, one skipped, no conflicts, no reviews, and no deployment. The hashes cannot be treated as one platform-verified artifact.
- Rechecked `/private/tmp/pigauto_0.11.0.tar.gz`: SHA-256 `1bbc2a7b45395c5c4dd9ef7769cab95e3e56e90535020b2ab10f123d74ec6c30`, 5,050,355 bytes, modified 2026-10-08 01:43:38 MDT. The release ledger has no source or receipt binding for this loose archive, so it does not count as release evidence.
- An independent read-only artifact-gate review found no factual mismatch in the refreshed claims and flagged this loose archive as a separate unbound file. The review did not approve the exact release artifact or close G9.

## 6. Tests of the Tests

The planted retired Markdown-route control was rejected by the site crawler, followed by a clean rerun. `git diff --check` and JSON parsing detect formatting and serialization errors in the edited records; neither validates package behavior.

## 7a. Issue Ledger

- PR #228: three merge conflicts; draft, no reviews or deployment.
- PR #231: draft; one skipped check; no review or deployment.
- Release ledger: `NOT_READY`; no artifact or evidence recorded.
- BirdTree redistribution evidence, local visual inspection, post-correction deployment checks, checksum-bound Windows results, final artifact checks, and independent exact-artifact review remain open.

## 8. Consistency Audit

Compared both ledger addenda with the current Chrome pages, the official BirdTree pages, and the independent site reviewer’s retained crawler run. Claims are scoped to those sources. Research-use citation instructions are distinguished from redistribution permission. Local crawler success is not presented as visual or deployment verification.

## 9. What Did Not Go Smoothly

The first lease attempts could not write to the brain’s registry under the workspace sandbox. A narrow permission request allowed the required scoped lease. The brain closeout generator resolves its repository root to the brain vault even when invoked from this worktree, so it was not used to create this repository report; the generated path would have pointed outside the project.

## 10. Known Residuals

PR #228 remains conflicted and both PRs remain drafts. The 10573… candidate is described in PR #228 and has a passing macOS check, but its Windows results are not checksum-bound and the release ledger is still empty. The locally retired route has not been verified on a deployment after the correction. BirdTree pages support attributed research use but do not establish redistribution rights for the bundled tree objects. Rights and visual-review gates are open.

## 11. Team Learning

For release receipts, record PR state, artifact identity, platform binding, and site evidence as separate facts. A green local crawler and a green source check do not establish deployment or exact-artifact status.

Memory receipt: loaded the pigauto LOAD-FIRST manifest from `route.py pigauto`, the repository instructions, and the pigauto CRAN 0.11 memory entry. The routed prediction-path and uncertainty guards did not change this documentation-only refresh. Golden Set: not in scope for this ledger update.

## 12. Cross-Product Coverage

Covers: current GitHub PR state, local site-crawler evidence, and release-ledger status.

Does NOT cover: package source behavior, a frozen post-merge tarball, checksum-bound Windows results, visual review, deployed-page verification, BirdTree redistribution rights, CRAN submission, or scientific validity beyond the audited workflows.
