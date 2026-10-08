# After-task report: site evidence reconciliation

## 1. Goal

Make the CRAN 0.11 evidence PR independently verify the source PR's retired-route control and distinguish the stale failing site snapshot from the fresh passing build.

## 2. Implemented

Added the retired-Markdown-route negative-control receipt to the evidence checkout, where `GATES.md` already referenced it. Corrected the 3,548-reference snapshot entry to state that the current source crawler finds `articles/simulation-study.md` and exits 1. The fresh cleaned site at PR #231 source commit `2f7305c7ce4d13e03f676f4e48818f3a259d5473` remains the valid local result: 62 HTML pages, 3,443 local references, 34 retired routes, and zero errors.

## 3a. Decisions and Rejected Alternatives

Kept the route-checker implementation in source PR #231 and copied only its reproducible negative-control output into evidence PR #228. The older site snapshot was not counted as passing evidence. No package code, site source, or public deployment was changed.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/site-crawler-simulation-retirement-control.txt`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-site-evidence-reconciliation.md`

## 5. Checks Run

- `lane_preflight.sh . --file ...` reported one pigauto lane. The primary checkout's dirty GNN-attribution files were left untouched.
- `route.py pigauto` loaded the pigauto LOAD-FIRST manifest. A brain search confirmed the recorded CRAN rights issue remains unresolved.
- The current PR #231 crawler returned exit 1 for `/private/tmp/pigauto-cran011-site-20261008` because it contains `articles/simulation-study.md`.
- The same crawler returned `SITE_CRAWL_OK` for `/private/tmp/pigauto-cran011-followup/_site`: 62 HTML pages, 3,443 local references, 34 retired routes, and zero errors. That build's source checkout is `2f7305c7ce4d13e03f676f4e48818f3a259d5473`.
- Chrome rechecked [PR #231](https://github.com/itchyshin/pigauto/pull/231): Draft, 14 commits, three of four checks successful with the pkgdown PR job skipped, no review, and no deployment. [PR #228](https://github.com/itchyshin/pigauto/pull/228) remains Draft and unmerged.
- `git diff --check` passed. No package tests ran because this slice changes evidence records only.

## 6. Tests of the Tests

The retained negative-control output shows the crawler rejecting an injected retired Markdown route with exit 1 and the exact `Retired historical page served` error. The older saved site independently reproduces the same rejection. The fresh cleaned site passes the same crawler.

## 7a. Issue Ledger

G5 and G5b remain verified for the fresh local build at PR #231 commit `2f7305c`. The old 3,548-reference snapshot is superseded and fails the retired-route check. G0 remains open for BirdTree redistribution rights and final source provenance; G6 remains open for local visual review; G7 remains open for post-deployment route checks; G8 and G9 remain open for a final checksum-bound artifact and independent exact-artifact review.

## 8. Consistency Audit

The evidence PR now contains the negative-control receipt cited by its gate ledger. The source PR retains the checker and cleanup implementation. The 3,443-reference fresh-site count agrees across the crawler output, after-task record, and current gate entry. The task remains the sole pigauto lane; other checkouts were treated as worktrees, not separate lanes.

## 9. What Did Not Go Smoothly

The first comparison used the older evidence checkout's crawler and the current source checkout's crawler as if they were the same revision. Reading the exact source PR code and rerunning both site directories resolved the mismatch. No source-code change was needed.

## 10. Known Residuals

The local pages have not received visual review. PR #231 has not been merged or deployed, so live routes are not verified against this candidate. The final source tarball and platform checks are not bound to one checksum. BirdTree redistribution rights remain undocumented in the reviewed public terms.

## 11. Team Learning

Keep crawler results paired with both the exact crawler revision and the exact generated-site directory. A newer crawler may correctly reject an older site snapshot; that failure must be preserved as a failed historical result, not confused with the current clean build.

Memory receipt: Loaded the pigauto LOAD-FIRST manifest with `route.py pigauto`, ran the lane preflight, and checked the brain status record. The recovery-to-truth and lane-ownership guidance shaped this slice. Golden Set: not run; this was an evidence-only update and changed no package behavior.

## 12. Cross-Product Coverage

Covers the local pkgdown retirement checker, its negative control, and the retained evidence PR.

Does NOT cover local visual review, deployed-site behavior, redistribution clearance, final tarball checks, or CRAN submission.

## Editorial self-check

Style: 2/10, medium confidence; the report uses concrete paths, counts, and command outcomes, with one explicit correction rather than generic readiness language. Genre/coverage: repository after-task report, whole draft. Evidence/repair: the stale snapshot and fresh build are named separately; no further prose change needed. Gates: science = not applicable, facts = pass for the listed live checks and local command output, references = pass for the linked PR pages and recorded source commit. Provenance: self-review, draft `after-task-2026-10-08-site-evidence-reconciliation.md`, 2026-10-08; the route control was already observed and was not treated as a held-out test.
