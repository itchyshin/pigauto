# After-task report: local site gates

## 1. Goal

Bind the local site build and route checks to the current PR #231 source candidate in the CRAN 0.11 evidence ledger.

## 2. Implemented

Updated G4 evidence and G5/G5b commands for source commit `2f7305c7ce4d13e03f676f4e48818f3a259d5473`, then approved and reran both site gates through unlazy. Recorded the fresh build and crawl results, removed the stale lane claim from the gate evidence, and refreshed PR #228's status description. Tried the local rendered page in Chrome for G6; browser policy rejected its `file:` URL and prohibited workarounds.

## 3a. Decisions and Rejected Alternatives

Kept the site checks bound to the exact source commit. Did not start a local server or use another browser surface after Chrome's policy rejection. Did not treat structural HTML and link checks as a visual review.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-local-site-gates.md`

## 5. Checks Run

- Unlazy G5 and G5b passed with `--reverify --approve --timeout 600` against PR #231 commit `2f7305c7ce4d13e03f676f4e48818f3a259d5473`.
- The first attempt with the default 120-second limit timed out while rendering vignettes. The approved 600-second rerun completed successfully.
- The crawler reported 62 HTML pages, 3,443 local references, 34 retired routes, zero errors, and `SITE_CRAWL_OK`.
- Rendered-content inspection found the updated BirdTree download guidance, tree provenance, and legacy-only inverse-Wishart wording.
- Chrome refused the local `file:` URL under browser policy. G6 remains unmet.
- PR #228's description was updated and rechecked in Chrome; it remains Draft and unmerged.
- `gate-check --status` reports 6 met and 5 unmet gates: G0, G6, G7, G8, and G9.
- No package tests ran; this was an evidence and site-verification update.

## 6. Tests of the Tests

The crawler returned `SITE_CRAWL_OK` on the clean generated site. The separate route-retirement negative control remains recorded in `provenance/site-crawler-simulation-retirement-control.txt`; it was not rerun in this slice.

## 7a. Issue Ledger

G4, G5, and G5b are verified for the pinned PR #231 source candidate. G6 is still open because visual inspection could not proceed through Chrome's allowed URL protocols. G7 remains open until the candidate is merged and its deployment is checked. G0, G8, and G9 remain open in the release ledger.

## 8. Consistency Audit

The generated Getting Started and tree-uncertainty pages describe `tree300` as a posterior sample member and direct users to BirdTree retrieval guidance. The multiple-imputation article identifies the inverse-Wishart option as legacy behavior. The old lease note in G4 was stale; this task is the sole pigauto lane. Local route retirement does not establish that the deployed site has changed.

## 9. What Did Not Go Smoothly

The default 120-second unlazy limit was too short for the complete pkgdown build. The build had progressed into vignette rendering when it was stopped. Increasing the approved per-check limit to 600 seconds allowed both site gates to pass.

## 10. Known Residuals

The local site has not received visual review. The deployed pages and their search and sitemap routes have not been checked after this candidate's changes. The release tarball and its platform results are not bound to this source candidate.

## 11. Team Learning

The candidate site's current full build needs more than the runner's default 120 seconds; use a measured 600-second ceiling for this exact offline build. The repo has one active pigauto lane, this task.

Golden Set: not run because this update changes no package behavior.

## 12. Cross-Product Coverage

Covers the local generated website and route cleanup for PR #231 commit `2f7305c7ce4d13e03f676f4e48818f3a259d5473`.

Does NOT cover local visual inspection, post-merge deployment, deployed routes, the final tarball, CRAN readiness, or submission.
