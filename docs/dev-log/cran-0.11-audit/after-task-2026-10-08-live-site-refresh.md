## 1. Goal

Refresh direct Chrome evidence for the live pigauto homepage, article index, and retired simulation-study route. Keep deployment verification unmet until the source changes are merged and deployed.

## 2. Implemented

Added a dated, cache-busted live-route observation to `site-review-2026-10-08.md`. No package code, generated help, or website source was changed in this slice.

## 3a. Decisions and Rejected Alternatives

Used Chrome with unique query parameters to check live responses. Did not use a local server to preview candidate files because the browser policy recorded in the existing site review prohibits that route. Did not infer sitemap state after Chrome blocked its URL.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/site-review-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-live-site-refresh.md`

## 5. Checks Run

- `lane_preflight.sh` reported one pigauto lane and confirmed no other ref had unmerged work on the site-review file.
- Chrome opened the live homepage, article index, and retired route with `?audit=20261008codexlead`.
- Chrome showed the historical article linked from the live article index and served at its direct route; its body says it is excluded from the 0.11.0 site.
- Chrome showed the homepage's 0.11.0 badge and current default and MI guidance.
- Chrome blocked the sitemap request with `net::ERR_BLOCKED_BY_CLIENT`; sitemap status remains unknown.
- `check-after-task.R` passed its structural check; `slop_check.py` found 0 hits in both changed prose files; `git diff --check` passed.
- `closeout.py check` failed after structural validation because the brain-wide evidence gate reports five unrelated open imputation-simulation gates.

## 6. Tests of the Tests

No code test was added or run. The index observation establishes that the route is discoverable, and the separate direct request establishes that it is served. The unique query parameter reduces reliance on a previously cached page. This does not test post-deployment retirement.

## 7a. Issue Ledger

- G7 remains open: the retired simulation-study route is still linked in the deployed article index and served directly.
- The live sitemap check is unresolved because Chrome blocked the URL.
- The source/docs PR remains Draft and undeployed; this slice does not authorize merge or deployment.

## 8. Consistency Audit

The homepage's visible 0.11.0 guidance agrees with the candidate defaults and MI documentation described in the existing release evidence. The article index contradicts the candidate retirement result because this deployment predates those source changes. No other retained pages were assessed in this slice.

## 9. What Did Not Go Smoothly

The sitemap URL was blocked by the browser. The relative-path invocation of `closeout.py new` resolved under the brain repository and failed before creating a report. Its absolute-path retry was denied by the filesystem sandbox, so the report was written through an approved scoped escalation. `closeout.py check` then failed on five unrelated open imputation-simulation gates after the structural check passed. No out-of-repository file was created.

## 10. Known Residuals

The live article remains public. The Pages deployment commit, live sitemap, live search index, local visual rendering, post-deployment route behavior, final source tarball, and platform checks are not established here. The BirdTree redistribution basis remains unresolved.

## 11. Team Learning

For live-site retirement checks, inspect both the discovery page and the direct retired URL with cache-busting parameters. A source crawler cannot establish that GitHub Pages has deployed the cleaned output.

Memory receipt: loaded the pigauto LOAD-FIRST manifest with `route.py pigauto` and the repo AGENTS.md; these directed the live-route audit and separation of local from deployed evidence. Golden Set: not run; this slice did not modify package behavior.

## 12. Cross-Product Coverage

Covers the live homepage text, article index entry, and direct simulation-study route.

Does NOT cover the sitemap or search index, other retained articles or reference pages, local visual inspection, post-deployment behavior, package tests, or the exact CRAN tarball.

Style: 2/10, medium confidence; concise technical closeout. Genre/coverage: this report and the added site-review section. Evidence/repair: “still lists” and “returned the full historical article” state the two observed outcomes directly; no change needed. Gates: science = not applicable, no scientific result assessed; facts = pass for the listed browser observations; references = not applicable, no external citations added. Provenance: Codex self-review, live-site-refresh report, 2026-10-08; prior route-control evidence was already known, and this cache-busted recheck is new.
