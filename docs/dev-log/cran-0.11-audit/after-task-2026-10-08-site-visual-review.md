## 1. Goal

Close local site review for the pigauto 0.11 candidate and refresh the release ledger against the current source commit.

## 2. Implemented

Updated `GATES.md` with the fresh site build and crawl evidence for commit `4a724399ada88eff6a4b347b94e9047c2e0163d4`, Chrome visual review, and current PR #231 CI state. The exact-head output has 62 HTML pages, 613 search entries, 3,443 checked local references, and 34 retired routes, with no crawl errors. Chrome reviewed Getting Started and `trees300` directly; the other three pages checked earlier are byte-identical to this build. Corrected old “current” snapshots in G0 and G7 to label them historical. The lane census reports one active pigauto lane. A later exact-head matrix on `dcdf0366cc6b99729cb7dc85045d224edab7d4b8` passed Ubuntu R release, Ubuntu R-devel, and macOS R release in 19m42s total. It is candidate-source CI only. G6 is met; G0, G7, G8, and G9 remain unmet.

## 3a. Decisions and Rejected Alternatives

Used Chrome with the existing loopback-only site server after direct `file://` navigation failed. Kept local visual review distinct from live deployment verification. Did not treat candidate-source checks as frozen-artifact evidence. Did not merge, deploy, freeze a post-merge tarball, or submit to CRAN.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-site-visual-review.md`

## 5. Checks Run

- `bash ~/shinichi-brain/tools/lane_preflight.sh /Users/z3437171/Dropbox/Github\ Local/pigauto` reported one active lane and 58 total worktrees.
- Checked the retained exact-head pkgdown build and crawl logs. Build output records 62 HTML pages, 613 search entries, 3,443 local references, 34 retired routes, zero crawl errors, `SITE_CRAWL_OK`, and a clean `pkgdown::check_pkgdown()` result. Both log checksums are recorded in `GATES.md`.
- Compared rendered HTML hashes against the preceding visually inspected output: `index.html` `63f96fd4e932a45b8789c6844a00b792c94e868254108c68388cc9f7df131b99`; `articles/getting-started.html` `33c9ea9cbbe4fca0c230ca082f94e89cab3c612dac13c58bbaff58639d0e24e6`; `articles/multiple-imputation.html` `4124a2e181083812e00fc4f9790f30f4024af111fcccc68fcb3b2cf0a602e4b5`; `articles/tree-uncertainty.html` `e92e5c69fbb5b34bd8de39e03b487150a7f62074b82d28b0a4e6ac57e403bffc`; `reference/trees300.html` `82aa732d24ea848a4bce76da8716f16271f39113a382d03d0ca7e6bcccae54ea`.
- Reviewed the local rendered home page, Getting Started, multiple-imputation, tree-uncertainty, and `trees300` reference pages in Chrome.
- Chrome showed PR #231 remains Draft at head `dcdf036`. Run #37870408831 passed all three platform jobs in 19m42s: Ubuntu R release 10m24s, Ubuntu R-devel 12m33s, and macOS R release 19m10s (MPS-focused tests 4m52s; `R CMD check` 10m43s). The macOS job displayed one runner-capacity notice; Ubuntu jobs displayed runner-image migration notices. This does not validate a frozen tarball or force-Suggests check.
- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md` reports 7 met and 4 unmet gates. It parses the ledger; it does not rerun the gate commands.
- `git diff --check` passed.

## 6. Tests of the Tests

The site crawler has a retained planted retirement negative control at `provenance/site-crawler-simulation-retirement-control.txt`; it was not rerun in this slice. The visual review is manual and has no automated layout oracle.

## 7a. Issue Ledger

- G0 remains open: rights for redistributing bundled BirdTree-derived data are not established by the external pages reviewed, and no final post-merge source revision or artifact exists.
- G7 remains open: the corrected pages have not been deployed and checked against the full route, sitemap, search, and retirement inventory.
- G8 remains open: no exact post-merge tarball is frozen and checked across required platforms.
- G9 remains open: independent reviewers have not voted on that final artifact and deployment evidence.

## 8. Consistency Audit

Compared G6 and G5/G5b claims with the 4a72439 build logs, crawler output, rendered Chrome pages, and page hashes. The exact-head local output contains the corrected tree provenance and retrieval guidance, fixed-effect MI boundaries, and current reader-facing navigation. G0 rights and G7 deployment remain separate gates. The old deployed-site record and old PR snapshot are now explicitly marked historical in GATES.md.

## 9. What Did Not Go Smoothly

Chrome refused local `file://` navigation. The temporary loopback-only server allowed the requested Chrome review and resolved that limitation. The closeout wrapper resolved its root to the brain repository rather than this pigauto worktree, so the direct after-task structure validator was used with its unrelated ledger recheck disabled.

## 10. Known Residuals

This slice does NOT establish BirdTree redistribution rights, deployed-page correctness, frozen-tarball validity, checksum-bound platform results, or CRAN submission readiness. The site crawler does not check CSS URLs, `srcset`, or JavaScript-generated URLs.

## 11. Team Learning

A browser's refusal of `file://` does not by itself prevent local visual review; a loopback-only HTTP server can provide a normal rendered-page review when the browser permits it. The routed `pigauto` LOAD-FIRST manifest and repository instructions were consulted. Their prediction-path and uncertainty guards did not affect this documentation-only slice. Golden Set regression was not in scope.

## 12. Cross-Product Coverage

Covers: exact-commit local site build evidence, route and search retirement checks, and Chrome visual review of five rendered pages.

Does NOT cover: deployed website state, package runtime behavior, a post-merge tarball, checksum-bound external artifact checks, full CSS/srcset/JavaScript link behavior, redistribution rights, or CRAN submission.
