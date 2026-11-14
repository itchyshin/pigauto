# Remaining-gate preparation checkpoint

## 1. Goal

Continue PR #236's remaining release gates and make the approval sequence concrete. This checkpoint preserves the approved unmerged-evidence-PR deliverable.

## 2. Implemented

Created source-only draft PR #237 from main `780c084` with 11 independently reviewed files. It carries the already-tested portability examples and new installation-version guidance. CRAN metadata still lists 0.10.0; the guide now explains how to obtain the GitHub 0.11.0 workflow. Updated generated R and paired help stay consistent.

Captured all 161 warnings in a corrected local CRAN-layout replay of the unchanged c2 archive. Independent review matched every warning object and classified its consequence. Rebuilt the guide pages and verified their rendered text and layout in Codex's own browser.

## 3a. Decisions and Rejected Alternatives

Keep PR #236 unmerged. A separate source PR #237 allows a main deployment and clean post-merge freeze without merging governance evidence. Reject silently reusing c2 platform results for the changed README/vignette bytes. Reject suppressing warnings or weakening tests to obtain a cleaner report.

Assumptions made without asking: existing implementation/check authority covers the focused reader repair, source-only draft PR and bounded diagnostic replay. It does not cover source merging, deployment authorization or CRAN submission. No new compute campaign was inferred.

## 4. Files Touched

Source PR #237: README.md; NEWS.md; R/evaluate.R, R/plot.R and R/simulate_benchmark.R; their three paired help files; tests/testthat/test-new-features.R; vignettes/getting-started.Rmd and its generated R script. The extraction manifest lists all 11 paths and hashes.

Evidence PR #236: cran-comments.md; release-ledger.json; current GATES.md and its verbatim historical predecessor; this report; and provenance/submission-preparation-2026-10-17/, whose manifest lists the warning objects, rows, classifier, sessions, failed diagnostic, primary-record refresh, source extraction, local site logs and reviews. Task-only execution scratch remains under /private/tmp/pigauto-final-closure-20261017.

## 5. Checks Run

Exact c2 replay with NOT_CRAN unset, force-Suggests, a working torch tensor and one-thread caps: 3,195 pass, 0 fail, 161 warnings, 86 skips; 183.448 seconds for tests, 205.568 seconds including install. All 161 rows classified; unknown-warning control rejected. Source-surface tests, R/Rd parsing, git diff checks and pkgdown::check_pkgdown passed.

Fresh local site: 62 HTML, 3,587 references, 34 archived retired pages, zero errors; 44 retired-route checks passed; sitemap 62/search 614/path values 568 measured; six crawler controls passed. Source PR #237's fresh CI and final evidence-head CI are checked after publication and recorded in PR descriptions. Earlier exact c2 platform checks remain bounded predecessor evidence.

## 6. Tests of the Tests

An injected unknown warning is refused by the classifier. Independent review matched all 161 RDS messages, unique IDs and callers to CSV rows. The bad installed-directory replay produced four source-only failures and was retained separately; the corrected native layout restored the exact CRAN-mode totals. The first bad site copy failed the rendered-note assertion before a corrected build passed.

## 7a. Issue Ledger

Fixed: install-version ambiguity and stale CRAN-comments prose. Closed locally: 161 unclassified warning messages, with their real consequences retained. Preserved: unclassified Windows warning identities, measured timing, 86 skips, source-merge/deployment approval, final artifact freezing and final reviewer votes.

## 8. Consistency Audit

The source reviewer verified exactly 11 changed source files and 216 matching shared package files. NEWS is shipped; it is not build-excluded. The generated Getting Started R script agrees with its Rmd. Historical gates were archived verbatim so the current ledger has the approved nine definitions. Plan-versus-actual reconciliation is in the preparation receipt.

## 9. What Did Not Go Smoothly

The first warning diagnostic used installed-package layout and falsely activated source-only checks. The corrected replay changed the harness layout, not package tests. The first site scratch copy misplaced the Rmd and was caught by its rendered-content assertion. Relative slop-check paths resolved against the brain rather than the worktree; absolute paths were then used. No failed attempt was relabelled as successful.

## 10. Known Residuals

Overall NOT_READY: six approved-plan gates met, three open, none abandoned. Updated source awaits post-merge freezing and tarball/platform evidence; c2 remains platform-clean only for its own identity. Local warning identities do not prove Windows identities. Final timing approval, four conditional ledger rows, fresh Grace/Rose/Pat verdicts and exact submission text remain pending.

## 11. Team Learning

Audit installation instructions against the actually published version before deploying new-default documentation. Preserve original archive identity when doing diagnostics. CheckReporter captures warning objects even when CRAN-mode printing hides them. Source tests need the native check-directory layout for a faithful replay.

Memory receipt: used the repository LOAD-FIRST route, the brain-first search with repo fallback, CRAN release protocol, ultra-plan routing, unlazy and after-task disciplines. The known partial-green release mistake remains covered by the previously rerun validator controls; no new validator algorithm change or fresh campaign was made. Golden-Set regression was not rerun in this documentation-only continuation. No brain or Codex memory files were written.

## 12. Cross-Product Coverage

Covers: focused source examples, CRAN/GitHub installation guidance, local c2 warning classification, bounded local site rendering, primary version/reverse-dependency refresh and the concrete merge/submission approval route.

This checkpoint does NOT cover actual source merge or deployment, the new post-merge archive, final Windows warning messages, new scientific default decisions, broad inferential validity, final release-panel READY votes, CRAN upload, confirmation or acceptance. Those remain explicit gates.
