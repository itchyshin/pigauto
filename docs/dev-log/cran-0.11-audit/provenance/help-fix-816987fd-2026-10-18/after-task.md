# Pigauto 0.11 source-help approval checkpoint

## 1. Goal

Finish the authorized checks for evidence PR #236 and identify the next approval needed before CRAN submission. This report closes preparation of the help-correction slice. The release goal remains open.

## 2. Implemented

PR #237 merged with explicit approval as source 8e275ca; matching source CI and website deployment passed. Its frozen archive 14ee7e5d completed three direct remote native checks, an independent Mac check, the optional-model installation matrix and a four-platform warning audit. The final reader review then found two installed-help descriptions that contradict the baseline-only default. That archive is withdrawn as a submission candidate. Its bytes, checks and original panel verdicts remain preserved.

PR #238 at 816987fd changes the descriptions for impute() and pigauto_report() and their generated Rd pages. Both now distinguish default baseline-only behavior from explicit GNN training/comparison. All 49 executable R syntax trees match merged 8e; no runtime, scientific default, dependency or example change is included. Fresh local rendered-page review and all three source CI jobs pass. PR #238 is ready for review and cleanly mergeable; explicit merge/deployment approval remains pending.

## 3a. Decisions and Rejected Alternatives

Withdraw the defective candidate rather than declaring its help fixed by an unmerged source patch. Preserve successful and failed evidence, including the original Mac timing NOTE and strict-driver exit1. Keep the three original 14ee votes unchanged: Grace supports technical preparation with conditions; Rose requires record/timing reconciliation; Pat refuses the help contradiction. Pat's separate proposed-source PASS does not replace the archive verdict.

Assumptions made without asking: the existing implementation/check/commit/push/PR authority covers this bounded four-file repair and publication of its evidence to236. Merge authority remains separate under the approved plan. The source is comment-only by executable syntax-tree comparison, so prior scientific campaign evidence can be retained within its regime. That comparison does not waive fresh checks for changed installed archive bytes.

## 4. Files Touched

Source #238 contains only R/impute.R, R/report.R, man/impute.Rd and man/pigauto_report.Rd. Evidence236 updates cran-comments.md with a historical-draft warning, the canonical GATES.md and release-ledger.json, and adds provenance/final-14ee7e5d-2026-10-17/ and provenance/help-fix-816987fd-2026-10-18/. Per-tree manifests identify published evidence files and hashes. A staged Git-blob readback verified all 339 predecessor files against their working bytes, including preserved raw CRLF and whitespace. Those raw/derived findings were retained rather than normalized; authored-wrapper whitespace checks pass. All preceding ledger/wrapper versions are retained. The canonical dirty checkout and foreign worktrees remain untouched.

## 5. Checks Run

Source help was regenerated. Existing report tests passed 17 expectations, 0 failures, 0 skips with 2 expected tiny-validation warnings. The provisional install exited0 in 56.077 seconds; the installed-help/novice check exited0 in 11.712 seconds. It confirms configured/default GNN false, 2 missing cells filled, observed values/classes preserved, and corrected installed help. Browser launch was suppressed only in that headless R process. No browser permissions or persistent settings changed. Independent Pat review passed the proposed source and provisional install. Source CI 39087559368 passed Linux release/devel and macOS release at exact head 816987fd, in 10m53s/12m11s/15m53s respectively. The PR readback shows OPEN, ready and CLEAN, without a merge. Its deployment job is deliberately skipped on pull requests. Fresh site build and crawl pass: 62 HTML pages, 3,588 local references, 34 retired pages, zero errors, search 614 and sitemap 62. Rose reviewed both changed reference pages and adjacent Getting Started in Chrome at 1200 pixels; the default/opt-in text agrees and titles remain legible when wrapped. The temporary local server was stopped. The first build attempt failed because its restricted PATH omitted Rscript; its log and the successful corrected-PATH run remain retained.

Predecessor 14ee direct run 39082452558 reports native StatusOK for Windows release/devel and Linux release. Totals 645.813/622.432/336.330 seconds are retained. Mac native exit 0 reports one compare_methods timing NOTE: user 2.532, system 0.992, elapsed 18.406 seconds; native total 3564.464 seconds. Its strict driver exits 1. No native rerun was used to replace the NOTE or select faster timings. All four native suites report3,195 pass / 0 fail / 161 warnings / 86 skips. Actual warning objects from separate replays match the native totals and pass the four-cell identity/consequence audit. The native reporter retains totals, not those individual objects.

The predecessor optional-model matrix passed 12 steps over neither backend, each separately, and both. These are fresh pigauto installs with recorded dependency symlinks and three-fit real-object fixtures. They do not imply clean-room dependency compilation. Existing n = 300, M = 20, ten-seed Gaussian recovery remains bounded to 20% MCAR, zero true phylogenetic covariance and the recorded posterior settings. Missing-rgbif/terra tests pass 60 expectations with 8 expected skips, without proving live services. Six vignette diagnostics pass ; 404 spelling flags were reviewed in 1,886 contexts across 43 archive-bound files. A retained DOI403 is distinct from successful primary Crossref identity confirmation.

## 6. Tests of the Tests

The warning reader's first actual audit failed on successful R NULL serialized as JSON{}. Minimal serializer and four-receipt reproduction preceded a strict reader correction. Independent review exercised 142 cases including 136 refusals and checked 138 protected files unchanged. Both first failure and successful fresh four-cell audit remain available. Warning/skip identity controls reject unknown messages, changed callers or source locations, changed skip identities and mismatched native/replay counts. Timing parser controls and raw-stage comparisons preserve the distinction between CPU, process wall time and enclosing-loop time.

Curated-copy review verified 202 indexed files, 7 bundles, 228 original members and 108 views, with planted hash/size/path refusals. The final release-checkpoint verifier now explicitly refuses a withdrawn current_candidate=false. The generic validator's source-clean result is a source-only claim, not release readiness. The source-approval checkpoint verifier exercises one positive and seven refusal controls; it separately checks source CI, source/site identities, provisional installed behavior, original panel preservation and published-file hashes. Its full readback is run after consolidation.

## 7a. Issue Ledger

Current plan: 6 met, 3 unmet, 0 abandoned. G0-G5 have scoped evidence; G6 awaits the approved merge/deployment, G7 a successor archive and exact checks, G8 final review/timing/submission-ready verification. The local checkpoint counts preparation separately from human/source and release gates. A pending merge or submission decision is not a passed gate.

## 8. Consistency Audit

The new source descriptions match the observed baseline-only novice report. Neighboring README, Getting Started, mixed-type and MI text were checked; dedicated comparison functions intentionally document GNN-enabled defaults. The report title describes an imputation report. The full 49-file executable comparison supports runtime equivalence while the four changed installed files require a new archive. Canonical records explicitly identify 14ee as predecessor and preserve old point-in-time pending statements as snapshots with their original dates.

## 9. What Did Not Go Smoothly

Final reader review caught a material documentation defect after costly artifact checks had completed. The fix is small, but the corresponding installed bytes must be checked again. The warning audit initially stopped on the NULL/{} envelope mismatch; strict reproduction and review prevented an unsafe broad error bypass. A Linux timing parser initially omitted CPU/elapsed pairs and a vignette summary conflated process and enclosing time; raw-line comparisons corrected the derived summaries. A copied reproduction instruction named an absent helper; its actual helper and original scratch-input requirement were corrected without altering raw evidence.

The shared-host Mac check exceeded its initial estimate; revised bounds and actual times were reported. Other users' processes were preserved. The cause of the shared-host CPU/wall-time gap remains unestablished. One command used the wrong agent-mention-check argument shape and failed on a literal paragraph interpreted as a filename; rerunning with the actual body-file path passed. No PR body shipped before that corrected check. A broad after-task invocation began checking inherited imputation-sim ledgers; its own controller tree was stopped, preserving unrelated processes and receipts. The report structure passes separately, while the whole-goal acceptance result remains open. The first handoff check was inconclusive because its restricted PATH omitted Node; the corrected check returns 1 with inherited ledgers refused and current human gates still open. Full historical branch inventories remain local; the public packet contains a bounded disposition and receipt hashes. Subsequent handoff checking refuses unowned historical ledgers without executing their commands.

## 10. Known Residuals

Source #238 merge and automatic deployment are not approved. No successor archive exists. Its corresponding native/installed/documentation checks and fresh final review remain necessary. Predecessor Windows totals exceed 600 seconds and Mac has one timing NOTE; the final maintainer timing decision must use the successor's actual results. The DOI resolver still returns403 to automated checking. Local/native timing, live services and platform-specific installed size do not imply CRAN acceptance.

PR #236 stays draft and unmerged. No CRAN upload, confirmation, tag or announcement occurred. Only the next source merge/deployment approval should be requested at this checkpoint. Separate final CRAN authorization must identify the new archive and exact submission text after its gates pass.

## 11. Team Learning

A final reader must inspect installed help as well as the website and novice example. Source-only comments can still invalidate an exact installed-artifact claim. Retain source, checker, deployment and gzip identities separately. Keep successful technical checks, reader refusals and human authorization as distinct records. Brain and Codex memory writes were not authorized; durable receipts remain in repository evidence and task scratch.

## 12. Cross-Product Coverage

This slice covers the proposed help correction, its provisional installed novice workflow, and preservation/reconciliation of predecessor evidence. It does NOT cover a corrected frozen archive, new deployment, CRAN submission/acceptance, general interval coverage, joint covariance/loading/variance-component pooling, untested missingness or model families, live service availability, or GPU hardware validity. Source CI does not validate a frozen gzip; R-hub rebuilt-source checks do not prove original-gzip execution.

Golden Set: shipped report tests, novice example, existing source/default checks and retained bounded recovery fixtures; no new campaign. Memory receipt: prior boundaries retrieved from MEMORY.md lines 217-218 and rechecked against current source and retained primary evidence. No memory writes.

Assessment: separate exact-draft record accompanies this report. Highest current source rung: source-clean for proposed 816987fd. Next unproven rung: tarball-clean for a successor after approved merge and deployment. Overall release remains NOT_READY.
