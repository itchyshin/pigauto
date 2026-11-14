# Gates: pigauto CRAN0.11 audit

OWNS: docs/dev-log/cran-0.11-audit/, cran-comments.md

This is the current nine-gate ledger from the approved plan. Its preceding version is preserved verbatim in historical-gates-through-direct-platform-2026-10-17.md. Source-only PR #237 merged and deployed with explicit approval; PR #236 stays unmerged as the evidence PR.

The c2 archive remains predecessor platform evidence for source `9a457eb`. Final archive `14ee7e5d` is frozen from clean merged source `8e275ca`. All native and installed checks are complete on14ee, with one retained Mac timingNOTE and strict-driver exit1. The actual four-cell warning audit passes. The current source-clean claim awaits final independent timing disposition before promotion; overall submission verdict remains NOT_READY.

- [x] G0: Scope, ownership, publication history and provenance are verified.
  EVIDENCE: Dirty canonical checkout preserved; primary CRAN index refreshed 2026-10-17 lists 0.10.0; remote v0.11 tags absent. PR #237 is 259d7e98 from main `780c084`. See provenance/submission-preparation-2026-10-17/.
- [x] G1: Effective-default inventory and forwarding/configuration checks pass within their recorded scope.
  EVIDENCE: defaults-inventory.md and the 202-expectation current-source audit in the c2 evidence bundle. Subsequent source changes concern examples and installation guidance, without scientific default or runtime changes. The 49 R files match the tested c2 source in the site snapshot.
- [x] G2: Optional-backend dependency independence is established within its recorded installation scope.
  EVIDENCE: DESCRIPTION unchanged; drmTMB/gllvmTMB absent from dependency fields. Built c2 scans and four-cell installations are in the bundle. G7 repeats final-artifact dependency scanning.
- [x] G3: Bounded analysis workflows, independent arithmetic and recovery outcomes are recorded.
  EVIDENCE: recovery/RESULTS.md retains n = 300, M = 20, ten-seed drmTMB/gllvmTMB research-source checks and bounded bias results; these do not certify the final tarball or nominal coverage. The c2 exact-artifact availability matrix passed 12 steps. G7 repeats the final installed matrix.
- [x] G4: Audited documentation mismatches are resolved or explicitly bounded.
  EVIDENCE: PR #237 contains portable examples and CRAN/GitHub version guidance, with generated Rd/R counterparts. Source-surface tests, R/Rd parsing and independent extraction review pass. Fresh source-CI readback is recorded in the PR descriptions. Existing defaults, MI and rights boundaries remain.
- [x] G5: Local website build, structural checks and bounded rendered-page review pass.
  EVIDENCE: Recorded snapshot: 62 HTML, 3,587 local references, 34 archived retired pages, zero errors; 44 retired-route checks pass; sitemap 62/search 614/path values 568 measured. Cleaner/six crawler controls pass. Codex browser reviewed both changed guide pages; prior five-page sample remains bounded. Font/news overrides are recorded and do not prove production deployment.
- [x] G6: Actual matching-main deployment and live retained/redirected/retired URLs pass.
  EVIDENCE: Explicit user approval received; PR237 merged at8e275ca738129c8d608f06ffc407b8141c1a2636. Pages39080327369 succeeded on that commit. Fresh live verifier passes44 retired routes,62 retained routes, sitemap62/search614, no retired targets. Own-browser text/pixel review covers homepage, Getting Started, NEWS and the three changed reference pages. Narrow-viewport code scrolling and long NEWS sidebar truncation are retained limits. See provenance/post-merge-8e275ca-2026-10-17/.
- [ ] G7: A clean merged-source final artifact passes every affected exact-artifact gate.
  EVIDENCE: Final archive14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb is frozen from clean merged8e275ca:5,137,648 bytes,250 entries. Independent source/rights/inventory binding PASS; exact public-object readback matches. Direct Windows release/devel/Linux checks report StatusOK; Mac native exit0/Status1NOTE and strict-driver exit1 are retained. Final installed-MI matrix passed12steps. Four-cell separate-replay warning audit passed with input hashes reverified. URL/spelling/vignette/size/service proofs are complete within stated limits. Only final timing disposition remains open for this gate. No c2 artifact check is counted for these new bytes.
- [ ] G8: Final independent closure, timing decision and submission-ready verification pass.
  EVIDENCE: Pending final exact-artifact checks, fresh Grace/Rose/Pat READY votes, four conditional evidence rows, final cran-comments, maintainer timing decision, validator controls and unlazy reverification. Prior bounded PASS prose is not a final READY vote. CRAN upload/confirmation need separate explicit approval afterward.

Current plan: 7 met, 2 unmet, 0 abandoned. Final14ee actual separate-replay warnings and skips are audited on all4platforms; native counts match. The earlier c2 classification is a comparison contract, not final-artifact evidence. Source PR237 merged and deployed with explicit approval. PR236 remains unmerged. No CRAN submission occurred. See the preparation receipt for evidence and limits.
