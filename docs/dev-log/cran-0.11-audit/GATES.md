# Gates: pigauto CRAN0.11 audit

OWNS: docs/dev-log/cran-0.11-audit/, cran-comments.md

This is the current nine-gate ledger from the approved plan. Its preceding version is preserved verbatim in historical-gates-through-direct-platform-2026-10-17.md. Source-only PR #237 needs merge/deployment approval; PR #236 stays unmerged as the evidence PR.

The c2 archive remains platform-clean for source `9a457eb`. The new installation guidance changes shipped README/vignette bytes. Updated source is conservatively source-clean, next tarball-clean, pending clean post-merge freezing and new checks. Overall NOT_READY.

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
- [ ] G6: Actual matching-main deployment and live retained/redirected/retired URLs pass.
  EVIDENCE: Pending explicit PR #237 source merge and automatic Pages deployment approval. Latest verified deployment remains 780c084. Require matching successful deployment, 44 retired/62 retained live routes, discovery checks and changed-page delta browser review, including the Getting Started guide.
- [ ] G7: A clean merged-source final artifact passes every affected exact-artifact gate.
  EVIDENCE: Pending source merge. Rebuild without silently normalizing Packaged timestamps; freeze hash/size/inventory, then local/platform/installed-MI/rights/dependency/URL/spelling/timing checks. c2 is predecessor evidence for the changed installed bytes.
- [ ] G8: Final independent closure, timing decision and submission-ready verification pass.
  EVIDENCE: Pending final artifact/site, fresh Grace/Rose/Pat READY votes, four conditional evidence rows, final cran-comments, maintainer timing decision, validator controls and unlazy reverification. Prior bounded PASS prose is not a final READY vote. CRAN upload/confirmation need separate explicit approval afterward.

Current plan: 6 met, 3 unmet, 0 abandoned. Local c2 warning classification covers all 161 warnings; Windows warning identities remain unclassified. No merge, deployment or CRAN submission occurred in this continuation. See the preparation receipt for evidence and limits.
