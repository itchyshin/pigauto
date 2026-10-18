# Gates: pigauto CRAN 0.11 audit and authorized submission

OWNS: docs/dev-log/cran-0.11-audit/, cran-comments.md, script/cran-0.11-platform/, .github/workflows/frozen-tarball-check.yaml

The approved nine preparation gates are closed for frozen a93fece7 from clean merged b2e6054. User subsequently authorized CRAN submission once these gates passed; G9 records that separate external action. PR236 stays draft and unmerged. Preserved historical ledgers and predecessor votes are not current claims.

- [x] G0: Scope, ownership, publication history and provenance are verified.
  EVIDENCE: Owned release/cran-0.11-gate-780c084 lane; dirty canonical checkout preserved. Primary CRAN index2026-10-18 reports0.10.0 and zero reverse dependencies. Actual archive and workflow identities are separate in final-a93fece7 provenance.
- [x] G1: Effective-default inventory and forwarding/configuration checks pass within their recorded scope.
  EVIDENCE: Existing defaults-inventory.md and202-expectation audit; all49 executable R syntax trees of source238 match the prior runtime. No scientific default was changed. New installed help/configuration checked on a93fece7.
- [x] G2: Optional-backend dependency independence is established within its recorded installation scope.
  EVIDENCE: Frozen DESCRIPTION excludes drmTMB/gllvmTMB from dependency fields. New exact-archive matrix passes neither,drm-only,gllvm-only,both cells, with installed-unattached and reloaded real-object checks. Dependency symlinks are explicitly recorded.
- [x] G3: Bounded analysis workflows, independent arithmetic and recovery outcomes are recorded.
  EVIDENCE: recovery/RESULTS.md retains n300,M20,10seed Gaussian zero-phylogenetic-covariance results and bounded bias outcomes. New archive real-object matrix passes12steps; no campaign rerun or general coverage claim.
- [x] G4: Audited documentation mismatches have a tested resolution or explicit bound.
  EVIDENCE: Merged source238 corrected impute/report descriptions and generated help. Actual a93fece7 installed help/novice checks pass, defaultGNNfalse,2missingto0, observedvalues/classes preserved. DOI403 remains disclosed; fresh spelling review finds no ordinary typo among404flags.
- [x] G5: Fresh local website build and bounded rendered-page review pass.
  EVIDENCE: Source238 committed local site and independent Chrome review pass, retained in help-fix provenance. No local-only result is substituted for actual deployment.
- [x] G6: Corrected source was approved, merged, deployed and live pages and route controls pass.
  EVIDENCE: User Merge238 exercised; exact mainCI39090594763 and Pages39090594794 at b2e6054 pass. Live check44retired/62retained, sitemap/search and Chrome changed-reference/GettingStarted review pass.
- [x] G7: Clean merged-source frozen artifact passes all corresponding exact-artifact gates.
  CHECK: python3 /private/tmp/pigauto-final-closure-20261017/successor-b2e6054/verify-successor.py preparation
  EXPECT: SUCCESSOR_EXACT_PREPARATION_PROOF_VERIFIED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-pr236-validator-20261017; path=6df32789272f/34 entries; output=SUCCESSOR_EXACT_PREPARATION_PROOF_VERIFIED
- [x] G8: Fresh independent closure, timing disposition and submission-ready verification pass.
  CHECK: python3 /private/tmp/pigauto-final-closure-20261017/successor-b2e6054/verify-successor.py panel
  EXPECT: SUCCESSOR_FINAL_PANEL_AND_EXACT_PROOF_VERIFIED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-pr236-validator-20261017; path=6df32789272f/34 entries; output=SUCCESSOR_FINAL_PANEL_AND_EXACT_PROOF_VERIFIED
- [ ] G9: Authorized exact archive submitted and maintainer confirmation retained.
  EVIDENCE: pending; submission-ready does not prove upload, confirmation, incoming pass, acceptance or live CRAN availability.

Current extended-goal count:9met,1unmet,0abandoned. The original nine preparation gates are met; submission remains pending. Nothing has been submitted to CRAN.
