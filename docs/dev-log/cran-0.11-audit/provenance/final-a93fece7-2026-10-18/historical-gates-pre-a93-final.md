# Gates: pigauto CRAN 0.11 audit

OWNS: docs/dev-log/cran-0.11-audit/, cran-comments.md

The approved nine-gate plan remains in force. PR #237 merged and deployed with explicit approval. Evidence PR #236 remains draft and unmerged. Earlier gate records are preserved verbatim under provenance, including the record immediately before the new help finding.

Archive `14ee7e5d` is withdrawn as a submission candidate. Pat found installed help promising GNN training or comparison despite the default `gnn = FALSE`. Its native checks, installed matrix and warning audit remain valid predecessor evidence. They cannot validate the corrected installed bytes proposed in PR #238. There is no frozen successor yet.

- [x] G0: Scope, ownership, publication history and provenance are verified.
  EVIDENCE: Canonical dirty checkout preserved. Primary CRAN index refreshed 17 October 2026 lists0.10.0;0.11.0 was unpublished. PR237 approved merge8e275ca and deployed. Scoped help-fix branch816987fd is based on8e275ca. See provenance/help-fix-816987fd-2026-10-18/withdrawal-and-successor.json.
- [x] G1: Effective-default inventory and forwarding/configuration checks pass within their recorded scope.
  EVIDENCE: defaults-inventory.md and202-expectation source audit; all49 executable R syntax trees of proposed source816987fd match8e275ca. The help change adds no runtime/default decision. Candidate-specific installed verification still belongs toG7.
- [x] G2: Optional-backend dependency independence is established within its recorded installation scope.
  EVIDENCE: DESCRIPTION remains unchanged; drmTMB/gllvmTMB are absent from dependency fields. Predecessor14ee built scans and four availability cells pass. A successor scan and installed matrix remain required atG7.
- [x] G3: Bounded analysis workflows, independent arithmetic and recovery outcomes are recorded.
  EVIDENCE: recovery/RESULTS.md retains n=300, M=20, ten-seed Gaussian zero-phylogenetic-covariance checks and bounded bias results. They do not certify nominal coverage. Predecessor14ee installed matrix passes12steps with three-fit fixtures; no new campaign is justified by the comment-only correction.
- [x] G4: Audited documentation mismatches have a tested proposed-source resolution or an explicit bound.
  EVIDENCE: PR238 changes two roxygen descriptions and generated help only. Provisional install and installed novice workflow pass; default/configured GNN false,2missing cells filled, observed values/classes preserved. Existing report expectations17pass/0fail/0skip with2expected tiny-validation warnings. Pat's separate source-fix review passes. This does not fix the withdrawn archive or deployed pages; those gates remain open.
- [x] G5: Fresh local website build and bounded rendered-page review pass for source238.
  EVIDENCE: Fresh committed816987fd site build exits0;62HTML,3588local references,34retired pages,0crawler errors; search614/sitemap62. Independent Chrome review of both changed reference pages and adjacent Getting Started agrees with baselineFALSE/opt-inTRUE; readable at1200px with wrapped titles. Temporary server stopped. First missing-Rscript PATH failure and corrected successful run both retained. This is local-source evidence only.
- [ ] G6: The corrected source is approved, merged, deployed and its live pages and route controls pass.
  EVIDENCE: Source238 merge and automatic deployment need explicit approval. Earlier8e deployment passed, but its two incorrect reference descriptions remain until the source correction deploys. Verify the actual new deployment commit and both changed pages, plus retirement/discovery controls.
- [ ] G7: A new clean merged-source frozen artifact passes all corresponding exact-artifact gates.
  EVIDENCE: Successor not frozen. Predecessor14ee has three direct remote StatusOK results, Mac exit0/one timingNOTE and strict-driver exit1, installedMI12stepsPASS and warning audit4/4PASS. Grace technically accepts its disclosed NOTE for preparation; human timing disposition remains separate. Preserve all these results without counting them for changed installed bytes.
- [ ] G8: Final independent closure, maintainer timing decision and submission-ready verification pass.
  EVIDENCE: Original14ee panel retained: Grace technical-preparation READY with conditions; Rose NOT_READY record/timing; Pat NOT_READY help. New candidate needs fresh appropriate final review, reconciled comments/ledger, validator controls and unlazy readback. CRAN upload and maintainer confirmation require separate explicit approval.

Current approved-plan count:6met,3unmet,0abandoned. Highest historical technical evidence is predecessor14ee platform checking with its disclosed NOTE, not a current submission candidate. Corrected source238 has passed all three source CI jobs and local-site review; merge/deployment approval is pending. Nothing has been submitted to CRAN.
