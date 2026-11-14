# Rose independent final release review

Verdict: **NOT_READY pending final record consolidation and timing disposition.** This is an independent review of the evidence for the exact archive below. It does not authorize upload.

## Identity and evidence reviewed

- Archive SHA-256: `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`; 5,137,648 bytes; 250 entries.
- Clean source commit: `8e275ca738129c8d608f06ffc407b8141c1a2636`.
- Direct platform run: `39082452558`; native checks identify the same archive on Windows release, Windows devel, Linux release, and local macOS.
- Fresh four-cell warning audit: `/private/tmp/pigauto-final-closure-20261017/final-warning-audit/actual-results-null-error-fix/audit.json` and `report.md`; all four warning-object audits pass and raw native counts agree. The root readback receipt rehashes the audit inputs. The first failed audit and the JSON-NULL reader correction remain preserved. The warning objects are separate replay captures, not the original native reporter objects.
- Reviewed the exact-artifact platform bundles, installed-MI manifest and verification, final provenance and rights carry-forward record, bounded recovery results, CRAN comments draft, closeout draft, current GATES and release ledger.

## Findings

The source/archive identity, source-to-archive binding, inventory, unchanged package R source, data and NOTICE bytes, and recorded consent/tree-use provenance agree. The rights record is a documented carry-forward for unchanged material, not an independent legal opinion or a new redistribution licence. The source PR is merged and the matching Pages deployment passed; PR236 remains open and unmerged.

The direct Windows release, Windows devel, and Linux release checks report Status OK and exit 0. The macOS check completed with exit 0, zero errors and zero warnings, but Status NOTE for `compare_methods` (3.524 CPU seconds; 18.406 elapsed seconds). The note is explicitly still **UNADJUDICATED**. The warning audit correctly leaves that native verdict untouched. This needs the maintainer's timing disposition before a submission decision; it is not evidence of a failed check or a warning-audit defect.

The fresh warning audit resolves the platform warning-identity question: 3,195 pass, zero fail, 161 warnings, and 86 skips on each platform, with message/object identities reconciled. Its consequences explicitly limit small calibration samples and conformal sample size, tree-draw stability, selected covariate routes, ordinal imbalance, documented fallback, unmatched species, constant traits, and undefined correlations. It does not certify general coverage or inference.

The installed optional-backend matrix establishes installation and bounded interoperability across neither backend, drmTMB only, gllvmTMB only, and both. Its real-object checks use three fits per backend. The separate ten-seed, M=20 recovery work supports only the specified n=300 Gaussian fixed-effect regimes; gllvmTMB's one-seed mechanism check is not a campaign result. No evidence here supports general inference, coverage, joint covariance, loadings, variance-component pooling, live external-service availability, or GPU hardware execution. Optional model adapters remain available while neither backend is in pigauto's dependency fields; fresh libraries use recorded dependency symlinks and do not prove clean-room recompilation.

The DOI resolver returned 403, while Crossref returned 200 and matched the cited Hackett et al. 2008 record. The result supports citation identity and records access denial; it does not mean the URL was fixed. The spelling review found no ordinary typo among the 404 flagged tokens in the reviewed contexts. The size finding is INFO only. The inherited `_problems` fixture is inventoried and byte-identical to its predecessor, not treated as an executed test or a newly added fixture.

## Concrete blocker before upload

The current `task/final-closeout-draft.md` still says the four-platform warning audit is pending in Sections 5 and 10, although the fresh audit now passes. The current canonical `GATES.md` and `release-ledger.json` also retain stale pending prose for direct-platform, warning, URL/spelling/timing and final-artifact work. The underlying proofs are present; this is a records contradiction, not a missing native run. Consolidate those records against the fresh actual audit and the remaining evidence, preserve the macOS NOTE as unadjudicated until a maintainer decision, and rerun the final closure validator/readback before upload could be considered. The current drafts themselves are not upload-ready records.

## Highest proven rung and next unproven rung

Highest proven rung: exact-artifact platform checks completed on the recorded archive, with three Status OK results and one macOS Status NOTE, plus the passed four-cell warning-identity audit and bounded installed-MI checks. Source merge and deployed-site verification are also proven at `8e275ca`.

Next unproven rung: a reconciled, validator-passing final release ledger and closeout with the macOS timing NOTE disposition recorded. Beyond that, CRAN upload and maintainer confirmation remain separate and unauthorized steps; no CRAN acceptance or live publication claim follows.
