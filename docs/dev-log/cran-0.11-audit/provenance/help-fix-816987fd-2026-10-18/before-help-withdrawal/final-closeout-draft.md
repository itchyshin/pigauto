# Final pigauto0.11 release-audit closeout draft

## 1. Goal

Finish the remaining PR236 audit gates for one immutable pigauto0.11.0 archive and present the exact CRAN submission packet for separate approval. This is a working draft; final check results and panel votes are pending.

## 2. Implemented

The separately approved source PR237 is merged as main8e275ca738129c8d608f06ffc407b8141c1a2636. Its full source matrix and matching Pages deployment passed. The final source archive was built from clean merged main and frozen at SHA25614ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb,5137648 bytes,250 entries. PR236 retains the independently reviewed checker and partial freeze/site receipts at6cd5f65b06d0c3ee32ce4b2d5a67a22c06028116 and remains unmerged.

## 3a. Decisions and Rejected Alternatives

Keep the selected final archive immutable. Preserve predecessor logs and the entire preceding ledger verbatim. Reuse bounded scientific evidence only after measuring unchanged package R/data bytes. Collect actual warning objects on each platform rather than inferring their identity from counts. Keep the original native R CMD check verdict and timing separate from the diagnostic replay.

Assumptions made without asking: prior explicit public-upload authority covers the same release candidate on the approved package builder; this was verified against session history before upload. No new account token or verification request was sent. The release is an update from published0.10.0, supported by primary CRAN metadata. No scientific default decisions were made in this closeout.

## 4. Files Touched

Working draft. Final scoped path manifest will be attached after the final evidence publication. Source changes were confined to the eleven reviewed files in PR237; the canonical old dirty checkout was preserved.

## 5. Checks Run

Merged-main source matrix39080327341 passed all three legs, including its MPS-related step. That step is not evidence of GPU hardware execution; the separate local torch probe found CPU available and CUDA/MPS unavailable. Matching Pages39080327369 passed. Live verifier checked44 retired/62 retained routes, sitemap62 and search614. Six changed pages were reviewed as rendered text and pixels in Codex's browser. Independent source/archive verification, public object readback and checker publication review passed.

Direct native run39082452558 passed Windows release/devel and Linux release with StatusOK on the exact14ee gzip. Native totals were645.8126/622.4320/336.3305 seconds; Windows tests415/404 seconds and vignette rebuilds22/22 seconds. The local Mac native check completed with exit0 and one compare_methods timingNOTE, elapsed18.406 seconds, CPU3.524 seconds. Native total3564.464 seconds and driver3687.675 seconds are separate measurements. The strict driver returned1, preserving the NOTE for adjudication. All four native suites reported3195pass/0fail/161warnings/86skips. The separate local warning replay completed with exit0 in687.196driver seconds and669.225reporter seconds. The actual four-cell warning audit passed, retaining the nativeNOTE as unadjudicated. All actual input hashes were reverified.

The final installed optional-model matrix passed12steps across4availability cells. Missing-rgbif/terra fixtures passed60expectations with0fail/0warn/8expectedskips. All6vignette diagnostics passed with74.8865seconds render time,308.2825seconds summed R-process elapsed and310.1696seconds enclosing-loop elapsed. All404spelling flags were reviewed in1886contexts across43archive-bound files; no ordinary typo was identified among those flags. The URL checker retains one distinct DOI403, while primary Crossref200 confirms citation identity. These bounds remain explicit.

## 6. Tests of the Tests

The identity helper passed one positive and five refusal controls. The warning parser accepts a complete native summary and refuses a missing summary; the local-NOTE extension rejects14altered cases and leaves the nativeNOTE unadjudicated. Remote timing parsing passed10known-answer controls and15independent stage comparisons. The synopsis passed5positive timing controls and9planted-error refusals. Curated-copy review verified202indexed files,7bundles/228original members and108views, with4planted hash/size/path refusals. Unknown-warning/caller/source/skip refusals passed. Final actual-input readback passed for all4cells. The JSON-envelope guard fix passed142independent cases, including136refusals, with138protected inputs unchanged. No altered archive was submitted to a builder.

## 7a. Issue Ledger

G0 through G6 met. G7 exact-artifact verification and G8 independent closure remain unmet. No approved-plan gate was abandoned. The source/examples/install-guide issues addressed by PR237 are merged. CRAN submission remains unapproved and unperformed.

## 8. Consistency Audit

Verified clean merged-source binding, all49 package R files unchanged from the predecessor, six data files and NOTICE unchanged, and direct/public artifact identities. New installed guidance necessitated a new frozen archive despite unchanged runtime R code. Navigation, retained pages, retired URLs, sitemap and search were checked on the actual deployment. No runtime-only or universal inferential claim is inferred from documentation checks.

## 9. What Did Not Go Smoothly

The shallow evidence clone initially lacked enough common ancestry for a normal main merge. Fetching missing history restored the actual merge base and allowed a normal merge. A PR-description update hit an automatic permission-review timeout; the permitted single retry succeeded. Local R checks exceeded their initial estimates, with measured CPU/wall disparity; revised bounds were reported. A one-second process sample showed active R evaluation/file access. Available memory was53percent of128GiB, so memory pressure is not asserted as the cause.

The first derived remote summary omitted Linux timings because its parser did not handle CPU/elapsed pairs. Independent raw-line comparison caught this; the parser, report and controls were corrected without changing raw evidence. Individual-vignette process totals were also separated from enclosing-loop overhead. A copied reproduction instruction named an absent helper; the retained helper name and original-input requirement were corrected, indexed and reviewed while preserving prior versions and raw bundles. Earlier DNS and export failures remain retained rather than being silently replaced with success-only records. The first actual warning audit also failed because a successful RNULL serialized as JSON{} while its reader expected null. Minimal serializer and four-receipt reproduction preceded a narrow reader fix; independent AST/controls review confirmed all downstream checks unchanged. A separate readback helper initially compared the string run ID with an integer; the assertion was corrected to the actual strict string schema. Original failures remain visible.

## 10. Known Residuals

Actual platform warning/skip review is complete within its separate-replay scope. Native timing adjudication, final panel reconciliation, final claimed-rung validator and maintainer approval remain open in this draft. The valid DOI still returns403 to the URL checker. Windows native totals exceed600seconds and Mac has one timingNOTE. These observations do not establish CRAN acceptance. Narrow-viewport code blocks require horizontal scrolling and long NEWS sidebar labels truncate. No exported browser screenshot files are claimed. No CRAN submission, confirmation, tag, announcement or PR236 merge has occurred.

## 11. Team Learning

Bind source commit, checker commit, deployed commit and frozen gzip independently. New guidance can change installed bytes even when runtime code stays the same. Warning totals are a reconciliation check; message identity and consequences require their own retained evidence. Brain and Codex memory writes were not authorized, so durable task receipts remain in the repository and task scratch only.

## 12. Cross-Product Coverage

This closeout covers exact artifact identity, reviewed source/site changes, installation independence and bounded termwise fixed-effect interoperability when the final tests pass. It does NOT cover joint covariance, loadings or variance-component pooling, new scientific default choices, general inferential validity outside the retained recovery regimes, live external-service availability from mock/cache tests, or CRAN acceptance. R-hub rebuilt-source checks do NOT cover original-gzip execution. Source CI does NOT cover frozen-artifact Windows checks.

Golden Set: existing shipped tests and prior bounded recovery fixtures are retained within their stated regimes; no new campaign. Memory receipt: prior release boundaries retrieved from MEMORY.md217-218 and checked against current repository/primary evidence. No memory writes.

Assessment: working draft NOTASSESSED. It is not a final verdict or submission document.
