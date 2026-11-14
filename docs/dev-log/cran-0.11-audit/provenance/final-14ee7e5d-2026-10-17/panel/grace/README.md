# Grace: frozen-artifact technical review

**READY for technical submission preparation within the scope below.** No concrete technical blocker was found for this frozen archive. This is not an upload authorization or an overall submission-ready declaration. The highest technical rung supported by this review is **platform-clean, with the local timing NOTE explicitly adjudicated below**. The next unproven rung is **submission-ready**: final evidence consolidation, all fresh panel verdicts, validator readback and the separate human timing/upload decision remain necessary.

Reviewed independently on 18 October 2026 UTC (17 October in Edmonton), in independent judgment mode. No previous READY vote was counted. Package source, archive, native checks, workflow settings and other reviewers' verdicts were not changed. No fits, full checks, merges, submissions or external messages were run.

## Archive and source

The reviewed gzip is `/private/tmp/pigauto-release-candidate/sha256-14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb/pigauto_0.11.0.tar.gz`, SHA-256 `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`, 5,137,648 bytes. Fresh readback confirmed 250 entries, 236 regular files, 14 directories, no links or unsafe paths, and the identical persistent retained copy. Independent `git archive` comparison with source `8e275ca738129c8d608f06ffc407b8141c1a2636` found 215 byte-identical source files, only DESCRIPTION normalization, and 20 expected build-generated files. Whitespace-normalized original DESCRIPTION fields agree; Author, Maintainer, NeedsCompilation and Packaged are generated additions. Clean pre-build source evidence is retained separately.

The archive contains no package-owned src, native binaries or excluded BACE/docs/script/dev/.github payload. The inherited tracked `_problems/test-new-features-526.R` is explicitly included, source-identical, and an isolated historical extraction rather than an executed test in this suite. Its content contains an old within-species expectation; it must not be cited as a current prediction contract. This is not a new archive delta or a reason to rebuild this candidate.

The workflow identity is separate: direct run 39082452558 used synthetic merge `e94e40a2f66a7f27dcc69b515d762d69c6eadf66`, whose parents and tree bind to checker head `6cd5f65b06d0c3ee32ce4b2d5a67a22c06028116`. Neither is substituted for the archive source. R-hub run 39081694892 is supplementary rebuilt-source evidence; it is not original-gzip execution proof. No Win-builder server-side digest binding is asserted.

## Native platform evidence

I opened all three remote raw check logs from the retained ZIPs and independently checked download/pre/post archive identities, exact commands, CRAN environment, 29/29 declared dependencies, incoming feasibility, install/load, examples, tests, vignette rebuilding and PDF/HTML manuals. Windows release/devel and Linux release have native exit 0 and Status OK. makeindex initially absent was repaired before the check; retained after-canaries succeed and manuals rebuild. Courier/V8 tool receipts are present. All four native suites report 3,195 passes, zero failures, 161 test warnings and 86 skips. These test warnings are distinct from native R CMD check WARNING status.

The Mac native check also has exit 0, but Status 1 NOTE. Its strict driver correctly exits 1 because it refuses to accept a NOTE without adjudication. That distinction is retained. The Mac's incoming feasibility, dependency/install/load checks and both manual formats pass. A separate CPU tensor probe succeeds; no GPU/MPS hardware or clean-room dependency compilation claim follows.

## Timing adjudication

| Native platform | Total seconds | Tests | Vignette rebuild |
|---|---:|---:|---:|
| Windows release | 645.813 | 415 seconds | 22 seconds |
| Windows devel | 622.432 | 404 seconds | 22 seconds |
| Linux release | 336.330 | 206/206 seconds | 11/11 seconds |
| Mac arm64 | 3564.464 | 5/31 minutes | 17/120 seconds |

Linux and Mac pairs are rounded CPU/elapsed observations. The sole Mac NOTE is `compare_methods`: user 2.532, system 0.992, elapsed 18.406 seconds, hence CPU 3.524 seconds. The broad CPU/elapsed differences support shared-host delay as an explanation, but do not prove its cause. I accept this sole disclosed NOTE for technical preparation: it shows elapsed-time cost, not a failed example, excess measured CPU for this example or a correctness failure. It stays in the release record and comments. No native run was repeated to replace it.

The two Windows totals exceed the protocol's 600-second observed incoming signal. This is an update from published 0.10.0, not a first submission; native statuses are OK and vignettes take 22 seconds. The current [CRAN policy](https://cran.r-project.org/web/packages/policies.html) asks for minimal CPU/check time and examples of a few seconds; it does not give an immutable 600-second total limit. This is a material submission risk for the maintainer to accept or decline explicitly, not proof of CRAN acceptance and not a technical fault waived by approval. The strict first-submission timing prohibition does not describe this release profile. Individual vignette diagnostics all exit 0; their 74.8865-second render sum and 308.2825-second R-process sum are separate measurements, not replacement native timings. Disabled teaching chunks remain disabled.

## Warnings and skipped coverage

The actual four-cell audit passes after the narrow successful-R-NULL serialization correction. I reviewed the fix, strict identity contract, preserved failed output, and independent 142-case/138-protected-file review, then rehashed the actual audit inputs. A fresh raw-RDS decode compared messages, test/source/line/caller identities and all 86 skip identities against that strict contract. Unknown warnings, changed callers and missing replay errors remain rejected.

The 161 warning identities comprise 119 small-validation diagnostics, 9 undefined constant correlations, 4 constant matched traits, 5 few-tree-draw diagnostics, 5 ignored joint-baseline covariate diagnostics, 3 GNN-off covariate diagnostics, 12 ordinal-imbalance diagnostics, 1 documented signal-gate fallback, 2 conformal sample-limit diagnostics and 1 deliberately unmatched-species diagnostic. These are fixture limitations with real scientific consequences; they do not certify coverage or broad inferential validity. Warning objects originate in separate native-layout replays on the checked installations. Native logs retain totals, not these saved replay objects. This disclosure is required and appears in the comments draft.

## Conditional gates and ancillary findings

This is an update, with zero direct reverse dependencies measured over all five dependency fields in the retained primary CRAN index. No reverse-package set exists to test. The change-impact boundary is preserved: all 49 R files and all six data streams/NOTICE match the tested predecessor; source documentation changes do not constitute a new recovery campaign. NeedsCompilation is no; sanitizer/native-registration gates for package-owned compiled code are inapplicable. Compiled dependencies still receive runtime checks.

Archive regular bytes are 7,150,500; compressed bytes are 5,137,648. Native installed-size INFO is 5.9 MB remotely and 7.0 MB locally, with the largest data directory reported as 3.2 MB locally. These observations are below current CRAN's general 5 MB data/documentation and preferred 10 MB source-tarball guidance; installed total alone is not a per-component violation. No size WARNING occurs.

Installed absence evidence gives 60 passes, zero failures/warnings and 8 expected skips without rgbif/terra. Cache/helper/failure paths are bounded proof; skipped live/stored-point/raster integrations and live GBIF/WorldClim availability remain unproven. All 12 optional-MI install/fixture/real steps in four dependency states exit 0; their logs were freshly rehashed. Both measured model versions are 0.7.1; real fixtures use three fits per backend and recorded dependency symlinks, not clean-room compilation or a research recovery design. Both backends are absent from all pigauto dependency fields.

Spelling review adjudicates 404 flags/1,886 contexts/43 source files. I independently re-bound all 43 files to the archive; no ordinary typo remains identified among those flags. This is not certification of unflagged prose or bibliography. The one distinct DOI resolver 403 remains, repeated in Getting Started source/HTML. Retained primary Crossref metadata matches Hackett's Science 2008 article and its DOI. Current [CRAN URL guidance](https://cran.r-project.org/web/packages/URL_checks.html) permits explanation of access-denial exceptions; the comments draft explicitly describes this finding without calling it fixed. I accept that disclosure for preparation, without predicting CRAN's decision.

## Rights, authority and remaining conditions

All six data objects and NOTICE match the carried-forward streams. The component ledger covers shipped data, compatibility/teaching fixtures, artwork and excluded development/comparator assets. Later rights records explicitly supersede the early ledger's overbroad MIT wording: the megatrees software licence is not a BirdTree data grant. The existing direct maintainer warranty with citation/credit remains the recorded inclusion basis. Sole-author/maintainer metadata is unchanged. This technical review gives no independent legal guarantee and asserts no newly published BirdTree redistribution licence.

Before any submission-ready claim, consolidate these completed receipts and this scoped verdict into the final durable wrapper/comments/conditional ledger, add the fresh Rose and Pat verdicts, and pass the exact claimed-rung validator. Historical pending snapshots should be labelled or superseded without being treated as missing raw proof. The draft comments accurately disclose native status versus test warnings, replay provenance, timing, DOI denial and conditional limits. The maintainer must separately decide the timing risk and approve the exact archive's CRAN upload/confirmation. Evidence PR236 merge, tagging and announcement remain unauthorized. None of those decisions may waive a newly found concrete technical, identity or rights blocker.

Review receipts: `independent-readback.json`, `raw-warning-readback.json`, and `verdict.json` in this directory. This report is a fresh Grace vote only; it cannot substitute for the whole panel or external CRAN states.
