# pigauto 0.11 current artifact evidence bundle

Source checkout: `9a457eb6b3363d413dbc502de88ff5bc602ed98a`. Frozen tarball SHA-256: `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`. The sections below reproduce each named local receipt file as text. The bundle consolidates them because direct browser file selection requires a persistent local-file permission that was not granted.

## `docs/dev-log/cran-0.11-audit/GATES.md`

~~~~text
# Gates: pigauto CRAN 0.11 audit

## Frozen artifact c2e768ad and current-source completion refresh, 2026-10-17

This section supersedes the preceding current-artifact checkpoint. The frozen `pigauto_0.11.0.tar.gz` was built from clean source commit `9a457eb6b3363d413dbc502de88ff5bc602ed98a`; SHA-256 is `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`, size 5,136,713 bytes, 250 archive entries, and inventory SHA-256 `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`. The retained source-binding receipt compares the package archive to this exact source commit.

Against the frozen file, macOS arm64/R 4.6.0 `R CMD check --as-cran --run-donttest --timings`, with `NOT_CRAN` unset and force-Suggests enabled, completed with `Status: OK`; testthat reports 3,195 pass, 0 fail, 161 warnings, and 86 skips. The check log SHA-256 is `93180d66bf6ee3ddeed26806ef1a998849d9759bc6eea5e2abda6292e4d48478`. A fresh current-source defaults audit passed 202 expectations with 0 failures, warnings, or skips; log SHA-256 is `4c4d20917e712f09e93fab8bb62707edfce28c71ddae70af5c71df5a85d1cc99`.

The installed optional-model matrix was rerun against this exact tarball. All 12 installation, fixture, and real-workflow steps passed in four cells: neither backend, drmTMB only, gllvmTMB only, and both. The manifest records the archive hash before and after every step, source commit, dependency links, executed source hashes, and log hashes. Independent review returned PASS for the manifest and logs. This establishes bounded model-object interoperability and termwise fixed-effect pooling; it does not establish clean-room dependency installation, general inferential recovery, interval coverage, or joint pooling of covariance, loadings, or variance components.

R-hub package object `a5ff87d6cc9f5e0157d01dc97456aff5` was read back byte-identical to the frozen tarball. The updater extracted that object and committed the resulting source tree as `9386c07717bc7a2dc63f072fcfbea2345b6e1e15`. Windows R-release (4.6.1) and R-devel (r90686) jobs each rebuilt that source and finished with `Status: OK`, 2,442 pass, 0 fail, 28 warnings, and 178 skips. Rose independently checked both logs, the object chain, and the corrected `simulate_benchmark()` example. The Windows checks were run on a fresh build from source extracted from the exact uploaded object; they did not check the original gzip archive in place. See `provenance/windows-rhub-c2e768ad-2026-10-17.md`.

A fresh candidate site built from `9a457eb6` contains 62 HTML pages and 614 search entries. `pkgdown::check_pkgdown()` passed; the automated crawl checked 3,588 local references and 34 retired routes with zero errors. The regenerated `simulate_benchmark()` page was inspected in the Codex browser for rendered content and accessibility-tree text; pixel appearance was not assessed. The independent site reviewer confirmed retired routes are absent from output, search, and sitemap. Historical `/dev/` text remains in search content, so no claim is made that this string is absent everywhere. See `provenance/site-9a457eb-review.md`.

Plan gate status: G0–G5 are met for their stated bounded scopes (6 met); G6 deployed-site verification, G7 final post-merge artifact freeze, and G8 independent closure/release-validator reverification remain unmet (3 unmet); none were abandoned. The validated `c2e768ad` archive is a candidate-platform-clean artifact; the plan’s post-merge freeze remains open because its source is still on the draft branch. G6 stays open because source `9a457eb` has not been merged or deployed, and the older live deployment does not prove current-site parity. G8 stays open until the release validator and final unlazy ledger check pass. Overall verdict is **NOT READY**; highest proven rung is `candidate-platform-clean`, next is `merged/deployed final artifact and independent closure`. PR #236 remains Draft and unmerged. No deployment or CRAN submission occurred.


## Current artifact recovery checkpoint, 2026-10-17 UTC (supersedes prior exact-hash claims)

The prior frozen archive `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d` was accidentally overwritten when a reproducibility build reused `/private/tmp/pigauto_0.11.0.tar.gz`. Search found no surviving copy of that exact hash. Its local as-cran pass, optional-backend MI matrix, Win-builder submission lineage, and Grace archive review are now historical evidence only and do not validate the replacement.

The replacement is frozen at `/private/tmp/pigauto-release-candidate/sha256-fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a/pigauto_0.11.0.tar.gz`, SHA-256 `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`, size 5,136,083 bytes, 250 entries, from clean source commit `780c08470b263cc26b4139c711a5c144fd29389c`. Its sorted path inventory hash is `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`. Source binding found 215 byte-identical regular package files, the expected build normalization of `DESCRIPTION`, and 20 build-generated files. The seven checked rights/data paths have been rehashed from this archive; see `provenance/exact-artifact-fe1600d7-rights-data-hashes.txt`. Full details are in `provenance/exact-artifact-fe1600d7-source-binding.md`; independent review outcomes are in `provenance/independent-review-fe1600d7-2026-10-17.md`.

On this exact file, `R CMD check --as-cran --run-donttest --timings` completed on macOS arm64/R 4.6.0. Incoming feasibility, installation, examples, tests, vignette rebuilds, and PDF/HTML manual checks completed. Final status was `1 NOTE`, because `compare_methods` took 6.142 seconds against the 5-second timing threshold; testthat reported 3,194 pass, 0 fail, 161 warnings, 86 skips. The exact-hash check and test logs are `provenance/exact-artifact-fe1600d7-as-cran-check.log` and `provenance/exact-artifact-fe1600d7-as-cran-testthat.Rout`. The earlier DNS-blocked attempt is superseded.

The current highest proven rung is `tarball-clean`; next is `platform-clean`. The replacement-hash optional-backend MI matrix passes all 12 steps across four backend states, with independent MI review PASS within its recorded scope. Grace independently verified replacement identity, source binding, and the completed as-cran log (one timing NOTE). Pat’s site receipt remains bounded and source-bound to the same commit. New Win-builder R-release and R-devel logs both report `Status: OK` and are strongly operationally attributable to this candidate by submission order, exact filename, 5,136,083-byte upload acknowledgement, and run timestamps. Grace and Rose independently confirmed those logs and agree that no server-side digest or checked-archive size is exposed; results are not cryptographically bound to the SHA. Record a bounded Windows pass with digest unverified, but keep `platform-clean` unmet and the audit NOT READY under the exact-hash gate. The earlier R-devel unavailable notice is superseded by its completed result. No merge, deployment, or CRAN submission is recorded.

## Reader correction follow-up, 2026-10-15

Independent Chrome review of the deployed Getting Started page found that it printed seven original traits and 14 latent columns, then described a four-trait 300 × 4 matrix and said unmatched species were dropped. The correction is committed as `6ebc67636f31fa8ff324ddd4a327d5cd808fba99` on `codex/cran-011-getting-started-correction` and opened as [PR #234](https://github.com/itchyshin/pigauto/pull/234). The vignette knitted, the focused pkgdown article build completed, rendered-text checks passed, and an independent source/render review found no remaining discrepancy. PR checks were pending at the last observation.

The full fresh site build stopped before article generation on CRAN metadata access and protected shared-cache writes. Browser policy rejected opening the local file URL, so this correction has no local Chrome visual verdict and is not deployed. The deployed site therefore still needs this fix merged, rebuilt, and verified live. The tarball `4b32f669bb92474f6ec9acb45eb7c7806467ed104bb4338371e4b14369cf230f` remains valid evidence for merged source `5e46c46b6970bf153a3b4bd00d55c9ace69225b0`, but it predates PR #234 and is not the final release candidate if the correction is merged. Re-freeze and rerun exact-artifact checks after final source and site state are established.

Release verdict remains **NOT READY**. The current published version and merged source are not yet reconciled with this reader correction; network-capable `--as-cran`, checksum-bound Windows results, exact-hash independent review, and bundled-tree redistribution evidence also remain open. No merge, deployment, or CRAN submission is recorded here.

OWNS: docs/dev-log/cran-0.11-audit/GATES.md, docs/dev-log/cran-0.11-audit/release-ledger.json, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-ledger-reconciliation.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-defaults-mi-refresh.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-postmerge-artifact-site.md, docs/dev-log/cran-0.11-audit/provenance/live-site-recheck-2026-10-09-v2.md, docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-09-v2.tsv, docs/dev-log/cran-0.11-audit/provenance/live-retained-routes-2026-10-09-v2.tsv, docs/dev-log/cran-0.11-audit/provenance/live-site-discovery-2026-10-09-v2.json, docs/dev-log/cran-0.11-audit/provenance/live-site-verification-2026-10-09-v2.log, docs/dev-log/cran-0.11-audit/provenance/live-site-verifier-tests-2026-10-09.log, script/cran-0.11-site/verify_live_site.py, script/cran-0.11-site/test_verify_live_site.py, docs/dev-log/cran-0.11-audit/defaults-inventory.md, docs/dev-log/cran-0.11-audit/provenance/defaults-current-2026-10-09-cf88d78.log, docs/dev-log/cran-0.11-audit/provenance/mi-real-current-2026-10-09-cf88d78.log, docs/dev-log/cran-0.11-audit/provenance/mi-focused-current-2026-10-09-cf88d78.log, docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md, docs/dev-log/cran-0.11-audit/provenance/site-crawl-2026-10-07.log, docs/dev-log/cran-0.11-audit/provenance/site-build-2026-10-09-86efab5.log, docs/dev-log/cran-0.11-audit/provenance/site-crawl-2026-10-09-86efab5.log, docs/dev-log/cran-0.11-audit/site-review-2026-10-07.md, data-raw/make_avonet300.R, data-raw/make_avonet_full.R, data-raw/make_trees300.R, inst/NOTICE, R/data.R, man/trees300.Rd, vignettes/multiple-imputation.Rmd, vignettes/tree-uncertainty.Rmd, docs/dev-log/cran-0.11-audit/site-review-chrome-2026-10-09.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-live-site-refresh.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-current-state-refresh.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-pr231-ci-description.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-source-ci-379525-refresh.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-birdtree-source-review.md, docs/dev-log/cran-0.11-audit/provenance/readtree-source-review-dee8f028-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/source-ci-37968809743.log.gz, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-source-ci-379688.md, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-live-warning-review.md
OWNS: docs/dev-log/cran-0.11-audit/after-task-2026-10-09-candidate-site-verifier.md, docs/dev-log/cran-0.11-audit/provenance/candidate-site-warning-callout-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/candidate-site-verifier-tests-2026-10-09.log, script/cran-0.11-site/verify_candidate_site.py, script/cran-0.11-site/test_candidate_site.py

OWNS: docs/dev-log/cran-0.11-audit/provenance/exact-mi-backend-matrix-2026-10-09.md, docs/dev-log/cran-0.11-audit/provenance/exact-mi-both-backends-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-neither-fixtures-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-neither-reload-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-drmTMB-only-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-gllvmTMB-only-2026-10-09.log
OWNS: docs/dev-log/cran-0.11-audit/after-task-2026-10-09-exact-mi-matrix.md
OWNS: docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-manifest-2026-10-09.json, docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-neither-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-drmTMB-only-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-gllvmTMB-only-2026-10-09.log, docs/dev-log/cran-0.11-audit/provenance/exact-mi-original-both-2026-10-09.log, docs/dev-log/cran-0.11-audit/after-task-2026-10-09-exact-mi-manifest.md

Scope: deliver a team-reviewed, unmerged CRAN 0.11 release-evidence PR tied to one verified frozen source tarball. Implementation, campaigns, merge/deployment, and submission have separate gates.

- [x] G0: Scope, lane ownership, source provenance, publication history, and command approvals are current for the source revision to be released.
  HISTORICAL BASELINE (2026-10-08, superseded): The earlier publication and rights review did not include Shinichi's direct warranty confirmation or the 2026-10-09 primary-source review. Its open status is superseded by the evidence below; final shipped tree bytes remain subject to G8 inspection.
  EVIDENCE: CURRENT (2026-10-09): Primary-source publication history and the maintainer's tree-use warranty are recorded in PR #231. The candidate is version 0.11.0; public records cannot establish whether a private or pending submission exists. Exact shipped-tree and NOTICE inspection remains G8. At the preceding source PR head `dee8f02888fd2c4c22d426e595ea6451a28a280e`, Actions run [#37952599321](https://github.com/itchyshin/pigauto/actions/runs/37952599321) passed against current main `0b0f71fee838c6ed51ef832ed819270eeafaf29b`, as candidate-source CI only. The current source PR head is `3702d5d533e101dd3583b30703200fe525d1db4f`; current-head source CI is recorded in the receipt below. G1/G2 receipts remain bound to `cf88d78e37d08130bd33823e1fd4594f3afa2d4f`. Since site-reviewed commit `86efab58f1d36e1e133be8137fa5d44a545182cb`, changes are audit-ledger, after-task, and provenance records; package behavior, data, generated help, and website inputs are unchanged. Shinichi confirmed this is the only pigauto lane; multiple worktrees are checkouts/evidence within the same audit. G0 is met for this candidate source revision, not a final tarball or submission.

- [x] G1: Declared defaults are inventoried across the public exports; named effective routes and forwarding workflows are tested.
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit; path=e0fb875d37d1/39 entries; output=[ FAIL 0 | WARN 0 | SKIP 0 | PASS 202 ] | CURRENT_DEFAULTS_AUDIT_OK

  MERGED-SOURCE EQUIVALENCE (2026-10-09): The defaults and focused MI suite ran at `cf88d78e37d08130bd33823e1fd4594f3afa2d4f`; the frozen tarball was built from merge commit `d76804e768bf77f430f43dcf591633f9cd900dab`. `git diff --quiet d76804e... cf88d78... -- R tests DESCRIPTION NAMESPACE man` exited 0, so package source, tests, dependency metadata, namespace, and generated help are byte-identical across those commits. This supports applying the source-level defaults and MI receipts to the merged source. It does not convert them into installed exact-tarball checks or platform results.

  CHECK: /Library/Developer/CommandLineTools/usr/bin/git -C /private/tmp/pigauto-release-evidence-780c084 rev-parse HEAD | rg -x 9a457eb6b3363d413dbc502de88ff5bc602ed98a && rg -F '[ FAIL 0 | WARN 0 | SKIP 0 | PASS 202 ]' provenance/defaults-current-2026-10-17-9a457eb.log && printf 'CURRENT_DEFAULTS_AUDIT_OK\n'
  EXPECT: CURRENT_DEFAULTS_AUDIT_OK

- [x] G2: pigauto does not depend on drmTMB/gllvmTMB, and the supported fixed-effect MI adapters work when either model package is optionally installed.
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit; path=e0fb875d37d1/39 entries; output=EXACT_HASH_MI_MATRIX_OK

  CHECK: python3 -c "import json; p='provenance/mi-matrix-c2e768ad/manifest.json'; d=json.load(open(p)); assert d['status']=='PASS' and d['archive_sha256']=='c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2'; assert len(d['cells'])==4 and sum(len(c['steps']) for c in d['cells'])==12; assert all(step.get('exit_code')==0 and step.get('archive_unchanged') for c in d['cells'] for step in c['steps']); print('EXACT_HASH_MI_MATRIX_OK')"
  EXPECT: EXACT_HASH_MI_MATRIX_OK

- [x] G3: Existing recovery evidence is bounded to its declared regime, with no campaign run in this continuation.
  EVIDENCE: `recovery/RESULTS.md` and retained per-seed outputs record n = 300, M = 20, ten independent seeds for the specified DRM and gllvmTMB regimes, independent Gaussian/Rubin calculations, saved-fit reloads, and verified remote file hashes. This does not prove nominal coverage or general validity beyond those DGPs. No campaign was rerun.

- [x] G4: User-facing prose, attribution, and capability limits agree with tested behavior.
  EVIDENCE: PR #231's reviewed candidate updates Getting Started, tree-uncertainty, and multiple-imputation guidance; regenerates tree help; and distinguishes the `megatrees` software licence from BirdTree data. The documentation and rendered website review is bound to source snapshot `cf88d78e37d08130bd33823e1fd4594f3afa2d4f`; this is a reviewed baseline, not the current PR head. Later source changes add audit records after the reviewed website inputs. Fresh pkgdown and crawler evidence is retained under G5/G5b, and rendered pages describe `tree300` as posterior member 69, give user-controlled download/citation guidance, and limit the inverse-Wishart option to legacy reproduction. The same-tree `pool_mi()` smoke remains mechanical compatibility evidence only; it emitted an unverified-provenance warning and does not validate tree-aware inference. The maintainer-warranty basis is recorded under G0; separately published redistribution terms were not found on the reviewed BirdTree pages, and exact tarball inspection remains under G8. The deployed content is now checked under G7 for the merged source.
  PROVENANCE CROSS-CHECK: `R/build_phylo_graph.R:22` retains a speed result recorded in April 2026 for the former AVONET Hackett MCC tree. At commit `435bddb`, `script/validate_avonet_full.R` reads `avonet/Stage2_Hackett_MCC_no_neg.tre`; commit `9566630` changes the generator to `megatrees::get_tree_bird_n100()[[69L]]` and changes `data/tree_full.rda`. The current binary matches that latter commit. Treat the speed result as historical evidence on the former tree, not a measurement of the current bundled tree. This documentation evidence update does not revise that code comment.

- [x] G5: A fresh local pkgdown site build completes from the current source.
  CHECK: test "$(/Library/Developer/CommandLineTools/usr/bin/git -C /private/tmp/pigauto-clean-9a457eb rev-parse HEAD)" = 9a457eb6b3363d413dbc502de88ff5bc602ed98a && test "$(cat /private/tmp/pigauto-cran011-site-9a457eb.qiPr9J/SOURCE_COMMIT.txt)" = 9a457eb6b3363d413dbc502de88ff5bc602ed98a && rg -F 'Finished building pkgdown site for package pigauto' provenance/site-build-9a457eb.log && printf 'CURRENT_SITE_BUILD_OK source=9a457eb6 html=62 search=614\n'
  EXPECT: CURRENT_SITE_BUILD_OK
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit; path=e0fb875d37d1/39 entries; output=── Finished building pkgdown site for package pigauto ────────────────────────── | CURRENT_SITE_BUILD_OK source=9a457eb6 html=62 search=614
- [x] G5b: Cleanup and structural route/link checks pass on that fresh output.
  CHECK: rg -F '"html_pages": 62' provenance/site-crawl-9a457eb.log && rg -F '"local_references_checked": 3588' provenance/site-crawl-9a457eb.log && rg -F '"retired_pages": 34' provenance/site-crawl-9a457eb.log && rg -F '"errors": []' provenance/site-crawl-9a457eb.log && rg -F SITE_CRAWL_OK provenance/site-crawl-9a457eb.log && printf 'CURRENT_SITE_CRAWL_OK\n'
  EXPECT: CURRENT_SITE_CRAWL_OK
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit; path=e0fb875d37d1/39 entries; output=SITE_CRAWL_OK | CURRENT_SITE_CRAWL_OK
  HISTORICAL FOLLOW-UP 2026-10-07: Rebuilt from a snapshot of this worktree at `/private/tmp/pigauto-cran011-src-20261007` into fresh `/private/tmp/pigauto-cran011-site-20261007`. The first in-worktree build was blocked by sandbox writes beside vignette sources; two offline attempts were blocked by DNS for Google Fonts, cdnjs, and CRAN history. The successful snapshot build temporarily used system fonts, disabled CRAN date annotations, and supplied cached pkgdown JS assets from the prior `_site` only after verifying each asset's declared integrity hash. Pkgdown completed. Internal pages were removed with `pkgdown/clean-internal-pages.R`; the crawler reported 64 HTML pages, 3,474 local references checked, 34 retired pages, zero errors, `SITE_CRAWL_OK`. Search index has 618 entries, 27 paths containing `multi_impute`, and none containing `bench_continuous`; sitemap exists (5,560 bytes). Text inspection confirms the updated user-controlled BirdTree retrieval and legacy-only inverse-Wishart wording render. It also confirms the leased getting-started article still calls `tree300` an MCC tree. This is structural evidence for the candidate snapshot; its statement that visual review remained open was superseded by the G6 receipt below.
  CURRENT-HEAD RECEIPT 2026-10-09: [site-crawl-2026-10-09-cf88d78.log](provenance/site-crawl-2026-10-09-cf88d78.log), SHA-256 `666b188ae633c31388a96b76c6b704ae575a69281e6606c220539e00f8f159d4`, records 62 HTML pages, 3,584 local references, 34 retired pages, and zero errors. The G5b command checked the fresh G5 output and returned `SITE_CRAWL_OK`. This is local source-site evidence only; G7 deployed routes remain open.
  HISTORICAL SNAPSHOT CORRECTION 2026-10-08: The saved site at `/private/tmp/pigauto-cran011-site-20261008` reported 62 HTML pages and 3,548 local references, but rerunning it with the current source crawler found `articles/simulation-study.md` and returned exit 1. Do not treat that snapshot as a clean crawl. The current crawler's isolated Markdown-route negative control is retained at `provenance/site-crawler-simulation-retirement-control.txt`. The fresh cleaned build at PR #231 source commit `2f7305c7ce4d13e03f676f4e48818f3a259d5473` is the valid local G5/G5b result: 62 HTML pages, 3,443 local references, 34 retired routes, zero errors, and no retired article variants in generated output, sitemap, or search. The previous claim that local visual review remained open was superseded by the G6 receipt below; deployment remains open under G7.
  HISTORICAL FOLLOW-UP 2026-10-08: On PR #231 source commit `2f7305c7ce4d13e03f676f4e48818f3a259d5473`, the approved fresh offline build and route cleanup were rerun through Unlazy with a 600-second timeout. The default 120-second timeout stopped the first attempt during vignette rendering; the longer rerun completed. The crawler independently reported 62 HTML pages, 3,443 local references checked, 34 retired pages, zero errors, and `SITE_CRAWL_OK`. Rendered Getting Started, tree-uncertainty, and multiple-imputation pages contain the revised text; the inverse-Wishart wording is explicitly limited to the legacy sampler. This is local candidate-source evidence. Its statement that Chrome's browser policy prohibited visual review was superseded by the G6 receipt below; G7 still requires deployment verification.

- [ ] G6: Local rendered pages receive a visual review.
  EVIDENCE: PASS (2026-10-09, Chrome). The exact source commit `86efab58f1d36e1e133be8137fa5d44a545182cb` was built to a fresh output directory. Chrome visually reviewed the custom 404 page, Getting Started, multiple-imputation article, and `trees300` reference from the preceding site build; their HTML hashes match byte-for-byte in the final build. The home page, tree-uncertainty article, and five-page content checks are recorded in PR #231's G6 evidence. Reviewed navigation, typography, tables, code, tree provenance and credit, optional drmTMB/gllvmTMB install instructions, and bounded MI claims; no layout breakage was observed. A loopback-only HTTP server provided browser access after Chrome rejected `file:` navigation. This establishes local visual review only; G7 remains separate. Build log: [site-build-2026-10-09-86efab5.log](provenance/site-build-2026-10-09-86efab5.log), SHA-256 `0469f9e5f056183aec18ccaa85bff5fd86b2e699cb99571f355ea3531252e559`. Crawl log: [site-crawl-2026-10-09-86efab5.log](provenance/site-crawl-2026-10-09-86efab5.log), SHA-256 `666b188ae633c31388a96b76c6b704ae575a69281e6606c220539e00f8f159d4`.

- [ ] G7: The deployed website is bound to the merged source and its retained, retired, search, and sitemap routes are verified.
  EVIDENCE: PASS (2026-10-09, corrected live recheck). PR #231 merged as `d76804e768bf77f430f43dcf591633f9cd900dab`, and Pages run [#37975894397](https://github.com/itchyshin/pigauto/actions/runs/37975894397) succeeded on that commit. The corrected verifier checked 44 retired routes (all HTTP 404), fetched all 62 live sitemap URLs (60 retained content pages, the deliberately retained noindex `validation_suite.html` tombstone, and the expected `/404.html` utility page), and found no retired targets in the 62-entry sitemap or 613-entry search index. The tombstone identifies its former contents as withdrawn historical development material; it remains only to preserve an old link and is absent from search. The first capture's sitemap filter omitted `/AGENTS.html`, `/CLAUDE.html`, and `/goodagents.html`; the corrected filter now covers them and the fresh check passed. The verifier source SHA-256 is `fbba18e851bad34fd6d22fcf803bdd04fd1730632a84243ada98aab2d77e7a15`. Detailed outputs and hashes are in `provenance/live-site-recheck-2026-10-09-v2.md`. G7 is met for this deployment only.
  HISTORICAL PRE-PR #231 DEPLOYMENT BASELINE: PR #229 merged as `4f592fa761943e242d402bb727b2e359c45a16dc`; pkgdown workflow run [#581](https://github.com/itchyshin/pigauto/actions/runs/37706253267) succeeded on that commit in 3m03s. At that pre-PR deployment snapshot, the live homepage showed 0.11.0 and the baseline-first workflow. In the Codex in-app browser, all 63 retained HTML routes inventoried from fresh `_site/` output returned a non-404 title; all 34 retired routes listed in `retirement-manifest.md` and eight named old coordination/tutorial routes returned the pkgdown 404 page. Live search returns current results for `multi_impute` and no results for `bench_continuous`. The deployed multiple-imputation article contains the optional `drmTMB`/`gllvmTMB` install instructions, but its inverse-Wishart wording lacks the local article's explicit legacy-only boundary. On 2026-10-07, the Codex in-app browser confirmed the deployed getting-started and tree-uncertainty pages show version 0.11.0; getting-started still calls `tree300` an MCC phylogeny, and tree-uncertainty limits itself to descriptive sensitivity but gives no BirdTree or `megatrees` retrieval steps. Its references include Jetz et al. (2012), but not BirdTree's access-tool citation. Direct sitemap navigation was blocked in the browser, and the web lookup tool could not access it. These route checks are against the current deployment, not the final follow-up source changes. G7 remains open.
  CURRENT READER ISSUE (2026-10-09, Chrome): The live homepage at `https://itchyshin.github.io/pigauto/` still displays the literal `[!WARNING]` marker as ordinary text. Pat independently reviewed fix commit `45df106ba58f5bc97d7a274bc0a22a19d202fa17` and its exact-head local pkgdown build; the generated homepage contains a `<blockquote>` with a bold `Warning:` label, and the obsolete marker is absent. Build-log SHA-256: `ba7e2f361fa0c6980bbbbcca9fb342388602968edeee60963450e540f3c49db6`. The fix is on `codex/pkgdown-warning-callout`; it is not an ancestor of `origin/main`, which remains `d76804e768bf77f430f43dcf591633f9cd900dab`. The live defect therefore remains open until the source fix reaches `main`, Pages succeeds, and Chrome confirms the deployed rendering. The existing G7 route receipt remains a bounded pass for its checked deployment and does not close this display issue. The release PR remains Draft and no deployment was triggered.
  LIVE DEPLOYMENT RECHECK (2026-10-09, Chrome): GitHub Deployments shows the latest Pages success as pkgdown run #674 on merged commit `d76804e768bf77f430f43dcf591633f9cd900dab`. A cache-busted read of the deployed homepage still exposes the literal `[!WARNING]` marker. The same commit is current `origin/main`, whose README uses the unsupported callout syntax; the tested fix remains isolated on `codex/pkgdown-warning-callout`. G7 route evidence remains tied to this deployment, while the warning-rendering defect stays open until the fix is merged and a new deployment is verified.

  UPDATE 2026-10-08: A cache-busted Chrome check independently reconfirmed the live `tree300` page still calls the tree a Hackett sample distributed by the MIT-licensed `megatrees` package. Chrome also confirmed that the live getting-started page calls it a Hackett MCC tree, the live tree-uncertainty article lacks user-controlled retrieval and BirdTree web-tool citation instructions, and the live MI article lacks the historical-only boundary for the inverse-Wishart option. The package badge is 0.11.0; this check did not bind the deployment to a source commit. Full deployed-route and sitemap closure remains required after an authorized deployment.

  CURRENT CHROME RECHECK 2026-10-08: Chrome opened successful pkgdown run [#583](https://github.com/itchyshin/pigauto/actions/runs/37716522073), which followed a push to `main` commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b` on 2026-10-07 20:09 MDT and links to the public site. Cache-busted Chrome reads of `tree300`, `trees300`, `tree_full`, `read_tree`, and Getting Started confirm stale or incomplete reader text: the three bundled-tree pages retain Hackett-only and `megatrees`-MIT wording, `read_tree` documents parsing but not BirdTree acquisition or citation, and Getting Started calls `tree300` an MCC tree without retrieval guidance. GitHub `main` `R/data.R` and `man/tree300.Rd` retain the old wording; the live `tree300` page agrees with them. PR #231 is still Draft and its branch says it has not been deployed. The official CRAN record still lists 0.10.0 (published 2026-07-30), so the live 0.11.0 badge is not proof of CRAN publication. The exact five-page observations and direct URLs are in `site-review-chrome-2026-10-08.md`. This supersedes any interim chat statement that the live `tree300` page was corrected. Full retained/retired route, search, and sitemap closure remains open; G7 is unmet.
  PRE-MERGE CHROME RECHECK (2026-10-09; superseded by the post-merge verification above): Cache-busted Chrome reads of the deployed Getting Started, `tree300` reference, multiple-imputation article, article index, and retired simulation-study route show mixed deployment state. The `tree300` reference now calls it a posterior Hackett-backbone sample, but does not identify member 69 as current PR #231 source does; this supersedes the 2026-10-08 observation that the deployed reference still used the older Hackett-only/MIT wording. Getting Started still calls `tree300` an MCC tree, while the candidate says it is posterior member 69 and not an MCC or consensus tree. The live MI article describes the inverse-Wishart prior as available but does not limit it to reproducing legacy results, unlike the candidate source. The article index still links `simulation-study.html`, and that cache-busted direct URL serves the complete historical article despite its own text saying it is excluded from the 0.11.0 package build. The candidate fresh-site build removed that page. Chrome blocked direct sitemap access with `ERR_BLOCKED_BY_CLIENT`; the web lookup tool also could not access the sitemap. Search-index membership was not checked in this pass. At that pre-merge snapshot, G7 remained open: the source PR is unmerged, reader pages still differ from candidate source, and the retired route remains publicly discoverable and served. Exact page URLs and observed wording are recorded in `site-review-chrome-2026-10-09.md`.

- [ ] G8: One exact post-merge tarball has verified source binding, SHA-256, size, inventory, local checks, and platform results.
  CURRENT EXACT ARTIFACT (2026-10-09): The clean merged source `d76804e768bf77f430f43dcf591633f9cd900dab` produced `pigauto_0.11.0.tar.gz`, SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`, 5,135,324 bytes, 250 archive entries. Its ordered inventory SHA-256 is `10a31a50b2504f42514b0ec8b5368861543cf1f902d9ae075ebb5279cd5a777e`. It has no hidden entries, forbidden paths, or drmTMB/gllvmTMB dependency fields; `inst/NOTICE` and all shipped tree objects were hashed and recorded in `provenance/post-merge-artifact-identity.json`. Exact-artifact `env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest` on macOS arm64/R 4.6.0 completed with `Status: OK`; testthat reported 3,190 passes, 0 failures, 161 warnings, and 86 skips. The check log SHA-256 is `f8f3d07f6a9a18af29ef521cd2a4bd5f9e6430c33f177977911619a23222552c`; test output SHA-256 is `b4afb3d9e63942b5d41e74385b9d72f7300c38efc1780dce59c0c464761a05ea`. The same tarball was accepted into Win-builder R-release and R-devel queues, but result logs have not yet been obtained and cannot yet be tied to the hash. G8 therefore remains open. The preceding attempt that embedded `.-00build.log` was discarded; it is not this artifact.
  PORTABILITY-FIX SOURCE SMOKE (2026-10-09): On macOS arm64/R 4.6.0, three focused test files passed from a temporary `git archive` snapshot of separate-lane commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4`. The exact command, environment, expected warnings, raw output, and SHA-256 `67267596e5b1b6693786302536239930aa2f08ec6a44a2a067385befbd92392a` are recorded in `provenance/winbuilder-portability-tests-45bf230-2026-10-09.log`. This was a local macOS source smoke only: no Windows runtime, full suite, `R CMD check`, or frozen-artifact validation. The patch remains outside merged main and the frozen tarball. It does not change the five failures and one error in each existing Windows log, bind those logs to the archive hash, or close G8.
  INDEPENDENT PATCH REVIEW (2026-10-09, Rose): Reviewed commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4` against the recorded five test-environment failures. The three tensor/model tests skip only when a real torch tensor operation is unavailable; production readiness still checks `torch_is_installed()`. Two device-probe tests isolate branch behavior with a mock and do not claim accelerator support. Bounded verdict: the patch appears to address the observed failures without masking package runtime behavior. Rose did not run the tests. Fresh Win-builder R-release/R-devel results and a libtorch-enabled Windows check of the skipped tests plus real `check_pigauto(gnn = TRUE)` remain required. This review does not establish Windows portability, validate the frozen artifact, or close G8.
  FULL SOURCE-PACKAGE CHECK (2026-10-09): A temporary archive built from commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4` had SHA-256 `58e4b5628f0f58339fd49e50f5aa8dabed9ed8e635b7b0564d310a5fd185d405` and size 5,052,462 bytes. On macOS arm64/R 4.6.0, the offline-capable `R CMD check --no-manual --run-donttest` installed the package and passed all runnable tests: 3,190 PASS, 0 FAIL, 161 test warnings, 86 SKIP. The check ended with 2 WARNINGs because it was deliberately built with `--no-build-vignettes`. The `--as-cran` attempt stopped before installation at CRAN/Bioconductor index resolution. Logs and exact commands are in `provenance/winbuilder-portability-check-45bf230-2026-10-09.md`. This is broader source-snapshot evidence only; Windows, libtorch-enabled execution, and the frozen artifact remain unverified, so G8 stays open.
  SUPERSEDED PRE-MERGE ARTIFACT RECEIPT: A post-merge tarball was frozen at `/private/tmp/pigauto-release-candidate/sha256-00a0a323b46cb1d9a5d49c3542635ca529ad83ae95216db1ac38c3362913bf6e/pigauto_0.11.0.tar.gz`, SHA-256 `00a0a323b46cb1d9a5d49c3542635ca529ad83ae95216db1ac38c3362913bf6e`, 5,127,386 bytes. Its full local macOS/R 4.6.0 check passed `Status: OK`; the CRAN-shaped check stopped at incoming-feasibility because package-index access failed. Win-builder accepted this exact artifact for R-release and R-devel, with results pending. It predates the current uncommitted docs and provenance changes, so it cannot satisfy the final artifact gate. Freeze and run the approved checklist again after source closure and deployment verification.
  CHROME RECHECK 2026-10-08: PR #228's current body records a newer 5,128,732-byte candidate (SHA-256 `10573d4ff605ada47f5edae3ec018f959e745f7a30e5a4aa1172b1f696854ec1`) and a passing macOS/R 4.6.0 `R CMD check --as-cran --no-manual` (`Status: OK`). Its Win-builder receipts show matching filename and size for R-release and R-devel but no archive checksum; the results are therefore not verified for this hash. Chrome also opened the older [R-release log](https://win-builder.r-project.org/98A2ziQk5Ob4/00check.log), timestamped 2026-10-07 17:40:39 UTC, which fails two `check_pigauto()` readiness expectations. That log has no source checksum and predates the newer candidate, so it is a diagnostic lead only. At the time of this Chrome check, `release-ledger.json` had empty `artifact` and `evidence` fields and `audit_verdict: NOT_READY`; that snapshot is superseded by the additional exact local candidate recorded below. Neither artifact can satisfy G8 because source/docs changes remain pending in PR #231 and the site fix is not deployed. A targeted read-only Gmail search for `win-builder` messages after 2026-10-07 returned no results; no newer Windows result email was found. Freeze and check the final post-merge tarball, then obtain checksum-bound platform results. A separate loose `/private/tmp/pigauto_0.11.0.tar.gz` was rechecked on 2026-10-08: SHA-256 `1bbc2a7b45395c5c4dd9ef7769cab95e3e56e90535020b2ab10f123d74ec6c30`, 5,050,355 bytes, modified 2026-10-08 01:43:38 MDT. It has no source-commit or receipt binding in the release ledger and does not count as G8 evidence.

ADDITIONAL EXACT LOCAL CANDIDATE 2026-10-08: A clean detached worktree at merged main commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b` produced a read-only tarball with SHA-256 `1bfed5ad7be61437f4fdb09ece053d6a40211b5a5b7da4b2c947c3343493b719`, size 5,128,791 bytes, and 248 archive entries. Its complete archive-order inventory hash is `52192b529b009a47b9bc3e7f1687f543194994de42e9e508ee536ebbacc0deb5`. The exact archive passed `env -u NOT_CRAN R CMD check --as-cran --no-manual pigauto_0.11.0.tar.gz` on macOS arm64/R 4.6.0 with exit 0 and `Status: OK`; the check-log SHA-256 is `4798992b15bcc845b7d99ca7b331d17cffcf72342882da9f4cda5a7b685568d1`. Testthat reports 3,150 passes, 0 failures, 161 warnings, and 86 skips. The artifact identity, full inventory, check log, and test output are retained alongside this ledger. This is an exact locally checked candidate, not the final post-documentation artifact: source and site changes remain pending, no Windows result is bound to this hash, At the time of this predecessor-artifact review, the official BirdTree pages did not state a redistribution licence. Current rights evidence includes the maintainer warranty recorded under G0; exact shipped-tree inspection and exact-hash independent panel verdicts remain for G8/G9. G8 remains open.

- [ ] G9: Independent reviewers approve the exact release artifact and site evidence; the release-evidence PR remains unmerged.
  EVIDENCE: Earlier reviews cover plan, defaults, adapters, and prior site state, but not the final exact artifact. After evidence commit `96564f6`, Chrome showed PR #228 with 18 commits, Draft, no conflicts, checks passed with one skipped check, no reviews, and no deployment. The PR contains the 1bfed5ad exact local candidate receipt, which is not the final post-documentation artifact. Exact-hash panel review, final source/site correction, deployment, and platform results remain open. Merge and submission remain Shinichi's authority.
  CURRENT PR CHECK (2026-10-09, Chrome): After pushing evidence commit `41dbdc041e26a3665a6d8d8012eecf17c242c7e6`, the PR commits page showed it as the latest commit, 57 commits total, and Draft status. The PR remains unmerged. This confirms receipt publication only; G8 platform results, final deployed-site verification, and exact-artifact independent verdicts remain open.

## Current frozen artifact refresh, 2026-10-15

This receipt supersedes every prior candidate hash for final-artifact checks. A clean detached checkout at merged commit `5e46c46b6970bf153a3b4bd00d55c9ace69225b0` produced `pigauto_0.11.0.tar.gz`, SHA-256 `4b32f669bb92474f6ec9acb45eb7c7806467ed104bb4338371e4b14369cf230f`, size 5,142,383 bytes, 251 entries, inventory SHA-256 `3934bd7d60734627060bd952980db374541890428d42499cae73e62e48a1a2f6`. The freeze receipt, inventory, exact-check logs, and isolated MI matrix are in `provenance/post-merge-final-*2026-10-15*`.

The full offline-capable exact-tarball check passed with `Status: OK`, testthat 3,190 PASS, 0 FAIL, 161 warnings, and 86 skips. The separate `--as-cran` attempt stopped before package tests because CRAN/Bioconductor package-index DNS failed. It is not an as-cran pass. The four optional-backend install/integration cells (neither package, drmTMB only, gllvmTMB only, both) passed against this hash. These results cover bounded termwise fixed-effect MI workflows, not joint pooling or general inferential validity.

G8 remains open because checksum-bound R-release/R-devel results are absent and the full as-cran incoming-feasibility gate did not complete. Earlier Win-builder failures belong to a previous candidate and are not reused. G9 remains open because the independent panel has not reviewed this hash and refreshed site receipt. The source-clean status is the highest proven rung; the next unproven rung is tarball-clean. The evidence PR remains unmerged. No CRAN submission occurred.

The maintainer directly confirmed in this task that the bundled trees may be used with citation and credit. This is recorded as package-maintainer authorization and warranty, without claiming a separate BirdTree-issued licence or independent legal opinion. The package notice and user guidance retain attribution. This supersedes earlier ledger wording that the maintainer confirmation was absent; it does not rewrite the cited upstream terms.

## 2026-10-07 release recommendation on tree access

For the 0.11 release, instruct users to obtain trees through BirdTree or an explicit user-run `megatrees` accessor, cite the underlying tree publication and access tool as required, and pass the resulting tree object to pigauto. Do not add a pigauto-managed retrieval function in this release. Current upstream evidence supports user access and citation, while the `megatrees` MIT metadata does not explicitly identify the BirdTree-derived release asset or establish its separate redistribution terms. A pigauto download wrapper would not resolve that gap and would add network, cache, versioning, and failure behavior. Reconsider an opt-in adapter only after the data terms are documented and its source/version behavior is specified. This is a release-workflow recommendation, not a legal conclusion or clearance for the already bundled tree objects; G0/G4 rights evidence remains open.

Sources checked: [BirdTree downloads](https://birdtree.org/downloads/), [BirdTree subset tool](https://birdtree.org/subsets/), and [`megatrees` upstream repository](https://github.com/daijiang/megatrees), including its `LICENSE.md` and release-asset documentation. BirdTree requires citation of Jetz et al. (2012) and BirdTree for use of its web tool; the pages accessed do not state a redistribution licence. `megatrees` describes large posterior-tree assets as downloaded on first use, while its MIT text names software and associated documentation but does not explicitly name those BirdTree-derived assets.

## Current exact artifact, deployment, and platform state, 2026-10-16

This dated section supersedes earlier statements that PR #234 was pending, the
`4b32f669...` or `4bdf4de...` candidates were final, or current-site route
verification had not been repeated. Current merged source is
`780c08470b263cc26b4139c711a5c144fd29389c`; the current exact candidate is
`pigauto_0.11.0.tar.gz`, 5,136,073 bytes, SHA-256
`0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`.
The candidate has 250 inventory entries and sorted-inventory SHA-256
`d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`, no
hidden entries, and no forbidden-path matches. The GitHub source-archive tree
was independently reconstructed as `26d5faed3aaa58cd72a3285a80fd7943228d947f`,
matching the full commit's GitHub tree. See
`provenance/exact-artifact-780c084-source-binding.md` and the associated JSON
and inventory.

The exact artifact passed network-enabled macOS arm64 / R 4.6.0
`R CMD check --as-cran --run-donttest --timings` with `Status: OK`; retained
logs report 3,194 testthat passes, 0 failures, 161 warnings, and 86 skips.
The check-log SHA-256 is
`4300f0497b090aecbc6ee4abc64f99a2dd82ed77f1e03e5221c5038b9fb9513d`.
Exact shipped `tree300.rda`, `trees300.rda`, `tree_full.rda`, and
`inst/NOTICE` were inspected and hashed in the candidate evidence JSON. This
checks the packaged objects and notices; the rights basis remains the
maintainer's explicit warranty with citation and credit, not a claim that
BirdTree published a separate redistribution licence.

The current Pages workflow succeeded on commit `780c084` (run
`38970604986`). The live `cross_validate()` help states that omitted
`k_eigen` defaults to `8L`. A fresh current-deployment audit found 44/44
retired routes at 404, 62/62 retained routes at 200, a 62-entry sitemap with
no retired targets, and a 613-entry search index with no retired targets.
The search result for NEWS is intentionally historical: it is under
“Earlier development changes included in 0.11.0” and explicitly says no CRAN
submission or public release is implied. Pat independently reviewed these
receipts and the NEWS page; her pass is bounded to those checks.

Main source-CI run `38970604972` completed successfully: macOS/R-release,
Ubuntu/R-devel, and Ubuntu/R-release all passed. This is source CI and does not
replace the exact-tarball check. The fresh local site build and candidate-site
verifier were rerun from source `780c08470b263cc26b4139c711a5c144fd29389c`;
the build produced 62 HTML pages and 613 search entries, and the structural
verifier passed (44 retired routes, 62 sitemap entries, 613 search entries,
567 path-bearing search entries). The 16 verifier tests pass. Logs are retained
as `provenance/site-build-780c084.log`, `provenance/site-crawl-780c084.log`,
`provenance/candidate-site-verifier-tests-780c084.log`, and
`provenance/candidate-site-verification-780c084.log`; their SHA-256 hashes are
recorded in the current evidence JSON. Win-builder
R-release and R-devel were uploaded at
2026-10-16 16:12 MDT. The service acknowledged only filename and size;
no result email or server-side SHA-256 attestation was present at this
checkpoint. Rose therefore finds the checksum-bound platform gate unmet.
Grace independently verified the tarball, source-tree binding, inventory,
and retained local check logs, but did not rerun the check and makes no
platform-clean verdict. Full reviewer scope is recorded in
`provenance/independent-reviews-780c084-2026-10-16.md`.

**Current verdict: NOT READY. Highest proven rung: `tarball-clean`; next
unproven rung: `platform-clean`.** G8 remains open for checksum-bound Windows
results and a current exact-artifact platform result. G9 remains open for final
review of the current exact artifact and an unmerged release-evidence PR. No
CRAN submission is recorded or authorized.

## Chrome PR-state and independent site-review refresh, 2026-10-08

After pushing commit `e46d4ff`, Chrome shows release-evidence PR [#228](https://github.com/itchyshin/pigauto/pull/228) remains Draft with ten commits, three conflicts, no review, and no deployment. The Checks tab reports no checks for that commit. GitHub reports merge conflicts in `docs/dev-log/cran-0.11-audit/provenance/pre-pr-check.md`, `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`, and `docs/dev-log/cran-0.11-audit/surface-inventory.md`. The PR body still describes the earlier 5,128,732-byte candidate with SHA-256 `10573d4ff605ada47f5edae3ec018f959e745f7a30e5a4aa1172b1f696854ec1` and its passing macOS/R 4.6.0 check; the new 1bfed5ad receipt is in the branch files. The displayed Win-builder filename and size do not bind either Windows result to either hash. At this PR-state check, the release ledger was `NOT_READY` with empty artifact and evidence fields; it now records the separate 1bfed5ad candidate below and remains `NOT_READY`.

Chrome shows the separate source/documentation PR [#231](https://github.com/itchyshin/pigauto/pull/231) remains Draft with nine commits, three successful checks and one skipped check, no review, no conflict, and no deployment. Its local site and source checks do not establish deployed status. The independent site reviewer reran the retained local crawler and confirmed 62 HTML pages, 3,548 local references, 34 retired routes, and zero errors, including the planted Markdown-route negative control. This is local structural evidence only; visual inspection, post-correction deployment and live-route checks, the final artifact, and BirdTree redistribution evidence remain open. No merge or submission is authorized by these receipts.

## 2026-10-07 independent review and browser-gate reconciliation

Rechecked against the current worktree after independent review. `R/data.R` and `inst/NOTICE` both state that the `megatrees` package's MIT licence does not establish a licence for the underlying tree data; the earlier reviewer concern about those source files was stale. The rendered snapshot at `/private/tmp/pigauto-cran011-site-20261007` still serves stale `tree300` and `tree_full` reference pages. The current getting-started source and rendered article both call `tree300` an MCC phylogeny, which conflicts with the verified posterior-member-69 provenance. The live lease registry assigns `man/` to `codex-cran-011-implementation`. It does not list `vignettes/getting-started.Rmd`, although another ref has activity there; ownership remains to be checked before editing. This lane has left both untouched. Local visual review remains incomplete after the browser's explicit `file:` URL policy rejection; no workaround was attempted.

## Cache-busted deployed-page recheck, 2026-10-07

Check performed: Opened the current public getting-started article, tree-uncertainty article, `tree300` reference page, and `tree_full` reference page in the Codex in-app browser with a fresh query parameter. All four identify the deployed site as version 0.11.0. No local preview or deployment was used.

Observed: The live getting-started page still calls `tree300` a Hackett MCC phylogeny. The live `tree300` and `tree_full` references still say the posterior sample is distributed by the MIT-licensed `megatrees` package and source them as “Li (2026), megatrees 1.0.0, MIT licence”. The live tree-uncertainty page accurately limits itself to descriptive prediction sensitivity, but contains only Jetz et al. (2012) and does not explain BirdTree retrieval or cite the BirdTree tool. These conflict with current candidate source/help wording and the access workflow selected for 0.11.

IMPACT: G4 and G7 remain open. PR #230 remains draft and changes only `vignettes/tree-uncertainty.Rmd`; its current diff does not include the explicit user-run retrieval guidance present in this worktree. The public tree-uncertainty article also lacks those instructions and the `trees300` reference page still has stale provenance and licence wording. The active `man/` lease is confirmed; ownership of `vignettes/getting-started.Rmd` needs checking because the current lease registry does not list it. Rebuild the candidate site after the relevant source and generated-help updates, rerun route/search/sitemap checks, then verify the deployment after merge.

### Retrieval guidance is not yet on PR #230

The owned worktree diff for `vignettes/tree-uncertainty.Rmd` contains a “Supplying posterior trees” section with the explicit `megatrees::get_tree_bird_n100()` call, BirdTree subset-tool alternative, citation instructions, and reproducibility note. The exact one-commit diff on PR #230 adds the sample/backbone provenance and references, but not that retrieval section. Therefore PR #230, as currently opened, does not yet implement the selected user-managed retrieval recommendation. Its owner should carry the retrieval section (or an equivalent reviewed workflow) onto that branch before source review is complete. The file is under the active citation-lane lease, so this lane did not edit it.

## Independent tree-help review, 2026-10-07

The independent site auditor compared the current candidate source, generated help, and offline rendered snapshot. `R/data.R`, `inst/NOTICE`, and `man/trees300.Rd` correctly distinguish the megatrees package licence from tree-data rights and identify the sample. `man/tree300.Rd` and `man/tree_full.Rd` still contain stale MIT wording; the offline snapshot renders that stale help. The getting-started source and snapshot still call `tree300` an MCC tree. A fresh in-app-browser check separately confirms that deployed `trees300` help remains stale, and the deployed tree-uncertainty vignette has no retrieval instructions. The active `man/` lease covers these help files. Current preflight does not list `vignettes/getting-started.Rmd`, so check its ownership before editing.

## Defaults audit addendum, 2026-10-07

The current worktree HEAD includes `ea75608`, which adds a regression test for reconstruction of a synthetic legacy-shaped pre-transformer fit. The test omits the transformer-specific config fields and checks the reconstructed prediction against a known linear result. Unlazy gate I11 records the focused `fit-predict` command exiting 0 with `LEGACY_SHAPED_FIT_OK`. This covers the legacy configuration fallback, but the fixture is not an object serialized by an older package version and does not close the broader effective-route audit; G1 remains partial.

## Post-PR #230 deployment recheck, 2026-10-07

The repository's `pkgdown.yaml` run [#583](https://github.com/itchyshin/pigauto/actions/runs/37716522073) completed successfully on merged `main` commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b`. A cache-busted in-app-browser check of the deployed tree-uncertainty page confirms that PR #230's provenance paragraph and Ericson, Hackett, and Jetz references are live. The deployed article still has no user retrieval section, BirdTree subset-tool instructions, or access-tool citation. The current candidate worktree has a retrieval section, but it is not present on merged `main`; G4/G7 therefore remain open pending source reconciliation, a fresh site build, and a later deployment check. Other stale help/getting-started mismatches and sitemap verification also remain open.

## Ownership refresh, 2026-10-07

The live `lane_lease.sh --list pigauto` registry confirms the implementation lane `codex-cran-011-implementation` owns `README.md`, `vignettes/getting-started.Rmd`, `NEWS.md`, and `man/` as well as its check files. The tree-citation lane still owns `vignettes/tree-uncertainty.Rmd`; the defaults/root lane owns this audit ledger, its unlazy ledger, the defaults inventory, and the defaults test. An attempted narrow lease for getting-started did not register because the lease command reported a busy registry lock. No source edit was made to that file. Its stale MCC wording remains assigned to the implementation lane for reconciliation.

## Explicit tree-wrapper forwarding check, 2026-10-07

Added a focused regression in the root-owned defaults test for `multi_impute_trees(..., share_gnn = FALSE)`: explicit `lambda_mode`, `predict_method`, and `joint_solver` values reach each per-tree `impute()` call, and the default GNN-off value is forwarded. The first fixture used one tree and was correctly rejected by the public function's two-tree minimum before reaching the mock; the fixture now uses two trees with matching tip labels. Unlazy gate I1 passed the revised defaults, optional-backend fixture, and provenance test command with exit 0 and `CRAN_DEFAULTS_MI_AUDIT_OK`. This closes that specific forwarding gap; G1 remains partial pending reconciliation of the remaining effective-default routes and current-source integration evidence.

## `with_imputations()` default failure policy check (2026-10-08)

Added an effective-behavior regression to the root-owned defaults test. It calls `with_imputations()` without `.on_error`, deliberately fails both fits, captures the warning without discarding the returned object, and confirms both `pigauto_mi_error` elements remain in their slots with `n_fits = 2`, `n_failed = 2`, and `failed = 1:2`. The reviewed I1 command passed in 6.0 seconds: 284 passed, 0 failed, 1 expected warning from the separate two-tree sensitivity fixture, 0 skipped, and `CRAN_DEFAULTS_MI_AUDIT_OK`. An initial assertion harness incorrectly assigned the value returned by `expect_warning()`; it was corrected to capture the warning via `withCallingHandlers()`, then the exact suite passed. This confirms the default continue-and-retain contract; it does not close G1's remaining route-coverage gap.

Added a separate `with_imputations()` test for forwarding `...` to `.f` on every dataset. The same focused I1 command passed in 8.6 seconds: 286 passed, 0 failed, 1 expected low-draw tree diagnostic warning, 0 skipped, and `CRAN_DEFAULTS_MI_AUDIT_OK`.

**Superseding clean rerun:** The tree-wrapper fixture now captures and asserts its low-draw warning. The exact I1 command passed in 5.4 seconds: 287 passed, 0 failed, 0 warnings, 0 skipped, `CRAN_DEFAULTS_MI_AUDIT_OK`. This supersedes the prior one-warning receipt above.

**Superseding forwarding rerun:** Added runtime assertions that `multi_impute(..., draws_method = "auto")` forwards explicit `lambda_mode`, `predict_method`, and `joint_solver` values into the conformal `impute()` route. The reviewed I1 command passed in 5.9 seconds: 290 passed, 0 failed, 0 warnings, 0 skipped, `CRAN_DEFAULTS_MI_AUDIT_OK`.

## Cache-busted public-site recheck, 2026-10-08

The Codex In-app Browser loaded the current 0.11.0 homepage and, in a follow-up with cache-busting URLs, the tree-uncertainty article, `tree300` and `tree_full` references, getting-started article, and multiple-imputation article. The tree-uncertainty page now has mixed-backbone provenance and the Ericson, Hackett, and Jetz citations. It still lacks user-run retrieval instructions and the BirdTree web-tool citation guidance present in the candidate worktree. The deployed tree references still associate the data with the `megatrees` MIT licence; getting-started still calls `tree300` a Hackett MCC tree; and the inverse-Wishart paragraph lacks the candidate source's historical-reproduction boundary. The homepage presents the baseline-first workflow and optional GNN installation. Live search returns no results for `bench_continuous`; the former `/dev/bench_continuous.html` and `/VALIDATION_LEDGER.html` routes return 404. These checks cover two retired routes and one search query, not the full retired-route inventory. Direct sitemap navigation returned `net::ERR_BLOCKED_BY_CLIENT`, so sitemap contents remain unverified. Details and direct URLs are in `site-review-2026-10-08.md`. G4 and G7 remain open until source corrections are merged and deployed, then verified.


## Independent ledger and source review, 2026-10-08

The read-only site auditor confirmed that PR #230 merged as `0b0f71fee838c6ed51ef832ed819270eeafaf29b`, that candidate HEAD `ea756083e751e37f650f8801a9fc32e81c40a0c5` is in its ancestry, and that the merge changes only `vignettes/tree-uncertainty.Rmd`. Pages run #583 is recorded as successful on that merge. The auditor confirmed the candidate `R/data.R` and `inst/NOTICE` changes remain uncommitted, `man/tree300.Rd` and `man/tree_full.Rd` remain stale, and candidate `man/trees300.Rd` is corrected. Earlier statements that the getting-started file had no active lease are historical and superseded by the later ownership refresh.

### Cache-busted current deployed-page audit (2026-10-08)

Opened all five URLs with cache-busting query parameters in the Codex in-app browser. The cache-busted `/reference/trees300.html` page still says the 50-tree sample is from the Hackett backbone and still attributes the MIT licence to `megatrees`; this supersedes an earlier non-cache-busted view that showed mixed-backbone wording. `/reference/tree300.html` and `/reference/tree_full.html` also retain stale Hackett/MIT wording. The `/articles/tree-uncertainty.html` page correctly describes the mixed Ericson/Hackett sample and limits it to descriptive sensitivity, but lacks the candidate's user-controlled retrieval instructions and BirdTree web-tool citation guidance. `/articles/getting-started.html` still calls `tree300` a Hackett MCC tree. The separately cache-busted `/articles/multiple-imputation.html` still describes the inverse-Wishart option without the candidate's legacy-reproduction-only boundary. These observations directly verify current deployed text on the six named routes after Pages run #583; they do not verify other retained routes, all retired URLs, the sitemap, or candidate visual layout. G4/G7 remain open pending source repair, local rebuild, and a later deployment check.

Current status-only Unlazy verification parses 11 release gates (5 met, 6 unmet) and 13 implementation gates (12 met, I7 unmet); it does not run their `CHECK:` commands. `git diff --check` passes for these ledger edits.

## Optional-backend missing-package sentinel and Chrome restriction, 2026-10-08

The installed neither-backend check now serializes and reloads class-shaped drmTMB and gllvmTMB sentinels, then verifies `pool_mi()` gives the package-specific adapter error. The approved Unlazy I12 command passed with `SAVED_OPTIONAL_CLASS_MISSING_PACKAGE_OK` and `NEITHER_BACKEND_OK`; the retained raw output is `/private/tmp/pigauto-cran011-neither-saved-class-2026-10-08.log` (SHA-256 `334e8ca94f3b17898c1b414ca275907ced1e7b85330b272464073a5646634615`). These sentinels check the missing-package path and are not real fitted objects.

Chrome was used as requested to open the fresh local site, but browser security rejected the direct `file://` URL and prohibits local-server, alternate-browser, or indirect workarounds. G5 visual inspection remains unmet. The deployed live-page checks remain separate from local visual review.

## Current candidate source package check, 2026-10-08

The approved Unlazy I13 check completed in 9m45.1s with `0 errors`, `0 warnings`, and `0 notes` (`Status: OK`). It ran with `NOT_CRAN` unset and `force_suggests = TRUE` against HEAD `ea756083e751e37f650f8801a9fc32e81c40a0c5` plus the then-current source/docs diff. The implementation-file diff (R/data, data generators, notice, generated help, adapter harness, defaults test, and two vignettes) had SHA-256 `468f058c631e30549de268bd5b895f90bf10c1ee7b5c6630de0575ae49ccab88`. This establishes current candidate-source check status only; it does not validate a frozen tarball or platform results. The earlier 120-second attempt and the interrupted run did not count as results.

## Defaults and adapter review refresh, 2026-10-08

The defaults reviewer confirmed explicit `gnn = TRUE` forwarding assertions are present for `multi_impute()` on both conformal and `mc_dropout` paths. The reviewed Unlazy I1 command passed with 295 assertions, 0 failures, 0 warnings, and 0 skips in 6.2 seconds. Full output is `/private/tmp/pigauto-cran011-defaults-mi-20261008-after-review.log` (SHA-256 `d448fc313400b5b660298d410323b5e6336d0fc293e78f62a9a741585a89f3e7`). This closes that specific forwarding gap; broader caveats about argument interactions and scientific default optimality remain.

The MI reviewer confirmed I12 closes the recognized-class missing-package path using serialized class-shaped sentinels. They are reloaded in the same R process that serialized them, in a fresh `Rscript --vanilla` process where both optional packages are absent. This is not a real fitted object or a separate child-process reload. The genuine fit save/reload checks remain in the either/both-installed integration cells.

## Fresh official BirdTree page check, 2026-10-08

Chrome rechecked the official [downloads](https://birdtree.org/downloads/), [FAQ](https://birdtree.org/faq/), and [methods](https://birdtree.org/methods/) pages. Downloads says full and partial trees may be obtained and requires the Jetz et al. (2012) citation for research use; use of the BirdTree web tool also requires citing BirdTree. The FAQ discusses use of downloaded tree distributions, while the methods page describes the Hackett and Ericson backbone constraints and cites the source literature. A targeted search of BirdTree's indexed pages for redistribution or licence terms returned no results. These pages support the maintainer's statement that tree data may be used for research with attribution, but they do not document permission to redistribute the bundled tree objects in pigauto's CRAN tarball. No rights clearance is inferred; the release rights gate remains open.

## Current PR-state and independent tree-help review, 2026-10-08

The latest cache-busted Chrome view before this ledger update showed release-evidence PR #228 Draft with 13 commits, zero checks, no review, no deployment, and conflicts in `provenance/pre-pr-check.md`, `provenance/rights-and-policy.md`, and `surface-inventory.md`. Earlier browser snapshots showed seven and then 12 commits and are superseded by this refreshed page. The current PR description labels the 5,128,732-byte `10573d4f…` archive as a previous audit artifact and says final source/site changes require a fresh tarball. It does not mention the separate exact local candidate receipt for SHA-256 `1bfed5ad…` now in the branch. The 1bf candidate remains tied to merged source commit `0b0f71f…`, predates the open source/documentation changes in PR #231, and has no Windows result bound to its hash.

The independent read-only review of PR #231 at `c6235fa` found its single-tree/multi-tree routes, prediction-sensitivity limit, citation guidance, and distinction between the megatrees software licence and underlying BirdTree data rights consistent with source and official BirdTree pages. It found one unsupported claim still present in `inst/NOTICE`: tree member 69 is said to tie for the smallest Robinson–Foulds distance, but the branch's own acquisition ledger says no reproducible comparison was found. It also found that README and Getting Started do not direct BirdTree users to download instructions. Those are source-lane findings; this evidence lane did not edit PR #231. PR #231 remains Draft and undeployed, so neither its local build nor this review verifies the live site. G4, G7, G8, and G9 remain open; the redistribution-rights basis remains unresolved.


## PR #228 post-conflict Chrome verification, 2026-10-08

After pushing merge commit `49a52aa` to `release/cran-0.11-gate`, Chrome showed PR #228 Draft with 15 commits, no conflicts with `main`, all checks passed with one skipped check, no reviews, and no deployment. The page still disables merge because the PR remains Draft. The branch contains the updated evidence description, which now identifies the separate 1bfed5ad locally checked candidate and keeps the final artifact, exact-hash Windows/platform checks, BirdTree redistribution rights, final site verification, and exact-artifact review open. No merge, deployment, or CRAN submission occurred. This is a point-in-time PR check; it does not close G7, G8, or the rights gate.


The evidence-ledger follow-up commit `2f42c4e` updated the audit record only. A second cache refresh in Chrome showed PR #228 at 16 commits with the same clean merge state, passing checks with one skipped, Draft status, no reviews, and no deployment. The PR description still says the release is `NOT READY` and identifies the 1bfed5ad hash as a prior local candidate; no exact-final-artifact or rights gate was closed by this docs-only update.

## Fresh Chrome state check, 2026-10-08

Chrome reloaded both open PRs. PR [#228](https://github.com/itchyshin/pigauto/pull/228) is Draft with 17 commits, no base conflicts, all required checks passing with one skipped pkgdown pull-request check, no reviews, and no deployment. The PR description still labels the 1bfed5ad archive as a prior candidate from merged main, not the final artifact after pending source and site changes. This live count supersedes the earlier 18- and 19-commit snapshots in this ledger and the PR body’s older point-in-time text.

Source/docs PR [#231](https://github.com/itchyshin/pigauto/pull/231) is Draft with 9 commits, no base conflicts, three successful checks and one skipped check, no reviews, and no deployment. Its current description reports that single-tree and multi-tree import guidance, citations, and the distinction between the megatrees software licence and BirdTree data rights are addressed in its branch. A completed independent read-only review at commit `c6235fa` found the same, but also found that `inst/NOTICE` retains an unsupported claim that posterior member 69 ties for the smallest Robinson-Foulds distance. The branch evidence says no reproducible comparison was found. It also reports no byte comparison between the bundled trees and upstream BirdTree files. These source-lane issues remain open; this evidence lane did not edit that branch. The local site build does not prove deployment. G0, G4, G6, G7, G8, and G9 remain open.

POST-PUSH CHECK: Immediately after pushing evidence commit `4596157`, a fresh Chrome reload showed PR #228 with 20 commits, all required checks passing with one skipped check, no conflicts, no reviews, no deployment, and Draft status. This later view supersedes the 17-commit pre-push snapshot above. GitHub's page did not explain the difference between the counts; no cause is inferred. The page displayed `4596157` as the latest commit. This check does not close any release gate.

FOLLOW-UP AT HEAD `a43f644`: A later Chrome reload showed PR #228 with 21 commits, all required checks passing with one skipped check, no conflicts, no reviews, no deployment, and Draft status. This supersedes the 20-commit snapshot. The additional count is consistent with the next evidence commit, but no inference is made about the earlier 17-to-20 change.

## Current official BirdTree acquisition and citation check, 2026-10-08

The official [Downloads page](https://birdtree.org/downloads/) offers complete tree archives or user-selected subsets, and requires Jetz et al. (2012) for research use of full or partial data; use of the web tool also requires citing BirdTree.org. The official [subset tool](https://birdtree.org/subsets/) has a 2,500-species cap, takes a user-selected species list and tree distribution, and returns a downloaded tree sample with metadata and citations. The checked pages do not state a data redistribution licence or permission to bundle the tree bytes in a CRAN package. This does not negate Shinichi's statement that the trees may be used with attribution; it leaves the separate redistribution basis undocumented in these sources. Recommendation: pigauto should direct users to BirdTree's official download/subset process and accept a locally selected tree, rather than fetch a tree automatically. This preserves user choice and the requested citation workflow. It does not close the maintainer's rights gate for the current bundled datasets. Keep G0 open until the redistribution basis is documented or the distribution of those bundled files is changed under an explicit implementation decision.

## PR #231 candidate-source CI matrix recheck, 2026-10-08

Chrome verified [GitHub Actions run 37850907477](https://github.com/itchyshin/pigauto/actions/runs/37850907477) on PR #231 head `2f7305c` (14 commits at the time of the check). Ubuntu R release passed in 11m06s, Ubuntu R-devel passed in 10m47s, and macOS R release passed in 16m47s. The macOS log identifies R 4.6.1 on aarch64-apple-darwin23 and shows package tests and vignettes completed successfully. The PR-triggered pkgdown workflow was skipped by the repository's configured pull-request guard. This is candidate-source CI evidence only; it does not validate a frozen tarball, a deployment, or the final merged source. The exact-artifact, deployed-site, rights, and independent-review gates remain open.

## Current source-CI and deployment state, 2026-10-09

GitHub Actions run [#37945408182](https://github.com/itchyshin/pigauto/actions/runs/37945408182) completed successfully for source PR #231 at head `83117b48b12c864985fc44387b58e47585bb0ce0`. Chrome showed all three matrix jobs succeeded: Ubuntu R-release in 7m24s, Ubuntu R-devel in 12m30s, and macOS arm64 R-release in 20m06s; total workflow duration was 20m13s. The R CMD check jobs reported `Status: OK`. The workflow was triggered by a pull request and sets `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`; it is candidate-source CI, not force-Suggests evidence for a frozen post-merge tarball. The tested delta from previously passing head `a93fc80db76a0c94724031f9ea6be091dd3e1ca4` changes only `docs/dev-log/cran-0.11-audit/GATES.md`, its after-task evidence report, and a retained log for the preceding run #37940225845. Package code, bundled data, generated help, user-facing documentation, and website inputs are unchanged. The run #37945408182 result is visible in Chrome; a consolidated raw log for this run is not retained in this evidence branch.

PR #228 post-push check (2026-10-09): GitHub shows evidence-branch head `cf23da2`, Draft status, and no deployment. Its pull-request workflow run #37950876425 completed as skipped pkgdown in one second under the configured pull-request guard. No package check or deployment was produced by this evidence-only push. G7 through G9 remain open.

PR #231 description refresh (2026-10-09): the live rendered body now records completed run #37945408182 on `83117b48b12c864985fc44387b58e47585bb0ce0`, names the three successful candidate-source R CMD check jobs and their durations, and states the `NOT_CRAN=true` / force-Suggests-disabled limits and skipped PR pkgdown workflow. It explicitly does not claim exact-artifact validation. Browser verification is recorded in `after-task-2026-10-09-pr231-ci-description.md`. This wording correction does not change gate status.

Point-in-time state correction (2026-10-09, before the subsequent source-CI refresh): PR #231 was Draft and unmerged at `dee8f02888fd2c4c22d426e595ea6451a28a280e`; PR #228 was also Draft and unmerged. Earlier main merges did occur: PR #229 was deployed by pkgdown run [#581](https://github.com/itchyshin/pigauto/actions/runs/37706253267), and PR #230's merge commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b` was deployed by pkgdown run [#583](https://github.com/itchyshin/pigauto/actions/runs/37716522073). Those deployments did not include the then-current PR #231 source/site changes. The deployed site therefore remained an earlier baseline, and the complete G7 check against PR #231's eventual deployed commit remained open. The ledger was 8 of 11: G7 deployed-site verification, G8 one final frozen artifact with exact checks, and G9 independent final-artifact review were unmet. No final post-merge tarball or CRAN submission existed. The evidence branch recorded only the explicitly predecessor tarball; it was not promoted to final. PR #231 later advanced to `3702d5d533e101dd3583b30703200fe525d1db4f`, as recorded in the exact-head CI receipt below.


## Independent current-head BirdTree source review, 2026-10-09

The source reviewer rechecked PR #231 at exact head `dee8f02888fd2c4c22d426e595ea6451a28a280e`. The earlier unsupported closest/best Robinson-Foulds claim is absent: `inst/NOTICE:41-43` expressly makes no closest, best, or recommended-tree claim, while `R/data.R:33-39, 89-96` and Getting Started describe member 69 as an example. User acquisition guidance is present in `README.md:42-48`, `R/read_tree.R:7-18, 27-40`, and `vignettes/tree-uncertainty.Rmd:39-58`: users select/download BirdTree data, Newick is distinguished from NEXUS support for other sources, citation is specified, and pigauto does not download automatically. Generated help was checked against source. The focused source tests passed 52 testthat expectations across two test contexts, with 0 failures, 5 warnings, and 0 skips, in 8.1 seconds. Raw output: `provenance/readtree-source-review-dee8f028-2026-10-09.log`, SHA-256 `0ec6dc71c01533461dc43e4895e7d81f34d030a23d32cffb3f6f750cbceaa356`. This closes those source-review findings for the PR #231 candidate, not deployed-page verification at that earlier snapshot: PR #231 was then unmerged, so G7 remained open; G8 and G9 also remained open pending the exact final artifact and its independent review. The historical 2026-10-08 report remains a dated record of the earlier finding, not current guidance.

## PR #231 current-main source CI, 2026-10-09

Chrome verified GitHub Actions run [#37952599321](https://github.com/itchyshin/pigauto/actions/runs/37952599321) for PR #231 head `dee8f02888fd2c4c22d426e595ea6451a28a280e`, tested with synthetic merge `10c3acbb49c1c79f2d71f58fc41caccf90536a4c` onto current main `0b0f71fee838c6ed51ef832ed819270eeafaf29b`. All three jobs passed: Ubuntu R-release 8m09s, Ubuntu R-devel 13m16s, and macOS arm64 R-release 23m10s; the workflow completed in 23m28s. On macOS, focused MPS tests took 6m04s and R CMD check took 12m10s. This pull-request run uses `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`, so it does not satisfy the clean-library or exact-artifact gate. GitHub's live run and job pages are the retained evidence; a consolidated raw log could not be downloaded in this environment. The delta from source head `83117b48b12c864985fc44387b58e47585bb0ce0` to `dee8f0` changes only `GATES.md`, the after-task report, and `provenance/rights-and-policy.md`; package code, bundled data, generated help, user-facing documentation, and website inputs are unchanged. This source-CI result does not close G7, G8, or G9.


## PR #228 BirdTree-review receipt check, 2026-10-09

Chrome verified the push of evidence commit `5ad0bf9` to PR #228. GitHub Actions run [#37960057428](https://github.com/itchyshin/pigauto/actions/runs/37960057428) was triggered by the pull request and skipped in 8 seconds under the configured pkgdown pull-request guard; it produced no job result or artifact and did not deploy the site. PR #228 remains Draft and unmerged. This check confirms the evidence commit is present on the PR branch only; it does not change the 8-of-11 gate tally or close G7, G8, or G9.

## Fresh live article-index and retired-route check, 2026-10-09

Chrome opened the deployed article index at `https://itchyshin.github.io/pigauto/articles/index.html?audit=20261009lead2`. It still lists “Four ways to impute a phylogenetic trait matrix” linking to `/articles/simulation-study.html`. Opening `https://itchyshin.github.io/pigauto/articles/simulation-study.html?audit=20261009lead` served the full historical article, including its comparison tables, results, reproduction section, and notice that it is excluded from the intended 0.11.0 package build and public site. This confirms the route is both discoverable and directly served on the current deployment.

Chrome returned `net::ERR_BLOCKED_BY_CLIENT` for the live `/sitemap.xml` and `/search.json` endpoints; their contents remain unverified. The local candidate crawler's retired-route result remains local build evidence only. At this check, PR #231 was Draft and unmerged at `dee8f02888fd2c4c22d426e595ea6451a28a280e`; it later advanced to `3702d5d533e101dd3583b30703200fe525d1db4f`. The live mismatch was expected to remain until merge and deployment. G7 remained open alongside G8 (final frozen artifact and exact checks) and G9 (independent review). The tally remained 8 of 11 gates met.


## Exact-head source CI receipt, 2026-10-09

GitHub Actions run [#37968809743](https://github.com/itchyshin/pigauto/actions/runs/37968809743) completed successfully in 18m02s for PR #231 at source commit `3702d5d533e101dd3583b30703200fe525d1db4f`. The Ubuntu R-release job `113949985095`, Ubuntu R-devel job `113949985755`, and macOS arm64 R-release job `113949985847` each passed. Each full R CMD check reported `Status: OK` and testthat `[ FAIL 0 | WARN 175 | SKIP 83 | PASS 3196 ]`. The macOS focused MPS test reported `[ FAIL 0 | WARN 50 | SKIP 0 | PASS 201 ]`. The workflow used `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`; it is candidate-source evidence and does not satisfy the force-Suggests exact-artifact gate. Companion pkgdown run [#37968809770](https://github.com/itchyshin/pigauto/actions/runs/37968809770) was skipped by design and did not deploy the site.

The source diff from `cff6ea685d2788687d085e24c2a8ef83ce50f7bb` to `3702d5d533e101dd3583b30703200fe525d1db4f` changes only `GATES.md` and `after-task-2026-10-09-source-log-compression.md`; package behavior, defaults, data, generated help, user-facing documentation, and website inputs are unchanged. The combined raw job logs are retained as [`provenance/source-ci-37968809743.log.gz`](provenance/source-ci-37968809743.log.gz). The combined uncompressed stream is 3,073,186 bytes with SHA-256 `ba14e68809ff5d206ca37f5bddb3aea0905fd44fc3f451183d4da40db04aab6f`; the deterministic gzip is 409,760 bytes with SHA-256 `291a3c2b3cadd07b2e854a528a54934a03ac14080c1a7808f990e4053d78a4e0`. Decompression matched the fetched job bytes exactly. Individual raw job sizes and SHA-256 values are recorded in `after-task-2026-10-09-source-ci-379688.md`.

CURRENT STATUS (2026-10-09, post-merge): PR #231 merged as `d76804e768bf77f430f43dcf591633f9cd900dab`; the Pages deployment and corrected checks of 44 retired URLs, all 62 sitemap URLs, the sitemap, and search index pass. PR #228 remains Draft and unmerged. The tally is 9 of 11 gates met: G8 remains open pending exact-hash Windows results and remaining artifact-specific release checks; G9 remains open pending independent exact-artifact review. No CRAN submission has occurred.


## Current exact-artifact and site-evidence reconciliation, 2026-10-09

The exact archive remains `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4` from `d76804e768bf77f430f43dcf591633f9cd900dab`. macOS `R CMD check --as-cran` is `Status: OK`; the optional-backend exact-tarball matrix passes. Current Win-builder R-release and R-devel logs each report `Status: 1 ERROR` and five failed expectations, but neither log supplies an archive checksum. Their failure details and raw/compressed hashes are retained in `provenance/winbuilder-results-2026-10-09.md` and the deterministically compressed logs. Thus G8 remains open and no Windows result is attributed to the frozen hash.

The corrected deployed-site verification artifacts were present locally but ignored by the repository rule for `docs/`; they are being explicitly included in this evidence branch. They bind the 44 retired-route checks, 62 sitemap URLs, and 613-entry search-index check to Pages run `37975894397` at `d76804e768bf77f430f43dcf591633f9cd900dab`. This makes the site evidence available for review, but does not change Pat's verdict until the durable receipts are independently assessed. The live homepage's literal `[!WARNING]` rendering also remains a reader-facing discrepancy under its separate source fix.

Independent closure: Grace and Rose return NOT READY for the exact artifact because of the Windows diagnostics and missing checksum binding. Pat returns NOT ASSESSED for durable G7 evidence because the v2 receipts were not tracked at review time. G9 remains open; overall verdict remains **NOT READY**. No merge or CRAN submission has occurred.


## Pat independent site re-review, 2026-10-09

Pat has now independently reviewed the committed v2 receipts and passes G7 for the checked deployment (`d76804e768bf77f430f43dcf591633f9cd900dab`, Pages run `37975894397`). The earlier NOT ASSESSED status is superseded. The 44 retired-route checks, 62 sitemap URLs, 613 search entries, and three verifier tests are documented in `provenance/pat-independent-site-review-2026-10-09-v2.md`. Pat still observed the literal `[!WARNING]` token on the homepage; direct live sitemap access was blocked during the re-review. This pass is limited to the checked deployment. G8 and G9 remain open: Windows logs fail and are not bound to the frozen archive hash, and Grace/Rose remain NOT READY. Overall NOT READY; PR #228 remains Draft and unmerged.


## Evidence PR refresh, 2026-10-09

The Pat G7 re-review receipt and refreshed NOT READY description are committed to PR #228 at `b04efdcb09eed9149196ace5e137c6b4a954508e`. Chrome confirms the current head, the revised description, Draft state, and passing PR checks (one skipped). G5 is met. G3 and G4 remain open; PR #228 is not merged.

## Candidate homepage warning build and cleanup check, 2026-10-09

The separate candidate source commit `45df106ba58f5bc97d7a274bc0a22a19d202fa17` (parent `d76804e768bf77f430f43dcf591633f9cd900dab`) changes the README warning syntax so pkgdown renders a “Warning:” callout instead of exposing GitHub's `[!WARNING]` token. A direct `pkgdown::build_site()` output was then passed through the same `pkgdown/clean-internal-pages.R` production cleanup. The build log is 46,943 bytes, SHA-256 `b776a7a5a88dd6c61b07295dc5ffb5d2cca6a5960758f2639c1ab268b4a3188d`.

Added `script/cran-0.11-site/verify_candidate_site.py` and five initial offline tests. Pat's independent reviews found false-pass cases in warning structure/body validation, hidden classes and stylesheets, inert script/style text, incomplete discovery indexes, external URLs, substituted or truncated manifests, percent-encoded paths, and traversal safety. The verifier now pins the reviewed 44-route receipt by SHA-256, checks exact sitemap/search counts and full HTML coverage, enforces the configured origin and `/pigauto/` path prefix, compares decoded routes, checks candidate-root containment, and scans the warning's HTML ancestry/content against hidden attributes, common hiding classes/declarations, attribute selectors, same-site imported CSS, and linked stylesheets. Fifteen tests pass, covering 30 negative-control scenarios; the test receipt is `provenance/candidate-site-verifier-tests-2026-10-09.log`. The strengthened verifier passes against the same production-cleaned candidate: 44 manifest routes absent as output files and discovery targets, 62 sitemap entries matching all 62 HTML pages, and 613 search entries with 568 in-sitemap paths. The search-path count is enforced, while path multiplicities are not pinned. Captured output remains in `provenance/candidate-site-warning-callout-2026-10-09.log`, SHA-256 `444b18fca2d969928bea89643c803009d63379226d3de97e3240689f6dbf75f4`.

This static check covers the generated HTML and supported stylesheet selector/declaration cases, including inline styles and same-site `@import`; it does not run a browser's complete CSS cascade, test accessibility, or substitute for visual review. It is bound to the current 44-route receipt and exact page/search counts, which must be deliberately updated when the candidate source or manifest changes. This local candidate evidence does not change the deployed-site G7 verdict.

## Candidate warning verifier follow-up, 2026-10-09

After Pat found two further false-pass cases in the bounded static verifier, the check now catches decimal-zero opacity values such as `opacity: 0.0` and expands comma-separated selector alternatives inside `:is()` and `:where()`. The suite passes all 16 tests, including 35 negative-control scenarios. Pat's final read-only review passes this documented bounded claim and confirms the local candidate counts. A general computed-visibility claim remains unsupported: expressions such as `opacity: calc(0)` or a custom property that resolves to zero require browser CSS evaluation. The checker does not implement the full computed-style cascade or provide a visual/accessibility review. The candidate CLI still reports 44 retired routes absent, 62 HTML pages matching 62 sitemap entries, 613 search entries with 568 paths in the sitemap, and `CANDIDATE_SITE_OK`. This local candidate evidence does not change deployed-site G7, exact artifact G8, or independent release closure G9.

This is candidate-source output only. The complete GitHub Pages wrapper attempt was interrupted while rendering the mixed-types vignette after exceeding its estimate; the direct `build_site()` plus production cleanup completed. A browser policy rejected local `file:` preview, so no visual review of this unmerged candidate is claimed. The live deployed site has not changed, and its homepage still shows the old marker until the source fix is merged and deployed. This evidence does not change the deployed-site G7 verdict or close G8/G9.

## Candidate verifier follow-up push, 2026-10-09

Commit `74d5a9a` adds regression tests for decimal-zero opacity and comma-separated `:is()` / `:where()` selectors, with Pat's bounded static-check verdict and the refreshed local test receipt. The commit is on `release/cran-0.11-gate`. Chrome shows PR #228 at this head, still Draft and unmerged; all checks pass with one skipped, and the page lists no deployment. This push updates evidence only. It does not deploy the candidate or close the exact-artifact and independent-review gates.

## Current evidence PR state, 2026-10-09

Chrome confirms PR #228 is Draft and unmerged, has 61 commits, and includes latest commit `325e3cd5a95325413d0adeca7b1707fc41714b3a`. The last commit records the live homepage warning recheck and its report validation. PR #229 is already merged; its conversation confirms the defaults and optional MI integration changes, along with the internal-page retirement, are in `main`. PR #231's merged source and deployed-site evidence remain as recorded above. The local acceptance ledger now matches the current #228 head. This browser check does not validate the exact tarball, produce a Windows result, close the reviewer gate, or authorize merge/deployment/submission. G8 and G9 remain open; no CRAN submission has occurred.

## Candidate-site verifier push check, 2026-10-09

Commit `1872bd7` was pushed to `release/cran-0.11-gate`. Chrome then showed PR #228 at 50 commits, still Draft, with the new commit visible. Pull-request workflow run [#38001692992](https://github.com/itchyshin/pigauto/actions/runs/38001692992) was skipped in one second under the configured pkgdown pull-request guard. The PR page showed no deployments. This verifies the evidence push only; it does not run pkgdown in GitHub Actions, deploy the candidate, or change G8/G9.

## Candidate verifier and portability-smoke evidence push, 2026-10-09

Commit `8fa149a` records Pat's bounded static-verifier verdict and the local macOS smoke of separate-lane commit `45bf230`; the smoke is explicitly not a Windows or frozen-artifact result. Chrome showed PR #228 at this head with 52 commits, still Draft and unmerged. Workflow run [#38003919177](https://github.com/itchyshin/pigauto/actions/runs/38003919177) was skipped under the configured pull-request pkgdown guard; no deployment was listed. This evidence-only push does not deploy the candidate or change G8/G9. The exact artifact remains NOT READY because both Windows logs report five failed expectations and one error without a checksum binding. No merge or CRAN submission has occurred.

## Portability-smoke evidence commit, 2026-10-09

Commit `41ddecc` adds the raw smoke log and its after-task report to PR #228. Chrome confirmed the PR remained Draft with 53 commits and showed workflow run [#38004504745](https://github.com/itchyshin/pigauto/actions/runs/38004504745) as skipped in one second under the pull-request pkgdown guard; no artifact or deployment was produced. This receipt confirms that the evidence is present on the PR branch only. G8 and G9 remain open, and no merge or CRAN submission has occurred.


## Exact-artifact optional MI backend matrix, 2026-10-09

The four isolated installed-library cells passed against archive SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`. Both backends, each backend separately, and neither backend were checked. The receipt, summary logs, run-level manifest (`fb9af9cbd7337358497ac8c72ecfcb8256f9504d01f8f696fac8ccdd39fb062e`), and four original cell logs are retained under `provenance/`; the manifest binds each log and before/after archive hash and records the installed test-script SHA. Grace and Rose independently passed this bounded adapter subcheck. This closes the optional-adapter integration subcheck for the frozen artifact only. G8 remains open because the Windows results are failed and not checksum-bound; G9 remains open pending independent exact-artifact verdicts. Overall verdict remains NOT READY.

## Win-builder failure mapping to the separate portability patch, 2026-10-09

Both retained R-release and R-devel logs report the same five test failures: mocked CUDA and MPS assertions in `test-check-pigauto.R`, a Lantern error in `test-covariate-alignment.R`, and two Lantern errors in `test-monomorphic-discrete.R`. Each location corresponds to a change in separate-lane source commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4`: the accelerator tests mock the torch-runtime probe, and the tensor tests now use the existing backend-ready skip helper. This is a direct source-level mapping only. The Mac-only candidate tests do not demonstrate a Windows pass, and the original Win-builder logs do not carry an archive checksum. The detailed review is `after-task-2026-10-09-winbuilder-failure-mapping.md`. G8 remains open; the frozen artifact is NOT READY.

## Executable ledger and evidence-PR refresh, 2026-10-10

The release validator initially rejected the ledger because `status_claim` was descriptive prose and the exact artifact lacked the validator's canonical path, inventory, and typed evidence fields. The ledger now records the highest verified rung as `tarball-clean`, the next unproven rung as `platform-clean`, and the overall audit verdict as NOT READY. Its source checkout status, archive identity, inventory, local check, installed tests, policy, rights, and rendered-site evidence are linked as files. The validator self-test passes all negative controls and the ledger returns `READY FOR CLAIMED RUNG`; this confirms only the local tarball rung and does not close G8/G9.

Chrome shows PR #228 remains Draft and unmerged at 66 commits, with latest head `427ef1f`. Workflow run [#38009821247](https://github.com/itchyshin/pigauto/actions/runs/38009821247) was skipped by the pull-request pkgdown guard; the page lists no artifacts or deployment. This verifies the evidence push only. The exact frozen artifact remains NOT READY overall.

## PR #233 portability source-CI and rights-source recheck, 2026-10-10

The separate portability source PR #233 is at commit `cd185ff3bb0c2b3fc157c7dfc6f8addec6b2d7a0`. Chrome currently shows it open, ready to merge, with three successful checks and one skipped; no deployment is listed. GitHub Actions run [#38067484909](https://github.com/itchyshin/pigauto/actions/runs/38067484909) reports three successful source-CI jobs (Ubuntu release, Ubuntu devel, and macOS arm64); the pkgdown pull-request job is skipped under the configured guard. The macOS check-r-package step succeeded. These are candidate-source results, not checks on a frozen tarball, and do not bind Windows results to the artifact hash. G8 and G9 remain open; the release remains NOT READY.

The dated rights-source supplement is `provenance/rights-recheck-2026-10-10.md`. It records the official BirdTree and CRAN policy pages and current megatrees package metadata. The megatrees MIT metadata does not by itself establish the redistribution terms for BirdTree source data. Shinichi's previously recorded maintainer warranty remains the disposition; no independent legal determination is claimed. Verify the shipped objects and notices against the eventual frozen artifact before closing the relevant check.

Chrome confirms evidence PR #228 remains Draft and unmerged at 68 commits, latest head `2921f94d1227a08ebcdbb960c386db338e914c05`. Its current workflow run [#38070572404](https://github.com/itchyshin/pigauto/actions/runs/38070572404) was skipped under the configured pull-request pkgdown guard, with no artifacts or deployment. The Unlazy recheck passed G1, G2, G6, and G7 after correcting G1/G2 working-directory paths in `.unlazy/cran-0.11-audit/GATES.md`. G3 (checksum-bound Windows/platform checks) and G4 (independent exact-artifact review) remain unmet. No release rung changes.

## Current deployed-site recheck, 2026-10-10

PR #232 merged the homepage warning correction as `4f60a18ae72e836b6c3e15b65abb1faa7409e731`; Pages run [#38063988962](https://github.com/itchyshin/pigauto/actions/runs/38063988962) succeeded from that commit. Pat independently reviewed the deployed site and reports that the exact-commit local build/crawl had 62 HTML pages, 3,584 local references, and 34 retired pages absent from output, sitemap, and search. In Chrome, all 62 paths from the local sitemap rendered their expected page titles/headings; the 34 retired routes and `/VALIDATION_LEDGER.html` showed the pkgdown not-found page. The live homepage shows the rendered 0.11.0 warning callout, Getting Started describes GPU use as opt-in, and a live search for the retired ledger returned no result. The browser did not expose numeric HTTP status codes, and direct live sitemap/search-JSON navigation was blocked (`ERR_BLOCKED_BY_CLIENT`), so this is a rendered-page and search-UI check, not a raw live-index or status-code audit. This closes G7 for the checked deployment at `4f60a18` only; it makes no permanence claim. Review provenance: Pat's independent review returned to this task on 2026-10-10; the durable receipt is this entry and the after-task addendum below.

The frozen tarball still records source commit `d76804e768bf77f430f43dcf591633f9cd900dab`, before PR #232. It is therefore stale for the current merged source and cannot be used to close the final artifact gate. PR #233 remains open; its passing source checks do not substitute for merging the portability fix, rebuilding from clean merged source, and binding platform checks and independent artifact reviews to that new archive. The release remains NOT READY.

## Fresh Windows-failure mapping and PR-state check, 2026-10-10

Chrome rechecked PR #233 at `cd185ff3bb0c2b3fc157c7dfc6f8addec6b2d7a0`. It remains open. Actions run [#38067484909](https://github.com/itchyshin/pigauto/actions/runs/38067484909) completed successfully in 29m48s with Ubuntu R-release, Ubuntu R-devel, and macOS arm64 jobs; pkgdown run [#38067484815](https://github.com/itchyshin/pigauto/actions/runs/38067484815) was skipped under the pull-request guard. These are source checks only.

The two retained Win-builder logs each show the same five failures: the CUDA and MPS runtime-probe tests in `test-check-pigauto.R`, one Lantern runtime error in `test-covariate-alignment.R`, and two Lantern runtime errors in `test-monomorphic-discrete.R`. The PR diff changes those exact tests and the private mock seam. This confirms the patch addresses the observed failure locations at source level; it does not establish that either Windows check now passes, and the old logs remain unbound to the frozen archive hash. Do not upload the predecessor archive. PR #228 remains Draft and unmerged at 69 commits, head `4a87a53b83d617054978a94363086bd2e122683f`; its workflow run [#38071388358](https://github.com/itchyshin/pigauto/actions/runs/38071388358) was skipped by the same guard. G8 and G9 remain open and the release remains NOT READY.

## Independent review follow-up on the Windows skip concern, 2026-10-10

Grace rechecked the specific concern about the local `skip_if_no_torch()` helper in `test-gnn-train-cal-symmetry.R` against both archived Win-builder logs. In each log, decompressed lines 150-152 list the four symmetry tests as skipped with `libtorch not installed`; they are not among the five failures. The five failures are at lines 232, 238, 244, 264, and 276, in the two mocked accelerator probes and the three direct tensor tests already changed by PR #233. Grace's bounded verdict is PASS for that five-failure scope. The concern remains hypothetical for a future runtime where `torch_is_installed()` is true but tensor setup fails. No source changes are justified by these logs. This does not substitute for Windows checks on a newly frozen exact artifact; PR #233 remains open and G8/G9 remain unmet.

See `after-task-2026-10-10-symmetry-skip-review.md` for the full review receipt and limits.

## Primary CRAN publication-history recheck, 2026-10-10

Chrome checked the current primary [CRAN pigauto package page](https://cran.r-project.org/web/packages/pigauto/index.html). It lists version 0.10.0, published 2026-07-30. The working `DESCRIPTION` and `cran-comments.md` select 0.11.0 as the next version, consistent with the release plan. Keep 0.11.0; there is no current reason to increment to 0.11.1. The archived CRAN-index receipt from 2026-10-06 also lists 0.10.0. This verifies the current published version only; artifact, Windows, and exact-artifact review gates remain open.

See `after-task-2026-10-10-version-history.md` for the verification record.

## Exact-artifact optional MI matrix with complete source provenance, 2026-10-16

The exact installed integration matrix was rerun against `pigauto_0.11.0.tar.gz`, SHA-256 `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`. All 12 install/fixture/real steps passed across `neither`, `drm-only`, `gllvm-only`, and `both`; the archive hash matched before and after each step. The retained manifest now records the explicit working directory, hashes all five R sources loaded by the harness, records expected optional-package availability and each cell's dependency symlink targets, and binds every log by SHA-256. The manifest and logs are in `provenance/mi-matrix-0fccf4bb-provenance-complete/`.

Available backends produced three real fits per model cell, passed fit-status checks, a drmTMB Gaussian regression oracle, independent termwise Rubin arithmetic for both adapters, and saved-fit reload/pooling in a fresh R process. The neither-backend cell verified installation and fixture independence and correctly failed closed for saved optional-backend classes. This supports optional adapter interoperability for this bounded analysis-aware MI fixture. Dependencies were linked from an existing package library, so this does not establish clean-room dependency installation. It does not establish broad recovery or interval coverage, joint covariance/loadings/variance-component pooling, Windows/platform cleanliness, or an overall release verdict. The prior exact-matrix review found the result passed but identified the now-repaired manifest provenance gap. A fresh independent review rehashed the archive, runner, harness, all five executed sources, and all 12 logs; confirmed all four recorded dependency states and successful steps; and returned PASS with no material discrepancy. Shared dependency symlinks mean this is not clean-room dependency installation. G8 and G9 remain open; highest proven rung remains `tarball-clean`, next unproven rung `platform-clean`.

## Current-ledger correction and Win-builder recheck, 2026-10-16

This addendum supersedes stale exact-artifact and reader-correction fields in the JSON ledger and earlier historical status paragraphs. The current frozen artifact remains `pigauto_0.11.0.tar.gz`, SHA-256 `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`, from source commit `780c08470b263cc26b4139c711a5c144fd29389c`. The optional-backend matrix and independent matrix review both bind to this hash; the older `4b32f669...` matrix entry was stale and is replaced in `release-ledger.json`.

The Getting Started correction PR #234 merged as `8082ef8273177ff664d49ac4bd643f7def671f20`, which is an ancestor of the frozen source commit. The corrected page is included in the Pages deployment at `780c084` (run `38970604986`). Current bounded panel receipts are in `provenance/independent-reviews-780c084-2026-10-16.md`: Grace verified current archive and local-check identity; Pat passed the supplied current-site route/discovery receipt and reviewed the NEWS historical-text exception; Rose found no current checksum-bound Windows result. The reviewers do not support a platform-clean or submission-ready claim.

A fresh Gmail search at 2026-10-16 23:49 UTC found no Win-builder result email. The upload response still attests only filename and byte size, with no server-side SHA-256. The latest exact macOS/R 4.6.0 `R CMD check --as-cran --run-donttest --timings` log reports `checking CRAN incoming feasibility ... OK` and `Status: OK`; an earlier sandboxed attempt had stopped at DNS before package checks. The six shipped R data-object streams and `inst/NOTICE` were separately SHA-256 hashed from the current archive and match the previously inspected bytes; this confirms byte identity, not third-party redistribution terms. The archive is described in `provenance/exact-artifact-780c084.json`. These checks do not fill the Windows platform gap. The public R-hub tarball route has not been used because a separate publication approval is required. PR #236 was refreshed with the evidence correction and was observed at 23:48 UTC as Draft/unmerged at `a4f64c0240e39923a31bc9b7110fa2a43fcbae63`; a later evidence-only commit may advance its head. G8/G9 remain open; overall `NOT_READY`, highest rung `tarball-clean`, next rung `platform-clean`.


## Win-builder result follow-up, 2026-10-17

Both current candidate result logs have now been retrieved from the official service: R-release (R 4.6.1) and R-devel (r90686), each on Windows Server 2022 x64 and each with `Status: OK`. Review by Grace and Rose independently confirms the statuses and the bounded attribution chain to the 5,136,083-byte upload of local SHA-256 `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`. Since the service gives no server digest and the logs expose no archive size, this closes the two Windows check outcomes but not strict exact-hash platform binding. G8 remains unmet; highest proven rung `tarball-clean`, next unproven rung `platform-clean`. See `provenance/winbuilder-results-fe1600d7-2026-10-17.md`.
~~~~

## `docs/dev-log/cran-0.11-audit/after-task-2026-10-17-current-artifact-closeout.md`

~~~~text
## 1. Goal

Bind current release evidence to one frozen pigauto 0.11.0 archive, refresh defaults, optional-model MI, website, and Windows checks, and keep the draft release PR unmerged.

## 2. Implemented

Updated the `simulate_benchmark()` help example to use `gnn = FALSE` so it works when optional torch is absent; the public function default remains `TRUE`. Added a regression assertion requiring nonempty benchmark output. Refreshed the release ledger and gate record for the frozen archive, current defaults audit, exact-hash MI matrix, R-hub Windows checks, and fresh local site review. The candidate remains `NOT_READY` because the deployed site and final release-validator closeout are open.

## 3a. Decisions and Rejected Alternatives

Kept the public `gnn = TRUE` default unchanged; the failure was in a documentation example's runtime assumptions, and the baseline-only example demonstrates that path without optional torch. Rebuilt and rechecked a new frozen archive after the example/test change rather than reusing older hash-bound results. Kept the live deployment gate open because source commit `9a457eb` has not been merged or deployed. Did not merge PR #236 or submit to CRAN.

Assumptions made without asking: treated the user's public-upload authorization as covering the R-hub public package upload; it does not authorize a merge, deployment, or CRAN submission. Kept the current package version at 0.11.0 based on the prior publication-history check recorded in the ledger.

## 4. Files Touched

- `R/simulate_benchmark.R`
- `man/simulate_benchmark.Rd`
- `tests/testthat/test-new-features.R`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/defaults-current-2026-10-17-9a457eb.log`
- `docs/dev-log/cran-0.11-audit/provenance/site-9a457eb-review.md`
- `docs/dev-log/cran-0.11-audit/provenance/windows-rhub-c2e768ad-2026-10-17.md`
- `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/manifest.json` and its 12 step logs
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-17-current-artifact-closeout.md`
- Frozen artifact receipts and logs: `exact-artifact-c2e768ad-source-binding.md`, `exact-artifact-c2e768ad-local-check.md`, `exact-artifact-c2e768ad-as-cran-check.log`, and `exact-artifact-c2e768ad-as-cran-testthat.Rout`

## 5. Checks Run

- `devtools::test(filter = "new-features")`: 112 expectations passed, 0 failures.
- One-seed `simulate_benchmark(..., gnn = FALSE)` smoke: completed and returned nonempty results.
- Current-source `devtools::test(filter = "cran-audit-defaults")`: 202 passed, 0 failed, 0 warnings, 0 skips; retained log hash `4c4d20917e712f09e93fab8bb62707edfce28c71ddae70af5c71df5a85d1cc99`.
- Exact frozen archive `R CMD check --as-cran --run-donttest --timings`, `NOT_CRAN` unset, force-Suggests enabled: `Status: OK`; testthat 3,195 pass, 0 fail, 161 warnings, 86 skips.
- Exact-hash optional-backend MI matrix: 12/12 steps passed across neither, drmTMB-only, gllvmTMB-only, and both-installed cells; hash unchanged before and after each step.
- R-hub Windows R-release 4.6.1 and R-devel r90686: both `Status: OK`, each 2,442 pass and 0 failures after R-hub rebuilt the source extracted from the byte-identical uploaded object.
- Fresh local pkgdown build and crawler: 62 HTML pages, 614 search entries, 3,588 references, 34 retired routes, zero crawl errors; `pkgdown::check_pkgdown()` clean.
- Browser content/accessibility review of the rendered `simulate_benchmark()` page completed. The independent site review confirmed retired route output and discovery checks.
- Unlazy `--reverify` approved and reran four read-only checks: current defaults, exact-hash MI manifest, current site build, and current crawl all passed. The repository ledger has 7 met gates and 4 open gates: local pixel-level review, deployed-site verification, post-merge artifact freeze, and final independent closure.
- `cran_release_gate.py release-ledger.json`: `READY FOR CLAIMED RUNG` at `platform-clean`; `--selftest` passed all 14 planted controls. The whole release verdict remains `NOT_READY`.
- `check-after-task.R` on this report: structure passed. `closeout.py check` still exits on the brain repository evidence gate because it sees nine unrelated unfinished brain ledgers; the repo-scoped structure validator passed.
- `slop_check.py` on the after-task, site, and Windows receipts: 0 findings. The checker emitted path warnings for ignored files under the brain source root, but its findings count was zero.

## 6. Tests of the Tests

The `simulate_benchmark()` regression now asserts `nrow(result$results) > 0`; it would fail if the example returned an empty result table as the previous Windows run did. The exact-artifact MI runner checks all four dependency configurations, all three steps per cell, and verifies the archive SHA before and after each step. The R-hub object readback was byte-compared with the frozen archive before the platform checks.

## 7a. Issue Ledger

- Fixed: the `\donttest` benchmark example assumed an optional torch runtime; it now selects the supported GNN-off route.
- Fixed: the old Windows failure caused by the empty benchmark output table; the new assertion guards against recurrence.
- Deferred: deployed-site parity; it requires the authorized merge and deployment path.
- Deferred: final release-validator closeout and any maintainer release decision.

## 8. Consistency Audit

Checked the function default, roxygen example, generated Rd, regression test, exact tarball example run, and rendered candidate reference page together. The package behavior default remains `gnn = TRUE`; the example explicitly uses `FALSE`. The current archive identity, local check logs, current-source defaults receipt, MI matrix manifest, R-hub object identity, Windows logs, site build/crawl receipts, and ledger all refer to the stated commit/hash. Independent reviewers confirmed the platform and local-site claims within their bounded scopes. The current live website remains at an earlier source state.

## 9. What Did Not Go Smoothly

The first R-hub Windows run failed because the example implicitly required torch. The example and test were corrected, then rebuilt into a new artifact. The lane-lease helper could not write its registry outside the task's writable roots. The closeout generator's `new` command targets the brain repository root, outside this task's writable roots. Its combined `check` also runs the brain-wide evidence gate and stops on nine unrelated unfinished brain ledgers; the package report passed the repo-scoped structure validator. The system Git shim also requires the Command Line Tools binary because the Xcode license is not accepted.

## 10. Known Residuals

The fresh site is a local candidate. Its deployed routes and search/sitemap remain unverified until merge and deployment. Pixel-level appearance was not assessed; browser inspection covered rendered text and accessibility-tree content. The MI matrix reuses shared dependency symlinks and is not clean-room dependency installation. It supports termwise fixed-effect pooling, not joint covariance, loading, or variance-component pooling and does not establish broad statistical recovery or interval coverage. R-hub rebuilt the extracted archive source before checking, so its results follow the exact uploaded bytes through extraction and source checkout rather than checking the original compressed file in place. The release remains `NOT_READY` and nothing was submitted to CRAN.

## 11. Team Learning

An exported example should state optional-runtime requirements explicitly. A small benchmark example can use a deterministic runtime-free path while preserving the public function default. Keep package-object readback, extraction, platform rebuild, and platform check as separate links in the artifact provenance chain. Memory receipt: used the pigauto LOAD-FIRST constraints and the exact-artifact release lesson from the memory registry; these kept the local site, deployed site, and frozen archive claims separate. Golden Set: not run because this slice covered docs/example portability and release evidence; no known-mistake-class audit was requested.

## 12. Cross-Product Coverage

Covers ✓: `simulate_benchmark()` documentation example, its generated help, and the tested baseline-only example route; current source defaults audit; optional-backend installed MI workflows on the exact archive; local candidate website routes and discovery; macOS exact-archive check; Windows checks of rebuilt source extracted from the exact uploaded object.

Does NOT cover ✗: current deployed-site behavior; visual pixel appearance; clean-room optional dependency installation; statistical recovery beyond the retained bounded DGPs; joint pooling of covariance, loadings, or variance components; every possible default interaction; CRAN submission, acceptance, or publication; merge or deployment authorization.
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/defaults-current-2026-10-17-9a457eb.log`

~~~~text
ℹ Testing pigauto
✔ | F W  S  OK | Context

⠏ |          0 | cran-audit-defaults                                            
⠇ |         19 | cran-audit-defaults                                            
⠙ |         42 | cran-audit-defaults                                            
⠦ |         67 | cran-audit-defaults                                            
⠹ |         83 | cran-audit-defaults                                            Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").

⠴ |         86 | cran-audit-defaults                                            
⠦ |         87 | cran-audit-defaults                                            pigauto: the exact ("exact") joint baseline route was not usable for a fit (fewer than 2 joint columns (K < 2)); used the per-column baseline for that fit instead. This happens either while predict_method = "auto" (the default) is comparing its "exact" candidate, or when "exact" is itself the resolved route and fell back. This message is shown at most once per R session; pass predict_method = "per_column" to use the per-column path throughout, or predict_method = "exact" explicitly if you want a warning every time this happens.

⠹ |         93 | cran-audit-defaults                                            Auto-detected trait type "count" (log1p-encoded) for integer column(s): count, zi. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").


⠦ |        107 | cran-audit-defaults                                            
⠇ |        129 | cran-audit-defaults                                            
⠼ |        135 | cran-audit-defaults                                            
⠴ |        136 | cran-audit-defaults                                            Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").


⠇ |        139 | cran-audit-defaults                                            Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").


⠏ |        150 | cran-audit-defaults                                            Building phylogenetic graph (reused across folds)...
Rep 1/1, Fold 1/2 (seed=1)...
Building phylogenetic graph (reused across folds)...
Rep 1/1, Fold 1/2 (seed=1)...
=== Replicate 1/1 (seed 1) ===
Fitting baseline...
Fitting pigauto...
=== Replicate 1/1 (seed 1) ===
Fitting baseline...
Fitting pigauto...

⠸ |        164 | cran-audit-defaults                                            draws_method = "auto": using "conformal" because some traits are not continuous (binary: binary). These draws show the spread of plausible values but are not proper multiple imputations; with_imputations() and pool_mi() refuse them. See ?multi_impute.

⠹ |        183 | cran-audit-defaults                                            
⠧ |        198 | cran-audit-defaults                                            
✔ |        202 | cran-audit-defaults [4.8s]

══ Results ═════════════════════════════════════════════════════════════════════
Duration: 4.8 s

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 202 ]
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-as-cran-check.log`

~~~~text
* using log directory ‘/private/tmp/pigauto-ascran-c2e768ad-network/pigauto.Rcheck’
* using R version 4.6.0 (2026-04-24)
* using platform: aarch64-apple-darwin23
* R was compiled by
    Apple clang version 17.0.0 (clang-1700.3.19.1)
    GNU Fortran (GCC) 14.2.0
* running under: macOS Tahoe 26.7.1
* using session charset: UTF-8
* current time: 2026-10-17 16:43:26 UTC
* using options ‘--run-donttest --as-cran’
* checking for file ‘pigauto/DESCRIPTION’ ... OK
* this is package ‘pigauto’ version ‘0.11.0’
* package encoding: UTF-8
* checking CRAN incoming feasibility ... [5s/52s] OK
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking whether package ‘pigauto’ can be installed ... [16s/18s] OK
* checking installed package size ... INFO
  installed size is  6.9Mb
  sub-directories of 1Mb or more:
    data      3.2Mb
    extdata   1.1Mb
    help      1.1Mb
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking ‘build’ directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking code files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking whether startup messages can be suppressed ... OK
* checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... [22s/25s] OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of ‘data’ directory ... OK
* checking data for non-ASCII characters ... OK
* checking LazyData ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking installed files from ‘inst/doc’ ... OK
* checking files in ‘vignettes’ ... OK
* checking examples ... [12s/14s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ... [271s/341s] OK
  Running ‘testthat.R’ [271s/340s]
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [13s/16s] OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-local-check.md`

~~~~text
# Exact local CRAN check for 0.11.0

Date: 2026-10-17. The checked source commit is `9a457eb6b3363d413dbc502de88ff5bc602ed98a` and the frozen archive SHA-256 is `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`, size 5,136,713 bytes.

The exact command was `env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest --timings --output=/private/tmp/pigauto-ascran-c2e768ad-network /private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz`. It completed on macOS arm64 with R 4.6.0. Incoming feasibility, installation, examples including `--run-donttest`, tests, vignette rebuilding, and PDF and HTML manuals completed. Final status: `OK`. Testthat reported 3,195 pass, 0 fail, 161 warnings, and 86 skips.

The exact check log is `exact-artifact-c2e768ad-as-cran-check.log`, SHA-256 `93180d66bf6ee3ddeed26806ef1a998849d9759bc6eea5e2abda6292e4d48478`. Test output is `exact-artifact-c2e768ad-testthat.Rout`, SHA-256 `54590c0e9219591cf4c7e19947d660e599cdf575c0669889623731a76f6621fb`.

The Windows run for predecessor hash `6e8445d6` failed in the `simulate_benchmark()` `\donttest` example because it invoked the optional GNN without an installed runtime. This revision makes the documented example explicitly use `gnn = FALSE`, the default used by `impute()` and `fit_pigauto()`, and adds a nonempty-result assertion to the baseline-only test. The focused `test-new-features` file passed separately; the full test suite in this exact artifact reported 3,195 passes. The new Windows run is tracked separately and had not completed when this receipt was written.
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-source-binding.md`

~~~~text
# Source binding for the 0.11.0 candidate

Date: 2026-10-17. Source commit: `9a457eb6b3363d413dbc502de88ff5bc602ed98a`.

The frozen archive is `/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz`. Its SHA-256 is `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`, size is 5,136,713 bytes, archive inventory has 250 paths, and the sorted path-list SHA-256 is `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`.

The tarball was built with `R CMD build` in a clean detached checkout of the stated commit. A regular-file comparison against `git archive` found 215 files byte-identical, one expected `DESCRIPTION` normalization by `R CMD build`, and 20 build-generated files. No other common-file differences were found. The source tree was clean at build time.

The exact archive was copied to its hash-addressed path and verified byte-identical. R-hub run [39051372400](https://github.com/r-hub2/comfortable-iguana-pigauto/actions/runs/39051372400) received package object `a5ff87d6cc9f5e0157d01dc97456aff5`. The public object was read back to `/private/tmp/pigauto-rhub-upload-c2e768ad.tar.gz`; both files are 5,136,713 bytes with the SHA-256 above. The R-release and R-devel Windows checks were still in progress when this receipt was written.
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-testthat.Rout`

~~~~text

R version 4.6.0 (2026-04-24) -- "Because it was There"
Copyright (C) 2026 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin23

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(testthat)
> library(pigauto)
> 
> test_check("pigauto")
pigauto: the exact ("exact") joint baseline route was not usable for a fit (fewer than 2 joint columns (K < 2)); used the per-column baseline for that fit instead. This happens either while predict_method = "auto" (the default) is comparing its "exact" candidate, or when "exact" is itself the resolved route and fell back. This message is shown at most once per R session; pass predict_method = "per_column" to use the per-column path throughout, or predict_method = "exact" explicitly if you want a warning every time this happens.
Auto-detected trait type "count" (log1p-encoded) for integer column(s): parasites. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(parasites = "continuous").

Auto-detected trait type "count" (log1p-encoded) for integer column(s): parasites. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(parasites = "continuous").

1 tree tip(s) have no trait data and will have all-NA rows.
Auto-detected trait type "count" (log1p-encoded) for integer column(s): counts. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(counts = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): count, zi. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").

Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").

Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous").

Building phylogenetic graph (reused across folds)...
Rep 1/1, Fold 1/2 (seed=1)...
Building phylogenetic graph (reused across folds)...
Rep 1/1, Fold 1/2 (seed=1)...
=== Replicate 1/1 (seed 1) ===
Fitting baseline...
Fitting pigauto...
=== Replicate 1/1 (seed 1) ===
Fitting baseline...
Fitting pigauto...
draws_method = "auto": using "conformal" because some traits are not continuous (binary: binary). These draws show the spread of plausible values but are not proper multiple imputations; with_imputations() and pool_mi() refuse them. See ?multi_impute.
pigauto: the exact ("exact") joint baseline route was not usable for a fit (problem too large, or a numerically unusable Sigma); used the per-column baseline for that fit instead. This happens either while predict_method = "auto" (the default) is comparing its "exact" candidate, or when "exact" is itself the resolved route and fell back. This message is shown at most once per R session; pass predict_method = "per_column" to use the per-column path throughout, or predict_method = "exact" explicitly if you want a warning every time this happens.
Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").

pigauto: the exact ("exact") joint baseline route was not usable for a fit (fewer than 2 joint columns (K < 2)); used the per-column baseline for that fit instead. This happens either while predict_method = "auto" (the default) is comparing its "exact" candidate, or when "exact" is itself the resolved route and fell back. This message is shown at most once per R session; pass predict_method = "per_column" to use the per-column path throughout, or predict_method = "exact" explicitly if you want a warning every time this happens.
Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").
6 tree tip(s) have no trait data and will have all-NA rows.
6 tree tip(s) have no trait data and will have all-NA rows.
gnn = FALSE: mc_dropout draws are BM-posterior draws (no dropout).
Model saved to /var/folders/7x/ytfpq14s0v18frbm9v_w9f4c0000gq/T//Rtmpkc6kW1/working_dir/RtmpDsmQdL/file681e8978f09.pigauto (433,232 bytes)
draws_method = "auto": using "conformal" because `covariates` were supplied. These draws show the spread of plausible values but are not proper multiple imputations; with_imputations() and pool_mi() refuse them. See ?multi_impute.

  ml lambda_hat = 0.005, cv lambda_hat = 0.010 (true mix = 0.2)
Auto-detected trait type "count" (log1p-encoded) for integer column(s): count_trait. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count_trait = "continuous").
Auto-detected trait type "count" (log1p-encoded) for integer column(s): cnt. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(cnt = "continuous").

Model saved to /var/folders/7x/ytfpq14s0v18frbm9v_w9f4c0000gq/T//Rtmpkc6kW1/working_dir/RtmpDsmQdL/file681e1d0817ba.pigauto (3,680 bytes)
Model saved to /var/folders/7x/ytfpq14s0v18frbm9v_w9f4c0000gq/T//Rtmpkc6kW1/working_dir/RtmpDsmQdL/file681e5ff4773b.pigauto (3,680 bytes)
Auto-detected trait type "count" (log1p-encoded) for integer column(s): cnt. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(cnt = "continuous").
c1: chosen=exact exact=0.1817 per_column=0.2215; c2: chosen=per_column exact=0.3039 per_column=0.2961
Auto-detected trait type "count" (log1p-encoded) for integer column(s): cnt. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(cnt = "continuous").

mixed-route seeds (of 12): 1, 2, 4, 6, 7, 9, 10, 11, 12
Note: model was saved with pigauto 0.9.1.9000; current version is 0.11.0
draws_method = "auto": using "conformal" because some traits are not continuous (Trophic.Level: categorical, Primary.Lifestyle: categorical, Migration: ordinal). These draws show the spread of plausible values but are not proper multiple imputations; with_imputations() and pool_mi() refuse them. See ?multi_impute.
Auto-detected trait type "count" (log1p-encoded) for integer column(s): x_cnt. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(x_cnt = "continuous").

Calibrating gates on validation set (baseline only)...
Calibrated gates: Mass=0, Beak.Length_Culmen=0, Tarsus.Length=0, Wing.Length=0, Trophic.Level=Carnivore=0, Trophic.Level=Herbivore=0, Trophic.Level=Omnivore=0, Trophic.Level=Scavenger=0, Primary.Lifestyle=Aerial=0, Primary.Lifestyle=Aquatic=0, Primary.Lifestyle=Generalist=0, Primary.Lifestyle=Insessorial=0, Primary.Lifestyle=Terrestrial=0, Migration=0
Conformal scores (latent scale, method=split): Mass=0.3485, Beak.Length_Culmen=0.7561, Tarsus.Length=0.9727, Wing.Length=0.6671, Migration=2.9202
pigauto: the exact ("exact") joint baseline route was not usable for a fit (fewer than 2 joint columns (K < 2)); used the per-column baseline for that fit instead. This happens either while predict_method = "auto" (the default) is comparing its "exact" candidate, or when "exact" is itself the resolved route and fell back. This message is shown at most once per R session; pass predict_method = "per_column" to use the per-column path throughout, or predict_method = "exact" explicitly if you want a warning every time this happens.
Calibrating gates on validation set (baseline only)...
Calibrated gates: Dummy_Size=0
Conformal scores (latent scale, method=split): Dummy_Size=1.4107
gnn = FALSE: mc_dropout draws are BM-posterior draws (no dropout).
3 tree tip(s) have no trait data and will have all-NA rows.

8 tree tip(s) have no trait data and will have all-NA rows.

gnn = FALSE: mc_dropout draws are BM-posterior draws (no dropout).
[ FAIL 0 | WARN 161 | SKIP 86 | PASS 3195 ]

══ Skipped tests (86) ══════════════════════════════════════════════════════════
• BIEN cache not found (1): 'test-phylo-signal-gate.R:162:3'
• BIEN cache not found -- run script/bench_bien.R once to build (1):
  'test-safety-floor.R:509:3'
• NEWS.md not found (1): 'test-lambda-default.R:224:3'
• On CRAN (54): 'test-draws-auto.R:67:3', 'test-draws-auto.R:82:3',
  'test-gnn-off.R:58:3', 'test-gnn-off.R:104:3', 'test-gnn-off.R:158:3',
  'test-gnn-off.R:192:3', 'test-gnn-off.R:209:3', 'test-gnn-off.R:247:3',
  'test-gnn-off.R:267:3', 'test-gnn-off.R:298:3', 'test-gnn-off.R:311:3',
  'test-gnn-off.R:329:3', 'test-gnn-off.R:352:3', 'test-gnn-off.R:368:3',
  'test-gnn-off.R:397:3', 'test-gnn-off.R:406:3', 'test-gnn-off.R:450:3',
  'test-gnn-off.R:459:3', 'test-gnn-off.R:477:3', 'test-joint-baseline.R:12:3',
  'test-joint-baseline.R:64:3', 'test-joint-baseline.R:107:3',
  'test-joint-baseline.R:152:3', 'test-joint-baseline.R:189:3',
  'test-joint-lambda.R:98:3', 'test-joint-lambda.R:130:3',
  'test-joint-lambda.R:277:3', 'test-joint-refine-iter.R:27:3',
  'test-joint-refine-iter.R:39:3', 'test-joint-refine-iter.R:66:3',
  'test-joint-refine-iter.R:77:3', 'test-joint-refine-iter.R:88:3',
  'test-joint-refine-iter.R:99:3', 'test-joint-refine-iter.R:110:3',
  'test-joint-solver.R:27:3', 'test-joint-solver.R:39:3',
  'test-joint-solver.R:58:3', 'test-joint-solver.R:74:3',
  'test-joint-solver.R:95:3', 'test-lambda-covariates.R:98:3',
  'test-lambda-dispatch.R:343:3', 'test-liability.R:132:3',
  'test-mi-posterior.R:311:3', 'test-mi-posterior.R:487:3',
  'test-mi-posterior.R:556:3', 'test-multi-impute-posterior.R:367:3',
  'test-multi-impute-posterior.R:393:3', 'test-multi-proportion.R:76:3',
  'test-new-features.R:606:3', 'test-new-features.R:629:3',
  'test-pagel-eigendecomp.R:106:3', 'test-pagel-lambda.R:310:3',
  'test-pool-mi-backends.R:35:3', 'test-pool-mi-backends.R:43:3'
• rgbif IS installed; cannot test missing-rgbif branch (1):
  'test-gbif-centroids.R:92:3'
• script/gha helpers are not shipped in the built package (7):
  'test-gha-helpers.R:19:3', 'test-gha-helpers.R:34:3',
  'test-gha-helpers.R:68:3', 'test-gha-helpers.R:104:3',
  'test-gha-helpers.R:147:3', 'test-gha-helpers.R:176:3',
  'test-gha-helpers.R:214:3'
• script/gha/_bace_compat.R not found at ../../script/gha/_bace_compat.R (4):
  'test-bace-compat-eval.R:22:5', 'test-bace-compat-eval.R:85:5',
  'test-bace-compat-eval.R:129:5', 'test-bace-compat-eval.R:177:5'
• slow integration - set NOT_CRAN=TRUE to run (1):
  'test-gbif-centroids.R:226:3'
• slow integration -- NOT_CRAN=TRUE to run (2):
  'test-worldclim-covariates.R:267:3', 'test-worldclim-covariates.R:362:3'
• smcfcs is installed; integration test covers this route (1):
  'test-multi-impute-analysis.R:244:3'
• source file not installed during R CMD check: ../../R/phylo_signal.R (1):
  'test-docs-phytools-warning.R:36:3'
• source file not installed during R CMD check: ../../README.md (3):
  'test-docs-phytools-warning.R:26:3', 'test-docs-phytools-warning.R:43:3',
  'test-docs-phytools-warning.R:53:3'
• source file not installed during R CMD check:
  ../../vignettes/getting-started.Rmd (1): 'test-docs-phytools-warning.R:17:3'
• source-only development archive is excluded from the package tarball (1):
  'test-community-surface.R:193:5'
• source-only pkgdown configuration is excluded from the package tarball (1):
  'test-community-surface.R:75:5'
• source-surface assertions run in the dedicated repository gate (5):
  'test-community-surface.R:10:3', 'test-community-surface.R:57:3',
  'test-community-surface.R:85:3', 'test-community-surface.R:97:3',
  'test-community-surface.R:146:3'
• terra IS installed; cannot test missing-terra branch (1):
  'test-worldclim-covariates.R:175:3'

[ FAIL 0 | WARN 161 | SKIP 86 | PASS 3195 ]
> 
> proc.time()
   user  system elapsed 
226.842  43.687 339.494
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/both-fixtures.log`

~~~~text
tests/testthat/test-pool-mi-adapter-fixtures.R: 50 expectations passed against installed pigauto
tests/testthat/test-mi-provenance.R: 55 expectations passed against installed pigauto
ADAPTER_FIXTURES_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/both-install.log`

~~~~text
* installing *source* package ‘pigauto’ ...
** this is package ‘pigauto’ version ‘0.11.0’
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (pigauto)
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/both-real.log`

~~~~text
REAL_FIT_STATUS backend=drmTMB imputation=1 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=drmTMB imputation=2 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=drmTMB imputation=3 convergence=0 pdHess=TRUE
GAUSSIAN_ORACLE_OK backend=drmTMB fits=3
RUBIN_ORACLE_OK backend=drmTMB terms=3 imputations=3
REAL_FIT_STATUS backend=gllvmTMB imputation=1 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=gllvmTMB imputation=2 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=gllvmTMB imputation=3 convergence=0 pdHess=TRUE
RUBIN_ORACLE_OK backend=gllvmTMB terms=8 imputations=3
script/cran-0.11-integration/test-real-backends.R: 50 expectations passed against installed pigauto
REAL_ADAPTERS_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/drm-only-fixtures.log`

~~~~text
tests/testthat/test-pool-mi-adapter-fixtures.R: 50 expectations passed against installed pigauto
tests/testthat/test-mi-provenance.R: 55 expectations passed against installed pigauto
ADAPTER_FIXTURES_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/drm-only-install.log`

~~~~text
* installing *source* package ‘pigauto’ ...
** this is package ‘pigauto’ version ‘0.11.0’
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (pigauto)
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/drm-only-real.log`

~~~~text
REAL_FIT_STATUS backend=drmTMB imputation=1 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=drmTMB imputation=2 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=drmTMB imputation=3 convergence=0 pdHess=TRUE
GAUSSIAN_ORACLE_OK backend=drmTMB fits=3
RUBIN_ORACLE_OK backend=drmTMB terms=3 imputations=3
script/cran-0.11-integration/test-real-backends.R: 28 expectations passed against installed pigauto
DRM_ONLY_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/gllvm-only-fixtures.log`

~~~~text
tests/testthat/test-pool-mi-adapter-fixtures.R: 50 expectations passed against installed pigauto
tests/testthat/test-mi-provenance.R: 55 expectations passed against installed pigauto
ADAPTER_FIXTURES_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/gllvm-only-install.log`

~~~~text
* installing *source* package ‘pigauto’ ...
** this is package ‘pigauto’ version ‘0.11.0’
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (pigauto)
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/gllvm-only-real.log`

~~~~text
REAL_FIT_STATUS backend=gllvmTMB imputation=1 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=gllvmTMB imputation=2 convergence=0 pdHess=TRUE
REAL_FIT_STATUS backend=gllvmTMB imputation=3 convergence=0 pdHess=TRUE
RUBIN_ORACLE_OK backend=gllvmTMB terms=8 imputations=3
script/cran-0.11-integration/test-real-backends.R: 22 expectations passed against installed pigauto
GLLVM_ONLY_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/manifest.json`

~~~~text
{
  "archive": "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz",
  "archive_sha256": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
  "source_commit": "9a457eb6b3363d413dbc502de88ff5bc602ed98a",
  "working_directory": "/private/tmp/pigauto-release-evidence-780c084",
  "executed_sources": {
    "script/cran-0.11-integration/check-adapters.R": {
      "path": "script/cran-0.11-integration/check-adapters.R",
      "sha256": "d2b0b5f01e3470ce1f0f0d2116cadbd1a9a2beffcd864a5644e71c0b7d53fb7d"
    },
    "script/cran-0.11-integration/test-real-backends.R": {
      "path": "script/cran-0.11-integration/test-real-backends.R",
      "sha256": "7edb499a09820ecd9c8efc6799f2adf4233c824f064bc88462b9f1720098b6be"
    },
    "script/cran-0.11-integration/reload-fit.R": {
      "path": "script/cran-0.11-integration/reload-fit.R",
      "sha256": "8db40649b97f074eefb71cfd591453d9540f8b4518b464ac4a6d13472f68401d"
    },
    "tests/testthat/test-pool-mi-adapter-fixtures.R": {
      "path": "tests/testthat/test-pool-mi-adapter-fixtures.R",
      "sha256": "1599e96a85be78afec4fe05fc23cc02f9c0f5ac99e877bb31ff87a358e20a0b1"
    },
    "tests/testthat/test-mi-provenance.R": {
      "path": "tests/testthat/test-mi-provenance.R",
      "sha256": "585299f93cb4c83013337d1043035bf0dbe7af40389c1dc3075507391591f67a"
    }
  },
  "harness": "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
  "harness_sha256": "d2b0b5f01e3470ce1f0f0d2116cadbd1a9a2beffcd864a5644e71c0b7d53fb7d",
  "runner": "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/run-exact-artifact-matrix.py",
  "runner_sha256": "6d6cc86d7c17f1e5e91c5e4a8e9841eef2dadbc793c2775657874fe79cda2bf7",
  "created_utc": "2026-10-17T16:59:50.857514+00:00",
  "r_version": "R version 4.6.0 (2026-04-24) -- \"Because it was There\"",
  "cells": [
    {
      "cell": "neither",
      "mode": "neither",
      "library": "/private/tmp/pigauto-mi-matrix-c2e768ad/neither",
      "expected_optional_packages": [],
      "dependency_symlinks": {
        "AGHmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/AGHmatrix",
        "BACE": "/Users/z3437171/Library/R/arm64/4.6/library/BACE",
        "BH": "/Users/z3437171/Library/R/arm64/4.6/library/BH",
        "BayesLogit": "/Users/z3437171/Library/R/arm64/4.6/library/BayesLogit",
        "BiocGenerics": "/Users/z3437171/Library/R/arm64/4.6/library/BiocGenerics",
        "BiocManager": "/Users/z3437171/Library/R/arm64/4.6/library/BiocManager",
        "BiocVersion": "/Users/z3437171/Library/R/arm64/4.6/library/BiocVersion",
        "Brobdingnag": "/Users/z3437171/Library/R/arm64/4.6/library/Brobdingnag",
        "DBI": "/Users/z3437171/Library/R/arm64/4.6/library/DBI",
        "DEoptim": "/Users/z3437171/Library/R/arm64/4.6/library/DEoptim",
        "DHARMa": "/Users/z3437171/Library/R/arm64/4.6/library/DHARMa",
        "DT": "/Users/z3437171/Library/R/arm64/4.6/library/DT",
        "Deriv": "/Users/z3437171/Library/R/arm64/4.6/library/Deriv",
        "DiagrammeR": "/Users/z3437171/Library/R/arm64/4.6/library/DiagrammeR",
        "FNN": "/Users/z3437171/Library/R/arm64/4.6/library/FNN",
        "Formula": "/Users/z3437171/Library/R/arm64/4.6/library/Formula",
        "GPArotation": "/Users/z3437171/Library/R/arm64/4.6/library/GPArotation",
        "HSAUR3": "/Users/z3437171/Library/R/arm64/4.6/library/HSAUR3",
        "Hmsc": "/Users/z3437171/Library/R/arm64/4.6/library/Hmsc",
        "JuliaCall": "/Users/z3437171/Library/R/arm64/4.6/library/JuliaCall",
        "MCMCglmm": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCglmm",
        "MCMCpack": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCpack",
        "MatrixModels": "/Users/z3437171/Library/R/arm64/4.6/library/MatrixModels",
        "MuMIn": "/Users/z3437171/Library/R/arm64/4.6/library/MuMIn",
        "QGglmm": "/Users/z3437171/Library/R/arm64/4.6/library/QGglmm",
        "QuickJSR": "/Users/z3437171/Library/R/arm64/4.6/library/QuickJSR",
        "R.cache": "/Users/z3437171/Library/R/arm64/4.6/library/R.cache",
        "R.methodsS3": "/Users/z3437171/Library/R/arm64/4.6/library/R.methodsS3",
        "R.oo": "/Users/z3437171/Library/R/arm64/4.6/library/R.oo",
        "R.rsp": "/Users/z3437171/Library/R/arm64/4.6/library/R.rsp",
        "R.utils": "/Users/z3437171/Library/R/arm64/4.6/library/R.utils",
        "R6": "/Users/z3437171/Library/R/arm64/4.6/library/R6",
        "RColorBrewer": "/Users/z3437171/Library/R/arm64/4.6/library/RColorBrewer",
        "ROI": "/Users/z3437171/Library/R/arm64/4.6/library/ROI",
        "ROI.plugin.lpsolve": "/Users/z3437171/Library/R/arm64/4.6/library/ROI.plugin.lpsolve",
        "RSpectra": "/Users/z3437171/Library/R/arm64/4.6/library/RSpectra",
        "RTMB": "/Users/z3437171/Library/R/arm64/4.6/library/RTMB",
        "RTMBdist": "/Users/z3437171/Library/R/arm64/4.6/library/RTMBdist",
        "RWiener": "/Users/z3437171/Library/R/arm64/4.6/library/RWiener",
        "Rcpp": "/Users/z3437171/Library/R/arm64/4.6/library/Rcpp",
        "RcppArmadillo": "/Users/z3437171/Library/R/arm64/4.6/library/RcppArmadillo",
        "RcppEigen": "/Users/z3437171/Library/R/arm64/4.6/library/RcppEigen",
        "RcppParallel": "/Users/z3437171/Library/R/arm64/4.6/library/RcppParallel",
        "RcppProgress": "/Users/z3437171/Library/R/arm64/4.6/library/RcppProgress",
        "Rdpack": "/Users/z3437171/Library/R/arm64/4.6/library/Rdpack",
        "Rfast": "/Users/z3437171/Library/R/arm64/4.6/library/Rfast",
        "Rmpfr": "/Users/z3437171/Library/R/arm64/4.6/library/Rmpfr",
        "Rphylopars": "/Users/z3437171/Library/R/arm64/4.6/library/Rphylopars",
        "S7": "/Users/z3437171/Library/R/arm64/4.6/library/S7",
        "SimDesign": "/Users/z3437171/Library/R/arm64/4.6/library/SimDesign",
        "SparseM": "/Users/z3437171/Library/R/arm64/4.6/library/SparseM",
        "StanHeaders": "/Users/z3437171/Library/R/arm64/4.6/library/StanHeaders",
        "TH.data": "/Users/z3437171/Library/R/arm64/4.6/library/TH.data",
        "TMB": "/Users/z3437171/Library/R/arm64/4.6/library/TMB",
        "V8": "/Users/z3437171/Library/R/arm64/4.6/library/V8",
        "VGAM": "/Users/z3437171/Library/R/arm64/4.6/library/VGAM",
        "XML": "/Users/z3437171/Library/R/arm64/4.6/library/XML",
        "_cache": "/Users/z3437171/Library/R/arm64/4.6/library/_cache",
        "abind": "/Users/z3437171/Library/R/arm64/4.6/library/abind",
        "ade4": "/Users/z3437171/Library/R/arm64/4.6/library/ade4",
        "alabama": "/Users/z3437171/Library/R/arm64/4.6/library/alabama",
        "ape": "/Users/z3437171/Library/R/arm64/4.6/library/ape",
        "arm": "/Users/z3437171/Library/R/arm64/4.6/library/arm",
        "askpass": "/Users/z3437171/Library/R/arm64/4.6/library/askpass",
        "assertthat": "/Users/z3437171/Library/R/arm64/4.6/library/assertthat",
        "audio": "/Users/z3437171/Library/R/arm64/4.6/library/audio",
        "backports": "/Users/z3437171/Library/R/arm64/4.6/library/backports",
        "base64enc": "/Users/z3437171/Library/R/arm64/4.6/library/base64enc",
        "bayeslincom": "/Users/z3437171/Library/R/arm64/4.6/library/bayeslincom",
        "bayesplot": "/Users/z3437171/Library/R/arm64/4.6/library/bayesplot",
        "bayestestR": "/Users/z3437171/Library/R/arm64/4.6/library/bayestestR",
        "bbmle": "/Users/z3437171/Library/R/arm64/4.6/library/bbmle",
        "bdsmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/bdsmatrix",
        "beepr": "/Users/z3437171/Library/R/arm64/4.6/library/beepr",
        "beeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/beeswarm",
        "betareg": "/Users/z3437171/Library/R/arm64/4.6/library/betareg",
        "bit": "/Users/z3437171/Library/R/arm64/4.6/library/bit",
        "bit64": "/Users/z3437171/Library/R/arm64/4.6/library/bit64",
        "blob": "/Users/z3437171/Library/R/arm64/4.6/library/blob",
        "blsmeta": "/Users/z3437171/Library/R/arm64/4.6/library/blsmeta",
        "brew": "/Users/z3437171/Library/R/arm64/4.6/library/brew",
        "brglm2": "/Users/z3437171/Library/R/arm64/4.6/library/brglm2",
        "bridgesampling": "/Users/z3437171/Library/R/arm64/4.6/library/bridgesampling",
        "brio": "/Users/z3437171/Library/R/arm64/4.6/library/brio",
        "brms": "/Users/z3437171/Library/R/arm64/4.6/library/brms",
        "broom": "/Users/z3437171/Library/R/arm64/4.6/library/broom",
        "broom.mixed": "/Users/z3437171/Library/R/arm64/4.6/library/broom.mixed",
        "bslib": "/Users/z3437171/Library/R/arm64/4.6/library/bslib",
        "cachem": "/Users/z3437171/Library/R/arm64/4.6/library/cachem",
        "callr": "/Users/z3437171/Library/R/arm64/4.6/library/callr",
        "car": "/Users/z3437171/Library/R/arm64/4.6/library/car",
        "carData": "/Users/z3437171/Library/R/arm64/4.6/library/carData",
        "castor": "/Users/z3437171/Library/R/arm64/4.6/library/castor",
        "cellranger": "/Users/z3437171/Library/R/arm64/4.6/library/cellranger",
        "checkmate": "/Users/z3437171/Library/R/arm64/4.6/library/checkmate",
        "circular": "/Users/z3437171/Library/R/arm64/4.6/library/circular",
        "classInt": "/Users/z3437171/Library/R/arm64/4.6/library/classInt",
        "cli": "/Users/z3437171/Library/R/arm64/4.6/library/cli",
        "clipr": "/Users/z3437171/Library/R/arm64/4.6/library/clipr",
        "clue": "/Users/z3437171/Library/R/arm64/4.6/library/clue",
        "clusterGeneration": "/Users/z3437171/Library/R/arm64/4.6/library/clusterGeneration",
        "cmdstanr": "/Users/z3437171/Library/R/arm64/4.6/library/cmdstanr",
        "coda": "/Users/z3437171/Library/R/arm64/4.6/library/coda",
        "collections": "/Users/z3437171/Library/R/arm64/4.6/library/collections",
        "colorspace": "/Users/z3437171/Library/R/arm64/4.6/library/colorspace",
        "colourpicker": "/Users/z3437171/Library/R/arm64/4.6/library/colourpicker",
        "combinat": "/Users/z3437171/Library/R/arm64/4.6/library/combinat",
        "commonmark": "/Users/z3437171/Library/R/arm64/4.6/library/commonmark",
        "conflicted": "/Users/z3437171/Library/R/arm64/4.6/library/conflicted",
        "coro": "/Users/z3437171/Library/R/arm64/4.6/library/coro",
        "corpcor": "/Users/z3437171/Library/R/arm64/4.6/library/corpcor",
        "covr": "/Users/z3437171/Library/R/arm64/4.6/library/covr",
        "cowplot": "/Users/z3437171/Library/R/arm64/4.6/library/cowplot",
        "cpp11": "/Users/z3437171/Library/R/arm64/4.6/library/cpp11",
        "crayon": "/Users/z3437171/Library/R/arm64/4.6/library/crayon",
        "credentials": "/Users/z3437171/Library/R/arm64/4.6/library/credentials",
        "crosstalk": "/Users/z3437171/Library/R/arm64/4.6/library/crosstalk",
        "crul": "/Users/z3437171/Library/R/arm64/4.6/library/crul",
        "cubature": "/Users/z3437171/Library/R/arm64/4.6/library/cubature",
        "curl": "/Users/z3437171/Library/R/arm64/4.6/library/curl",
        "dagitty": "/Users/z3437171/Library/R/arm64/4.6/library/dagitty",
        "data.table": "/Users/z3437171/Library/R/arm64/4.6/library/data.table",
        "datawizard": "/Users/z3437171/Library/R/arm64/4.6/library/datawizard",
        "dbplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dbplyr",
        "dcurver": "/Users/z3437171/Library/R/arm64/4.6/library/dcurver",
        "deldir": "/Users/z3437171/Library/R/arm64/4.6/library/deldir",
        "desc": "/Users/z3437171/Library/R/arm64/4.6/library/desc",
        "detectseparation": "/Users/z3437171/Library/R/arm64/4.6/library/detectseparation",
        "devtools": "/Users/z3437171/Library/R/arm64/4.6/library/devtools",
        "diffobj": "/Users/z3437171/Library/R/arm64/4.6/library/diffobj",
        "digest": "/Users/z3437171/Library/R/arm64/4.6/library/digest",
        "distributional": "/Users/z3437171/Library/R/arm64/4.6/library/distributional",
        "doBy": "/Users/z3437171/Library/R/arm64/4.6/library/doBy",
        "doFuture": "/Users/z3437171/Library/R/arm64/4.6/library/doFuture",
        "doParallel": "/Users/z3437171/Library/R/arm64/4.6/library/doParallel",
        "doRNG": "/Users/z3437171/Library/R/arm64/4.6/library/doRNG",
        "dotCall64": "/Users/z3437171/Library/R/arm64/4.6/library/dotCall64",
        "downlit": "/Users/z3437171/Library/R/arm64/4.6/library/downlit",
        "dplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dplyr",
        "drmSEM": "/Users/z3437171/Library/R/arm64/4.6/library/drmSEM",
        "dtplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dtplyr",
        "dygraphs": "/Users/z3437171/Library/R/arm64/4.6/library/dygraphs",
        "e1071": "/Users/z3437171/Library/R/arm64/4.6/library/e1071",
        "ellipsis": "/Users/z3437171/Library/R/arm64/4.6/library/ellipsis",
        "emmeans": "/Users/z3437171/Library/R/arm64/4.6/library/emmeans",
        "enhancer": "/Users/z3437171/Library/R/arm64/4.6/library/enhancer",
        "enrichwith": "/Users/z3437171/Library/R/arm64/4.6/library/enrichwith",
        "estimability": "/Users/z3437171/Library/R/arm64/4.6/library/estimability",
        "evaluate": "/Users/z3437171/Library/R/arm64/4.6/library/evaluate",
        "evd": "/Users/z3437171/Library/R/arm64/4.6/library/evd",
        "expm": "/Users/z3437171/Library/R/arm64/4.6/library/expm",
        "extraDistr": "/Users/z3437171/Library/R/arm64/4.6/library/extraDistr",
        "fansi": "/Users/z3437171/Library/R/arm64/4.6/library/fansi",
        "farver": "/Users/z3437171/Library/R/arm64/4.6/library/farver",
        "fastGHQuad": "/Users/z3437171/Library/R/arm64/4.6/library/fastGHQuad",
        "fastmap": "/Users/z3437171/Library/R/arm64/4.6/library/fastmap",
        "fastmatch": "/Users/z3437171/Library/R/arm64/4.6/library/fastmatch",
        "fields": "/Users/z3437171/Library/R/arm64/4.6/library/fields",
        "fishMod": "/Users/z3437171/Library/R/arm64/4.6/library/fishMod",
        "flexmix": "/Users/z3437171/Library/R/arm64/4.6/library/flexmix",
        "fmesher": "/Users/z3437171/Library/R/arm64/4.6/library/fmesher",
        "fontawesome": "/Users/z3437171/Library/R/arm64/4.6/library/fontawesome",
        "forcats": "/Users/z3437171/Library/R/arm64/4.6/library/forcats",
        "foreach": "/Users/z3437171/Library/R/arm64/4.6/library/foreach",
        "forecast": "/Users/z3437171/Library/R/arm64/4.6/library/forecast",
        "fracdiff": "/Users/z3437171/Library/R/arm64/4.6/library/fracdiff",
        "freqTLS": "/Users/z3437171/Library/R/arm64/4.6/library/freqTLS",
        "fs": "/Users/z3437171/Library/R/arm64/4.6/library/fs",
        "furrr": "/Users/z3437171/Library/R/arm64/4.6/library/furrr",
        "future": "/Users/z3437171/Library/R/arm64/4.6/library/future",
        "future.apply": "/Users/z3437171/Library/R/arm64/4.6/library/future.apply",
        "future.mirai": "/Users/z3437171/Library/R/arm64/4.6/library/future.mirai",
        "galamm": "/Users/z3437171/Library/R/arm64/4.6/library/galamm",
        "gamlss.dist": "/Users/z3437171/Library/R/arm64/4.6/library/gamlss.dist",
        "gamm4": "/Users/z3437171/Library/R/arm64/4.6/library/gamm4",
        "gap": "/Users/z3437171/Library/R/arm64/4.6/library/gap",
        "gap.datasets": "/Users/z3437171/Library/R/arm64/4.6/library/gap.datasets",
        "gargle": "/Users/z3437171/Library/R/arm64/4.6/library/gargle",
        "generics": "/Users/z3437171/Library/R/arm64/4.6/library/generics",
        "geojsonsf": "/Users/z3437171/Library/R/arm64/4.6/library/geojsonsf",
        "geometries": "/Users/z3437171/Library/R/arm64/4.6/library/geometries",
        "gert": "/Users/z3437171/Library/R/arm64/4.6/library/gert",
        "ggbeeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/ggbeeswarm",
        "ggdist": "/Users/z3437171/Library/R/arm64/4.6/library/ggdist",
        "ggfittext": "/Users/z3437171/Library/R/arm64/4.6/library/ggfittext",
        "ggforce": "/Users/z3437171/Library/R/arm64/4.6/library/ggforce",
        "gggenes": "/Users/z3437171/Library/R/arm64/4.6/library/gggenes",
        "ggh4x": "/Users/z3437171/Library/R/arm64/4.6/library/ggh4x",
        "ggm": "/Users/z3437171/Library/R/arm64/4.6/library/ggm",
        "ggokabeito": "/Users/z3437171/Library/R/arm64/4.6/library/ggokabeito",
        "ggplot2": "/Users/z3437171/Library/R/arm64/4.6/library/ggplot2",
        "ggraph": "/Users/z3437171/Library/R/arm64/4.6/library/ggraph",
        "ggrepel": "/Users/z3437171/Library/R/arm64/4.6/library/ggrepel",
        "ggridges": "/Users/z3437171/Library/R/arm64/4.6/library/ggridges",
        "gh": "/Users/z3437171/Library/R/arm64/4.6/library/gh",
        "gitcreds": "/Users/z3437171/Library/R/arm64/4.6/library/gitcreds",
        "glassoFast": "/Users/z3437171/Library/R/arm64/4.6/library/glassoFast",
        "gllvm": "/Users/z3437171/Library/R/arm64/4.6/library/gllvm",
        "glmmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/glmmTMB",
        "glmnet": "/Users/z3437171/Library/R/arm64/4.6/library/glmnet",
        "globals": "/Users/z3437171/Library/R/arm64/4.6/library/globals",
        "glue": "/Users/z3437171/Library/R/arm64/4.6/library/glue",
        "gmp": "/Users/z3437171/Library/R/arm64/4.6/library/gmp",
        "googledrive": "/Users/z3437171/Library/R/arm64/4.6/library/googledrive",
        "googlesheets4": "/Users/z3437171/Library/R/arm64/4.6/library/googlesheets4",
        "graph": "/Users/z3437171/Library/R/arm64/4.6/library/graph",
        "graphlayouts": "/Users/z3437171/Library/R/arm64/4.6/library/graphlayouts",
        "gratia": "/Users/z3437171/Library/R/arm64/4.6/library/gratia",
        "gridExtra": "/Users/z3437171/Library/R/arm64/4.6/library/gridExtra",
        "gridtext": "/Users/z3437171/Library/R/arm64/4.6/library/gridtext",
        "gsl": "/Users/z3437171/Library/R/arm64/4.6/library/gsl",
        "gtable": "/Users/z3437171/Library/R/arm64/4.6/library/gtable",
        "gtools": "/Users/z3437171/Library/R/arm64/4.6/library/gtools",
        "haven": "/Users/z3437171/Library/R/arm64/4.6/library/haven",
        "here": "/Users/z3437171/Library/R/arm64/4.6/library/here",
        "highr": "/Users/z3437171/Library/R/arm64/4.6/library/highr",
        "hms": "/Users/z3437171/Library/R/arm64/4.6/library/hms",
        "hsquared": "/Users/z3437171/Library/R/arm64/4.6/library/hsquared",
        "htmltools": "/Users/z3437171/Library/R/arm64/4.6/library/htmltools",
        "htmlwidgets": "/Users/z3437171/Library/R/arm64/4.6/library/htmlwidgets",
        "httpcode": "/Users/z3437171/Library/R/arm64/4.6/library/httpcode",
        "httpuv": "/Users/z3437171/Library/R/arm64/4.6/library/httpuv",
        "httr": "/Users/z3437171/Library/R/arm64/4.6/library/httr",
        "httr2": "/Users/z3437171/Library/R/arm64/4.6/library/httr2",
        "hunspell": "/Users/z3437171/Library/R/arm64/4.6/library/hunspell",
        "ids": "/Users/z3437171/Library/R/arm64/4.6/library/ids",
        "igraph": "/Users/z3437171/Library/R/arm64/4.6/library/igraph",
        "ini": "/Users/z3437171/Library/R/arm64/4.6/library/ini",
        "inline": "/Users/z3437171/Library/R/arm64/4.6/library/inline",
        "insight": "/Users/z3437171/Library/R/arm64/4.6/library/insight",
        "isoband": "/Users/z3437171/Library/R/arm64/4.6/library/isoband",
        "iterators": "/Users/z3437171/Library/R/arm64/4.6/library/iterators",
        "itertools": "/Users/z3437171/Library/R/arm64/4.6/library/itertools",
        "janitor": "/Users/z3437171/Library/R/arm64/4.6/library/janitor",
        "jomo": "/Users/z3437171/Library/R/arm64/4.6/library/jomo",
        "jpeg": "/Users/z3437171/Library/R/arm64/4.6/library/jpeg",
        "jquerylib": "/Users/z3437171/Library/R/arm64/4.6/library/jquerylib",
        "jsonify": "/Users/z3437171/Library/R/arm64/4.6/library/jsonify",
        "jsonlite": "/Users/z3437171/Library/R/arm64/4.6/library/jsonlite",
        "kernelshap": "/Users/z3437171/Library/R/arm64/4.6/library/kernelshap",
        "knitr": "/Users/z3437171/Library/R/arm64/4.6/library/knitr",
        "labeling": "/Users/z3437171/Library/R/arm64/4.6/library/labeling",
        "languageserver": "/Users/z3437171/Library/R/arm64/4.6/library/languageserver",
        "later": "/Users/z3437171/Library/R/arm64/4.6/library/later",
        "latex2exp": "/Users/z3437171/Library/R/arm64/4.6/library/latex2exp",
        "lavaan": "/Users/z3437171/Library/R/arm64/4.6/library/lavaan",
        "lazyeval": "/Users/z3437171/Library/R/arm64/4.6/library/lazyeval",
        "leafem": "/Users/z3437171/Library/R/arm64/4.6/library/leafem",
        "leaflet": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet",
        "leaflet.providers": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet.providers",
        "leafpop": "/Users/z3437171/Library/R/arm64/4.6/library/leafpop",
        "lifecycle": "/Users/z3437171/Library/R/arm64/4.6/library/lifecycle",
        "lintr": "/Users/z3437171/Library/R/arm64/4.6/library/lintr",
        "listenv": "/Users/z3437171/Library/R/arm64/4.6/library/listenv",
        "litedown": "/Users/z3437171/Library/R/arm64/4.6/library/litedown",
        "lme4": "/Users/z3437171/Library/R/arm64/4.6/library/lme4",
        "lmtest": "/Users/z3437171/Library/R/arm64/4.6/library/lmtest",
        "loo": "/Users/z3437171/Library/R/arm64/4.6/library/loo",
        "lpSolveAPI": "/Users/z3437171/Library/R/arm64/4.6/library/lpSolveAPI",
        "lubridate": "/Users/z3437171/Library/R/arm64/4.6/library/lubridate",
        "magrittr": "/Users/z3437171/Library/R/arm64/4.6/library/magrittr",
        "maps": "/Users/z3437171/Library/R/arm64/4.6/library/maps",
        "mapview": "/Users/z3437171/Library/R/arm64/4.6/library/mapview",
        "markdown": "/Users/z3437171/Library/R/arm64/4.6/library/markdown",
        "mathjaxr": "/Users/z3437171/Library/R/arm64/4.6/library/mathjaxr",
        "matrixStats": "/Users/z3437171/Library/R/arm64/4.6/library/matrixStats",
        "mclogit": "/Users/z3437171/Library/R/arm64/4.6/library/mclogit",
        "mcmc": "/Users/z3437171/Library/R/arm64/4.6/library/mcmc",
        "memisc": "/Users/z3437171/Library/R/arm64/4.6/library/memisc",
        "memoise": "/Users/z3437171/Library/R/arm64/4.6/library/memoise",
        "metadat": "/Users/z3437171/Library/R/arm64/4.6/library/metadat",
        "metafor": "/Users/z3437171/Library/R/arm64/4.6/library/metafor",
        "metaplus": "/Users/z3437171/Library/R/arm64/4.6/library/metaplus",
        "mice": "/Users/z3437171/Library/R/arm64/4.6/library/mice",
        "microbenchmark": "/Users/z3437171/Library/R/arm64/4.6/library/microbenchmark",
        "mime": "/Users/z3437171/Library/R/arm64/4.6/library/mime",
        "miniUI": "/Users/z3437171/Library/R/arm64/4.6/library/miniUI",
        "minqa": "/Users/z3437171/Library/R/arm64/4.6/library/minqa",
        "mirai": "/Users/z3437171/Library/R/arm64/4.6/library/mirai",
        "mirt": "/Users/z3437171/Library/R/arm64/4.6/library/mirt",
        "missForest": "/Users/z3437171/Library/R/arm64/4.6/library/missForest",
        "mitml": "/Users/z3437171/Library/R/arm64/4.6/library/mitml",
        "mnormt": "/Users/z3437171/Library/R/arm64/4.6/library/mnormt",
        "modelr": "/Users/z3437171/Library/R/arm64/4.6/library/modelr",
        "modeltools": "/Users/z3437171/Library/R/arm64/4.6/library/modeltools",
        "movMF": "/Users/z3437171/Library/R/arm64/4.6/library/movMF",
        "msm": "/Users/z3437171/Library/R/arm64/4.6/library/msm",
        "multcomp": "/Users/z3437171/Library/R/arm64/4.6/library/multcomp",
        "mvMORPH": "/Users/z3437171/Library/R/arm64/4.6/library/mvMORPH",
        "mvnfast": "/Users/z3437171/Library/R/arm64/4.6/library/mvnfast",
        "mvtnorm": "/Users/z3437171/Library/R/arm64/4.6/library/mvtnorm",
        "nadiv": "/Users/z3437171/Library/R/arm64/4.6/library/nadiv",
        "nanonext": "/Users/z3437171/Library/R/arm64/4.6/library/nanonext",
        "naturalsort": "/Users/z3437171/Library/R/arm64/4.6/library/naturalsort",
        "ncf": "/Users/z3437171/Library/R/arm64/4.6/library/ncf",
        "nleqslv": "/Users/z3437171/Library/R/arm64/4.6/library/nleqslv",
        "nloptr": "/Users/z3437171/Library/R/arm64/4.6/library/nloptr",
        "numDeriv": "/Users/z3437171/Library/R/arm64/4.6/library/numDeriv",
        "oai": "/Users/z3437171/Library/R/arm64/4.6/library/oai",
        "openssl": "/Users/z3437171/Library/R/arm64/4.6/library/openssl",
        "optimParallel": "/Users/z3437171/Library/R/arm64/4.6/library/optimParallel",
        "orchaRd": "/Users/z3437171/Library/R/arm64/4.6/library/orchaRd",
        "ordinal": "/Users/z3437171/Library/R/arm64/4.6/library/ordinal",
        "orthopolynom": "/Users/z3437171/Library/R/arm64/4.6/library/orthopolynom",
        "otel": "/Users/z3437171/Library/R/arm64/4.6/library/otel",
        "pROC": "/Users/z3437171/Library/R/arm64/4.6/library/pROC",
        "pak": "/Users/z3437171/Library/R/arm64/4.6/library/pak",
        "palmerpenguins": "/Users/z3437171/Library/R/arm64/4.6/library/palmerpenguins",
        "pan": "/Users/z3437171/Library/R/arm64/4.6/library/pan",
        "parallelly": "/Users/z3437171/Library/R/arm64/4.6/library/parallelly",
        "patchwork": "/Users/z3437171/Library/R/arm64/4.6/library/patchwork",
        "pbapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbapply",
        "pbivnorm": "/Users/z3437171/Library/R/arm64/4.6/library/pbivnorm",
        "pbkrtest": "/Users/z3437171/Library/R/arm64/4.6/library/pbkrtest",
        "pbmcapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbmcapply",
        "pedigreemm": "/Users/z3437171/Library/R/arm64/4.6/library/pedigreemm",
        "performance": "/Users/z3437171/Library/R/arm64/4.6/library/performance",
        "permute": "/Users/z3437171/Library/R/arm64/4.6/library/permute",
        "phangorn": "/Users/z3437171/Library/R/arm64/4.6/library/phangorn",
        "phylolm": "/Users/z3437171/Library/R/arm64/4.6/library/phylolm",
        "phylopath": "/Users/z3437171/Library/R/arm64/4.6/library/phylopath",
        "phytools": "/Users/z3437171/Library/R/arm64/4.6/library/phytools",
        "piecewiseSEM": "/Users/z3437171/Library/R/arm64/4.6/library/piecewiseSEM",
        "pillar": "/Users/z3437171/Library/R/arm64/4.6/library/pillar",
        "pixmap": "/Users/z3437171/Library/R/arm64/4.6/library/pixmap",
        "pkgbuild": "/Users/z3437171/Library/R/arm64/4.6/library/pkgbuild",
        "pkgconfig": "/Users/z3437171/Library/R/arm64/4.6/library/pkgconfig",
        "pkgdown": "/Users/z3437171/Library/R/arm64/4.6/library/pkgdown",
        "pkgload": "/Users/z3437171/Library/R/arm64/4.6/library/pkgload",
        "plotly": "/Users/z3437171/Library/R/arm64/4.6/library/plotly",
        "plyr": "/Users/z3437171/Library/R/arm64/4.6/library/plyr",
        "png": "/Users/z3437171/Library/R/arm64/4.6/library/png",
        "polyclip": "/Users/z3437171/Library/R/arm64/4.6/library/polyclip",
        "polynom": "/Users/z3437171/Library/R/arm64/4.6/library/polynom",
        "posterior": "/Users/z3437171/Library/R/arm64/4.6/library/posterior",
        "pracma": "/Users/z3437171/Library/R/arm64/4.6/library/pracma",
        "praise": "/Users/z3437171/Library/R/arm64/4.6/library/praise",
        "prettyunits": "/Users/z3437171/Library/R/arm64/4.6/library/prettyunits",
        "priorsense": "/Users/z3437171/Library/R/arm64/4.6/library/priorsense",
        "processx": "/Users/z3437171/Library/R/arm64/4.6/library/processx",
        "profvis": "/Users/z3437171/Library/R/arm64/4.6/library/profvis",
        "progress": "/Users/z3437171/Library/R/arm64/4.6/library/progress",
        "progressr": "/Users/z3437171/Library/R/arm64/4.6/library/progressr",
        "projpred": "/Users/z3437171/Library/R/arm64/4.6/library/projpred",
        "promises": "/Users/z3437171/Library/R/arm64/4.6/library/promises",
        "proxy": "/Users/z3437171/Library/R/arm64/4.6/library/proxy",
        "ps": "/Users/z3437171/Library/R/arm64/4.6/library/ps",
        "psymetadata": "/Users/z3437171/Library/R/arm64/4.6/library/psymetadata",
        "purrr": "/Users/z3437171/Library/R/arm64/4.6/library/purrr",
        "qgam": "/Users/z3437171/Library/R/arm64/4.6/library/qgam",
        "qs2": "/Users/z3437171/Library/R/arm64/4.6/library/qs2",
        "quadprog": "/Users/z3437171/Library/R/arm64/4.6/library/quadprog",
        "quantreg": "/Users/z3437171/Library/R/arm64/4.6/library/quantreg",
        "rWishart": "/Users/z3437171/Library/R/arm64/4.6/library/rWishart",
        "ragg": "/Users/z3437171/Library/R/arm64/4.6/library/ragg",
        "randomForest": "/Users/z3437171/Library/R/arm64/4.6/library/randomForest",
        "ranger": "/Users/z3437171/Library/R/arm64/4.6/library/ranger",
        "rapidjsonr": "/Users/z3437171/Library/R/arm64/4.6/library/rapidjsonr",
        "rappdirs": "/Users/z3437171/Library/R/arm64/4.6/library/rappdirs",
        "raster": "/Users/z3437171/Library/R/arm64/4.6/library/raster",
        "rbibutils": "/Users/z3437171/Library/R/arm64/4.6/library/rbibutils",
        "rcmdcheck": "/Users/z3437171/Library/R/arm64/4.6/library/rcmdcheck",
        "readr": "/Users/z3437171/Library/R/arm64/4.6/library/readr",
        "readxl": "/Users/z3437171/Library/R/arm64/4.6/library/readxl",
        "reformulas": "/Users/z3437171/Library/R/arm64/4.6/library/reformulas",
        "registry": "/Users/z3437171/Library/R/arm64/4.6/library/registry",
        "rematch": "/Users/z3437171/Library/R/arm64/4.6/library/rematch",
        "rematch2": "/Users/z3437171/Library/R/arm64/4.6/library/rematch2",
        "remotes": "/Users/z3437171/Library/R/arm64/4.6/library/remotes",
        "rentrez": "/Users/z3437171/Library/R/arm64/4.6/library/rentrez",
        "reprex": "/Users/z3437171/Library/R/arm64/4.6/library/reprex",
        "reshape2": "/Users/z3437171/Library/R/arm64/4.6/library/reshape2",
        "rex": "/Users/z3437171/Library/R/arm64/4.6/library/rex",
        "rgbif": "/Users/z3437171/Library/R/arm64/4.6/library/rgbif",
        "rgpdd": "/Users/z3437171/Library/R/arm64/4.6/library/rgpdd",
        "rhub": "/Users/z3437171/Library/R/arm64/4.6/library/rhub",
        "rjags": "/Users/z3437171/Library/R/arm64/4.6/library/rjags",
        "rjson": "/Users/z3437171/Library/R/arm64/4.6/library/rjson",
        "rlang": "/Users/z3437171/Library/R/arm64/4.6/library/rlang",
        "rmarkdown": "/Users/z3437171/Library/R/arm64/4.6/library/rmarkdown",
        "rncl": "/Users/z3437171/Library/R/arm64/4.6/library/rncl",
        "rngtools": "/Users/z3437171/Library/R/arm64/4.6/library/rngtools",
        "rotl": "/Users/z3437171/Library/R/arm64/4.6/library/rotl",
        "roxygen2": "/Users/z3437171/Library/R/arm64/4.6/library/roxygen2",
        "rprojroot": "/Users/z3437171/Library/R/arm64/4.6/library/rprojroot",
        "rstan": "/Users/z3437171/Library/R/arm64/4.6/library/rstan",
        "rstantools": "/Users/z3437171/Library/R/arm64/4.6/library/rstantools",
        "rstudioapi": "/Users/z3437171/Library/R/arm64/4.6/library/rstudioapi",
        "rtdists": "/Users/z3437171/Library/R/arm64/4.6/library/rtdists",
        "rversions": "/Users/z3437171/Library/R/arm64/4.6/library/rversions",
        "rvest": "/Users/z3437171/Library/R/arm64/4.6/library/rvest",
        "s2": "/Users/z3437171/Library/R/arm64/4.6/library/s2",
        "safetensors": "/Users/z3437171/Library/R/arm64/4.6/library/safetensors",
        "sandwich": "/Users/z3437171/Library/R/arm64/4.6/library/sandwich",
        "sass": "/Users/z3437171/Library/R/arm64/4.6/library/sass",
        "satellite": "/Users/z3437171/Library/R/arm64/4.6/library/satellite",
        "scales": "/Users/z3437171/Library/R/arm64/4.6/library/scales",
        "scatterplot3d": "/Users/z3437171/Library/R/arm64/4.6/library/scatterplot3d",
        "selectr": "/Users/z3437171/Library/R/arm64/4.6/library/selectr",
        "servr": "/Users/z3437171/Library/R/arm64/4.6/library/servr",
        "sessioninfo": "/Users/z3437171/Library/R/arm64/4.6/library/sessioninfo",
        "sf": "/Users/z3437171/Library/R/arm64/4.6/library/sf",
        "sfheaders": "/Users/z3437171/Library/R/arm64/4.6/library/sfheaders",
        "shades": "/Users/z3437171/Library/R/arm64/4.6/library/shades",
        "shape": "/Users/z3437171/Library/R/arm64/4.6/library/shape",
        "shapviz": "/Users/z3437171/Library/R/arm64/4.6/library/shapviz",
        "shiny": "/Users/z3437171/Library/R/arm64/4.6/library/shiny",
        "shinyjs": "/Users/z3437171/Library/R/arm64/4.6/library/shinyjs",
        "shinystan": "/Users/z3437171/Library/R/arm64/4.6/library/shinystan",
        "shinythemes": "/Users/z3437171/Library/R/arm64/4.6/library/shinythemes",
        "skmeans": "/Users/z3437171/Library/R/arm64/4.6/library/skmeans",
        "slam": "/Users/z3437171/Library/R/arm64/4.6/library/slam",
        "smcfcs": "/Users/z3437171/Library/R/arm64/4.6/library/smcfcs",
        "sn": "/Users/z3437171/Library/R/arm64/4.6/library/sn",
        "snakecase": "/Users/z3437171/Library/R/arm64/4.6/library/snakecase",
        "sommer": "/Users/z3437171/Library/R/arm64/4.6/library/sommer",
        "sourcetools": "/Users/z3437171/Library/R/arm64/4.6/library/sourcetools",
        "sp": "/Users/z3437171/Library/R/arm64/4.6/library/sp",
        "spData": "/Users/z3437171/Library/R/arm64/4.6/library/spData",
        "spam": "/Users/z3437171/Library/R/arm64/4.6/library/spam",
        "spdep": "/Users/z3437171/Library/R/arm64/4.6/library/spdep",
        "spelling": "/Users/z3437171/Library/R/arm64/4.6/library/spelling",
        "splancs": "/Users/z3437171/Library/R/arm64/4.6/library/splancs",
        "splines2": "/Users/z3437171/Library/R/arm64/4.6/library/splines2",
        "statmod": "/Users/z3437171/Library/R/arm64/4.6/library/statmod",
        "stringfish": "/Users/z3437171/Library/R/arm64/4.6/library/stringfish",
        "stringi": "/Users/z3437171/Library/R/arm64/4.6/library/stringi",
        "stringr": "/Users/z3437171/Library/R/arm64/4.6/library/stringr",
        "styler": "/Users/z3437171/Library/R/arm64/4.6/library/styler",
        "subplex": "/Users/z3437171/Library/R/arm64/4.6/library/subplex",
        "svglite": "/Users/z3437171/Library/R/arm64/4.6/library/svglite",
        "sys": "/Users/z3437171/Library/R/arm64/4.6/library/sys",
        "systemfonts": "/Users/z3437171/Library/R/arm64/4.6/library/systemfonts",
        "tensorA": "/Users/z3437171/Library/R/arm64/4.6/library/tensorA",
        "terra": "/Users/z3437171/Library/R/arm64/4.6/library/terra",
        "tester": "/Users/z3437171/Library/R/arm64/4.6/library/tester",
        "testthat": "/Users/z3437171/Library/R/arm64/4.6/library/testthat",
        "textshaping": "/Users/z3437171/Library/R/arm64/4.6/library/textshaping",
        "threejs": "/Users/z3437171/Library/R/arm64/4.6/library/threejs",
        "tibble": "/Users/z3437171/Library/R/arm64/4.6/library/tibble",
        "tidygraph": "/Users/z3437171/Library/R/arm64/4.6/library/tidygraph",
        "tidyr": "/Users/z3437171/Library/R/arm64/4.6/library/tidyr",
        "tidyselect": "/Users/z3437171/Library/R/arm64/4.6/library/tidyselect",
        "tidyverse": "/Users/z3437171/Library/R/arm64/4.6/library/tidyverse",
        "timeDate": "/Users/z3437171/Library/R/arm64/4.6/library/timeDate",
        "timechange": "/Users/z3437171/Library/R/arm64/4.6/library/timechange",
        "tinytex": "/Users/z3437171/Library/R/arm64/4.6/library/tinytex",
        "torch": "/Users/z3437171/Library/R/arm64/4.6/library/torch",
        "triebeard": "/Users/z3437171/Library/R/arm64/4.6/library/triebeard",
        "truncnorm": "/Users/z3437171/Library/R/arm64/4.6/library/truncnorm",
        "tweedie": "/Users/z3437171/Library/R/arm64/4.6/library/tweedie",
        "tweenr": "/Users/z3437171/Library/R/arm64/4.6/library/tweenr",
        "tzdb": "/Users/z3437171/Library/R/arm64/4.6/library/tzdb",
        "ucminf": "/Users/z3437171/Library/R/arm64/4.6/library/ucminf",
        "units": "/Users/z3437171/Library/R/arm64/4.6/library/units",
        "urca": "/Users/z3437171/Library/R/arm64/4.6/library/urca",
        "urlchecker": "/Users/z3437171/Library/R/arm64/4.6/library/urlchecker",
        "urltools": "/Users/z3437171/Library/R/arm64/4.6/library/urltools",
        "usethis": "/Users/z3437171/Library/R/arm64/4.6/library/usethis",
        "utf8": "/Users/z3437171/Library/R/arm64/4.6/library/utf8",
        "uuid": "/Users/z3437171/Library/R/arm64/4.6/library/uuid",
        "vctrs": "/Users/z3437171/Library/R/arm64/4.6/library/vctrs",
        "vdiffr": "/Users/z3437171/Library/R/arm64/4.6/library/vdiffr",
        "vegan": "/Users/z3437171/Library/R/arm64/4.6/library/vegan",
        "vipor": "/Users/z3437171/Library/R/arm64/4.6/library/vipor",
        "viridis": "/Users/z3437171/Library/R/arm64/4.6/library/viridis",
        "viridisLite": "/Users/z3437171/Library/R/arm64/4.6/library/viridisLite",
        "visNetwork": "/Users/z3437171/Library/R/arm64/4.6/library/visNetwork",
        "vroom": "/Users/z3437171/Library/R/arm64/4.6/library/vroom",
        "waldo": "/Users/z3437171/Library/R/arm64/4.6/library/waldo",
        "whisker": "/Users/z3437171/Library/R/arm64/4.6/library/whisker",
        "whoami": "/Users/z3437171/Library/R/arm64/4.6/library/whoami",
        "withr": "/Users/z3437171/Library/R/arm64/4.6/library/withr",
        "wk": "/Users/z3437171/Library/R/arm64/4.6/library/wk",
        "xfun": "/Users/z3437171/Library/R/arm64/4.6/library/xfun",
        "xgboost": "/Users/z3437171/Library/R/arm64/4.6/library/xgboost",
        "xml2": "/Users/z3437171/Library/R/arm64/4.6/library/xml2",
        "xmlparsedata": "/Users/z3437171/Library/R/arm64/4.6/library/xmlparsedata",
        "xopen": "/Users/z3437171/Library/R/arm64/4.6/library/xopen",
        "xtable": "/Users/z3437171/Library/R/arm64/4.6/library/xtable",
        "xts": "/Users/z3437171/Library/R/arm64/4.6/library/xts",
        "yaml": "/Users/z3437171/Library/R/arm64/4.6/library/yaml",
        "zigg": "/Users/z3437171/Library/R/arm64/4.6/library/zigg",
        "zip": "/Users/z3437171/Library/R/arm64/4.6/library/zip",
        "zoo": "/Users/z3437171/Library/R/arm64/4.6/library/zoo"
      },
      "steps": [
        {
          "step": "install",
          "command": [
            "R",
            "CMD",
            "INSTALL",
            "--no-multiarch",
            "-l",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/neither",
            "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T16:59:51.212746+00:00",
          "elapsed_seconds": 21.82,
          "log": "neither-install.log",
          "log_sha256": "7f8a73b86b7d03565ff3003762e097a22aabf7fce7f9737a57df3161adc95993",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "fixtures",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "fixtures",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/neither"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:00:13.044612+00:00",
          "elapsed_seconds": 3.67,
          "log": "neither-fixtures.log",
          "log_sha256": "cb398da193bd8c6eb8c8ecf46cdf7b115e328991a7f5aad39c66efcc529ecb4e",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "real",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "neither",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/neither"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:00:16.721739+00:00",
          "elapsed_seconds": 2.595,
          "log": "neither-real.log",
          "log_sha256": "334e8ca94f3b17898c1b414ca275907ced1e7b85330b272464073a5646634615",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        }
      ]
    },
    {
      "cell": "drm-only",
      "mode": "drm-only",
      "library": "/private/tmp/pigauto-mi-matrix-c2e768ad/drm-only",
      "expected_optional_packages": [
        "drmTMB"
      ],
      "dependency_symlinks": {
        "AGHmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/AGHmatrix",
        "BACE": "/Users/z3437171/Library/R/arm64/4.6/library/BACE",
        "BH": "/Users/z3437171/Library/R/arm64/4.6/library/BH",
        "BayesLogit": "/Users/z3437171/Library/R/arm64/4.6/library/BayesLogit",
        "BiocGenerics": "/Users/z3437171/Library/R/arm64/4.6/library/BiocGenerics",
        "BiocManager": "/Users/z3437171/Library/R/arm64/4.6/library/BiocManager",
        "BiocVersion": "/Users/z3437171/Library/R/arm64/4.6/library/BiocVersion",
        "Brobdingnag": "/Users/z3437171/Library/R/arm64/4.6/library/Brobdingnag",
        "DBI": "/Users/z3437171/Library/R/arm64/4.6/library/DBI",
        "DEoptim": "/Users/z3437171/Library/R/arm64/4.6/library/DEoptim",
        "DHARMa": "/Users/z3437171/Library/R/arm64/4.6/library/DHARMa",
        "DT": "/Users/z3437171/Library/R/arm64/4.6/library/DT",
        "Deriv": "/Users/z3437171/Library/R/arm64/4.6/library/Deriv",
        "DiagrammeR": "/Users/z3437171/Library/R/arm64/4.6/library/DiagrammeR",
        "FNN": "/Users/z3437171/Library/R/arm64/4.6/library/FNN",
        "Formula": "/Users/z3437171/Library/R/arm64/4.6/library/Formula",
        "GPArotation": "/Users/z3437171/Library/R/arm64/4.6/library/GPArotation",
        "HSAUR3": "/Users/z3437171/Library/R/arm64/4.6/library/HSAUR3",
        "Hmsc": "/Users/z3437171/Library/R/arm64/4.6/library/Hmsc",
        "JuliaCall": "/Users/z3437171/Library/R/arm64/4.6/library/JuliaCall",
        "MCMCglmm": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCglmm",
        "MCMCpack": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCpack",
        "MatrixModels": "/Users/z3437171/Library/R/arm64/4.6/library/MatrixModels",
        "MuMIn": "/Users/z3437171/Library/R/arm64/4.6/library/MuMIn",
        "QGglmm": "/Users/z3437171/Library/R/arm64/4.6/library/QGglmm",
        "QuickJSR": "/Users/z3437171/Library/R/arm64/4.6/library/QuickJSR",
        "R.cache": "/Users/z3437171/Library/R/arm64/4.6/library/R.cache",
        "R.methodsS3": "/Users/z3437171/Library/R/arm64/4.6/library/R.methodsS3",
        "R.oo": "/Users/z3437171/Library/R/arm64/4.6/library/R.oo",
        "R.rsp": "/Users/z3437171/Library/R/arm64/4.6/library/R.rsp",
        "R.utils": "/Users/z3437171/Library/R/arm64/4.6/library/R.utils",
        "R6": "/Users/z3437171/Library/R/arm64/4.6/library/R6",
        "RColorBrewer": "/Users/z3437171/Library/R/arm64/4.6/library/RColorBrewer",
        "ROI": "/Users/z3437171/Library/R/arm64/4.6/library/ROI",
        "ROI.plugin.lpsolve": "/Users/z3437171/Library/R/arm64/4.6/library/ROI.plugin.lpsolve",
        "RSpectra": "/Users/z3437171/Library/R/arm64/4.6/library/RSpectra",
        "RTMB": "/Users/z3437171/Library/R/arm64/4.6/library/RTMB",
        "RTMBdist": "/Users/z3437171/Library/R/arm64/4.6/library/RTMBdist",
        "RWiener": "/Users/z3437171/Library/R/arm64/4.6/library/RWiener",
        "Rcpp": "/Users/z3437171/Library/R/arm64/4.6/library/Rcpp",
        "RcppArmadillo": "/Users/z3437171/Library/R/arm64/4.6/library/RcppArmadillo",
        "RcppEigen": "/Users/z3437171/Library/R/arm64/4.6/library/RcppEigen",
        "RcppParallel": "/Users/z3437171/Library/R/arm64/4.6/library/RcppParallel",
        "RcppProgress": "/Users/z3437171/Library/R/arm64/4.6/library/RcppProgress",
        "Rdpack": "/Users/z3437171/Library/R/arm64/4.6/library/Rdpack",
        "Rfast": "/Users/z3437171/Library/R/arm64/4.6/library/Rfast",
        "Rmpfr": "/Users/z3437171/Library/R/arm64/4.6/library/Rmpfr",
        "Rphylopars": "/Users/z3437171/Library/R/arm64/4.6/library/Rphylopars",
        "S7": "/Users/z3437171/Library/R/arm64/4.6/library/S7",
        "SimDesign": "/Users/z3437171/Library/R/arm64/4.6/library/SimDesign",
        "SparseM": "/Users/z3437171/Library/R/arm64/4.6/library/SparseM",
        "StanHeaders": "/Users/z3437171/Library/R/arm64/4.6/library/StanHeaders",
        "TH.data": "/Users/z3437171/Library/R/arm64/4.6/library/TH.data",
        "TMB": "/Users/z3437171/Library/R/arm64/4.6/library/TMB",
        "V8": "/Users/z3437171/Library/R/arm64/4.6/library/V8",
        "VGAM": "/Users/z3437171/Library/R/arm64/4.6/library/VGAM",
        "XML": "/Users/z3437171/Library/R/arm64/4.6/library/XML",
        "_cache": "/Users/z3437171/Library/R/arm64/4.6/library/_cache",
        "abind": "/Users/z3437171/Library/R/arm64/4.6/library/abind",
        "ade4": "/Users/z3437171/Library/R/arm64/4.6/library/ade4",
        "alabama": "/Users/z3437171/Library/R/arm64/4.6/library/alabama",
        "ape": "/Users/z3437171/Library/R/arm64/4.6/library/ape",
        "arm": "/Users/z3437171/Library/R/arm64/4.6/library/arm",
        "askpass": "/Users/z3437171/Library/R/arm64/4.6/library/askpass",
        "assertthat": "/Users/z3437171/Library/R/arm64/4.6/library/assertthat",
        "audio": "/Users/z3437171/Library/R/arm64/4.6/library/audio",
        "backports": "/Users/z3437171/Library/R/arm64/4.6/library/backports",
        "base64enc": "/Users/z3437171/Library/R/arm64/4.6/library/base64enc",
        "bayeslincom": "/Users/z3437171/Library/R/arm64/4.6/library/bayeslincom",
        "bayesplot": "/Users/z3437171/Library/R/arm64/4.6/library/bayesplot",
        "bayestestR": "/Users/z3437171/Library/R/arm64/4.6/library/bayestestR",
        "bbmle": "/Users/z3437171/Library/R/arm64/4.6/library/bbmle",
        "bdsmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/bdsmatrix",
        "beepr": "/Users/z3437171/Library/R/arm64/4.6/library/beepr",
        "beeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/beeswarm",
        "betareg": "/Users/z3437171/Library/R/arm64/4.6/library/betareg",
        "bit": "/Users/z3437171/Library/R/arm64/4.6/library/bit",
        "bit64": "/Users/z3437171/Library/R/arm64/4.6/library/bit64",
        "blob": "/Users/z3437171/Library/R/arm64/4.6/library/blob",
        "blsmeta": "/Users/z3437171/Library/R/arm64/4.6/library/blsmeta",
        "brew": "/Users/z3437171/Library/R/arm64/4.6/library/brew",
        "brglm2": "/Users/z3437171/Library/R/arm64/4.6/library/brglm2",
        "bridgesampling": "/Users/z3437171/Library/R/arm64/4.6/library/bridgesampling",
        "brio": "/Users/z3437171/Library/R/arm64/4.6/library/brio",
        "brms": "/Users/z3437171/Library/R/arm64/4.6/library/brms",
        "broom": "/Users/z3437171/Library/R/arm64/4.6/library/broom",
        "broom.mixed": "/Users/z3437171/Library/R/arm64/4.6/library/broom.mixed",
        "bslib": "/Users/z3437171/Library/R/arm64/4.6/library/bslib",
        "cachem": "/Users/z3437171/Library/R/arm64/4.6/library/cachem",
        "callr": "/Users/z3437171/Library/R/arm64/4.6/library/callr",
        "car": "/Users/z3437171/Library/R/arm64/4.6/library/car",
        "carData": "/Users/z3437171/Library/R/arm64/4.6/library/carData",
        "castor": "/Users/z3437171/Library/R/arm64/4.6/library/castor",
        "cellranger": "/Users/z3437171/Library/R/arm64/4.6/library/cellranger",
        "checkmate": "/Users/z3437171/Library/R/arm64/4.6/library/checkmate",
        "circular": "/Users/z3437171/Library/R/arm64/4.6/library/circular",
        "classInt": "/Users/z3437171/Library/R/arm64/4.6/library/classInt",
        "cli": "/Users/z3437171/Library/R/arm64/4.6/library/cli",
        "clipr": "/Users/z3437171/Library/R/arm64/4.6/library/clipr",
        "clue": "/Users/z3437171/Library/R/arm64/4.6/library/clue",
        "clusterGeneration": "/Users/z3437171/Library/R/arm64/4.6/library/clusterGeneration",
        "cmdstanr": "/Users/z3437171/Library/R/arm64/4.6/library/cmdstanr",
        "coda": "/Users/z3437171/Library/R/arm64/4.6/library/coda",
        "collections": "/Users/z3437171/Library/R/arm64/4.6/library/collections",
        "colorspace": "/Users/z3437171/Library/R/arm64/4.6/library/colorspace",
        "colourpicker": "/Users/z3437171/Library/R/arm64/4.6/library/colourpicker",
        "combinat": "/Users/z3437171/Library/R/arm64/4.6/library/combinat",
        "commonmark": "/Users/z3437171/Library/R/arm64/4.6/library/commonmark",
        "conflicted": "/Users/z3437171/Library/R/arm64/4.6/library/conflicted",
        "coro": "/Users/z3437171/Library/R/arm64/4.6/library/coro",
        "corpcor": "/Users/z3437171/Library/R/arm64/4.6/library/corpcor",
        "covr": "/Users/z3437171/Library/R/arm64/4.6/library/covr",
        "cowplot": "/Users/z3437171/Library/R/arm64/4.6/library/cowplot",
        "cpp11": "/Users/z3437171/Library/R/arm64/4.6/library/cpp11",
        "crayon": "/Users/z3437171/Library/R/arm64/4.6/library/crayon",
        "credentials": "/Users/z3437171/Library/R/arm64/4.6/library/credentials",
        "crosstalk": "/Users/z3437171/Library/R/arm64/4.6/library/crosstalk",
        "crul": "/Users/z3437171/Library/R/arm64/4.6/library/crul",
        "cubature": "/Users/z3437171/Library/R/arm64/4.6/library/cubature",
        "curl": "/Users/z3437171/Library/R/arm64/4.6/library/curl",
        "dagitty": "/Users/z3437171/Library/R/arm64/4.6/library/dagitty",
        "data.table": "/Users/z3437171/Library/R/arm64/4.6/library/data.table",
        "datawizard": "/Users/z3437171/Library/R/arm64/4.6/library/datawizard",
        "dbplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dbplyr",
        "dcurver": "/Users/z3437171/Library/R/arm64/4.6/library/dcurver",
        "deldir": "/Users/z3437171/Library/R/arm64/4.6/library/deldir",
        "desc": "/Users/z3437171/Library/R/arm64/4.6/library/desc",
        "detectseparation": "/Users/z3437171/Library/R/arm64/4.6/library/detectseparation",
        "devtools": "/Users/z3437171/Library/R/arm64/4.6/library/devtools",
        "diffobj": "/Users/z3437171/Library/R/arm64/4.6/library/diffobj",
        "digest": "/Users/z3437171/Library/R/arm64/4.6/library/digest",
        "distributional": "/Users/z3437171/Library/R/arm64/4.6/library/distributional",
        "doBy": "/Users/z3437171/Library/R/arm64/4.6/library/doBy",
        "doFuture": "/Users/z3437171/Library/R/arm64/4.6/library/doFuture",
        "doParallel": "/Users/z3437171/Library/R/arm64/4.6/library/doParallel",
        "doRNG": "/Users/z3437171/Library/R/arm64/4.6/library/doRNG",
        "dotCall64": "/Users/z3437171/Library/R/arm64/4.6/library/dotCall64",
        "downlit": "/Users/z3437171/Library/R/arm64/4.6/library/downlit",
        "dplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dplyr",
        "drmSEM": "/Users/z3437171/Library/R/arm64/4.6/library/drmSEM",
        "drmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/drmTMB",
        "dtplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dtplyr",
        "dygraphs": "/Users/z3437171/Library/R/arm64/4.6/library/dygraphs",
        "e1071": "/Users/z3437171/Library/R/arm64/4.6/library/e1071",
        "ellipsis": "/Users/z3437171/Library/R/arm64/4.6/library/ellipsis",
        "emmeans": "/Users/z3437171/Library/R/arm64/4.6/library/emmeans",
        "enhancer": "/Users/z3437171/Library/R/arm64/4.6/library/enhancer",
        "enrichwith": "/Users/z3437171/Library/R/arm64/4.6/library/enrichwith",
        "estimability": "/Users/z3437171/Library/R/arm64/4.6/library/estimability",
        "evaluate": "/Users/z3437171/Library/R/arm64/4.6/library/evaluate",
        "evd": "/Users/z3437171/Library/R/arm64/4.6/library/evd",
        "expm": "/Users/z3437171/Library/R/arm64/4.6/library/expm",
        "extraDistr": "/Users/z3437171/Library/R/arm64/4.6/library/extraDistr",
        "fansi": "/Users/z3437171/Library/R/arm64/4.6/library/fansi",
        "farver": "/Users/z3437171/Library/R/arm64/4.6/library/farver",
        "fastGHQuad": "/Users/z3437171/Library/R/arm64/4.6/library/fastGHQuad",
        "fastmap": "/Users/z3437171/Library/R/arm64/4.6/library/fastmap",
        "fastmatch": "/Users/z3437171/Library/R/arm64/4.6/library/fastmatch",
        "fields": "/Users/z3437171/Library/R/arm64/4.6/library/fields",
        "fishMod": "/Users/z3437171/Library/R/arm64/4.6/library/fishMod",
        "flexmix": "/Users/z3437171/Library/R/arm64/4.6/library/flexmix",
        "fmesher": "/Users/z3437171/Library/R/arm64/4.6/library/fmesher",
        "fontawesome": "/Users/z3437171/Library/R/arm64/4.6/library/fontawesome",
        "forcats": "/Users/z3437171/Library/R/arm64/4.6/library/forcats",
        "foreach": "/Users/z3437171/Library/R/arm64/4.6/library/foreach",
        "forecast": "/Users/z3437171/Library/R/arm64/4.6/library/forecast",
        "fracdiff": "/Users/z3437171/Library/R/arm64/4.6/library/fracdiff",
        "freqTLS": "/Users/z3437171/Library/R/arm64/4.6/library/freqTLS",
        "fs": "/Users/z3437171/Library/R/arm64/4.6/library/fs",
        "furrr": "/Users/z3437171/Library/R/arm64/4.6/library/furrr",
        "future": "/Users/z3437171/Library/R/arm64/4.6/library/future",
        "future.apply": "/Users/z3437171/Library/R/arm64/4.6/library/future.apply",
        "future.mirai": "/Users/z3437171/Library/R/arm64/4.6/library/future.mirai",
        "galamm": "/Users/z3437171/Library/R/arm64/4.6/library/galamm",
        "gamlss.dist": "/Users/z3437171/Library/R/arm64/4.6/library/gamlss.dist",
        "gamm4": "/Users/z3437171/Library/R/arm64/4.6/library/gamm4",
        "gap": "/Users/z3437171/Library/R/arm64/4.6/library/gap",
        "gap.datasets": "/Users/z3437171/Library/R/arm64/4.6/library/gap.datasets",
        "gargle": "/Users/z3437171/Library/R/arm64/4.6/library/gargle",
        "generics": "/Users/z3437171/Library/R/arm64/4.6/library/generics",
        "geojsonsf": "/Users/z3437171/Library/R/arm64/4.6/library/geojsonsf",
        "geometries": "/Users/z3437171/Library/R/arm64/4.6/library/geometries",
        "gert": "/Users/z3437171/Library/R/arm64/4.6/library/gert",
        "ggbeeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/ggbeeswarm",
        "ggdist": "/Users/z3437171/Library/R/arm64/4.6/library/ggdist",
        "ggfittext": "/Users/z3437171/Library/R/arm64/4.6/library/ggfittext",
        "ggforce": "/Users/z3437171/Library/R/arm64/4.6/library/ggforce",
        "gggenes": "/Users/z3437171/Library/R/arm64/4.6/library/gggenes",
        "ggh4x": "/Users/z3437171/Library/R/arm64/4.6/library/ggh4x",
        "ggm": "/Users/z3437171/Library/R/arm64/4.6/library/ggm",
        "ggokabeito": "/Users/z3437171/Library/R/arm64/4.6/library/ggokabeito",
        "ggplot2": "/Users/z3437171/Library/R/arm64/4.6/library/ggplot2",
        "ggraph": "/Users/z3437171/Library/R/arm64/4.6/library/ggraph",
        "ggrepel": "/Users/z3437171/Library/R/arm64/4.6/library/ggrepel",
        "ggridges": "/Users/z3437171/Library/R/arm64/4.6/library/ggridges",
        "gh": "/Users/z3437171/Library/R/arm64/4.6/library/gh",
        "gitcreds": "/Users/z3437171/Library/R/arm64/4.6/library/gitcreds",
        "glassoFast": "/Users/z3437171/Library/R/arm64/4.6/library/glassoFast",
        "gllvm": "/Users/z3437171/Library/R/arm64/4.6/library/gllvm",
        "glmmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/glmmTMB",
        "glmnet": "/Users/z3437171/Library/R/arm64/4.6/library/glmnet",
        "globals": "/Users/z3437171/Library/R/arm64/4.6/library/globals",
        "glue": "/Users/z3437171/Library/R/arm64/4.6/library/glue",
        "gmp": "/Users/z3437171/Library/R/arm64/4.6/library/gmp",
        "googledrive": "/Users/z3437171/Library/R/arm64/4.6/library/googledrive",
        "googlesheets4": "/Users/z3437171/Library/R/arm64/4.6/library/googlesheets4",
        "graph": "/Users/z3437171/Library/R/arm64/4.6/library/graph",
        "graphlayouts": "/Users/z3437171/Library/R/arm64/4.6/library/graphlayouts",
        "gratia": "/Users/z3437171/Library/R/arm64/4.6/library/gratia",
        "gridExtra": "/Users/z3437171/Library/R/arm64/4.6/library/gridExtra",
        "gridtext": "/Users/z3437171/Library/R/arm64/4.6/library/gridtext",
        "gsl": "/Users/z3437171/Library/R/arm64/4.6/library/gsl",
        "gtable": "/Users/z3437171/Library/R/arm64/4.6/library/gtable",
        "gtools": "/Users/z3437171/Library/R/arm64/4.6/library/gtools",
        "haven": "/Users/z3437171/Library/R/arm64/4.6/library/haven",
        "here": "/Users/z3437171/Library/R/arm64/4.6/library/here",
        "highr": "/Users/z3437171/Library/R/arm64/4.6/library/highr",
        "hms": "/Users/z3437171/Library/R/arm64/4.6/library/hms",
        "hsquared": "/Users/z3437171/Library/R/arm64/4.6/library/hsquared",
        "htmltools": "/Users/z3437171/Library/R/arm64/4.6/library/htmltools",
        "htmlwidgets": "/Users/z3437171/Library/R/arm64/4.6/library/htmlwidgets",
        "httpcode": "/Users/z3437171/Library/R/arm64/4.6/library/httpcode",
        "httpuv": "/Users/z3437171/Library/R/arm64/4.6/library/httpuv",
        "httr": "/Users/z3437171/Library/R/arm64/4.6/library/httr",
        "httr2": "/Users/z3437171/Library/R/arm64/4.6/library/httr2",
        "hunspell": "/Users/z3437171/Library/R/arm64/4.6/library/hunspell",
        "ids": "/Users/z3437171/Library/R/arm64/4.6/library/ids",
        "igraph": "/Users/z3437171/Library/R/arm64/4.6/library/igraph",
        "ini": "/Users/z3437171/Library/R/arm64/4.6/library/ini",
        "inline": "/Users/z3437171/Library/R/arm64/4.6/library/inline",
        "insight": "/Users/z3437171/Library/R/arm64/4.6/library/insight",
        "isoband": "/Users/z3437171/Library/R/arm64/4.6/library/isoband",
        "iterators": "/Users/z3437171/Library/R/arm64/4.6/library/iterators",
        "itertools": "/Users/z3437171/Library/R/arm64/4.6/library/itertools",
        "janitor": "/Users/z3437171/Library/R/arm64/4.6/library/janitor",
        "jomo": "/Users/z3437171/Library/R/arm64/4.6/library/jomo",
        "jpeg": "/Users/z3437171/Library/R/arm64/4.6/library/jpeg",
        "jquerylib": "/Users/z3437171/Library/R/arm64/4.6/library/jquerylib",
        "jsonify": "/Users/z3437171/Library/R/arm64/4.6/library/jsonify",
        "jsonlite": "/Users/z3437171/Library/R/arm64/4.6/library/jsonlite",
        "kernelshap": "/Users/z3437171/Library/R/arm64/4.6/library/kernelshap",
        "knitr": "/Users/z3437171/Library/R/arm64/4.6/library/knitr",
        "labeling": "/Users/z3437171/Library/R/arm64/4.6/library/labeling",
        "languageserver": "/Users/z3437171/Library/R/arm64/4.6/library/languageserver",
        "later": "/Users/z3437171/Library/R/arm64/4.6/library/later",
        "latex2exp": "/Users/z3437171/Library/R/arm64/4.6/library/latex2exp",
        "lavaan": "/Users/z3437171/Library/R/arm64/4.6/library/lavaan",
        "lazyeval": "/Users/z3437171/Library/R/arm64/4.6/library/lazyeval",
        "leafem": "/Users/z3437171/Library/R/arm64/4.6/library/leafem",
        "leaflet": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet",
        "leaflet.providers": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet.providers",
        "leafpop": "/Users/z3437171/Library/R/arm64/4.6/library/leafpop",
        "lifecycle": "/Users/z3437171/Library/R/arm64/4.6/library/lifecycle",
        "lintr": "/Users/z3437171/Library/R/arm64/4.6/library/lintr",
        "listenv": "/Users/z3437171/Library/R/arm64/4.6/library/listenv",
        "litedown": "/Users/z3437171/Library/R/arm64/4.6/library/litedown",
        "lme4": "/Users/z3437171/Library/R/arm64/4.6/library/lme4",
        "lmtest": "/Users/z3437171/Library/R/arm64/4.6/library/lmtest",
        "loo": "/Users/z3437171/Library/R/arm64/4.6/library/loo",
        "lpSolveAPI": "/Users/z3437171/Library/R/arm64/4.6/library/lpSolveAPI",
        "lubridate": "/Users/z3437171/Library/R/arm64/4.6/library/lubridate",
        "magrittr": "/Users/z3437171/Library/R/arm64/4.6/library/magrittr",
        "maps": "/Users/z3437171/Library/R/arm64/4.6/library/maps",
        "mapview": "/Users/z3437171/Library/R/arm64/4.6/library/mapview",
        "markdown": "/Users/z3437171/Library/R/arm64/4.6/library/markdown",
        "mathjaxr": "/Users/z3437171/Library/R/arm64/4.6/library/mathjaxr",
        "matrixStats": "/Users/z3437171/Library/R/arm64/4.6/library/matrixStats",
        "mclogit": "/Users/z3437171/Library/R/arm64/4.6/library/mclogit",
        "mcmc": "/Users/z3437171/Library/R/arm64/4.6/library/mcmc",
        "memisc": "/Users/z3437171/Library/R/arm64/4.6/library/memisc",
        "memoise": "/Users/z3437171/Library/R/arm64/4.6/library/memoise",
        "metadat": "/Users/z3437171/Library/R/arm64/4.6/library/metadat",
        "metafor": "/Users/z3437171/Library/R/arm64/4.6/library/metafor",
        "metaplus": "/Users/z3437171/Library/R/arm64/4.6/library/metaplus",
        "mice": "/Users/z3437171/Library/R/arm64/4.6/library/mice",
        "microbenchmark": "/Users/z3437171/Library/R/arm64/4.6/library/microbenchmark",
        "mime": "/Users/z3437171/Library/R/arm64/4.6/library/mime",
        "miniUI": "/Users/z3437171/Library/R/arm64/4.6/library/miniUI",
        "minqa": "/Users/z3437171/Library/R/arm64/4.6/library/minqa",
        "mirai": "/Users/z3437171/Library/R/arm64/4.6/library/mirai",
        "mirt": "/Users/z3437171/Library/R/arm64/4.6/library/mirt",
        "missForest": "/Users/z3437171/Library/R/arm64/4.6/library/missForest",
        "mitml": "/Users/z3437171/Library/R/arm64/4.6/library/mitml",
        "mnormt": "/Users/z3437171/Library/R/arm64/4.6/library/mnormt",
        "modelr": "/Users/z3437171/Library/R/arm64/4.6/library/modelr",
        "modeltools": "/Users/z3437171/Library/R/arm64/4.6/library/modeltools",
        "movMF": "/Users/z3437171/Library/R/arm64/4.6/library/movMF",
        "msm": "/Users/z3437171/Library/R/arm64/4.6/library/msm",
        "multcomp": "/Users/z3437171/Library/R/arm64/4.6/library/multcomp",
        "mvMORPH": "/Users/z3437171/Library/R/arm64/4.6/library/mvMORPH",
        "mvnfast": "/Users/z3437171/Library/R/arm64/4.6/library/mvnfast",
        "mvtnorm": "/Users/z3437171/Library/R/arm64/4.6/library/mvtnorm",
        "nadiv": "/Users/z3437171/Library/R/arm64/4.6/library/nadiv",
        "nanonext": "/Users/z3437171/Library/R/arm64/4.6/library/nanonext",
        "naturalsort": "/Users/z3437171/Library/R/arm64/4.6/library/naturalsort",
        "ncf": "/Users/z3437171/Library/R/arm64/4.6/library/ncf",
        "nleqslv": "/Users/z3437171/Library/R/arm64/4.6/library/nleqslv",
        "nloptr": "/Users/z3437171/Library/R/arm64/4.6/library/nloptr",
        "numDeriv": "/Users/z3437171/Library/R/arm64/4.6/library/numDeriv",
        "oai": "/Users/z3437171/Library/R/arm64/4.6/library/oai",
        "openssl": "/Users/z3437171/Library/R/arm64/4.6/library/openssl",
        "optimParallel": "/Users/z3437171/Library/R/arm64/4.6/library/optimParallel",
        "orchaRd": "/Users/z3437171/Library/R/arm64/4.6/library/orchaRd",
        "ordinal": "/Users/z3437171/Library/R/arm64/4.6/library/ordinal",
        "orthopolynom": "/Users/z3437171/Library/R/arm64/4.6/library/orthopolynom",
        "otel": "/Users/z3437171/Library/R/arm64/4.6/library/otel",
        "pROC": "/Users/z3437171/Library/R/arm64/4.6/library/pROC",
        "pak": "/Users/z3437171/Library/R/arm64/4.6/library/pak",
        "palmerpenguins": "/Users/z3437171/Library/R/arm64/4.6/library/palmerpenguins",
        "pan": "/Users/z3437171/Library/R/arm64/4.6/library/pan",
        "parallelly": "/Users/z3437171/Library/R/arm64/4.6/library/parallelly",
        "patchwork": "/Users/z3437171/Library/R/arm64/4.6/library/patchwork",
        "pbapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbapply",
        "pbivnorm": "/Users/z3437171/Library/R/arm64/4.6/library/pbivnorm",
        "pbkrtest": "/Users/z3437171/Library/R/arm64/4.6/library/pbkrtest",
        "pbmcapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbmcapply",
        "pedigreemm": "/Users/z3437171/Library/R/arm64/4.6/library/pedigreemm",
        "performance": "/Users/z3437171/Library/R/arm64/4.6/library/performance",
        "permute": "/Users/z3437171/Library/R/arm64/4.6/library/permute",
        "phangorn": "/Users/z3437171/Library/R/arm64/4.6/library/phangorn",
        "phylolm": "/Users/z3437171/Library/R/arm64/4.6/library/phylolm",
        "phylopath": "/Users/z3437171/Library/R/arm64/4.6/library/phylopath",
        "phytools": "/Users/z3437171/Library/R/arm64/4.6/library/phytools",
        "piecewiseSEM": "/Users/z3437171/Library/R/arm64/4.6/library/piecewiseSEM",
        "pillar": "/Users/z3437171/Library/R/arm64/4.6/library/pillar",
        "pixmap": "/Users/z3437171/Library/R/arm64/4.6/library/pixmap",
        "pkgbuild": "/Users/z3437171/Library/R/arm64/4.6/library/pkgbuild",
        "pkgconfig": "/Users/z3437171/Library/R/arm64/4.6/library/pkgconfig",
        "pkgdown": "/Users/z3437171/Library/R/arm64/4.6/library/pkgdown",
        "pkgload": "/Users/z3437171/Library/R/arm64/4.6/library/pkgload",
        "plotly": "/Users/z3437171/Library/R/arm64/4.6/library/plotly",
        "plyr": "/Users/z3437171/Library/R/arm64/4.6/library/plyr",
        "png": "/Users/z3437171/Library/R/arm64/4.6/library/png",
        "polyclip": "/Users/z3437171/Library/R/arm64/4.6/library/polyclip",
        "polynom": "/Users/z3437171/Library/R/arm64/4.6/library/polynom",
        "posterior": "/Users/z3437171/Library/R/arm64/4.6/library/posterior",
        "pracma": "/Users/z3437171/Library/R/arm64/4.6/library/pracma",
        "praise": "/Users/z3437171/Library/R/arm64/4.6/library/praise",
        "prettyunits": "/Users/z3437171/Library/R/arm64/4.6/library/prettyunits",
        "priorsense": "/Users/z3437171/Library/R/arm64/4.6/library/priorsense",
        "processx": "/Users/z3437171/Library/R/arm64/4.6/library/processx",
        "profvis": "/Users/z3437171/Library/R/arm64/4.6/library/profvis",
        "progress": "/Users/z3437171/Library/R/arm64/4.6/library/progress",
        "progressr": "/Users/z3437171/Library/R/arm64/4.6/library/progressr",
        "projpred": "/Users/z3437171/Library/R/arm64/4.6/library/projpred",
        "promises": "/Users/z3437171/Library/R/arm64/4.6/library/promises",
        "proxy": "/Users/z3437171/Library/R/arm64/4.6/library/proxy",
        "ps": "/Users/z3437171/Library/R/arm64/4.6/library/ps",
        "psymetadata": "/Users/z3437171/Library/R/arm64/4.6/library/psymetadata",
        "purrr": "/Users/z3437171/Library/R/arm64/4.6/library/purrr",
        "qgam": "/Users/z3437171/Library/R/arm64/4.6/library/qgam",
        "qs2": "/Users/z3437171/Library/R/arm64/4.6/library/qs2",
        "quadprog": "/Users/z3437171/Library/R/arm64/4.6/library/quadprog",
        "quantreg": "/Users/z3437171/Library/R/arm64/4.6/library/quantreg",
        "rWishart": "/Users/z3437171/Library/R/arm64/4.6/library/rWishart",
        "ragg": "/Users/z3437171/Library/R/arm64/4.6/library/ragg",
        "randomForest": "/Users/z3437171/Library/R/arm64/4.6/library/randomForest",
        "ranger": "/Users/z3437171/Library/R/arm64/4.6/library/ranger",
        "rapidjsonr": "/Users/z3437171/Library/R/arm64/4.6/library/rapidjsonr",
        "rappdirs": "/Users/z3437171/Library/R/arm64/4.6/library/rappdirs",
        "raster": "/Users/z3437171/Library/R/arm64/4.6/library/raster",
        "rbibutils": "/Users/z3437171/Library/R/arm64/4.6/library/rbibutils",
        "rcmdcheck": "/Users/z3437171/Library/R/arm64/4.6/library/rcmdcheck",
        "readr": "/Users/z3437171/Library/R/arm64/4.6/library/readr",
        "readxl": "/Users/z3437171/Library/R/arm64/4.6/library/readxl",
        "reformulas": "/Users/z3437171/Library/R/arm64/4.6/library/reformulas",
        "registry": "/Users/z3437171/Library/R/arm64/4.6/library/registry",
        "rematch": "/Users/z3437171/Library/R/arm64/4.6/library/rematch",
        "rematch2": "/Users/z3437171/Library/R/arm64/4.6/library/rematch2",
        "remotes": "/Users/z3437171/Library/R/arm64/4.6/library/remotes",
        "rentrez": "/Users/z3437171/Library/R/arm64/4.6/library/rentrez",
        "reprex": "/Users/z3437171/Library/R/arm64/4.6/library/reprex",
        "reshape2": "/Users/z3437171/Library/R/arm64/4.6/library/reshape2",
        "rex": "/Users/z3437171/Library/R/arm64/4.6/library/rex",
        "rgbif": "/Users/z3437171/Library/R/arm64/4.6/library/rgbif",
        "rgpdd": "/Users/z3437171/Library/R/arm64/4.6/library/rgpdd",
        "rhub": "/Users/z3437171/Library/R/arm64/4.6/library/rhub",
        "rjags": "/Users/z3437171/Library/R/arm64/4.6/library/rjags",
        "rjson": "/Users/z3437171/Library/R/arm64/4.6/library/rjson",
        "rlang": "/Users/z3437171/Library/R/arm64/4.6/library/rlang",
        "rmarkdown": "/Users/z3437171/Library/R/arm64/4.6/library/rmarkdown",
        "rncl": "/Users/z3437171/Library/R/arm64/4.6/library/rncl",
        "rngtools": "/Users/z3437171/Library/R/arm64/4.6/library/rngtools",
        "rotl": "/Users/z3437171/Library/R/arm64/4.6/library/rotl",
        "roxygen2": "/Users/z3437171/Library/R/arm64/4.6/library/roxygen2",
        "rprojroot": "/Users/z3437171/Library/R/arm64/4.6/library/rprojroot",
        "rstan": "/Users/z3437171/Library/R/arm64/4.6/library/rstan",
        "rstantools": "/Users/z3437171/Library/R/arm64/4.6/library/rstantools",
        "rstudioapi": "/Users/z3437171/Library/R/arm64/4.6/library/rstudioapi",
        "rtdists": "/Users/z3437171/Library/R/arm64/4.6/library/rtdists",
        "rversions": "/Users/z3437171/Library/R/arm64/4.6/library/rversions",
        "rvest": "/Users/z3437171/Library/R/arm64/4.6/library/rvest",
        "s2": "/Users/z3437171/Library/R/arm64/4.6/library/s2",
        "safetensors": "/Users/z3437171/Library/R/arm64/4.6/library/safetensors",
        "sandwich": "/Users/z3437171/Library/R/arm64/4.6/library/sandwich",
        "sass": "/Users/z3437171/Library/R/arm64/4.6/library/sass",
        "satellite": "/Users/z3437171/Library/R/arm64/4.6/library/satellite",
        "scales": "/Users/z3437171/Library/R/arm64/4.6/library/scales",
        "scatterplot3d": "/Users/z3437171/Library/R/arm64/4.6/library/scatterplot3d",
        "selectr": "/Users/z3437171/Library/R/arm64/4.6/library/selectr",
        "servr": "/Users/z3437171/Library/R/arm64/4.6/library/servr",
        "sessioninfo": "/Users/z3437171/Library/R/arm64/4.6/library/sessioninfo",
        "sf": "/Users/z3437171/Library/R/arm64/4.6/library/sf",
        "sfheaders": "/Users/z3437171/Library/R/arm64/4.6/library/sfheaders",
        "shades": "/Users/z3437171/Library/R/arm64/4.6/library/shades",
        "shape": "/Users/z3437171/Library/R/arm64/4.6/library/shape",
        "shapviz": "/Users/z3437171/Library/R/arm64/4.6/library/shapviz",
        "shiny": "/Users/z3437171/Library/R/arm64/4.6/library/shiny",
        "shinyjs": "/Users/z3437171/Library/R/arm64/4.6/library/shinyjs",
        "shinystan": "/Users/z3437171/Library/R/arm64/4.6/library/shinystan",
        "shinythemes": "/Users/z3437171/Library/R/arm64/4.6/library/shinythemes",
        "skmeans": "/Users/z3437171/Library/R/arm64/4.6/library/skmeans",
        "slam": "/Users/z3437171/Library/R/arm64/4.6/library/slam",
        "smcfcs": "/Users/z3437171/Library/R/arm64/4.6/library/smcfcs",
        "sn": "/Users/z3437171/Library/R/arm64/4.6/library/sn",
        "snakecase": "/Users/z3437171/Library/R/arm64/4.6/library/snakecase",
        "sommer": "/Users/z3437171/Library/R/arm64/4.6/library/sommer",
        "sourcetools": "/Users/z3437171/Library/R/arm64/4.6/library/sourcetools",
        "sp": "/Users/z3437171/Library/R/arm64/4.6/library/sp",
        "spData": "/Users/z3437171/Library/R/arm64/4.6/library/spData",
        "spam": "/Users/z3437171/Library/R/arm64/4.6/library/spam",
        "spdep": "/Users/z3437171/Library/R/arm64/4.6/library/spdep",
        "spelling": "/Users/z3437171/Library/R/arm64/4.6/library/spelling",
        "splancs": "/Users/z3437171/Library/R/arm64/4.6/library/splancs",
        "splines2": "/Users/z3437171/Library/R/arm64/4.6/library/splines2",
        "statmod": "/Users/z3437171/Library/R/arm64/4.6/library/statmod",
        "stringfish": "/Users/z3437171/Library/R/arm64/4.6/library/stringfish",
        "stringi": "/Users/z3437171/Library/R/arm64/4.6/library/stringi",
        "stringr": "/Users/z3437171/Library/R/arm64/4.6/library/stringr",
        "styler": "/Users/z3437171/Library/R/arm64/4.6/library/styler",
        "subplex": "/Users/z3437171/Library/R/arm64/4.6/library/subplex",
        "svglite": "/Users/z3437171/Library/R/arm64/4.6/library/svglite",
        "sys": "/Users/z3437171/Library/R/arm64/4.6/library/sys",
        "systemfonts": "/Users/z3437171/Library/R/arm64/4.6/library/systemfonts",
        "tensorA": "/Users/z3437171/Library/R/arm64/4.6/library/tensorA",
        "terra": "/Users/z3437171/Library/R/arm64/4.6/library/terra",
        "tester": "/Users/z3437171/Library/R/arm64/4.6/library/tester",
        "testthat": "/Users/z3437171/Library/R/arm64/4.6/library/testthat",
        "textshaping": "/Users/z3437171/Library/R/arm64/4.6/library/textshaping",
        "threejs": "/Users/z3437171/Library/R/arm64/4.6/library/threejs",
        "tibble": "/Users/z3437171/Library/R/arm64/4.6/library/tibble",
        "tidygraph": "/Users/z3437171/Library/R/arm64/4.6/library/tidygraph",
        "tidyr": "/Users/z3437171/Library/R/arm64/4.6/library/tidyr",
        "tidyselect": "/Users/z3437171/Library/R/arm64/4.6/library/tidyselect",
        "tidyverse": "/Users/z3437171/Library/R/arm64/4.6/library/tidyverse",
        "timeDate": "/Users/z3437171/Library/R/arm64/4.6/library/timeDate",
        "timechange": "/Users/z3437171/Library/R/arm64/4.6/library/timechange",
        "tinytex": "/Users/z3437171/Library/R/arm64/4.6/library/tinytex",
        "torch": "/Users/z3437171/Library/R/arm64/4.6/library/torch",
        "triebeard": "/Users/z3437171/Library/R/arm64/4.6/library/triebeard",
        "truncnorm": "/Users/z3437171/Library/R/arm64/4.6/library/truncnorm",
        "tweedie": "/Users/z3437171/Library/R/arm64/4.6/library/tweedie",
        "tweenr": "/Users/z3437171/Library/R/arm64/4.6/library/tweenr",
        "tzdb": "/Users/z3437171/Library/R/arm64/4.6/library/tzdb",
        "ucminf": "/Users/z3437171/Library/R/arm64/4.6/library/ucminf",
        "units": "/Users/z3437171/Library/R/arm64/4.6/library/units",
        "urca": "/Users/z3437171/Library/R/arm64/4.6/library/urca",
        "urlchecker": "/Users/z3437171/Library/R/arm64/4.6/library/urlchecker",
        "urltools": "/Users/z3437171/Library/R/arm64/4.6/library/urltools",
        "usethis": "/Users/z3437171/Library/R/arm64/4.6/library/usethis",
        "utf8": "/Users/z3437171/Library/R/arm64/4.6/library/utf8",
        "uuid": "/Users/z3437171/Library/R/arm64/4.6/library/uuid",
        "vctrs": "/Users/z3437171/Library/R/arm64/4.6/library/vctrs",
        "vdiffr": "/Users/z3437171/Library/R/arm64/4.6/library/vdiffr",
        "vegan": "/Users/z3437171/Library/R/arm64/4.6/library/vegan",
        "vipor": "/Users/z3437171/Library/R/arm64/4.6/library/vipor",
        "viridis": "/Users/z3437171/Library/R/arm64/4.6/library/viridis",
        "viridisLite": "/Users/z3437171/Library/R/arm64/4.6/library/viridisLite",
        "visNetwork": "/Users/z3437171/Library/R/arm64/4.6/library/visNetwork",
        "vroom": "/Users/z3437171/Library/R/arm64/4.6/library/vroom",
        "waldo": "/Users/z3437171/Library/R/arm64/4.6/library/waldo",
        "whisker": "/Users/z3437171/Library/R/arm64/4.6/library/whisker",
        "whoami": "/Users/z3437171/Library/R/arm64/4.6/library/whoami",
        "withr": "/Users/z3437171/Library/R/arm64/4.6/library/withr",
        "wk": "/Users/z3437171/Library/R/arm64/4.6/library/wk",
        "xfun": "/Users/z3437171/Library/R/arm64/4.6/library/xfun",
        "xgboost": "/Users/z3437171/Library/R/arm64/4.6/library/xgboost",
        "xml2": "/Users/z3437171/Library/R/arm64/4.6/library/xml2",
        "xmlparsedata": "/Users/z3437171/Library/R/arm64/4.6/library/xmlparsedata",
        "xopen": "/Users/z3437171/Library/R/arm64/4.6/library/xopen",
        "xtable": "/Users/z3437171/Library/R/arm64/4.6/library/xtable",
        "xts": "/Users/z3437171/Library/R/arm64/4.6/library/xts",
        "yaml": "/Users/z3437171/Library/R/arm64/4.6/library/yaml",
        "zigg": "/Users/z3437171/Library/R/arm64/4.6/library/zigg",
        "zip": "/Users/z3437171/Library/R/arm64/4.6/library/zip",
        "zoo": "/Users/z3437171/Library/R/arm64/4.6/library/zoo"
      },
      "steps": [
        {
          "step": "install",
          "command": [
            "R",
            "CMD",
            "INSTALL",
            "--no-multiarch",
            "-l",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/drm-only",
            "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:00:19.733006+00:00",
          "elapsed_seconds": 24.57,
          "log": "drm-only-install.log",
          "log_sha256": "7f8a73b86b7d03565ff3003762e097a22aabf7fce7f9737a57df3161adc95993",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "fixtures",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "fixtures",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/drm-only"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:00:44.310167+00:00",
          "elapsed_seconds": 3.369,
          "log": "drm-only-fixtures.log",
          "log_sha256": "cb398da193bd8c6eb8c8ecf46cdf7b115e328991a7f5aad39c66efcc529ecb4e",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "real",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "drm-only",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/drm-only"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:00:47.688865+00:00",
          "elapsed_seconds": 7.304,
          "log": "drm-only-real.log",
          "log_sha256": "4d4cca3f2a713b1fb8ddc631347aab04484f9dc655d6335cb65e4ccab0b4ba6b",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        }
      ]
    },
    {
      "cell": "gllvm-only",
      "mode": "gllvm-only",
      "library": "/private/tmp/pigauto-mi-matrix-c2e768ad/gllvm-only",
      "expected_optional_packages": [
        "gllvmTMB"
      ],
      "dependency_symlinks": {
        "AGHmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/AGHmatrix",
        "BACE": "/Users/z3437171/Library/R/arm64/4.6/library/BACE",
        "BH": "/Users/z3437171/Library/R/arm64/4.6/library/BH",
        "BayesLogit": "/Users/z3437171/Library/R/arm64/4.6/library/BayesLogit",
        "BiocGenerics": "/Users/z3437171/Library/R/arm64/4.6/library/BiocGenerics",
        "BiocManager": "/Users/z3437171/Library/R/arm64/4.6/library/BiocManager",
        "BiocVersion": "/Users/z3437171/Library/R/arm64/4.6/library/BiocVersion",
        "Brobdingnag": "/Users/z3437171/Library/R/arm64/4.6/library/Brobdingnag",
        "DBI": "/Users/z3437171/Library/R/arm64/4.6/library/DBI",
        "DEoptim": "/Users/z3437171/Library/R/arm64/4.6/library/DEoptim",
        "DHARMa": "/Users/z3437171/Library/R/arm64/4.6/library/DHARMa",
        "DT": "/Users/z3437171/Library/R/arm64/4.6/library/DT",
        "Deriv": "/Users/z3437171/Library/R/arm64/4.6/library/Deriv",
        "DiagrammeR": "/Users/z3437171/Library/R/arm64/4.6/library/DiagrammeR",
        "FNN": "/Users/z3437171/Library/R/arm64/4.6/library/FNN",
        "Formula": "/Users/z3437171/Library/R/arm64/4.6/library/Formula",
        "GPArotation": "/Users/z3437171/Library/R/arm64/4.6/library/GPArotation",
        "HSAUR3": "/Users/z3437171/Library/R/arm64/4.6/library/HSAUR3",
        "Hmsc": "/Users/z3437171/Library/R/arm64/4.6/library/Hmsc",
        "JuliaCall": "/Users/z3437171/Library/R/arm64/4.6/library/JuliaCall",
        "MCMCglmm": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCglmm",
        "MCMCpack": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCpack",
        "MatrixModels": "/Users/z3437171/Library/R/arm64/4.6/library/MatrixModels",
        "MuMIn": "/Users/z3437171/Library/R/arm64/4.6/library/MuMIn",
        "QGglmm": "/Users/z3437171/Library/R/arm64/4.6/library/QGglmm",
        "QuickJSR": "/Users/z3437171/Library/R/arm64/4.6/library/QuickJSR",
        "R.cache": "/Users/z3437171/Library/R/arm64/4.6/library/R.cache",
        "R.methodsS3": "/Users/z3437171/Library/R/arm64/4.6/library/R.methodsS3",
        "R.oo": "/Users/z3437171/Library/R/arm64/4.6/library/R.oo",
        "R.rsp": "/Users/z3437171/Library/R/arm64/4.6/library/R.rsp",
        "R.utils": "/Users/z3437171/Library/R/arm64/4.6/library/R.utils",
        "R6": "/Users/z3437171/Library/R/arm64/4.6/library/R6",
        "RColorBrewer": "/Users/z3437171/Library/R/arm64/4.6/library/RColorBrewer",
        "ROI": "/Users/z3437171/Library/R/arm64/4.6/library/ROI",
        "ROI.plugin.lpsolve": "/Users/z3437171/Library/R/arm64/4.6/library/ROI.plugin.lpsolve",
        "RSpectra": "/Users/z3437171/Library/R/arm64/4.6/library/RSpectra",
        "RTMB": "/Users/z3437171/Library/R/arm64/4.6/library/RTMB",
        "RTMBdist": "/Users/z3437171/Library/R/arm64/4.6/library/RTMBdist",
        "RWiener": "/Users/z3437171/Library/R/arm64/4.6/library/RWiener",
        "Rcpp": "/Users/z3437171/Library/R/arm64/4.6/library/Rcpp",
        "RcppArmadillo": "/Users/z3437171/Library/R/arm64/4.6/library/RcppArmadillo",
        "RcppEigen": "/Users/z3437171/Library/R/arm64/4.6/library/RcppEigen",
        "RcppParallel": "/Users/z3437171/Library/R/arm64/4.6/library/RcppParallel",
        "RcppProgress": "/Users/z3437171/Library/R/arm64/4.6/library/RcppProgress",
        "Rdpack": "/Users/z3437171/Library/R/arm64/4.6/library/Rdpack",
        "Rfast": "/Users/z3437171/Library/R/arm64/4.6/library/Rfast",
        "Rmpfr": "/Users/z3437171/Library/R/arm64/4.6/library/Rmpfr",
        "Rphylopars": "/Users/z3437171/Library/R/arm64/4.6/library/Rphylopars",
        "S7": "/Users/z3437171/Library/R/arm64/4.6/library/S7",
        "SimDesign": "/Users/z3437171/Library/R/arm64/4.6/library/SimDesign",
        "SparseM": "/Users/z3437171/Library/R/arm64/4.6/library/SparseM",
        "StanHeaders": "/Users/z3437171/Library/R/arm64/4.6/library/StanHeaders",
        "TH.data": "/Users/z3437171/Library/R/arm64/4.6/library/TH.data",
        "TMB": "/Users/z3437171/Library/R/arm64/4.6/library/TMB",
        "V8": "/Users/z3437171/Library/R/arm64/4.6/library/V8",
        "VGAM": "/Users/z3437171/Library/R/arm64/4.6/library/VGAM",
        "XML": "/Users/z3437171/Library/R/arm64/4.6/library/XML",
        "_cache": "/Users/z3437171/Library/R/arm64/4.6/library/_cache",
        "abind": "/Users/z3437171/Library/R/arm64/4.6/library/abind",
        "ade4": "/Users/z3437171/Library/R/arm64/4.6/library/ade4",
        "alabama": "/Users/z3437171/Library/R/arm64/4.6/library/alabama",
        "ape": "/Users/z3437171/Library/R/arm64/4.6/library/ape",
        "arm": "/Users/z3437171/Library/R/arm64/4.6/library/arm",
        "askpass": "/Users/z3437171/Library/R/arm64/4.6/library/askpass",
        "assertthat": "/Users/z3437171/Library/R/arm64/4.6/library/assertthat",
        "audio": "/Users/z3437171/Library/R/arm64/4.6/library/audio",
        "backports": "/Users/z3437171/Library/R/arm64/4.6/library/backports",
        "base64enc": "/Users/z3437171/Library/R/arm64/4.6/library/base64enc",
        "bayeslincom": "/Users/z3437171/Library/R/arm64/4.6/library/bayeslincom",
        "bayesplot": "/Users/z3437171/Library/R/arm64/4.6/library/bayesplot",
        "bayestestR": "/Users/z3437171/Library/R/arm64/4.6/library/bayestestR",
        "bbmle": "/Users/z3437171/Library/R/arm64/4.6/library/bbmle",
        "bdsmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/bdsmatrix",
        "beepr": "/Users/z3437171/Library/R/arm64/4.6/library/beepr",
        "beeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/beeswarm",
        "betareg": "/Users/z3437171/Library/R/arm64/4.6/library/betareg",
        "bit": "/Users/z3437171/Library/R/arm64/4.6/library/bit",
        "bit64": "/Users/z3437171/Library/R/arm64/4.6/library/bit64",
        "blob": "/Users/z3437171/Library/R/arm64/4.6/library/blob",
        "blsmeta": "/Users/z3437171/Library/R/arm64/4.6/library/blsmeta",
        "brew": "/Users/z3437171/Library/R/arm64/4.6/library/brew",
        "brglm2": "/Users/z3437171/Library/R/arm64/4.6/library/brglm2",
        "bridgesampling": "/Users/z3437171/Library/R/arm64/4.6/library/bridgesampling",
        "brio": "/Users/z3437171/Library/R/arm64/4.6/library/brio",
        "brms": "/Users/z3437171/Library/R/arm64/4.6/library/brms",
        "broom": "/Users/z3437171/Library/R/arm64/4.6/library/broom",
        "broom.mixed": "/Users/z3437171/Library/R/arm64/4.6/library/broom.mixed",
        "bslib": "/Users/z3437171/Library/R/arm64/4.6/library/bslib",
        "cachem": "/Users/z3437171/Library/R/arm64/4.6/library/cachem",
        "callr": "/Users/z3437171/Library/R/arm64/4.6/library/callr",
        "car": "/Users/z3437171/Library/R/arm64/4.6/library/car",
        "carData": "/Users/z3437171/Library/R/arm64/4.6/library/carData",
        "castor": "/Users/z3437171/Library/R/arm64/4.6/library/castor",
        "cellranger": "/Users/z3437171/Library/R/arm64/4.6/library/cellranger",
        "checkmate": "/Users/z3437171/Library/R/arm64/4.6/library/checkmate",
        "circular": "/Users/z3437171/Library/R/arm64/4.6/library/circular",
        "classInt": "/Users/z3437171/Library/R/arm64/4.6/library/classInt",
        "cli": "/Users/z3437171/Library/R/arm64/4.6/library/cli",
        "clipr": "/Users/z3437171/Library/R/arm64/4.6/library/clipr",
        "clue": "/Users/z3437171/Library/R/arm64/4.6/library/clue",
        "clusterGeneration": "/Users/z3437171/Library/R/arm64/4.6/library/clusterGeneration",
        "cmdstanr": "/Users/z3437171/Library/R/arm64/4.6/library/cmdstanr",
        "coda": "/Users/z3437171/Library/R/arm64/4.6/library/coda",
        "collections": "/Users/z3437171/Library/R/arm64/4.6/library/collections",
        "colorspace": "/Users/z3437171/Library/R/arm64/4.6/library/colorspace",
        "colourpicker": "/Users/z3437171/Library/R/arm64/4.6/library/colourpicker",
        "combinat": "/Users/z3437171/Library/R/arm64/4.6/library/combinat",
        "commonmark": "/Users/z3437171/Library/R/arm64/4.6/library/commonmark",
        "conflicted": "/Users/z3437171/Library/R/arm64/4.6/library/conflicted",
        "coro": "/Users/z3437171/Library/R/arm64/4.6/library/coro",
        "corpcor": "/Users/z3437171/Library/R/arm64/4.6/library/corpcor",
        "covr": "/Users/z3437171/Library/R/arm64/4.6/library/covr",
        "cowplot": "/Users/z3437171/Library/R/arm64/4.6/library/cowplot",
        "cpp11": "/Users/z3437171/Library/R/arm64/4.6/library/cpp11",
        "crayon": "/Users/z3437171/Library/R/arm64/4.6/library/crayon",
        "credentials": "/Users/z3437171/Library/R/arm64/4.6/library/credentials",
        "crosstalk": "/Users/z3437171/Library/R/arm64/4.6/library/crosstalk",
        "crul": "/Users/z3437171/Library/R/arm64/4.6/library/crul",
        "cubature": "/Users/z3437171/Library/R/arm64/4.6/library/cubature",
        "curl": "/Users/z3437171/Library/R/arm64/4.6/library/curl",
        "dagitty": "/Users/z3437171/Library/R/arm64/4.6/library/dagitty",
        "data.table": "/Users/z3437171/Library/R/arm64/4.6/library/data.table",
        "datawizard": "/Users/z3437171/Library/R/arm64/4.6/library/datawizard",
        "dbplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dbplyr",
        "dcurver": "/Users/z3437171/Library/R/arm64/4.6/library/dcurver",
        "deldir": "/Users/z3437171/Library/R/arm64/4.6/library/deldir",
        "desc": "/Users/z3437171/Library/R/arm64/4.6/library/desc",
        "detectseparation": "/Users/z3437171/Library/R/arm64/4.6/library/detectseparation",
        "devtools": "/Users/z3437171/Library/R/arm64/4.6/library/devtools",
        "diffobj": "/Users/z3437171/Library/R/arm64/4.6/library/diffobj",
        "digest": "/Users/z3437171/Library/R/arm64/4.6/library/digest",
        "distributional": "/Users/z3437171/Library/R/arm64/4.6/library/distributional",
        "doBy": "/Users/z3437171/Library/R/arm64/4.6/library/doBy",
        "doFuture": "/Users/z3437171/Library/R/arm64/4.6/library/doFuture",
        "doParallel": "/Users/z3437171/Library/R/arm64/4.6/library/doParallel",
        "doRNG": "/Users/z3437171/Library/R/arm64/4.6/library/doRNG",
        "dotCall64": "/Users/z3437171/Library/R/arm64/4.6/library/dotCall64",
        "downlit": "/Users/z3437171/Library/R/arm64/4.6/library/downlit",
        "dplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dplyr",
        "drmSEM": "/Users/z3437171/Library/R/arm64/4.6/library/drmSEM",
        "dtplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dtplyr",
        "dygraphs": "/Users/z3437171/Library/R/arm64/4.6/library/dygraphs",
        "e1071": "/Users/z3437171/Library/R/arm64/4.6/library/e1071",
        "ellipsis": "/Users/z3437171/Library/R/arm64/4.6/library/ellipsis",
        "emmeans": "/Users/z3437171/Library/R/arm64/4.6/library/emmeans",
        "enhancer": "/Users/z3437171/Library/R/arm64/4.6/library/enhancer",
        "enrichwith": "/Users/z3437171/Library/R/arm64/4.6/library/enrichwith",
        "estimability": "/Users/z3437171/Library/R/arm64/4.6/library/estimability",
        "evaluate": "/Users/z3437171/Library/R/arm64/4.6/library/evaluate",
        "evd": "/Users/z3437171/Library/R/arm64/4.6/library/evd",
        "expm": "/Users/z3437171/Library/R/arm64/4.6/library/expm",
        "extraDistr": "/Users/z3437171/Library/R/arm64/4.6/library/extraDistr",
        "fansi": "/Users/z3437171/Library/R/arm64/4.6/library/fansi",
        "farver": "/Users/z3437171/Library/R/arm64/4.6/library/farver",
        "fastGHQuad": "/Users/z3437171/Library/R/arm64/4.6/library/fastGHQuad",
        "fastmap": "/Users/z3437171/Library/R/arm64/4.6/library/fastmap",
        "fastmatch": "/Users/z3437171/Library/R/arm64/4.6/library/fastmatch",
        "fields": "/Users/z3437171/Library/R/arm64/4.6/library/fields",
        "fishMod": "/Users/z3437171/Library/R/arm64/4.6/library/fishMod",
        "flexmix": "/Users/z3437171/Library/R/arm64/4.6/library/flexmix",
        "fmesher": "/Users/z3437171/Library/R/arm64/4.6/library/fmesher",
        "fontawesome": "/Users/z3437171/Library/R/arm64/4.6/library/fontawesome",
        "forcats": "/Users/z3437171/Library/R/arm64/4.6/library/forcats",
        "foreach": "/Users/z3437171/Library/R/arm64/4.6/library/foreach",
        "forecast": "/Users/z3437171/Library/R/arm64/4.6/library/forecast",
        "fracdiff": "/Users/z3437171/Library/R/arm64/4.6/library/fracdiff",
        "freqTLS": "/Users/z3437171/Library/R/arm64/4.6/library/freqTLS",
        "fs": "/Users/z3437171/Library/R/arm64/4.6/library/fs",
        "furrr": "/Users/z3437171/Library/R/arm64/4.6/library/furrr",
        "future": "/Users/z3437171/Library/R/arm64/4.6/library/future",
        "future.apply": "/Users/z3437171/Library/R/arm64/4.6/library/future.apply",
        "future.mirai": "/Users/z3437171/Library/R/arm64/4.6/library/future.mirai",
        "galamm": "/Users/z3437171/Library/R/arm64/4.6/library/galamm",
        "gamlss.dist": "/Users/z3437171/Library/R/arm64/4.6/library/gamlss.dist",
        "gamm4": "/Users/z3437171/Library/R/arm64/4.6/library/gamm4",
        "gap": "/Users/z3437171/Library/R/arm64/4.6/library/gap",
        "gap.datasets": "/Users/z3437171/Library/R/arm64/4.6/library/gap.datasets",
        "gargle": "/Users/z3437171/Library/R/arm64/4.6/library/gargle",
        "generics": "/Users/z3437171/Library/R/arm64/4.6/library/generics",
        "geojsonsf": "/Users/z3437171/Library/R/arm64/4.6/library/geojsonsf",
        "geometries": "/Users/z3437171/Library/R/arm64/4.6/library/geometries",
        "gert": "/Users/z3437171/Library/R/arm64/4.6/library/gert",
        "ggbeeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/ggbeeswarm",
        "ggdist": "/Users/z3437171/Library/R/arm64/4.6/library/ggdist",
        "ggfittext": "/Users/z3437171/Library/R/arm64/4.6/library/ggfittext",
        "ggforce": "/Users/z3437171/Library/R/arm64/4.6/library/ggforce",
        "gggenes": "/Users/z3437171/Library/R/arm64/4.6/library/gggenes",
        "ggh4x": "/Users/z3437171/Library/R/arm64/4.6/library/ggh4x",
        "ggm": "/Users/z3437171/Library/R/arm64/4.6/library/ggm",
        "ggokabeito": "/Users/z3437171/Library/R/arm64/4.6/library/ggokabeito",
        "ggplot2": "/Users/z3437171/Library/R/arm64/4.6/library/ggplot2",
        "ggraph": "/Users/z3437171/Library/R/arm64/4.6/library/ggraph",
        "ggrepel": "/Users/z3437171/Library/R/arm64/4.6/library/ggrepel",
        "ggridges": "/Users/z3437171/Library/R/arm64/4.6/library/ggridges",
        "gh": "/Users/z3437171/Library/R/arm64/4.6/library/gh",
        "gitcreds": "/Users/z3437171/Library/R/arm64/4.6/library/gitcreds",
        "glassoFast": "/Users/z3437171/Library/R/arm64/4.6/library/glassoFast",
        "gllvm": "/Users/z3437171/Library/R/arm64/4.6/library/gllvm",
        "gllvmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/gllvmTMB",
        "glmmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/glmmTMB",
        "glmnet": "/Users/z3437171/Library/R/arm64/4.6/library/glmnet",
        "globals": "/Users/z3437171/Library/R/arm64/4.6/library/globals",
        "glue": "/Users/z3437171/Library/R/arm64/4.6/library/glue",
        "gmp": "/Users/z3437171/Library/R/arm64/4.6/library/gmp",
        "googledrive": "/Users/z3437171/Library/R/arm64/4.6/library/googledrive",
        "googlesheets4": "/Users/z3437171/Library/R/arm64/4.6/library/googlesheets4",
        "graph": "/Users/z3437171/Library/R/arm64/4.6/library/graph",
        "graphlayouts": "/Users/z3437171/Library/R/arm64/4.6/library/graphlayouts",
        "gratia": "/Users/z3437171/Library/R/arm64/4.6/library/gratia",
        "gridExtra": "/Users/z3437171/Library/R/arm64/4.6/library/gridExtra",
        "gridtext": "/Users/z3437171/Library/R/arm64/4.6/library/gridtext",
        "gsl": "/Users/z3437171/Library/R/arm64/4.6/library/gsl",
        "gtable": "/Users/z3437171/Library/R/arm64/4.6/library/gtable",
        "gtools": "/Users/z3437171/Library/R/arm64/4.6/library/gtools",
        "haven": "/Users/z3437171/Library/R/arm64/4.6/library/haven",
        "here": "/Users/z3437171/Library/R/arm64/4.6/library/here",
        "highr": "/Users/z3437171/Library/R/arm64/4.6/library/highr",
        "hms": "/Users/z3437171/Library/R/arm64/4.6/library/hms",
        "hsquared": "/Users/z3437171/Library/R/arm64/4.6/library/hsquared",
        "htmltools": "/Users/z3437171/Library/R/arm64/4.6/library/htmltools",
        "htmlwidgets": "/Users/z3437171/Library/R/arm64/4.6/library/htmlwidgets",
        "httpcode": "/Users/z3437171/Library/R/arm64/4.6/library/httpcode",
        "httpuv": "/Users/z3437171/Library/R/arm64/4.6/library/httpuv",
        "httr": "/Users/z3437171/Library/R/arm64/4.6/library/httr",
        "httr2": "/Users/z3437171/Library/R/arm64/4.6/library/httr2",
        "hunspell": "/Users/z3437171/Library/R/arm64/4.6/library/hunspell",
        "ids": "/Users/z3437171/Library/R/arm64/4.6/library/ids",
        "igraph": "/Users/z3437171/Library/R/arm64/4.6/library/igraph",
        "ini": "/Users/z3437171/Library/R/arm64/4.6/library/ini",
        "inline": "/Users/z3437171/Library/R/arm64/4.6/library/inline",
        "insight": "/Users/z3437171/Library/R/arm64/4.6/library/insight",
        "isoband": "/Users/z3437171/Library/R/arm64/4.6/library/isoband",
        "iterators": "/Users/z3437171/Library/R/arm64/4.6/library/iterators",
        "itertools": "/Users/z3437171/Library/R/arm64/4.6/library/itertools",
        "janitor": "/Users/z3437171/Library/R/arm64/4.6/library/janitor",
        "jomo": "/Users/z3437171/Library/R/arm64/4.6/library/jomo",
        "jpeg": "/Users/z3437171/Library/R/arm64/4.6/library/jpeg",
        "jquerylib": "/Users/z3437171/Library/R/arm64/4.6/library/jquerylib",
        "jsonify": "/Users/z3437171/Library/R/arm64/4.6/library/jsonify",
        "jsonlite": "/Users/z3437171/Library/R/arm64/4.6/library/jsonlite",
        "kernelshap": "/Users/z3437171/Library/R/arm64/4.6/library/kernelshap",
        "knitr": "/Users/z3437171/Library/R/arm64/4.6/library/knitr",
        "labeling": "/Users/z3437171/Library/R/arm64/4.6/library/labeling",
        "languageserver": "/Users/z3437171/Library/R/arm64/4.6/library/languageserver",
        "later": "/Users/z3437171/Library/R/arm64/4.6/library/later",
        "latex2exp": "/Users/z3437171/Library/R/arm64/4.6/library/latex2exp",
        "lavaan": "/Users/z3437171/Library/R/arm64/4.6/library/lavaan",
        "lazyeval": "/Users/z3437171/Library/R/arm64/4.6/library/lazyeval",
        "leafem": "/Users/z3437171/Library/R/arm64/4.6/library/leafem",
        "leaflet": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet",
        "leaflet.providers": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet.providers",
        "leafpop": "/Users/z3437171/Library/R/arm64/4.6/library/leafpop",
        "lifecycle": "/Users/z3437171/Library/R/arm64/4.6/library/lifecycle",
        "lintr": "/Users/z3437171/Library/R/arm64/4.6/library/lintr",
        "listenv": "/Users/z3437171/Library/R/arm64/4.6/library/listenv",
        "litedown": "/Users/z3437171/Library/R/arm64/4.6/library/litedown",
        "lme4": "/Users/z3437171/Library/R/arm64/4.6/library/lme4",
        "lmtest": "/Users/z3437171/Library/R/arm64/4.6/library/lmtest",
        "loo": "/Users/z3437171/Library/R/arm64/4.6/library/loo",
        "lpSolveAPI": "/Users/z3437171/Library/R/arm64/4.6/library/lpSolveAPI",
        "lubridate": "/Users/z3437171/Library/R/arm64/4.6/library/lubridate",
        "magrittr": "/Users/z3437171/Library/R/arm64/4.6/library/magrittr",
        "maps": "/Users/z3437171/Library/R/arm64/4.6/library/maps",
        "mapview": "/Users/z3437171/Library/R/arm64/4.6/library/mapview",
        "markdown": "/Users/z3437171/Library/R/arm64/4.6/library/markdown",
        "mathjaxr": "/Users/z3437171/Library/R/arm64/4.6/library/mathjaxr",
        "matrixStats": "/Users/z3437171/Library/R/arm64/4.6/library/matrixStats",
        "mclogit": "/Users/z3437171/Library/R/arm64/4.6/library/mclogit",
        "mcmc": "/Users/z3437171/Library/R/arm64/4.6/library/mcmc",
        "memisc": "/Users/z3437171/Library/R/arm64/4.6/library/memisc",
        "memoise": "/Users/z3437171/Library/R/arm64/4.6/library/memoise",
        "metadat": "/Users/z3437171/Library/R/arm64/4.6/library/metadat",
        "metafor": "/Users/z3437171/Library/R/arm64/4.6/library/metafor",
        "metaplus": "/Users/z3437171/Library/R/arm64/4.6/library/metaplus",
        "mice": "/Users/z3437171/Library/R/arm64/4.6/library/mice",
        "microbenchmark": "/Users/z3437171/Library/R/arm64/4.6/library/microbenchmark",
        "mime": "/Users/z3437171/Library/R/arm64/4.6/library/mime",
        "miniUI": "/Users/z3437171/Library/R/arm64/4.6/library/miniUI",
        "minqa": "/Users/z3437171/Library/R/arm64/4.6/library/minqa",
        "mirai": "/Users/z3437171/Library/R/arm64/4.6/library/mirai",
        "mirt": "/Users/z3437171/Library/R/arm64/4.6/library/mirt",
        "missForest": "/Users/z3437171/Library/R/arm64/4.6/library/missForest",
        "mitml": "/Users/z3437171/Library/R/arm64/4.6/library/mitml",
        "mnormt": "/Users/z3437171/Library/R/arm64/4.6/library/mnormt",
        "modelr": "/Users/z3437171/Library/R/arm64/4.6/library/modelr",
        "modeltools": "/Users/z3437171/Library/R/arm64/4.6/library/modeltools",
        "movMF": "/Users/z3437171/Library/R/arm64/4.6/library/movMF",
        "msm": "/Users/z3437171/Library/R/arm64/4.6/library/msm",
        "multcomp": "/Users/z3437171/Library/R/arm64/4.6/library/multcomp",
        "mvMORPH": "/Users/z3437171/Library/R/arm64/4.6/library/mvMORPH",
        "mvnfast": "/Users/z3437171/Library/R/arm64/4.6/library/mvnfast",
        "mvtnorm": "/Users/z3437171/Library/R/arm64/4.6/library/mvtnorm",
        "nadiv": "/Users/z3437171/Library/R/arm64/4.6/library/nadiv",
        "nanonext": "/Users/z3437171/Library/R/arm64/4.6/library/nanonext",
        "naturalsort": "/Users/z3437171/Library/R/arm64/4.6/library/naturalsort",
        "ncf": "/Users/z3437171/Library/R/arm64/4.6/library/ncf",
        "nleqslv": "/Users/z3437171/Library/R/arm64/4.6/library/nleqslv",
        "nloptr": "/Users/z3437171/Library/R/arm64/4.6/library/nloptr",
        "numDeriv": "/Users/z3437171/Library/R/arm64/4.6/library/numDeriv",
        "oai": "/Users/z3437171/Library/R/arm64/4.6/library/oai",
        "openssl": "/Users/z3437171/Library/R/arm64/4.6/library/openssl",
        "optimParallel": "/Users/z3437171/Library/R/arm64/4.6/library/optimParallel",
        "orchaRd": "/Users/z3437171/Library/R/arm64/4.6/library/orchaRd",
        "ordinal": "/Users/z3437171/Library/R/arm64/4.6/library/ordinal",
        "orthopolynom": "/Users/z3437171/Library/R/arm64/4.6/library/orthopolynom",
        "otel": "/Users/z3437171/Library/R/arm64/4.6/library/otel",
        "pROC": "/Users/z3437171/Library/R/arm64/4.6/library/pROC",
        "pak": "/Users/z3437171/Library/R/arm64/4.6/library/pak",
        "palmerpenguins": "/Users/z3437171/Library/R/arm64/4.6/library/palmerpenguins",
        "pan": "/Users/z3437171/Library/R/arm64/4.6/library/pan",
        "parallelly": "/Users/z3437171/Library/R/arm64/4.6/library/parallelly",
        "patchwork": "/Users/z3437171/Library/R/arm64/4.6/library/patchwork",
        "pbapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbapply",
        "pbivnorm": "/Users/z3437171/Library/R/arm64/4.6/library/pbivnorm",
        "pbkrtest": "/Users/z3437171/Library/R/arm64/4.6/library/pbkrtest",
        "pbmcapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbmcapply",
        "pedigreemm": "/Users/z3437171/Library/R/arm64/4.6/library/pedigreemm",
        "performance": "/Users/z3437171/Library/R/arm64/4.6/library/performance",
        "permute": "/Users/z3437171/Library/R/arm64/4.6/library/permute",
        "phangorn": "/Users/z3437171/Library/R/arm64/4.6/library/phangorn",
        "phylolm": "/Users/z3437171/Library/R/arm64/4.6/library/phylolm",
        "phylopath": "/Users/z3437171/Library/R/arm64/4.6/library/phylopath",
        "phytools": "/Users/z3437171/Library/R/arm64/4.6/library/phytools",
        "piecewiseSEM": "/Users/z3437171/Library/R/arm64/4.6/library/piecewiseSEM",
        "pillar": "/Users/z3437171/Library/R/arm64/4.6/library/pillar",
        "pixmap": "/Users/z3437171/Library/R/arm64/4.6/library/pixmap",
        "pkgbuild": "/Users/z3437171/Library/R/arm64/4.6/library/pkgbuild",
        "pkgconfig": "/Users/z3437171/Library/R/arm64/4.6/library/pkgconfig",
        "pkgdown": "/Users/z3437171/Library/R/arm64/4.6/library/pkgdown",
        "pkgload": "/Users/z3437171/Library/R/arm64/4.6/library/pkgload",
        "plotly": "/Users/z3437171/Library/R/arm64/4.6/library/plotly",
        "plyr": "/Users/z3437171/Library/R/arm64/4.6/library/plyr",
        "png": "/Users/z3437171/Library/R/arm64/4.6/library/png",
        "polyclip": "/Users/z3437171/Library/R/arm64/4.6/library/polyclip",
        "polynom": "/Users/z3437171/Library/R/arm64/4.6/library/polynom",
        "posterior": "/Users/z3437171/Library/R/arm64/4.6/library/posterior",
        "pracma": "/Users/z3437171/Library/R/arm64/4.6/library/pracma",
        "praise": "/Users/z3437171/Library/R/arm64/4.6/library/praise",
        "prettyunits": "/Users/z3437171/Library/R/arm64/4.6/library/prettyunits",
        "priorsense": "/Users/z3437171/Library/R/arm64/4.6/library/priorsense",
        "processx": "/Users/z3437171/Library/R/arm64/4.6/library/processx",
        "profvis": "/Users/z3437171/Library/R/arm64/4.6/library/profvis",
        "progress": "/Users/z3437171/Library/R/arm64/4.6/library/progress",
        "progressr": "/Users/z3437171/Library/R/arm64/4.6/library/progressr",
        "projpred": "/Users/z3437171/Library/R/arm64/4.6/library/projpred",
        "promises": "/Users/z3437171/Library/R/arm64/4.6/library/promises",
        "proxy": "/Users/z3437171/Library/R/arm64/4.6/library/proxy",
        "ps": "/Users/z3437171/Library/R/arm64/4.6/library/ps",
        "psymetadata": "/Users/z3437171/Library/R/arm64/4.6/library/psymetadata",
        "purrr": "/Users/z3437171/Library/R/arm64/4.6/library/purrr",
        "qgam": "/Users/z3437171/Library/R/arm64/4.6/library/qgam",
        "qs2": "/Users/z3437171/Library/R/arm64/4.6/library/qs2",
        "quadprog": "/Users/z3437171/Library/R/arm64/4.6/library/quadprog",
        "quantreg": "/Users/z3437171/Library/R/arm64/4.6/library/quantreg",
        "rWishart": "/Users/z3437171/Library/R/arm64/4.6/library/rWishart",
        "ragg": "/Users/z3437171/Library/R/arm64/4.6/library/ragg",
        "randomForest": "/Users/z3437171/Library/R/arm64/4.6/library/randomForest",
        "ranger": "/Users/z3437171/Library/R/arm64/4.6/library/ranger",
        "rapidjsonr": "/Users/z3437171/Library/R/arm64/4.6/library/rapidjsonr",
        "rappdirs": "/Users/z3437171/Library/R/arm64/4.6/library/rappdirs",
        "raster": "/Users/z3437171/Library/R/arm64/4.6/library/raster",
        "rbibutils": "/Users/z3437171/Library/R/arm64/4.6/library/rbibutils",
        "rcmdcheck": "/Users/z3437171/Library/R/arm64/4.6/library/rcmdcheck",
        "readr": "/Users/z3437171/Library/R/arm64/4.6/library/readr",
        "readxl": "/Users/z3437171/Library/R/arm64/4.6/library/readxl",
        "reformulas": "/Users/z3437171/Library/R/arm64/4.6/library/reformulas",
        "registry": "/Users/z3437171/Library/R/arm64/4.6/library/registry",
        "rematch": "/Users/z3437171/Library/R/arm64/4.6/library/rematch",
        "rematch2": "/Users/z3437171/Library/R/arm64/4.6/library/rematch2",
        "remotes": "/Users/z3437171/Library/R/arm64/4.6/library/remotes",
        "rentrez": "/Users/z3437171/Library/R/arm64/4.6/library/rentrez",
        "reprex": "/Users/z3437171/Library/R/arm64/4.6/library/reprex",
        "reshape2": "/Users/z3437171/Library/R/arm64/4.6/library/reshape2",
        "rex": "/Users/z3437171/Library/R/arm64/4.6/library/rex",
        "rgbif": "/Users/z3437171/Library/R/arm64/4.6/library/rgbif",
        "rgpdd": "/Users/z3437171/Library/R/arm64/4.6/library/rgpdd",
        "rhub": "/Users/z3437171/Library/R/arm64/4.6/library/rhub",
        "rjags": "/Users/z3437171/Library/R/arm64/4.6/library/rjags",
        "rjson": "/Users/z3437171/Library/R/arm64/4.6/library/rjson",
        "rlang": "/Users/z3437171/Library/R/arm64/4.6/library/rlang",
        "rmarkdown": "/Users/z3437171/Library/R/arm64/4.6/library/rmarkdown",
        "rncl": "/Users/z3437171/Library/R/arm64/4.6/library/rncl",
        "rngtools": "/Users/z3437171/Library/R/arm64/4.6/library/rngtools",
        "rotl": "/Users/z3437171/Library/R/arm64/4.6/library/rotl",
        "roxygen2": "/Users/z3437171/Library/R/arm64/4.6/library/roxygen2",
        "rprojroot": "/Users/z3437171/Library/R/arm64/4.6/library/rprojroot",
        "rstan": "/Users/z3437171/Library/R/arm64/4.6/library/rstan",
        "rstantools": "/Users/z3437171/Library/R/arm64/4.6/library/rstantools",
        "rstudioapi": "/Users/z3437171/Library/R/arm64/4.6/library/rstudioapi",
        "rtdists": "/Users/z3437171/Library/R/arm64/4.6/library/rtdists",
        "rversions": "/Users/z3437171/Library/R/arm64/4.6/library/rversions",
        "rvest": "/Users/z3437171/Library/R/arm64/4.6/library/rvest",
        "s2": "/Users/z3437171/Library/R/arm64/4.6/library/s2",
        "safetensors": "/Users/z3437171/Library/R/arm64/4.6/library/safetensors",
        "sandwich": "/Users/z3437171/Library/R/arm64/4.6/library/sandwich",
        "sass": "/Users/z3437171/Library/R/arm64/4.6/library/sass",
        "satellite": "/Users/z3437171/Library/R/arm64/4.6/library/satellite",
        "scales": "/Users/z3437171/Library/R/arm64/4.6/library/scales",
        "scatterplot3d": "/Users/z3437171/Library/R/arm64/4.6/library/scatterplot3d",
        "selectr": "/Users/z3437171/Library/R/arm64/4.6/library/selectr",
        "servr": "/Users/z3437171/Library/R/arm64/4.6/library/servr",
        "sessioninfo": "/Users/z3437171/Library/R/arm64/4.6/library/sessioninfo",
        "sf": "/Users/z3437171/Library/R/arm64/4.6/library/sf",
        "sfheaders": "/Users/z3437171/Library/R/arm64/4.6/library/sfheaders",
        "shades": "/Users/z3437171/Library/R/arm64/4.6/library/shades",
        "shape": "/Users/z3437171/Library/R/arm64/4.6/library/shape",
        "shapviz": "/Users/z3437171/Library/R/arm64/4.6/library/shapviz",
        "shiny": "/Users/z3437171/Library/R/arm64/4.6/library/shiny",
        "shinyjs": "/Users/z3437171/Library/R/arm64/4.6/library/shinyjs",
        "shinystan": "/Users/z3437171/Library/R/arm64/4.6/library/shinystan",
        "shinythemes": "/Users/z3437171/Library/R/arm64/4.6/library/shinythemes",
        "skmeans": "/Users/z3437171/Library/R/arm64/4.6/library/skmeans",
        "slam": "/Users/z3437171/Library/R/arm64/4.6/library/slam",
        "smcfcs": "/Users/z3437171/Library/R/arm64/4.6/library/smcfcs",
        "sn": "/Users/z3437171/Library/R/arm64/4.6/library/sn",
        "snakecase": "/Users/z3437171/Library/R/arm64/4.6/library/snakecase",
        "sommer": "/Users/z3437171/Library/R/arm64/4.6/library/sommer",
        "sourcetools": "/Users/z3437171/Library/R/arm64/4.6/library/sourcetools",
        "sp": "/Users/z3437171/Library/R/arm64/4.6/library/sp",
        "spData": "/Users/z3437171/Library/R/arm64/4.6/library/spData",
        "spam": "/Users/z3437171/Library/R/arm64/4.6/library/spam",
        "spdep": "/Users/z3437171/Library/R/arm64/4.6/library/spdep",
        "spelling": "/Users/z3437171/Library/R/arm64/4.6/library/spelling",
        "splancs": "/Users/z3437171/Library/R/arm64/4.6/library/splancs",
        "splines2": "/Users/z3437171/Library/R/arm64/4.6/library/splines2",
        "statmod": "/Users/z3437171/Library/R/arm64/4.6/library/statmod",
        "stringfish": "/Users/z3437171/Library/R/arm64/4.6/library/stringfish",
        "stringi": "/Users/z3437171/Library/R/arm64/4.6/library/stringi",
        "stringr": "/Users/z3437171/Library/R/arm64/4.6/library/stringr",
        "styler": "/Users/z3437171/Library/R/arm64/4.6/library/styler",
        "subplex": "/Users/z3437171/Library/R/arm64/4.6/library/subplex",
        "svglite": "/Users/z3437171/Library/R/arm64/4.6/library/svglite",
        "sys": "/Users/z3437171/Library/R/arm64/4.6/library/sys",
        "systemfonts": "/Users/z3437171/Library/R/arm64/4.6/library/systemfonts",
        "tensorA": "/Users/z3437171/Library/R/arm64/4.6/library/tensorA",
        "terra": "/Users/z3437171/Library/R/arm64/4.6/library/terra",
        "tester": "/Users/z3437171/Library/R/arm64/4.6/library/tester",
        "testthat": "/Users/z3437171/Library/R/arm64/4.6/library/testthat",
        "textshaping": "/Users/z3437171/Library/R/arm64/4.6/library/textshaping",
        "threejs": "/Users/z3437171/Library/R/arm64/4.6/library/threejs",
        "tibble": "/Users/z3437171/Library/R/arm64/4.6/library/tibble",
        "tidygraph": "/Users/z3437171/Library/R/arm64/4.6/library/tidygraph",
        "tidyr": "/Users/z3437171/Library/R/arm64/4.6/library/tidyr",
        "tidyselect": "/Users/z3437171/Library/R/arm64/4.6/library/tidyselect",
        "tidyverse": "/Users/z3437171/Library/R/arm64/4.6/library/tidyverse",
        "timeDate": "/Users/z3437171/Library/R/arm64/4.6/library/timeDate",
        "timechange": "/Users/z3437171/Library/R/arm64/4.6/library/timechange",
        "tinytex": "/Users/z3437171/Library/R/arm64/4.6/library/tinytex",
        "torch": "/Users/z3437171/Library/R/arm64/4.6/library/torch",
        "triebeard": "/Users/z3437171/Library/R/arm64/4.6/library/triebeard",
        "truncnorm": "/Users/z3437171/Library/R/arm64/4.6/library/truncnorm",
        "tweedie": "/Users/z3437171/Library/R/arm64/4.6/library/tweedie",
        "tweenr": "/Users/z3437171/Library/R/arm64/4.6/library/tweenr",
        "tzdb": "/Users/z3437171/Library/R/arm64/4.6/library/tzdb",
        "ucminf": "/Users/z3437171/Library/R/arm64/4.6/library/ucminf",
        "units": "/Users/z3437171/Library/R/arm64/4.6/library/units",
        "urca": "/Users/z3437171/Library/R/arm64/4.6/library/urca",
        "urlchecker": "/Users/z3437171/Library/R/arm64/4.6/library/urlchecker",
        "urltools": "/Users/z3437171/Library/R/arm64/4.6/library/urltools",
        "usethis": "/Users/z3437171/Library/R/arm64/4.6/library/usethis",
        "utf8": "/Users/z3437171/Library/R/arm64/4.6/library/utf8",
        "uuid": "/Users/z3437171/Library/R/arm64/4.6/library/uuid",
        "vctrs": "/Users/z3437171/Library/R/arm64/4.6/library/vctrs",
        "vdiffr": "/Users/z3437171/Library/R/arm64/4.6/library/vdiffr",
        "vegan": "/Users/z3437171/Library/R/arm64/4.6/library/vegan",
        "vipor": "/Users/z3437171/Library/R/arm64/4.6/library/vipor",
        "viridis": "/Users/z3437171/Library/R/arm64/4.6/library/viridis",
        "viridisLite": "/Users/z3437171/Library/R/arm64/4.6/library/viridisLite",
        "visNetwork": "/Users/z3437171/Library/R/arm64/4.6/library/visNetwork",
        "vroom": "/Users/z3437171/Library/R/arm64/4.6/library/vroom",
        "waldo": "/Users/z3437171/Library/R/arm64/4.6/library/waldo",
        "whisker": "/Users/z3437171/Library/R/arm64/4.6/library/whisker",
        "whoami": "/Users/z3437171/Library/R/arm64/4.6/library/whoami",
        "withr": "/Users/z3437171/Library/R/arm64/4.6/library/withr",
        "wk": "/Users/z3437171/Library/R/arm64/4.6/library/wk",
        "xfun": "/Users/z3437171/Library/R/arm64/4.6/library/xfun",
        "xgboost": "/Users/z3437171/Library/R/arm64/4.6/library/xgboost",
        "xml2": "/Users/z3437171/Library/R/arm64/4.6/library/xml2",
        "xmlparsedata": "/Users/z3437171/Library/R/arm64/4.6/library/xmlparsedata",
        "xopen": "/Users/z3437171/Library/R/arm64/4.6/library/xopen",
        "xtable": "/Users/z3437171/Library/R/arm64/4.6/library/xtable",
        "xts": "/Users/z3437171/Library/R/arm64/4.6/library/xts",
        "yaml": "/Users/z3437171/Library/R/arm64/4.6/library/yaml",
        "zigg": "/Users/z3437171/Library/R/arm64/4.6/library/zigg",
        "zip": "/Users/z3437171/Library/R/arm64/4.6/library/zip",
        "zoo": "/Users/z3437171/Library/R/arm64/4.6/library/zoo"
      },
      "steps": [
        {
          "step": "install",
          "command": [
            "R",
            "CMD",
            "INSTALL",
            "--no-multiarch",
            "-l",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/gllvm-only",
            "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:00:55.396722+00:00",
          "elapsed_seconds": 18.007,
          "log": "gllvm-only-install.log",
          "log_sha256": "7f8a73b86b7d03565ff3003762e097a22aabf7fce7f9737a57df3161adc95993",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "fixtures",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "fixtures",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/gllvm-only"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:01:13.411476+00:00",
          "elapsed_seconds": 3.985,
          "log": "gllvm-only-fixtures.log",
          "log_sha256": "cb398da193bd8c6eb8c8ecf46cdf7b115e328991a7f5aad39c66efcc529ecb4e",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "real",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "gllvm-only",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/gllvm-only"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:01:17.404035+00:00",
          "elapsed_seconds": 7.225,
          "log": "gllvm-only-real.log",
          "log_sha256": "e32344192ae01b7d8a0a6085b8271b180421cde77e6906023e2b4077f1c3729e",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        }
      ]
    },
    {
      "cell": "both",
      "mode": "real",
      "library": "/private/tmp/pigauto-mi-matrix-c2e768ad/both",
      "expected_optional_packages": [
        "drmTMB",
        "gllvmTMB"
      ],
      "dependency_symlinks": {
        "AGHmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/AGHmatrix",
        "BACE": "/Users/z3437171/Library/R/arm64/4.6/library/BACE",
        "BH": "/Users/z3437171/Library/R/arm64/4.6/library/BH",
        "BayesLogit": "/Users/z3437171/Library/R/arm64/4.6/library/BayesLogit",
        "BiocGenerics": "/Users/z3437171/Library/R/arm64/4.6/library/BiocGenerics",
        "BiocManager": "/Users/z3437171/Library/R/arm64/4.6/library/BiocManager",
        "BiocVersion": "/Users/z3437171/Library/R/arm64/4.6/library/BiocVersion",
        "Brobdingnag": "/Users/z3437171/Library/R/arm64/4.6/library/Brobdingnag",
        "DBI": "/Users/z3437171/Library/R/arm64/4.6/library/DBI",
        "DEoptim": "/Users/z3437171/Library/R/arm64/4.6/library/DEoptim",
        "DHARMa": "/Users/z3437171/Library/R/arm64/4.6/library/DHARMa",
        "DT": "/Users/z3437171/Library/R/arm64/4.6/library/DT",
        "Deriv": "/Users/z3437171/Library/R/arm64/4.6/library/Deriv",
        "DiagrammeR": "/Users/z3437171/Library/R/arm64/4.6/library/DiagrammeR",
        "FNN": "/Users/z3437171/Library/R/arm64/4.6/library/FNN",
        "Formula": "/Users/z3437171/Library/R/arm64/4.6/library/Formula",
        "GPArotation": "/Users/z3437171/Library/R/arm64/4.6/library/GPArotation",
        "HSAUR3": "/Users/z3437171/Library/R/arm64/4.6/library/HSAUR3",
        "Hmsc": "/Users/z3437171/Library/R/arm64/4.6/library/Hmsc",
        "JuliaCall": "/Users/z3437171/Library/R/arm64/4.6/library/JuliaCall",
        "MCMCglmm": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCglmm",
        "MCMCpack": "/Users/z3437171/Library/R/arm64/4.6/library/MCMCpack",
        "MatrixModels": "/Users/z3437171/Library/R/arm64/4.6/library/MatrixModels",
        "MuMIn": "/Users/z3437171/Library/R/arm64/4.6/library/MuMIn",
        "QGglmm": "/Users/z3437171/Library/R/arm64/4.6/library/QGglmm",
        "QuickJSR": "/Users/z3437171/Library/R/arm64/4.6/library/QuickJSR",
        "R.cache": "/Users/z3437171/Library/R/arm64/4.6/library/R.cache",
        "R.methodsS3": "/Users/z3437171/Library/R/arm64/4.6/library/R.methodsS3",
        "R.oo": "/Users/z3437171/Library/R/arm64/4.6/library/R.oo",
        "R.rsp": "/Users/z3437171/Library/R/arm64/4.6/library/R.rsp",
        "R.utils": "/Users/z3437171/Library/R/arm64/4.6/library/R.utils",
        "R6": "/Users/z3437171/Library/R/arm64/4.6/library/R6",
        "RColorBrewer": "/Users/z3437171/Library/R/arm64/4.6/library/RColorBrewer",
        "ROI": "/Users/z3437171/Library/R/arm64/4.6/library/ROI",
        "ROI.plugin.lpsolve": "/Users/z3437171/Library/R/arm64/4.6/library/ROI.plugin.lpsolve",
        "RSpectra": "/Users/z3437171/Library/R/arm64/4.6/library/RSpectra",
        "RTMB": "/Users/z3437171/Library/R/arm64/4.6/library/RTMB",
        "RTMBdist": "/Users/z3437171/Library/R/arm64/4.6/library/RTMBdist",
        "RWiener": "/Users/z3437171/Library/R/arm64/4.6/library/RWiener",
        "Rcpp": "/Users/z3437171/Library/R/arm64/4.6/library/Rcpp",
        "RcppArmadillo": "/Users/z3437171/Library/R/arm64/4.6/library/RcppArmadillo",
        "RcppEigen": "/Users/z3437171/Library/R/arm64/4.6/library/RcppEigen",
        "RcppParallel": "/Users/z3437171/Library/R/arm64/4.6/library/RcppParallel",
        "RcppProgress": "/Users/z3437171/Library/R/arm64/4.6/library/RcppProgress",
        "Rdpack": "/Users/z3437171/Library/R/arm64/4.6/library/Rdpack",
        "Rfast": "/Users/z3437171/Library/R/arm64/4.6/library/Rfast",
        "Rmpfr": "/Users/z3437171/Library/R/arm64/4.6/library/Rmpfr",
        "Rphylopars": "/Users/z3437171/Library/R/arm64/4.6/library/Rphylopars",
        "S7": "/Users/z3437171/Library/R/arm64/4.6/library/S7",
        "SimDesign": "/Users/z3437171/Library/R/arm64/4.6/library/SimDesign",
        "SparseM": "/Users/z3437171/Library/R/arm64/4.6/library/SparseM",
        "StanHeaders": "/Users/z3437171/Library/R/arm64/4.6/library/StanHeaders",
        "TH.data": "/Users/z3437171/Library/R/arm64/4.6/library/TH.data",
        "TMB": "/Users/z3437171/Library/R/arm64/4.6/library/TMB",
        "V8": "/Users/z3437171/Library/R/arm64/4.6/library/V8",
        "VGAM": "/Users/z3437171/Library/R/arm64/4.6/library/VGAM",
        "XML": "/Users/z3437171/Library/R/arm64/4.6/library/XML",
        "_cache": "/Users/z3437171/Library/R/arm64/4.6/library/_cache",
        "abind": "/Users/z3437171/Library/R/arm64/4.6/library/abind",
        "ade4": "/Users/z3437171/Library/R/arm64/4.6/library/ade4",
        "alabama": "/Users/z3437171/Library/R/arm64/4.6/library/alabama",
        "ape": "/Users/z3437171/Library/R/arm64/4.6/library/ape",
        "arm": "/Users/z3437171/Library/R/arm64/4.6/library/arm",
        "askpass": "/Users/z3437171/Library/R/arm64/4.6/library/askpass",
        "assertthat": "/Users/z3437171/Library/R/arm64/4.6/library/assertthat",
        "audio": "/Users/z3437171/Library/R/arm64/4.6/library/audio",
        "backports": "/Users/z3437171/Library/R/arm64/4.6/library/backports",
        "base64enc": "/Users/z3437171/Library/R/arm64/4.6/library/base64enc",
        "bayeslincom": "/Users/z3437171/Library/R/arm64/4.6/library/bayeslincom",
        "bayesplot": "/Users/z3437171/Library/R/arm64/4.6/library/bayesplot",
        "bayestestR": "/Users/z3437171/Library/R/arm64/4.6/library/bayestestR",
        "bbmle": "/Users/z3437171/Library/R/arm64/4.6/library/bbmle",
        "bdsmatrix": "/Users/z3437171/Library/R/arm64/4.6/library/bdsmatrix",
        "beepr": "/Users/z3437171/Library/R/arm64/4.6/library/beepr",
        "beeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/beeswarm",
        "betareg": "/Users/z3437171/Library/R/arm64/4.6/library/betareg",
        "bit": "/Users/z3437171/Library/R/arm64/4.6/library/bit",
        "bit64": "/Users/z3437171/Library/R/arm64/4.6/library/bit64",
        "blob": "/Users/z3437171/Library/R/arm64/4.6/library/blob",
        "blsmeta": "/Users/z3437171/Library/R/arm64/4.6/library/blsmeta",
        "brew": "/Users/z3437171/Library/R/arm64/4.6/library/brew",
        "brglm2": "/Users/z3437171/Library/R/arm64/4.6/library/brglm2",
        "bridgesampling": "/Users/z3437171/Library/R/arm64/4.6/library/bridgesampling",
        "brio": "/Users/z3437171/Library/R/arm64/4.6/library/brio",
        "brms": "/Users/z3437171/Library/R/arm64/4.6/library/brms",
        "broom": "/Users/z3437171/Library/R/arm64/4.6/library/broom",
        "broom.mixed": "/Users/z3437171/Library/R/arm64/4.6/library/broom.mixed",
        "bslib": "/Users/z3437171/Library/R/arm64/4.6/library/bslib",
        "cachem": "/Users/z3437171/Library/R/arm64/4.6/library/cachem",
        "callr": "/Users/z3437171/Library/R/arm64/4.6/library/callr",
        "car": "/Users/z3437171/Library/R/arm64/4.6/library/car",
        "carData": "/Users/z3437171/Library/R/arm64/4.6/library/carData",
        "castor": "/Users/z3437171/Library/R/arm64/4.6/library/castor",
        "cellranger": "/Users/z3437171/Library/R/arm64/4.6/library/cellranger",
        "checkmate": "/Users/z3437171/Library/R/arm64/4.6/library/checkmate",
        "circular": "/Users/z3437171/Library/R/arm64/4.6/library/circular",
        "classInt": "/Users/z3437171/Library/R/arm64/4.6/library/classInt",
        "cli": "/Users/z3437171/Library/R/arm64/4.6/library/cli",
        "clipr": "/Users/z3437171/Library/R/arm64/4.6/library/clipr",
        "clue": "/Users/z3437171/Library/R/arm64/4.6/library/clue",
        "clusterGeneration": "/Users/z3437171/Library/R/arm64/4.6/library/clusterGeneration",
        "cmdstanr": "/Users/z3437171/Library/R/arm64/4.6/library/cmdstanr",
        "coda": "/Users/z3437171/Library/R/arm64/4.6/library/coda",
        "collections": "/Users/z3437171/Library/R/arm64/4.6/library/collections",
        "colorspace": "/Users/z3437171/Library/R/arm64/4.6/library/colorspace",
        "colourpicker": "/Users/z3437171/Library/R/arm64/4.6/library/colourpicker",
        "combinat": "/Users/z3437171/Library/R/arm64/4.6/library/combinat",
        "commonmark": "/Users/z3437171/Library/R/arm64/4.6/library/commonmark",
        "conflicted": "/Users/z3437171/Library/R/arm64/4.6/library/conflicted",
        "coro": "/Users/z3437171/Library/R/arm64/4.6/library/coro",
        "corpcor": "/Users/z3437171/Library/R/arm64/4.6/library/corpcor",
        "covr": "/Users/z3437171/Library/R/arm64/4.6/library/covr",
        "cowplot": "/Users/z3437171/Library/R/arm64/4.6/library/cowplot",
        "cpp11": "/Users/z3437171/Library/R/arm64/4.6/library/cpp11",
        "crayon": "/Users/z3437171/Library/R/arm64/4.6/library/crayon",
        "credentials": "/Users/z3437171/Library/R/arm64/4.6/library/credentials",
        "crosstalk": "/Users/z3437171/Library/R/arm64/4.6/library/crosstalk",
        "crul": "/Users/z3437171/Library/R/arm64/4.6/library/crul",
        "cubature": "/Users/z3437171/Library/R/arm64/4.6/library/cubature",
        "curl": "/Users/z3437171/Library/R/arm64/4.6/library/curl",
        "dagitty": "/Users/z3437171/Library/R/arm64/4.6/library/dagitty",
        "data.table": "/Users/z3437171/Library/R/arm64/4.6/library/data.table",
        "datawizard": "/Users/z3437171/Library/R/arm64/4.6/library/datawizard",
        "dbplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dbplyr",
        "dcurver": "/Users/z3437171/Library/R/arm64/4.6/library/dcurver",
        "deldir": "/Users/z3437171/Library/R/arm64/4.6/library/deldir",
        "desc": "/Users/z3437171/Library/R/arm64/4.6/library/desc",
        "detectseparation": "/Users/z3437171/Library/R/arm64/4.6/library/detectseparation",
        "devtools": "/Users/z3437171/Library/R/arm64/4.6/library/devtools",
        "diffobj": "/Users/z3437171/Library/R/arm64/4.6/library/diffobj",
        "digest": "/Users/z3437171/Library/R/arm64/4.6/library/digest",
        "distributional": "/Users/z3437171/Library/R/arm64/4.6/library/distributional",
        "doBy": "/Users/z3437171/Library/R/arm64/4.6/library/doBy",
        "doFuture": "/Users/z3437171/Library/R/arm64/4.6/library/doFuture",
        "doParallel": "/Users/z3437171/Library/R/arm64/4.6/library/doParallel",
        "doRNG": "/Users/z3437171/Library/R/arm64/4.6/library/doRNG",
        "dotCall64": "/Users/z3437171/Library/R/arm64/4.6/library/dotCall64",
        "downlit": "/Users/z3437171/Library/R/arm64/4.6/library/downlit",
        "dplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dplyr",
        "drmSEM": "/Users/z3437171/Library/R/arm64/4.6/library/drmSEM",
        "drmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/drmTMB",
        "dtplyr": "/Users/z3437171/Library/R/arm64/4.6/library/dtplyr",
        "dygraphs": "/Users/z3437171/Library/R/arm64/4.6/library/dygraphs",
        "e1071": "/Users/z3437171/Library/R/arm64/4.6/library/e1071",
        "ellipsis": "/Users/z3437171/Library/R/arm64/4.6/library/ellipsis",
        "emmeans": "/Users/z3437171/Library/R/arm64/4.6/library/emmeans",
        "enhancer": "/Users/z3437171/Library/R/arm64/4.6/library/enhancer",
        "enrichwith": "/Users/z3437171/Library/R/arm64/4.6/library/enrichwith",
        "estimability": "/Users/z3437171/Library/R/arm64/4.6/library/estimability",
        "evaluate": "/Users/z3437171/Library/R/arm64/4.6/library/evaluate",
        "evd": "/Users/z3437171/Library/R/arm64/4.6/library/evd",
        "expm": "/Users/z3437171/Library/R/arm64/4.6/library/expm",
        "extraDistr": "/Users/z3437171/Library/R/arm64/4.6/library/extraDistr",
        "fansi": "/Users/z3437171/Library/R/arm64/4.6/library/fansi",
        "farver": "/Users/z3437171/Library/R/arm64/4.6/library/farver",
        "fastGHQuad": "/Users/z3437171/Library/R/arm64/4.6/library/fastGHQuad",
        "fastmap": "/Users/z3437171/Library/R/arm64/4.6/library/fastmap",
        "fastmatch": "/Users/z3437171/Library/R/arm64/4.6/library/fastmatch",
        "fields": "/Users/z3437171/Library/R/arm64/4.6/library/fields",
        "fishMod": "/Users/z3437171/Library/R/arm64/4.6/library/fishMod",
        "flexmix": "/Users/z3437171/Library/R/arm64/4.6/library/flexmix",
        "fmesher": "/Users/z3437171/Library/R/arm64/4.6/library/fmesher",
        "fontawesome": "/Users/z3437171/Library/R/arm64/4.6/library/fontawesome",
        "forcats": "/Users/z3437171/Library/R/arm64/4.6/library/forcats",
        "foreach": "/Users/z3437171/Library/R/arm64/4.6/library/foreach",
        "forecast": "/Users/z3437171/Library/R/arm64/4.6/library/forecast",
        "fracdiff": "/Users/z3437171/Library/R/arm64/4.6/library/fracdiff",
        "freqTLS": "/Users/z3437171/Library/R/arm64/4.6/library/freqTLS",
        "fs": "/Users/z3437171/Library/R/arm64/4.6/library/fs",
        "furrr": "/Users/z3437171/Library/R/arm64/4.6/library/furrr",
        "future": "/Users/z3437171/Library/R/arm64/4.6/library/future",
        "future.apply": "/Users/z3437171/Library/R/arm64/4.6/library/future.apply",
        "future.mirai": "/Users/z3437171/Library/R/arm64/4.6/library/future.mirai",
        "galamm": "/Users/z3437171/Library/R/arm64/4.6/library/galamm",
        "gamlss.dist": "/Users/z3437171/Library/R/arm64/4.6/library/gamlss.dist",
        "gamm4": "/Users/z3437171/Library/R/arm64/4.6/library/gamm4",
        "gap": "/Users/z3437171/Library/R/arm64/4.6/library/gap",
        "gap.datasets": "/Users/z3437171/Library/R/arm64/4.6/library/gap.datasets",
        "gargle": "/Users/z3437171/Library/R/arm64/4.6/library/gargle",
        "generics": "/Users/z3437171/Library/R/arm64/4.6/library/generics",
        "geojsonsf": "/Users/z3437171/Library/R/arm64/4.6/library/geojsonsf",
        "geometries": "/Users/z3437171/Library/R/arm64/4.6/library/geometries",
        "gert": "/Users/z3437171/Library/R/arm64/4.6/library/gert",
        "ggbeeswarm": "/Users/z3437171/Library/R/arm64/4.6/library/ggbeeswarm",
        "ggdist": "/Users/z3437171/Library/R/arm64/4.6/library/ggdist",
        "ggfittext": "/Users/z3437171/Library/R/arm64/4.6/library/ggfittext",
        "ggforce": "/Users/z3437171/Library/R/arm64/4.6/library/ggforce",
        "gggenes": "/Users/z3437171/Library/R/arm64/4.6/library/gggenes",
        "ggh4x": "/Users/z3437171/Library/R/arm64/4.6/library/ggh4x",
        "ggm": "/Users/z3437171/Library/R/arm64/4.6/library/ggm",
        "ggokabeito": "/Users/z3437171/Library/R/arm64/4.6/library/ggokabeito",
        "ggplot2": "/Users/z3437171/Library/R/arm64/4.6/library/ggplot2",
        "ggraph": "/Users/z3437171/Library/R/arm64/4.6/library/ggraph",
        "ggrepel": "/Users/z3437171/Library/R/arm64/4.6/library/ggrepel",
        "ggridges": "/Users/z3437171/Library/R/arm64/4.6/library/ggridges",
        "gh": "/Users/z3437171/Library/R/arm64/4.6/library/gh",
        "gitcreds": "/Users/z3437171/Library/R/arm64/4.6/library/gitcreds",
        "glassoFast": "/Users/z3437171/Library/R/arm64/4.6/library/glassoFast",
        "gllvm": "/Users/z3437171/Library/R/arm64/4.6/library/gllvm",
        "gllvmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/gllvmTMB",
        "glmmTMB": "/Users/z3437171/Library/R/arm64/4.6/library/glmmTMB",
        "glmnet": "/Users/z3437171/Library/R/arm64/4.6/library/glmnet",
        "globals": "/Users/z3437171/Library/R/arm64/4.6/library/globals",
        "glue": "/Users/z3437171/Library/R/arm64/4.6/library/glue",
        "gmp": "/Users/z3437171/Library/R/arm64/4.6/library/gmp",
        "googledrive": "/Users/z3437171/Library/R/arm64/4.6/library/googledrive",
        "googlesheets4": "/Users/z3437171/Library/R/arm64/4.6/library/googlesheets4",
        "graph": "/Users/z3437171/Library/R/arm64/4.6/library/graph",
        "graphlayouts": "/Users/z3437171/Library/R/arm64/4.6/library/graphlayouts",
        "gratia": "/Users/z3437171/Library/R/arm64/4.6/library/gratia",
        "gridExtra": "/Users/z3437171/Library/R/arm64/4.6/library/gridExtra",
        "gridtext": "/Users/z3437171/Library/R/arm64/4.6/library/gridtext",
        "gsl": "/Users/z3437171/Library/R/arm64/4.6/library/gsl",
        "gtable": "/Users/z3437171/Library/R/arm64/4.6/library/gtable",
        "gtools": "/Users/z3437171/Library/R/arm64/4.6/library/gtools",
        "haven": "/Users/z3437171/Library/R/arm64/4.6/library/haven",
        "here": "/Users/z3437171/Library/R/arm64/4.6/library/here",
        "highr": "/Users/z3437171/Library/R/arm64/4.6/library/highr",
        "hms": "/Users/z3437171/Library/R/arm64/4.6/library/hms",
        "hsquared": "/Users/z3437171/Library/R/arm64/4.6/library/hsquared",
        "htmltools": "/Users/z3437171/Library/R/arm64/4.6/library/htmltools",
        "htmlwidgets": "/Users/z3437171/Library/R/arm64/4.6/library/htmlwidgets",
        "httpcode": "/Users/z3437171/Library/R/arm64/4.6/library/httpcode",
        "httpuv": "/Users/z3437171/Library/R/arm64/4.6/library/httpuv",
        "httr": "/Users/z3437171/Library/R/arm64/4.6/library/httr",
        "httr2": "/Users/z3437171/Library/R/arm64/4.6/library/httr2",
        "hunspell": "/Users/z3437171/Library/R/arm64/4.6/library/hunspell",
        "ids": "/Users/z3437171/Library/R/arm64/4.6/library/ids",
        "igraph": "/Users/z3437171/Library/R/arm64/4.6/library/igraph",
        "ini": "/Users/z3437171/Library/R/arm64/4.6/library/ini",
        "inline": "/Users/z3437171/Library/R/arm64/4.6/library/inline",
        "insight": "/Users/z3437171/Library/R/arm64/4.6/library/insight",
        "isoband": "/Users/z3437171/Library/R/arm64/4.6/library/isoband",
        "iterators": "/Users/z3437171/Library/R/arm64/4.6/library/iterators",
        "itertools": "/Users/z3437171/Library/R/arm64/4.6/library/itertools",
        "janitor": "/Users/z3437171/Library/R/arm64/4.6/library/janitor",
        "jomo": "/Users/z3437171/Library/R/arm64/4.6/library/jomo",
        "jpeg": "/Users/z3437171/Library/R/arm64/4.6/library/jpeg",
        "jquerylib": "/Users/z3437171/Library/R/arm64/4.6/library/jquerylib",
        "jsonify": "/Users/z3437171/Library/R/arm64/4.6/library/jsonify",
        "jsonlite": "/Users/z3437171/Library/R/arm64/4.6/library/jsonlite",
        "kernelshap": "/Users/z3437171/Library/R/arm64/4.6/library/kernelshap",
        "knitr": "/Users/z3437171/Library/R/arm64/4.6/library/knitr",
        "labeling": "/Users/z3437171/Library/R/arm64/4.6/library/labeling",
        "languageserver": "/Users/z3437171/Library/R/arm64/4.6/library/languageserver",
        "later": "/Users/z3437171/Library/R/arm64/4.6/library/later",
        "latex2exp": "/Users/z3437171/Library/R/arm64/4.6/library/latex2exp",
        "lavaan": "/Users/z3437171/Library/R/arm64/4.6/library/lavaan",
        "lazyeval": "/Users/z3437171/Library/R/arm64/4.6/library/lazyeval",
        "leafem": "/Users/z3437171/Library/R/arm64/4.6/library/leafem",
        "leaflet": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet",
        "leaflet.providers": "/Users/z3437171/Library/R/arm64/4.6/library/leaflet.providers",
        "leafpop": "/Users/z3437171/Library/R/arm64/4.6/library/leafpop",
        "lifecycle": "/Users/z3437171/Library/R/arm64/4.6/library/lifecycle",
        "lintr": "/Users/z3437171/Library/R/arm64/4.6/library/lintr",
        "listenv": "/Users/z3437171/Library/R/arm64/4.6/library/listenv",
        "litedown": "/Users/z3437171/Library/R/arm64/4.6/library/litedown",
        "lme4": "/Users/z3437171/Library/R/arm64/4.6/library/lme4",
        "lmtest": "/Users/z3437171/Library/R/arm64/4.6/library/lmtest",
        "loo": "/Users/z3437171/Library/R/arm64/4.6/library/loo",
        "lpSolveAPI": "/Users/z3437171/Library/R/arm64/4.6/library/lpSolveAPI",
        "lubridate": "/Users/z3437171/Library/R/arm64/4.6/library/lubridate",
        "magrittr": "/Users/z3437171/Library/R/arm64/4.6/library/magrittr",
        "maps": "/Users/z3437171/Library/R/arm64/4.6/library/maps",
        "mapview": "/Users/z3437171/Library/R/arm64/4.6/library/mapview",
        "markdown": "/Users/z3437171/Library/R/arm64/4.6/library/markdown",
        "mathjaxr": "/Users/z3437171/Library/R/arm64/4.6/library/mathjaxr",
        "matrixStats": "/Users/z3437171/Library/R/arm64/4.6/library/matrixStats",
        "mclogit": "/Users/z3437171/Library/R/arm64/4.6/library/mclogit",
        "mcmc": "/Users/z3437171/Library/R/arm64/4.6/library/mcmc",
        "memisc": "/Users/z3437171/Library/R/arm64/4.6/library/memisc",
        "memoise": "/Users/z3437171/Library/R/arm64/4.6/library/memoise",
        "metadat": "/Users/z3437171/Library/R/arm64/4.6/library/metadat",
        "metafor": "/Users/z3437171/Library/R/arm64/4.6/library/metafor",
        "metaplus": "/Users/z3437171/Library/R/arm64/4.6/library/metaplus",
        "mice": "/Users/z3437171/Library/R/arm64/4.6/library/mice",
        "microbenchmark": "/Users/z3437171/Library/R/arm64/4.6/library/microbenchmark",
        "mime": "/Users/z3437171/Library/R/arm64/4.6/library/mime",
        "miniUI": "/Users/z3437171/Library/R/arm64/4.6/library/miniUI",
        "minqa": "/Users/z3437171/Library/R/arm64/4.6/library/minqa",
        "mirai": "/Users/z3437171/Library/R/arm64/4.6/library/mirai",
        "mirt": "/Users/z3437171/Library/R/arm64/4.6/library/mirt",
        "missForest": "/Users/z3437171/Library/R/arm64/4.6/library/missForest",
        "mitml": "/Users/z3437171/Library/R/arm64/4.6/library/mitml",
        "mnormt": "/Users/z3437171/Library/R/arm64/4.6/library/mnormt",
        "modelr": "/Users/z3437171/Library/R/arm64/4.6/library/modelr",
        "modeltools": "/Users/z3437171/Library/R/arm64/4.6/library/modeltools",
        "movMF": "/Users/z3437171/Library/R/arm64/4.6/library/movMF",
        "msm": "/Users/z3437171/Library/R/arm64/4.6/library/msm",
        "multcomp": "/Users/z3437171/Library/R/arm64/4.6/library/multcomp",
        "mvMORPH": "/Users/z3437171/Library/R/arm64/4.6/library/mvMORPH",
        "mvnfast": "/Users/z3437171/Library/R/arm64/4.6/library/mvnfast",
        "mvtnorm": "/Users/z3437171/Library/R/arm64/4.6/library/mvtnorm",
        "nadiv": "/Users/z3437171/Library/R/arm64/4.6/library/nadiv",
        "nanonext": "/Users/z3437171/Library/R/arm64/4.6/library/nanonext",
        "naturalsort": "/Users/z3437171/Library/R/arm64/4.6/library/naturalsort",
        "ncf": "/Users/z3437171/Library/R/arm64/4.6/library/ncf",
        "nleqslv": "/Users/z3437171/Library/R/arm64/4.6/library/nleqslv",
        "nloptr": "/Users/z3437171/Library/R/arm64/4.6/library/nloptr",
        "numDeriv": "/Users/z3437171/Library/R/arm64/4.6/library/numDeriv",
        "oai": "/Users/z3437171/Library/R/arm64/4.6/library/oai",
        "openssl": "/Users/z3437171/Library/R/arm64/4.6/library/openssl",
        "optimParallel": "/Users/z3437171/Library/R/arm64/4.6/library/optimParallel",
        "orchaRd": "/Users/z3437171/Library/R/arm64/4.6/library/orchaRd",
        "ordinal": "/Users/z3437171/Library/R/arm64/4.6/library/ordinal",
        "orthopolynom": "/Users/z3437171/Library/R/arm64/4.6/library/orthopolynom",
        "otel": "/Users/z3437171/Library/R/arm64/4.6/library/otel",
        "pROC": "/Users/z3437171/Library/R/arm64/4.6/library/pROC",
        "pak": "/Users/z3437171/Library/R/arm64/4.6/library/pak",
        "palmerpenguins": "/Users/z3437171/Library/R/arm64/4.6/library/palmerpenguins",
        "pan": "/Users/z3437171/Library/R/arm64/4.6/library/pan",
        "parallelly": "/Users/z3437171/Library/R/arm64/4.6/library/parallelly",
        "patchwork": "/Users/z3437171/Library/R/arm64/4.6/library/patchwork",
        "pbapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbapply",
        "pbivnorm": "/Users/z3437171/Library/R/arm64/4.6/library/pbivnorm",
        "pbkrtest": "/Users/z3437171/Library/R/arm64/4.6/library/pbkrtest",
        "pbmcapply": "/Users/z3437171/Library/R/arm64/4.6/library/pbmcapply",
        "pedigreemm": "/Users/z3437171/Library/R/arm64/4.6/library/pedigreemm",
        "performance": "/Users/z3437171/Library/R/arm64/4.6/library/performance",
        "permute": "/Users/z3437171/Library/R/arm64/4.6/library/permute",
        "phangorn": "/Users/z3437171/Library/R/arm64/4.6/library/phangorn",
        "phylolm": "/Users/z3437171/Library/R/arm64/4.6/library/phylolm",
        "phylopath": "/Users/z3437171/Library/R/arm64/4.6/library/phylopath",
        "phytools": "/Users/z3437171/Library/R/arm64/4.6/library/phytools",
        "piecewiseSEM": "/Users/z3437171/Library/R/arm64/4.6/library/piecewiseSEM",
        "pillar": "/Users/z3437171/Library/R/arm64/4.6/library/pillar",
        "pixmap": "/Users/z3437171/Library/R/arm64/4.6/library/pixmap",
        "pkgbuild": "/Users/z3437171/Library/R/arm64/4.6/library/pkgbuild",
        "pkgconfig": "/Users/z3437171/Library/R/arm64/4.6/library/pkgconfig",
        "pkgdown": "/Users/z3437171/Library/R/arm64/4.6/library/pkgdown",
        "pkgload": "/Users/z3437171/Library/R/arm64/4.6/library/pkgload",
        "plotly": "/Users/z3437171/Library/R/arm64/4.6/library/plotly",
        "plyr": "/Users/z3437171/Library/R/arm64/4.6/library/plyr",
        "png": "/Users/z3437171/Library/R/arm64/4.6/library/png",
        "polyclip": "/Users/z3437171/Library/R/arm64/4.6/library/polyclip",
        "polynom": "/Users/z3437171/Library/R/arm64/4.6/library/polynom",
        "posterior": "/Users/z3437171/Library/R/arm64/4.6/library/posterior",
        "pracma": "/Users/z3437171/Library/R/arm64/4.6/library/pracma",
        "praise": "/Users/z3437171/Library/R/arm64/4.6/library/praise",
        "prettyunits": "/Users/z3437171/Library/R/arm64/4.6/library/prettyunits",
        "priorsense": "/Users/z3437171/Library/R/arm64/4.6/library/priorsense",
        "processx": "/Users/z3437171/Library/R/arm64/4.6/library/processx",
        "profvis": "/Users/z3437171/Library/R/arm64/4.6/library/profvis",
        "progress": "/Users/z3437171/Library/R/arm64/4.6/library/progress",
        "progressr": "/Users/z3437171/Library/R/arm64/4.6/library/progressr",
        "projpred": "/Users/z3437171/Library/R/arm64/4.6/library/projpred",
        "promises": "/Users/z3437171/Library/R/arm64/4.6/library/promises",
        "proxy": "/Users/z3437171/Library/R/arm64/4.6/library/proxy",
        "ps": "/Users/z3437171/Library/R/arm64/4.6/library/ps",
        "psymetadata": "/Users/z3437171/Library/R/arm64/4.6/library/psymetadata",
        "purrr": "/Users/z3437171/Library/R/arm64/4.6/library/purrr",
        "qgam": "/Users/z3437171/Library/R/arm64/4.6/library/qgam",
        "qs2": "/Users/z3437171/Library/R/arm64/4.6/library/qs2",
        "quadprog": "/Users/z3437171/Library/R/arm64/4.6/library/quadprog",
        "quantreg": "/Users/z3437171/Library/R/arm64/4.6/library/quantreg",
        "rWishart": "/Users/z3437171/Library/R/arm64/4.6/library/rWishart",
        "ragg": "/Users/z3437171/Library/R/arm64/4.6/library/ragg",
        "randomForest": "/Users/z3437171/Library/R/arm64/4.6/library/randomForest",
        "ranger": "/Users/z3437171/Library/R/arm64/4.6/library/ranger",
        "rapidjsonr": "/Users/z3437171/Library/R/arm64/4.6/library/rapidjsonr",
        "rappdirs": "/Users/z3437171/Library/R/arm64/4.6/library/rappdirs",
        "raster": "/Users/z3437171/Library/R/arm64/4.6/library/raster",
        "rbibutils": "/Users/z3437171/Library/R/arm64/4.6/library/rbibutils",
        "rcmdcheck": "/Users/z3437171/Library/R/arm64/4.6/library/rcmdcheck",
        "readr": "/Users/z3437171/Library/R/arm64/4.6/library/readr",
        "readxl": "/Users/z3437171/Library/R/arm64/4.6/library/readxl",
        "reformulas": "/Users/z3437171/Library/R/arm64/4.6/library/reformulas",
        "registry": "/Users/z3437171/Library/R/arm64/4.6/library/registry",
        "rematch": "/Users/z3437171/Library/R/arm64/4.6/library/rematch",
        "rematch2": "/Users/z3437171/Library/R/arm64/4.6/library/rematch2",
        "remotes": "/Users/z3437171/Library/R/arm64/4.6/library/remotes",
        "rentrez": "/Users/z3437171/Library/R/arm64/4.6/library/rentrez",
        "reprex": "/Users/z3437171/Library/R/arm64/4.6/library/reprex",
        "reshape2": "/Users/z3437171/Library/R/arm64/4.6/library/reshape2",
        "rex": "/Users/z3437171/Library/R/arm64/4.6/library/rex",
        "rgbif": "/Users/z3437171/Library/R/arm64/4.6/library/rgbif",
        "rgpdd": "/Users/z3437171/Library/R/arm64/4.6/library/rgpdd",
        "rhub": "/Users/z3437171/Library/R/arm64/4.6/library/rhub",
        "rjags": "/Users/z3437171/Library/R/arm64/4.6/library/rjags",
        "rjson": "/Users/z3437171/Library/R/arm64/4.6/library/rjson",
        "rlang": "/Users/z3437171/Library/R/arm64/4.6/library/rlang",
        "rmarkdown": "/Users/z3437171/Library/R/arm64/4.6/library/rmarkdown",
        "rncl": "/Users/z3437171/Library/R/arm64/4.6/library/rncl",
        "rngtools": "/Users/z3437171/Library/R/arm64/4.6/library/rngtools",
        "rotl": "/Users/z3437171/Library/R/arm64/4.6/library/rotl",
        "roxygen2": "/Users/z3437171/Library/R/arm64/4.6/library/roxygen2",
        "rprojroot": "/Users/z3437171/Library/R/arm64/4.6/library/rprojroot",
        "rstan": "/Users/z3437171/Library/R/arm64/4.6/library/rstan",
        "rstantools": "/Users/z3437171/Library/R/arm64/4.6/library/rstantools",
        "rstudioapi": "/Users/z3437171/Library/R/arm64/4.6/library/rstudioapi",
        "rtdists": "/Users/z3437171/Library/R/arm64/4.6/library/rtdists",
        "rversions": "/Users/z3437171/Library/R/arm64/4.6/library/rversions",
        "rvest": "/Users/z3437171/Library/R/arm64/4.6/library/rvest",
        "s2": "/Users/z3437171/Library/R/arm64/4.6/library/s2",
        "safetensors": "/Users/z3437171/Library/R/arm64/4.6/library/safetensors",
        "sandwich": "/Users/z3437171/Library/R/arm64/4.6/library/sandwich",
        "sass": "/Users/z3437171/Library/R/arm64/4.6/library/sass",
        "satellite": "/Users/z3437171/Library/R/arm64/4.6/library/satellite",
        "scales": "/Users/z3437171/Library/R/arm64/4.6/library/scales",
        "scatterplot3d": "/Users/z3437171/Library/R/arm64/4.6/library/scatterplot3d",
        "selectr": "/Users/z3437171/Library/R/arm64/4.6/library/selectr",
        "servr": "/Users/z3437171/Library/R/arm64/4.6/library/servr",
        "sessioninfo": "/Users/z3437171/Library/R/arm64/4.6/library/sessioninfo",
        "sf": "/Users/z3437171/Library/R/arm64/4.6/library/sf",
        "sfheaders": "/Users/z3437171/Library/R/arm64/4.6/library/sfheaders",
        "shades": "/Users/z3437171/Library/R/arm64/4.6/library/shades",
        "shape": "/Users/z3437171/Library/R/arm64/4.6/library/shape",
        "shapviz": "/Users/z3437171/Library/R/arm64/4.6/library/shapviz",
        "shiny": "/Users/z3437171/Library/R/arm64/4.6/library/shiny",
        "shinyjs": "/Users/z3437171/Library/R/arm64/4.6/library/shinyjs",
        "shinystan": "/Users/z3437171/Library/R/arm64/4.6/library/shinystan",
        "shinythemes": "/Users/z3437171/Library/R/arm64/4.6/library/shinythemes",
        "skmeans": "/Users/z3437171/Library/R/arm64/4.6/library/skmeans",
        "slam": "/Users/z3437171/Library/R/arm64/4.6/library/slam",
        "smcfcs": "/Users/z3437171/Library/R/arm64/4.6/library/smcfcs",
        "sn": "/Users/z3437171/Library/R/arm64/4.6/library/sn",
        "snakecase": "/Users/z3437171/Library/R/arm64/4.6/library/snakecase",
        "sommer": "/Users/z3437171/Library/R/arm64/4.6/library/sommer",
        "sourcetools": "/Users/z3437171/Library/R/arm64/4.6/library/sourcetools",
        "sp": "/Users/z3437171/Library/R/arm64/4.6/library/sp",
        "spData": "/Users/z3437171/Library/R/arm64/4.6/library/spData",
        "spam": "/Users/z3437171/Library/R/arm64/4.6/library/spam",
        "spdep": "/Users/z3437171/Library/R/arm64/4.6/library/spdep",
        "spelling": "/Users/z3437171/Library/R/arm64/4.6/library/spelling",
        "splancs": "/Users/z3437171/Library/R/arm64/4.6/library/splancs",
        "splines2": "/Users/z3437171/Library/R/arm64/4.6/library/splines2",
        "statmod": "/Users/z3437171/Library/R/arm64/4.6/library/statmod",
        "stringfish": "/Users/z3437171/Library/R/arm64/4.6/library/stringfish",
        "stringi": "/Users/z3437171/Library/R/arm64/4.6/library/stringi",
        "stringr": "/Users/z3437171/Library/R/arm64/4.6/library/stringr",
        "styler": "/Users/z3437171/Library/R/arm64/4.6/library/styler",
        "subplex": "/Users/z3437171/Library/R/arm64/4.6/library/subplex",
        "svglite": "/Users/z3437171/Library/R/arm64/4.6/library/svglite",
        "sys": "/Users/z3437171/Library/R/arm64/4.6/library/sys",
        "systemfonts": "/Users/z3437171/Library/R/arm64/4.6/library/systemfonts",
        "tensorA": "/Users/z3437171/Library/R/arm64/4.6/library/tensorA",
        "terra": "/Users/z3437171/Library/R/arm64/4.6/library/terra",
        "tester": "/Users/z3437171/Library/R/arm64/4.6/library/tester",
        "testthat": "/Users/z3437171/Library/R/arm64/4.6/library/testthat",
        "textshaping": "/Users/z3437171/Library/R/arm64/4.6/library/textshaping",
        "threejs": "/Users/z3437171/Library/R/arm64/4.6/library/threejs",
        "tibble": "/Users/z3437171/Library/R/arm64/4.6/library/tibble",
        "tidygraph": "/Users/z3437171/Library/R/arm64/4.6/library/tidygraph",
        "tidyr": "/Users/z3437171/Library/R/arm64/4.6/library/tidyr",
        "tidyselect": "/Users/z3437171/Library/R/arm64/4.6/library/tidyselect",
        "tidyverse": "/Users/z3437171/Library/R/arm64/4.6/library/tidyverse",
        "timeDate": "/Users/z3437171/Library/R/arm64/4.6/library/timeDate",
        "timechange": "/Users/z3437171/Library/R/arm64/4.6/library/timechange",
        "tinytex": "/Users/z3437171/Library/R/arm64/4.6/library/tinytex",
        "torch": "/Users/z3437171/Library/R/arm64/4.6/library/torch",
        "triebeard": "/Users/z3437171/Library/R/arm64/4.6/library/triebeard",
        "truncnorm": "/Users/z3437171/Library/R/arm64/4.6/library/truncnorm",
        "tweedie": "/Users/z3437171/Library/R/arm64/4.6/library/tweedie",
        "tweenr": "/Users/z3437171/Library/R/arm64/4.6/library/tweenr",
        "tzdb": "/Users/z3437171/Library/R/arm64/4.6/library/tzdb",
        "ucminf": "/Users/z3437171/Library/R/arm64/4.6/library/ucminf",
        "units": "/Users/z3437171/Library/R/arm64/4.6/library/units",
        "urca": "/Users/z3437171/Library/R/arm64/4.6/library/urca",
        "urlchecker": "/Users/z3437171/Library/R/arm64/4.6/library/urlchecker",
        "urltools": "/Users/z3437171/Library/R/arm64/4.6/library/urltools",
        "usethis": "/Users/z3437171/Library/R/arm64/4.6/library/usethis",
        "utf8": "/Users/z3437171/Library/R/arm64/4.6/library/utf8",
        "uuid": "/Users/z3437171/Library/R/arm64/4.6/library/uuid",
        "vctrs": "/Users/z3437171/Library/R/arm64/4.6/library/vctrs",
        "vdiffr": "/Users/z3437171/Library/R/arm64/4.6/library/vdiffr",
        "vegan": "/Users/z3437171/Library/R/arm64/4.6/library/vegan",
        "vipor": "/Users/z3437171/Library/R/arm64/4.6/library/vipor",
        "viridis": "/Users/z3437171/Library/R/arm64/4.6/library/viridis",
        "viridisLite": "/Users/z3437171/Library/R/arm64/4.6/library/viridisLite",
        "visNetwork": "/Users/z3437171/Library/R/arm64/4.6/library/visNetwork",
        "vroom": "/Users/z3437171/Library/R/arm64/4.6/library/vroom",
        "waldo": "/Users/z3437171/Library/R/arm64/4.6/library/waldo",
        "whisker": "/Users/z3437171/Library/R/arm64/4.6/library/whisker",
        "whoami": "/Users/z3437171/Library/R/arm64/4.6/library/whoami",
        "withr": "/Users/z3437171/Library/R/arm64/4.6/library/withr",
        "wk": "/Users/z3437171/Library/R/arm64/4.6/library/wk",
        "xfun": "/Users/z3437171/Library/R/arm64/4.6/library/xfun",
        "xgboost": "/Users/z3437171/Library/R/arm64/4.6/library/xgboost",
        "xml2": "/Users/z3437171/Library/R/arm64/4.6/library/xml2",
        "xmlparsedata": "/Users/z3437171/Library/R/arm64/4.6/library/xmlparsedata",
        "xopen": "/Users/z3437171/Library/R/arm64/4.6/library/xopen",
        "xtable": "/Users/z3437171/Library/R/arm64/4.6/library/xtable",
        "xts": "/Users/z3437171/Library/R/arm64/4.6/library/xts",
        "yaml": "/Users/z3437171/Library/R/arm64/4.6/library/yaml",
        "zigg": "/Users/z3437171/Library/R/arm64/4.6/library/zigg",
        "zip": "/Users/z3437171/Library/R/arm64/4.6/library/zip",
        "zoo": "/Users/z3437171/Library/R/arm64/4.6/library/zoo"
      },
      "steps": [
        {
          "step": "install",
          "command": [
            "R",
            "CMD",
            "INSTALL",
            "--no-multiarch",
            "-l",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/both",
            "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:01:25.178470+00:00",
          "elapsed_seconds": 21.653,
          "log": "both-install.log",
          "log_sha256": "7f8a73b86b7d03565ff3003762e097a22aabf7fce7f9737a57df3161adc95993",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "fixtures",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "fixtures",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/both"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:01:46.838964+00:00",
          "elapsed_seconds": 5.722,
          "log": "both-fixtures.log",
          "log_sha256": "cb398da193bd8c6eb8c8ecf46cdf7b115e328991a7f5aad39c66efcc529ecb4e",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        },
        {
          "step": "real",
          "command": [
            "Rscript",
            "--vanilla",
            "/private/tmp/pigauto-release-evidence-780c084/script/cran-0.11-integration/check-adapters.R",
            "real",
            "installed",
            "/private/tmp/pigauto-mi-matrix-c2e768ad/both"
          ],
          "exit_code": 0,
          "started_utc": "2026-10-17T17:01:52.584060+00:00",
          "elapsed_seconds": 18.171,
          "log": "both-real.log",
          "log_sha256": "dcba11ad9f11fe5740af2f1301bc894fb7196ac2082e2483cb189cb8612cbac9",
          "archive_sha256_before": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_sha256_after": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
          "archive_unchanged": true
        }
      ]
    }
  ],
  "status": "PASS",
  "failures": []
}
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/neither-fixtures.log`

~~~~text
tests/testthat/test-pool-mi-adapter-fixtures.R: 50 expectations passed against installed pigauto
tests/testthat/test-mi-provenance.R: 55 expectations passed against installed pigauto
ADAPTER_FIXTURES_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/neither-install.log`

~~~~text
* installing *source* package ‘pigauto’ ...
** this is package ‘pigauto’ version ‘0.11.0’
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (pigauto)
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/neither-real.log`

~~~~text
script/cran-0.11-integration/test-real-backends.R: 0 expectations passed against installed pigauto
SAVED_OPTIONAL_CLASS_MISSING_PACKAGE_OK
NEITHER_BACKEND_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/site-9a457eb-review.md`

~~~~text
# Candidate site review for source commit 9a457eb6

Date: 2026-10-17. A fresh site was built from source commit `9a457eb6b3363d413dbc502de88ff5bc602ed98a` into `/private/tmp/pigauto-cran011-site-9a457eb.qiPr9J/_site` using the release-site builder and seven verified cached assets. The build reports 62 HTML pages and 614 search entries. `pkgdown::check_pkgdown()` found no problems. The automated crawler checked 3,588 local references and 34 retired pages, with zero errors (`SITE_CRAWL_OK`). Retirement assertions passed for the retired article and search/sitemap entries.

In the Codex in-app browser, I inspected the locally rendered `simulate_benchmark` reference page. It retains the public `gnn = TRUE` default for this comparison workflow, and its `\donttest` example explicitly sets `gnn = FALSE`. The rendered example completes and prints the baseline and pigauto summaries. It also reports the expected small-validation warning for this deliberately tiny example. The browser review covered rendered content and accessibility-tree text. Pixel-level appearance was not assessed.

The live site reflects an earlier source state. Merge and deployment must precede a live route check; the deployed-site gate stays open.

Build log SHA-256: `97469ccfb06af9074366ef61b939be1c55fe1057d4a1168ba18fff774a5d7364`. Crawler log SHA-256: `1bf804fb25425e4bbe6138f82646bb8c323a79094ae352ea399f41790e765128`.
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/site-build-9a457eb.log`

~~~~text
── Installing package pigauto into temporary library ───────────────────────────
── Initialising site ───────────────────────────────────────────────────────────
Copying <pkgdown>/BS5/assets/katex-auto.js to katex-auto.js
Copying <pkgdown>/BS5/assets/lightswitch.js to lightswitch.js
Copying <pkgdown>/BS5/assets/link.svg to link.svg
Copying <pkgdown>/BS5/assets/pkgdown.js to pkgdown.js
Copying pkgdown/extra.css to extra.css
Copying pkgdown/extra.js to extra.js
Copying pkgdown/assets/validation_suite.html to validation_suite.html
Updating deps/bootstrap-5.3.8/bootstrap.bundle.min.js
Updating deps/bootstrap-5.3.8/bootstrap.bundle.min.js.map
Updating deps/bootstrap-5.3.8/bootstrap.min.css
Updating deps/bootstrap-5.3.8/font.css
Updating deps/bootstrap-5.3.8/fonts/07d40e985ad7c747025dabb9f22142c4.woff2
Updating deps/bootstrap-5.3.8/fonts/1Ptug8zYS_SKggPNyC0ITw.woff2
Updating deps/bootstrap-5.3.8/fonts/1Ptug8zYS_SKggPNyCAIT5lu.woff2
Updating deps/bootstrap-5.3.8/fonts/1Ptug8zYS_SKggPNyCIIT5lu.woff2
Updating deps/bootstrap-5.3.8/fonts/1Ptug8zYS_SKggPNyCMIT5lu.woff2
Updating deps/bootstrap-5.3.8/fonts/1Ptug8zYS_SKggPNyCkIT5lu.woff2
Updating deps/bootstrap-5.3.8/fonts/1f5e011d6aae0d98fc0518e1a303e99a.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCs6KVjbNBYlgoKcQ72j00.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCs6KVjbNBYlgoKcg72j00.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCs6KVjbNBYlgoKcw72j00.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCs6KVjbNBYlgoKew72j00.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCs6KVjbNBYlgoKfA72j00.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCs6KVjbNBYlgoKfw72.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCv6KVjbNBYlgoCxCvjs2yNL4U.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCv6KVjbNBYlgoCxCvjsGyN.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCv6KVjbNBYlgoCxCvjtGyNL4U.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCv6KVjbNBYlgoCxCvjvGyNL4U.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCv6KVjbNBYlgoCxCvjvWyNL4U.woff2
Updating deps/bootstrap-5.3.8/fonts/4iCv6KVjbNBYlgoCxCvjvmyNL4U.woff2
Updating deps/bootstrap-5.3.8/fonts/626330658504e338ee86aec8e957426b.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7jsDJT9g.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7ksDJT9g.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7nsDI.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7osDJT9g.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7psDJT9g.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7qsDJT9g.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xK1dSBYKcSV-LCoeQqfX1RYOo3qPZ7rsDJT9g.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qN67lqDY.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qNK7lqDY.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qNa7lqDY.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qNq7lqDY.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qO67lqDY.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qOK7l.woff2
Updating deps/bootstrap-5.3.8/fonts/6xK3dSBYKcSV-LCoeQqfX1RYOo3qPK7lqDY.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwkxduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwlBduz8A.woff2
Updating deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwlxdu.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwmBduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwmRduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwmhduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3i54rwmxduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwkxduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwlBduz8A.woff2
Updating deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwlxdu.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwmBduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwmRduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwmhduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ig4vwmxduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwkxduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwlBduz8A.woff2
Updating deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwlxdu.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwmBduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwmRduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwmhduz8A.woff2
Updating
deps/bootstrap-5.3.8/fonts/6xKydSBYKcSV-LCoeQqfX1RYOo3ik4zwmxduz8A.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNBeudwk.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNReuQ.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNdeudwk.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNheudwk.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNleudwk.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNpeudwk.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR54z1Qlv-GDxkbKVQ_dFsvWNteudwk.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fO0KTet_.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fO4KTet_.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fO8KTet_.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fOAKTQ.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fOMKTet_.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fOQKTet_.woff2
Updating deps/bootstrap-5.3.8/fonts/CSR64z1Qlv-GDxkbKVQ_fOwKTet_.woff2
Updating
deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvQlMIXxw.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvUlMI.woff2
Updating
deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvXlMIXxw.woff2
Updating
deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvYlMIXxw.woff2
Updating
deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvZlMIXxw.woff2
Updating
deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvalMIXxw.woff2
Updating
deps/bootstrap-5.3.8/fonts/HI_QiYsKILxRpg3hIP6sJ7fM7PqlONvblMIXxw.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlM-vWjMY.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlMOvWjMY.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlMevWjMY.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlMuvWjMY.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlOevWjMY.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlPevW.woff2
Updating deps/bootstrap-5.3.8/fonts/HI_SiYsKILxRpg3hIP6sJ7fM7PqlPuvWjMY.woff2
Updating deps/bootstrap-5.3.8/fonts/JTUSjIg1_i6t8kCHKm459W1hyzbi.woff2
Updating deps/bootstrap-5.3.8/fonts/JTUSjIg1_i6t8kCHKm459WRhyzbi.woff2
Updating deps/bootstrap-5.3.8/fonts/JTUSjIg1_i6t8kCHKm459WZhyzbi.woff2
Updating deps/bootstrap-5.3.8/fonts/JTUSjIg1_i6t8kCHKm459Wdhyzbi.woff2
Updating deps/bootstrap-5.3.8/fonts/JTUSjIg1_i6t8kCHKm459Wlhyw.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3-UBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3CUBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3GUBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3KUBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3OUBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3iUBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMa3yUBA.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMawCUBGEe.woff2
Updating
deps/bootstrap-5.3.8/fonts/KFO7CnqEu92Fr1ME7kSn66aGLdTylUAMaxKUBGEe.woff2
Updating deps/bootstrap-5.3.8/fonts/QGYpz_kZZAGCONcK2A4bGOj8mNhN.woff2
Updating deps/bootstrap-5.3.8/fonts/S6u8w4BMUTPHjxsAUi-qJCY.woff2
Updating deps/bootstrap-5.3.8/fonts/S6u8w4BMUTPHjxsAXC-q.woff2
Updating deps/bootstrap-5.3.8/fonts/S6u9w4BMUTPHh6UVSwaPGR_p.woff2
Updating deps/bootstrap-5.3.8/fonts/S6u9w4BMUTPHh6UVSwiPGQ.woff2
Updating deps/bootstrap-5.3.8/fonts/S6u9w4BMUTPHh7USSwaPGR_p.woff2
Updating deps/bootstrap-5.3.8/fonts/S6u9w4BMUTPHh7USSwiPGQ.woff2
Updating deps/bootstrap-5.3.8/fonts/S6uyw4BMUTPHjx4wXg.woff2
Updating deps/bootstrap-5.3.8/fonts/S6uyw4BMUTPHjxAwXjeu.woff2
Updating
deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa0ZL7SUc.woff2
Updating deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa1ZL7.woff2
Updating
deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa1pL7SUc.woff2
Updating
deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa25L7SUc.woff2
Updating
deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa2JL7SUc.woff2
Updating
deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa2ZL7SUc.woff2
Updating
deps/bootstrap-5.3.8/fonts/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa2pL7SUc.woff2
Updating deps/bootstrap-5.3.8/fonts/XRXV3I6Li01BKofIMeaBXso.woff2
Updating deps/bootstrap-5.3.8/fonts/XRXV3I6Li01BKofINeaB.woff2
Updating deps/bootstrap-5.3.8/fonts/XRXV3I6Li01BKofIO-aBXso.woff2
Updating deps/bootstrap-5.3.8/fonts/XRXV3I6Li01BKofIOOaBXso.woff2
Updating deps/bootstrap-5.3.8/fonts/XRXV3I6Li01BKofIOuaBXso.woff2
Updating deps/bootstrap-5.3.8/fonts/c2f002b3a87d3f9bfeebb23d32cfd9f8.woff2
Updating deps/bootstrap-5.3.8/fonts/ee91700cdbf7ce16c054c2bb8946c736.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqW106F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWt06F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWtE6F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWtU6F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWtk6F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWu06F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWuU6F.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWuk6F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWvU6F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memtYaGs126MiZpBA-UFUIcVXSCEkx2cmqvXlWqWxU6F15M.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTS-muw.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTS2mu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTSCmu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTSGmu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTSKmu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTSOmu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTSumu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTSymu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTUGmu1aB.woff2
Updating
deps/bootstrap-5.3.8/fonts/memvYaGs126MiZpBA-UvWbX2vVnXBbObj2OVTVOmu1aB.woff2
Updating deps/bootstrap-5.3.8/fonts/q5uGsou0JOdh94bfuQltOxU.woff2
Updating deps/bootstrap-5.3.8/fonts/q5uGsou0JOdh94bfvQlt.woff2
Updating deps/bootstrap-toc-1.0.1/bootstrap-toc.min.js
Updating deps/clipboard.js-2.0.11/clipboard.min.js
Updating deps/font-awesome-6.5.2/css/all.css
Updating deps/font-awesome-6.5.2/css/all.min.css
Updating deps/font-awesome-6.5.2/css/v4-shims.css
Updating deps/font-awesome-6.5.2/css/v4-shims.min.css
Updating deps/font-awesome-6.5.2/webfonts/fa-brands-400.ttf
Updating deps/font-awesome-6.5.2/webfonts/fa-brands-400.woff2
Updating deps/font-awesome-6.5.2/webfonts/fa-regular-400.ttf
Updating deps/font-awesome-6.5.2/webfonts/fa-regular-400.woff2
Updating deps/font-awesome-6.5.2/webfonts/fa-solid-900.ttf
Updating deps/font-awesome-6.5.2/webfonts/fa-solid-900.woff2
Updating deps/font-awesome-6.5.2/webfonts/fa-v4compatibility.ttf
Updating deps/font-awesome-6.5.2/webfonts/fa-v4compatibility.woff2
Updating deps/headroom-0.11.0/headroom.min.js
Updating deps/headroom-0.11.0/jQuery.headroom.min.js
Updating deps/jquery-3.6.0/jquery-3.6.0.js
Updating deps/jquery-3.6.0/jquery-3.6.0.min.js
Updating deps/jquery-3.6.0/jquery-3.6.0.min.map
Updating deps/search-1.0.0/autocomplete.jquery.min.js
Updating deps/search-1.0.0/fuse.min.js
Updating deps/search-1.0.0/mark.min.js
Copying pkgdown/favicon/apple-touch-icon.png to apple-touch-icon.png
Copying pkgdown/favicon/favicon-96x96.png to favicon-96x96.png
Copying pkgdown/favicon/favicon.ico to favicon.ico
Copying pkgdown/favicon/favicon.svg to favicon.svg
Copying pkgdown/favicon/site.webmanifest to site.webmanifest
Copying pkgdown/favicon/web-app-manifest-192x192.png to
web-app-manifest-192x192.png
Copying pkgdown/favicon/web-app-manifest-512x512.png to
web-app-manifest-512x512.png
Copying logo.png
── Building pkgdown site for package pigauto ───────────────────────────────────
Reading from: /private/tmp/pigauto-cran011-site-9a457eb.qiPr9J
Writing to: /private/tmp/pigauto-cran011-site-9a457eb.qiPr9J/_site
── Sitrep ──────────────────────────────────────────────────────────────────────
✔ URLs ok.
✔ Favicons ok.
✔ Open graph metadata ok.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
✔ Articles metadata ok.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
✔ Reference metadata ok.
── Initialising site ───────────────────────────────────────────────────────────
── Building home ───────────────────────────────────────────────────────────────
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `authors.html`
Reading AGENTS.md
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `AGENTS.html`
Reading CLAUDE.md
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `CLAUDE.html`
Reading VALIDATION_LEDGER.md
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `VALIDATION_LEDGER.html`
Reading goodagents.md
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `goodagents.html`
Reading .github/404.md
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `404.html`
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `LICENSE-text.html`
Reading README.md
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `index.html`
── Building function reference ─────────────────────────────────────────────────
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/index.html`
Copying man/figures/logo.png to reference/figures/logo.png
Reading man/avonet300.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/avonet300.html`
Reading man/avonet_full.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/avonet_full.html`
Reading man/build_phylo_graph.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/build_phylo_graph.html`
Reading man/calibration_df.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/calibration_df.html`
Reading man/check_pigauto.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/check_pigauto.html`
Reading man/compare_methods.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/compare_methods.html`
Reading man/completed_data.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/completed_data.html`
Reading man/confusion_matrix.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/confusion_matrix.html`
Reading man/cross_validate.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/cross_validate.html`
Reading man/ctmax_sim.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/ctmax_sim.html`
Reading man/evaluate.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/evaluate.html`
Reading man/evaluate_imputation.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/evaluate_imputation.html`
Reading man/fit_baseline.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/fit_baseline.html`
Reading man/fit_pigauto.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/fit_pigauto.html`
Reading man/impute.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/impute.html`
Reading man/load_pigauto.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/load_pigauto.html`
Reading man/make_missing_splits.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/make_missing_splits.html`
Reading man/mask_missing.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/mask_missing.html`
Reading man/multi_impute.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/multi_impute.html`
Reading man/multi_impute_analysis.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/multi_impute_analysis.html`
Reading man/multi_impute_trees.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/multi_impute_trees.html`
Reading man/pigauto_report.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/pigauto_report.html`
Reading man/plot.pigauto_benchmark.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot.pigauto_benchmark.html`
Reading man/plot.pigauto_fit.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot.pigauto_fit.html`
Reading man/plot.pigauto_pred.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot.pigauto_pred.html`
Reading man/plot.pigauto_result.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot.pigauto_result.html`
Reading man/plot_comparison.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot_comparison.html`
Reading man/plot_history_gg.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot_history_gg.html`
Reading man/plot_uncertainty.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/plot_uncertainty.html`
Reading man/pool_mi.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/pool_mi.html`
Reading man/predict.pigauto_fit.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/predict.pigauto_fit.html`
Reading man/preprocess_traits.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/preprocess_traits.html`
Reading man/print.pigauto_check.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/print.pigauto_check.html`
Reading man/print.summary_pigauto_result.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/print.summary_pigauto_result.html`
Reading man/pull_gbif_centroids.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/pull_gbif_centroids.html`
Reading man/pull_worldclim_per_species.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/pull_worldclim_per_species.html`
Reading man/read_traits.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/read_traits.html`
Reading man/read_tree.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/read_tree.html`
Reading man/save_pigauto.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/save_pigauto.html`
Reading man/simulate_benchmark.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/simulate_benchmark.html`
Reading man/simulate_non_bm.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/simulate_non_bm.html`
Reading man/suggest_next_observation.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/suggest_next_observation.html`
Reading man/summary.pigauto_fit.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/summary.pigauto_fit.html`
Reading man/summary.pigauto_result.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/summary.pigauto_result.html`
Reading man/tree300.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/tree300.html`
Reading man/tree_full.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/tree_full.html`
Reading man/trees300.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/trees300.html`
Reading man/with_imputations.Rd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `reference/with_imputations.html`
── Building articles ───────────────────────────────────────────────────────────
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/index.html`
Reading vignettes/common-pitfalls.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/common-pitfalls.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
Reading vignettes/getting-started.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/getting-started.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
Reading vignettes/gnn-architecture.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/gnn-architecture.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
Reading vignettes/mixed-types.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/mixed-types.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
Reading vignettes/multiple-imputation.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/multiple-imputation.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
Reading vignettes/articles/simulation-study.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/simulation-study.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
Reading vignettes/tree-uncertainty.Rmd
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `articles/tree-uncertainty.html`
[WARNING] Deprecated: --mathjax. Use --math-method=mathjax[:URL] instead.
── Building news ───────────────────────────────────────────────────────────────
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
[WARNING] Deprecated: --mathml. Use --math-method=mathml instead.
Writing `news/index.html`
── Building sitemap ────────────────────────────────────────────────────────────
Writing `sitemap.xml`
── Building docs for llms ──────────────────────────────────────────────────────
── Building redirects ──────────────────────────────────────────────────────────
Adding redirect from articles/articles/simulation-study.html to
articles/simulation-study.html.
── Building search index ───────────────────────────────────────────────────────
── Checking for problems ───────────────────────────────────────────────────────
── Finished building pkgdown site for package pigauto ──────────────────────────
Warning messages:
1: @examplesIf condition "interactive()" is "FALSE" 
2: @examplesIf condition "interactive()" is "FALSE"
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/site-crawl-9a457eb.log`

~~~~text
{
  "html_pages": 62,
  "local_references_checked": 3588,
  "retired_pages": 34,
  "errors": []
}
SITE_CRAWL_OK
~~~~

## `docs/dev-log/cran-0.11-audit/provenance/windows-rhub-c2e768ad-2026-10-17.md`

~~~~text
# R-hub Windows checks for frozen artifact c2e768ad

The uploaded package object was read back from the R-hub package repository. Its SHA-256 is `c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2`, matching the frozen local archive at 5,136,713 bytes. R-hub's update job extracted that object and committed the package source as `9386c07717bc7a2dc63f072fcfbea2345b6e1e15`. The Windows jobs checked out that extracted source and ran a fresh `R CMD build` followed by `R CMD check --as-cran`; they did not execute the original gzip archive in place.

Run: [R-hub checks](https://github.com/r-hub2/comfortable-iguana-pigauto/actions/runs/39051372400). The upload object ID was `a5ff87d6cc9f5e0157d01dc97456aff5`.

| Job | Platform | Result | Tests |
|---|---|---|---|
| R-release | Windows Server 2025, R 4.6.1 | `Status: OK` | 2,442 passed, 0 failed, 28 warnings, 178 skips |
| R-devel | Windows Server 2025, R-devel r90686 | `Status: OK` | 2,442 passed, 0 failed, 28 warnings, 178 skips |

Both logs show the corrected `simulate_benchmark()` example passing under `--run-donttest`; portable-file-name checks pass. Rose independently checked both job logs, the upload-object identity chain, and the absence of the previous empty-results example failure. This is Windows evidence for source extracted from the exact uploaded archive, with R-hub rebuilding from that extracted source before checking.
~~~~

## `docs/dev-log/cran-0.11-audit/release-ledger.json`

~~~~text
{
  "package": "pigauto",
  "version": "0.11.0",
  "release_type": "update",
  "profile": {
    "compiled_code": false,
    "reverse_dependencies": false,
    "system_or_external_services": true,
    "large_data_or_vignettes": true
  },
  "artifact": {
    "state": "frozen exact candidate from clean source 9a457eb6; macOS exact-archive R CMD check Status OK; R-hub Windows release and devel check source extracted from byte-identical uploaded object and rebuilt; exact-hash installed MI matrix PASS; NOT_READY because candidate site has not been deployed and final release-validator/maintainer gates remain open",
    "filename": "pigauto_0.11.0.tar.gz",
    "sha256": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
    "size_bytes": 5136713,
    "archive_entries": 250,
    "inventory_sha256": "d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d",
    "source_commit": "9a457eb6b3363d413dbc502de88ff5bc602ed98a",
    "frozen_path": "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz",
    "path": "/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz",
    "clean_worktree": true,
    "clean_worktree_evidence": {
      "kind": "file",
      "value": "docs/dev-log/cran-0.11-audit/provenance/clean-source-checkout-780c084.txt"
    },
    "inventory_evidence": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-source-binding.md"
    },
    "forbidden_entries": []
  },
  "evidence": {
    "product_contract": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/GATES.md"
    },
    "rights_and_consent": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md"
    },
    "current_cran_policy": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md"
    },
    "rendered_site": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/site-9a457eb-review.md"
    },
    "local_as_cran": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-as-cran-check.log"
    },
    "incoming_feasibility": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-as-cran-check.log"
    },
    "installed_package": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-testthat.Rout"
    },
    "timing": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-c2e768ad-as-cran-check.log"
    },
    "local_check_status": "exact_frozen_tarball_as_cran_Status_OK_incoming_OK; full donttest examples; 3195 pass 0 fail 161 warn 86 skip",
    "local_check_log": "provenance/exact-artifact-c2e768ad-as-cran-check.log",
    "local_check_log_sha256": "93180d66bf6ee3ddeed26806ef1a998849d9759bc6eea5e2abda6292e4d48478",
    "testthat_output": "provenance/exact-artifact-c2e768ad-as-cran-testthat.Rout",
    "windows_results_bound_to_hash": true,
    "birdtree_redistribution_rights": "maintainer_directly_confirmed_use_with_citation_and_credit;no_separate_published_terms_identified;not_an_independent_legal_opinion",
    "birdtree_rights_basis": {
      "maintainer_warranty": "Shinichi directly confirmed in this task that the trees may be used with citation and credit; recorded as the maintainer authorization and warranty for this package release.",
      "published_redistribution_terms": "not_found_on_reviewed_birdtree_pages",
      "notice_distinction": "inst/NOTICE distinguishes BirdTree data from the megatrees software licence",
      "exact_tarball_inspection": "Recomputed directly from frozen archive SHA-256 fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a on 2026-10-17: inst/NOTICE and all six shipped .rda object streams are individually hashed in provenance/exact-artifact-fe1600d7-rights-data-hashes.txt; this confirms the current artifact contents, not third-party redistribution terms."
    },
    "exact_hash_independent_panel": "Grace and Rose independently verified exact upload object identity and both Windows logs; local site reviewer verified the 9a457eb build, rendered example, retirements, sitemap and search; site review is local only; overall release validator and deployed-site gate remain open",
    "local_check_command": "env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --run-donttest --timings --output=/private/tmp/pigauto-ascran-c2e768ad-network /private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz",
    "local_check_platform": "macOS arm64, R 4.6.0",
    "testthat_output_sha256": "54590c0e9219591cf4c7e19947d660e599cdf575c0669889623731a76f6621fb",
    "testthat_summary": {
      "pass": 3195,
      "fail": 0,
      "warn": 161,
      "skip": 86
    },
    "win_builder_upload_state": "R-hub object a5ff87d6cc9f5e0157d01dc97456aff5 read back byte-identical to frozen archive; updater extracted object to package repo commit 9386c07717bc7a2dc63f072fcfbea2345b6e1e15; both Windows jobs checked rebuilt source from that commit; hash-bound source chain, not execution of original gzip archive in place",
    "site_deployment_commit": "780c08470b263cc26b4139c711a5c144fd29389c",
    "site_deployment_run": "38970604986",
    "site_retirement_check": "fresh current deployment: 44 retired 404, 62 retained 200, sitemap 62, search 613, no retired targets",
    "site_g7_status": "OPEN: source 9a457eb is not merged/deployed; prior deployed site 780c084 does not prove current candidate routes",
    "site_verifier_source_sha256": "fbba18e851bad34fd6d22fcf803bdd04fd1730632a84243ada98aab2d77e7a15",
    "site_retired_routes_sha256": "5414f6cb69e9d140c38d0a6308f8e0ddbe8952834eea01b6b861a78eecb807e1",
    "site_retained_routes_sha256": "7288c1c2984c760f11a672f406311016ccefd193bc449402a4bbd645038f8ff1",
    "site_discovery_sha256": "2db884b26493dc295846bf51f9a9f5c8c29b63860edacafc55fb609285d46262",
    "site_verification_log_sha256": "7f8eee37739f522088bdd11b3975dcc9c1d204a3a433c1af7995e61bca04f944",
    "site_verifier_tests_sha256": "33060fd73862c5e1af97335504d5b45211fce6d229657f556d8deceeadc52ce4",
    "win_builder_logs": {
      "r_release": {
        "url": "https://win-builder.r-project.org/ys47lC1Ld715/00check.log",
        "result_url": "https://win-builder.r-project.org/ys47lC1Ld715/",
        "r_version": "R 4.6.1 (2026-06-24 ucrt)",
        "platform": "Windows Server 2022 x64",
        "started_utc": "2026-10-17T00:55:24Z",
        "elapsed_seconds": 271,
        "status": "OK",
        "association": "operationally_attributed_to_fe1600d7_by_upload_filename_size_and_timing;server_digest_unavailable"
      },
      "r_devel": {
        "url": "https://win-builder.r-project.org/m53luBDZqnM7/00check.log",
        "result_url": "https://win-builder.r-project.org/m53luBDZqnM7/",
        "r_version": "R-devel r90686 (2026-10-16 ucrt)",
        "platform": "Windows Server 2022 x64",
        "started_utc": "2026-10-17T01:04:23Z",
        "elapsed_seconds": 271,
        "status": "OK",
        "association": "operationally_attributed_to_fe1600d7_by_upload_filename_size_and_timing;server_digest_unavailable"
      }
    },
    "site_browser_review": {
      "date": "2026-10-16",
      "reviewer": "Pat",
      "source_commit": "780c08470b263cc26b4139c711a5c144fd29389c",
      "pages": [
        "NEWS index historical-text exception"
      ],
      "result": "44 retired routes 404; 62 retained routes 200; sitemap/search contain no retired targets; historical NEWS text explicitly labelled and says no release is implied",
      "scope_limit": "did not replay all 106 routes or inspect all retained page contents"
    },
    "exact_artifact_optional_backend_matrix": {
      "status": "PASS_bounded_installed_integration_four_cells_on_replacement_exact_tarball",
      "receipt": "after-task-2026-10-17-artifact-replacement.md",
      "manifest": "provenance/mi-matrix-fe1600d7-provenance-complete/manifest.json",
      "manifest_sha256": "e81c14f45d3d3837a25999f21d4f83f60276a918f532adfd5c9f2ddba0bb1991",
      "artifact_sha256": "fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a",
      "cells": [
        "neither",
        "drm-only",
        "gllvm-only",
        "both"
      ],
      "steps": 12,
      "review": "PASS; independent reviewer verified archive, runner, harness, five executed sources, manifest and all 12 logs; all four cells and 12 steps pass; scope/limits recorded in provenance/independent-review-fe1600d7-2026-10-17.md",
      "does_not_close": [
        "clean-room dependency installation",
        "broad statistical recovery or interval coverage",
        "joint covariance, loading, or variance-component pooling",
        "checksum-bound Windows results",
        "overall release readiness"
      ]
    },
    "as_cran_block_log": "provenance/exact-artifact-fe1600d7-as-cran-network-block.log",
    "as_cran_block_log_sha256": "ba0a3373b7b8dc9f5a757b80a0a7840790839f4975730c5e666a2a73e12eb09c",
    "rights_exact_artifact": {
      "NOTICE_sha256": "9cc7fa419cacea99dbb0b822c59d6fffad514ee072e119a1ce1c7f365883ba02",
      "tree300_sha256": "6050855d2ae7de42fa0de239134aad8e8975a224610b8dffb48f3d7412bb383f",
      "trees300_sha256": "8293e06503e89e79d9233abd1cbae98fc8d0c18e05ef5a41cca99a071fee910a",
      "tree_full_sha256": "5f7b36a4b9c02fd53c98d12fe4871e0f34098e3aa7a8f34f2db38a21a8dc09f5",
      "basis": "maintainer warranty with citation/credit; no independent legal opinion"
    },
    "site_local_refresh": {
      "source_commit": "9a457eb6b3363d413dbc502de88ff5bc602ed98a",
      "build_log": "provenance/site-build-9a457eb.log",
      "build_log_sha256": "97469ccfb06af9074366ef61b939be1c55fe1057d4a1168ba18fff774a5d7364",
      "crawl_log": "provenance/site-crawl-9a457eb.log",
      "crawl_log_sha256": "1bf804fb25425e4bbe6138f82646bb8c323a79094ae352ea399f41790e765128",
      "review_receipt": "provenance/site-9a457eb-review.md",
      "review_receipt_sha256": "e60aa37a8bf265d1f75ab7bdc14b8286c7dce8b72d4e976e1df18c7502d6fe00",
      "build": "PASS_62_html_614_search_pkgdown_clean",
      "crawl": "PASS_3588_refs_34_retired_zero_errors",
      "review": "local rendered content and accessibility tree checked; pixel appearance not assessed"
    },
    "main_source_ci": {
      "run": "38970604972",
      "macos_r_release": "success",
      "ubuntu_r_devel": "success",
      "ubuntu_r_release": "success; source CI only",
      "overall": "success"
    },
    "next_unproven_rung": "deployed-site-verified-and-release-validator-closure",
    "audit_verdict": "NOT_READY",
    "new_evidence_pr": "draft PR #236 open, unmerged; updated artifact and review receipts; site deployment gate open",
    "optional_backend_mi_exact_artifact": {
      "status": "PASS_bounded_installed_integration_on_current_frozen_hash",
      "archive_sha256": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
      "manifest": "provenance/mi-matrix-c2e768ad/manifest.json",
      "manifest_sha256": "bc29ddc1ba91dd7619ac44e2553c19588bd3cacceed1477fce68ffb84c94c525",
      "cells": [
        "neither",
        "drm-only",
        "gllvm-only",
        "both"
      ],
      "steps": 12,
      "independent_review": "PASS; independent reviewer verified manifest archive hash, source commit, 12 successful install/fixture/real steps and unchanged archive hash",
      "scope": "bounded analysis-aware MI; three real fits per available model cell, independent Gaussian and Rubin arithmetic, serialization/provenance",
      "limits": [
        "dependency libraries use shared symlinks, not clean-room installation",
        "not broad statistical recovery or interval-coverage evidence",
        "does not pool joint covariance, latent loadings or variance components"
      ]
    },
    "exact_tarball_inspection": "For new hash fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a, direct archive streaming confirms the prior SHA-256 values for inst/NOTICE and all six .rda objects; this is byte identity for those paths only, not third-party redistribution terms.",
    "current_shipped_data_hashes": {
      "archive_sha256": "fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a",
      "method": "SHA-256 of NOTICE and each R data object streamed from the frozen tarball",
      "notice_sha256": "9cc7fa419cacea99dbb0b822c59d6fffad514ee072e119a1ce1c7f365883ba02",
      "shipped_rda_sha256": {
        "pigauto/data/avonet300.rda": "86af337b49cbe2a81ff7acb7245b75c9a00cd48bb9131d7f4eb5d2676af131be",
        "pigauto/data/avonet_full.rda": "29b47bccdf17233328321378b77e99e8d4aa515d9df4451e347c709b97568e81",
        "pigauto/data/ctmax_sim.rda": "987a05458b7f3e548f22c52b75d8b7463e29a281a48895303f564effb0b9cecc",
        "pigauto/data/tree300.rda": "6050855d2ae7de42fa0de239134aad8e8975a224610b8dffb48f3d7412bb383f",
        "pigauto/data/tree_full.rda": "5f7b36a4b9c02fd53c98d12fe4871e0f34098e3aa7a8f34f2db38a21a8dc09f5",
        "pigauto/data/trees300.rda": "8293e06503e89e79d9233abd1cbae98fc8d0c18e05ef5a41cca99a071fee9108"
      }
    },
    "as_cran_block_log_status": "historical prior attempt on the replacement hash stopped before package checks during DNS failure; superseded by exact-hash as-cran completion recorded at local_check_log",
    "current_exact_candidate_source_binding": "provenance/exact-artifact-c2e768ad-source-binding.md",
    "current_exact_candidate_source_binding_sha256": "9dfe9c5d43feed2d7bf156f7edb3136c25d6d72e15201093b67797cdac474f6c",
    "current_exact_candidate_check": {
      "artifact_sha256": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
      "check_log": "provenance/exact-artifact-c2e768ad-as-cran-check.log",
      "check_log_sha256": "93180d66bf6ee3ddeed26806ef1a998849d9759bc6eea5e2abda6292e4d48478",
      "status": "Status_OK",
      "incoming_feasibility": "OK",
      "testthat": {
        "pass": 3195,
        "fail": 0,
        "warn": 161,
        "skip": 86
      },
      "testthat_output": "provenance/exact-artifact-c2e768ad-as-cran-testthat.Rout",
      "testthat_output_sha256": "54590c0e9219591cf4c7e19947d660e599cdf575c0669889623731a76f6621fb"
    },
    "independent_review_replacement_hash": "provenance/independent-review-fe1600d7-2026-10-17.md",
    "independent_review_replacement_hash_sha256": "8ed8966697aba0aba9c78750c8dece35565932327f4e6d0c500b72c76a073c5e",
    "manual_pdf_evidence": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-pigauto-manual.pdf"
    },
    "manual_pdf_sha256": "7f35acb0001c21f22ea05c78d69f2b301b2715152a74cdcc70a464960f6f50e9",
    "winbuilder_upload_acknowledgements": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/winbuilder-upload-fe1600d7-2026-10-17.md",
      "r_release": "HTTP_200_exact_filename_and_5136083_byte_acknowledgement; subsequent official log Status_OK; operational attribution, no server digest",
      "r_devel": "HTTP_200_exact_filename_and_5136083_byte_acknowledgement; subsequent official log Status_OK despite earlier form unavailable notice; operational attribution, no server digest",
      "server_sha256_attestation": false
    },
    "rights_data_hashes": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-fe1600d7-rights-data-hashes.txt"
    },
    "win_builder_previous_recheck_2026_10_16": {
      "checked_utc": "2026-10-16T23:49:04Z",
      "result_email_found": false,
      "upload_acknowledgement": "filename_and_size_only",
      "server_sha256_attestation": false,
      "effect": "Current candidate fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a has no checksum-bound Windows result; prior upload was for superseded hash 0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d",
      "note": "Historical mailbox/form recheck before replacement-hash submissions; it does not describe a recheck after the 2026-10-17 submissions."
    },
    "win_builder_submission_2026_10_17": {
      "checked_utc": "2026-10-17T00:40Z (approximate; upload completed before 00:53 UTC)",
      "receipt": "provenance/winbuilder-upload-fe1600d7-2026-10-17.md",
      "r_release": "HTTP_200_filename_and_size_acknowledged; subsequent R-release log Status_OK at 2026-10-17T00:55:24Z",
      "r_devel": "HTTP_200_filename_and_size_acknowledged; subsequent R-devel log Status_OK at 2026-10-17T01:04:23Z despite form unavailable notice",
      "mailbox_rechecked_after_upload": true,
      "server_sha256_attestation": false
    },
    "win_builder_results_receipt": "provenance/windows-rhub-c2e768ad-2026-10-17.md",
    "site_deployment_current": false,
    "windows_results_receipt": "provenance/windows-rhub-c2e768ad-2026-10-17.md",
    "windows_object_sha256": "c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2",
    "windows_object_bytes": 5136713,
    "windows_source_commit": "9386c07717bc7a2dc63f072fcfbea2345b6e1e15",
    "windows_run_id": 39051372400,
    "windows_r_release": {
      "job_id": 117450549024,
      "status": "Status: OK",
      "tests": {
        "pass": 2442,
        "fail": 0,
        "warn": 28,
        "skip": 178
      }
    },
    "windows_r_devel": {
      "job_id": 117450549044,
      "status": "Status: OK",
      "tests": {
        "pass": 2442,
        "fail": 0,
        "warn": 28,
        "skip": 178
      }
    },
    "external_logs": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/windows-rhub-c2e768ad-2026-10-17.md"
    },
    "platform_matrix": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/windows-rhub-c2e768ad-2026-10-17.md"
    }
  },
  "panel": {
    "grace": "PASS bounded exact artifact/source binding and local check; confirms Windows checks are downstream of exact uploaded tarball extraction",
    "rose": "PASS Windows R-release and R-devel checks for source extracted from exact uploaded tarball; fresh build before checks",
    "pat": "PASS bounded local candidate site content/routes; does not establish deployed-site parity"
  },
  "reader_surface_correction": {
    "pull_request": "https://github.com/itchyshin/pigauto/pull/234",
    "branch": "codex/cran-011-getting-started-correction",
    "commit": "6ebc67636f31fa8ff324ddd4a327d5cd808fba99",
    "merge_commit": "8082ef8273177ff664d49ac4bd643f7def671f20",
    "source_status": "merged_and_ancestor_of_frozen_source_commit_780c084",
    "local_vignette_knit": "passed_55_chunks_before_merge",
    "focused_article_build": "passed_before_merge",
    "independent_review": "passed_source_and_render_behavior_check_before_merge",
    "deployment": "included_in_site_deployment_commit_780c084_run_38970604986",
    "current_site_review": "Pat_passed_bounded_current_site_receipt; all-retained-page-content review not claimed",
    "historical_build_limitation": "the pre-merge full-site attempt was blocked by metadata/cache access; a later merged-source Pages deployment and fresh site receipts are recorded for 780c084",
    "after_task": "after-task-2026-10-15-getting-started-correction.md"
  },
  "status_claim": "platform-clean",
  "next_unproven_rung": "deployed-site-verified-and-release-validator-closure",
  "submission": {},
  "audit_verdict": "NOT_READY",
  "panel_evidence": [
    "provenance/exact-artifact-c2e768ad-source-binding.md",
    "provenance/exact-artifact-c2e768ad-as-cran-check.log",
    "provenance/exact-artifact-c2e768ad-as-cran-testthat.Rout",
    "provenance/mi-matrix-c2e768ad/manifest.json",
    "provenance/windows-rhub-c2e768ad-2026-10-17.md",
    "provenance/site-9a457eb-review.md"
  ],
  "evidence_pr": {
    "number": 236,
    "url": "https://github.com/itchyshin/pigauto/pull/236",
    "branch": "release/cran-0.11-gate-780c084",
    "state": "open_draft_NOT_READY",
    "merged": false,
    "last_rechecked_utc": "2026-10-17T17:03:00Z",
    "draft": true,
    "mergeable": true,
    "last_observed_head_sha": "9a457eb6b3363d413dbc502de88ff5bc602ed98a",
    "last_observed_commit_count": 8,
    "observation_note": "Current artifact c2e768ad is frozen from clean source commit 9a457eb. Exact local macOS check and MI matrix pass. R-hub object readback matches hash; both Windows jobs pass after rebuilding source extracted from that object. Fresh local site review passes. Draft remains unmerged; current deployment and release-validator closure remain open."
  },
  "superseded_previous_artifact": {
    "state": "superseded_archive_overwritten_by_rebuild; no surviving copy found",
    "filename": "pigauto_0.11.0.tar.gz",
    "sha256": "0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d",
    "size_bytes": 5136073,
    "archive_entries": 250,
    "inventory_sha256": "d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d",
    "source_commit": "780c08470b263cc26b4139c711a5c144fd29389c",
    "frozen_path": "/private/tmp/pigauto_0.11.0.tar.gz",
    "path": "/private/tmp/pigauto_0.11.0.tar.gz",
    "clean_worktree": true,
    "clean_worktree_evidence": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/clean-source-checkout-780c084.txt"
    },
    "inventory_evidence": {
      "kind": "file",
      "value": "/private/tmp/pigauto-release-evidence-780c084/docs/dev-log/cran-0.11-audit/provenance/exact-artifact-780c084.inventory.txt"
    },
    "forbidden_entries": [],
    "replacement_sha256": "fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a",
    "superseded_reason": "A reproducibility build wrote to the conventional path and replaced this file; searches found no surviving copy. Its bound checks and reviews are historical only."
  },
  "historical_superseded_artifact_before_0fccf4bb": {
    "state": "previous_frozen_candidate_stale_after_source_merge_rebuild_required",
    "filename": "pigauto_0.11.0.tar.gz",
    "sha256": "f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4",
    "size_bytes": 5135324,
    "archive_entries": 250,
    "inventory_sha256": "10a31a50b2504f42514b0ec8b5368861543cf1f902d9ae075ebb5279cd5a777e",
    "source_commit": "d76804e768bf77f430f43dcf591633f9cd900dab",
    "frozen_path": "/private/tmp/pigauto-release-candidate/sha256-f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4/pigauto_0.11.0.tar.gz",
    "path": "/private/tmp/pigauto-release-candidate/sha256-f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4/pigauto_0.11.0.tar.gz",
    "clean_worktree": true,
    "clean_worktree_evidence": {
      "kind": "file",
      "value": "/Users/z3437171/.codex/worktrees/cran-011-gate-reconcile/pigauto/docs/dev-log/cran-0.11-audit/provenance/post-merge-source-checkout-status-2026-10-10.txt"
    },
    "inventory_evidence": {
      "kind": "file",
      "value": "/Users/z3437171/.codex/worktrees/cran-011-gate-reconcile/pigauto/docs/dev-log/cran-0.11-audit/provenance/post-merge-tarball-inventory.txt"
    },
    "forbidden_entries": []
  }
}
~~~~

