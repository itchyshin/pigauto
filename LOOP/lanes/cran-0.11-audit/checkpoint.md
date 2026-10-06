# Checkpoint, 2026-10-06

## Current post-merge status

PR #226 merged as `bb5835d1b214d783da7b0b414df99aa6ba926bc7`. Pages run #570
completed successfully and published that commit. The live homepage, getting-started
article, multiple-imputation article, and `multi_impute()` reference are current.
All 34 retired benchmark URLs and four retired walkthrough URLs return 404; the
64-location sitemap and 574-path search index contain no retired routes. The live
homepage has no per-trait benchmark links. Details are in
`docs/dev-log/cran-0.11-audit/post-merge-verification-2026-10-06.md`.

The exact merged-source tarball is `pigauto_0.11.0.tar.gz`, SHA-256
`d34d469981386546b4ad03b9c7aab81674556708eeaf524c1c7365b6c68dec26`
(5,122,013 bytes; 247 members; no forbidden development paths). Its full
macOS/R 4.6.0 `R CMD check --as-cran --no-manual` completed `Status: OK`.
The merged commit also passed Ubuntu R release, Ubuntu R-devel, and macOS R
release checks in Actions run #37526171411. The exact tarball has been submitted
to win-builder's R-release and R-devel forms; those results are pending. G7
remains partial until those exact-artifact platform results are recorded.
Retained logs and artifact metadata are under
`docs/dev-log/cran-0.11-audit/provenance/`.

The bounded G8 audit criterion is complete. Fresh Grace, Rose and Pat reviews
returned READY; 13 planted validator controls failed closed; cache-busted direct
browser checks confirmed the four current reader surfaces; artifact and check
hashes were re-verified. Receipts are in
`docs/dev-log/cran-0.11-audit/provenance/review-panel-2026-10-06.md` and adjacent
files.

Next unmet gates: the exact tarball's Windows R-release and R-devel results, then
BirdTree underlying-data redistribution rights. The release ledger and CRAN
authorization remain `NOT_READY` until platform evidence is complete and that
right is resolved. No upstream contact or CRAN submission has occurred.

## Historical pre-merge checkpoint (retained for sequence)


PLATFORM: Codex | ON BRANCH: `release/cran-0.11-audit` | LANE: pre-CRAN 0.11 audit. The direct-to-main lane uses another checkout. Work here is isolated, with no edits to that lane.

S0 to S2 are complete. Version 0.11.0, publication provenance, all 33 exports, and 126 bounded default-path assertions are recorded. Four isolated installed-library adapter cells passed. The final-adapter drmTMB and gllvmTMB recovery campaigns each completed ten registered seeds and passed all prespecified fixed-effect bias margins and per-seed mechanism checks. The seven-seed gllvmTMB continuation finished in 5,236 seconds under two one-core Totoro workers. Coverage counts remain descriptive at ten seeds. Totoro retains all raw RDS files; all 21 copied gllvmTMB CSV/log files matched remote SHA-256 digests.

S3 to S5 are complete in this worktree. Active reader surfaces and generated help agree with current defaults; 39 historical site assets were archived byte-for-byte. The fresh local site passed 65 pages and 3,621 local references, with 34 retired direct pages absent; desktop and mobile views were inspected. The public deployment remains unchanged until merge.

A follow-up Symbolizer comparison narrowed the README and getting-started downstream-inference wording. That rebuild and refreshed search index passed 65 pages and 3,622 local references. The rendered prose was checked locally.

The multiple-imputation article also now distinguishes fit/extraction examples from validated downstream inference. Its subsequent full local rebuild passed 65 pages, 3,621 local references, and 34 retired pages absent. The rendered prose was checked; the earlier viewport inspection was not repeated.

A pre-PR tarball built from the corrected worktree has SHA-256 `a9c7ce7b65aec751edd0c8246697ecb51eda9ffe19c874566ad246b4b67c4852`. `R CMD check --as-cran --no-manual` on macOS/R 4.6.0 finished `Status: OK`; the check log is in `docs/dev-log/cran-0.11-audit/provenance/pre-pr-check.log`. The exact release candidate must still be frozen from the merged source. The release ledger remains `NOT_READY`.

Draft audit PR [#226](https://github.com/itchyshin/pigauto/pull/226) is open against `main`. Commit `70089b9` passed all three R CMD check jobs (Ubuntu release/devel and macOS release). A subsequent reader check corrected the shipped getting-started R companion and one NEWS sentence. The rebuilt local site passed 65 pages and 3,622 local references. CI must rerun after those two corrections are pushed. The PR remains a draft; human merge is a separate gate. The other active lane is the direct-to-main checkout, which this worktree has not touched.

Next: commit and push the final reader corrections, then verify PR checks. After human merge, verify the deployed site, freeze the exact tarball, run platform checks, resolve the open BirdTree redistribution-rights point, and prepare the separate evidence PR. No CRAN submission is authorized.
