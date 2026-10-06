# Checkpoint, 2026-10-06

PLATFORM: Codex | ON BRANCH: `release/cran-0.11-audit` | LANE: pre-CRAN 0.11 audit. The direct-to-main lane uses another checkout. Work here is isolated, with no edits to that lane.

S0 to S2 are complete. Version 0.11.0, publication provenance, all 33 exports, and 126 bounded default-path assertions are recorded. Four isolated installed-library adapter cells passed. The final-adapter drmTMB and gllvmTMB recovery campaigns each completed ten registered seeds and passed all prespecified fixed-effect bias margins and per-seed mechanism checks. The seven-seed gllvmTMB continuation finished in 5,236 seconds under two one-core Totoro workers. Coverage counts remain descriptive at ten seeds. Totoro retains all raw RDS files; all 21 copied gllvmTMB CSV/log files matched remote SHA-256 digests.

S3 to S5 are complete in this worktree. Active reader surfaces and generated help agree with current defaults; 39 historical site assets were archived byte-for-byte. The fresh local site passed 65 pages and 3,621 local references, with 34 retired direct pages absent; desktop and mobile views were inspected. The public deployment remains unchanged until merge.

A follow-up Symbolizer comparison narrowed the README and getting-started downstream-inference wording. That rebuild and refreshed search index passed 65 pages and 3,622 local references. The rendered prose was checked locally.

The multiple-imputation article also now distinguishes fit/extraction examples from validated downstream inference. Its subsequent full local rebuild passed 65 pages, 3,621 local references, and 34 retired pages absent. The rendered prose was checked; the earlier viewport inspection was not repeated.

A pre-PR tarball built from the corrected worktree has SHA-256 `a9c7ce7b65aec751edd0c8246697ecb51eda9ffe19c874566ad246b4b67c4852`. `R CMD check --as-cran --no-manual` on macOS/R 4.6.0 finished `Status: OK`; the check log is in `docs/dev-log/cran-0.11-audit/provenance/pre-pr-check.log`. The exact release candidate must still be frozen from the merged source. The release ledger remains `NOT_READY`.

Draft audit PR [#226](https://github.com/itchyshin/pigauto/pull/226) is open against `main`. The last pushed head, `bd9a9f4`, passed all three R CMD check jobs (Ubuntu release/devel and macOS release). The completed gllvmTMB results and corrected BirdTree provenance comment are being added to the audit branch; CI will need to rerun after push. It remains a draft; human merge is a separate gate. The other active lane is the direct-to-main checkout, which this worktree has not touched.

Next: commit the final recovery receipts, push the audit PR, and verify its checks. After human merge, verify the deployed site, freeze the exact tarball, run platform checks, resolve the open BirdTree redistribution-rights point, and prepare the separate evidence PR. No CRAN submission is authorized.
