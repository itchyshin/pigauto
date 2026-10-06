# Checkpoint, 2026-10-06

PLATFORM: Codex | ON BRANCH: `release/cran-0.11-audit` | LANE: pre-CRAN 0.11 audit. The direct-to-main lane uses another checkout. Work here is isolated, with no edits to that lane.

S0 and S1 are complete. Version 0.11.0, publication provenance, all 33 exports, and 126 bounded default-path assertions are recorded. S2 has four isolated installed-library adapter cells passing, a ten-seed final-adapter drmTMB replay meeting the preset bias margin, and three completed gllvmTMB mechanism seeds. Seeds 2026100602 and 2026100603 ran concurrently in 1,309 seconds under two one-core workers. Seven seeds remain. At the measured pace, four two-worker waves would take about 87 minutes; this replaces the earlier 191.3-minute serial projection for planning. No further gllvmTMB seed has been launched.

S3 to S5 are complete in this worktree. Active reader surfaces and generated help agree with current defaults; 39 historical site assets were archived byte-for-byte. The fresh local site passed 65 pages and 3,621 local references, with 34 retired direct pages absent; desktop and mobile views were inspected. The public deployment remains unchanged until merge.

A follow-up Symbolizer comparison narrowed the README and getting-started downstream-inference wording. The rebuilt site and refreshed search index passed 65 pages and 3,622 local references; the previous 3,621-reference crawl remains the pre-PR snapshot. The rendered prose was checked locally.

A pre-PR tarball built from the corrected worktree has SHA-256 `a9c7ce7b65aec751edd0c8246697ecb51eda9ffe19c874566ad246b4b67c4852`. `R CMD check --as-cran --no-manual` on macOS/R 4.6.0 finished `Status: OK`; the check log is in `docs/dev-log/cran-0.11-audit/provenance/pre-pr-check.log`. The exact release candidate must still be frozen from the merged source. The release ledger remains `NOT_READY`.

Draft audit PR [#226](https://github.com/itchyshin/pigauto/pull/226) is open against `main` from commit `60195260c9ef2e74a60c3f9c41700b7b5dbb27da`. It remains a draft; human merge is a separate gate. The other active lane is the direct-to-main checkout, which this worktree has not touched.

Next: review PR checks and coordinate the seven remaining gllvmTMB seeds under a two-worker, bounded plan. After human merge, verify the deployed site, freeze the exact tarball, run platform checks, resolve the open BirdTree redistribution-rights point, and prepare the separate evidence PR. No CRAN submission is authorized.
