# GNN incremental-value campaign: plan versus actual

## Planned execution

The approved plan required a frozen ADEMP design and acceptance ledger, deterministic checks and same-fit ablation tests, a bounded pilot, a Totoro capacity check, and approval before any production estimate above three hours. Production was allocated at 500 datasets in each of five scenarios, with both gate settings. The approved production ceiling was 630 minutes, using no more than 48 single-thread workers and staying below the 150-core ownership cap. The article edit waited for the shared validation lane to release `vignettes/validation-results.Rmd`. Deployment was deferred at the initial closeout. The user then explicitly authorized pushing and deploying; the closeout record was updated after live verification. No API or method change was authorized.

The plan anticipated a pilot of up to 50 datasets and about one hour. The pilot established that production would exceed three hours, so the user explicitly approved Totoro production before it started. Initial throughput projections during production varied from about 11 to 12.5 hours and were treated as provisional. The plan did not specify a reliable aggregate wall-clock estimate for all implementation and closeout work; no single total estimate is claimed here.

## Actual execution

The implementation and its review finished in the isolated `codex/gnn-incremental-value` worktree. A source fix keeps internal validation and test outcomes out of the optional phylogenetic-signal estimate. Negative controls exposed issues in the comparator, audit and small-pilot classification; each issue was corrected and tested before production.

The 50-dataset pilot covered all five scenarios. A 10-dataset-per-scenario checkpoint was retained after correcting the pilot's confirmatory classification. Invalid earlier outputs remain preserved and labelled.

Totoro ran production with 48 single-thread workers from 2026-10-19 20:40 MDT to 2026-10-20 05:20 MDT. Elapsed time was about 520 minutes (8 h 40 min), within the approved 630-minute cap. The run completed 5,000 pigauto fits and 30,000 method attempts. Capacity snapshots showed substantial memory headroom and total owned CPU use below 150. DRAC and Kohaku were not needed for this CPU workload.

The strict audit passed on Totoro. We transferred the full archive and verified its SHA-256, then regenerated summaries locally. The local strict audit rejected the different runtime identity; the Totoro receipt remains the authoritative audit.

We waited for PR #239 to merge and release the article path. Then we copied its committed article into this worktree and edited the isolated copy. The older untracked article in the other worktree remains untouched. After PR #239 merged, we rebased the 30 local commits onto `b8153d8`; the only conflict was the article add/add, resolved by retaining the released article plus this GNN section. The branch now has no diff deletions against current main.

The full package tests passed on the rebased tree (3,655 passed, 0 failed, 8 skipped, 184 warnings), and pkgdown build/check completed locally. `devtools::check(error_on = "warning")` passed with zero errors, warnings or notes in 9m 21s. Unlazy re-verification and local commits completed the initial local closeout. The user subsequently authorized pushing and deployment. PR #240 merged as `78b1306d18d2a27b6f1caef574c0782c1cc3b4dc`; its required checks settled green. The `pkgdown.yaml` push run #758 completed successfully on `main`, including the Pages deploy step. The live validation article was checked and contains the measured table and scenario limits. See `DEPLOYMENT-20261020.md`.

## Deviations and unresolved items

The exact production source bundle had no Git metadata, so `git_head` is `NA`. Source identity is instead bound by the bundle SHA-256 and 66 per-file hashes. This is documented in the scientific report and run configuration. The total time for all planning, lane waiting, implementation, review, and closeout was not separately instrumented; only production elapsed time is measured. Optional-data tests were skipped because their external caches were unavailable, and these skips are not treated as passes.
