# GNN incremental predictive-value study: acceptance gates

Scope: determine whether pigauto's gated neural correction adds held-out observation-level predictive value in the frozen repeated-observation scenarios. This ledger governs only the runner, report, and local evidence paths. The validation-results pkgdown article is reserved until the active validation lane releases it.

The user approved testing both `phylo_signal_gate = FALSE` (the shipped default and confirmatory policy) and `TRUE` (policy sensitivity). They are separate paired fits on identical generated data, outer masks, and fit seed. The gate is evaluated before training; after calibration it can force the neural weight to zero, so a final-weight toggle is not an equivalent policy comparison.

OWNS: `script/gnn_incremental/**`, `docs/dev-log/gnn-incremental/**`, `.unlazy/gnn-incremental/**`

- [x] G1: The generator, estimand, and phylogenetic GLS predictions match the frozen ADEMP contract and an independently calculated small covariance oracle; the deliberately broken control is detected.
  CHECK: `Rscript --vanilla script/gnn_incremental/test_harness.R`
  EXPECT: `HARNESS_VALIDATION_PASS`
  EVIDENCE: `Rscript --vanilla script/gnn_incremental/test_harness.R` returned `HARNESS_VALIDATION_PASS` on 2026-10-19; local worktree source is bound by the final bundle SHA recorded in the pilot receipt.

- [x] G2: The same-fit ablation changes only the gated neural contribution, preserves the fitted linear covariate term and mean-baseline weight, and keeps predictions row-aligned.
  CHECK: `Rscript --vanilla script/gnn_incremental/test_ablation.R`
  EXPECT: `ABLATION_VALIDATION_PASS`
  EVIDENCE: `Rscript --vanilla script/gnn_incremental/test_ablation.R` returned `ABLATION_VALIDATION_PASS` on 2026-10-19.

- [x] G3: An expanded 50-dataset pilot from the corrected source bundle has valid outputs, exact source/environment identity, and reconciled requested/completed/failed/nonconverged attempts across all five scenarios. Earlier invalid pilots remain retained: the quadratic-basis error archive and the comparator-row-selector archive at `script/gnn_incremental/results/archive/invalid-comparator-pilot-20261019.tar.gz` (SHA-256 `761c7fb9c80812f288590172a830b9ae462ac7fe74f230c86bfc36eac5e4db7c`). The selector fix has a regression fixture. Pilot v3 used source bundle SHA-256 `1d65ad064fb01bd6c297d0bf6cfb81a65d7fc78a848607f30a6c59ebb5590217`; all 50 dataset records, 100 pigauto gate fits, and 600 method attempts reconciled, with zero failures, nonconvergence, or skips. Totoro audit returned `RUN_ACCOUNTING_PASS`; local re-summarization from the verified retained archive returned `SUMMARY_REGENERATION_PASS`. An added negative-first regression check caught and fixed pilot-size confirmatory classification; with 10 replicates, confirmatory outcomes now remain indeterminate below the planned 500. Receipt: `docs/dev-log/gnn-incremental/PILOT-20261019.md`.
  CHECK: `Rscript --vanilla script/gnn_incremental/audit_run.R --run-dir results/pilot-20261019-v3-n10`
  EXPECT: `RUN_ACCOUNTING_PASS`
  EVIDENCE: Totoro audit retained in `script/gnn_incremental/results/pilot-20261019-v3-n10.tar.gz`; pilot receipt `docs/dev-log/gnn-incremental/PILOT-20261019.md`.

- [ ] G4: The guarded Totoro production launcher refuses to run without explicit approval and a fresh capacity-check marker, accepts only a unique persistent run directory and the exact source-bundle hash, limits concurrency to at most 48 single-thread workers, and initializes exactly 500 datasets per frozen scenario. Production then contains all attempts with no duplicates, silent skips, or stale reuse, and accounts for comparator-specific failures.
  CHECK: `Rscript --vanilla script/gnn_incremental/audit_run.R --run-dir results/production`
  EXPECT: `RUN_ACCOUNTING_PASS`
  EVIDENCE: The exact extracted compact production bundle from HEAD `a9540510a642a21584cbd1f798bbd718d5c62baa`, SHA-256 `0b50049ce8a1f6f1aaa81b9bb8a3005dee9a5620a9bf7dd0d68f0a660316fffb`, passes package installation, harness, ablation, accounting, and launcher guard checks both locally and on Totoro in an isolated library. Production remains pending user approval. Expanded-pilot estimate: 8 h 39 min mean raw compute, 9 h 43 min conservative raw-compute budget, and 10–10.5 h including run overhead.

- [ ] G5: Production scientific summaries regenerate from retained replicate-level results and report paired uncertainty across datasets, both missingness strata, gate behavior, failures, and convergence. The 50-dataset pilot summary reproduces locally from its retained archive and remains explicitly indeterminate because it is below the planned 500 replicates per scenario.
  CHECK: `Rscript --vanilla script/gnn_incremental/summarise.R --verify --run-dir results/production`
  EXPECT: `SUMMARY_REGENERATION_PASS`
  EVIDENCE: Pilot `SUMMARY_REGENERATION_PASS` on local retained archive; production remains pending approval and run.

- [ ] G6: The pkgdown results article renders and its package-site structure passes the package's pkgdown check after the active validation lane releases the shared article path.
  CHECK: `Rscript --vanilla -e 'pkgdown::build_site(); pkgdown::check_pkgdown()'`
  EXPECT: `No problems detected`
  EVIDENCE: pending

- [ ] G7: Independent statistical review confirms each claim matches its measured regime and that limitations, deviations, and unsupported outcomes remain visible.
  REVIEWER: Fresh statistical-review lens; review the frozen manifest, retained result tables, and final report.
  EVIDENCE: pending

- [ ] G8: The after-task report and plan-versus-actual reconciliation are complete; all gates are reverified and local changes are committed.
  CHECK: `Rscript --vanilla script/gnn_incremental/check_closeout.R`
  EXPECT: `CLOSEOUT_VALIDATION_PASS`
  EVIDENCE: pending

Independent review update (2026-10-19): the reviewer found the original audit trusted the stored configuration fingerprint and accepted any unique task indices. Added verification of the fingerprint against the serialized settings and exact task sequence checks in both the runner and audit, with negative controls in `test_accounting.R`. The positive and tampered-run checks pass, and follow-up statistical review passed. An independent direct `phytools::phylosig` oracle now checks split masking. The first pilot also detected the missing `train_idx` assumption in the comparator harness; an observed-value plus split-mask selector and regression fixture now pass, and one short route check confirms OLS/GLS complete. The fresh five-scenario pilot remains pending.
