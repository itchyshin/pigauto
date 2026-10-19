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

- [x] G3: A five-dataset pilot from the corrected source bundle has valid outputs, exact source/environment identity, and reconciled requested/completed/failed/nonconverged attempts across all five scenarios. The earlier 25-dataset attempt remains invalid because the quadratic-basis matrix had the wrong shape. The first corrected-source pilot is also invalid for comparator assessment because the runner referenced nonexistent `splits$train_idx`, leaving OLS/GLS with no training rows; this run is retained at `script/gnn_incremental/results/archive/invalid-comparator-pilot-20261019.tar.gz` (SHA-256 `761c7fb9c80812f288590172a830b9ae462ac7fe74f230c86bfc36eac5e4db7c`). The row selector fix is exercised by a regression fixture and the fresh route check. Corrected pilot v2 used source bundle SHA-256 `0619d51b6abe99a08a07cb5b4da65e4be2349940c5f7c1ef4d9ea02f83ca91bc`; 60/60 method attempts succeeded with zero failures or nonconvergence. Totoro audit returned `RUN_ACCOUNTING_PASS`; summary regeneration returned `SUMMARY_REGENERATION_PASS`. Local re-summarization of the retained archive also returned `SUMMARY_REGENERATION_PASS`. Receipt: `docs/dev-log/gnn-incremental/PILOT-20261019.md`.
  CHECK: `Rscript --vanilla script/gnn_incremental/audit_run.R --run-dir results/pilot-20261019-v2`
  EXPECT: `RUN_ACCOUNTING_PASS`
  EVIDENCE: Totoro `logs/audit.log`; pilot receipt `docs/dev-log/gnn-incremental/PILOT-20261019.md`.

- [ ] G4: Production contains exactly 500 requested datasets per frozen scenario only after the post-pilot runtime estimate checkpoint is explicitly approved; no duplicates, silent skips, or stale-result reuse are allowed, and all attempts are accounted for, including comparator-specific failures.
  CHECK: `Rscript --vanilla script/gnn_incremental/audit_run.R --run-dir results/production`
  EXPECT: `RUN_ACCOUNTING_PASS`
  EVIDENCE: Pending valid pilot and user approval if the measured production estimate exceeds three hours.

- [ ] G5: Production scientific summaries regenerate from retained replicate-level results and report paired uncertainty across datasets, both missingness strata, gate behavior, failures, and convergence. Pilot summary reproduction passed, but five datasets cannot support efficacy inference.
  CHECK: `Rscript --vanilla script/gnn_incremental/summarise.R --verify --run-dir results/production`
  EXPECT: `SUMMARY_REGENERATION_PASS`
  EVIDENCE: Pilot `SUMMARY_REGENERATION_PASS` on Totoro and local retained archive; production remains pending approval and run.

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
