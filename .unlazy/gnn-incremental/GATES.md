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

- [ ] G4: The guarded Totoro production launcher refuses to run without explicit approval and a fresh capacity-check marker, accepts only a unique persistent run directory and the exact source-bundle hash, limits concurrency to at most 48 single-thread workers, and initializes exactly 500 datasets per frozen scenario. Its audit also requires an explicit recognized terminal status for every method in every requested gate fit, including `not_attempted` with a visible count, and rejects missing or unknown statuses. Production then contains all attempts with no duplicates, silent skips, or stale reuse, and accounts for comparator-specific failures.
  CHECK: `Rscript --vanilla script/gnn_incremental/audit_run.R --run-dir results/production`
  EXPECT: `RUN_ACCOUNTING_PASS`
  EVIDENCE: Current compact source from HEAD `bbc7e8b1e0886f8236550b688e9a1211315f653c`, SHA-256 `76c594de2c063c559516652151ea63e26a13ed4d2097ced3e25dbaaacb1d4695`, passes package installation, harness, ablation, accounting, and launcher-guard checks locally and on Totoro in isolated libraries; Totoro verified the transferred checksum. User approved Totoro production on 2026-10-19. Fresh pre-launch capacity: load 63.23/63.38/64.26 on 384 cores, 414 GiB available, about 18 Shinichi-owned CPU equivalents before launch; 48 single-thread workers project about 66 total, under the 150-core cap. Production run `production-20261019-v1` was initialized on Totoro at 2026-10-19 20:40 MDT with 2,500 datasets / 5,000 gate-specific fits, using the exact verified bundle and isolated R library. Authorization and capacity snapshot are retained at `/home/snakagaw/hsq_work/pigauto-gnn-incremental-20261019-corrected/results/production-20261019-v1.approval.txt`; launcher PID 3351516. First startup check showed 48 active R workers, no task errors, and load 76.59; production outcomes remain pending. Expanded-pilot estimate: 8 h 39 min mean raw compute, 9 h 43 min conservative raw-compute budget, and 10–10.5 h including run overhead.

- [ ] G5: Production scientific summaries regenerate from retained replicate-level results and report paired uncertainty across datasets, both missingness strata, gate behavior, failures, and convergence. The 50-dataset pilot summary reproduces locally from its retained archive and remains explicitly indeterminate because it is below the planned 500 replicates per scenario.
  CHECK: `Rscript --vanilla script/gnn_incremental/summarise.R --verify --run-dir results/production`
  EXPECT: `SUMMARY_REGENERATION_PASS`
  EVIDENCE: Pilot `SUMMARY_REGENERATION_PASS` on local retained archive; production run `production-20261019-v1` is in progress on Totoro.

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

Independent review update (2026-10-19): an earlier review found the original audit trusted the stored configuration fingerprint and accepted any unique task indices. Added verification of the fingerprint against the serialized settings and exact task sequence checks in both the runner and audit, with negative controls in `test_accounting.R`. The positive and tampered-run checks pass. An independent direct `phytools::phylosig` oracle now checks split masking. The first pilot also detected the missing `train_idx` assumption in the comparator harness; an observed-value plus split-mask selector and regression fixture now pass, and one short route check confirms OLS/GLS complete. The expanded five-scenario pilot is retained below.

Independent statistical review update (2026-10-19): the reviewer confirmed the estimand, matched ablation, multiplicity plan, and pilot classifications. It found that the audit summarized method-level outcomes without rejecting a missing or unknown method status. Added missing/unknown-status negative controls and fail-closed terminal-status accounting. The new negative control failed against the old auditor and `test_accounting.R` now passes with the strengthened auditor. The retained pilot itself records all 600 method attempts successful; this audit gap did not alter the recorded pilot predictions, but the updated verifier still requires validation against the exact refreshed production bundle before G4 can pass. Public `impute()` defaults to `phylo_signal_gate = FALSE`; the campaign keeps `lambda_mode = "estimate"` in both gate arms.

Follow-up independent verification (2026-10-19, HEAD `edd9d6f62382c48af110be8bac08a5a42a4b9269`): PASS. The reviewer checked that the audit's terminal-status patterns cover all method statuses emitted by `campaign.R`, and ran `test_accounting.R` with positive, missing-status, and unknown-status fixtures. Output: `ACCOUNTING_AND_SUMMARY_VALIDATION_PASS`. This validates the local audit change; it does not substitute for testing the refreshed production bundle on Totoro.

Compact-bundle verification (2026-10-19): the runner-ready bundle omits prior results and unrelated documentation, while including package sources, tests, campaign scripts, manifest, and ADEMP specification. Its SHA-256 is `76c594de2c063c559516652151ea63e26a13ed4d2097ced3e25dbaaacb1d4695`. It installed as pigauto 0.11.0 and passed `test_harness.R`, `test_ablation.R`, `test_accounting.R`, and `test_production_launcher.sh` both locally and on Totoro. Outputs: `HARNESS_VALIDATION_PASS`, `ABLATION_VALIDATION_PASS`, `ACCOUNTING_AND_SUMMARY_VALIDATION_PASS`, and `PRODUCTION_LAUNCHER_GUARD_PASS`. These results verify the exact current source bundle and guarded launch path, not production results. G4 remains pending explicit approval and the 500-replicate run.
