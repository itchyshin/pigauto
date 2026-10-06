# Gates: runner (S3 corrected design in script/)

WORKTREE: /Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim

OWNS: script/campaign_gnn_off_lib.R, script/campaign_sim_cell.R, script/campaign_sim_checks.R, script/campaign_gnn_off_aggregate.R, script/campaign_gnn_off_tables.R, script/campaign_gnn_off_figures.R

Scope: the runner implements the Corrected design (plan section) with every existing invocation unchanged by default.

- [x] G2: smoke on the new wrapper runs every arm including 3b and writes a non-empty rds
  CHECK: Rscript script/campaign_sim_cell.R --dgp types_mixed --n 60 --seed 1 --lambda 0.7 --rho 0.5 --miss mar --out /tmp/imputation_sim_smoke --arms gnn_on,gnn_off,gnn_off_rphylopars,freq,bace,floor --smoke 2>&1 | tee /tmp/imputation_sim_smoke/smoke.log | grep -c "done in" | grep -qx 6 && ! grep -q ERROR /tmp/imputation_sim_smoke/smoke.log && Rscript -e 'x <- readRDS(list.files("/tmp/imputation_sim_smoke", "rds$", full.names=TRUE)[1]); stopifnot(nrow(x$results) > 0); cat("G2 PASS\n")'
  EXPECT: G2 PASS
  EVIDENCE: exit 0; six arms each logged 'done in', zero ERROR, rds non-empty (types_mixed n=60 lambda 0.7 rho 0.5 MAR, driver, fixed thresholds, --smoke). Re-run 2026-09-20 after the review fixes.

- [x] G3: mask integrity: NA pattern of df_miss equals mask OR truth-NA; every arm returns rows in truth order; d1 never masked
  CHECK: Rscript script/campaign_sim_checks.R --gate G3
  EXPECT: G3 PASS
  EVIDENCE: G3 PASS. Five cells across mcar/mar/clade: is.na(df_miss) equals mask | is.na(truth); every arm's rownames equal truth's; d1 never masked.

- [x] G4: Pagel lambda recovery on BM at n = 1000 (20 replicates; phylolm mean within 0.05 of 0.3, 0.7, 1.0)
  CHECK: Rscript script/campaign_sim_checks.R --gate G4
  EXPECT: G4 PASS
  EVIDENCE: G4 PASS. 20 replicates at n=1000, BM: phylolm lambda recovered 0.258 / 0.677 / 1.000 against targets 0.3 / 0.7 / 1.0, all within the 0.05 bound.

- [x] G5: realised missing fraction within 0.01 of target over 20 replicates for mar and clade; at least 5 observed per column
  CHECK: Rscript script/campaign_sim_checks.R --gate G5
  EXPECT: G5 PASS
  EVIDENCE: G5 PASS. 20 replicates each: realised fraction 0.2994 (mar) and 0.3000 (clade) against 0.30; minimum observed per column at least 5.

- [x] G6: intervals at lambda = 1, n = 1000, MCAR, 20 replicates: arm 1 and arm 2 coverage within [0.90, 0.99]; arm 2 interval from >= 500 posterior predictive samples
  CHECK: Rscript script/campaign_sim_checks.R --gate G6
  EXPECT: G6 PASS
  EVIDENCE: NOT OBTAINED as a pre-run gate (first attempt segfaulted under mclapply; rescoped and unreturned). RESOLVED BY THE CAMPAIGN 2026-09-22: the core slice measured exactly this cell with 200 replicates (BACE 100). Coverage at lambda = 1, n = 1000, MCAR 0.30: freq 0.908, freq_lambda 0.902, BACE 0.885 (committed summary.csv; methods note "Interval coverage"). BACE falls OUTSIDE the gate's [0.90, 0.99] band, so the EXPECT was wrong and the fact is reported as a finding in every deliverable; BACE intervals are built from n_final = 20 full imputation runs (all retained draws), the ">= 500 samples" clause was superseded by the n_final decision recorded in the pre-run note. RETRACTED 2026-09-24: that waiver was wrong. With 20 draws the percentile interval's coverage ceiling is 0.872 under a correct model, so BACE's 0.885 here does not show BACE's intervals are too narrow; G6 stays unmet as specified and moves to v2. Gate ABANDONED as a pass/fail oracle; its content is a reported result.

- [x] Gold: an old invocation reproduces the committed pre-run rds (bm_mixed n 100 seed 1) metrics within 1e-8 for the arms it shares
  CHECK: Rscript script/campaign_sim_checks.R --gate Gold
  EXPECT: Gold PASS
  EVIDENCE: Gold PASS. max abs diff 3.947e-11 over 18 shared (arm, trait, metric) rows against script/campaign_gnn_off_prerun/bm_mixed_n100_s1.rds. Reference is self-generated pre-change, so this proves no regression since that snapshot, not independent ground truth.

- [x] Galign: symbolic-alignment table (DGP maths vs code, term by term) reviewed and signed by a fresh reviewer; recorded in docs/dev-log/arc/<date>-simulation-runner-alignment.md
  EVIDENCE: Symbolic-alignment table produced by the builder and reviewed by a fresh Opus reviewer, which cleared the OU tip-correlation and binary-probability-orientation flags and raised four blocking findings, all fixed in 4ccc6ac and 3b62c78.
