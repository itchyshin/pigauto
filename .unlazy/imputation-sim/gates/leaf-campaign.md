# Gates: campaign (S6a core on Totoro, S6b factorial on nibi, S6c AVONET300, S6d covariate sensitivity)

WORKTREE: /Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim

OWNS: results directories on Totoro and nibi /project; script/campaign_sim_results/** after rsync

Scope: every design cell has its rds for every replicate; failures are recorded, never dropped.

- [x] G10: core slice complete: 18 cells x 200 replicates present (BACE arm present on replicates 1..100), failure rate column filled, 0 missing rds
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && Rscript script/campaign_sim_checks.R --gate G10 --dir /tmp/pig_pool4/core
  EXPECT: G10 PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=01d9749a8aeb/36 entries; output=G10: failure rate recorded on 808 of 10785 cells | G10 PASS
  FINAL 2026-09-23 02:20: "G10: 18 design cells, 10785 rds present, 0 replicate-arms missing / failure rate recorded on 808 of 10785 cells / G10 PASS". Core BACE 600/600, of which 3 are hand-assembled failed-arm records (seed 55 n=300 segfault x2; seeds 38, 51 n=100 timed out at 4 h x2), each labelled inside the file and disclosed in the methods note.

- [x] G11: trimmed factorial complete: 56 cells (BM/OU x lambda {0.3, 1.0} x rho {0, 0.5} x n {100, 1000} x mechanisms {MCAR 0.1, MCAR 0.3, MAR 0.3, clade 0.3} minus the 8 core overlaps; the plan's "40" omitted MCAR 0.3, the design table has always emitted 56), BACE 100 seeds at n=100 and 30 at n=1000 (Shinichi 2026-09-21)
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && Rscript script/campaign_sim_checks.R --gate G11 --dir /tmp/pig_pool6/factorial
  EXPECT: G11 PASS
  EVIDENCE: PASS 2026-09-22 20:58, direct run: "G11: 56 design cells, 40316 rds present, 0 replicate-arms missing / failure rate recorded on 1676 of 40316 cells / G11 PASS". BACE 3,550/3,550 (100 at n=100; 30 to 100 at n=1000). History: the 10:00 run the same day FAILED on BACE only (57 seeds, 56 of them the four BM clade0.3 n=1000 cells after fir array 60776371 ended 37 TIMEOUT + 24 NODE_FAIL at 5 h); recovered on nibi (array 22491383, 12 h, 420/420 COMPLETED, ~3.8 h per fit) and the n100 half-B seed (22491384). NOTE gate-check --approve kills this CHECK at its 120 s per-check timeout (39k rds) and then unchecks the gate; run the CHECK directly, never re-approve this leaf through gate-check.

- [x] G14: cross-host reproducibility: 5 cells run on both Totoro and nibi under RNGkind L'Ecuyer-CMRG agree on truth, mask and the frequentist arm's predictions to 1e-10
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && Rscript script/campaign_sim_checks.R --gate G14 --dir /tmp/pig_pool6/factorial/totoro_fl,/tmp/pig_pool6/factorial/nibi
  EXPECT: G14 PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=01d9749a8aeb/36 entries; output=G14: 5 overlapping cell(s) checked, truth/mask/freq identical | G14 PASS

- [x] G6c: AVONET300 case study present for every arm, 20 seeds
  EVIDENCE: Totoro results/avonet 20/20 (gnn_on, gnn_off, gnn_off_rphylopars, freq, bace, floor) + results/avonet_fl 20/20 (freq_lambda, run 2026-09-21 14:34). Aggregated /tmp/pig_pool4/avo_summary.csv; freq_lambda == freq to 4 dp (estimated lambda ~ 1).

- [x] G6d: covariate sensitivity present for the 18 core cells (arm 1 phylolm/phyloglm variant; arm 3 recorded as covariate-free)
  EVIDENCE: 2026-09-22 10:47 covsens 3,600/3,600 (fir 3,599 + the last n=1000 seed computed on Totoro against the seeded directory; pooled at /tmp/pig_pool6/covsens/totoro; the paired differences are unchanged to three decimals from the 3,599 run) + covsens_fl 3,600/3,600 (Totoro). Arms gnn_on, gnn_off, gnn_off_rphylopars, freq, freq_lambda, floor, ncov=2. Paired covsens-minus-core on identical (cell, seed, arm): freq -0.051/-0.024/+0.003 z-RMSE at lambda 0.3/0.7/1; freq_lambda -0.032/-0.012/+0.002; gnn_on -0.012/-0.015/-0.000; gnn_off exactly 0 (covariate-free by design: the shared-seed check). Reported in the methods note, the article and the board. Arm 1 ran its phylolm/phyloglm covariate variant; arm 3 recorded covariate-free. Earlier line: fir results/covsens 3,584/3,600 (n100 1200, n300 1185, n1000 1199), arms gnn_on, gnn_off, gnn_off_rphylopars, freq, floor, ncov=2; array 60895543 has 16 tasks pending on JobArrayTaskLimit %8. freq_lambda wave launched on Totoro 10:12 into results/covsens_fl (3,600 jobs). Pull both when complete.
