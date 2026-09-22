# Gates: campaign (S6a core on Totoro, S6b factorial on nibi, S6c AVONET300, S6d covariate sensitivity)

OWNS: results directories on Totoro and nibi /project; script/campaign_sim_results/** after rsync

Scope: every design cell has its rds for every replicate; failures are recorded, never dropped.

- [ ] G10: core slice complete: 18 cells x 200 replicates present (BACE arm present on replicates 1..100), failure rate column filled, 0 missing rds
  CHECK: Rscript script/campaign_sim_checks.R --gate G10 --dir script/campaign_sim_results/raw_index.csv
  EXPECT: G10 PASS
  EVIDENCE: pending

- [ ] G11: trimmed factorial complete: 40 cells (BM/OU x lambda {0.3, 1.0} x rho {0, 0.5} x n {100, 1000} x mechanisms {MCAR 0.1, MAR 0.3, clade 0.3} minus the core), same replicate conditions
  CHECK: Rscript script/campaign_sim_checks.R --gate G11 --dir script/campaign_sim_results/raw_index.csv
  EXPECT: G11 PASS
  EVIDENCE: pending

- [ ] G14: cross-host reproducibility: 5 cells run on both Totoro and nibi under RNGkind L'Ecuyer-CMRG agree on truth, mask and the frequentist arm's predictions to 1e-10
  CHECK: Rscript script/campaign_sim_checks.R --gate G14
  EXPECT: G14 PASS
  EVIDENCE: pending

- [ ] G6c: AVONET300 case study present for every arm, 20 seeds
  EVIDENCE: pending

- [ ] G6d: covariate sensitivity present for the 18 core cells (arm 1 phylolm/phyloglm variant; arm 3 recorded as covariate-free)
  EVIDENCE: pending
