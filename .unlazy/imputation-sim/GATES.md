# Gates: imputation-sim (phylogenetic imputation simulation study, arc/imputation-sim)

OWNS: script/campaign_sim*, script/campaign_gnn_off_lib.R, script/campaign_gnn_off_aggregate.R, script/campaign_gnn_off_tables.R, script/campaign_gnn_off_figures.R, docs/dev-log/arc/**, vignettes/articles/**, .Rbuildignore, _pkgdown.yml

Scope: run the four-arm imputation simulation (frequentist, BACE, pigauto GNN off 3a/3b, GNN on) to completion on the corrected design, with pre-run before launch, paired MCSE on every number, and three deliverables (Artifact, BACE methods note, pkgdown article).

Plan: /Users/z3437171/.claude/plans/distributed-snacking-turtle.md

Leaves:
- gates/leaf-env.md      G1 G7 G8            (S2a, S2b, S2c)
- gates/leaf-runner.md   G2 G3 G4 G5 G6 Gold (S3, S3-verify)
- gates/leaf-prerun.md   G9                  (S4)  -> G0 Shinichi
- gates/leaf-campaign.md G10 G11 G14         (S6, MECHANICAL-VERIFY)
- gates/leaf-results.md  G12 G13             (S7, S7a-c)
