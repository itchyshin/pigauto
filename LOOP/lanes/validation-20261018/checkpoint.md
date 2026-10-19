# Checkpoint: production assessment and closeout

Branch `codex/comprehensive-validation-20261018`; base `b2e6054685e6bdedb073ba46e058af2c0a619a18`. Managed worktree `/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto`. Original dirty checkout remains untouched. No push, merge, release or brain writes.

Read `docs/dev-log/validation-20261018/scientific-report.md`, `plan-versus-actual.md`, `validation-ledger.md`, and `review-production.md`. The full 151-scenario, 32,240-job manifest ran on Totoro. Scientific assessment: 12 inferential rows failed, 32 are indeterminate and none is supported.

Prior evidence remains applicable within scope: base suite PASS 3554, FAIL 0, WARN 182, SKIP 9; contract fixtures 4 tests and 86 expectations, with no failures, errors or skips; available-dependency package check Status OK. The strict package check remains blocked by unavailable brms and rgbif. Pilot and follow-up results remain source-qualified historical stages.

Production approval is retained in `docs/dev-log/validation-20261018/receipts/production-approval.rds` and the resume receipt. Compressed result tables and execution log are in the same receipt folder. Totoro raw results remain at `/home/snakagaw/hsq_work/pigauto_validation_20261018/source-5e25a197/script/validation_20261018/results/production/`; checksums for the 32,384 worker RDS files are in `production-raw-rds.sha256`.

The expected-roster audit is complete: 5,697 scenario/method/estimand groups and 1,122,060 planned job-level keys are all present; zero absent, partial, unexpected or unrostered rows. Strict reliability-bin and posterior-interval payload checks found zero missing or failed groups. The five BIEN trait names and hashes were rederived read-only on Totoro and exactly compared with the local roster. Independent regeneration of 4,000 posterior DGP cells matched all 8,000 x/y interval denominators. Final independent review of the corrected status validation and receipt binding is pending. Preserve runtime inventory variation and scientific failures. Do not change package methods or defaults without a separate approved decision.
