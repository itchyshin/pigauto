# Comprehensive scientific validation

This campaign tests existing pigauto behaviour. It does not change methods or defaults.

## Start

Use the isolated worktree, read docs/dev-log/validation-20261018/PLAN.md and inspect .unlazy/validation-20261018/GATES.md. Run the harness, runner and contract checks first.

`Rscript script/validation_20261018/run.R --mode manifest` freezes the scenario table and expanded jobs. `bash script/validation_20261018/pilot.sh` runs a serial pilot, limited to30minutes overall and10minutes per job. It selects the first job per adapter plus large-n prediction sentinels and uses production settings. It is a feasibility/timing sample, not calibration evidence. Missing attempts remain missing.

`Rscript script/validation_20261018/summarise.R --pilot` verifies result accounting. A scientific failure does not make evidence disappear. Production requires an explicit RDS approval receipt with approved=TRUE, current source_identity, runtime_identity and jobs_hash bound to the frozen full manifest, plus exact workers (integer), budget_sec, job_timeout_sec and execution_context (totoro or drac). No receipt is created automatically.

## Interpretation

Supported inference requires independent truth, suitable imputation and analysis models, adequate convergence and reported Monte Carlo uncertainty. The scientific criteria are coverage0.925-0.975, absolute meanSE/empiricalSD0.9-1.1 and bias/empiricalSD<=0.1. Insufficient precision is indeterminate. Stress regimes are descriptive. No package-wide conformal coverage guarantee is asserted.

Source identity covers package R files, DESCRIPTION/NAMESPACE and this campaign. Results carry host, R/session packages, controls, warnings, diagnostics and runtime. Resume rejects stale identity or seeds. Archives/receipts and raw RDS outputs are retained separately. Production jobs must be immutable; any code fix creates a new identity and affected evidence must be rerun.

## Compute

Recheck Totoro CPU, memory and existing user jobs before launch. This pilot uses one worker and single-thread libraries. Production starts at no more than48workers, memory-limited, with total user usage below150cores. Use existing DRAC connections if headroom is inadequate, with explicit account/time/memory and persistent project storage. Never use GitHub Actions for campaigns.

## Review and delivery

Keep executable run-state ledgers untracked; force-add the durable plan/report under ignored docs/. Record every failed, unavailable, timed-out and missing attempt. Independent reviewers assess the harness, code and final claims. Full-production scientific results remain pending until approved and executed.

References: Morris, White and Crowther (2019), https://doi.org/10.1002/sim.8086; Williams et al.(2024), https://doi.org/10.1111/2041-210X.14415.

## Production dispatch

The production launcher partitions the exact approved job manifest across at most48single-thread workers, each with separate output and logs. Its approval preflight runs before workers start. Totoro context is allowed only on the Totoro host; DRAC context requires an active scheduler allocation with at least the declared worker CPUs. The launcher binds the exact total budget and per-job timeout to the receipt. Global summary reconciles worker outputs against the full manifest, not just files that happened to finish.

Prediction uses current-default GNN-off and mean/mode controls. The older broad comparator roster is not rerun. GNN comparisons are explicit in the covariate and nonlinear workflow slice. No unavailable dataset is silently replaced by another one.

## Follow-up pilot and reporting

`Rscript script/validation_20261018/run.R --mode pilot --supplement --out <directory>` selects seventeen full-setting cases: all twelve analysis size/mechanism cells, AVONET and configured BIEN, tree and posterior n=300, and a MAR prediction sentinel. Use the same owned GNU timeout wrapper and thread caps; the retained follow-up uses a fifteen-minute global cap and five-minute job caps.

Run `summarise.R --pilot --out <directory>` separately for each frozen identity. Raw identities are never rewritten. The evaluator records its own file hashes when reanalysing older pilot records. It exports attempts, complete result rows, scalar summaries, paired method differences, reliability bins and posterior interval tables. Legacy intervals without expected-cell counts remain unavailable.

`estimate.R <primary-pilot-directory> <other-pilot-directory> ...` combines timing observations only, retaining each raw identity. Successful matched-size and mechanism routes replace failed startup timings, while separate failure-history flags remain visible. It never pools scientific results across source identities.

The runtime fingerprint includes RNG kind. BIEN requires both explicit input-path environment variables, whose bytes are hash-bound. The current calibration metric is fixed-bin confidence ECE, with full masked-cell accounting; it is descriptive and does not assert a nominal-coverage guarantee.
