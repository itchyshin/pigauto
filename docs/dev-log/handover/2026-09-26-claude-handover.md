# Session Handoff: continuous freq-vs-BACE campaign done; start the DISCRETE equivalent

Meta: 2026-09-26 · from Claude Code (Opus 5.5) to a fresh Claude Code session · lane `claude:pigauto:rubin`
(released) · worktree `/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace` · branch `arc/rubin-freq-bace`

You are Claude, picking up the pigauto Rubin lane. The continuous campaign is finished and published. Your job is
to **plan** the discrete-trait equivalent (binary, ordinal, categorical) and bring Shinichi a pre-run plan. Do not
launch a campaign without his approval: it will exceed the 3-hour line (D-287).

## Mission control

| Repo / worktree | Branch | State | What shipped | Next, by leverage |
|---|---|---|---|---|
| pigauto (worktree above) | `arc/rubin-freq-bace` @ `c05d880`, pushed, no PR | LANDED on the branch; not merged | Continuous campaign: 6,600 fits, report + summary + accuracy pages | 1. Discrete study design (this handover). 2. Shinichi decides whether the branch gets a PR. |

## Critical context

1. **The result files do NOT contain the imputed datasets.** Each `rds` keeps summary scores, the true data
   (`truth`, including `bin`, `ord`, `cat3`) and the `mask`, but not the 20 imputations. The discrete traits
   therefore **cannot be scored from the finished runs**; the discrete study needs new runs. (An earlier chat
   message said the opposite; it was wrong and was corrected.) New runs must save at least the imputed
   discrete columns.
2. **The continuous campaign was the frequentist method's home ground.** Only c1, c2 (and prp) were scored,
   under Brownian motion, where Rphylopars is the correct model. The frequentist arms impute only the
   continuous block; BACE imputes every trait. This campaign did not measure BACE's claimed advantage (one
   model for mixed trait types). Simulation v1 measured discrete per-value accuracy only, against castor (see
   Q1), without multiple imputation, Rubin pooling or a downstream estimand. The discrete study fills that gap.
3. **Hold Dan.** Shinichi's standing instruction: do not contact Dan Noble. The pages are private.

## What was accomplished (continuous campaign, 2026-09-24 to 26)

- 18 cells (n 100/300/1000 × λ 0.3/0.7/1 × ρ 0/0.5), M = 20, arms freq A (proper: parametric bootstrap then
  joint conditional draw), freq B (improper: one fit), BACE as shipped, chained BACE, BACE + residual (negative
  control), plus complete data. BACE 3,000 fits, frequentist 3,600, 0 missing.
- Verdict: BACE as shipped under-covers at λ 0.3/0.7 (about 0.91 per value; downstream correlation 0.90 at
  n = 300 and 1000) because its 20 final runs share one starting dataset. Chaining the final runs fixes it:
  chained BACE matches freq A on coverage and accuracy at n ≥ 300; at n = 100 it is less accurate and pulls
  the slope toward zero (about −0.02 against 0.001).
- Pages (private, Shinichi owns):
  - report https://claude.ai/artifact/PFkoRFtTjox4tndEbuBtPQ
  - summary https://claude.ai/artifact/9vNrmz2PUfjW9vZwuKmvFc
  - accuracy https://claude.ai/artifact/CGGmi4Km6unGPhxEiv94qk
  - pre-run https://claude.ai/artifact/WMHxQ8idafD5sLAVEJtg9Q
- Full detail, including every fix and residual: `docs/dev-log/after-task/2026-09-25-rubin-campaign.md`.

## Current working state

- Working: all runners, tests and build scripts; lane suite `testthat::test_dir("script/tests-rubin")` FAIL 0
  PASS 234; all 16 lane gates met (re-approved and re-run 2026-09-26).
- In progress: nothing. No jobs on any cluster.
- Not working / blocked: nothing in this lane.

## Key decisions and rationale (still binding)

- Never edit `R/` or `BACE/` in this lane. BACE is imported as installed.
- BACE settings for campaigns: runs 5, nitt 50,000, burnin 10,000, thin 25 (Shinichi; BACE's own convergence
  check does not measure mixing).
- Coverage target is the population value; correlation is GLS-whitened on Fisher z; scoring PGLS uses the fast
  eigenbasis REML estimator (`est_pgls_slope_fast`), which is more exact than `nlme::gls` at boundary λ.
- Compute: D-287 (estimate first; over 3 h needs a plan, a pre-run with results, and approval), D-143 (Totoro
  ≤150 cores for this user across lanes), D-64 (attach only via `~/.ssh/cm-*` sockets, never trigger Duo),
  never compute on a DRAC login node. Proven hosts: fir, rorqual, nibi, Totoro (narval's CPUs cannot run the
  fir-built library).
- `ScheduleWakeup` is forbidden in this project (see `CLAUDE.md`); use Monitor or background Bash.

## Landing state

`tools/handoff_gate.sh` output, annotated:

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| pigauto `arc/rubin-freq-bace` `c05d880` (+ this handover commit) | y | y | none | LANDED on branch |
| `.unlazy/rubin-freq-bace/gates/*` (16 gates) | untracked by design | n/a | n/a | ALL MET |
| `.unlazy/imputation-sim/gates/*` (5 ledgers UNMET) | n/a | n/a | n/a | PROTECTED: another lane's (`arc/imputation-sim`); not ours to close |
| "89 UNPUSHED on other branch(es)" | n/a | n/a | n/a | PROTECTED: other branches in the shared pigauto repo, not this lane |
| Result pool `~/pigauto_rubin_pool` (Mac, outside Dropbox; 361 MB) | n/a | n/a | n/a | local data of record; `agg/report.json` feeds all pages |

FINDING-OF-RECORD: BACE as shipped gives improper MI (shared final-run anchor); chaining fixes coverage;
chained BACE ties proper frequentist MI at n ≥ 300 on continuous traits. vault-note: NOT WRITTEN. Brain writes
need Shinichi's explicit approval; canonical source is this branch
(`docs/dev-log/after-task/2026-09-25-rubin-campaign.md`). Propose the note to him; do not write it unasked.

## Files created / modified (this lane, `git diff --stat 3f28e26..HEAD`)

- Runners: `script/rubin_lib.R`, `script/rubin_freq.R`, `script/rubin_bace.R`, `script/rubin_cell.R`,
  `script/rubin_checks.R`, `script/rubin_convcheck.R`, `script/rubin_prerun_summary.R`
- Cluster tooling: `script/rubin_campaign_nibi.sh`, `script/rubin_prerun_nibi.sh`, `script/rubin_prerun_totoro.sh`,
  `script/rubin_totoro.sh`, `script/rubin_status.sh`, `script/rubin_retry_prep.R`, `script/rubin_pool.sh`,
  `script/rubin_missing.R`
- Aggregation and pages: `script/rubin_campaign_aggregate.R`, `script/rubin_report_build.sh`,
  `script/rubin_report_template.html`, `script/rubin_summary_template.html`, `script/rubin_accuracy_template.html`
- Tests: `script/tests-rubin/test-*.R` (13 files)
- Docs (force-added; `docs/` is gitignored): `docs/dev-log/arc/2026-09-24-rubin-*.md`,
  `docs/dev-log/arc/2026-09-25-rubin-{prerun-results,campaign-report,summary,accuracy}.html`,
  `docs/dev-log/after-task/2026-09-2{4,5}-rubin-*.md`, `docs/dev-log/plan-actual/2026-09-24-rubin-steps1-2.md`,
  `.unlazy/rubin-freq-bace/{GOAL,arcs,checkpoint}.md`, and this file.

## Next immediate steps (OWED)

1. Run `~/shinichi-brain/tools/lane_preflight.sh "/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace"`
   and claim a lease for the discrete lane. Classify this handover's items OWED / DONE / RETRACTED / PROTECTED.
2. Open an ultra-plan for the discrete study (Shinichi asked for this lane to "start the discrete one"). Say
   back to him in one sentence what you will build before spending compute.
3. Settle the three design questions below with Shinichi, each with your recommendation drafted.
4. Feasibility smoke on the Mac (under 3 h): one cell, n = 100, the chosen frequentist discrete arm plus BACE,
   saving imputations, scoring the chosen estimands.
5. Pre-run plan with a measured time estimate, then STOP for approval (D-287).

## Open questions for Shinichi (draft a recommendation for each)

**Q1. What is the frequentist equivalent for discrete traits?** No off-the-shelf phylogenetic frequentist method
imputes ordinal or categorical traits jointly with continuous ones; that gap is itself part of the finding.
Candidates:
- (a) Liability-scale joint MVN: plug-in liabilities for the discrete cells, joint Rphylopars fit, draw
  liabilities, threshold back. pigauto already has this idea (`R/joint_threshold_baseline.R`, OVR for
  categorical), but it is an approximation, not a likelihood fit, and must be re-implemented in `script/`
  (no `R/` edits in this lane).
- (b) Frequentist chained equations with a phylogenetic GLMM per discrete trait (`glmmTMB` binomial/probit with
  a phylogenetic random effect; installed), parametric bootstrap for properness. This mirrors BACE's
  chained-equation structure, so it isolates "frequentist vs Bayesian". But `glmmTMB` has no ordinal family and
  no multinomial, and getting the random effect for species with a missing response needs a feasibility check.
- (c) **castor Mk hidden-state prediction, already used in simulation v1** (`run_freq()` in
  `script/campaign_gnn_off_lib.R`: `castor::hsp_mk_model`, equal rates for binary and categorical, stepwise
  SUEDE for ordinal, tip likelihoods kept as class probabilities). It covers all three discrete types, but it
  is single-trait: it cannot use c1, c2 or the driver. As used in v1 it is a plug-in (fixed fitted rates),
  so it is the analogue of freq B. A proper version (freq A analogue) needs a parametric bootstrap of the Mk
  rates (simulate states on the tree at the fitted rates, refit, recompute tip probabilities) before each
  draw. Draws from castor's per-tip marginal probabilities ignore dependence between tips, so a joint draw
  needs checking. v1 failure to reuse: castor errors on a monomorphic trait (common at λ = 1 with fixed
  thresholds), and v1 scored it at the floor.
- v1 already measured discrete ACCURACY (no MI, no Rubin, no downstream estimand), per
  `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md` and the board
  https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24 (v10). At λ = 0.3 BACE leads by 11 to 17 accuracy points
  (0.587 against 0.422 at n = 1000), and every other arm is at or below the mode floor. At λ = 1 castor and
  pigauto lead (0.96 to 0.97 against BACE's 0.85). A plausible reading, AGENT-INFERRED and untested: at low
  signal the tree carries little information, so BACE's use of the correlated continuous traits is what
  helps, and single-trait castor cannot do that. Note that v1's BACE ran before the chaining fix and the input
  cleaning.
- Suggested recommendation: use castor as the frequentist arm for all three discrete types (plug-in = freq
  B analogue; bootstrap = freq A analogue), plus (b) for binary only if the smoke shows the random effect for
  missing species can be recovered. Its lack of cross-trait information is the comparison the BACE paper
  needs.

**Q2. What do we score?** Per value: accuracy of the imputed class, Brier score and calibration of the pooled
class probabilities. There is no natural "coverage" for a class, so a coverage analogue would be the rate
at which the true class falls in the prediction set built from the 20 draws. Downstream: something with a
discrete trait in it. Options are the PGLS slope of c1 on `bin` (binary predictor) or a phylogenetic logistic
regression `bin ~ c1` (`phylolm::phyloglm`). Its true value is not simply ρ, so the truth must be computed
(analytically on the liability scale, or as a large Monte Carlo mean of the complete-data estimate).

**Q3. Scale.** The continuous design had 18 cells with 200 reps (100 for BACE at n = 1000). BACE cost is
dominated by n = 1000 (about 6 h per fit on fir, 19 GB). Suggested recommendation: n = 100 and 300 first
(λ 0.3/0.7/1, ρ 0/0.5, 200 reps), with n = 1000 decided after, as in the continuous campaign.

## Gotchas and failed approaches (do not repeat)

- `RNGkind("L'Ecuyer-CMRG")` leaks between test files unless restored; `make_cell()` datasets depend on it.
- BACE stops when a discrete level is unobserved (common at λ = 1): `fit_bace_mi()` drops unused levels. For the
  discrete study this matters more: a trait with an empty level changes what can be scored.
- Seeded BACE failures on fixed code are deterministic; retry only after a code change. MCMCglmm segfaults are
  not deterministic; rerun them.
- Request cluster memory from measured peak plus margin (40 GB blocked fir scheduling; 24 GB ran).
- `docs/` and `.unlazy/` are gitignored: `git add --force`, as a separate command (a failing plain add aborts an
  `&&` chain).
- The destructive-command hook blocks `rm -rf`; move to the scratchpad instead.
- Shipped prose: no em dashes; run `python3 ~/shinichi-brain/tools/slop_check.py <absolute path>`.
- Commit attribution line: `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>`.

## How to resume

Environment: R with BACE (installed), MCMCglmm, Rphylopars, glmmTMB, phylolm, nlme, ape, testthat. Safe check:
`Rscript -e 'testthat::test_dir("script/tests-rubin")'` from the worktree (about 5 minutes; expect FAIL 0).
Cluster status: `bash script/rubin_status.sh`. Never stage `.unlazy/imputation-sim/` or files outside this lane.

Read in order: this file, `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md` (v1's discrete results and
castor), `docs/dev-log/after-task/2026-09-25-rubin-campaign.md`,
`docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md`, then the repo `AGENTS.md`.

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-26-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```
