# Handover to Codex: pigauto, 2026-10-06

You are Codex, picking up pigauto from a Claude Code session (Opus 5.5, 2026-10-01 to 2026-10-06). Read
`AGENTS.md` first (it is native to you), then this file. Nothing below depends on that Claude chat.

## Critical context

- Repo: `itchyshin/pigauto` (R package). `main` is at `c9dc465` (PR #223 merged 2026-10-06 12:24 UTC). No
  open PRs. Version in DESCRIPTION: `0.11.0.9002`; CRAN has `0.10.0`.
- Your next job is the CRAN 0.11 release gate, approved by Shinichi in D-318 (2026-10-05): "CRAN 0.11.x
  after the bug-fix lane closes". The trigger is met: issues #196 to #215 are closed, #221 and #223 are merged.
  Shinichi chose Codex (not Cursor) for this (2026-10-06).
- Known blocker carried from D-315: `gllvmTMB` is in `Suggests` but is not on CRAN. The submission waits
  for it or drops it from `Suggests` (with the code paths guarded). Shinichi decides; bring him the two
  options with their costs.
- Submission itself is Shinichi's. You produce evidence and report the highest proven rung of the ladder in
  `~/shinichi-brain/protocols/cran-release-gate.md` (never "CRAN ready").

## Goals and plans

- Mission (AGENTS.md): mixed-type phylogenetic imputation with honest uncertainty; usability first
  (D-139). Preserve `r_cal = 0` as a valid fallback.
- Roadmap, most important first:
  1. CRAN 0.11 release (now).
  2. Follow-ups from #223 (below), smallest first.
  3. Ordinal accuracy against BACE at low phylogenetic signal (open research question).

## What was accomplished (this Claude session)

| PR | What | Status |
|---|---|---|
| #202, #203, #204 | Posterior MI: caveat numbers; Cholesky crash fix; `residual_prior = "sep"` default (fixes the lambda = 1, rho = 0.5 under-coverage: 0.750 to 0.965 at n = 1000) | merged |
| #222 | Docs: `log_transform = TRUE` re-logs already-logged traits (PanTHERIA slope shifts up to 3.9 SE; 1.8 with `FALSE`); opt-in `MI_REALDATA_LOG_TRANSFORM` in the real-data runner | merged |
| #223 | New defaults: `discrete_lambda = "estimate"` (binary, ordinal, zi-count gate and one-vs-rest categorical columns get their own Pagel lambda); `safety_floor = FALSE`, `phylo_signal_gate = FALSE` (decision D-320) | merged c9dc465 |
| #175 | GNN sentinel evidence PR | closed (branch kept) |

Evidence for #223: `docs/dev-log/discrete-lambda/README.md` (8 screens on Totoro, about 10,000 paired runs,
including the four-arm study's datasets paired with stored BACE results; `fisher-review.md` and
`rose-review.md` are the independent reviews). The shipped code reproduces the prototype on 5,475 datasets
exactly. Key numbers: discrete accuracy +0.04 to +0.12 at lambda <= 0.3 (discrete lambda alone about +0.01 to
+0.04); binary and categorical match or beat BACE; continuous error never worse; costs Brier +0.010 to +0.016 at
lambda near 1 with n >= 1000 and +0.014 to +0.020 on AVONET's categorical traits.

## Current working state

- Working: `main` passes `NOT_CRAN=true devtools::test()` (FAILED 0, ERRORS 0, PASSED 3312, SKIPPED 8 on
  the merged #223 branch) and `rcmdcheck --as-cran` (0 errors, 0 warnings, 1 NOTE: development version number;
  `gllvmTMB` not in mainstream repositories). CI green on Ubuntu release, Ubuntu devel, macOS.
- In progress: nothing. No runs on Totoro.
- Blocked: CRAN submission on the `gllvmTMB` decision.

## Key decisions (brain `~/shinichi-brain/memory/DECISIONS.md`)

- D-315 (2026-10-02): GNN off by default; `draws_method = "auto"`; gllvmTMB blocker noted.
- D-318 (2026-10-05): CRAN 0.11.x after the bug-fix lane; then the release gate.
- D-320 (2026-10-06): discrete lambda + gate/floor off by default; into 0.11 if merged before the gate (it was).

## Files created or modified (this session, on main)

- R: `R/joint_threshold_baseline.R`, `R/ovr_categorical.R`, `R/fit_baseline.R`, `R/fit_pigauto.R`,
  `R/impute.R`, `R/multi_impute.R`, `R/multi_impute_trees.R`, `R/mi_posterior.R`, `R/preprocess_traits.R`;
  `man/*.Rd` for those.
- Tests: `tests/testthat/test-discrete-lambda.R` (new), `test-lambda-per-type.R`, `test-lambda-dispatch.R`,
  `test-exact-default.R`, `test-gnn-off.R`, `test-phylo-signal-gate.R`, `test-mi-posterior.R`.
- Docs: `NEWS.md`, `vignettes/getting-started.Rmd` (+ `.R`), `vignettes/articles/simulation-study.Rmd`,
  `vignettes/multiple-imputation.Rmd`, `vignettes/common-pitfalls*`.
- Dev-log: `docs/dev-log/discrete-lambda/` (README, screens 1 to 8, reviews, scripts),
  `docs/dev-log/mi-posterior/sep_validation/`, `docs/dev-log/mi-posterior/pantheria_logtf/`.
- Real-data runner: `script/mi_realdata/01_run.R` (opt-in `MI_REALDATA_LOG_TRANSFORM`).
- This handover: `docs/dev-log/handover/2026-10-06-codex-handover.md` (branch `handover/2026-10-06-codex`).

## Branches left on origin (not merged, by design)

- `arc/rubin-freq-bace` (176 ahead): the Rubin freq-vs-BACE study, its report `script/rubin_study/study.qmd`
  and the campaign write-up `docs/dev-log/arc/2026-10-02-rubin-pigauto-campaign.md`. Long-lived study branch.
- `feat/discrete-lambda-ordinal` (8e99670): cumulative ordinal prototype, measured and NOT adopted.
- `feat/discrete-lambda-auto` (6969a20): validation-chosen lambda prototype, NOT adopted (has a known bug).
- `analysis/pantheria-logtransform-sensitivity`: content already on main via #222; safe to delete on
  Shinichi's word.

## Next immediate steps (yours)

1. Run the CRAN release gate (`~/shinichi-brain/protocols/cran-release-gate.md`; also the `cran-release-gate`
   skill if you have it): update profile, so a reverse-dependency check, a clean temporary-library install,
   examples, Windows vignette timing (win-builder), and a frozen tarball. Report the highest proven rung.
2. Bring Shinichi the `gllvmTMB` decision with options: (a) wait for gllvmTMB on CRAN; (b) drop it from
   `Suggests` and guard the code paths that use it (find them with `grep -rn gllvmTMB R/ tests/ vignettes/`).
3. Release version number and NEWS heading (0.11.0.9002 is a dev version): propose, do not decide.
4. After the release (or in parallel if it waits): the #223 follow-ups below.

## Follow-ups from #223 (open, small)

- The exact route's discrete accuracy under the new default: on the no-signal n = 60 fixture in
  `tests/testthat/test-exact-default.R`, 50 seeds give -0.014 (SE 0.010) against `per_column`. Not
  significant; worth a check on a real-size fixture.
- Categorical one-vs-rest lambdas are not reported in `$lambda_per_trait` (the rebuild is still deterministic).
- Ordinal accuracy vs BACE at low signal (0.38 vs 0.44 at n = 300, lambda = 0.3). Two ideas were tried and
  rejected (see the README, screens 5 and 7).

## Blockers and open questions

- `gllvmTMB` on CRAN (Shinichi's call).
- Nothing else blocks.

## Gotchas and failed approaches

- Totoro thread cap: set BOTH `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` (and `MKL_NUM_THREADS=1`).
  A smoke run with only OpenBLAS capped ran about 72 threads per R process (about 360 cores, over the 150-core
  cap, D-143) on 2026-10-05.
- Long jobs on Totoro must be detached (`setsid nohup ... < /dev/null &`); a job run inside an ssh command dies
  or hangs the session.
- `docs/` is git-ignored except dev-log files; add dev-log files with `git add -f`.
- The private R libraries on Totoro (`~/R/lib-dlam`, `lib-dlord`, `lib-dlauto`, `lib-dlfinal`, and `~/pigauto_rubin/rlib_sepf`) are
  prototype builds for the screens; do not use them for release checks.
- Never edit `BACE/` (separate package).

## How to resume (Codex)

From the repo root:

```bash
export NOT_CRAN=true OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
git fetch origin && git switch main && git pull --ff-only
Rscript -e 'devtools::test()'                          # expect FAILED 0
Rscript -e 'rcmdcheck::rcmdcheck(args = "--as-cran")'  # expect 0 errors, 0 warnings, 1 NOTE (gllvmTMB)
```

Live-toolchain work (real fits, `R CMD check`, win-builder, rendering) is yours. Totoro
(`snakagaw@totoro.biology.ualberta.ca`, ControlMaster socket `~/.ssh/cm-*totoro*`, at most 150 cores) for any
campaign; never GitHub Actions for campaigns. Any run over 3 hours needs Shinichi's approval first (D-287).
Before any public claim, run a fresh Rose audit.

Resume prompt to paste into Codex: `Rehydrate from docs/dev-log/handover/2026-10-06-codex-handover.md, then
continue with the Next Immediate Steps (CRAN 0.11 release gate).`

## Mission control

| repo | main / CI | what shipped this session | next, most important first |
|---|---|---|---|
| itchyshin/pigauto | c9dc465, CI green | #202 #203 #204 (posterior MI sep prior), #222 (log_transform docs), #223 (discrete lambda + gate/floor off) | CRAN 0.11 gate; gllvmTMB decision; #223 follow-ups; ordinal vs BACE |

Results page (private): https://claude.ai/artifact/Y644K1sLxbrWsAmTQvZspC ·
Mission Control: `~/shinichi-brain/Shinichi/Dashboards/mission-control/live/status/pigauto.json`.

## Full Codex prompt (added 2026-10-06 at Shinichi's request)

Paste this into a fresh Codex session started in the pigauto repo. It stands on its own.

```text
You are Codex, taking over the pigauto R package (GitHub: itchyshin/pigauto) from a Claude Code session that ran 2026-10-01 to 2026-10-06. You have none of that session's chat; everything you need is in the repo and the files named below.

READ FIRST, IN ORDER
1. AGENTS.md (repo root): project rules. Hard rules: preserve r_cal = 0 as a valid fallback; never edit BACE/; usability over contextualised accuracy (D-139).
2. docs/dev-log/handover/2026-10-06-codex-handover.md: the full handover (state, decisions, file list, gotchas).
3. ~/shinichi-brain/protocols/cran-release-gate.md: the release gate you will run (also the cran-release-gate skill if you have it).
4. docs/dev-log/discrete-lambda/README.md: evidence for the last big change (#223), only if you need it.

WHERE THINGS STAND
- main is at d99f9fa (handover doc, #224) or later; the last code change is c9dc465 (#223). No open PRs.
- DESCRIPTION: Version 0.11.0.9002 (development). CRAN currently has 0.10.0.
- On the #223 branch merged with main: NOT_CRAN=true devtools::test() gave FAILED 0, ERRORS 0, PASSED 3312, SKIPPED 8. rcmdcheck --as-cran gave 0 errors, 0 warnings and 1 NOTE (dev version number, plus "Suggests or Enhances not in mainstream repositories: gllvmTMB"). CI is green on Ubuntu R release, Ubuntu R devel and macOS.
- Nothing is running on Totoro.

WHAT SHIPPED RECENTLY (so the release notes are right)
- #223 (D-320): new defaults. discrete_lambda = "estimate", which gives binary, ordinal, zi-count gate and one-vs-rest categorical columns their own Pagel lambda. safety_floor = FALSE and phylo_signal_gate = FALSE. discrete_lambda = "fixed_1", safety_floor = TRUE, phylo_signal_gate = TRUE restores the old behaviour exactly.
- #204: posterior multiple imputation uses residual_prior = "sep" by default.
- #222: docs. log_transform = TRUE re-logs traits that are already on a log scale.
- #201 and #200: draws_method = "auto" and gnn = FALSE are the defaults.
- #216 to #221: the usability bug-fix lane (issues #196 to #215, all closed).
NEWS.md (section "# pigauto 0.11.0.9002 (dev)") already describes all of these. Check it reads as release notes before submission.

YOUR JOB: THE CRAN 0.11 RELEASE GATE (decision D-318, Shinichi, 2026-10-05; trigger met)
Follow cran-release-gate.md. It is fail-closed: the default verdict is NOT READY. Never say "CRAN ready". Report the highest rung of its ladder you have proven, and the next unproven rung.
This is an UPDATE submission, so the conditional gates apply:
  a. Reverse-dependency check of packages that depend on pigauto (revdepcheck or tools::package_dependencies on CRAN). Record the result even if there are none.
  b. Build one source tarball (R CMD build), record its sha256, and run every later check on that exact tarball.
  c. R CMD check --as-cran on the tarball in a clean temporary library. Examples and vignettes must run.
  d. Platform checks on the same source: win-builder (R-release and R-devel), and macOS builder or R-hub if available. Record vignette build time on Windows; incoming has a roughly 10-minute signal.
  e. Package size and timing budgets (tarball size, per-vignette time, test time under NOT_CRAN unset).
  f. Refresh the current CRAN Repository Policy and submission checklist from cran.r-project.org before judging any threshold, and label each threshold as current policy, observed incoming behaviour, or our own margin.
  g. Check DESCRIPTION (Title case, Description wording and spelling, Authors@R, URL and BugReports, License), cran-comments.md, and inst/WORDLIST / spelling.
Write all evidence to docs/dev-log/cran-0.11/ (one README with the ladder, plus raw logs), on a branch named release/cran-0.11-gate. Open a PR. Do not merge it.

DECISIONS YOU MUST BRING TO SHINICHI (do not make them)
1. gllvmTMB is in Suggests but not on CRAN (carried from D-315). Option (a): wait until gllvmTMB is on CRAN. Option (b): drop it from Suggests, and guard or remove every use (find them with: grep -rn gllvmTMB R/ tests/ vignettes/ man/ DESCRIPTION). For each option give the cost: what users lose, files touched, and whether --as-cran becomes NOTE-free. If he picks (b), implement it on the release branch with tests.
2. The release version number. Propose 0.11.0 or 0.11.1 with a reason (check NEWS for any 0.11.0 heading). Do not change DESCRIPTION's version until he decides.
3. Submission itself. You never submit; you stop at "submission-ready" or lower.

FOLLOW-UPS FROM #223 (smaller; do them after the gate, or while waiting on win-builder)
- The exact route's discrete accuracy under the new default. On the no-signal n = 60 fixture in tests/testthat/test-exact-default.R, 50 seeds give -0.014 (SE 0.010) against per_column. Measure on a realistic fixture (n = 300, lambda 0.3 and 1, about 50 seeds) and report. Change code only if there is a real, significant loss, and then as a PR.
- Categorical one-vs-rest lambdas are not reported in $lambda_per_trait (binary and ordinal are). Add them in a backward-compatible way, with a test.
- Ordinal accuracy against BACE at low phylogenetic signal (0.38 vs 0.44 at n = 300, lambda = 0.3). Research only. Two ideas were already tried and rejected (README screens 5 and 7, branches feat/discrete-lambda-ordinal and feat/discrete-lambda-auto); do not repeat them.

ENVIRONMENT (you run the live toolchain)
export NOT_CRAN=true OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
git fetch origin && git switch main && git pull --ff-only
Rscript -e 'devtools::test()'                          # expect FAILED 0
Rscript -e 'rcmdcheck::rcmdcheck(args = "--as-cran")'  # expect 0 errors, 0 warnings, 1 NOTE (gllvmTMB)

COMPUTE RULES (binding)
- Heavy runs go on Totoro: snakagaw@totoro.biology.ualberta.ca, using the existing ControlMaster socket ~/.ssh/cm-*totoro*. Never trigger Duo. Use at most 150 cores in total (D-143).
- Cap BOTH OpenBLAS and OpenMP threads on every R process. Capping only OpenBLAS ran about 72 threads per process (about 360 cores) on 2026-10-05.
- Detach long jobs: setsid nohup ... < /dev/null &. A job inside an ssh command dies or hangs.
- State a time estimate before any run. Anything over 3 hours needs Shinichi's approval first (D-287). A run that overruns its estimate by 30% stops and is re-reported.
- Never use GitHub Actions for campaigns. Never compute on DRAC login nodes.
- The private Totoro libraries ~/R/lib-dlam, lib-dlord, lib-dlauto, lib-dlfinal and ~/pigauto_rubin/rlib_sepf are prototype builds. Never use them for release checks.

REPO AND GIT RULES
- Never commit to main; work on branches and open PRs; never auto-merge. Shinichi merges, or says "merge #N".
- docs/ is git-ignored except dev-log files: add those with git add -f.
- Stage explicit paths only. Never stage files you did not create (other lanes' untracked files exist in some checkouts).
- Before any PR body, comment or public claim: run python3 ~/shinichi-brain/tools/agent_mention_check.py --text <file>. Never @-mention an agent name. Shipped prose: no em dashes; run python3 ~/shinichi-brain/tools/slop_check.py <file>.
- Run a fresh Rose audit (.codex/agents if present; otherwise an independent reviewer pass) before any public claim, for example before saying the release gate reached a rung.

BRANCHES LEFT ON ORIGIN (leave them alone unless told)
- arc/rubin-freq-bace: the Rubin freq-vs-BACE study. Long-lived and unmerged by design.
- feat/discrete-lambda-ordinal and feat/discrete-lambda-auto: rejected prototypes, kept as a record.
- analysis/pantheria-logtransform-sensitivity: already on main via #222. Delete only when Shinichi says so.

STATUS REPORTING
- Update Mission Control at milestones: ~/shinichi-brain/Shinichi/Dashboards/mission-control/live/status/pigauto.json (now.focus, now.next_safe_action, now.active_lane). Commit only that file in the vault, which is local-only with no remote.
- Do not write elsewhere in ~/shinichi-brain without Shinichi's approval.
- When you report to Shinichi, put the answer first in plain language: the ladder rung reached, what blocks the next rung, and the one decision you need from him.

DONE MEANS
The release-gate PR is open, with docs/dev-log/cran-0.11/README.md showing every gate's evidence on one frozen tarball. The highest proven rung is stated. The gllvmTMB and version decisions are put to Shinichi with options. Mission Control is updated. Nothing is submitted.
```
