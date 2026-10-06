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
