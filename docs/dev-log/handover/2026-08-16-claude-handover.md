# Session Handoff: pigauto — Tier 0 + Tier 1 COMPLETE, two campaigns landed (Claude → Claude)

Meta: 2026-08-16 ~08:20 MDT · from Claude (`AUTHOR = claude`) · `TARGET = claude`
Authoritative copy: this file. Chat is disposable. You inherit **no** authoring chat.

You are **Claude**, picking up **pigauto only**. A long overnight lane closed Tier 0, merged
9 PRs, re-ran the whole bench suite, and completed both Tier-1 Totoro campaigns. **Nothing is
running.** Nothing is mid-flight. Do not re-merge, do not relaunch, do not CRAN-submit.

**This supersedes `2026-08-16-cursor-handover.md`** — Shinichi decided to stay on Claude, so
no Cursor lane was started. That file remains accurate about repo state; use this one.

## Critical Context

1. **`origin/main` = `175ebdc`**, Version `0.10.0.9000`, Suggests `BACE`. Nine PRs merged
   (#158–#166). CI green throughout. `--as-cran`: 0 errors / 0 warnings / 1 known note.
2. **Both Totoro campaigns are COMPLETE, zero failures.** Regime map 5,400 jobs
   (`DONE ok=5400`), coverage 1,920 jobs (`DONE wall=2152s ok=1920 skip=0 fail=0`). Results
   are summarised, pulled back, committed. **Nothing to resume.**
3. **Two new findings that constrain what the paper may claim** — read both before writing
   any prose about the GNN or about uncertainty:
   - `docs/dev-log/2026-08-16-regime-map-results.md`
   - `docs/dev-log/2026-08-16-coverage-results.md`
4. **Bench suite re-run complete, 8/8.** The training leak is undetectable in every trait
   type (`2026-08-16-bench-rerun-results.md`).
5. **PROTECTED / do not touch:** standalone BACE `@ce8bc87`; in-tree `pigauto/BACE`
   `@de87d8c`; the 15 dirty carried-over files; ~20 foreign unpushed branches; EM
   (`max_iter > 0`); DRM.jl; GLLVM.jl (Shinichi said leave it alone). Never `git add -A`.
   Never archive this project.

## The two findings, in one paragraph each

**Regime map (#135).** At λ = 1, where pigauto's BM baseline is correctly specified, the
linear specificity control **F1 is null** (mean |Δ|/MCSE 1.51, 11% of cells detectable) while
the nonlinear family **F2 is not** (3.29, 58%). So the GNN recovers structure the joint-MVN
baseline structurally cannot represent, and adds nothing where it can. That is the evidence
#135 was opened to get. **Separately — do not merge the two —** the much larger gains at
λ < 1 are the baseline's λ = 1 assumption being violated, not a second win; the effect shows
up in the *linear* family, which identifies it as misspecification. pigauto already ships a
Pagel-λ baseline, and **this campaign did not test whether fitting λ in the baseline recovers
the same gain**. No low-λ claim should ship until it does. That is the highest-value
follow-up in the repo right now.

**Coverage.** Conformal intervals **undercover at n = 100** (0.853–0.861 vs nominal 0.95,
6/6 cells inadequate) and are well calibrated at n ≥ 300 (0.948–0.963, 6/6 adequate).
**See the CORRECTION at the foot of that doc before using this paragraph** -- the campaign
was MCAR-only, and the repo's real-data benches show 0.87-0.91 coverage at n up to 10,654,
which the arithmetic mechanism below cannot explain (a second, likely-exchangeability
defect). `conformal_split_val = TRUE` does **not** repair it, which rules out the documented
gate/conformal-reuse explanation. The mechanism is arithmetic: split conformal needs
n_val ≥ 19 for 95%, and n = 100 with default splits gives n_val ≈ 10, where the ceiling is
n_val/(n_val+1) ≈ 0.909 — **95% is unreachable before any noise**. This also retires the old
0.884–0.887 number, which was dismissed as leak-tainted and was in fact measuring this.
`DESCRIPTION`, `README` and `gnn-architecture.Rmd` §7 all present ≥95% coverage as the
guarantee. **Same defect class as #157.**

## What Was Accomplished

| PR | Item |
|---|---|
| #158 | P0 honesty fixes + two-layer #157 gate floor + 9-edit doc bundle (**closes #157**) |
| #159 | pkgdown destination fix (post-merge fallout) |
| #160 | P1-7, P1-6 |
| #161 | P1-13 + P1-10 fail-closed guards |
| #162 | P1-5 joint-solver attribution |
| #163 | P1-12 multi-obs phylo-signal gate (silent no-op) |
| #164 | P1-9 zi_count zeros in splits (silent) |
| #165 | P1-11 documented, deliberately not fixed |
| #166 | all 8 bench drivers leaked their PSOCK cluster on error (found by the drmTMB lane) |

Plus: Rose P1 queue **8/9 addressed**; brain lessons filed (vault `0cf87d1`); the
brain-write rule's false D-37 citation corrected in 5 repos on Shinichi's instruction.

## Current Working State

- **Working:** everything. `main` green, both campaigns complete, bench suite 8/8.
- **In progress:** **nothing.** No background jobs, no Totoro jobs, no open PRs.
- **Blocked:** CRAN (Suggests BACE; CRAN BACE 404). P1-8 and the brain-write-policy question
  need Shinichi.

## Landing State

`handoff_gate.sh` → FAIL: 15 uncommitted on this branch + ~40 unpushed commits on ~20 other
branches + a stale `index.lock`. All declared; nothing silently dropped.

| Artifact | Committed | Pushed | State |
|---|---|---|---|
| `origin/main` `175ebdc` | y | y | **LANDED.** 9 PRs. Do not re-merge. |
| `origin/handover/2026-08-09-cursor` (all docs) | y | y | **LANDED.** No PR into main — do not open one. |
| Totoro `~/pigauto_regime_map/` | n/a | n/a | **COMPLETE.** 5,400 + 1,920 RDS + summaries retained on Totoro. |
| GLLVM.jl `fix/brain-write-citation` | y | y | **CARRIED-OVER.** Unmerged by Shinichi's instruction — the rule is not on GLLVM.jl's `main`, so merging would *add* a rule, not fix one. |
| 15 dirty carried-over files | n | n | **CARRIED-OVER — do not stage.** Incl. `AGENTS.md`, whose brain-write line is corrected but uncommitted. |
| ~40 unpushed commits, ~20 foreign branches | mixed | n | **CARRIED-OVER.** Not this lane. |
| Ghost `D`/`MM` entries under `docs/` in `git status` | — | — | **ARTEFACT** of the `GIT_INDEX_FILE` workaround; files verified on disk and origin. `git reset` clears them once the lock is gone. |
| stale `.git/index.lock` | — | — | **REPORT-ONLY.** Cursor gitWorker artefact, `lsof`-verified unheld; harness blocks `rm`. |

## Next Immediate Steps (OWED only)

Classify against **current** git/gh before acting.

1. **Rehydrate.** `~/shinichi-brain/tools/lane_preflight.sh "<repo>"`; `git fetch origin`;
   `git status`; read this file + the two results docs + the coordination board (**both**
   Active Lane Split rows). Expect `main-direct` in the census — that is this lane's own PR
   merges under Shinichi's GitHub identity, not a foreign agent.
2. **Nothing is owed automatically.** Both campaigns are done and written up. The next move
   is Shinichi's, from this menu:
   - **The Pagel-λ comparison** (highest value — it decides whether the low-λ regime-map
     gains survive a correctly-specified baseline).
   - **One working external comparison (arguably now the top item).** Every BACE
     head-to-head in `script/` reports "**BACE skipped** (not installed or failed)" --
     `bench_bace_avonet_head_to_head.md`, `bench_pantheria_bace_head_to_head.md`,
     `bench_avonet_bace.md`. Rphylopars is a *dependency* (the joint-MVN solver), not a
     benchmarked comparator. **The repo currently has no working comparison against any
     external package**, so no comparative claim can ship. Rphylopars standalone on
     continuous traits is the cheapest credible comparator.
   - **The real-data coverage defect** (see the CORRECTION): run
     `bench_missingness_mechanism.R` against the exchangeability hypothesis. The n=100
     arithmetic fix (warn at `n_val < 19`) is a one-line honesty fix, not the priority.
   - **P1-8** — covariates ignored by the joint/threshold-joint baseline; a real feature with
     a design choice in it.
   - **Paper framing** — Tier 1.2 said the choice waits for the regime map and must not be
     made retroactively by it. The evidence is now in.
3. **Before any public claim**, spawn the Rose lens (repo rule) and re-read
   `2026-08-16-regime-map-results.md` §"What this licenses" — it states explicitly what the
   data does and does not support.

## Gotchas & Failed Approaches (this lane paid for all of these)

- **Never `renice` a process that will fork PSOCK workers** — niceness is inherited at fork;
  workers miss the connect window. Cost ~2 h. Bound cluster load by *worker count*.
- **`parLapply` logs nothing until every cell finishes.** A master at 0% CPU is normal, not
  hung. Killed two healthy drivers on this.
- **torch hangs the R process AFTER work is written.** Drive chains off the driver's own
  completion marker, not process exit.
- **The bench drivers resume from their own output `.rds`.** A stale copy makes a "re-run"
  compare old numbers to themselves — it logs `Running 0 cells` and `0.0s`. Always verify
  non-zero cells and non-trivial wall time. This bit once (`multi_proportion`).
- **`git switch -c` in a shared checkout moves HEAD under whoever is standing in it.** Did
  this to GLLVM.jl and CBIC; caught by a peer lane and repaired. Use a temp `git worktree`.
- **Price the JOB, not the FIT.** The D-139 pre-run timed `impute()` only; the campaign ran
  ~2.5 h against a 1.5 h estimate.
- **Shell parsing lied about campaign health four times** — multi-line `ssh` output collapsed
  into wrong counts, twice reporting a healthy campaign as failed/finished. Verify state with
  a direct, single-value query before believing a monitor.
- `docs/` is `.gitignore`d — `git add -f` specific paths only.

## How to Resume (Claude)

Working directory `/Users/z3437171/Dropbox/Github Local/pigauto`, branch
`handover/2026-08-09-cursor` for docs; code is `origin/main` `175ebdc`.
Toolchain: R + `devtools` + `torch`; `NOT_CRAN=true` for tests. Totoro reachable via the
standing ControlMaster socket; campaigns live in `~/pigauto_regime_map/`.

Safe verification (~1 min, touches nothing):
```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git fetch origin && git log -1 --oneline origin/main      # expect 175ebdc
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-type-detection.R")'
```

Housekeeping only Shinichi can do (harness blocks agent `rm` on `.git`):
```bash
rm "/Users/z3437171/Dropbox/Github Local/pigauto/.git/index.lock"
git -C "/Users/z3437171/Dropbox/Github Local/pigauto" reset
```

Read order: `AGENTS.md` (**note: uncommitted and differs from every committed ref**) → this
file → `2026-08-16-regime-map-results.md` → `2026-08-16-coverage-results.md` →
`2026-08-16-overnight-report.md` → `docs/dev-log/coordination-board.md` →
`docs/dev-log/2026-08-15-publication-checklist.md`.

Rehydrate (human starts a **fresh Claude** session in the repo and pastes):

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-16-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

Interactive equivalent: `claude "Rehydrate from docs/dev-log/handover/2026-08-16-claude-handover.md, then continue with the Next Immediate Steps."`

## Mission control

| | |
|---|---|
| Repo | `/Users/z3437171/Dropbox/Github Local/pigauto` |
| main | `175ebdc` · `0.10.0.9000` · Suggests `BACE` |
| Merged | PRs #158–#166 (9) · CI green · `--as-cran` 0E/0W/1N |
| Campaigns | regime map 5,400 ✅ · coverage 1,920 ✅ · **0 failures** · summarised + committed |
| Benches | 8/8 re-run; leak undetectable in every type |
| Rose P1 | 8/9 addressed; **P1-8 open** |
| Open for Shinichi | Pagel-λ comparison · coverage remedy · P1-8 · paper framing · brain-write policy |
| Do not | re-merge · CRAN-submit · relaunch campaigns · touch GLLVM.jl · `git add -A` · archive |

> Related: `docs/dev-log/2026-08-16-regime-map-results.md` ·
> `docs/dev-log/2026-08-16-coverage-results.md` ·
> `docs/dev-log/2026-08-16-bench-rerun-results.md` ·
> `docs/dev-log/2026-08-16-overnight-report.md` ·
> `docs/dev-log/2026-08-15-publication-checklist.md` · `docs/dev-log/coordination-board.md`
