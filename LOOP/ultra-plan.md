```
🎯 GOAL
Solo platform: Cursor (this session; read, not inferred — R/compute in scope does NOT mean Codex)
Deliverable: An approved G0 decision on pigauto's BACE wrap API (`bace_imp` chain averages vs
  `bace_imp` + `bace_final_imp` proper MI), plus a runnable Phase-3 plan scoped to
  `R/fit_baseline_bace.R` + tests + NEWS, executed later via `/goal` — not in this chat.
HEADLINE: Decide the wrap API. Everything else in this lane (re-bench vs BACE@ce8bc87, OVR
  footnote, snapshot refresh) is downstream of that one choice.
IN PARALLEL: OVR naming footnote (docs-only, no estimator change); standalone BACE @ce8bc87
  install check; coordination-board status line.
DEFER (fenced — NOT in this plan):
  · P0 / PR #155 merge work (sibling lane, do not redo)
  · vendor-sync of in-tree `pigauto/BACE` @de87d8c (PROTECTED)
  · any edit to standalone BACE @ce8bc87 (PROTECTED)
  · DESCRIPTION joint-Σ / 95%-coverage claim gate (P1, other lane)
  · EM restore (`max_iter > 0`) — forbidden without a new G0
  · rebase onto `main`; merge anything to `main`
  · covariates into joint baseline (P1-8); DRAC GNN attribution ladder
DISCIPLINE: verify=no wrap claim without a BACE-installed run of the T4 smoke + a re-bench vs
  @ce8bc87 · compute=local (wrap smoke) then Totoro if a full re-bench ladder is approved ·
  closure=Shinichi answers the single Phase-0.4 question below; plan then moves to `/goal`.
```

Plan artifact · 2026-08-09 · author Ada (Cursor) · Phases 0–2 read-only, then **G0 LOCKED**.

---

## ✅ G0 LOCKED — 2026-08-09

**WHAT SHINICHI TOLD US** (verbatim): *"it should be final yes — OK opt-in cool /ultra-plan it"*.

Locked decisions:

| Item | Locked value | Provenance |
|---|---|---|
| Wrap API | **Option B-minus.** `fit_baseline_bace()` gains `final_imp = FALSE` (default) + `n_final`. When `final_imp = TRUE`, call `BACE::bace_final_imp()` on the existing `bace_imp` object and derive `mu`/`se` from `all_datasets`. **Default path must stay byte-identical.** | Shinichi, approved |
| `n_final` default | **15L** (matches Dan's Study B). | **Ada's IF-YOU-DO-NOT-MIND default — Shinichi did not choose between 15 and 50.** Do **not** flip to 50 without a new ask. |
| Slack to Dan | **Parked — not requested.** Do not draft or send. G0 is sufficient. | Shinichi, by omission |
| OVR footnote (row S2) | In scope, docs-only. | same G0 |

Phase 3 executes from here. Everything in the DEFER fence above remains fenced.

---

## Phase 0 — Orient

Lane: **pigauto BACE-wrap / re-bench**, branch `handover/2026-08-09-cursor` @ `8b5037f`.
Sibling lane: P0 / PR [#155](https://github.com/itchyshin/pigauto/pull/155) — do not touch.

`python3 ~/shinichi-brain/tools/route.py pigauto` LOAD-FIRST loaded. Relevant invariants for this
lane: *"diff main before building"*, *"prediction-path correctness is the first audit lane"*,
*"keep uncertainty propagation first-class"*. The last of those is what makes the wrap's `se`
semantics the substantive part of this decision, not the API cosmetics.

Repo `AGENTS.md` rule 3 ("BM baseline is a kernel, not a method… scope the change as *add a
parameter, default off* rather than *replace*") is directly load-bearing on the recommendation.

---

## Phase 0.25 — Prior-work sweep (RECEIPT — hard gate, evidence cited)

### Surface 1 — repo git

| Probe | Result |
|---|---|
| `git status -sb` | `## handover/2026-08-09-cursor...origin/handover/2026-08-09-cursor` — in sync. 6 modified (uinit/banners/LOOP), 10 untracked (DRAC GNN attribution). **All left unstaged**, per handover. |
| `git log --oneline -20` | HEAD `8b5037f docs: Cursor handover…`; below it `5800b01` merge-from-parent, `0ed6a73` hygiene, `58708ac` P0 blockers — i.e. this branch already carries the P0 tip. |
| `git branch -a` | ~90 local + ~70 remote branches. **No branch matching `bace`/`wrap`/`final-imp` exists** — no prior wrap attempt to rebuild. |
| `git worktree list` | 17 worktrees (2 prunable under `/private/tmp`). **None** on a BACE/wrap branch. |
| `git stash list` | 1 stash, `stash@{0}: On codex/portable-benchmark-paths: agents/claude wip` — unrelated lane, not touched. |
| `branch_drift_check.sh` | `handover/2026-08-09-cursor vs origin/main: 9 ahead, 43 behind ⚠ DRIFTED`. Rebase is **forbidden** this lane; mitigation is a paths-scoped diff (see MECHANICAL-VERIFY row M1). |
| `lane_preflight.sh` | `no codex lane detected in the last 12h`. Per D-87 this is **weak** evidence, not proof of sole ownership. Open PRs: only #155. |
| `gh pr view 155` | **OPEN**, base `fix/ci-install-libtorch`, head `fix/p0-review-blockers`, `MERGEABLE` / `UNSTABLE`, updated 14:08 UTC. |
| `gh pr checks 155` | `ubuntu-latest (release)` **pending**; `pkgdown` skipping. **CI has not concluded** → merge acceptance not met → not merged by me. |

### Surface 2 — sister packages (read-only, both PROTECTED)

| Clone | HEAD | State |
|---|---|---|
| `/Users/z3437171/Dropbox/Github Local/BACE` | `ce8bc87` | `main...origin/main` in sync. `R/` includes `bace.R`, `bace_imp.R`, `bace_final_imp.R`, `pool_mi.R`, `with_imputations.R`, `accessors.R`, `phylo_signal_summary.R`. |
| `/Users/z3437171/Dropbox/Github Local/pigauto/BACE` | `de87d8c` | `main...origin/main [behind 86]`. Untouched. |

Signatures read (read-only, no edits):

- `bace(fixformula, ran_phylo_form, phylo, data, nitt=6000, thin=5, burnin=1000, runs=10,`
  `n_final=50, …, ovr_categorical=FALSE, …)` — chain **plus** final phase in one call.
- `bace_final_imp(bace_object, fixformula, ran_phylo_form, phylo, nitt=6000, thin=5, burnin=1000,`
  `n_final=50, …, ovr_categorical=FALSE, …)` — **consumes a `bace` object produced by `bace_imp`**.
  Returns `class "bace_final"` with `all_datasets` = list of length `n_final`.
- `bace_imp(…, runs=10, …, ovr_categorical=FALSE, …)` — the chain phase pigauto already calls.

**Load-bearing consequence:** because `bace_final_imp()` takes the `bace_imp()` object as input,
Option B is **additive** — pigauto keeps its existing chain call and appends a final phase. It is
not a rewrite of the wrap.

BACE roxygen on `n_final` (quoted): default 50 *"large enough that the empirical 2.5 and 97.5
percent quantiles of per-cell imputed values give a reasonable estimate of the 95 percent
posterior predictive interval. Smaller values undercover: n_final = 5 gives roughly 67 percent
effective coverage, n_final = 20 …"*. Dan's Study B ran `n_final = 15`.

### Surface 3 — pigauto wrap + its consumers

- `R/fit_baseline_bace.R:92` calls `BACE::bace_imp()` only (`runs=5, nitt=4000, burnin=1000, thin=10`).
  Lines 104–107 drop `Initial_Data` and take `fit$data[2:(n_iter+1)]`; lines 165–177 set
  `mu` = mean over those datasets and `se` = between-dataset `sd`.
- Exported (`NAMESPACE:29`), documented (`man/fit_baseline_bace.Rd`), listed on the pkgdown
  reference index (`_pkgdown.yml:147`). So the API is public — a signature change needs NEWS.
- Test surface is thin: **one** smoke, `test-shipping-coverage.R:288` `[T4]`, guarded by
  `skip_if_not_installed("BACE")` **and** a `tryCatch → skip`. It asserts only `list` / `mu` /
  `se` / `nrow`. `test-bace-compat-eval.R` (222 lines) tests `script/gha/_bace_compat.R` eval
  plumbing, **not** the wrap.
- Head-to-head benches call `BACE::bace()` (not the wrap): `script/bench_avonet_bace.R:219`,
  `bench_bace_avonet_head_to_head.R:115`, `bench_pantheria_bace_head_to_head.R:156`,
  `bench_fishbase.R:317`.

### Surface 4 — the environment (NEW finding, not in the handover)

```
Rscript -e 'requireNamespace("BACE", quietly=TRUE)'  →  FALSE
```

**BACE is not installed on this machine at all.** The handover's open question was "installed BACE
may not be `@ce8bc87`"; the measured truth is stronger — it is **absent**. This explains, without
appeal to a code defect, why the T4 smoke skips and why all three local benches printed
*"BACE skipped (not installed or failed)"*. Nothing in this lane — status quo or change — can be
verified locally until standalone `@ce8bc87` is installed. That is Phase 3 step zero.

### Surface 5 — brain

MCP `search_notes` (`search_all_projects: true`, never `project:`) — D-82 walk not needed, MCP
answered on the first call.

- Query *"pigauto BACE wrap bace_final_imp baseline"* → top hits are indexed **BACE source files**
  (`bace_final_imp.R`, `.Rd`, `demo_bace_wrapper.R`) plus two design notes:
  `2026-04-20-phase-8-benchmarks-design` (*"`bace_default`: `BACE::bace()` with OVR enabled"*) and
  `2026-05-16-bace-headtohead-ci-design` (*"`BACE/` — the vendored sister package, never touched by
  pigauto work"*; 6 datasets, 30% MCAR, fixed seed).
- Query *"pigauto invariants WHAT-WORKS MODEL-ROUTING WORKING-STYLE OPEN-LOOPS"* → `59-missing-data-layer`
  (two-path strategy: gllvmTMB `with_imputations()` vs pigauto for rich GNN / mixed types / tree
  uncertainty), `71-pigauto-mi-sister-path`, `pigauto_state_and_next_steps`.
- Deterministic vault greps: `memory/AGENT_LOG.md` — **no** `bace` hits. `memory/DECISIONS.md` —
  **no** `bace` hits. `memory/OPEN_QUESTIONS.md` — two hits, both the *"`BACE` is missing — add it"*
  web-page/repo-listing item, not a wrap decision. `journal/` — one hit, same web-page item.
  `projects/deep-research/README.md` — none.

**Sweep verdict:** there is **no prior wrap attempt** anywhere — no branch, no worktree, no stash,
no brain decision. The two brain design notes cover **bench** design (`BACE::bace()` as a *rival*),
not the wrap-as-baseline API. This G0 is genuinely unresolved, and the decision, once made, is a
durable one that belongs in `memory/DECISIONS.md` (with Shinichi's approval — D-37 boundary).

---

## Phase 0.3 — Routing

Per `MODEL-ROUTING` / the Cursor adapter: scout and mechanical sweep on **Cursor Models**
(Composer 2.5 / Grok 4.5); judgment, synthesis and adversarial review on **Other Models**
(Auto Cost or pinned Claude/GPT); live R, `devtools::check()`, BACE install and any re-bench
**hand off** to Codex — those are toolchain work, not Cursor-parent work.

### Phase 0.3b — two-bar usage

**UNVERIFIED.** Cursor Settings → Usage is not readable from this session (no shell/MCP surface
exposes the two-bar figures). Per instruction this does not block. Owner heuristic stands: ~3%/day
average on the Cursor Models bar; on-demand disabled. Recommend Shinichi glances the bars before
approving the Phase-3 route, since row S3 below is deliberately parked on Codex rather than either
Cursor bar.

---

## Phase 0.4 — G0 (ONE question, homework already done)

**QUESTION.** Should `fit_baseline_bace()` stay on cheap `bace_imp()` chain averages, or gain an
opt-in `final_imp = TRUE` / `n_final` path that appends `BACE::bace_final_imp()` and derives
`mu` / `se` from proper MI draws — default **off** in v0.10, flipped only after a measured
re-bench against `daniel1noble/BACE@ce8bc87`?

**WHY NOW.** Dan published Study B (~10k sims, 9,997 ok) on the `bace()` + `bace_final_imp` path
on 8 Aug. pigauto's wrap has never called that path. Until the wrap can reach it, no pigauto
sentence may cite Study B, and the wrap's `se` has no documented coverage property. The sweep also
found BACE is **not installed here**, so this decision determines what we install and re-bench
against — answering it first avoids a wasted install-and-measure cycle.

**TEAM VIEW.**

- **Ada (orchestration).** Additive, not replacement. `bace_final_imp()` consumes the object
  `bace_imp()` already returns, and its `all_datasets` slot has the same shape as the
  `fit$data[-1]` list the wrap already re-encodes — so the code delta is one new branch plus an
  argument, not a rewrite. `AGENTS.md` rule 3 says add-a-parameter-default-off.
- **Rose (adversarial / claim gate).** Two traps. (1) **Do not double-count Dan's fix.** His
  imputed-as-observed bug was in `bace_final_imp`; `bace_imp` already reinstated response NAs.
  pigauto never called the buggy function, so moving to it repairs **nothing** in pigauto's current
  output — NEWS must not say "fixes imputed-as-observed". What Option B buys is MI propriety and
  Study-B alignment pigauto *lacks*, which is a different and weaker claim. (2) Any wrap change
  changes `se`, and `se` is a **public** output feeding `fit_pigauto(baseline=)` and MI draws —
  so "no behaviour change" is false unless the new path is default-off.
- **Fisher/Gauss (inference).** The status quo `se` is the SD across `runs` **successive sweeps of
  a chained-equations chain**. Those are autocorrelated by construction and include convergence
  transient; they are not independent posterior-predictive draws, so that SD is not an imputation
  standard error and has no coverage guarantee. `bace_final_imp`'s `n_final` draws are the object
  with the documented interval property (and BACE's own roxygen warns `n_final = 5` ≈ 67%
  effective coverage). Expect the new `se` to be **larger** — that is the correction, not a
  regression, but it will move every downstream interval.

**RECOMMENDATION.** **Option B-minus** — add `final_imp = FALSE` (default) + `n_final = 15L` to
`fit_baseline_bace()`; when `TRUE`, append `BACE::bace_final_imp()` and compute `mu` / `se` from
`all_datasets`. Ship default-off in v0.10 with a NEWS note that states the propriety gap honestly
and does **not** claim a bug fix. Flip the default only after row S3's re-bench. Rationale: it
makes the Study-B-aligned object reachable (unblocking every future pigauto-vs-BACE sentence),
costs nothing to existing users, respects the add-a-parameter rule, and keeps `r_cal = 0`-style
fallback logic — the cheap path remains a valid choice for small compute budgets.

**IF YOU DO NOT MIND.** Two sub-answers would sharpen Phase 3: (a) `n_final` default — 15 (matches
Dan's Study B, cheap) or 50 (BACE's own default, better tails)? (b) Do you want the single Slack
to Dan sent, or is your G0 sufficient? The handover says do **not** ask him to re-run 10k; the only
question worth sending is the wrap-API one.

**WHAT CONTINUES REGARDLESS.** The docs-only OVR footnote (row S2) and the coordination-board
status line are unblocked by this answer. The BACE `@ce8bc87` install (row S0) is required under
*both* options — status quo cannot be verified either, since BACE is absent.

---

## Phase 0.5 — NotebookLM

Offered once, not run. A NotebookLM pass over Dan's `simulation-report.html` +
`.agents/investigation-2026-08.md` could give a grounded, citable summary of Study B's coverage
columns for the eventual NEWS / manuscript wording. **This is not a novelty or priority claim
against the literature**, and the G0 decision is **not** gated on it. Say the word and it runs;
otherwise it stays parked.

---

## Phase 1 — Decomposition

The wrap change is small and well-fenced. The expensive, uncertain part is measurement, and it is
blocked on an install that does not exist yet. So the lane decomposes as:
**install → recon the delta → implement behind a default-off flag → mechanically verify →
measure → reconcile.** Rows S3+ only begin if S0 succeeds.

Scope (if approved): `R/fit_baseline_bace.R`, `man/fit_baseline_bace.Rd` (regenerated),
`tests/testthat/` (new focused file), `NEWS.md`.
Out of scope: BACE source (either clone), vendor-sync, `DESCRIPTION`, EM, GNN, DRAC, PR #155.

---

## Phase 2 — Runnable plan

| # | Slice | Bar | Model / mode | Depends | Done when |
|---|---|---|---|---|---|
| **S0** | Install standalone BACE `@ce8bc87` (`devtools::install("/Users/z3437171/Dropbox/Github Local/BACE")`) — **never** the in-tree clone. Record `packageVersion` + SHA. | **hand off** | Codex Sol | — | `requireNamespace("BACE")` TRUE **and** recorded SHA == `ce8bc87`. |
| **R1** | **RECON.** Read-only: exact shape of `bace_final_imp()$all_datasets` vs `bace_imp()$data[-1]`; which columns come back as factor vs numeric; whether the wrap's existing re-encode loop (lines 109–163) works unchanged on `all_datasets`. Note any `ovr_categorical` / `nitt_cat_mult` args the wrap would need to expose. **No edits.** | Cursor Models | Composer 2.5 (Agent) | S0 | A written delta list: what changes in the wrap, what does not. |
| **S1** | Implement `final_imp = FALSE` + `n_final = 15L` in `fit_baseline_bace()`; when TRUE, call `BACE::bace_final_imp()` on the existing `bace_imp` object and derive `mu`/`se` from `all_datasets`. Roxygen: document the propriety difference and the runtime cost. `devtools::document()`. | Cursor Models | Composer 2.5 (Agent) | R1 | Both paths run; default path is byte-identical to today. |
| **S2** | **Docs-only, parallel.** `R/ovr_categorical.R` comment/roxygen: pigauto OVR is no longer "the OVR strategy BACE uses" as BACE's *current default* (BACE default is now multinomial). Footnote only — pigauto's reason (Rphylopars full-rank) is unchanged. **Not an estimator change.** | Other Models | Auto Cost | — | Wording states the divergence without implying pigauto should follow. |
| **S3** | Re-bench wrap (both paths) vs `BACE::bace()` @ `ce8bc87`. Reuse `script/bench_avonet_bace.R` / `bench_bace_avonet_head_to_head.R`. Report RMSE/accuracy **and** interval coverage — coverage is the whole point of the change. | **hand off** | Codex Sol; Totoro if the ladder grows | S0, S1 | Numbers exist for both wrap paths on the same seed/split. No "pigauto vs BACE" sentence before this. |
| **M1** | **MECHANICAL-VERIFY.** (a) `git log --oneline HEAD..origin/main -- R/fit_baseline_bace.R tests/ NEWS.md` — the 43-behind drift check on *our* paths only (no rebase). (b) `NOT_CRAN=true` focused testthat on the new file + `[T4]`. (c) `devtools::check()` clean on the wrap surface. (d) Confirm default path output is unchanged vs a pre-change run. | **hand off** | Codex Sol | S1 | All four green, output pasted — not asserted. |
| **M2** | New focused test file `tests/testthat/test-fit-baseline-bace-final-imp.R`: default-off back-compat; `final_imp=TRUE` returns same `mu`/`se` dims; `se` is finite and non-negative; both guarded `skip_if_not_installed("BACE")`. | Cursor Models | Grok 4.5 (Agent) | S1 | Tests pass with BACE installed; skip cleanly without it. |
| **N1** | NEWS entry. Must state: new opt-in path, default unchanged, `se` semantics differ, **and must not claim it fixes "imputed-as-observed"** (that bug was in BACE's `bace_final_imp`, which pigauto never called). | Other Models | Auto Cost / pinned Claude | S1, S3 | Rose reads it and finds no overclaim. |
| **X1** | **Melissa RECONCILE.** Plan-vs-actual: every row above vs git reality (`git log`, `git diff --stat`, test output, bench artefacts). Flag any row claimed done without a receipt. | Other Models | Auto Cost | M1, M2, N1 | Reconciliation table written; ≥2 NOT-DONE verdicts withhold the completion claim (D-43). |
| **X2** | After-task report + brain write **proposal** (not a write). Draft the `memory/DECISIONS.md` entry for the wrap-API choice; stage it for Shinichi's approval (D-37). | Other Models | Auto Cost | X1 | Report in `docs/dev-log/`; brain entry drafted, **not** committed to the vault. |

**Gates.** S0 fails → the whole lane stops and reports; do not fake a bench.
S3 shows the proper-MI path is *worse* on point accuracy but better on coverage → that is an
expected trade, report both, do **not** silently pick the winner.
Never `git add -A`; stage by explicit path; the uinit and DRAC files stay unstaged throughout.

---

## Classification of handover items

| # | Item | Class | Evidence |
|---|---|---|---|
| 1 | Lane preflight + `git status` + `gh pr view 155` + classify | **DONE** | Phase 0.25 receipt above. |
| 2 | Read `AGENTS.md`, both handovers, `R/fit_baseline_bace.R`, coordination board | **DONE** | All five read this session. |
| 3 | Ultra-plan G0, STOP for approval | **DONE (Phases 0–2)** → **OWED (Shinichi's answer)** | This file; question in Phase 0.4. |
| 4 | Install standalone BACE `@ce8bc87`; re-bench wrap | **OWED — after G0** | `requireNamespace("BACE")` → FALSE. Row S0/S3. |
| 5 | Optional docs-only OVR footnote | **OWED — approved-in-principle, unstarted** | Row S2; no `R/ovr_categorical.R` edit made. |
| 6 | PR #155 merge to parent | **OWED — acceptance NOT met** | `ubuntu-latest (release)` **pending**; CI has not concluded. Not merged. |
| 7 | P0 blockers + B1–B3 | **DONE (sibling lane)** | `58708ac`, `0ed6a73`, `5800b01` on this branch's history. |
| 8 | Standalone BACE `@ce8bc87` | **PROTECTED** | Read-only inspection only; `main...origin/main` clean. |
| 9 | In-tree `pigauto/BACE` `@de87d8c` | **PROTECTED** | `[behind 86]`, untouched, no vendor-sync. |
| 10 | Dirty uinit files (6) | **PROTECTED (do not stage)** | Still modified-unstaged in `git status -sb`. |
| 11 | DRAC `dev/gnn_attribution_*`, `script/*gnn_attribution*`, `script/returned_gnn_attr/` | **PROTECTED (do not stage)** | Still untracked. |
| 12 | Rebase onto `main` | **RETRACTED / forbidden** | 43 behind, locked; mitigation is paths-scoped diff M1. |
| 13 | EM restore (`max_iter > 0`) | **RETRACTED / forbidden** | Needs a new G0. |
| 14 | Ask Dan to re-run 10k | **RETRACTED** | Dan note: do not ask. Only the wrap question is worth sending. |

---

## Post-G0 handoff

Execution does **not** happen in this chat. On approval, paste the `/goal` prompt recorded in the
return message into a fresh Cursor session (or hand rows S0/S3/M1 to Codex directly).

> Related: `docs/dev-log/handover/2026-08-09-cursor-handover.md` ·
> `docs/dev-log/handover/2026-08-09-bace-dan-progress.md` · `docs/dev-log/coordination-board.md`
