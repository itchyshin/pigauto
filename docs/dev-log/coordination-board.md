# Coordination Board — pigauto

## Active Lane Split

| Lane | Owner | Branch | Handover / truth | Status |
|---|---|---|---|---|
| **P0 review blockers + B1–B3** | Claude (2026-08-15 arc) | `arc/p0-onto-main` `f5e8416` (from `fix/ci-install-libtorch` `21d2ea6`) | `docs/dev-log/handover/2026-08-15-p0-arc-handover.md` · after-task `docs/dev-log/after-task/2026-08-15-p0-arc.md` · PR [#155](https://github.com/itchyshin/pigauto/pull/155) → `fix/ci-install-libtorch` | **MERGED ONTO A BRANCH, NOT LANDED.** 7/7 conflicts resolved, suite run, branch pushed, **no PR**. **BLOCKER:** `test-safety-floor.R:712` strict val-floor passes on `main`, fails on `main`+P0 (controlled). `main` unchanged at `416561b`. |
| **BACE wrap / re-bench (pigauto remit)** | **LANDED on main** | `origin/main` `416561b` (PR [#156](https://github.com/itchyshin/pigauto/pull/156) **MERGED**) · Version `0.10.0.9000` | Fresh Claude: `docs/dev-log/handover/2026-08-14-claude-handover.md` · prior Cursor: `docs/dev-log/handover/2026-08-13-cursor-handover.md` · plan `docs/dev-log/handover/2026-08-13-wrap-merge-g0-plan.md` · after-task `docs/dev-log/after-task/2026-08-13-bace-wrap-merge.md` | **DONE.** GitHub-dev wrap on `main`. Post-merge CI green. Not a CRAN tarball. BACE still 404. Do not re-merge. Do not CRAN-submit from this `main`. |

Rehydrate must read **both** rows. A single AGENTS.md snapshot pointer would orphan a sibling — there is no Live Phase Snapshot in AGENTS.md; **this table is the split**.

## Current Rule

- Wrap **LANDED** on `origin/main` `416561b` (PR #156 MERGED, 2026-08-13). Version `0.10.0.9000`. Suggests `BACE`. **Do not re-merge #156. Do not CRAN-submit from this `main`** until BACE is on CRAN or Suggests BACE is dropped again.
- P0 sibling #155 remains on `fix/ci-install-libtorch` only. **Do not merge that parent to `main` unless Shinichi locks it as a new G0 that day.**
- Next Claude reads `docs/dev-log/handover/2026-08-14-claude-handover.md`. Rehydrate both Active Lane Split rows. If wrap still landed and no new G0: **STOP**.
- Do not rebase `handover/2026-08-09-cursor` onto main. Leave its dirty uinit / gnn tree unstaged.
- Do not modify standalone BACE or in-tree `pigauto/BACE`.
- No EM restore (`max_iter>0`) without a new G0. No DRM.jl from this board.

## Status

- 2026-08-09: P0 landed on origin; #155 merge-from-parent done; R-CMD-check running. BACE wrap lane opened via Cursor handover. Dirty uinit / DRAC scripts / README banners stay unstaged.
- 2026-08-09 (later): Wrap lane ultra-plan Phases 0–2 complete, **stopped at G0**. Sweep found no
  prior wrap attempt (no branch / worktree / stash / brain decision) and that **BACE is not
  installed on this machine** (`requireNamespace("BACE")` → FALSE) — so the T4 smoke and all three
  local benches skip for environment reasons, not code reasons. PR #155 CI still **pending**;
  acceptance not met, **not merged**. No `R/`, BACE, or commit changes made.
- 2026-08-09 (G0): Shinichi locked **Option B-minus** (`final_imp = FALSE` default + `n_final = 15L`).
  Slack later parked (`/goal`: "No Slack to Dan") — do not draft or send.
- 2026-08-09 (G0 LOCKED): Shinichi approved Option B-minus — *"it should be final yes — OK opt-in
  cool /ultra-plan it"*, then confirmed in a structured AskQuestion: B-minus, **`n_final = 15`
  explicit**. Phase 3 ran on `handover/2026-08-09-cursor`. Do not start a second Phase 3 / ultra-plan.
  PR #155 untouched.
- 2026-08-09 (KEEP wrap): Shinichi chose **keep the wrap** and asked for a landing-path design
  onto current main (`b615579` deleted `R/fit_baseline_bace.R` for the v0.10.0 CRAN surface).
  **No merge this turn.** Slack to Dan stays parked (`/goal`: "No Slack to Dan"). No public
  pigauto-vs-BACE sentence. No vendor-sync. M1c wrap Rd `\usage` **OK**; focused tests
  18 pass / 1 skip; full check still 2 WARN / 3 NOTE (pre-existing). Receipt
  `/tmp/pigauto_m1c3/pigauto.Rcheck/00check.log`.
- 2026-08-09 (landing timing): Shinichi locked **A after #155** — new branch from
  `origin/main`, restore wrap, **GitHub-dev-only** until BACE is on CRAN / v0.10.0 ships.
  Do not rebase `handover/2026-08-09-cursor`. Wrap lane waits on #155 clearing `R/`.
- 2026-08-09 (landing executed): #155 is **MERGED** to `fix/ci-install-libtorch`. Wrap
  restored on isolated worktree `/tmp/pigauto-bace-wrap-restore`, branch
  `feat/bace-wrap-restore` @ `b180555` (from `origin/main` `bf46991`). Pushed
  `origin/feat/bace-wrap-restore`. Focused tests `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]`.
  **No PR into `main`.** Handover checkout left dirty/untouched. Neither BACE tree touched.
- 2026-08-09 (wrap closeout): Focused tests re-run on
  `/tmp/pigauto-bace-wrap-restore` → `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]`.
  Draft DO-NOT-MERGE PR [#156](https://github.com/itchyshin/pigauto/pull/156)
  (`feat/bace-wrap-restore` → `main`). After-task
  `docs/dev-log/after-task/2026-08-09-bace-wrap-closeout.md`. Melissa LIGHT
  `docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md`.
  **STOP until CRAN.** Do not merge. Do not start DRM.jl / DESCRIPTION
  claim-gate / EM / Slack / vendor-sync / public pigauto-vs-BACE.
- 2026-08-11 (wrap handover): Wrap row marked **CLOSED / wait CRAN / PR 156**.
  Durable Cursor handover
  `docs/dev-log/handover/2026-08-11-cursor-handover.md`. `origin/main` still
  Version 0.10.0 (`bf46991`). Do not re-restore. Do not merge #156. DRM.jl is
  a sibling repo with its own handover — this board stays pigauto-only.
- 2026-08-13 (CRAN live, plan first): CRAN `pigauto` 0.10.0
  Date/Publication 2026-07-30 17:10:22 UTC. CRAN `BACE` still 404.
  Shinichi picked **plan first**, not merge. Plan
  `docs/dev-log/handover/2026-08-13-wrap-merge-g0-plan.md`. #156 remains
  draft. Do not undraft or merge until G0 lock.
- 2026-08-13 (G0 locked + merged): Shinichi locked merge G0. NEWS
  `a54e6a4` on wrap branch. PR [#156](https://github.com/itchyshin/pigauto/pull/156)
  MERGED (`416561b` on `origin/main`). Version `0.10.0.9000`. Do not
  CRAN-submit. P0 parent still not on `main`.
- 2026-08-14 (Cursor → Claude handover): Wrap still LANDED `416561b`.
  Post-merge CI green (R-CMD-check 31755795655, pkgdown 31755795738).
  CRAN pigauto 0.10.0 live; CRAN BACE 404. Memo verdict
  **agree-with-corrections**. Recommended next G0 *if named*: land P0
  on main (parent 46 behind), then claim-gate. **Not locked — STOP.**
  Fresh handover `docs/dev-log/handover/2026-08-14-claude-handover.md`.
  Dirty uinit / GNN / 44 foreign unpushed commits stay CARRIED-OVER.
- 2026-08-14 (Claude rehydration): Claude rehydrated from
  `docs/dev-log/handover/2026-08-14-claude-handover.md`. All handover claims
  reconcile against live git/gh/CRAN (wrap `416561b` landed, CI green,
  #155 parent-only 46/9, CRAN pigauto 200 / BACE 404, #135 OPEN) — zero
  discrepancies. Brain sweep found **no decision locking "land P0 on main"**.
  **No G0 locked — STOP stands.** P0 remains off `main`. Dirty uinit / GNN and
  the 44 foreign commits stay CARRIED-OVER; no R code ran. Receipt
  `docs/dev-log/after-task/2026-08-14-claude-rehydration.md`.
- 2026-08-15 (P0 arc — merged, verified, **NOT landed**): Shinichi locked
  "land P0 + fix the false doc claim" **conditional on a clean check**. Merge
  of `fix/ci-install-libtorch` onto `416561b` done on `arc/p0-onto-main`
  `f5e8416` (pushed, **no PR**): 7/7 conflicts resolved exactly as the
  `git merge-tree` probe predicted, `document()` clean. **Check not clean →
  did not land.** `test-safety-floor.R:712` strict val-floor **passes on
  pristine `main`, fails on `main`+P0** (0.3434 / 0.3442 vs 0.3410),
  controlled comparison, test provably ran not skipped. Hypothesis: P0 fix #4
  removes a held-out-context leak, so the GNN is honestly worse on val and the
  5% threshold was calibrated on the leaky pipeline — **REFUTED by experiment
  2026-08-15**: restoring the leak shows the gate CLOSES pre-fix (`r_cal_gnn`
  0.0, blend == bm) but OPENS to 0.1 under P0 and degrades below pure BM. That
  is a **safety-invariant failure, not a stale threshold** — `r_cal = 0` stopped
  being a protective fallback. Next question is why calibration selects a
  harmful gate; do NOT loosen the threshold. `--as-cran` deliberately NOT run
  (no value on a branch with an unexplained failure). Claim-gate: **1 BLOCKING**
  — `gnn-architecture.Rmd` §5 says `pred$se` is "Exact under BM" while `main`'s
  default path broadcasts one `sd(observed)` to every tip (verified in code);
  P0 is the fix. `main` UNCHANGED at `416561b`. Shinichi stepping away from
  pigauto; **do not resume this lane unasked.** Handover
  `docs/dev-log/handover/2026-08-15-p0-arc-handover.md`.
