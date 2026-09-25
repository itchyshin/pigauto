# Handover: posterior multiple imputation (Claude, 2026-09-25)

## Start here

- Lane: `claude:pigauto-mi-posterior`.
- Worktree: `/Users/z3437171/Dropbox/Github Local/pigauto-mi-posterior`.
- Branch: `arc/mi-posterior`, pushed.
- Draft PR: itchyshin/pigauto#189 (base `main`; #187 merged). It stays a draft and is not merged.
- Read, in order:
  1. `docs/dev-log/mi-posterior/results.md`
  2. `docs/dev-log/after-task/2026-09-24-mi-posterior.md`
  3. `docs/dev-log/mi-posterior/design.md` (sections 5d and 5e record every decision)
  4. `docs/dev-log/plan-actual/2026-09-24-mi-posterior.md`

## What exists

`multi_impute(draws_method = "posterior")` for continuous traits. It gives proper Bayesian multiple
imputation from Sigma_P %x% R + Sigma_E %x% I, with automatic chain extension, and ships with tests,
documentation and NEWS.

Simulation: 40 regimes x 200 reps.
- In-model G7 (per-cell coverage) passes.
- In-model G6 passes on the pooled relative SE-ratio rule Shinichi chose on 2026-09-25 (mean
  1.059 in [0.95, 1.10]; `design.md` 5f). Three rows (twins 35, 36, 38) sit above the per-row band;
  the M2 review showed this is real heterogeneity, not noise. They are reported, not gated; the MI
  SE ratio there is 0.98 to 1.06 in absolute terms.

Real data: all 10 cells done; G8 met (`REALDATA_COMPLETE`).

## Decisions owed by Shinichi

1. ~~The G6 SE-ratio rule~~: decided 2026-09-25 (pooled mean, option 3) and applied.
2. **Default draws method.** It is still `"conformal"`; changing it was never in scope.
3. **Merge of PR #189**, after 1 and 2.
4. **Optional follow-ups:**
   - a `log_transform = FALSE` PanTHERIA sensitivity run;
   - the Sigma_E prior's pull on the residual correlation at lambda near 1;
   - tip-variance heterogeneity on non-ultrametric trees (pigauto-wide).
5. **Brain proposals** staged in `docs/dev-log/mi-posterior/brain-proposals.md` (not written to the
   vault).

## Landing State

| Item | State | Why | Resume |
|---|---|---|---|
| `arc/mi-posterior` | pushed; draft PR #189 open | CARRIED-OVER: awaiting Shinichi's decisions (default, merge) | `cd <worktree>; git pull; gh pr view 189` |
| `.unlazy/mi-posterior/GATES.md` | local only (`.unlazy/` in `info/exclude`) | the acceptance ledger; all 12 gates met (2026-09-25) | `node ~/.claude/skills/unlazy/scripts/gate-check.mjs --status .unlazy/mi-posterior/GATES.md` from the worktree |
| Totoro run folders `/home/snakagaw/pigauto_mi_posterior/{69670d44f9,9597e18b79,7a0f47030f,diag_*}` | on Totoro (home, not scratch) | raw campaign outputs; summaries and evidence are committed | ssh via `~/.ssh/cm-snakagaw@totoro.biology.ualberta.ca:22` |
| Unpushed commits on other branches (reported by `handoff_gate.sh`: `spec/vulcan-gpu-avonet9993`, several `codex/*`, `experiment/*`, `handover/*` and others) | not this lane's | other lanes' state; not touched (D-88) | owners of those lanes |

The durable findings are listed below. Each names the vault note it should become. None of these
notes is written yet: the vault is approval-gated for this repo's agents, so every one is a PROPOSAL
awaiting Shinichi (full text in `docs/dev-log/mi-posterior/brain-proposals.md`). The canonical source
for each is `results.md` or `diagnosis.md` on `arc/mi-posterior`.

FINDING-OF-RECORD: proper posterior MI fixes the downstream slope attenuation of conformal MI (paired bias within 0.02 in-model; per-cell coverage about 95%)  vault-note: [[pigauto posterior multiple imputation gives honest per-cell intervals and near-unbiased downstream slopes]] (PROPOSED, not yet written)
FINDING-OF-RECORD: a gate script must exit non-zero on failure and print its token only on the pass line  vault-note: [[A gate script must fail closed]] (PROPOSED, not yet written)
FINDING-OF-RECORD: BLAS thread caps must be exported in the launching shell  vault-note: [[Export BLAS thread caps in the launching shell, never only inside R]] (PROPOSED, not yet written)
FINDING-OF-RECORD: an MCSE for a ratio of SDs from the same replicates must include their correlation  vault-note: [[A Monte Carlo SE for a ratio of paired SDs must account for their correlation]] (PROPOSED, not yet written)
FINDING-OF-RECORD: in-model twins separate model-family mismatch from estimator error  vault-note: [[In-model twins separate model-family mismatch from estimator error]] (PROPOSED, not yet written)

## Compute notes

- Totoro was heavily contended all day (load 300 to 500 from other users).
- Other lanes of this account launched uncapped jobs: drmTMB at 16:41 and `exact_prerun_cell.R` at
  about 17:52, taking the account above 150 cores. Shinichi chose to keep this lane's campaign
  running.
- The FishBase real-data cell (10,484 tips, 5 traits) ran for more than 16 h at 69670d4.
- All of this lane's processes are listed under Final state.

## Final state (2026-09-25, about 04:30)

Ledger (local `.unlazy/mi-posterior/GATES.md`, re-verified with `gate-check --reverify` at the final
`R/` and `tests/` trees):

| Gate | State |
|---|---|
| G1 to G5c | met |
| G6 | not met: 3 relative SE-ratio rows (decision owed) |
| G7 | met |
| G8 | not met yet: the FishBase cell was still running (see below) |
| M1 | met |
| M2 | findings fixed; the D-43 panel withheld nothing (`docs/dev-log/mi-posterior/d43-panel.md`); the formal PROCEED waits for the final real-data section |

Evidence: `docs/dev-log/mi-posterior/evidence/ledger/`.

The FishBase cell, 69670d4 code, is the last real-data cell. It had run for more than 19.5 h. It uses
1 core on Totoro (pid in `ps -u snakagaw`, command `01_run.R --args fishbase structured 20260818`). A
tiny timing run at 02:16 measured about 1.1 s per sweep on a quiet machine, so a full cell needs
about 7.3 h unloaded; the day's heavy contention slowed it several-fold. To finish G8 once
`real_run.log` shows 10 `CELL_OK`:

1. `rsync -a -e "ssh -o ControlPath=$HOME/.ssh/cm-snakagaw@totoro.biology.ualberta.ca:22" --exclude logs snakagaw@totoro.biology.ualberta.ca:/home/snakagaw/pigauto_mi_posterior/69670d44f9/campaign/real/ <worktree>/script/mi_realdata/returned/`
2. From the worktree: `Rscript script/mi_realdata/02_summarise.R`, then `Rscript script/mi_realdata/03_acceptance.R`
   (G8 prints `REALDATA_COMPLETE` only if all 10 cells are ok, fail-closed).
3. Replace `docs/dev-log/mi-posterior/real_preview/` with the final tables. Update the real-data section
   of `results.md`, recomputing every number. Then record M2's PROCEED in `GATES.md` and run
   `gate-check --reverify`.

If the cell is killed or fails, G8 cannot pass on one code SHA without rerunning FishBase. Rerunning
it at the new code would also require the other 9 cells at that code (G8 requires one SHA). That
decision is Shinichi's.

Processes of this lane still running on Totoro: only the FishBase cell. Everything else finished or
was killed by this lane. Leases: released at the end of the run (`lane_lease.sh --list pigauto`).

## Update (2026-09-25, after the G6 decision)

- FishBase finished and converged; G8 is met with all 10 cells at 69670d4 (commit 2a89591). No
  process of this lane runs on Totoro now.
- Shinichi chose G6 option 3. `04_acceptance.R` gates the mean relative SE ratio over the 24 gated
  phylolm rows in [0.95, 1.10] (`MI_SE_RULE=pooled_relative`, the new default) and reports each
  row. Result: 1.0591, `SIM_ACCEPT_PASS`. Recorded in `design.md` 5f.
- Ledger: all 12 gates met. G6 was re-run after the rule change. A full `--reverify` was started
  and stopped at Shinichi's request; the last complete one is `evidence/ledger/gate_reverify_2026-09-25.log`
  (the `R/` and `tests/` trees are unchanged since). To redo it:
  `node ~/.claude/skills/unlazy/scripts/gate-check.mjs --reverify --root <worktree> --cwd <worktree> --timeout 3600 .unlazy/mi-posterior/GATES.md`.
- The "Final state" section above is the record as of 04:30 and is kept as written.
- Still owed by Shinichi: the default draws method, the merge of #189, and the brain proposals.
