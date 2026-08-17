# After-task report — Arc 3 (solver convergence · per-type λ · UQ tail)

Date: 2026-08-17 · Lane: Claude (solo) · Scope authorised by Shinichi ("go ahead 1-6"),
with the #169 keep/revert decision left at my stated default (keep as yardstick).

## 1. What was asked

Six ranked items after Arc 2 closed: (1) make the in-house Σ solver converge, (2) decide
#169's fate, (3) multi-seed external comparison, (4) per-type λ dispatch, (5) real-data
Mondrian re-run, (6) the UQ tail. Items 3 and 5 were sequenced last by design.

## 2. What was delivered

- **Slice 1 — CLOSED BY MEASUREMENT, verdict: do not ship the port.** A pre-registered
  known-Σ recovery sim (12 cells × 50 reps × 6 arms, 0 failures) found the actual
  mechanism: **the in-house solver estimates Σ and then discards it.** E0 (single-pass) and
  E1 (Fisher-ML) produce *numerically identical* predictions in all 12 cells — paired
  difference exactly 0.0000 — despite Σ̂ differing by nearly 2× in Frobenius error, because
  `max_iter = 0` short-circuits to per-column BM. No Σ estimator can help while that holds.
  Fisher-ML is additionally worse at λ=1 (off-diagonal sign accuracy 0.993 → 0.912), the
  precise failure its iid-species assumption predicts. Docs:
  `2026-08-17-sigma-recovery-design.md` (pre-registered) + `...-results.md`.
- **Slice 4 — MERGED (PR #170).** Per-type λ dispatch: discrete traits keep their
  joint/OVR baselines when `lambda_mode != "fixed_1"`. Removes the measured 19-pp
  Trophic.Level regression (0.600 → 0.789, exactly the joint-path value) while retaining
  the continuous lift. `fixed_1` byte-identical, verified against origin/main.
- **Slice 6 — PR #171 at the gate.** Tree-MI conformal draws (P1-11) via an extracted
  `.conformal_draws()` helper with `multi_impute()` verified byte-identical; zi_count
  conditional-on-non-zero conformal intervals (explicitly *not* an E[X] interval); one
  definitive `$se`-by-type doc section (P1-6/7). Caught a latent bug in passing: the
  Mondrian locality lookup always indexed `latent_cols[1]`, silently wrong for zi_count.
- **Slice 2 — decided by default, stated:** `joint_solver = "rphylopars"` (#169) stays as
  the measured yardstick until the dependency-free fix reaches parity, then is removable.

## 3. What was verified, and how

Sim design pre-registered and committed *before* any estimator comparison ran. Every
merged PR: targeted tests → full suite → `--as-cran` (0E/0W/1 known NOTE) → 3-platform CI.
Suite grew 2,057 → 2,127 passes across the arc, 0 failures throughout. Slice 4's real-data
claim was verified on the same AVONET mask as the bench it cites, not a fresh one.

## 4. What is NOT covered

- Recovery sim: simulated matrix-normal BM (the model all arms assume), K=4, MCAR 30%,
  single-obs, n ≤ 300, λ ∈ {1.0, 0.2}, solver-level only — no GNN, no mixed types.
- Coverage measured there is the BM analytic SE, not conformal intervals.
- Slice 4's verification is single-seed (Mass swings materially rep-to-rep at n=300).
- Slices 3 (multi-seed external comparison) and 5 (real-data Mondrian) NOT started — both
  were sequenced to follow the solver verdict, which redirected the work.

## 5. Surprises

1. **The headline was the wrong question answered rightly.** The slice set out to improve a
   Σ estimator; the measurement showed Σ quality is irrelevant as wired. Had the sim not
   included both estimators at `max_iter = 0`, the identical-prediction signature — the
   entire finding — would have been invisible.
2. **My own pre-registered criterion was ill-posed.** Rule 1 compared Frobenius error
   against phylopars, whose `phylocov` is on a different scale; its error sits at ~1.1 in
   every cell regardless of design. A tree-height normalisation did not fix it. The rule is
   marked unevaluable in the results doc rather than quietly reinterpreted.
3. **Refinement trades accuracy for calibration.** Turning Σ on recovers ~85% of the gap at
   λ=0.2 but drives SE coverage 0.925 → 0.618 — the posterior variance is not updated when
   the mean moves. This is almost certainly the 2026-05-17 "EM diverged" note seen from
   another angle, and it converts a vague history into a testable defect.
4. A third silent-parameter-drop class was already fixed in Arc 2 (#167); the Mission
   Control board's standing recruit item predicted exactly this recurrence.

## 6. Decisions taken (and whose)

- Do not ship Fisher-ML (mine, forced by measurement; port left unmerged as opt-in internal).
- Keep single-pass Σ — better where it matters (high signal).
- Slice 4's hybrid: liability fit still *sees* continuous columns, but their output comes
  from the λ path via the existing fallthrough (agent's call, reviewed and endorsed).
- **NOT taken — needs your G0:** enabling refinement is the EM cell-restore that the Arc 2
  goal explicitly fenced ("EM restore needs its own G0"). Slice 1b is scoped and stopped.

## 7. Cost

Local only — no Totoro this arc (deliberate: the account was carrying two foreign lanes).
Three Sonnet agents (port, λ dispatch, UQ tail), ~72 min of simulation, four full suites,
three `--as-cran` runs.

## 8. State of the tree

- `origin/main` = `584ea81` (#170 merged). #171 open, gates nearly clear.
- `handover/2026-08-09-cursor`: +3 doc commits (sim design pre-registered, sim results,
  this report).
- Branches pushed: `arc/lambda-per-type` (merged), `arc/uq-tail` (PR #171),
  `arc/inhouse-sigma-convergence` (**deliberately unmerged** — the port does not ship, but
  the sim harness and the opt-in `sigma_method` live there for slice 1b).
- 15 carried-over files: untouched (verified). No `BACE/` commits. No CRAN action.

## 9. Open items, ranked

1. **Slice 1b — needs your G0** (crosses the EM fence): enable refinement, fix the
   posterior variance so coverage returns to ~0.93–0.95, add a stopping rule (more
   iterations helped at λ=0.2, hurt at λ=1.0). Gate = re-run the recovery sim, then the
   5-seed AVONET bench.
2. Merge #171 once CI/`--as-cran` clear (mechanical).
3. Slice 3 — multi-seed external comparison (after 1b, so it benchmarks the fixed solver).
4. Slice 5 — real-data Mondrian re-run (fishbase/pantheria).
5. Capability-surface row updates: the "joint baseline / optional dependencies" wording is
   stale twice over, and the conformal row now has two distinct remedies to name.
6. Paper framing — yours; evidence complete and strengthened again this arc.

## 10. Claims discipline check

No claim exceeds its regime. The one negative result (Fisher-ML) is reported as such with
its mechanism named. My own methodological error is recorded in the results doc rather than
buried. The identical-prediction finding is stated with the exact evidence (paired
difference 0.0000 in 12/12 cells) and the code line that explains it.

## 11. Handover

Snapshot, coordination board, and the Mission Control status file are refreshed in the same
close-out. Mission Control had been 15 days stale ("LONG HOLD — do not open work here")
through two active arcs; corrected and committed to the vault at `bbac784`.
