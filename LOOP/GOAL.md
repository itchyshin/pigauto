# GOAL — pigauto GNN evidence sentinel pre-run (IMMUTABLE)

## Mission
Bounded Totoro pre-run (12 fits) measuring per-fit wall time and paired
blend-vs-baseline held-out loss for regime-map sentinels F1/F2/F3 at
n ∈ {100, 300, 1000}, 30% MCAR, candidate SHA `6fddd79`.

## Slices
- S4: commit-pinned driver `script/gnn_evidence_sentinel_prerun.R`
- S5: Totoro execution (12 parallel jobs, ≤100 workers)
- S6: timing table + revised full-campaign estimate

## Invariants
- Evidence-only lane: do NOT touch `codex/pigauto-0-11-trust-usability` / PR #174
- Branch from candidate `6fddd79`, not `handover/2026-08-09-cursor`
- `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`
- Preserve `r_cal=0` fallback; record floor interventions and failures
- No push, merge, release, or public claims from this lane

## Definition of done
S4–S6 complete with gate audit G1–G8; timing evidence returned to parent.
