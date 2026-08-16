# Campaign checkpoint — resume point for this lane (survives compaction)

Written 2026-08-16 ~07:20 MDT, immediately before a context compaction, so nothing
about the in-flight campaign lives only in a chat window. **Claude lane holds pigauto
until this is done; Shinichi is not starting the Cursor lane until then.**

## Where things stand

- `origin/main` = `175ebdc` · 9 PRs merged (#158–#166) · CI green ·
  `--as-cran` 0 errors / 0 warnings / 1 known note.
- **Regime-map campaign running on Totoro.** `~/pigauto_regime_map/`,
  `snakagaw@totoro.biology.ualberta.ca`. At 07:15: **4,791 / 5,400, 0 failures**, ~20 min out.
- A background poll is running locally; it exits when `pgrep -fc run_campaign.sh` hits 0.
- **Summariser already staged and syntax-checked on Totoro**:
  `~/pigauto_regime_map/summarise_regime_map.R`.
- `zi_count` bench retry still running locally in `.worktrees/bench-rerun` — low value,
  nothing depends on it, safe to abandon.

## The remaining sequence (this is the part that was only in chat)

1. **Confirm the campaign ended cleanly.**
   ```bash
   ssh snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_regime_map; \
     echo "results=$(ls results/regime_map | wc -l)/5400"; \
     echo "xargs=$(pgrep -fc run_campaign.sh || echo 0)"; \
     echo "fails=$(grep -c "^FAIL " campaign_regime_map.log || echo 0)"'
   ```
   Failures > 0 → **stop and report**; do not summarise a partial grid as if complete.
2. **Summarise it.**
   ```bash
   ssh snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_regime_map && Rscript summarise_regime_map.R'
   scp snakagaw@totoro.biology.ualberta.ca:~/pigauto_regime_map/regime_map_summary.md /tmp/
   ```
   The script emits the pre-registered ADEMP quantities: paired Δ (blend − baseline) with
   MCSE per cell, the `|Δ|/MCSE ≥ 3` detection rule, P(gate open), P(floor fired), the
   >5% safety tail, and an explicit failed-job count.
3. **Launch the coverage campaign** (24 cells, ~0.5 h; Shinichi approved both):
   ```bash
   ssh snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_regime_map && nohup bash run_campaign.sh coverage 100 > launch_coverage.log 2>&1 &'
   ```
   Design: `docs/dev-log/2026-08-15-coverage-campaign-design.md`.
4. **Summarise coverage** (needs a small variant of the summariser: `cover95` and `width`
   are already recorded per job by `regime_cell.R`).
5. **Write both up** into `docs/dev-log/2026-08-16-regime-map-results.md`, commit on
   `handover/2026-08-09-cursor` via the `GIT_INDEX_FILE` route (stale `index.lock`),
   push, and refresh `docs/dev-log/handover/2026-08-16-cursor-handover.md`.
6. **Tell Shinichi it is done** — that is his signal to start the Cursor lane.

## Interpreting the result (decided in advance, so it is not decided by the answer)

The design's honest expectations, stated before seeing the numbers:

- **F1 linear** is the *specificity control*. The joint baseline is near-optimal there by
  construction, so Δ ≈ 0 and a closed gate is the correct outcome. **A GNN "win" on F1 is a
  red flag, not a success** — it would suggest the comparison is leaking.
- **F2 nonlinear** is the GNN's claimed niche. The pre-run saw it beat the baseline 6/6 with
  gates wide open, which is the one encouraging signal so far.
- **F3 mixed-type** exercises the discrete calibration path.
- Either outcome is publishable: lift in mapped regimes → claim those regimes; OFF ≥ ON on
  the corrected pipeline → the headline becomes the safety-gated architecture. **Do not let
  the result choose the framing retroactively** — that decision is Shinichi's, and the
  roadmap (Tier 1.2) says so.

## Standing constraints (unchanged)

Never `git add -A`. The 15 carried-over dirty files stay unstaged — including `AGENTS.md`,
whose brain-write line is corrected but uncommitted. Do not relaunch the campaign from zero
(it is per-cell resumable). No CRAN submission. P1-8 and the brain-write-policy question are
Shinichi's, not this lane's. GLLVM.jl is deliberately untouched.
