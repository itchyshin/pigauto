# Session Handoff (to Claude): Mondrian conformal, from shipped code to a defensible default

Meta: 2026-09-22 evening · from Claude (Fable 5.1, vault session, lease `claude:shinichi-brain:76566` on
`docs/dev-log/handover/,useful/` only) · platform-agnostic; either tool resumes from this file.

## Critical Context

1. **`conformal_method = "mondrian"` is already on `main`.** PR #168 merged (`0e7a34a` feature,
   `4a3e456` floor 10 → 19, `295c36c` back-compat test), then #171. Do not re-implement it. The work
   that remains is evidence and prose, not code.
2. **The default is still `"split"`, and that is correct today.** Mondrian has been verified only on
   simulated mechanisms (`docs/dev-log/2026-08-16-mechanism-coverage-results.md`). The real-data
   benches where undercoverage was first seen (fishbase 0.89 to 0.91 at n = 10,654; pantheria 0.87 to
   0.94) have never been re-run with it. Flipping the default before that is the thing this handover
   exists to prevent.
3. **This lane touched nothing in `R/`, `tests/`, or `script/`.** A foreign lane (codex/cursor) and
   four Claude lanes were live on pigauto at write time; the checkout is on
   `handover/2026-08-09-cursor` with 18 uncommitted files that are not this lane's (D-88).

## What Was Accomplished

- A step-by-step explainer of split conformal in pigauto for Shinichi and Szymek, with three SVG
  mechanism figures (data flow; the rank argument; clade failure vs Mondrian repair) and the split vs
  Mondrian comparison table: https://claude.ai/artifact/6m8UyPBoNHMV4LePeutMDV. Filed in the vault as
  `memory/pigauto-split-conformal-explainer.md`.
- A Methods paragraph for the paper, `useful/mondrian-methods-paragraph.md` (this branch), every
  number read from `R/fit_helpers.R`, `NEWS.md` and the 08-16 results doc. Placement: under the UQ
  Methods, after the split-conformal paragraph, before the MI draws. Five references to verify.
- Szymek's worry ("the cluster splits may not be done correctly") classified: **a measured finding,
  not a misunderstanding**. Calibration cells come from the observed complement; missing cells cluster
  in sparse clades; the two are not exchangeable. MAR_phylo coverage 0.923 (n = 300) / 0.927
  (n = 1000) vs MCAR 0.961 / 0.957; does not wash out with n.

## Current Working State

- Working: Mondrian on `main`, 7 tests in `tests/testthat/test-mondrian-conformal.R` (locality
  hand-check, near/far scores differ, 19-floor fallback, clear errors on missing inputs, isolated
  clade widens, split reproduces the formula, split carries no mondrian state).
- In progress: nothing in code.
- Not done: real-data re-run with Mondrian; the default flip decision; the paper paragraph is a draft
  in `useful/`, not merged into `useful/paper_section_draft.md`, which currently has no UQ section at
  all (headings 1 to 7 cover architecture, covariates, evidence, discussion, limitations).

## Key Decisions & Rationale

- Default stays `"split"` until real-data confirmation (this handover, and the roxygen on `main`:
  "These choices support nominal held-out diagnostics, not package-certified coverage").
- Stratum floor is 19 residuals, the smallest n with n/(n+1) ≥ 0.95, so ~38 val cells per trait; below
  it Mondrian is a no-op by design and records `fallback`. The remedy at small n is more held-out data,
  not more strata (08-16 results, B2 verification).
- Mondrian is single-obs only (locality is per species) and requires `gnn = TRUE`
  (`R/fit_pigauto.R:429–455`).
- D-139: the real-data re-run is a campaign. Estimate before running (below).

## Landing State

`tools/handoff_gate.sh /Users/z3437171/Dropbox/Github\ Local/pigauto` → **GATE FAIL**: 18
uncommitted files on `handover/2026-08-09-cursor`, 89 unpushed commits across 13 other branches. None
of it is this lane's; declared here so it is not invisible.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `pigauto` `handover/2026-09-22-claude-mondrian` (this file + `useful/mondrian-methods-paragraph.md`) | y | y | none, docs only | LANDED |
| `pigauto` `main` `ab02e31` (Mondrian feature, PR #168; UQ tail, PR #171) | y | y | merged | LANDED (prior lanes) |
| `pigauto` checkout `handover/2026-08-09-cursor`: `M .gitignore AGENTS.md CLAUDE.md R/utils_torch.R README.md _pkgdown.yml`; `?? .ignore dev/gnn_attribution_low_lambda_smoke.{R,log,md} script/bootstrap_gnn_attribution_fir.sh script/collect_gnn_attribution_low_lambda.R` | n | n | none | CARRIED-OVER, **another lane's** (GNN attribution work); not read, not staged (D-88) |
| 13 branches with unpushed commits (`feat/gnn-off` 12, `spec/vulcan-gpu-avonet9993` 8, `handover/2026-09-19-claude` 2, ten others with 1 each) | y | n | various | CARRIED-OVER, **other lanes'**; listed by the gate, owners unknown from here; `lane_rescue.py` backs them up nightly |
| Vault `memory/pigauto-split-conformal-explainer.md` `61b7f585` | y | local-only by design (D-37) | n/a | LANDED |

FINDING-OF-RECORD: the clade-split exchangeability failure and its Mondrian repair (measured 08-16;
explained 09-22) vault-note: [[pigauto-split-conformal-explainer]].

## Next Immediate Steps

1. **Real-data re-run with Mondrian** (the gate for the default flip). Same seeds and settings as the
   committed benches, `conformal_method = "mondrian"` added:
   `script/bench_fishbase.R`, `script/bench_pantheria_bace_head_to_head.R`,
   `script/bench_avonet_full_local.R`. Report per-trait coverage split vs mondrian, realised n_val,
   how many traits hit the fallback, median width change. Compute: Totoro, ≤150 cores (D-143).
   D-139 estimate: fishbase at n = 10,654 ran the 08-16 mechanism cells at 10 to 20 min each for
   n = 1000; a single full-data fit with Mondrian is one fit, so expect **1 to 3 h wall per dataset**,
   under the 30-minute line only if run as one job each; state the estimate in the run log before
   launching. Pre-run: AVONET300 with mondrian (minutes) to confirm the option threads through
   `impute()`.
2. **Decide the default** on that evidence: flip to `"mondrian"` only if every real-data dataset is
   within 3 MCSE of nominal without MCAR-style width inflation above ~10%, and only for the regime
   where it activates (single-obs, gnn on, ≥38 val cells per trait). Otherwise keep `"split"` and
   make the fallback warning louder.
3. **Merge the Methods paragraph** into `useful/paper_section_draft.md` as a new UQ section, verify
   the five references against DOIs (Garfield), and add the 08-16 coverage table as a supplementary
   table.
4. Optional, cheap: an `impute()` message when Mondrian falls back for a trait, naming the realised
   stratum sizes, so users at n ≈ 300 learn why their intervals did not widen.

## Blockers / Open Questions

- Who owns the 18 uncommitted files on the checkout and the 13 unpushed branches. Not this lane's to
  resolve; overlap is Shinichi's call (D-87).
- Whether MNAR (value-dependent) missingness, which Mondrian repaired only marginally (0.932 → 0.940),
  needs a weighted-conformal arm before the paper claims coverage under "realistic missingness".

## Gotchas & Failed Approaches

- The first B2 run used a 10-residual stratum floor and pinned MCAR coverage at exactly 13/14 = 0.929:
  the per-stratum n/(n+1) ceiling. That is why the floor is 19. Do not lower it.
- Mondrian with `gnn = FALSE` errors by design; the GNN-off campaign of 09-19 is a different interval.
- `make_missing_splits()` under a mechanism weights one sample then splits it into val and test, so val
  and test are exchangeable with each other even under MNAR; that design cannot reproduce the real
  seam. The 08-16 campaign applied mechanisms as genuine NAs for that reason. Reuse its runner
  (`mech_cell.R` on Totoro, `~/pigauto_regime_map/`) rather than the April bench.
- Do not quote the fishbase/pantheria numbers as "intermittent" or as ceiling effects; the 08-16 work
  showed the mechanism explains them at n where the ceiling cannot.

## How to Resume

```sh
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
~/shinichi-brain/tools/lane_preflight.sh .          # name the ONE lane you take
git fetch origin && git switch -c arc/mondrian-realdata origin/main
sed -n '1,60p' docs/dev-log/handover/2026-09-22-claude-handover-mondrian-conformal.md
sed -n '1,120p' docs/dev-log/2026-08-16-mechanism-coverage-results.md
grep -n "mondrian" R/fit_helpers.R | head -40
```
Then step 1 above. Read order: this file → 08-16 results → `NEWS.md` "mondrian" entry →
`tests/testthat/test-mondrian-conformal.R` → the three bench scripts.

**One-command resume (paste into a fresh Claude session opened in the pigauto repository):**

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-22-claude-handover-mondrian-conformal.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

Rehydration for Claude: run `tools/lane_preflight.sh` from the vault, name the one lane you take, then
classify every item above as OWED / DONE / RETRACTED / PROTECTED against `origin/main` before doing
anything. Expected as of 2026-09-22: Mondrian code DONE (PR #168); real-data re-run OWED; paper UQ
section OWED; the 18 dirty files and 13 unpushed branches PROTECTED (other lanes). No snapshot pointer
exists in this repo's `AGENTS.md` and there is no coordination board; with a foreign lane and four
Claude lanes live, none was added, so that a single pointer cannot orphan a sibling (handover-skill
Step 4). Environment: R with torch installed as this repo's `AGENTS.md` describes; safe verification
is `devtools::test(filter = "mondrian-conformal")`; campaigns go to Totoro under the 150-core cap, never
run locally at scale (D-200). Do not stage anything under `dev/` or `script/` you did not create.
