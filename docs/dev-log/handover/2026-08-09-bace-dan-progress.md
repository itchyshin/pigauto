# Handover — Dan Noble BACE progress vs pigauto in-tree BACE

Date: 2026-08-09  
Author: Cursor Auto (synthesis only; no BACE / pigauto `R/` edits)  
Sources: `/Users/z3437171/Downloads/simulation-report.html` (published 8 Aug 2026); standalone clone `/Users/z3437171/Dropbox/Github Local/BACE` @ `ce8bc87`; in-tree `pigauto/BACE` @ `de87d8c`; `https://github.com/daniel1noble/BACE`; `R/fit_baseline_bace.R`; pigauto BACE benches + `useful/bace_results_snapshot/`. Slack used as context only, not evidence.

Dan Noble (`daniel1noble`) is a close collaborator (ANU; ~94 joint papers; BACE `aut` alongside Shinichi `cre` and Szymek). He owns the public remote `daniel1noble/BACE`. Do not invent a “BACE lead” or lab-role beyond that.

---

## Where local BACE lives + HEAD vs origin

| Clone | Path | Remote | HEAD | vs `origin/main` |
|---|---|---|---|---|
| **Working copy (use this)** | `/Users/z3437171/Dropbox/Github Local/BACE` | `https://github.com/daniel1noble/BACE.git` | `ce8bc87` 2026-08-09 *Updated project with new findings for imrpovements* (Dan) | **in sync** (`main...origin/main`) |
| **In-tree reference** | `/Users/z3437171/Dropbox/Github Local/pigauto/BACE` | same origin | `de87d8c` 2026-04-01 *fixing categorical formula…* (Szymek) | **behind 86 commits** |
| GitHub | `daniel1noble/BACE` public `main` | pushed 2026-08-08 22:24 UTC | same `ce8bc87` | — |

No `NEWS.md` on BACE. README on standalone documents mixed types, `bace()` / `bace_final_imp()`, OVR as optional (`ovr_categorical = TRUE`), default now multinomial. In-tree README is the older “should not be used until further testing” draft.

Nested `pigauto/BACE` is a separate git repo (not a submodule pin in pigauto’s index). AGENTS.md: do not modify BACE as pigauto work. Syncing the vendored tree is a separate, explicit ask.

Local extra branch on standalone: `feature/categorical-improvements` is 2 commits ahead of its origin counterpart (gaussian ±5 SD clip). Not on `main`.

---

## What Dan actually changed (file-level, matches Slack)

Slack paraphrase → what is in the repo + HTML (not the paraphrase).

Four-agent investigation write-up: `BACE/.agents/investigation-2026-08.md` (poisson / MI propriety / categorical+OVR / evaluation). Report cites it.

| Slack item | Measured in git / report | Files (Aug 8–9 commits on `main`) |
|---|---|---|
| Imputed used as observed | **Yes.** Dominant bug: `bace_final_imp` fit final models on chain point-fills as observed. `bace_imp` already reinstated response NAs; final phase did not. | `R/bace_final_imp.R`, `R/bace_imp.R`, `R/bace.R` (`5709c2b` *Fixes to functions to improve models*) |
| Priors by response type | **Yes.** Gaussian R back to Hadfield (V=1, ν=0.002); 2026-04 ν=2 inflated residual var ~3×. Poisson log-link-scaled R (V=0.1, ν=1) + PX α.V=1. | `R/model_functions.R` (`5709c2b`, `2ef0ee0`) |
| Poisson point estimate | **Yes (not in Slack; load-bearing).** Mean-of-lognormals → median rate + guard `3·max_obs+5`. Pre-fix catastrophe (rep 587: imputation 1,039,177 for truth 5). | `R/model_functions.R` |
| Type detection | **Yes (not in Slack).** Integers: drop marginal AIC race; non-negative integers → poisson. | `R/prep_functions.R` (`2ef0ee0`) |
| `sim_bace` tweaks | **Yes.** Distinct per-liability slopes + REs; returns `response_probs` for Bayes ceiling. | `R/simulate_simBACE.R` (`e78e67b`) |
| OVR default | **Yes (Slack omitted).** `ovr_categorical` default **FALSE** (true multinomial). Investigation: OVR quasi-separation in 1-obs-per-species; multinomial better calibrated and ~2.6× faster in this sim. | `R/bace.R`, `R/bace_imp.R`, `R/bace_final_imp.R` |
| ~10k sims + push | **Yes.** Study B design 5 types × 2 mechanisms × 1000 = 10,000; **9,997/10,000 ok** (3 failures). Commits `1d7db1d` (rerun + report) + `ce8bc87` (`.agents/` state). Wall ~5.3 h (project-state). | `dev/13_response_type_simulation.R`, `dev/simulation_results/*`, `.agents/simulation-report.html` |

Tests: `8314c99`, plus `tests/testthat/test-performance-fixes.R` (claimed 1337 pass / 0 fail in `.agents/project-state.md` — not re-run here).

PRs: none for this pass (direct pushes to `main`). Last merged PR is #12 (2026-06-14, GHA sim-reference sharding).

---

## Report numbers (from the HTML tables, not Slack)

**Study A** — gaussian slope, 40 reps each mechanism. Parameter recovery vs complete-case and oracle.

| Method | Mech | Slope bias | Coverage | Mean bias |
|---|---|---:|---:|---:|
| BACE | MAR | 0.065 | **0.950** | 0.011 |
| complete-case | MAR | 0.098 | 0.950 | −0.259 |
| oracle | MAR | 0.071 | 0.950 | 0.000 |
| BACE | MCAR | 0.070 | 0.975 | −0.003 |
| complete-case | MCAR | 0.045 | 0.975 | −0.002 |
| oracle | MCAR | 0.057 | 1.000 | 0.000 |

**Study B** — value recovery / calibration. `n_final = 15` in the report machinery. Cover (q) = empirical 2.5–97.5% of 15 draws (perfect-imputer ceiling ≈ 0.83–0.845). Cover (t) = t-interval (gaussian/poisson). Cover (prob) = probability-set (discrete).

| Type | Mech | n | Recovery (Pearson / class) | Spearman | Bal.acc | Bayes ceiling | Median bias | Cover (q) | Cover (t) | Cover (prob) | Runtime (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gaussian | MAR | 1000 | 0.830 | 0.795 | — | — | 0.005 | 0.844 | **0.949** | — | 10.2 |
| gaussian | MCAR | 1000 | 0.881 | 0.852 | — | — | −0.005 | 0.845 | **0.947** | — | 10.0 |
| poisson | MAR | 1000 | 0.687 | 0.713 | — | — | 0.035 | 0.867 | **0.916** | — | 12.1 |
| poisson | MCAR | 1000 | 0.744 | 0.706 | — | — | 0.075 | 0.902 | **0.932** | — | 11.9 |
| binary | MAR | 999 | 0.859 | — | 0.785 | — | −0.041 | 0.983 | — | **0.941** | 11.3 |
| binary | MCAR | 1000 | 0.855 | — | 0.826 | — | 0.000 | 0.981 | — | **0.938** | 11.3 |
| ordinal | MAR | 999 | 0.877 | — | 0.820 | — | −0.014 | 0.991 | — | **0.959** | 11.8 |
| ordinal | MCAR | 1000 | 0.871 | — | 0.831 | — | 0.000 | 0.990 | — | **0.957** | 11.8 |
| categorical | MAR | 1000 | 0.807 | — | **0.628** | **0.709** | 0.101 | 0.984 | — | **0.984** | 21.7 |
| categorical | MCAR | 999 | 0.783 | — | **0.635** | **0.722** | 0.091 | 0.981 | — | **0.984** | 22.1 |

Before vs after (HTML §4.3 + project-state; categorical DGP changed so accuracy is not like-for-like):

- Poisson catastrophe P(|std bias| > 1): pre MCAR 0.088 / MAR 0.146 → post 0.028 / 0.065.
- Pre-fix Study B coverage (q): gaussian ~0.76, poisson ~0.69 (HTML “What changed”).
- Categorical 2026-07 balanced accuracy 0.43–0.46 sat next to a degenerate DGP Bayes ceiling ≈ 0.50.

### Slack vs measured

| Slack | HTML / project-state |
|---|---|
| Coverage near nominal across types | Gaussian t **0.949 / 0.947** ≈ nominal. Poisson t **0.916 / 0.932** under. Binary/ordinal prob **0.938–0.959**. Categorical prob **0.984** (conservative). |
| Poisson/gaussian imputed-value corr 0.75–0.88 | Gaussian Pearson **0.830 / 0.881**. Poisson Pearson **0.687 / 0.744** (not in 0.75–0.88). Poisson log-scale 0.75 / 0.78 is in `.agents/project-state.md`, not the HTML Pearson column. |
| Categorical still a problem | Remaining: bal-acc **0.628 / 0.635** vs Bayes **0.709 / 0.722** (~8–9 pp below ceiling); sets slightly conservative. Report frames the *old* near-chance result as mostly DGP, not a model failure. Slack is stronger than the report. |

Regime: `sim_bace()` mixed-type phylogenetic MI, MCAR+MAR, production MCMC (nitt 4–5k, burnin 0.8–1k, thin 4–5, `n_final` 15). Not AVONET / PanTHERIA. Not pigauto.

---

## In-tree BACE vs this pass

In-tree `de87d8c` (1 Apr 2026) predates: OVR, `nitt_cat_mult`, `phylo_signal_summary`, `pool_mi` / `with_imputations` / accessors, the 2026-04 ν=2 prior change **and** its 2026-08 revert, the final-imp NA fix, poisson median guard, deterministic type rule, non-degenerate `sim_bace` multinomial DGP, and the 10k re-run.

`R/` diff (standalone vs in-tree): every shared modelling file differs; standalone-only: `accessors.R`, `phylo_signal_summary.R`, `pool_mi.R`, `with_imputations.R`.

---

## pigauto wrap + what pigauto already claims vs BACE

`R/fit_baseline_bace.R` calls **`BACE::bace_imp()` only** (`runs=5`, `nitt=4000`, `burnin=1000`, `thin=10`), averages chain datasets onto latent scale, returns `mu`/`se`. It does **not** call `bace()` or `bace_final_imp()`. The wrap therefore:

- already used the chain phase that reinstated NAs (the “imputed-as-observed” bug was in **final** imp);
- will pick up prior / type / poisson-point / OVR-default changes **only if** the installed BACE is current (wrap does not pass `ovr_categorical` → new default FALSE);
- is not the workflow Dan just validated (chain → final proper MI → pool).

Head-to-head **scripts** (`script/bench_avonet_bace.R`, `bench_bace_avonet_head_to_head.R`, `bench_pantheria_bace_head_to_head.R`) call **`BACE::bace()`** (includes `n_final`). Local markdown outputs say **“BACE skipped (not installed or failed)”** — no pigauto-vs-BACE numbers on those runs. GHA snapshot `useful/bace_results_snapshot/` is BACE Actions run **25329857467** (2026-05-13), i.e. **pre** the Aug MI/prior/poisson/DGP fixes.

pigauto Level-C OVR (`R/ovr_categorical.R`) still says it is “the OVR strategy BACE uses” (K independent threshold-joint / Rphylopars fits so phylopars stays full-rank). That is a **different estimator** from BACE’s MCMCglmm OVR binaries. BACE’s new default is multinomial, with OVR documented as miscalibrated in the 1-obs regime. The shared slogan is stale; pigauto OVR’s reason (phylopars singularity) is unchanged.

Claude memory `project_bace.md` (typo `"categorcial"`, early Gelman-prior TODOs) is stale relative to `ce8bc87`. Do not load it as current BACE state.

---

## Implication for pigauto (2–3 sentences)

The wrap’s `bace_imp` call is still a valid optional baseline entry point, but it is not the object Dan’s 10k study evaluated, and the in-tree vendored BACE is 86 commits behind the package that study used. Any pigauto-vs-BACE claim (GHA snapshot, skipped local benches, wrap-as-baseline) needs a **re-bench against `daniel1noble/BACE@ce8bc87`**, using `bace()`+`bace_final_imp` if the comparison is meant to match Dan’s MI propriety fix — not a silent `devtools::install("pigauto/BACE")`. Categorical overlap with pigauto’s OVR story is naming-only: BACE now defaults off OVR for MCMCglmm 1-obs; pigauto OVR exists so Rphylopars does not go singular — do not copy BACE’s multinomial default into pigauto without a separate measurement.

---

## Honest next contact with Dan (if any)

Do **not** ask him to re-explain the six bugs or re-run 10k for pigauto. He already wrote `.agents/investigation-2026-08.md` + the HTML.

If pinging:

1. Confirm we read `simulation-report.html` (Study B table + §4.3) and the investigation note.
2. One decision: should pigauto’s wrap stay on cheap `bace_imp` chain averages, or move to `bace()` + `bace_final_imp` now that final imp is proper MI? That is the wrap question; it is not “sync in-tree BACE this week.”
3. One alignment sentence for the ms / docs: pigauto Level-C OVR ≠ current BACE default. Ask if he wants a shared footnote; not a pigauto rewrite.
4. Optional later: refresh `useful/bace_results_snapshot/` from a new BACE GHA run when he is ready. Not urgent.

Do not recommend a “highest leverage” next move. Do not vendor-sync `pigauto/BACE` unless Shinichi explicitly asks.

---

## Resume

```bash
cd "/Users/z3437171/Dropbox/Github Local/BACE" && git log -1 --oneline && git status -sb
open "/Users/z3437171/Downloads/simulation-report.html"
```
