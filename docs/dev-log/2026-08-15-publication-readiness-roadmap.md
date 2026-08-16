# pigauto — publication readiness: gaps and work plan

Date: 2026-08-15, **updated 2026-08-16 (overnight: Tier 0 CLOSED)** · Author: Claude
(Fable 5 → Opus 5), from the P0-arc evidence base
Audience: Shinichi, returning to pigauto after a break, deciding what to do before a methods/software paper.

Every item below is motivated by a **measurement or filed finding**, not a hunch. Sources:
`docs/dev-log/arc/2026-08-15-p0-onto-main-arc-notes.md` (F1–F6, three experiments) ·
`docs/dev-log/arc/2026-08-15-p0-claim-gate-findings.md` · GitHub
[#157](https://github.com/itchyshin/pigauto/issues/157) ·
[#135](https://github.com/itchyshin/pigauto/issues/135) · Rose P1 queue
(`docs/dev-log/handover/2026-08-08-rose-close.md`) · `NEWS.md` · `cran-comments.md`.
Regime discipline applies throughout: where a number is cited, its regime is named.

---

## 0. Honest baseline — where the package stands today

**Solid, and worth saying so:**

- CRAN release 0.10.0 (published 2026-07-30). 558+ tests. Eight trait types in one latent
  matrix. Bundled AVONET300 + trees300 + ctmax_sim. Vignette suite with an unusually
  well-hedged claims culture — the adversarial claim-gate audit found `_pkgdown.yml`,
  `tree-uncertainty.Rmd`, and `common-pitfalls.Rmd` completely clean.
- **One validated island exists and is the template for everything else:** the analysis-aware
  MI backend was gated by a 6,000-task campaign — all 24 fixed-effect cells, 93.9–96.3%
  coverage, pooled-SE/empirical-SD 0.942–1.030 (scope: ONE incomplete continuous covariate,
  MAR, Gaussian `lm` / binomial-logit `glm` / one-random-intercept `lmer`). That is what a
  publishable claim looks like: narrow, gated, multi-seed, recorded.
- The unique selling points remain intact and no other package combines them:
  (a) unified mixed-type imputation; (b) conformal intervals tied to validation residuals;
  (c) multi-tree Rubin pooling; (d) multi-obs `obs_refine`; (e) `suggest_next_observation()`.

**~~Not solid~~ — ALL FOUR CLOSED overnight 2026-08-15/16. Kept for the record:**

- ~~Four correctness fixes off `main`~~ → **LANDED** via PR #158 (`main` `3677a85`).
- ~~`main` ships a documented guarantee its code does not satisfy~~ → **FIXED** in the same
  PR: the `gnn-architecture.Rmd` §5 SE row is split by path, and the per-tip Henderson SE
  now makes the joint path's claim true.
- ~~#157 gate/predict surface mismatch~~ → **FIXED and CLOSED**: a two-layer margin-based BM
  floor (`calibrate_gates()` + a post-refine floor in `fit_pigauto()`), plus the val-floor
  test moved to the calibration-matched surface.
- ~~Every GNN benchmark predates the leak fix~~ → **RE-RUN**; the pre-fix numbers reproduce
  within 1.1 MCSE in every completed type (§1.1). The provenance cloud has lifted.

**What is genuinely not solid, as of 2026-08-16:**

- The GNN's *value* is still unmeasured across regimes (#135) — the re-run says the leak did
  not move the numbers, not that the GNN earns its place. That is Tier 1.2, designed and
  priced, awaiting a go/no-go.
- Conformal coverage vs n is unmeasured (Tier 1.3, same status).
- **CORRECTION (2026-08-16, same night):** an earlier draft of this line said the joint
  baseline "ships as a single-pass Henderson init; cross-trait Σ is never iterated". That
  was **wrong**, and it framed a deliberate evidence-based choice as a deficiency. Reading
  `R/joint_mvn_solver.R:291–305`: with `max_iter = 0` the solver still returns "the
  closed-form `L̂ᵀR⁻¹L̂/n` Σ estimate — consistent under matrix-normal BM". What is disabled
  is the **cross-trait EM cell-refinement**, switched off 2026-05-17 because it *diverged*
  on strong-signal data (Σ blow-up via `build_conditional_prior()` ignoring R-mediated
  cross-row correlation) and was empirically much worse (0.93 → 0.53 argmax accuracy,
  K=5 BM, n=100, miss=0.30). Σ estimation is not missing; a refinement loop that made
  things worse is switched off, with the reason recorded in the source.
  **The real Tier-2.1 question is therefore narrower:** is the closed-form Σ estimate
  adequate, and does any doc claim more than it delivers? (See P1-5 below.)
- No known-Σ recovery simulation exists for the multivariate machinery (Tier 2.2).
- Five Rose P1 items remain (P1-5, P1-8, P1-9, P1-11, P1-12); four closed overnight
  (P1-6, P1-7, P1-10, P1-13).

---

## ~~TIER 0~~ — CLOSED 2026-08-16. Kept for provenance; see `2026-08-15-p0-land.md`.

### 0.1 Decide #157, then land P0 — **DONE** (surface B locked; PR #158 merged; #157 closed)

The paper's central safety claim — the calibrated gate guarantees the blend is never worse
than the phylogenetic baseline — is currently **enforced on the wrong surface**. Calibration
evaluates val cells baseline-filled and unpinned (`R/fit_pigauto.R:819`, `:429`); `predict()`
pins them to their own truth (`R/predict_pigauto.R:339`). Measured gap on `main`:
max|A−B| = 0.029 at the calibrated gate, 0.043 forced open; on main+P0 it breaches the 5%
floor (`test-safety-floor.R:712`, 2/4 reps in isolation, 2/2 in suite).

The design decision only you can make: **which surface is the invariant about?** A genuinely
missing cell has no truth to pin, so surface B (truth hidden) is what real users receive —
that argues the invariant and the test belong on B, evaluated via
`predict(fit, .mask_observed_idx = splits$val_idx)`. But even on B, the P0-branch gate opened
for a GNN that did not improve on BM (0.3315 vs 0.3248) — so the gate-acceptance guards
(relative-gain floor, half-A/half-B cross-check) likely also need to assert on the delivered
surface, not just be re-measured.

Then land `arc/p0-onto-main` (`f5e8416` — merge done, 7/7 conflicts resolved, suite otherwise
green: covariate-alignment 13/13, henderson-s-inv 18/18, bm-internal 58/58). Without P0 on
main, the SE, covariate, and zi_count claims in any manuscript are false as shipped.

*Scope: the #157 decision is hours of thought; the fix is 1–2 focused days
(calibration/predict context + test surface + guard re-check); the land is already prepared.*

### 0.2 Fix the public claims the audit flagged — **DONE** (9-edit bundle in #158)

- **BLOCKING**: split the `gnn-architecture.Rmd` §5 SE table row by path (single-column
  `bm_internal.R` = "exact under BM" is defensible; joint/threshold-joint Henderson path is
  not, pre-P0).
- **5 SHOULD-FIX**: covariate row-alignment silence in `README.md` + `getting-started.Rmd`;
  vague SE claim in `mixed-types.Rmd`; leak-provenance caveat on the v0.9.3 ablation numbers
  duplicated in `NEWS.md` and `gnn-architecture.Rmd` §3.5.
  Full text: `docs/dev-log/arc/2026-08-15-p0-claim-gate-findings.md`.

*Scope: half a day, prose only, most wording already drafted in the findings file.*

### 0.3 `--as-cran` on the landed result — **DONE** (0 errors / 0 warnings / 1 known note)

Never run on the merged branch (deliberately skipped while the blocker was unexplained).
Expect the pre-existing 2 WARN / 3 NOTE and nothing new. ~45 min local.

---

## TIER 1 — The GNN: answer the question a reviewer will ask first

**"What does the GNN actually buy?"** The only multi-rep answer on record is unfavourable
AND unusable: the v0.9.3 ablation (60 reps, one regime) found `pigauto_OFF` beat
`pigauto_ON` (z-RMSE 1.038 vs 1.056) — and it was measured on the leaky training loop, so it
cannot be cited in either direction.

### 1.1 Re-run the per-type bench suite post-P0 — **DONE (5/8) 2026-08-16**

Ran locally (throttled), not Totoro — each driver is ~15–45 min, well under the D-139 line.
**Result: the leak had no detectable effect.** Paired per-rep deltas are below 1.1 MCSE in
every completed type (continuous 0.69, binary 1.05, ordinal 0.08, count 0.39, categorical
0.45); categorical reproduces 6/7 scenarios exactly, because its gate closes completely.
Full write-up + caveats (power bounds effects only above ~0.01–0.02; this is not "the leak
was harmless"): `docs/dev-log/2026-08-16-bench-rerun-results.md`.
Remaining: proportion, zi_count, multi_proportion.

### 1.2 Build the regime map (#135)

Where does the GNN add value, as a function of phylogenetic signal × n × missingness ×
trait type, with gate-firing diagnostics (how often, how wide, which traits). Multi-seed,
gated like the MI campaign. Two publishable outcomes, both fine:

- GNN lifts in identifiable regimes → the paper claims those regimes and shows the map.
- OFF ≥ ON survives the corrected pipeline → the honest headline becomes the **safety-gated
  architecture** ("an ML correction that provably cannot hurt you, and a gate that tells you
  when phylogeny alone suffices") plus the mixed-type/UQ/design contributions. That is still
  a paper; it is a different paper. Decide the framing AFTER 1.1–1.2, not before.

### 1.3 Conformal coverage campaign

The documented claim is ≥95% marginal coverage under split-conformal assumptions. The only
multi-rep number on record is **0.884–0.887** (v0.9.3 ablation — leak-tainted, and measured
with ~10–15 validation cells per trait, where the conformal quantile is extremely noisy).
Needed: coverage vs n (say 80 / 300 / 1000), per trait type, multi-seed, with n_val recorded.
If small-n undercoverage is real, the fix may be documentation (minimum-n guidance) rather
than code — but it must be measured. *Totoro; D-139.*

---

## TIER 2 — The latent multivariate machinery (Level C)

### 2.1 Decide: restore EM or scope the claim

`max_iter = 0` ships. The "joint multivariate baseline" is a **single-pass Henderson-init**
fit — cross-trait Σ is never iterated, and the threshold-model prior is the plug-in N(0,1)
that Phase 6 EM was supposed to replace. Two defensible routes:

- **Scope the claim** (cheap): the paper says "single-pass joint MVN with Henderson
  initialisation" and never says "EM". P1-5 (joint-Σ claim-gate) closes by wording.
- **Restore EM** (a dedicated numerical G0): only with a known-Σ recovery simulation gating
  it — the reason it was disabled matters and must be re-established before re-enabling.

### 2.2 Known-truth recovery simulation for the multivariate baseline

Nothing currently demonstrates the joint path recovers a known Σ, or that per-tip SEs are
calibrated (SE coverage vs nominal), or that the binary/ordinal liability decode
`Φ(μ/√(1+σ²))` produces calibrated probabilities once P0's per-tip SE feeds it. This is the
recovery-to-truth lane the repo doctrine already demands, applied to the newest machinery.
One ADEMP-style sim: BM-evolve K correlated traits with known Σ and λ ladder, mask, recover.
*Design first (simulation-design skill), then Totoro; D-139.*

### 2.3 Make the headline features compose — P1-8

Covariates are **ignored by the joint / threshold-joint baselines**. A user who supplies
climate covariates and has ≥2 continuous traits silently gets a covariate-free baseline.
Either thread covariates through (larger), or detect-and-message the fallback (small, honest,
publishable as a documented limitation).

### 2.4 The remaining P1 items that touch manuscript claims

- **P1-9** zi_count zeros excluded from val/test → zi_count evaluation numbers are computed
  on non-zeros only; fix before publishing any zi_count metric.
- **P1-11** tree-MI has no conformal draws → `multi_impute_trees()` draws are not the
  documented conformal method; fix or scope the tree-uncertainty claim.
- **P1-12** multi-obs phylo-signal gate no-op → the safety story has a hole in multi-obs
  mode specifically.
- **P1-7 / P1-6** `$se` is three different objects / conformal wording — documentation
  consolidation, half a day.
- **P1-13** type detection mislabelling counts — small correctness fix in `read_traits`.

### 2.5 OVR categorical head-to-head

The 72% Trophic.Level number in the docs is **BACE's**, not pigauto's. Since the in-house
solver replaced Rphylopars, the OVR path has no post-solver accuracy benchmark against LP,
against a true multinomial, or against BACE. One bench script, AVONET + one simulated grid.

---

## TIER 3 — Packaging and comparisons

- **Head-to-head, regime-labelled**: Rphylopars, BACE, missForest(+eigenvectors), on the
  corrected pipeline. Some exist (`bench_avonet_missingness`, TabPFN follow-up) — re-run
  post-P0 and label regimes. Reviewers will not accept self-benchmarks alone.
- **CRAN**: next cut must drop Suggests `BACE` or wait for BACE on CRAN (still 404 as of
  2026-08-15). A JOSS/MEE software paper wants the CRAN version citable.
- **Reproducibility pack**: the three #157 scripts under `docs/dev-log/arc/` are the pattern —
  every number in the manuscript should have a committed script that regenerates it.

---

## Suggested paper shape (opinion, clearly labelled as such)

Lead with what is defensible **today** and no competitor combines: unified mixed-type
imputation with distribution-free conformal UQ, tree-uncertainty propagation via Rubin
pooling, and active sampling design (`suggest_next_observation()` is the novel hook no other
phylogenetic imputer exposes). Present the GNN as a **safety-gated correction** whose
contribution is regime-mapped by Tier 1 — whatever that map says. Do not publish: any
pre-P0 benchmark number; "exact under BM" for the joint path; a pigauto-vs-BACE headline
from the wrap coverage numbers (S3 was wrap-config vs wrap-config); an EM claim while
`max_iter = 0`.

**Venue fit**: MEE application/methods paper is the natural target; JOSS as a fallback for
the software alone if the methods evidence takes longer.

---

## Sequencing (dependencies, not dates)

```
#157 decision ──→ P0 fix+land ──→ --as-cran ──→ claim fixes (0.2)
                        │
                        ├─→ 1.1 bench re-runs ──→ 1.2 regime map ──→ paper framing decision
                        │                    └──→ 1.3 coverage campaign
                        │
                        └─→ 2.2 recovery sim (needs per-tip SE landed)
2.1 EM decision, 2.3–2.5, Tier 3: parallel, unblocked
```

Compute: 1.1–1.3 and 2.2 are Totoro campaigns — **D-139: written estimate + pre-run test +
approval before each**. Everything in Tier 0 and the small P1 fixes are local Mac work.

The single highest-leverage next action when you return: **the #157 surface decision.** It
unblocks the land, the safety claim, and every benchmark that follows.
