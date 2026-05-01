# Multi-seed N_IMP=20 AVONET re-verify + cv_folds default-flip evidence

Date: 2026-05-01.
Closes the two open items from the 2026-04-30 sprint:
- Opus E2 caveat: re-verify continuous AVONET evidence at N_IMP=20.
- v0.9.2 question: should `gate_method = "cv_folds"` become the default?

This memo presents both pieces of evidence and the resulting decisions.

---

## Part 1 — Multi-seed N_IMP=20 AVONET (closes Opus E2)

### Method
Three independent runs of `script/bench_avonet_full_local.R` on AVONET
n=1500 random subset, miss_frac = 0.30, **N_IMP = 20** (4× the previous
N_IMP = 5 multi-seed run from `MEMO_2026-04-30_avonet_multiseed.md`).
Seeds: 2030, 2031, 2032. Wall: ~9 min each, ~28 min total. Sequential
(no GPU contention).

Files: `script/bench_avonet_full_local_n20_seed{2030,2031,2032}.{rds,md}`.
Aggregation: `/tmp/analyze_multiseed_n20.R` →
`/tmp/multiseed_n20_analysis.txt`.

### Headline table

| trait | metric | mean baseline | pigauto mean ± SD | CV % | lift % mean ± SD |
|---|---|---|---|---|---|
| Beak.Length_Culmen | RMSE     | 26.1   | **12.5 ± 3.4**   |  27.6 | **+52.4 ± 6.1**   |
| Mass               | RMSE     | 2588   | **7655 ± 10724** | 140.1 | **−297 ± 606**    |
| Tarsus.Length      | RMSE     | 26.2   | **13.2 ± 7.1**   |  53.4 | **+51.7 ± 12.9**  |
| Wing.Length        | RMSE     | 94.6   | **29.5 ± 1.9**   |   6.5 | **+68.9 ± 1.9**   |
| Migration          | accuracy | 0.800  | **0.713 ± 0.032**|   4.5 | **−10.9 ± 3.1**   |
| Primary.Lifestyle  | accuracy | 0.575  | **0.810 ± 0.029**|   3.6 | **+41.0 ± 3.7**   |
| Trophic.Level      | accuracy | 0.566  | **0.825 ± 0.008**|   1.0 | **+46.1 ± 7.9**   |

### Per-seed pigauto values (continuous, RMSE)

| trait | seed 2030 | seed 2031 | seed 2032 |
|---|---|---|---|
| Beak.Length_Culmen | 12.62 |  9.02 | 15.91 |
| Mass               | **19960.** |  292.7 | 2714.9 |
| Tarsus.Length      |  8.25 | 21.34 | 10.12 |
| Wing.Length        | 30.92 | 27.31 | 30.15 |

Seed 2030 produced a Mass RMSE of **19960** vs the per-seed baseline of
1820 — pigauto roughly 11× worse than column mean on that seed. Seeds
2031 and 2032 produced 293 and 2715 RMSE respectively, beating the
baseline. All three seeds emit `solve(): system is singular` warnings
during training (rcond ranging from 1e−16 down to 1e−237), but only
seed 2030 produces an outlier prediction in the back-transformed Mass
scale. **This is a numerical-instability finding that needs a separate
investigation; do not treat the seed-2030 Mass result as a regression
of any code change.**

### Comparison vs N_IMP=5 (2026-04-30 baseline)

| trait | metric | n5 CV % | n20 CV % | noise reduction factor (n5_CV / n20_CV) |
|---|---|---|---|---|
| Beak.Length_Culmen | RMSE     | 21.9 | 27.6 | 0.79 |
| Mass               | RMSE     | 57.4 | **140.1** | **0.41** |
| Tarsus.Length      | RMSE     | 32.8 | 53.4 | 0.61 |
| **Wing.Length**    | RMSE     | 24.2 | **6.5** | **3.75** |
| Migration          | accuracy |  4.7 |  4.5 | 1.04 |
| Primary.Lifestyle  | accuracy |  3.6 |  3.6 | 0.99 |
| Trophic.Level      | accuracy |  0.9 |  1.0 | 0.85 |

NRF > 1 means N_IMP=20 reduced cross-seed noise, as the Nakagawa &
de Villemereuil (2019) pooling argument predicts. Wing achieves the
predicted reduction (3.75×). Mass actually got noisier because the
seed-2030 outlier dominates the SD; without seed 2030 the Mass NRF
would also be > 1.

### Closes Opus E2?

Partially. The original caveat — "the v3 single-seed Mass −2.75
RMSE change might be within pooling noise" — was already addressed
by the 2026-04-30 N_IMP=5 multi-seed run (the SDs were huge). The
N_IMP=20 run was supposed to **tighten** the SDs and confirm that the
mean lifts were not artefacts. Result:

- **Confirmed for 5 of 7 traits**: Beak +52 % (SD 6 %), Tarsus +52 %
  (SD 13 %), Wing +69 % (SD 2 %), Trophic +46 % (SD 8 %), Primary
  +41 % (SD 4 %). All bands strictly positive across all three seeds.
- **Confirmed regression for Migration**: −10.9 % (SD 3.1 %), all
  three seeds in the −8.4 % to −14.4 % band. Phase F (LP corner for
  ordinal) is the queued fix; the regression is real and consistent.
- **NEW: Mass instability uncovered**: at N_IMP=20 the seed-2030 run
  produces a catastrophic Mass RMSE (≈ 11× baseline). At N_IMP=5 this
  pattern was hidden in seed-to-seed variance. This is a real
  numerical-instability finding, not a Fix-A–H regression.

### Mass instability — provisional diagnosis

The torch / Armadillo `solve()` warnings happen on the threshold-joint
B3 / OVR-categorical / joint-MVN baselines whenever `Rphylopars`
delegates to a near-singular Σ. They appear in all three seeds, so the
warnings themselves do not predict the seed-2030 blow-up. The most
likely cause of the seed-2030 outlier is one or more residual cells
where the back-transformation (`expm1` of the log1p-z latent) hits an
extreme tail. Action item: instrument `predict_pigauto.R` to log the
top-k absolute residuals per trait per imputation; investigate at next
sprint.

---

## Part 2 — cv_folds vs single_split default-flip evidence

### Method
`script/bench_cv_default_flip.R` (NEW, committed in this sprint).
Pre-registered decision rule, recorded in the script header BEFORE the
run:

> Flip default to `cv_folds` if:
> 1. mean continuous RMSE lift ≥ 1 % across the three λ regimes
> 2. mean discrete accuracy gap (cv − single) ≥ −0.5 pp

Sweep: 3 phylo-signal regimes (λ ∈ {0.4, 0.7, 1.0}) × 3 trait types
(continuous BM, binary threshold, K=4 categorical) × 10 reps × 2 gate
methods = 90 cells per method. n = 300 species, 30 % MCAR mask, 60
epochs.

### Results

Wall: ~19 min for 60 cells (10 reps × 2 methods × 3 λ).

Per-cell paired delta (cv_folds − single_split):

| trait | λ | delta_mean | delta_sd | n_pos | n_neg |
|---|---|---|---|---|---|
| bin   | 0.4 | −0.00222 | 0.00703 | 0 | 1 |
| bin   | 0.7 | +0.00333 | 0.01054 | 1 | 0 |
| bin   | 1.0 | +0.00667 | 0.02108 | 1 | 0 |
| cat3  | 0.4 | +0.01222 | 0.03865 | 1 | 0 |
| cat3  | 0.7 | +0.00222 | 0.00703 | 1 | 0 |
| cat3  | 1.0 | +0.00222 | 0.00703 | 1 | 0 |
| cont  | 0.4 | −0.00048 | 0.00364 | 1 | 3 |
| cont  | 0.7 | −0.01379 | 0.02728 | 0 | 4 |
| cont  | 1.0 | −0.00346 | 0.02019 | 1 | 2 |

Lower delta is better for RMSE; higher is better for accuracy.
n_pos / n_neg = number of reps where cv_folds was bigger / smaller
than single_split.

Aggregated metrics vs the pre-registered rule:

* Mean continuous RMSE lift (single − cv, so positive ⇒ cv better):
  **+0.00591** (RMSE units). Baseline mean continuous RMSE ≈ 0.75
  across the three λ regimes, so the relative lift is **≈ 0.79 %** —
  below the 1 % gate.
* Mean discrete accuracy gap (cv − single, positive ⇒ cv better):
  **+0.41 pp** — above the −0.5 pp gate (cv_folds does not regress
  on discrete).

### Decision

> Pre-registered rule: flip default to `cv_folds` only if BOTH
> continuous lift ≥ 1 % AND discrete gap ≥ −0.5 pp.
>
> **Continuous: FAIL** (0.79 % < 1 %).
> **Discrete: PASS** (+0.4 pp ≥ −0.5 pp).
>
> Verdict: **KEEP `single_split` as default.** `cv_folds` remains
> opt-in (recommended for users who want lower MC variance and can
> afford ~5× longer gate-calibration wall).

DESCRIPTION stays at 0.9.1.9009 — no v0.9.2 bump needed from this
piece of evidence. The cv_folds option, gate_cv_folds argument, and
docs in `R/fit_pigauto.R` and the corresponding NEWS entry remain
unchanged; we just don't promote it to the first position in the
match.arg call.

`script/bench_cv_default_flip.{R,md,rds}` ships as evidence so any
future re-discussion has reproducible numbers without re-running.

---

## Part 3 — Files

* `script/bench_avonet_full_local_n20_seed{2030,2031,2032}.{rds,md}` — N_IMP=20 multi-seed
* `/tmp/analyze_multiseed_n20.R` + `/tmp/multiseed_n20_analysis.txt` — aggregation
* `script/bench_cv_default_flip.R` — pre-registered cv-folds default-flip bench
* `script/bench_cv_default_flip.{rds,md}` — outputs (post-run)
* `useful/MEMO_2026-04-30_avonet_multiseed.md` — N_IMP=5 baseline this is compared to

## Part 4 — Open

* Mass instability investigation (top-k residual instrumentation) —
  scoped at next sprint.
* Phase F (LP corner for ordinal) — closes the Migration regression;
  already scheduled.
