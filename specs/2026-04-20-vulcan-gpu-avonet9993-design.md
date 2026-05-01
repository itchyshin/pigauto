# Vulcan GPU bundle — AVONET 9,993 × pigauto + BACE head-to-head

**Status:** Draft (awaiting user review)
**Author:** Claude Opus 4.7 (working with Shinichi Nakagawa)
**Date:** 2026-04-20

---

## Problem

All pigauto-vs-BACE comparisons to date are on **AVONET 300** (bundled
subset). That's tractable but underwhelming as a paper claim — the
Jones et al. 2009 and AVONET 2022 datasets have 9,993 species, which is
where pigauto's scaling story and where BACE's O(n²) MCMC chain start
diverging in practice. We need a **single apples-to-apples comparison
at full scale**.

CPU is infeasible: pigauto at n=9,993, 500 epochs on a Vulcan CPU node
is ~4–8 hr; BACE is ~10–20 hr of MCMCglmm chain. Running both on the
same SLURM job × CPU-only = a day per run, with GPU-capable nodes
sitting idle.

Vulcan has A100 and H100 partitions. Pigauto's GNN training is the
only compute-intensive stage (baseline fit is CPU-only Rphylopars),
and the GNN's O(n²) self-attention is precisely what GPUs accelerate
10–50× at n=10k. BACE stays on CPU in both paths (MCMCglmm has no GPU
support).

## Goal

Ship a **second Vulcan bundle** — `submit_v090_vulcan_gpu/` — that
submits one SLURM job requesting a GPU, runs pigauto on GPU for the
GNN stage, runs BACE on CPU (single chain), and produces a head-to-head
comparison RDS + report on the full AVONET 9,993 dataset.

Both coverage types captured (conformal + MC-dropout) per the standing
rule.

## Non-goals

- **No new R/ code.** Pigauto already routes to GPU automatically when
  `CUDA_VISIBLE_DEVICES` is set; the CPU-only Vulcan bundle explicitly
  unsets it. The GPU bundle simply doesn't.
- **No multi-seed.** Single seed for the MVP; multi-seed is a
  follow-up.
- **No distributed multi-GPU.** One GPU is enough for n=10k.
- **No BACE-on-GPU attempt.** MCMCglmm is C++/CPU only.
- **Not a replacement** for `submit_v090_vulcan/` (the calibration
  grid). Different question, different bundle.

## Design decisions

### 1. New bundle at `submit_v090_vulcan_gpu/`

Mirrors the existing CPU bundle structure:

```
submit_v090_vulcan_gpu/
├── README.md
├── submit_avonet9993_gpu.sh          — SLURM script, GPU job
├── Rscripts/
│   └── run_avonet9993_bace.R          — the actual bench driver
└── rsync_results_back.sh              — pull results back
```

### 2. SLURM parameters (`submit_avonet9993_gpu.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=pigauto_avonet9993_gpu
#SBATCH --account=aip-snakagaw
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=pigauto-avonet9993-%j.out

module load r/4.4.0 cuda/12 cudnn/9
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
# Leave CUDA_VISIBLE_DEVICES at SLURM's default (0)

Rscript Rscripts/run_avonet9993_bace.R
```

Runtime budget: 12 hr is comfortable. Expected actual ~4-8 hr.

### 3. Data acquisition

AVONET 3.0 + matching BirdTree (9,993 species) is already bundled in
pigauto as `avonet_full` + `tree_full` (see `script/bench_avonet_missingness.R`
which used it for the n=10k scaling claim). The GPU bench reuses those
bundled datasets directly — no download step.

### 4. Methods compared

| method | runs on | runtime (est.) |
|---|---|---|
| `mean_baseline` | CPU | seconds |
| `pigauto_default` | GPU | ~15–30 min training + ~5 min predict |
| `pigauto_em5` | GPU | ~30–60 min (Phase 6 EM × 5) |
| `bace_default` | CPU | ~4–8 hr MCMCglmm × 7 traits |

BACE run with `requireNamespace("BACE", quietly = TRUE)` + `tryCatch()`
— if it fails (memory, time), bench still produces pigauto-only
output, consistent with the CPU bundle.

### 5. Splits

Single seed 2026, MCAR 30% per trait. Same convention as Phase 8 AVONET
head-to-head for directly-comparable results.

### 6. Metrics (per trait, per method)

Continuous (`Mass`, `Beak.Length_Culmen`, `Tarsus.Length`, `Wing.Length`):
- RMSE
- Pearson r
- **coverage95_conformal** (from pigauto `pred$conformal_lower/upper`)
- **coverage95_mcdropout** (from `pred$imputed_datasets` with M=20)

Discrete (`Trophic.Level`, `Primary.Lifestyle`, `Migration`):
- Accuracy
- Log-loss
- Brier score
- Coverage via top-k credible set (mass ≥ 0.95) from MI draws

BACE outputs point predictions only (no interval) → coverage N/A for
that method.

### 7. Output + integration

- `submit_v090_vulcan_gpu/results/avonet9993_bace.rds` — tidy data.frame,
  same schema as the CPU bench RDS
- `avonet9993_bace.md` — human-readable summary
- `make_bench_avonet9993_bace_html.R` — HTML with winner-highlighted
  pivot table; lands at `pkgdown/assets/dev/bench_avonet9993_bace.html`
- Row in `make_validation_suite_html.R` for the new bench
- NEWS entry in the v0.9.1.9000 (dev) section

### 8. What we'll learn

Four paper-ready numbers:

1. **Does pigauto scale end-to-end at n=10k on real data?** Yes/no
   with timing evidence.
2. **Does pigauto beat BACE on continuous RMSE at this scale?** The
   v0.9.0 n=300 claim was 0.89–0.98 r (pigauto) vs 0.97 r (BACE); we'll
   see what happens at 33× scale.
3. **Does pigauto's categorical accuracy hold?** The v0.9.0 n=300
   claim was 77% on Trophic.Level (vs BACE 72%); n=9,993 should be
   stronger signal.
4. **Do both coverage types hold at large n?** Expect conformal ≈ 0.95
   (the guarantee kicks in cleanly when val cells ≫ 20 per trait);
   MC-dropout ≈ 0.90 (mildly under-covered based on prior benches).

### 9. BACE call at this scale

`BACE::bace()` on 9,993 species × 7 traits with the default 2000-iter
MCMCglmm chain is the slow point. Options:

- **Option A** (shipped default): single chain, 2000 iter, 500 burnin.
  Runtime: ~4–8 hr. Possibly hits memory ceiling on 64 GB — monitor.
- **Option B**: shorter chain (1000 iter, 200 burnin). ~2–4 hr. Less
  converged but fits within `--time=12:00:00`.
- **Option C**: drop BACE — pigauto-only at n=9,993. Save for a
  dedicated BACE run once bench infrastructure is proven.

Spec goes with A + a guard: `tryCatch()` around BACE so failure
doesn't abort the pigauto portion.

### 10. Fallback: if GPU allocation fails

GPU nodes can be scarce. If the job is stuck in the queue for > 1 hr,
user can:
- Cancel + resubmit to CPU partition (edit `--gres` off and bump
  `--time` to `24:00:00`)
- Or wait — the GPU job will eventually schedule.

Script prints clear status at each stage (`[STAGE] preprocess`,
`[STAGE] baseline`, `[STAGE] pigauto train`, `[STAGE] BACE`) so a
stuck stage is visible in the SLURM stdout.

## API surface

No user-facing R changes. New artefacts:

- `submit_v090_vulcan_gpu/README.md` — setup + submit walkthrough
- `submit_v090_vulcan_gpu/submit_avonet9993_gpu.sh` — SLURM script
- `submit_v090_vulcan_gpu/Rscripts/run_avonet9993_bace.R` — driver
- `submit_v090_vulcan_gpu/rsync_results_back.sh` — helper
- `script/make_bench_avonet9993_bace_html.R` — HTML gen (lives in the
  main tree, not the bundle — rendered from results pulled home)

## Testing

- Bundle shell scripts pass `bash -n`.
- R driver passes `parse()`.
- No unit tests — this is an integration bench.
- Smoke check: can run the driver script locally on avonet_full at a
  **tiny** subset (n=500) to verify end-to-end plumbing before pushing
  to Vulcan. Smoke mode triggered by env var `PIGAUTO_SMOKE=1` →
  subsets to first 500 species.

## Documentation

1. Bundle README with one-time setup (mirrors `submit_v090_vulcan/README.md`).
2. NEWS.md bullet:
   > **Vulcan GPU bundle for AVONET 9,993 × BACE head-to-head.** New
   > `submit_v090_vulcan_gpu/` directory ships a GPU-enabled SLURM
   > script that runs pigauto on GPU (GNN training ~15–30 min at n=10k)
   > and BACE on CPU (MCMCglmm ~4–8 hr) on the same splits. Coverage
   > and point metrics for all 7 AVONET traits. Run on PAICE
   > `aip-snakagaw` after the calibration grid (separate job, separate
   > bundle).

## Backward compatibility

Fully additive. No existing files touched.

## Success criteria

- [ ] `submit_avonet9993_gpu.sh` submits cleanly with `sbatch` and
      lands a GPU (`squeue` shows `R` with `gpu:1`).
- [ ] pigauto training stage completes on GPU in < 60 min.
- [ ] Bench writes `.rds` + `.md` even if BACE fails (graceful
      fallback).
- [ ] Report shows both coverage types for pigauto paths on all
      continuous traits.
- [ ] pkgdown validation suite gains a new row on merge.

## Out of scope (deferred)

- Multi-seed AVONET 9,993 (needs M × the compute — queue it separately)
- PanTHERIA × BACE at n=4,000 (same pattern, different taxon; queue
  after AVONET)
- Self-supervised transformer pretraining (weeks-scale research project)
- Distributed multi-GPU training (not needed at n=10k)
