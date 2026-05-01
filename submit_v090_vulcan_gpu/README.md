# pigauto Vulcan (PAICE) SLURM bundle — GPU head-to-head at n=9,993

Runs pigauto × BACE head-to-head on the full AVONET dataset (9,993 bird
species, 7 traits) in one SLURM job. Pigauto's GNN trains on GPU
(O(n²) attention benefits from parallel kernels at this scale); BACE
stays on CPU (no GPU support in MCMCglmm). Both coverage types
captured (conformal + MC-dropout).

## Who this is for

PAICE user with `aip-snakagaw` and Vulcan GPU partition access.

## Setup (one-time — assumes `submit_v090_vulcan/README.md` already
completed for the CPU bundle, i.e. pigauto + BACE + torch already in
`~/R/x86_64-pc-linux-gnu-library/4.4`)

If you skipped the CPU bundle and start here:

```bash
ssh snakagaw@vulcan.alliancecan.ca
module load r/4.4.0
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4
echo '.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")' >> ~/.Rprofile
R -e 'install.packages("pak", repos = "https://cran.r-project.org")'
R -e 'pak::pak("local::~/pigauto_src")'
R -e 'pak::pak("local::~/BACE_src")'
R -e 'torch::install_torch()'   # ~200 MB libtorch download
R -e 'library(pigauto); library(BACE); library(torch); cat("all loaded\n")'
```

**GPU-specific verification** (do this on a GPU compute node via
`salloc` before running the full job — saves a failed 12-hr reservation):

```bash
salloc --gres=gpu:1 --time=00:30:00 --mem=8G --account=aip-snakagaw
R -e 'library(torch); cat("cuda_is_available:", torch::cuda_is_available(), "\n"); cat("cuda_device_count:", torch::cuda_device_count(), "\n"); x <- torch::torch_randn(1000L, 1000L, device = "cuda"); cat("GPU tensor alloc + matmul wall (ms):", system.time(torch::torch_matmul(x, x))[["elapsed"]] * 1000, "\n")'
exit
```

Expect `cuda_is_available: TRUE`, `cuda_device_count: 1`, and a
sub-second matmul wall. If any of those fail, torch's libtorch install
isn't wired to CUDA — reinstall with `torch::install_torch(reinstall = TRUE)`.

## Per-campaign upload + submission

```bash
# From your Mac:
rsync -avz submit_v090_vulcan_gpu/ snakagaw@vulcan.alliancecan.ca:~/pigauto_vulcan_gpu/
```

On Vulcan:

```bash
cd ~/pigauto_vulcan_gpu
sbatch submit_avonet9993_gpu.sh
squeue -u snakagaw     # watch it queue + start
```

## Expected timeline

| stage | runs on | ~wall |
|---|---|---:|
| preprocess + build_phylo_graph (n² cophenetic) | CPU | ~5 min |
| fit_baseline (Rphylopars)                      | CPU | 5–15 min |
| pigauto GNN training (500 epochs)              | **GPU** | 30–60 min |
| predict + MI draws                             | GPU | ~5 min |
| BACE MCMCglmm on 7 traits                      | CPU | **4–8 hr** |
| Total                                          |    | **~5–10 hr** |

BACE dominates wall-clock. SLURM budget is 12 hr with safety margin.

## Smoke test first (recommended)

Before burning ~8 hr of GPU budget, verify the bench plumbing works on
a 500-species subset:

```bash
cd ~/pigauto_vulcan_gpu
PIGAUTO_SMOKE=1 Rscript Rscripts/run_avonet9993_bace.R
```

Takes ~5–10 min on a GPU node, much less on a login node. Produces
the same output files; RDS has `smoke = TRUE` attribute.

## What comes back

`~/pigauto_vulcan_gpu/bench_avonet9993_bace.{rds,md}` + the per-job
`run_avonet9993_bace_${SLURM_JOB_ID}.log`, archived as
`results-avonet9993-${SLURM_JOB_ID}.tar.gz` by the post-run `tar` step.

Pull home:

```bash
# From your Mac:
bash submit_v090_vulcan_gpu/rsync_results_back.sh snakagaw@vulcan.alliancecan.ca
```

Then render the HTML report locally:

```bash
Rscript script/make_bench_avonet9993_bace_html.R
```

Report lands at `pkgdown/assets/dev/bench_avonet9993_bace.html` and
linked from the pkgdown validation suite.

## If the GPU queue is slow

If `squeue -u snakagaw` shows the job stuck in `PD (Priority)` or
`PD (Resources)` for > 1 hr, options:

1. **Wait** — GPU queues on PAICE have wait times of hours to ~a day
   during peak.
2. **Cancel + resubmit to CPU** — edit the SLURM script:
   - Remove `#SBATCH --gres=gpu:1`
   - Bump `#SBATCH --time=24:00:00`
   - `sbatch submit_avonet9993_gpu.sh`
   - Expect total wall ~15–20 hr (pigauto GNN goes to CPU, BACE
     unchanged).

## Follow-ups (not in this bundle)

- **Multi-seed AVONET 9,993** — redo at seed ∈ {2026, 1, 2, ...} to get
  robustness bars. Add an array loop.
- **PanTHERIA × BACE** — same pattern, different taxon (mammals).
- **Self-supervised transformer pretraining** — separate bundle, weeks
  of compute.
