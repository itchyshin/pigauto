# pigauto v0.9.0 — Compute Canada submission bundle

This directory contains everything needed to run two heavy benchmarks on an
Alliance Canada (Compute Canada) cluster. As of the 2025-2026 infrastructure
renewal the general-purpose clusters are **Fir** (replaces Cedar), **Nibi**
(replaces Graham), **Rorqual** (replaces Béluga), **Trillium** (replaces
Niagara), and **Narval** (unchanged). Cedar / Graham / Béluga hostnames no
longer resolve — use the new names below.

The scripts and instructions in this directory use `fir.alliancecan.ca` as
the default; substitute any other active cluster the same way.

## Benchmarks

| Script | Job script | Description | Est. runtime |
|---|---|---|---|
| `Rscripts/bench_avonet_missingness.R` | `submit_avonet.sh` | AVONET full n=9,993 missingness sweep (20/50/80%, 3 methods) | 10-20 hr |
| `Rscripts/bench_scaling_v090_extended.R` | `submit_scaling_7500_10000.sh` | Scaling curve extension: n in {5000, 7500, 10000} | 6-10 hr |

Both scripts write `.rds`, `.md`, and (for scaling) `.png` result files to the
working directory, then the SLURM wrapper archives them into
`results-<jobname>-<JOBID>.tar.gz`.

---

## Prerequisites

- A valid Alliance Canada account with SSH access to a cluster.
- Your SLURM account code (e.g. `def-senak`). If unsure, run `sacctmgr show
  associations user=$USER` on the cluster.

---

## Step-by-step

### 1. Upload this bundle to the cluster

```bash
# From your local machine, inside the pigauto repo:
rsync -avz submit_v090_cloud/ \
  <CC_USERNAME>@fir.alliancecan.ca:~/pigauto_cloud/
```

Replace `<CC_USERNAME>` with your CC username and `fir` with your preferred
cluster (`nibi` / `rorqual` / `narval` / `trillium`).

### 2. SSH in and install pigauto

```bash
ssh <CC_USERNAME>@fir.alliancecan.ca
cd ~/pigauto_cloud
```

Load R and install pigauto from GitHub:

```bash
module load r/4.4.0

# Option A — pak (preferred, resolves all dependencies automatically)
R -e 'install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable/"); pak::pak("itchyshin/pigauto")'

# Option B — remotes (if pak is unavailable)
R -e 'install.packages("remotes"); remotes::install_github("itchyshin/pigauto")'
```

Install `torch` after pigauto (it bundles libtorch, ~500 MB download):

```bash
R -e 'torch::install_torch()'
```

**Important torch note**: CC compute nodes have no internet access. Run
`torch::install_torch()` from a login node BEFORE submitting jobs. The
libtorch binaries land in `~/.local/share/torch` by default and are accessible
from compute nodes via the shared filesystem.

Verify everything works:

```bash
R -e 'library(pigauto); library(torch); cat("OK\n")'
```

### 3. Edit the SLURM scripts

Open each SLURM script and replace `def-<PI>` with your actual account code:

```bash
# In ~/pigauto_cloud/
sed -i 's/def-<PI>/def-senak/g' submit_avonet.sh submit_scaling_7500_10000.sh
```

Or edit manually — the line is `#SBATCH --account=def-<PI>`.

### 4. Submit jobs

```bash
cd ~/pigauto_cloud
sbatch submit_avonet.sh
sbatch submit_scaling_7500_10000.sh
```

Both jobs can run concurrently if sufficient allocation is available.

### 5. Monitor progress

```bash
# Live queue status
squeue -u $USER

# Detailed accounting after completion
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS

# Live stdout of the running job (tailed log)
tail -f pigauto_avonet_missingness-<JOBID>.out
```

The SLURM `.out` file mirrors the Rscript stdout, which prints timestamped
stage checkpoints like `[  42.3s] [preprocess_20pct         ] wall = ...`.

### 6. Retrieve results

On your **local** machine:

```bash
cd /path/to/pigauto  # repo root
bash submit_v090_cloud/rsync_results_back.sh \
  <CC_USERNAME>@fir.alliancecan.ca
```

Archives land in `submit_v090_cloud/returned/`. Unpack:

```bash
cd submit_v090_cloud/returned
for f in *.tar.gz; do tar xzf "$f"; done
```

### 7. Regenerate HTML reports locally

Copy the `.rds` files to `script/` in your local pigauto repo and re-run the
HTML generators:

```bash
# Avonet missingness
cp submit_v090_cloud/returned/bench_avonet_missingness.rds script/
Rscript script/make_avonet_missingness_html.R

# Scaling (merge with local bench_scaling_v090.rds first if desired)
cp submit_v090_cloud/returned/bench_scaling_v090_extended.rds script/
# Then re-run whichever make_bench_scaling_html.R script applies.
```

---

## File layout

```
submit_v090_cloud/
├── README.md                           # this file
├── submit_avonet.sh                    # SLURM: AVONET missingness sweep
├── submit_scaling_7500_10000.sh        # SLURM: scaling curve n={5000,7500,10000}
├── Rscripts/
│   ├── bench_avonet_missingness.R      # self-contained; uses library(pigauto)
│   └── bench_scaling_v090_extended.R  # n_grid extended to 10000
├── data/
│   └── README.md                      # notes: data is bundled in the package
├── returned/                          # rsync_results_back.sh writes here
└── rsync_results_back.sh              # helper to pull .tar.gz files home
```

---

## Troubleshooting

**`torch::install_torch()` hangs or fails on a login node**
Some clusters rate-limit large downloads. Try: `R -e 'options(timeout = 600);
torch::install_torch()'`. If the download still fails, download libtorch
manually on your local machine and rsync `~/.local/share/torch` to the cluster.

**`module load r/4.4.0` not found**
Run `module spider r` to see available R versions and update the `module load`
line in the SLURM scripts accordingly.

**Job hits the time limit before finishing**
Increase `#SBATCH --time` (e.g. `36:00:00` for Avonet). The scripts checkpoint
after every stage, so a restarted run can be modified to skip already-completed
stages. The intermediate `.rds` written during the run already contains partial
results.

**`def-<PI>` account not found / job rejected**
Run `sacctmgr show associations user=$USER` on the cluster to list valid
account codes. Use the `Account` column value (e.g. `def-senak`).

**OOM kill at n=10000**
The scaling script allocates ~800 MB for the cophenetic distance matrix plus
large torch tensors. If 48 GB is insufficient, increase `#SBATCH --mem` to
`64G`. The script calls `graph$D <- NULL; gc()` before GNN training to free
the cophenetic matrix.
