# pigauto v0.9.1 — Vulcan (PAICE) SLURM bundle for the calibration grid

Runs the pigauto calibration grid across 4 phylogenetic-signal scenarios,
3 missingness mechanisms, 5 trait types, and N replicates. At the default
config (N=50 reps, ~3 min per fit on CPU) the total work is ~3000 fits
and parallelises cleanly as a SLURM array.

## Who this is for

You have PAICE access with `aip-snakagaw` on Vulcan.

## Setup (one-time)

```bash
# Login
ssh snakagaw@vulcan.alliancecan.ca

# Load R and install pigauto + torch (once)
module load r/4.4.0

# User-lib (CC default library is read-only)
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4
echo '.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")' >> ~/.Rprofile

# Install pigauto (needs to be accessible from login node)
# Either: from GitHub (need GitHub PAT) -- see Narval setup in ../submit_v090_cloud/README.md
# Or: rsync your local pigauto source into ~/pigauto_src/ and:
R -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4"); pak::pak("local::~/pigauto_src")'

# Install torch libtorch binaries on login node (compute nodes lack internet)
R -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4"); torch::install_torch()'

# Install BACE (needed by Rscripts/run_calibration_cell.R for sim_bace())
# BACE lives in-tree at pigauto/BACE; rsync it too or install from GitHub.
R -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4"); pak::pak("local::~/BACE_src")'
```

## Per-campaign upload + submission

```bash
# From your Mac:
rsync -avz submit_v090_vulcan/ snakagaw@vulcan.alliancecan.ca:~/pigauto_vulcan/

# On Vulcan:
cd ~/pigauto_vulcan
# Edit the array size and account if needed:
sed -i 's/aip-<PI>/aip-snakagaw/g' submit_calibration_grid.sh

# Submit
sbatch submit_calibration_grid.sh

# Monitor
squeue -u snakagaw
```

Each array task runs one (scenario, mechanism, trait_type) cell across all
`N_REPS` replicates and saves a per-cell RDS. When the whole array
completes, rsync the results home:

```bash
# From your Mac:
bash submit_v090_vulcan/rsync_results_back.sh snakagaw@vulcan.alliancecan.ca
```

## Expected runtime

- Per single fit (n=150, epochs=500): ~3 min CPU
- Per array task (50 reps × 1 cell): ~2.5 hr
- Array size = 4 scenarios × 3 mechanisms × 5 trait types = 60 tasks
- Wall time at 20-way parallelism: ~7.5 hours
- If Vulcan's queue has capacity for all 60 concurrent: ~2.5 hours

## What comes back

`~/pigauto_vulcan/results/cell_*.rds` — one per array task, tidy
`data.frame` with columns `metric`, `value`, `n_cells`, `scenario`,
`mechanism`, `trait_type`, `cell_id`, `rep`, `seed`, `wall_s`. The
`metric` column carries `coverage95`, `rmse`/`nrmse` (continuous/count),
or `accuracy` (discrete) per rep.

Merge locally after `rsync_results_back.sh`:

```r
files <- list.files("submit_v090_vulcan/returned",
                    pattern = "^cell_.*\\.rds$", full.names = TRUE)
grid  <- do.call(rbind, lapply(files, readRDS))
saveRDS(grid, "submit_v090_vulcan/returned/calibration_grid.rds")
```

Then feed the grid into a pkgdown calibration report (separate follow-up).

## GPU alternative (future)

Vulcan has A100/H100 partitions. For pigauto's GNN training, moving from
CPU to GPU gives ~10-50× speedup on the training stage (~5 min of each
fit). Not used in this bundle — adding `--gres=gpu:1` and the right
module loads would be a small follow-up once the CPU baseline works.
