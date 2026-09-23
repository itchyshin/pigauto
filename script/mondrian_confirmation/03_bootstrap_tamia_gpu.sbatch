#!/bin/bash
# One-off, project-local R/Torch runtime bootstrap for the Stage-C GPU smoke.
# Submit from Tamia after source is staged under $PIGAUTO_STAGEC_ROOT/source.
# This is deliberately a Slurm CPU job: no package compilation/download occurs
# on the login node.
#SBATCH --account=aip-snakagaw
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

set -euo pipefail
PIGAUTO_STAGEC_ROOT="${PIGAUTO_STAGEC_ROOT:-/project/aip-snakagaw/pigauto-stagec-gpu}"
module load httpproxy/1.0 r/4.5.0 cuda/12.6
export LD_LIBRARY_PATH="$EBROOTCUDA/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="$PIGAUTO_STAGEC_ROOT/Rlib"
mkdir -p "$R_LIBS_USER" "$PIGAUTO_STAGEC_ROOT/logs"

Rscript --vanilla - <<'RS'
repos <- c(CRAN = "https://cloud.r-project.org")
pkgs <- c("torch", "ape", "ggplot2", "rlang")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = repos,
                   dependencies = c("Depends", "Imports", "LinkingTo"))
}
stopifnot(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
if (!torch::torch_is_installed()) torch::install_torch()
stopifnot(torch::torch_is_installed())
cat("R/Torch bootstrap complete\n")
RS

R CMD INSTALL --library="$R_LIBS_USER" "$PIGAUTO_STAGEC_ROOT/source/pigauto"
