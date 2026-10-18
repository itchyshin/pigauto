#!/usr/bin/env bash
#SBATCH --job-name=pigauto-validation
#SBATCH --account=def-snakagaw_cpu
#SBATCH --time=00:35:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
# Pilot only. Production array resources must be bound to measured pilot approval.
: "${VALIDATION_ROOT:?Set the exact frozen source directory on persistent project storage}"
cd "$VALIDATION_ROOT"
bash script/validation_20261018/pilot.sh
