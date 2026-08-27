#!/bin/bash
# One-time bootstrap on Totoro for AVONET panel comparators.
# Installs remotes, Rphylopars, and in-tree BACE (if BACE/ synced).
#
# Usage (on Totoro):
#   cd ~/pigauto_gnn_evidence_campaign
#   bash script/bootstrap_gnn_evidence_avonet_panel.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

echo "[bootstrap-avonet] ROOT=${ROOT}"

Rscript -e '
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  if (!requireNamespace("Rphylopars", quietly = TRUE)) {
    message("Installing Rphylopars from CRAN ...")
    install.packages("Rphylopars", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  cat("Rphylopars:", requireNamespace("Rphylopars", quietly = TRUE), "\n")
'

if [[ -d "${ROOT}/BACE" && -f "${ROOT}/BACE/DESCRIPTION" ]]; then
  echo "[bootstrap-avonet] Installing in-tree BACE ..."
  Rscript -e "devtools::install('${ROOT}/BACE', upgrade=FALSE, quiet=TRUE)"
else
  echo "[bootstrap-avonet] BACE/ not synced — skip local install (panel will report BACE unavailable)"
fi

Rscript -e 'cat("BACE:", requireNamespace("BACE", quietly=TRUE),
                  " MCMCglmm:", requireNamespace("MCMCglmm", quietly=TRUE), "\n")'
echo "[bootstrap-avonet] done"
