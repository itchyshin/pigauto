#!/usr/bin/env bash
# Codespace setup for pigauto scaling benchmarks (n up to 10,000 species)
set -euo pipefail

echo "=== Installing pigauto dependencies ==="
Rscript -e 'install.packages("pak", repos = "https://pak.r-lib.org/stable")'
Rscript -e 'pak::pkg_install(c("ape", "phytools", "torch", "Rphylopars", "cli", "future.apply", "Matrix", "RSpectra", "igraph"))'

echo "=== Installing torch runtime (CPU) ==="
Rscript -e 'torch::install_torch()'

echo "=== Installing pigauto from local source ==="
Rscript -e 'pak::pkg_install(".", dependencies = TRUE)'

echo "=== Sanity check ==="
Rscript -e 'library(pigauto); library(torch); cat("pigauto:", as.character(packageVersion("pigauto")), "\n"); cat("torch tensor test:", sum(torch_ones(3)$to(dtype=torch_float())), "\n")'

echo "=== Setup complete ==="
echo "To run the scaling benchmark:"
echo "  Rscript script/bench_scaling_v031.R > script/bench_scaling_v031.log 2>&1 &"
echo "  tail -f script/bench_scaling_v031.log"
