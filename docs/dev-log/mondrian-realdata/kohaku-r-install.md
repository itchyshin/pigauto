# User-space R on kohaku (2026-09-23)

Host: `ssh -o BatchMode=yes snakagaw@kohaku.biology.ualberta.ca` (Ubuntu 24.04.4 LTS,
no sudo, RTX PRO 6000 Blackwell Max-Q Workstation Edition, driver 610.43.02, 32 vCPU).
Everything below was installed user-space under `/home/snakagaw`, capped at 8 CPU
threads (`OMP_NUM_THREADS=8 MAKEFLAGS=-j8 OPENBLAS_NUM_THREADS=1`).

## 1. micromamba + R 4.5 env

```bash
mkdir -p ~/.local/bin
cd /tmp
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest -o micromamba.tar.bz2
tar -xvjf micromamba.tar.bz2 bin/micromamba
mv bin/micromamba ~/.local/bin/micromamba
chmod +x ~/.local/bin/micromamba

export MAMBA_ROOT_PREFIX=~/micromamba
~/.local/bin/micromamba create -y -p ~/micromamba/envs/r45 -c conda-forge \
  "r-base=4.5" r-ape r-devtools r-rlang r-ggplot2 r-nlme r-matrix \
  compilers cmake make pkg-config
```

No edits were made to `~/.bashrc`; the env is called by absolute path
(`~/micromamba/envs/r45/bin/Rscript`) with `PATH=~/micromamba/envs/r45/bin:$PATH`
set per-command so the conda-forge compiler wrappers (`x86_64-conda-linux-gnu-cc`
etc.) are found by R's `Makeconf`. Without that `PATH` prefix, R package
compilation fails with `x86_64-conda-linux-gnu-cc: not found` even though
`compilers` is installed in the env.

- micromamba: 2.9.0
- R: 4.5.3 ("Reassured Reassurer")
- ape: 5.8.1
- devtools: 2.5.2

## 2. torch R package + libtorch (CUDA)

```bash
export PATH=~/micromamba/envs/r45/bin:$PATH
~/micromamba/envs/r45/bin/Rscript -e 'install.packages(c("bit","bit64","torch"), repos="https://cran.r-project.org")'
~/micromamba/envs/r45/bin/Rscript -e 'Sys.setenv(CUDA="12.8"); library(torch); torch::install_torch()'
```

- torch (R package): 0.17.0
- libtorch: `libtorch-shared-with-deps-2.8.0+cu128` (from
  `download.pytorch.org/libtorch/cu128/...`), lantern
  `lantern-0.17.0+cu128+x86_64-Linux` (from `torch-cdn.mlverse.org`)

### CUDA usable: yes, after one fix

First attempt failed at load time:

```
Error: Lantern is not loaded. Please use `install_torch()` to install
additional dependencies.
...
/home/snakagaw/micromamba/envs/r45/lib/R/library/torch/lib/liblantern.so -
libcudart.so.12: cannot open shared object file: No such file or directory
```

Root cause: the downloaded `libtorch-shared-with-deps-2.8.0+cu128` ships its
CUDA shared libraries with content-hash suffixes (e.g.
`libcudart-c3a75b33.so.12`, `libcublas-031ce6c2.so.12`,
`libnvrtc-2bb82d1a.so.12`, ...) whose `SONAME` (verified with `readelf -d`) is
the unsuffixed name (`libcudart.so.12`, `libcublas.so.12`, ...), but no
symlink from the `SONAME` to the hash-suffixed file was created on
extraction. `liblantern.so`'s `RUNPATH` is `$ORIGIN` (the same `lib/`
directory), so the dynamic linker needs a file literally named
`libcudart.so.12` there — the `SONAME` recorded *inside* the hash-suffixed
file is irrelevant to `dlopen`'s file-name search. This is a packaging gap in
that libtorch zip, not a CUDA-version incompatibility with Blackwell (sm_120).

Fix (entirely inside the package's own `lib/` directory, no sudo, no system
changes):

```bash
LIBDIR=/home/snakagaw/micromamba/envs/r45/lib/R/library/torch/lib
cd "$LIBDIR"
for f in *.so *.so.*; do
  [ -f "$f" ] || continue
  soname=$(readelf -d "$f" 2>/dev/null | grep -oP "(?<=Library soname: \[)[^]]+")
  if [ -n "$soname" ] && [ "$soname" != "$f" ] && [ ! -e "$soname" ]; then
    ln -s "$f" "$soname"
  fi
done
```

This created 9 symlinks (`libcudart.so.12`, `libcublas.so.12`,
`libcublasLt.so.12`, `libnvrtc.so.12`, `libnvrtc-builtins.so.12.8`,
`libcufile.so.0`, `libcufile_rdma.so.1`, `libcusparseLt.so.0`,
`libgomp.so.1`), after which `ldd liblantern.so` resolved cleanly and:

```r
library(torch)
torch::cuda_is_available()   # TRUE
torch::cuda_device_count()   # 1
```

### Correctness check

```r
x <- torch_randn(4096, 4096, device = "cuda")
y <- (x$matmul(x))$cpu()
xc <- as_array(x$cpu())
ref <- xc %*% xc
max(abs(as_array(y) - ref))   # 0.0003906971
```

Max abs difference ~3.9e-4 on values with magnitude on the order of tens
(sum of 4096 products of standard-normal draws), i.e. relative error
~1e-5 — the expected level of floating-point disagreement between cuBLAS's
GPU reduction order and R's CPU BLAS, not a correctness failure.

Device name (`nvidia-smi --query-gpu=name`): `NVIDIA RTX PRO 6000 Blackwell
Max-Q Workstation Edition`. `torch::cuda_get_device_properties(0)` did not
print a populated result in this torch build; `nvidia-smi` is the
authoritative device-name source used here.

GPU memory: a neighbour's `llama-server` held ~50.5 GiB of the 97.9 GiB total
throughout (`memory.used` read before starting any GPU work); this run's
peak usage stayed far under the 40 GiB cap (4096x4096 fp32 tensors are ~64
MiB each; the pigauto fit below trains on a 3,000-tip graph, also a small
fraction of the remaining ~47 GiB).

## 3. pigauto install + timing pre-run

```bash
export PATH=~/micromamba/envs/r45/bin:$PATH
~/micromamba/envs/r45/bin/Rscript -e 'install.packages("withr", repos="https://cran.r-project.org")'   # missing Import
mkdir -p ~/pigauto_mondrian_realdata/pigauto_src
tar -xzf ~/pigauto_mondrian_realdata/pigauto_7af133a.tar.gz -C ~/pigauto_mondrian_realdata/pigauto_src
~/micromamba/envs/r45/bin/Rscript -e 'install.packages("~/pigauto_mondrian_realdata/pigauto_src/pigauto", repos=NULL, type="source")'
```

- pigauto: 0.11.0 (tarball `pigauto_7af133a.tar.gz`, staged by Shinichi at
  `~/pigauto_mondrian_realdata/`)
- Extra dependency needed beyond the CAP list: `withr` (pigauto `Imports`).

**Input file note**: the task named
`~/pigauto_mondrian_realdata/input/fishbase-3000-input.rds`, but the file
actually staged there is `~/pigauto_mondrian_realdata/input/fb3000.rds`
(3,000 tips, 6 traits). Used `fb3000.rds` as the input for the pre-run below.

### Pre-run

```bash
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
cd ~/pigauto_mondrian_realdata/pigauto_src/pigauto
Rscript script/mondrian_confirmation/01_run_masked_confirmation.R \
  ~/pigauto_mondrian_realdata/input/fb3000.rds \
  ~/pigauto_mondrian_realdata/results/fishbase3000-kohaku-prerun \
  mcar 20260818 100 split
```

**Estimate written before running**: Tamia H100 did 1,500 tips x 10 epochs
(split) in 9.7 s. This run is 3,000 tips x 100 epochs, split-only. Scaling
epochs 10x and tips ~2x (conservative; attention cost may scale worse than
linearly in tips) gave an estimate of roughly 5-15 minutes, well under the
45-minute stop threshold, so the run proceeded directly.

**Result**: exit code 0, total wall time (script start to exit) 103 s.
`split.rds$elapsed_s` (the fit_pigauto() training-loop timer pigauto itself
records) = **43.086 s**. `run.log` contained only one warning
(`phylo_signal_gate requires the 'phytools' package; returning NA for all
traits.`) — no errors.

**Device**: the script does not pass `verbose = TRUE` to `fit_pigauto()`, so
the "Using device: ..." message was not printed to `run.log`. Device
selection (`get_device()` in `R/utils_torch.R`) picks `cuda` whenever
`torch::cuda_is_available()` is `TRUE`, which was confirmed `TRUE` in this
exact environment (Section 2) immediately before and is process-global for
the R session that ran the script — so the run used **cuda**.

Outputs written: `mask_receipt.rds` (227 KB), `split.rds` (55 KB), `run.log`
(98 bytes) in
`~/pigauto_mondrian_realdata/results/fishbase3000-kohaku-prerun/`.
