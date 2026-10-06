# Optional fixed-effect adapter integration checks

These tests are excluded from package builds (`script/` is in `.Rbuildignore`).
They exercise the real model adapters without making `drmTMB` or `gllvmTMB`
installation requirements for pigauto. The small analysis-aware imputation
fixture runs `multi_impute_analysis()` → `with_imputations()` → backend fits →
`pool_mi()`; it does not run the phylogenetic posterior sampler.

Run the source checks in a fresh R process after installing the desired
optional packages. The tests skip only the adapter whose package is absent:

```sh
NOT_CRAN=true OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla script/cran-0.11-integration/check-adapters.R real
```

The `real` check requires both packages and fails if either is absent. The
same script accepts `drm-only`, `gllvm-only`, and `neither`. Each mode checks
the packages visible to that R process before fitting and checks the exact
number of backend skips. Run `check-adapters.R fixtures` in every
configuration to keep the shipped tests independent of both packages.

To test the installed package, select an isolated library and pass its path
to the script. The parent and saved-fit reload child both verify that they
loaded pigauto from that library. For example, in the `both` cell:

```sh
AUDIT_LIB=/private/tmp/pigauto-cran-011-matrix/both
R_LIBS_USER="$AUDIT_LIB" R_LIBS_SITE=/nonexistent \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla script/cran-0.11-integration/check-adapters.R \
    fixtures installed "$AUDIT_LIB"
R_LIBS_USER="$AUDIT_LIB" R_LIBS_SITE=/nonexistent \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla script/cran-0.11-integration/check-adapters.R \
    real installed "$AUDIT_LIB"
```

Run `fixtures installed <cell>` in every cell. For the real checks, use
`neither`, `drm-only`, `gllvm-only`, or `real` to match the cell. Omitting
`installed <cell>` selects the source-check mode. A mismatched library path
fails before testing.

Use four isolated package-library configurations to demonstrate install and
runtime behavior:

1. **Neither package:** install pigauto into a clean library with neither
   package available, then run the command above. The core package install
   must succeed; optional adapter tests skip.
2. **drmTMB only:** add `drmTMB` and its own dependencies to that library,
   leave `gllvmTMB` absent, and run the command. The `drmTMB` test must pass
   and the `gllvmTMB` test skips.
3. **gllvmTMB only:** add `gllvmTMB` and its own dependencies, leave `drmTMB`
   absent, and run the command. The `gllvmTMB` test must pass and the
   `drmTMB` test skips.
4. **Both:** install both and run the command. Both real adapter tests must
   pass.

Installation independence requires an isolated library and a fresh `Rscript`
process in the neither-package cell. The real check saves
each backend's `pigauto_mi_fits` object, starts a fresh child process without
the backend attached, reloads it, and calls `pool_mi()`. Deserialising a fit
may itself load its backend namespace; that is distinct from attaching it.
The `gllvmTMB_multi` adapter resolves its registered `tidy.gllvmTMB_multi()`
method by loading the owning namespace at runtime. Users need that package
for gllvmTMB fits; pigauto can be installed without it. The
`drmTMB` adapter uses the fit's registered `coef()` / `vcov()` methods.

On the 2026-10-06 audit worktree, the four installation configurations were
created under `/private/tmp/pigauto-cran-011-matrix/`. Each contained links to
the same installed dependency library, omitting both optional backends or
including the indicated backend packages. Both `R_LIBS` and `R_LIBS_USER`
pointed at the cell, so fresh R processes saw only that cell and the base R
library. `R CMD INSTALL -l <cell> .` succeeded in all four cells with pigauto
0.11.0, including when neither backend was visible. Installed-mode fixtures
passed 20 adapter and 46 provenance expectations in every cell. Installed-mode
real checks passed eight drmTMB expectations in its cell, eight gllvmTMB
expectations in its cell, and 16 when both packages were present; the neither
cell had the expected two skips. Source-mode checks passed separately. This
checks installation independence and adapter routing on a small analysis-aware
fixture. It is separate from the n=300,
M=20 posterior recovery and interval-coverage campaign.
