# BACE data snapshot

Per-dataset truth `.rda` files (traits + trees) snapshotted from
`daniel1noble/BACE`. Used by pigauto's CI wrappers
(`script/gha/run_bench_*.R`) so pigauto sees the **exact same data
BACE saw** in its head-to-head run.

## Provenance

- Source: `daniel1noble/BACE/dev/testing_data/data/*.rda` at the
  commit underlying Actions run **25329857467** (manual
  `workflow_dispatch` on 2026-05-04).
- Datasets: AVONET (9,993 species, 16 cols), PanTHERIA, AmphiBIO,
  BIEN, GlobTherm, LepTraits (each at full species count).
- Snapshot date: 2026-05-16.

## Files

| file | rows × cols | size |
|---|---|---|
| `avonet_traits.rda` + `avonet_tree.rda` | 9993 × 16 | 242 + 211 KB |
| `pantheria_traits.rda` + `pantheria_tree.rda` | ~5400 × 9 | 62 + 42 KB |
| `amphibio_traits.rda` + `amphibio_tree.rda` | ~6800 × 5 | 61 + 77 KB |
| `bien_traits.rda` + `bien_tree.rda` | ~7800 × 8 | 239 + 281 KB |
| `globtherm_traits.rda` + `globtherm_tree.rda` | ~800 × 6 | 41 + 31 KB |
| `leptraits_traits.rda` + `leptraits_tree.rda` | ~5000 × 7 | 87 + 94 KB |

## How the wrappers use this

`script/gha/_bace_compat.R` ports BACE's `benchmark_engine.R`:

1. `set.seed(2026)` + `sample(rownames(traits_df), 2000)` — pick the
   same 2000 species BACE picked.
2. Apply BACE's exact `.apply_mask()` — per-trait 30% MCAR using the
   same RNG state, so the masked cells match BACE's masked cells
   row-for-row.
3. Log-transform the per-dataset `LOG_TRAITS` from BACE's
   `0[0-7]_benchmark_*.R` (e.g. AVONET: `mass_g`, `wing_length_mm`,
   `beak_length_culmen_mm`, `tarsus_length_mm`).
4. Fit pigauto's `multi_impute()` on the same masked data.
5. Back-transform log preds to raw scale and evaluate against
   BACE's truth long-format.

The result is a true apples-to-apples head-to-head: same species,
same masked cells, same log conventions.

## Re-snapshotting

When Dan adds/updates a dataset, fetch the new .rda files from
his repo:

```bash
for f in <new_files>; do
  gh api "repos/daniel1noble/BACE/contents/dev/testing_data/data/${f}.rda" \
    --jq '.content' | base64 -d > useful/bace_data_snapshot/data/${f}.rda
done
```

Then update this README's "Snapshot date" and the file list.
