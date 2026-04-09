# pigauto scaling diagnostic (Part 3 baseline, tag v0.3.0)

**Run**:    2026-04-09 10:05:13 -- 10:20:17 MDT (wall 15m 04s)
**Machine**: macOS Darwin 25.4.0 (arm64), R 4.5.2, CPU-only torch
**Commit**: `8bd867d` (tag `v0.3.0`)
**Dataset**: `ape::rcoal(n)` simulated trees; 3 continuous (BM) + 1 binary +
            1 ordinal traits; 15% MCAR missingness; one seed per n.
**Stages**:  preprocess -> graph -> splits -> baseline -> train (100 epochs)
            -> predict (5 imputations). 30-minute per-stage budget, 3-hour
            overall budget. No stage timed out.

## Headline numbers

| N       | wall (s) | graph wall (s) | graph share | R heap peak (MB) |
|---------|----------|----------------|-------------|------------------|
|    300  |    18.6  |     0.02       |      0 %    |     62           |
|  1,000  |    23.0  |     0.42       |      2 %    |    118           |
|  2,000  |    28.4  |     2.18       |      8 %    |    345           |
|  3,000  |    46.2  |     5.96       |     13 %    |    727           |
|  5,000  |    95.0  |    26.05       |     27 %    |  1,566           |
|  7,500  |   166.9  |    85.78       |     51 %    |  3,474           |
| 10,000  |   509.8  |   400.04       |   **78 %**  |  5,385           |

**pigauto v0.3.0 already completes end-to-end at N = 10,000** in 510
seconds (~8.5 minutes) on a CPU-only Apple Silicon laptop. Nothing
errors, nothing runs out of memory, nothing exceeds its per-stage
budget. The scaling wall is a *speed* problem, not a *feasibility*
problem.

## Per-stage scaling (fitted log-log exponents)

| Stage      | Fitted slope     | What this means                                    |
|------------|------------------|----------------------------------------------------|
| graph      | **O(N^2.7)**     | Dense Laplacian `eigen()` -- matches theory        |
| splits     | O(N^1.9)         | Mask matrix allocation; irrelevant in wall time    |
| baseline   | O(N^0.8)         | Rphylopars post-order traversal; very efficient    |
| train      | O(N^0.5)         | 100-epoch loop dominated by constant overhead      |
| predict    | O(N^0.6)         | Trivial on any realistic N                         |
| preprocess | ~0s at every N   | (excluded from fit)                                |

**At N = 10,000 the graph stage alone is 400 seconds, 78% of the total
wall time.** The other five stages together take 110 seconds.

## Full per-stage results

| N | stage | status | wall (s) | R heap max (MB) |
|---|---|---|---|---|
| 300 | preprocess | ok | 0.0 | 18 |
| 300 | graph | ok | 0.02 | 42 |
| 300 | splits | ok | 0.0 | 19 |
| 300 | baseline | ok | 1.5 | 62 |
| 300 | train | ok | 16.9 | 48 |
| 300 | predict | ok | 0.2 | 44 |
| 1000 | preprocess | ok | 0.0 | 40 |
| 1000 | graph | ok | 0.4 | 116 |
| 1000 | splits | ok | 0.0 | 49 |
| 1000 | baseline | ok | 3.2 | 118 |
| 1000 | train | ok | 19.2 | 54 |
| 1000 | predict | ok | 0.2 | 53 |
| 2000 | preprocess | ok | 0.0 | 40 |
| 2000 | graph | ok | 2.2 | 345 |
| 2000 | splits | ok | 0.0 | 78 |
| 2000 | baseline | ok | 4.6 | 223 |
| 2000 | train | ok | 21.4 | 84 |
| 2000 | predict | ok | 0.3 | 80 |
| 3000 | preprocess | ok | 0.0 | 41 |
| 3000 | graph | ok | 6.0 | 727 |
| 3000 | splits | ok | 0.0 | 126 |
| 3000 | baseline | ok | 8.6 | 383 |
| 3000 | train | ok | 31.3 | 121 |
| 3000 | predict | ok | 0.4 | 120 |
| 5000 | preprocess | ok | 0.0 | 42 |
| 5000 | graph | ok | 26.1 | 1566 |
| 5000 | splits | ok | 0.0 | 278 |
| 5000 | baseline | ok | 11.3 | 479 |
| 5000 | train | ok | 57.1 | 253 |
| 5000 | predict | ok | 0.5 | 247 |
| 7500 | preprocess | ok | 0.0 | 43 |
| 7500 | graph | ok | 85.8 | 3474 |
| 7500 | splits | ok | 0.0 | 576 |
| 7500 | baseline | ok | 18.4 | 918 |
| 7500 | train | ok | 61.6 | 497 |
| 7500 | predict | ok | 1.1 | 492 |
| 10000 | preprocess | ok | 0.0 | 45 |
| 10000 | graph | ok | 400.0 | 5385 |
| 10000 | splits | ok | 0.1 | 992 |
| 10000 | baseline | ok | 24.1 | 1590 |
| 10000 | train | ok | 83.9 | 836 |
| 10000 | predict | ok | 1.6 | 834 |

## Notes and caveats from the run

- **Rphylopars emits `solve(): system is singular (rcond: ...)` warnings**
  at N >= 2000 and retries an approximate solve. The warnings are
  non-fatal and the baseline matrix comes back fine. They grow in
  severity with N (smallest rcond ~1e-17 at N=2000, ~1e-39 at N=7500),
  suggesting Rphylopars is pushing a numerical boundary at large N. Not
  a blocker today but worth keeping an eye on if a user reports
  unstable baseline estimates at very large trees.
- **R heap memory at the graph stage** peaks at 5.4 GB for N=10,000.
  This is a worst-case *reported* figure (max-used since last reset),
  not resident RSS. The dense cophenetic matrix at N=10k is 10000^2 *
  8 bytes = 800 MB on its own; the rest is LAPACK workspace and
  intermediate copies inside `eigen()`. Fine on a 64 GB laptop, but
  sparse Lanczos would drop this to roughly 1 GB total including the
  sparse Laplacian.
- **Training at 100 epochs is a stand-in, not converged.** The real
  user workflow uses epochs = 2000 with early stopping. Extrapolating
  from 100 -> 2000 epochs at N=10k gives ~30 minutes of training, which
  is fine but not instant. After Fix A the eigendecomposition is no
  longer the rate limiter and training wall becomes the biggest
  remaining single cost -- but it is still sub-linear in N, so the
  user can train overnight at N = 10,000 without difficulty.
- **Multi-observation path is NOT exercised by this benchmark.** The
  benchmark uses one row per species. CLAUDE.md already flags the
  multi-obs code path as not having an end-to-end benchmark; that is
  still true after this run.

## The wall of first failure

**N = 10,000, stage `graph`, wall 400 seconds (6.7 minutes).** No
stage errors or times out at any N in the ladder. The wall is a wall
of *speed*, not correctness.

## Verdict on the Part 4 fixes

### Fix A (sparse Lanczos eigensolver via RSpectra) -- PROCEED

**Impact**: reduces the graph stage from O(N^2.7) to O(N*k*iters),
expected 2-5 seconds at N=10k instead of 400. Total pigauto wall at
N=10,000 drops from ~510 s to ~115 s (~80% wall reduction).

**Risk**: Low. RSpectra is a small stable C++ package. The dense
fallback path remains for edge cases and for sanity checks in tests.

**Recommendation**: implement and ship as v0.3.1.

### Fix B (cache cophenetic distances across pipeline calls) -- PROCEED (secondary)

**Impact**: 4 cophenetic calls -> 1 across the pipeline. At N = 10,000
each cophenetic computation is ~15-25 seconds (it is inside the graph
stage plus a fourth call in fit_baseline). Caching saves ~60 seconds
of wall at N=10k and ~2.4 GB of allocator churn. After Fix A the graph
stage no longer hides the cophenetic cost, so the savings become
visible.

**Risk**: Negligible. Pure refactor.

**Recommendation**: ship together with Fix A as v0.3.1.

### Fix C (kNN-sparse adjacency / attention) -- DEFER

**Rationale**: training at N = 10,000 is already only 84 seconds for
100 epochs, and its scaling exponent is O(N^0.5). There is no GPU
bottleneck to unlock because we are on CPU-only torch here. Dense
attention uses ~400 MB per layer at N=10k which is fine. Fix C is a
fix for N = 50,000+ scale, not for N = 10,000.

**Recommendation**: leave as an opt-in future feature; add
`k_neighbors = NULL` argument only when a concrete N >= 50,000 use
case arrives.

### Fix D (fast Felsenstein BM baseline) -- ABANDON

**Rationale**: the baseline stage scales **O(N^0.8)** and takes only
24 seconds at N=10,000. Rphylopars is already efficient enough that a
hand-rolled Felsenstein implementation would save a handful of seconds
at best, at the cost of reimplementing a well-tested phylogenetic
method. The original plan listed Fix D as conditional on baseline
being a bottleneck; the benchmark proves it is not.

**Recommendation**: do not implement. Keep Rphylopars as the canonical
BM baseline. If the singular-solve warnings ever become a real problem,
address that separately via Rphylopars options or a regularised solve,
not by rewriting the baseline from scratch.

## One-line conclusion

> **pigauto v0.3.0 is already usable at N = 10,000 (8.5 min wall,
> 5.4 GB peak R heap). Fix A alone reduces wall time by ~80%, Fix B
> saves another ~12%, Fixes C and D are unnecessary at this scale.**

## Next steps (not this PR)

1. Implement Fix A (`RSpectra::eigs_sym` replacement at
   `R/build_phylo_graph.R:13`, dense fallback for N < 500 or when
   RSpectra is unavailable) on a feature branch off this baseline.
2. Implement Fix B (cophenetic caching) on the same branch.
3. Rerun this exact benchmark and compare wall times per stage.
   Expected: graph drops from 400 s to ~5 s at N=10,000; baseline
   unchanged; total wall drops from ~510 s to ~115 s.
4. Tag the result as `v0.3.1` and document in `NEWS.md`.
