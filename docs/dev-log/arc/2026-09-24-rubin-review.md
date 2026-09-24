# Rubin MI review (Meng lens): freq vs BACE lane

Reviewer: Meng (missing-data / multiple-imputation lens), 2026-09-24. Read-only review of worktree
`pigauto-rubin-freq-bace` at `46f3217`. Probes live in `/private/tmp/claude-503/meng/` (p1 to p8, p6b); none
touched repo code.

## Verdict

**GO WITH FIXES for the BACE settings pre-run** (three cheap fixes below, none of which change the grid).
**The campaign is NOT ready**: two BLOCKING findings (B1, B2) must be resolved before it launches.

Pre-run fixes, before `CONFIRM=yes`:

1. Guard the scoring stage so every cell writes an rds (B1's fix; for the BACE-only pre-run the crash is
   unlikely but possible, and the pre-run plan's "failures are counted, not lost" depends on it).
2. Use only the `bace` arm in selection rule (c); label `bace_resid` a negative control (N8).
3. Correct the stale "posterior means" header in `script/rubin_cell.R:13-14`, and the plan's sentence
   "Each final dataset is therefore a posterior predictive draw" (B2).

Labels: **Evidence** = code read or probe output; **Inference** = follows from evidence but not measured in
the study's own DGP; **Speculation** = plausible, untested.

---

## BLOCKING

### B1. A failed freq A draw crashes the whole cell; no rds is written (BLOCKING for the campaign)

- **Evidence.** `mi_freq_A()` leaves `datasets[[m]] <- NULL` when two refits fail (`script/rubin_freq.R:195`,
  `next`). Scoring runs outside `run_arm()`'s `tryCatch` (`script/rubin_cell.R:138`). Probe p1:
  `vapply(list(NULL, df), function(s) as.numeric(s[["c1"]][idx]), numeric(3))` errors "values must be length
  3, but FUN(X[[1]]) result is length 0". So `score_cells()` aborts the script. Separately, if
  `est_pgls_slope()` ever returns `lambda_hat = NA`, `est_phylo_cor()` (`script/rubin_lib.R:165-169`) calls
  `chol()` on an NA matrix and errors ("'a' must have dims > 0" for the NULL case, p1). That path is also
  unguarded and applies to BACE datasets.
- **Why it matters.** The rerun is deterministic (same seed), so the cell is lost every time. Lost cells are
  exactly the hard ones (refit failures), so failure becomes informative missingness in the results table.
  A reviewer would ask how many cells were silently dropped.
- **Inference.** Refit failure is rare: 0 of 160 bootstrap refits failed in p2. But the campaign runs about
  200 reps x 12 cells x 20 refits.
- **Fix.** In `rubin_cell.R`, drop `NULL` entries before scoring (`sets <- Filter(Negate(is.null), sets)`)
  and record `m_used`. Wrap `score_cells()` and `score_estimands()` per arm in `tryCatch`, storing the error in
  `errors[[arm]]`. In `est_phylo_cor()`, return NA fields when `lambda` is not finite. `rubin_pool()` already
  needs `length(q) >= 2`, so pool on `m_ok` and report it (the G-S4a gate already checks `m_ok == M`).

### B2. BACE as shipped is not proper MI, and `bace_resid` does not test the fix (BLOCKING for campaign design and for anything sent to Dan)

- **Evidence (installed BACE, deparsed).** `bace_final_imp()` runs
  `lapply(seq_len(n_final), function(run) .one_final_run(run, last_data, ...))`. **Every final run starts
  from the same `last_data`**, the last convergence iterate. It does not continue from the previous final
  dataset. The in-code comment reads "Each run is independent (all start from the converged last_data),
  giving truly independent posterior draws suitable for Rubin's rules pooling." `last_data` comes from
  `bace_imp()`, whose convergence runs call `.predict_bace(..., sample = FALSE)`. Those are deterministic
  posterior-mean fills. Inside a final run, traits are refit in formula order (c1, c2, cnt, prp, bin, ord,
  cat3, d1). Trait v's missing cells go back to NA and are augmented, but its predictors hold (i) this run's
  draws for traits earlier in the order and (ii) **`last_data`'s deterministic fills, identical in all M
  runs**, for traits later in the order.
- **Evidence (probe p5, a real BACE fit, n = 60, seed 2, smoke chains).** c1's MCMCglmm design matrix
  `all_models[[i]][["c1"]]$X` is **identical across all 20 final runs** (TRUE). c2's is not, because c1 was
  redrawn first. So c1's between-imputation variance never includes the uncertainty in its missing
  predictors.
- **Evidence (toy, probe p6; iid tips, 4 traits at rho 0.5 plus a fully observed driver at 0.35, n = 100,
  30% MCAR, M = 20, 1000 reps; each imputation draws parameters from the flat-prior posterior and adds a
  residual draw).** BACE's scheme, 20 one-sweep runs from a common deterministic state:
  slope coverage **0.900** (MCSE 0.010), per-cell coverage **0.882 (c1) / 0.903 (c2)**, and mean B 0.75x
  that of proper independent FCS chains. The proper chains cover 0.955 / 0.949 / 0.950.
  Probe p6b (600 reps) tested a **chained** variant, where final run i starts from final dataset i-1.
  It restores coverage: slope 0.950, cells 0.943 / 0.947 at thin 1, and 0.957 / 0.950 / 0.946 at thin 2.
- **Inference.** The same deflation should appear in the phylogenetic DGP. How large it is there is
  unmeasured: it depends on the fraction of co-missing predictors and on rho. The direction is
  undercoverage, largest for traits early in the formula order. More `runs` or `nitt` will not remove it,
  because the anchor stays shared.
- **Consequences.**
  - The plan's correction block says "Each final dataset is therefore a posterior predictive draw". That is
    only true *conditional on a shared anchor*. The draws are not independent across m.
  - `bace_resid` was meant to answer "do proper predictive draws fix BACE's coverage?". With the installed
    build it adds a second residual (N8) and leaves the anchor untouched. So the campaign has no arm that
    answers its own question.
  - Nothing should tell Dan that BACE's final datasets are "truly independent posterior draws". His own
    code comment is the claim a reviewer would test.
- **Fix (design call for Shinichi).**
  - Keep `bace` (as shipped) as the headline BACE arm.
  - Replace `bace_resid` with a `bace_chain` arm that re-uses BACE's own code without editing it. Take
    `outb$initial_results` (class "bace") and call `BACE:::bace_final_imp(..., n_final = 1)` M times (or
    2M with thin 2), each time replacing `bace_object$data[[length(.)]]` with the previous run's dataset.
    This is feasible by inference, not yet tested.
  - Cost is about +M sweeps per fit: roughly x1.8 of the shipped wall time at `runs = 5`.
  - Tell Dan the mechanism and the toy numbers.
  - If Shinichi wants the pre-run to size `bace_chain` too, add it there; otherwise the pre-run can go
    ahead for the shipped arm alone.

---

## NON-BLOCKING

### N1. Pooling arithmetic is correct (claim upheld)

- **Evidence.** The installed `mice:::barnard.rubin` has no lambda floor (p1 prints its body), and
  `.barnard_rubin_df()` is algebraically identical to it. B = 0 gives the same df (56.098 at dfcom 58; Inf
  at dfcom Inf).
- `fmi = (riv + 2/(df+3))/(riv+1)` is mice's definition. Note it is evaluated at the Barnard-Rubin df,
  not Rubin's old df, as mice does.
- Per-cell interval: with W = 0, Rubin's df `(M-1)/lambda^2` has lambda = 1, so df = M - 1. The
  t_{M-1} x sqrt((1+1/M)B) interval is the exact predictive interval for a new exchangeable normal draw.
  The 0.9508 check validates that. It does not validate that any arm's draws are exchangeable with the
  truth, which is what the campaign measures.

### N2. Fisher-z pooling uses a t reference while the complete-data reference is normal

- **Evidence.** `pool_cor()` passes `df_com = n - 3` into Barnard-Rubin (`script/rubin_lib.R:195`), so even
  with B = 0 the CI uses t_{about 55} at n = 60. The complete-data gate uses `qnorm`
  (`script/rubin_checks.R:45`). Probe p8, complete data: at n = 60, coverage with z is 0.935 and with
  t_{n-3} is 0.943. At n = 100 the gains are +0.0 to +0.8 pp.
- **Fix.** Use `df_com = Inf` for Fisher z (normal complete-data reference, mice's default), or use
  t_{n-3} for the complete-data reference too. Either way, be consistent.

### N3. The complete-data reference is not 0.95; score MI arms against complete data on the same replicates

- **Evidence (p8, 600 reps each, driver on, MCSE about 0.010).**

  | cell | slope | cor (z) | cor, true lambda | mean lambda_hat | lambda_hat at 0 or 1 |
  |---|---|---|---|---|---|
  | n 60, lambda 0.7 | 0.923 | 0.935 | 0.957 | 0.622 | 2% |
  | n 100, lambda 0.3 | 0.947 | 0.948 | 0.952 | 0.280 | 10% |
  | n 100, lambda 0.7 | 0.942 | 0.938 | 0.950 | 0.652 | 1% |
  | n 100, lambda 1.0 | 0.947 | 0.935 | 0.940 | 1.000 | 61% |

- **Evidence.** The empirical sd(z) is 5 to 7% above sqrt(1/(n-3)).
- **Inference.** 1/(n-3) is exact for the GLS-whitened r when lambda is known. The whitened pairs are iid
  bivariate normal and GLS demeaning removes one df, as for Pearson's r. The shortfall comes from plugging in
  lambda_hat. corPagel's REML lambda_hat is biased low (0.62 at n = 60 for a true 0.7). The G-S4b n = 60
  numbers (0.930 / 0.920) agree with this.
- **Fix.** Add a `complete` row to `est_tab` per replicate, computed from `ref_slope` / `ref_cor` with the
  same df convention. Report MI coverage relative to it. Make no coverage claims at n = 60. `df_com = n - 2`
  for the slope is standard. Taking one more df for lambda is defensible but immaterial at n >= 100.

### N4. `est_pgls_slope()`'s bounded fallback is correct; mixing paths does not bias

- **Evidence (p4, 60 complete-data reps per cell).**
  - The fallback fires in 8/60 (n 60) and 5/60 (n 100) reps at lambda 0.3, and in 43/60 and 39/60 at
    lambda 1.0.
  - In all 240 reps, every free fit that was accepted sits at the bounded-REML maximum: the logLik gap is
    at most 0.01, and the maximum gap is 0.000.
  - So the estimator is bounded REML throughout. Both paths report vcov conditional on lambda_hat.
- **Minor.** `est_phylo_cor()` refits `est_pgls_slope()` internally, which doubles the gls cost per
  dataset. Pass `lambda` from the slope fit (performance only).

### N5. The Armadillo "system is singular" messages are benign; do not count them as failures

- **Evidence (p2).** They print to stderr (not as R warnings) in 38 to 40 of 40 bootstrap refits in every
  cell tried: n 60 / lambda 0.7, n 100 / lambda 0.3 (two seeds), and n 100 / lambda 1.0. All refits
  returned PD Sigma_p.
- **Evidence (p3).** For 10 refits, Rphylopars' `model = "lambda"` lambda* equals, to 1e-4 in lambda and
  1e-3 in logLik, the maximiser of an independent REML profile (`model = "BM"` on the lambda-transformed
  tree). The internal solves hit singular trial matrices during the optimisation, and the reported optimum
  is still the global one.

### N6. Freq A does what it claims; it is approximately proper, with caveats to state

- **Evidence.** The code simulates Y* ~ N(mu0, V0) from pars0, applies the original mask (`mask_block`),
  refits to get theta* (mu*, Sigma*, lambda*), and draws the **original** data's missing cells from the
  conditional under theta* (`Yt_orig` from `df_miss`, `script/rubin_freq.R:197-198`). It does not draw Y*'s
  cells and does not reuse theta_hat. mu* is used as well, so the root uncertainty is carried.
- **Caveats.**
  - **(a) Inference.** A parametric bootstrap approximates the posterior of theta only asymptotically.
    Its distribution centres on theta_hat + bias: lambda* means were 0.03 below lambda_hat0 at lambda 0.3
    and 0.015 below at 0.7 (p2), about 0.2 to 0.3 bootstrap SDs. Call it "approximately proper" in the paper.
  - **(b) Evidence.** At a true lambda of 1, lambda* has sd 0.002 (range 0.992 to 0.999), so lambda
    uncertainty is effectively nil there. That is plausible for rcoal trees with short tips, and harmless.
  - **(c) Evidence.** The block model is misspecified for prp. qlogis(prp) = L4 + N(0, 0.3^2) iid, but
    Rphylopars fits one shared lambda and no phenotypic error (`pheno_error` FALSE). prp's effective lambda
    is about 0.7/1.09. The effect on c1 / c2 is expected to be small (inference).
  - **(d) Speculation.** Under MAR cells, re-using the fixed mask while d1* is drawn independently of R
    approximates Var(theta_hat | R) under MCAR-like sampling (compare Kenward and Molenberghs 1998 on
    expected vs observed information). This does not affect the MCAR core. Check freq A vs B in the MAR
    cells before claiming anything there.
  - **(e) Evidence (p7).** With the driver block {c1, c2, prp, d1}, the conditional mean matches
    `anc_recon` to 6e-15. `diag(cond_cov)` is always at or below Rphylopars' `anc_var` (median relative
    gap 0.1%, max 5e-3), consistent with `anc_var` including root-estimation variance. Negligible, and it
    only affects freq B.

### N7. Freq B differs from A only by parameter uncertainty (confirmed)

- **Evidence.** It is the same `fit_block()` on `df_miss`, with M draws from the same conditional at
  theta_hat.
- **Note.** Freq B's `fmi` column is not a fraction of missing information, because B omits the parameter
  term. Do not tabulate it as one.

### N8. `bace_resid` is a double residual with a mis-reconstructed scale; report it only as a negative control, or drop it

- **Evidence.** The installed gaussian branch already adds `rnorm(0, sqrt(units))` (the plan's correction
  block, confirmed).
- **Evidence (p5).** `.bace_gaussian_sd()` (`script/rubin_bace.R:105-113`) assumes run i reads run i-1's
  dataset. BACE actually starts every run from `last_data` (B2), and it computes `sd_val` via
  `.extract_gaussian_attrs(data_i)` with the response's missing rows set back to NA. That makes `sd_val`
  the **observed-only** sd, the same in all runs. Against c1's observed-only sd of 0.623, the
  reconstruction ranged 0.570 to 0.773 (mean ratio 1.023). Similar errors for c2 and prp.
- **Evidence.** `set.seed(seed)` (`script/rubin_bace.R:139`) re-uses the DGP seed and resets the global
  stream mid-cell.
- **Fix.** Replace it with `bace_chain` (B2). If it is kept, use `sd(df_miss[[v]], na.rm = TRUE)` and a
  distinct seed (e.g. `seed + 7919L`).
- **cnt is fine.** The poisson `sample = TRUE` branch draws `rpois(exp(Liab))` from a posterior iteration's
  latent, which already includes the units overdispersion. It is a proper predictive draw; adding a residual
  would be wrong.

### N9. Both arms see d1, but BACE sees more auxiliaries; this is an efficiency confound, not a validity defect

- **Evidence (p1).** `cont_traits` = c1, c2, cnt, prp, d1, so `default_block_traits()` = **c1, c2, prp, d1**.
  Freq includes d1. BACE regresses every trait on all the others (`script/rubin_bace.R:55-56`), so it also
  uses cnt, bin, ord and cat3, whose liabilities correlate 0.5 with c1 / c2 at rho = 0.5.
- **Inference.** Omitting auxiliaries keeps freq's (c1, c2) imputation valid under MCAR and under MAR on d1,
  because d1 is in the block. So coverage comparisons are fair. Width, RMSE and interval score favour BACE
  partly by information, not method.
- **Fix.** State this in the methods ("stack vs stack", v1's freq design). Compare coverage as the validity
  axis, and read width only among arms that are valid.

### N10. prp and cnt per-cell scoring on the raw scale

- **Evidence (p5).** BACE models prp as gaussian on the raw proportion scale; 6 of 360 smoke imputations
  exceeded 1 (maximum 1.073). This is a property of BACE as shipped, worth telling Dan.
- **Inference.** Symmetric t intervals on a bounded, skewed scale distort coverage for both arms.
- **Fix.** Score prp on the logit scale (clipping BACE's draws for the transform, and reporting how many
  were out of range). cnt is BACE-only, so keep it out of cross-arm tables.

### N11. RNG pairing across arms and settings

- **Evidence.** All arms share one L'Ecuyer stream in sequence (`script/rubin_cell.R:135-160`). BACE's draws
  therefore depend on whether freq A / B ran first. The pre-run (`--arms bace,bace_resid`) will not
  reproduce campaign BACE results at the same seed.
- **Fix.** Call `set.seed(seed + <arm offset>)` before each arm.

### N12. Monte Carlo precision for BACE

- **Inference.** 100 BACE reps per cell give a coverage MCSE of about 0.022. A B2-sized shortfall (0.90 vs
  0.95) is only about 2.3 MCSE. Use 200 reps, or pre-specify pooling over lambda, before claiming
  undercoverage per cell. Per-cell coverage pools correlated cells within a replicate, so compute its MCSE
  by replicate.

### N13. Stale or mismatched text

- **Evidence.**
  - `script/rubin_cell.R:13-14` says "bace (as shipped: continuous imputations are posterior means)". The
    plan's Arms table repeats that.
  - Plan estimand 3 says "Pearson correlation", but the code uses the GLS-whitened correlation, which is the
    right choice.
  - `script/rubin_bace.R:92-99` has the wrong run-chaining premise.
- **Fix.** Correct all three before anything is shared.

### N14. `score_cells()` silently skips a trait if any draw is NA

- **Evidence.** `script/rubin_cell.R:87`. This is intended for traits outside the freq block, but it would
  also hide a partially failed BACE imputation.
- **Fix.** Skip only when *all* draws are NA; otherwise record `n_na`.

---

## What was not checked

- Freq A's coverage in the phylogenetic DGP. No multi-replicate freq A run was made here; it costs about
  18 s per replicate at n = 60.
- `bace_chain`'s feasibility inside BACE's internals.
- BACE's convergence diagnostic semantics. It diagnoses a deterministic-fill chain, which converges to a
  fixed point rather than to a distribution. That is relevant to how Dan describes "convergence", and it
  was not reviewed.
