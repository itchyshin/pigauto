# Review: posterior multiple imputation design (mi-posterior)

Reviewers: Gauss (Bayesian/statistical correctness), Rose (scope and claims).
Reviewed: `docs/dev-log/mi-posterior/design.md` (2026-09-24 version, 100 lines),
`.unlazy/mi-posterior/GATES.md`, `R/henderson_s_inv.R`, `R/draws_conditional.R`.

## BLOCKING (must change before building)

### B1. The precision used for `a` is underspecified and, read literally, gives the wrong scale

Design.md section 1 states the model's intent clearly: "A is the phylogenetic covariance
over nodes built so that its tip block equals pigauto's correlation matrix
R = cov2cor(vcv(tree))" and says to get there via "Q = A^{-1} over the extended tree with
the tip_sqrt_d correlation scaling, as draw_conditional_bm() already does." But section 2
step 1 then writes the precision of `vec(a)` as `Sigma_P^{-1} %x% Q + Sigma_E^{-1} %x% (Z'Z)`
using bare `Q`, and step 3's IW update uses bare `a' Q a`. `Q` as returned by
`build_henderson_S_inv()` (R/henderson_s_inv.R) is the RAW extended-tree precision: its own
header says plainly "For ultrametric trees diag(A) is constant and the two scales coincide up
to a global factor, but for coalescent/variable-depth trees the difference is substantial."
Concretely, the Schur complement of the raw Q over internal nodes gives `A_tip^{-1}`, i.e.
`Q_TT - Q_TI Q_II^{-1} Q_IT = A_tip^{-1}` where `A_tip = vcv(tree)` (non-unit diagonal for
non-ultrametric trees, and even for an ultrametric tree it is `d * R` for tip depth `d`, not
`R` itself unless the tree height happens to equal 1).

If the sampler literally draws `a` from `N(0, Sigma_P (x) Q^{-1})` using raw Q, then `a`'s
tip-marginal covariance is `Sigma_P (x) A_tip`, not `Sigma_P (x) R` as claimed in the
"Implied marginal" two lines later (section 1, line 23: `vec(Y) ~ N(mu, Sigma_P %x% R +
Sigma_E %x% I_n)`). These two statements are inconsistent as written. Concretely this means:
even on an ultrametric tree, the MCMC-sampled "Sigma_P" would really be `d * Sigma_P_true`
for tip depth `d`, so `lambda_k = Sigma_P[k,k] / (Sigma_P[k,k] + Sigma_E[k,k])` (design.md
line 24, and question 2 in the review brief) would be wrong unless the tree height is
normalised to 1. On a non-ultrametric tree the bias is tip-dependent and cannot be fixed by
a global rescale at all.

`draw_conditional_bm()` (R/draws_conditional.R) avoids this problem, but only because it
never needs `a` as a persistent state: it rescales the DATA vector going in
(`y[o] * henderson$tip_sqrt_d[o]`) and rescales the OUTPUT coming out
(`mu_tip / sd_r`), while using raw Q throughout for one single linear solve. That trick is
correct for a one-shot "solve for the unknowns given the knowns" computation, but the
mi-posterior sampler needs `a` itself as an explicit, persistent, correctly-scaled state
vector across sweeps, because step 3 needs the sufficient statistic `a' Q a` directly, and
step 2 needs `e = y - mu - a` at tips in the same units as `y`. "As draw_conditional_bm()
already does" is not a sufficient specification for that.

**The fix exists and is simple.** Let `D = diag(tip_sqrt_d at tip rows/cols, 1 at internal
rows/cols)` (N x N diagonal, so `D Q D` keeps Q's sparsity pattern exactly). Then
`Schur(D Q D, eliminate internal) = D_T (Q_TT - Q_TI Q_II^{-1} Q_IT) D_T = D_T A_tip^{-1} D_T
= R^{-1}` exactly (using `R^{-1} = D_T A_tip^{-1} D_T`, the standard covariance-to-correlation
identity, which is also the identity already used in `henderson_R_inv_apply`'s `cor_scale`
path). So `Q' = D Q D`, used consistently as the one precision matrix in both the `(a, mu)`
draw and the `a' Q' a` sufficient statistic, gives `a` a tip-block prior of exactly `R` and
leaves internal nodes on an arbitrary (irrelevant, since they never touch data) scale. This
also removes any need for separate rescale-in/rescale-out bookkeeping at each Gibbs step,
because `a` at tips is then directly in the same units as `y`.

Required change: rewrite section 2 steps 1 and 3 to say explicitly `Q'` (or `Q` rescaled by
`tip_sqrt_d` at tip rows/cols, identity at internal rows/cols), not bare `Q`, and add a
sentence stating this is a single fixed sparse matrix built once per tree (not per sweep).
Add a unit test (see R2 below) asserting the resulting tip-block equals
`cov2cor(ape::vcv(tree))` to numerical precision on both an ultrametric and a genuinely
non-ultrametric test tree.

### B2. The parameter-expansion (PX) scheme is asserted, not specified

Design.md section 2 (after step 4) says: "Parameter expansion for Sigma_P (Gelman 2006;
MCMCglmm's alpha.mu/alpha.V): a = diag(alpha) a*, with the IW on the working covariance...
Record the implementation choice in the code header." But step 3's formula, two paragraphs
above, is written as a plain IW update on the REAL `a` (`IW(nu_P + N, S_P + a' Q a)`), with no
indication of where PX enters. For question 3 in the review brief: is `a = diag(alpha) a*`
the right MCMCglmm-style scheme for a full K x K Sigma_P? Yes, this generalises correctly:
MCMCglmm expands a variance-component random effect with a per-trait diagonal alpha, updates
an IW on the WORKING covariance Sigma_P* using `a*' Q' a*` (not the real `a`), draws alpha
(commonly Gaussian, e.g. N(alpha.mu, alpha.V), by Gibbs or Metropolis depending on the
conditional), and recovers the real covariance as `Sigma_P = diag(alpha) Sigma_P*
diag(alpha)`. That is a real, well-established fix for the concentration problem: an IW prior
with minimal df (`nu = K + 1`, section 2's own choice) is genuinely weak on correlations but
still meaningfully informative on the marginal variances, more so as K grows or as the true
variance approaches the boundary (Sigma_P near singular, i.e. lambda near 0), which is exactly
the Gelman 2006 critique PX is meant to fix.

But design.md does not specify: (i) the prior on alpha (alpha.mu, alpha.V), (ii) whether alpha
is drawn jointly with `a*` or in a separate Gibbs/MH step, (iii) whether the `a`-drawing step
(section 2 step 1) draws the WORKING `a*` directly (with precision `Sigma_P*^{-1} (x) Q' +
Sigma_E^{-1} (x) (Z'Z)`, noting alpha multiplies into the likelihood term too since `a =
diag(alpha) a*` appears in `y - mu - Z a`) or draws real `a` and back-solves `a*`. This is not
a cosmetic gap: if an implementer instead uses `a' Q' a` (real a) in the IW update while ALSO
sampling alpha, the working covariance gets double-counted and the resulting prior on Sigma_P
is not what Gelman/MCMCglmm intended (nor is it clear it stays proper). Pin down the exact
equations (or point at the specific MCMCglmm formula being ported) before S1/S2 implement it,
and record it in design.md, not only "in the code header" after the fact, since G1/G2 need to
test against a written specification.

### B3. G4's convergence check covers 2 regimes; G6 relies on all of them being converged

GATES.md G4 checks split R-hat < 1.05 and bulk ESS > 400 "on regime 1 rep 1 and regime 23 rep
1" only. G6 ("simulation acceptance... every regime x analysis model") and G7 (per-cell
coverage per regime) both make claims across the FULL regime battery, with no per-regime
convergence gate of their own. If a regime other than 1 or 23 fails to converge (a likely
failure mode near the boundary regimes flagged in R3 below, or in especially small-n / large-K
regimes), G6/G7 could report misleading bias or coverage numbers attributable to
non-convergence rather than to the model, and nothing in the gate ledger would catch it. This
directly bears on question 5 ("are the gates able to fail honestly"): as written, no.

Required change: either (a) have `script/mi_gls/04_acceptance.R` and `05_cell_coverage.R`
compute and report (not necessarily gate on) R-hat/ESS for every regime they touch, with a
loud flag on any regime that fails G4's thresholds, or (b) extend G4 itself to run on all
regimes (or a representative superset larger than 2) before G6/G7 are allowed to run. Given
compute cost is a real constraint here, (a) is the cheaper fix and is sufficient as long as a
non-converged regime cannot silently pass G6/G7.

## REQUIRED (change during build)

**R1. G1's floor ("MI_POST_TESTS_OK" at just `passed >= 10`) is too weak for this sampler's
complexity.** A Gibbs sampler with a sparse joint (a, mu) draw, a PX scheme, and two IW
updates needs, at minimum: the Q' identity test from B1; a direct check that `a' Q' a` matches
a dense reference computation on a small tree; a check that the flat prior on mu gives a
proper (non-degenerate) full conditional; and the PX recovery formula `Sigma_P = diag(alpha)
Sigma_P* diag(alpha)` tested against a hand-computed small case. Ten tests total is not a
meaningful bar; ask S1/S2 to enumerate what G1 must cover, not just count assertions.

**R2. Clarify what G2 actually tests.** "with parameters fixed, draws match the exact dense
Gaussian conditional (n=40, K=2); Sigma_E = 0 reproduces the prototype conditional mean" is
read most naturally as two separate checks: (i) the general case (Sigma_P and Sigma_E both
positive, parameters fixed) matches the dense conditional of the FULL marginal `N(mu, Sigma_P
%x% R + Sigma_E %x% I_n)` for `y_mis | y_obs`, not just the `a`-posterior; and (ii) the
Sigma_E = 0 special case matches `draw_conditional_bm()`. Question 5 in the brief asks
specifically whether there is a check that per-cell draws are posterior predictive (i.e.
include Sigma_E, not just `a`). Design.md section 2 step 2 correctly specifies that per-cell
draws come from `y_mis | a, mu, y_obs, Sigma_E` (so Sigma_E's contribution is in scope by
design), but G2's one-line description does not make explicit that its dense reference
includes the Sigma_E term with Sigma_E > 0, only that Sigma_E = 0 is checked against the
prototype. Write G2 (or a G2b) so it unambiguously exercises Sigma_E > 0 against the true
dense conditional, since this is exactly the quantity feeding `mi$posterior$cell_interval`.

**R3. Confirm the regime table stresses the boundary, and that it is stressed by G3 as
written.** AVONET's own lambda is reported (design.md context, question 4) at 0.99; PX exists
specifically to keep mixing and point estimates sane in that regime. G3 as written checks
recovery at n = 1000 with 3 seeds, "within 0.1 of truth," with no regime named in GATES.md
itself (the regime table lives elsewhere, not shown to this review). Before relying on G3 as
evidence the boundary is handled, confirm: (a) at least one regime targets lambda >= 0.95 and
one targets lambda <= 0.05 for every trait; (b) a 0.1 absolute tolerance is tight enough to
catch prior-induced shrinkage at the boundary, since a naive (non-PX) IW(K+1, 0.01*I) posterior
mean for Sigma_E near a true value of 0 will not collapse to exactly 0, and the resulting
lambda bias could plausibly be small enough (a few hundredths) to pass a 0.1 gate while still
being a real, reportable prior-dominance effect worth describing in results.md even if it
passes the gate.

**R4. Confirm #187's lambda is the same quantity being compared against in G3.** Section 1
line 26 and G3's "within 2 posterior SD of #187's REML lambda" both assume #187's lambda is
computed on the same tip-level `R = cov2cor(vcv(tree))`, with no extra df correction and no
different variance normalisation. The algebra checks out under that assumption (marginal
correlation of `Sigma_P[k,k] R + Sigma_E[k,k] I` divided by its own diagonal is exactly
`lambda_k R + (1 - lambda_k) I`, Pagel's transform, for `lambda_k = Sigma_P[k,k] /
(Sigma_P[k,k] + Sigma_E[k,k])`), but this review did not have #187's code in scope to verify
the assumption itself. Have S1/S2 confirm it directly rather than inferring it from the name
"lambda."

**R5. Clade-biased mask coverage is named as an estimand (design.md section 5) but not gated
anywhere in GATES.md.** Either add an acceptance threshold for it (it is a legitimate test of
whether the tree-structured model is actually robust to phylogenetically clustered
missingness, which is a real selling point if true) or state explicitly in design.md that it
is descriptive only and is expected to run somewhat worse than the exchangeable-mask G7
threshold, so results.md does not overstate it later.

**R6. Test the reduction identities that section 1 asserts.** Line 26-27 claims the model
reduces to #187's single-lambda model when `Sigma_P = lambda * Sigma, Sigma_E = (1 - lambda) *
Sigma`, and to the `draw_conditional_bm()` prototype when `Sigma_E = 0`. The second is in scope
for G2. The first is not obviously tested anywhere in the current gate list; add it or note
it as untested in design.md so the claim is not left as an unverified assertion.

## SUGGESTIONS

- S1. Have `04_acceptance.R` / `05_cell_coverage.R` log per-regime R-hat/ESS in
  `docs/dev-log/mi-posterior/sim_summary.csv` even where not gated, for transparency per the
  after-task discipline (every number in results.md should trace to a file).
- S2. For G3's "phylogenetic correlation within 0.1," confirm this checks every off-diagonal
  entry of `cov2cor(Sigma_P)` when K > 2 in a regime, not a single summary scalar.
- S3. Given the tip/internal index-mapping complexity in the extended tree (B1), add a
  dimension/index assertion test that internal-node components of `a` are never accidentally
  used in place of tip values in the decode path (`build_completed`, `.dcb_decode_latent`-style
  helpers).

## VERDICT: REVISE

B1 and B2 are foundational: as literally specified, the sampler's core precision matrix for
`a` does not obviously reproduce the model design.md itself claims (the R-tip-block marginal),
and the PX scheme has no pinned-down equations for S1/S2 to implement against. B3 is a gating
gap that lets the main simulation claims (G6/G7) run without a convergence backstop across most
regimes. All three should be resolved in design.md (not left to implementation-time judgment)
before S1/S2 start building the sampler.
