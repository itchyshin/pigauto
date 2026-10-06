# Fisher review: should "none,est" become the default? (2026-10-05)

**Verdict: NOT SUPPORTED as a default change today; SUPPORTED as an opt-in.** The gains are real in the simulated
regime. The costs are understated, and the regime the safety machinery was built for was never tested.

## What holds

- Continuous zRMSE is never worse (every d <= 0.000; best -0.094, SE 0.005).
- Discrete accuracy at lambda = 0.3: +0.043 to +0.085 in all 14 types_mixed cells, each > 6 MCSE.
- MCSEs are sound: dataset is the unit, traits pooled within it.

## What the summary understates

1. **AVONET Brier is not small.** Default beats the mode by 0.019 (0.516 vs 0.535); none,est by 0.002. About 90%
   of the probability skill on the only real dataset is lost, for +0.006 accuracy. The cause is discrete lambda,
   not the gate: floor,est is also +0.017, none,l1 +0.001.
2. **The Brier cost grows with n.** bm_mixed: -0.018 (n100), +0.005 (n300), +0.013 (n1000). At lambda = 1,
   n = 1000 the +0.011 to +0.016 is a 7 to 11% relative rise (0.142 to 0.158). Real users sit above n = 1000,
   untested. Most of it is the floor's removal (none,l1 alone +0.012).
3. **"Everywhere" means types_mixed only.** bace_dgp gains +0.018 to +0.030, SE 0.011 to 0.020 (not resolved).
4. **ou_mixed duplicates bm_mixed.** All numbers match to 3 decimals, including continuous zRMSE, which an OU
   continuous DGP should change. OU evidence is the tm OU cells only. Check the job arguments.
5. **Coverage** has no MCSE; two cells are slightly lower (tm OU clade0.3 l0.3 n100 -0.003).

## Statistical validity

The one-sided 2-MCSE rule over about 228 comparisons expects about 5 false flags, but the 12 none,est flags form
one coherent pattern (lambda = 1, n = 1000, Brier), so they are real. The bigger problem is power: "not flagged"
is read as "not worse". tm BM clade0.3 l1 n1000 (k = 44) has Brier +0.005 (SE 0.007), an interval that covers
the losses seen elsewhere. The k = 21 and 44 cells (others 60) need explaining. bace_dgp (k = 20) cannot detect a 0.02 harm. A non-inferiority margin (say Brier +0.01)
is the right screen.

## Confounds and the missing regime

All simulated discrete traits are thresholded Gaussian liabilities, the model family of BACE and of pigauto's
threshold baseline. The floor and gate exist for real weak-signal data where the baseline was 15 to 101% worse
than the mean. Here no arm is ever worse than the mean (default zRMSE < 1 everywhere), so the floor never had a
job. Its value was assumed away, not tested. "Estimated lambda now goes towards 0" is plausible but unmeasured.
Also untested: the LP path (all-discrete data; gate off there gives none,l1, which put cat3 below the mode in
screen 3), multi-obs, n > 1000.

## Screen 3 attribution

Mechanically right: with the floor off, a fired gate sets r_cal_bm = 1 (`fit_pigauto.R` ~661), so "gate on, floor
off" equals "both off". But "the gate is the larger cost" depends on order: the gate acts only through the floor's
mean corner, so "floor off alone gives all of none,est" is equally true. The defensible reading: the gate misfires
on discrete traits with true lambda = 0.3 (default 0.476 vs mode 0.469), likely because Pagel's lambda on
thresholded data is biased low.

## Options for the owner

- **A. Adopt none,est now.** Gains as above. Cost: AVONET probability skill mostly gone; large-n Brier worse and
  trending worse; weak-signal real data unprotected.
- **B. Narrow change: gate off for discrete traits (or compute it on the liability scale), keep the floor.**
  Gains +0.029 to +0.066 at l0.3; n = 1000 Brier cost only +0.003 to +0.008. AVONET +0.017 remains until its per-trait cause is found.
- **C. Opt-in now; three cheap tests (under 3 h on Totoro) before any flip:** (i) types_mixed at lambda = 0 and
  0.1, n = 300 and 1000, 60 seeds; (ii) the BIEN plants data the floor was built for, 10 to 20 masks;
  (iii) AVONET per-trait Brier plus one n = 3000 high-signal cell. Flip only if non-inferior on all three.

Recommended: C, then A or B by the result.
