# GNN earnings v1 -- preliminary (BM-misspecified phylolm)

WARNING: v1 uses phylolm(model=BM) which assumes zero residual
noise (misspecified vs the DGP iid noise term).  phylolm is
actually HANDICAPPED here, yet it still crushes pigauto on
linear and nonlinear effects.  v2 with model=lambda will be
even harder on pigauto.

## Mean RMSE per (beta, f_type, method): 3 reps each


## tree300 (n=300)

```
 beta      f_type rmse.column_mean rmse.lm rmse.phylolm_blup
  0.2 interactive            0.981   0.989             1.007
  0.4 interactive            1.058   1.073             1.067
  0.6 interactive            0.886   0.892             0.910
  0.2      linear            0.973   0.913             0.803
  0.4      linear            1.026   0.797             0.625
  0.6      linear            1.026   0.633             0.372
  0.2   nonlinear            0.992   1.002             0.984
  0.4   nonlinear            1.021   0.980             0.963
  0.6   nonlinear            0.994   0.927             0.941
 rmse.pigauto_no_cov rmse.pigauto_cov_sfT ratio  better
               0.960                0.960 0.953 pigauto
               1.049                1.010 0.947 pigauto
               0.876                0.856 0.940 pigauto
               0.925                0.908 1.130 phylolm
               0.984                0.978 1.565 phylolm
               1.042                0.968 2.598 phylolm
               0.992                1.039 1.056 phylolm
               1.059                1.074 1.115 phylolm
               1.052                1.000 1.063 phylolm
```

## mammal_n1000 (n=1000)

```
 beta      f_type rmse.column_mean rmse.lm rmse.phylolm_blup
  0.2 interactive            1.025   1.033             0.926
  0.4 interactive            1.017   1.019             0.975
  0.6 interactive            0.975   0.986             0.969
  0.2      linear            1.037   0.950             0.741
  0.4      linear            1.006   0.792             0.543
  0.6      linear            0.982   0.639             0.144
  0.2   nonlinear            1.005   0.991             0.908
  0.4   nonlinear            0.996   0.968             0.879
  0.6   nonlinear            0.985   0.945             0.850
 rmse.pigauto_no_cov rmse.pigauto_cov_sfT ratio  better
               1.104                1.144 1.236 phylolm
               1.074                1.002 1.028 phylolm
               1.071                0.795 0.820 pigauto
               1.069                0.934 1.261 phylolm
               1.103                0.835 1.538 phylolm
               1.171                0.705 4.906 phylolm
               0.918                0.883 0.972 pigauto
               1.157                1.167 1.327 phylolm
               1.050                1.108 1.304 phylolm
```

## amphibio_n1500 (n=1500)

```
 beta      f_type rmse.column_mean rmse.lm rmse.phylolm_blup
  0.2 interactive            0.989   0.989             0.954
  0.4 interactive            0.984   0.986             0.944
  0.6 interactive            0.977   0.978             0.944
  0.2      linear            0.989   0.870             0.752
  0.4      linear            0.992   0.779             0.549
  0.6      linear            1.009   0.630             0.048
  0.2   nonlinear            1.003   0.992             0.960
  0.4   nonlinear            0.989   0.969             0.898
  0.6   nonlinear            0.995   0.951             0.876
 rmse.pigauto_no_cov rmse.pigauto_cov_sfT  ratio  better
               0.916                0.926  0.970 pigauto
               0.973                0.975  1.033 phylolm
               1.017                0.819  0.867 pigauto
               0.904                0.871  1.159 phylolm
               0.971                0.850  1.548 phylolm
               0.974                0.753 15.689 phylolm
               0.931                0.912  0.950 pigauto
               0.902                0.883  0.983 pigauto
               0.914                0.869  0.992 pigauto
```

## Observations from v1 (with BM-misspecified phylolm)

1. **Linear effects: pigauto LOSES badly across the board.**
   At beta=0.6, pigauto is 2.6x (tree300) to 4.9x (mammal n=1000) WORSE
   than phylolm.  The linear regime is exactly where phylolm should
   shine -- and it does.  pigauto cannot catch this purely-linear signal.

2. **Nonlinear effects: pigauto loses or ties.**
   tree300: pigauto 6-12% worse than phylolm at all beta levels.
   mammal: pigauto 30% worse at beta=0.4-0.6.
   The GNN is NOT extracting the sin*exp signal better than a linear
   model with BM-BLUP correction.

3. **Interactive effects: pigauto wins MARGINALLY.**
   tree300: 5-6% better at beta=0.4-0.6.
   mammal: 18% better at beta=0.6 (the deterministic cell).
   This is the ONLY regime where the GNN earns its keep.

## Why this is honest

v1 phylolm uses model=BM, which assumes zero residual noise.  Our DGP
has eps ~ N(0, 1-alpha-beta) -- so phylolm is misspecified.  Despite
this handicap, phylolm beats pigauto on most cells.  v2 (with the right
model=lambda) will be EVEN HARDER on pigauto.

## Implications for the paper

The honest claim from v1 is **narrower than I had hoped**:
- pigauto is NOT a generally better imputer than phylolm-BLUP on
  semi-synthetic BM-on-tree data with linear or nonlinear covariate
  effects.
- pigauto wins specifically on INTERACTIVE covariate effects, by 5-18%.
- The paper claim should be scoped accordingly OR the experiment
  should use a DGP that is not stacked in phylolm BLUP-favor (e.g.,
  non-BM phylogenetic process, multi-trait, mixed types, multi-obs).

