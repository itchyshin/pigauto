## Downstream slope (posterior_full vs complete data, same analysis model)

Paired bias = mean over reps of (MI pooled slope - complete-data slope in the same rep); MCSE in brackets.
SE ratio = mean pooled SE / empirical SD of the pooled slope; 'rel' = MI ratio / complete ratio (gated under phylolm, [0.90, 1.15]).
rel_mcse: approximate Monte Carlo SE of 'rel', treating the two empirical SDs as independent (each has relative SE about 1/sqrt(2(R-1))); reported only, not part of any gate.
Coverage truth: 0.7 in regimes 1-16; mean complete-data slope in 17-24.

| regime | analysis | paired_bias | se_ratio | complete_se_ratio | rel | rel_mcse | coverage | complete_coverage | plugin_se_ratio | plugin_coverage | converged |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 (lambda 1.0, n 300, MCAR, x_only) | gls | -0.013 (0.001) | 0.89 | 0.87 | 1.02 | 0.07 | 0.920 | 0.915 |  |  | 199/200 |
| 1 (lambda 1.0, n 300, MCAR, x_only) | phylolm | -0.026 (0.002) | 0.96 | 0.92 | 1.05 | 0.07 | 0.910 | 0.905 |  |  | 199/200 |
| 2 (lambda 0.5, n 300, MCAR, x_only) | gls | 0.000 (0.003) | 0.70 | 0.55 | 1.26 | 0.09 | 0.840 | 0.805 |  |  | 200/200 |
| 2 (lambda 0.5, n 300, MCAR, x_only) | phylolm | -0.014 (0.001) | 1.00 | 0.92 | 1.08 | 0.08 | 0.940 | 0.940 |  |  | 200/200 |
| 3 (lambda 1.0, n 1000, MCAR, x_only) | gls | -0.007 (0.001) | 0.91 | 0.89 | 1.02 | 0.07 | 0.899 | 0.915 |  |  | 198/200 |
| 3 (lambda 1.0, n 1000, MCAR, x_only) | phylolm | -0.019 (0.001) | 0.96 | 0.92 | 1.04 | 0.07 | 0.854 | 0.905 |  |  | 198/200 |
| 4 (lambda 0.5, n 1000, MCAR, x_only) | gls | 0.006 (0.002) | 0.71 | 0.55 | 1.29 | 0.09 | 0.865 | 0.765 |  |  | 200/200 |
| 4 (lambda 0.5, n 1000, MCAR, x_only) | phylolm | -0.010 (0.001) | 1.08 | 1.05 | 1.04 | 0.07 | 0.975 | 0.960 |  |  | 200/200 |
| 5 (lambda 1.0, n 300, MAR_phylo, x_only) | gls | -0.002 (0.002) | 0.94 | 0.89 | 1.06 | 0.08 | 0.918 | 0.930 |  |  | 195/200 |
| 5 (lambda 1.0, n 300, MAR_phylo, x_only) | phylolm | -0.012 (0.002) | 0.97 | 0.94 | 1.03 | 0.07 | 0.918 | 0.915 |  |  | 195/200 |
| 6 (lambda 0.5, n 300, MAR_phylo, x_only) | gls | 0.015 (0.003) | 0.64 | 0.56 | 1.16 | 0.08 | 0.770 | 0.750 |  |  | 200/200 |
| 6 (lambda 0.5, n 300, MAR_phylo, x_only) | phylolm | -0.003 (0.002) | 0.87 | 0.93 | 0.94 | 0.07 | 0.920 | 0.920 |  |  | 200/200 |
| 7 (lambda 1.0, n 1000, MAR_phylo, x_only) | gls | 0.002 (0.001) | 0.86 | 0.88 | 0.98 | 0.07 | 0.898 | 0.895 |  |  | 197/200 |
| 7 (lambda 1.0, n 1000, MAR_phylo, x_only) | phylolm | -0.009 (0.001) | 0.88 | 0.92 | 0.95 | 0.07 | 0.929 | 0.935 |  |  | 197/200 |
| 8 (lambda 0.5, n 1000, MAR_phylo, x_only) | gls | 0.015 (0.004) | 0.65 | 0.46 | 1.43 | 0.10 | 0.805 | 0.730 |  |  | 200/200 |
| 8 (lambda 0.5, n 1000, MAR_phylo, x_only) | phylolm | -0.004 (0.001) | 0.90 | 1.03 | 0.88 | 0.06 | 0.915 | 0.960 |  |  | 200/200 |
| 9 (lambda 1.0, n 300, MCAR, both) | gls | -0.013 (0.003) | 0.97 | 0.90 | 1.08 | 0.08 | 0.929 | 0.910 | 0.87 | 0.894 | 198/200 |
| 9 (lambda 1.0, n 300, MCAR, both) | phylolm | -0.021 (0.003) | 1.03 | 0.97 | 1.06 | 0.08 | 0.939 | 0.940 | 0.93 | 0.884 | 198/200 |
| 10 (lambda 0.5, n 300, MCAR, both) | gls | 0.011 (0.004) | 0.81 | 0.53 | 1.53 | 0.11 | 0.925 | 0.810 | 0.71 | 0.890 | 200/200 |
| 10 (lambda 0.5, n 300, MCAR, both) | phylolm | -0.009 (0.003) | 1.05 | 1.06 | 0.99 | 0.07 | 0.960 | 0.970 | 0.93 | 0.925 | 200/200 |
| 11 (lambda 1.0, n 1000, MCAR, both) | gls | -0.010 (0.001) | 1.03 | 0.91 | 1.14 | 0.08 | 0.950 | 0.930 | 0.91 | 0.925 | 199/200 |
| 11 (lambda 1.0, n 1000, MCAR, both) | phylolm | -0.016 (0.001) | 1.07 | 0.95 | 1.12 | 0.08 | 0.940 | 0.945 | 0.95 | 0.894 | 199/200 |
| 12 (lambda 0.5, n 1000, MCAR, both) | gls | 0.010 (0.003) | 0.88 | 0.58 | 1.54 | 0.11 | 0.900 | 0.740 | 0.79 | 0.890 | 200/200 |
| 12 (lambda 0.5, n 1000, MCAR, both) | phylolm | -0.006 (0.001) | 1.01 | 0.91 | 1.11 | 0.08 | 0.965 | 0.925 | 0.90 | 0.940 | 200/200 |
| 13 (lambda 1.0, n 300, MAR_phylo, both) | gls | -0.008 (0.003) | 0.98 | 1.01 | 0.97 | 0.07 | 0.954 | 0.960 | 0.87 | 0.919 | 197/200 |
| 13 (lambda 1.0, n 300, MAR_phylo, both) | phylolm | -0.013 (0.003) | 1.01 | 1.05 | 0.96 | 0.07 | 0.964 | 0.965 | 0.91 | 0.939 | 197/200 |
| 14 (lambda 0.5, n 300, MAR_phylo, both) | gls | 0.017 (0.004) | 0.92 | 0.58 | 1.58 | 0.11 | 0.920 | 0.725 | 0.81 | 0.870 | 200/200 |
| 14 (lambda 0.5, n 300, MAR_phylo, both) | phylolm | -0.000 (0.002) | 1.04 | 0.94 | 1.10 | 0.08 | 0.975 | 0.950 | 0.91 | 0.935 | 200/200 |
| 15 (lambda 1.0, n 1000, MAR_phylo, both) | gls | -0.002 (0.002) | 0.97 | 1.05 | 0.93 | 0.07 | 0.935 | 0.970 | 0.87 | 0.915 | 200/200 |
| 15 (lambda 1.0, n 1000, MAR_phylo, both) | phylolm | -0.008 (0.002) | 1.01 | 1.12 | 0.90 | 0.06 | 0.945 | 0.970 | 0.92 | 0.925 | 200/200 |
| 16 (lambda 0.5, n 1000, MAR_phylo, both) | gls | 0.016 (0.003) | 0.87 | 0.55 | 1.59 | 0.11 | 0.885 | 0.760 | 0.79 | 0.855 | 200/200 |
| 16 (lambda 0.5, n 1000, MAR_phylo, both) | phylolm | 0.001 (0.002) | 0.96 | 1.01 | 0.95 | 0.07 | 0.940 | 0.955 | 0.88 | 0.915 | 200/200 |
| 17 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MCAR, both) | gls | 0.006 (0.003) | 0.87 | 0.75 | 1.16 | 0.08 | 0.915 | 0.850 | 0.71 | 0.850 | 200/200 |
| 17 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MCAR, both) | phylolm | 0.000 (0.002) | 0.91 | 0.95 | 0.96 | 0.07 | 0.950 | 0.945 | 0.77 | 0.865 | 200/200 |
| 18 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MCAR, both) | gls | 0.002 (0.001) | 0.98 | 0.67 | 1.47 | 0.10 | 0.925 | 0.785 | 0.85 | 0.880 | 200/200 |
| 18 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MCAR, both) | phylolm | 0.001 (0.001) | 0.98 | 0.90 | 1.08 | 0.08 | 0.945 | 0.905 | 0.87 | 0.920 | 200/200 |
| 19 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MAR_phylo, both) | gls | 0.002 (0.002) | 0.93 | 0.74 | 1.26 | 0.09 | 0.940 | 0.860 | 0.78 | 0.885 | 200/200 |
| 19 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MAR_phylo, both) | phylolm | -0.001 (0.002) | 1.00 | 0.95 | 1.05 | 0.07 | 0.940 | 0.935 | 0.85 | 0.890 | 200/200 |
| 20 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MAR_phylo, both) | gls | 0.001 (0.001) | 0.88 | 0.64 | 1.37 | 0.10 | 0.925 | 0.810 | 0.75 | 0.875 | 200/200 |
| 20 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MAR_phylo, both) | phylolm | -0.001 (0.001) | 0.93 | 0.89 | 1.05 | 0.07 | 0.930 | 0.925 | 0.83 | 0.925 | 200/200 |
| 21 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MCAR, both) | gls | -0.009 (0.007) | 0.78 | 0.62 | 1.27 | 0.09 | 0.884 | 0.795 | 0.65 | 0.767 | 189/200 |
| 21 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MCAR, both) | phylolm | -0.015 (0.005) | 0.90 | 0.82 | 1.10 | 0.08 | 0.931 | 0.895 | 0.79 | 0.889 | 189/200 |
| 22 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MCAR, both) | gls | 0.004 (0.004) | 0.96 | 0.59 | 1.64 | 0.12 | 0.935 | 0.795 | 0.82 | 0.900 | 200/200 |
| 22 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MCAR, both) | phylolm | -0.004 (0.002) | 0.96 | 0.86 | 1.11 | 0.08 | 0.940 | 0.920 | 0.83 | 0.900 | 200/200 |
| 23 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MAR_phylo, both) | gls | -0.002 (0.008) | 0.82 | 0.69 | 1.19 | 0.09 | 0.875 | 0.820 | 0.65 | 0.799 | 184/200 |
| 23 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MAR_phylo, both) | phylolm | -0.010 (0.005) | 0.80 | 0.80 | 1.00 | 0.07 | 0.897 | 0.875 | 0.68 | 0.804 | 184/200 |
| 24 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MAR_phylo, both) | gls | -0.006 (0.004) | 0.90 | 0.56 | 1.59 | 0.11 | 0.925 | 0.760 | 0.76 | 0.880 | 200/200 |
| 24 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MAR_phylo, both) | phylolm | -0.004 (0.002) | 0.85 | 0.79 | 1.08 | 0.08 | 0.905 | 0.850 | 0.74 | 0.845 | 200/200 |

## Per-cell 95% predictive coverage (masked truth)

Coverage = covered cells / scored cells over all converged reps; width = mean interval width on the latent scale.
posterior_full vs conformal (impute(gnn = FALSE), same masked cells). MCAR is gated in [0.92, 0.98]; MAR_phylo (clade-biased) is reported only.

| regime | trait | mask | cells | posterior_cov | posterior_width | conformal_cov | conformal_width | width_ratio |
|---|---|---|---|---|---|---|---|---|
| 1 (lambda 1.0, n 300, MCAR, x_only) | x | MCAR | 17935 | 0.957 | 0.93 | 0.970 | 1.64 | 0.57 |
| 2 (lambda 0.5, n 300, MCAR, x_only) | x | MCAR | 18055 | 0.948 | 1.72 | 0.965 | 2.87 | 0.60 |
| 3 (lambda 1.0, n 1000, MCAR, x_only) | x | MCAR | 59539 | 0.949 | 0.81 | 0.969 | 1.33 | 0.60 |
| 4 (lambda 0.5, n 1000, MCAR, x_only) | x | MCAR | 60247 | 0.945 | 1.65 | 0.960 | 2.50 | 0.66 |
| 5 (lambda 1.0, n 300, MAR_phylo, x_only) | x | MAR_phylo | 17498 | 0.941 | 1.00 | 0.945 | 1.66 | 0.60 |
| 6 (lambda 0.5, n 300, MAR_phylo, x_only) | x | MAR_phylo | 17941 | 0.934 | 1.73 | 0.958 | 2.89 | 0.60 |
| 7 (lambda 1.0, n 1000, MAR_phylo, x_only) | x | MAR_phylo | 58984 | 0.939 | 0.89 | 0.937 | 1.30 | 0.68 |
| 8 (lambda 0.5, n 1000, MAR_phylo, x_only) | x | MAR_phylo | 60121 | 0.940 | 1.66 | 0.948 | 2.45 | 0.68 |
| 9 (lambda 1.0, n 300, MCAR, both) | x | MCAR | 17684 | 0.952 | 1.04 | 0.968 | 1.62 | 0.64 |
| 9 (lambda 1.0, n 300, MCAR, both) | y | MCAR | 17865 | 0.947 | 1.04 | 0.967 | 1.65 | 0.63 |
| 10 (lambda 0.5, n 300, MCAR, both) | x | MCAR | 17911 | 0.944 | 1.92 | 0.963 | 2.87 | 0.67 |
| 10 (lambda 0.5, n 300, MCAR, both) | y | MCAR | 18012 | 0.943 | 1.92 | 0.965 | 2.88 | 0.67 |
| 11 (lambda 1.0, n 1000, MCAR, both) | x | MCAR | 59719 | 0.951 | 0.92 | 0.966 | 1.32 | 0.69 |
| 11 (lambda 1.0, n 1000, MCAR, both) | y | MCAR | 59437 | 0.951 | 0.92 | 0.969 | 1.34 | 0.68 |
| 12 (lambda 0.5, n 1000, MCAR, both) | x | MCAR | 60235 | 0.945 | 1.85 | 0.959 | 2.52 | 0.73 |
| 12 (lambda 0.5, n 1000, MCAR, both) | y | MCAR | 60164 | 0.944 | 1.85 | 0.954 | 2.45 | 0.76 |
| 13 (lambda 1.0, n 300, MAR_phylo, both) | x | MAR_phylo | 17819 | 0.938 | 1.14 | 0.936 | 1.63 | 0.70 |
| 13 (lambda 1.0, n 300, MAR_phylo, both) | y | MAR_phylo | 17674 | 0.935 | 1.13 | 0.942 | 1.63 | 0.70 |
| 14 (lambda 0.5, n 300, MAR_phylo, both) | x | MAR_phylo | 17959 | 0.932 | 1.98 | 0.947 | 2.84 | 0.70 |
| 14 (lambda 0.5, n 300, MAR_phylo, both) | y | MAR_phylo | 17964 | 0.936 | 1.99 | 0.951 | 2.86 | 0.70 |
| 15 (lambda 1.0, n 1000, MAR_phylo, both) | x | MAR_phylo | 59980 | 0.938 | 1.00 | 0.936 | 1.31 | 0.77 |
| 15 (lambda 1.0, n 1000, MAR_phylo, both) | y | MAR_phylo | 59662 | 0.939 | 0.99 | 0.942 | 1.32 | 0.75 |
| 16 (lambda 0.5, n 1000, MAR_phylo, both) | x | MAR_phylo | 59827 | 0.936 | 1.90 | 0.940 | 2.39 | 0.79 |
| 16 (lambda 0.5, n 1000, MAR_phylo, both) | y | MAR_phylo | 60134 | 0.939 | 1.90 | 0.948 | 2.47 | 0.77 |
| 17 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MCAR, both) | x | MCAR | 17964 | 0.944 | 3.11 | 0.962 | 4.22 | 0.74 |
| 17 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MCAR, both) | y | MCAR | 17927 | 0.933 | 1.85 | 0.964 | 2.64 | 0.70 |
| 18 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MCAR, both) | x | MCAR | 60069 | 0.945 | 3.05 | 0.956 | 3.74 | 0.82 |
| 18 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MCAR, both) | y | MCAR | 59874 | 0.945 | 1.76 | 0.963 | 2.28 | 0.77 |
| 19 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MAR_phylo, both) | x | MAR_phylo | 18042 | 0.941 | 3.15 | 0.959 | 4.27 | 0.74 |
| 19 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 300, MAR_phylo, both) | y | MAR_phylo | 17930 | 0.944 | 2.04 | 0.949 | 2.64 | 0.77 |
| 20 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MAR_phylo, both) | x | MAR_phylo | 60007 | 0.947 | 3.13 | 0.955 | 3.78 | 0.83 |
| 20 (lambda 0.3/0.9 rP 0.7 rE 0.7, n 1000, MAR_phylo, both) | y | MAR_phylo | 60073 | 0.947 | 1.92 | 0.946 | 2.26 | 0.85 |
| 21 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MCAR, both) | x | MCAR | 16849 | 0.943 | 2.70 | 0.957 | 3.26 | 0.83 |
| 21 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MCAR, both) | y | MCAR | 16917 | 0.943 | 2.67 | 0.967 | 3.33 | 0.80 |
| 22 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MCAR, both) | x | MCAR | 59697 | 0.945 | 2.61 | 0.955 | 2.84 | 0.92 |
| 22 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MCAR, both) | y | MCAR | 60008 | 0.948 | 2.62 | 0.962 | 2.94 | 0.89 |
| 23 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MAR_phylo, both) | x | MAR_phylo | 16411 | 0.944 | 2.77 | 0.954 | 3.27 | 0.85 |
| 23 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 300, MAR_phylo, both) | y | MAR_phylo | 16656 | 0.944 | 2.76 | 0.958 | 3.32 | 0.83 |
| 24 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MAR_phylo, both) | x | MAR_phylo | 59635 | 0.947 | 2.68 | 0.954 | 2.93 | 0.92 |
| 24 (lambda 0.7/0.7 rP 0.7 rE 0.0, n 1000, MAR_phylo, both) | y | MAR_phylo | 60235 | 0.949 | 2.67 | 0.951 | 2.86 | 0.93 |

## Convergence (posterior_full)

Fits: 4800; converged: 4756 (99.08%). Median of max R-hat per regime: 1.001 to 1.003; median of min bulk ESS: 1,491 to 3,332.

Source: `docs/dev-log/mi-posterior/sim_summary.csv`, `docs/dev-log/mi-posterior/cell_coverage.csv`; code SHA 69670d44f96e5934d194abac8c6653ee82e11811.
