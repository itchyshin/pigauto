# TabPFN Reconciliation Benchmark

- Generated: 2026-06-11 07:30:36 MDT
- Git commit: `e0329d8503f2cbc8a79997a8182d8914f6db8dba`
- Scales: `300`
- Regimes: `same_row, same_row_lappe`
- Replicates: `5`
- Split mode: `cell`
- Dry run: `FALSE`
- Pigauto configs: `default, relaxed_gate`

## Aim

This benchmark tests whether Russell-style TabPFN gains are mostly explained by same-row cross-trait features, phylogenetic features, or their combination. The shuffled regime keeps the same feature distribution but breaks row-level cross-trait alignment.

## Status Counts

                method status Freq
              baseline     ok   20
               pigauto     ok   20
  pigauto_relaxed_gate     ok   20
       tabpfn_same_row     ok   20
 tabpfn_same_row_lappe     ok   20

## Mean RMSE By Method

                method scale_n   rmse pearson_r coverage_95 wall_sec
 tabpfn_same_row_lappe     300 0.3428    0.9361      0.9604    10.54
       tabpfn_same_row     300 0.3658    0.9254      0.9702    10.58
  pigauto_relaxed_gate     300 0.5126    0.8573      0.9900   179.23
               pigauto     300 0.5139    0.8573      0.9900   180.51
              baseline     300 0.5664    0.8303          NA     0.34

## Mean RMSE By Trait

                method scale_n              trait   rmse pearson_r coverage_95
 tabpfn_same_row_lappe     300 Beak.Length_Culmen 0.4692    0.8786      0.9840
       tabpfn_same_row     300 Beak.Length_Culmen 0.5364    0.8454      0.9629
  pigauto_relaxed_gate     300 Beak.Length_Culmen 0.6709    0.7422      0.9800
               pigauto     300 Beak.Length_Culmen 0.6780    0.7361      0.9800
              baseline     300 Beak.Length_Culmen 0.6784    0.7381          NA
 tabpfn_same_row_lappe     300               Mass 0.2242    0.9719      0.9539
       tabpfn_same_row     300               Mass 0.2247    0.9688      0.9742
  pigauto_relaxed_gate     300               Mass 0.3183    0.9455      1.0000
               pigauto     300               Mass 0.3205    0.9486      1.0000
              baseline     300               Mass 0.4333    0.8990          NA
 tabpfn_same_row_lappe     300      Tarsus.Length 0.3920    0.9296      0.9228
       tabpfn_same_row     300      Tarsus.Length 0.4166    0.9224      0.9482
  pigauto_relaxed_gate     300      Tarsus.Length 0.5543    0.8604      0.9895
               pigauto     300      Tarsus.Length 0.5571    0.8594      0.9895
              baseline     300      Tarsus.Length 0.5955    0.8290          NA
       tabpfn_same_row     300        Wing.Length 0.2853    0.9652      0.9956
 tabpfn_same_row_lappe     300        Wing.Length 0.2859    0.9643      0.9807
               pigauto     300        Wing.Length 0.5001    0.8852      0.9906
  pigauto_relaxed_gate     300        Wing.Length 0.5069    0.8810      0.9906
              baseline     300        Wing.Length 0.5586    0.8550          NA

## Pigauto Gate Audit

               method scale_n              trait   rmse r_cal_bm r_cal_gnn
              pigauto     300 Beak.Length_Culmen 0.6780   0.7000    0.3000
 pigauto_relaxed_gate     300 Beak.Length_Culmen 0.6709   0.6800    0.3200
              pigauto     300               Mass 0.3205   0.3932    0.5963
 pigauto_relaxed_gate     300               Mass 0.3183   0.3337    0.6663
              pigauto     300      Tarsus.Length 0.5571   0.6429    0.3571
 pigauto_relaxed_gate     300      Tarsus.Length 0.5543   0.8400    0.1600
              pigauto     300        Wing.Length 0.5001   0.5579    0.3705
 pigauto_relaxed_gate     300        Wing.Length 0.5069   0.5889    0.3383
 r_cal_mean
    0.00000
    0.00000
    0.01053
    0.00000
    0.00000
    0.00000
    0.07158
    0.07275

