# TabPFN Reconciliation Benchmark

- Generated: 2026-06-11 06:52:42 MDT
- Git commit: `e0329d8503f2cbc8a79997a8182d8914f6db8dba`
- Scales: `300`
- Regimes: `same_row, same_row_lappe, phylo_only`
- Replicates: `3`
- Split mode: `row_block_all`
- Dry run: `FALSE`
- Pigauto configs: `default, relaxed_gate`

## Aim

This benchmark uses row-block holdouts so validation/test species rows have their selected trait cells hidden together. In `row_block_all`, every latent column for held-out rows is masked, so same-row TabPFN features cannot borrow observed cells from the same species at prediction time.

## Status Counts

                method status Freq
              baseline     ok   12
               pigauto     ok   12
  pigauto_relaxed_gate     ok   12
     tabpfn_phylo_only     ok   12
       tabpfn_same_row     ok   12
 tabpfn_same_row_lappe     ok   12

## Mean RMSE By Method

                method scale_n   rmse pearson_r coverage_95 wall_sec
              baseline     300 0.6633    0.8307          NA   0.4673
               pigauto     300 0.6991    0.7956      0.9627 210.1267
  pigauto_relaxed_gate     300 0.7090    0.7833      0.9627 186.9593
     tabpfn_phylo_only     300 0.8656    0.6716      0.9605  10.5470
 tabpfn_same_row_lappe     300 1.1453    0.5270      0.9671  11.4098
       tabpfn_same_row     300 1.1634        NA      0.9825  11.5854

## Mean RMSE By Trait

                method scale_n              trait   rmse pearson_r coverage_95
              baseline     300 Beak.Length_Culmen 0.7030    0.7981          NA
               pigauto     300 Beak.Length_Culmen 0.7713    0.7100      0.9649
  pigauto_relaxed_gate     300 Beak.Length_Culmen 0.7835    0.6858      0.9649
     tabpfn_phylo_only     300 Beak.Length_Culmen 0.9051    0.5995      0.9912
 tabpfn_same_row_lappe     300 Beak.Length_Culmen 1.0730    0.4701      0.9912
       tabpfn_same_row     300 Beak.Length_Culmen 1.1462        NA      0.9912
               pigauto     300               Mass 0.6529    0.8379      0.9649
              baseline     300               Mass 0.6529    0.8379          NA
  pigauto_relaxed_gate     300               Mass 0.6730    0.8235      0.9649
     tabpfn_phylo_only     300               Mass 0.8288    0.7153      0.9561
       tabpfn_same_row     300               Mass 1.1218        NA      0.9737
 tabpfn_same_row_lappe     300               Mass 1.1402    0.6107      0.9386
              baseline     300      Tarsus.Length 0.6440    0.8478          NA
               pigauto     300      Tarsus.Length 0.6440    0.8478      0.9649
  pigauto_relaxed_gate     300      Tarsus.Length 0.6440    0.8478      0.9649
     tabpfn_phylo_only     300      Tarsus.Length 0.8719    0.6881      0.9649
 tabpfn_same_row_lappe     300      Tarsus.Length 1.1458    0.4745      0.9912
       tabpfn_same_row     300      Tarsus.Length 1.2098        NA      0.9912
              baseline     300        Wing.Length 0.6531    0.8393          NA
               pigauto     300        Wing.Length 0.7284    0.7867      0.9561
  pigauto_relaxed_gate     300        Wing.Length 0.7356    0.7763      0.9561
     tabpfn_phylo_only     300        Wing.Length 0.8566    0.6834      0.9298
       tabpfn_same_row     300        Wing.Length 1.1757        NA      0.9737
 tabpfn_same_row_lappe     300        Wing.Length 1.2221    0.5525      0.9474

## Pigauto Gate Audit

               method scale_n              trait   rmse r_cal_bm r_cal_gnn
              pigauto     300 Beak.Length_Culmen 0.7713   0.8000    0.2000
 pigauto_relaxed_gate     300 Beak.Length_Culmen 0.7835   0.7667    0.2333
              pigauto     300               Mass 0.6529   1.0000    0.0000
 pigauto_relaxed_gate     300               Mass 0.6730   0.9000    0.1000
              pigauto     300      Tarsus.Length 0.6440   1.0000    0.0000
 pigauto_relaxed_gate     300      Tarsus.Length 0.6440   1.0000    0.0000
              pigauto     300        Wing.Length 0.7284   0.6667    0.3333
 pigauto_relaxed_gate     300        Wing.Length 0.7356   0.6667    0.3333
 r_cal_mean
          0
          0
          0
          0
          0
          0
          0
          0

