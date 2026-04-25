# Day 2: multitrait GNN earnings (multitrait_FIXED_tree300, n=300, K=4)

Per (couple, trait, method): mean RMSE +/- SE on held-out cells.
Decisive trait: y4 (where the nonlinear sin*exp coupling lives).

```
             method couple trait_name mean_rmse     se_rmse
        column_mean    0.0         y1 0.8595405 0.063509355
      lm_crosstrait    0.0         y1 0.6971043 0.051777516
 pigauto_multitrait    0.0         y1 0.4810188 0.009931473
    rphylopars_blup    0.0         y1 0.4443872 0.011132327
        column_mean    0.5         y1 0.7736544 0.015949636
      lm_crosstrait    0.5         y1 0.6400236 0.050053976
 pigauto_multitrait    0.5         y1 0.5109677 0.019725882
    rphylopars_blup    0.5         y1 0.4637636 0.024586842
        column_mean    1.0         y1 0.8283159 0.041509791
      lm_crosstrait    1.0         y1 0.6484442 0.071140411
 pigauto_multitrait    1.0         y1 0.5428916 0.048360962
    rphylopars_blup    1.0         y1 0.4949758 0.053284317
        column_mean    0.0         y2 0.7512847 0.030574556
      lm_crosstrait    0.0         y2 0.6410057 0.044680516
 pigauto_multitrait    0.0         y2 0.5157049 0.009944141
    rphylopars_blup    0.0         y2 0.4693064 0.027800883
        column_mean    0.5         y2 0.7876725 0.030773028
      lm_crosstrait    0.5         y2 0.6354609 0.041882849
 pigauto_multitrait    0.5         y2 0.5308469 0.034471735
    rphylopars_blup    0.5         y2 0.4733233 0.020833532
        column_mean    1.0         y2 0.8457644 0.025238146
      lm_crosstrait    1.0         y2 0.7290208 0.048221500
 pigauto_multitrait    1.0         y2 0.5831074 0.041379339
    rphylopars_blup    1.0         y2 0.5128015 0.028550858
        column_mean    0.0         y3 0.8024527 0.044586505
      lm_crosstrait    0.0         y3 0.6836752 0.067809417
 pigauto_multitrait    0.0         y3 0.4815280 0.017045575
    rphylopars_blup    0.0         y3 0.4280574 0.027876647
        column_mean    0.5         y3 0.6997813 0.027282453
      lm_crosstrait    0.5         y3 0.6078317 0.026172175
 pigauto_multitrait    0.5         y3 0.4824307 0.012236626
    rphylopars_blup    0.5         y3 0.4416865 0.010916870
        column_mean    1.0         y3 0.7309627 0.051986518
      lm_crosstrait    1.0         y3 0.6849069 0.065036028
 pigauto_multitrait    1.0         y3 0.4981394 0.041765126
    rphylopars_blup    1.0         y3 0.4448077 0.027335781
        column_mean    0.0         y4 0.7432621 0.038522852
      lm_crosstrait    0.0         y4 0.6302083 0.028245502
 pigauto_multitrait    0.0         y4 0.5208874 0.028228045
    rphylopars_blup    0.0         y4 0.4792288 0.018232123
        column_mean    0.5         y4 0.8028828 0.117010001
      lm_crosstrait    0.5         y4 0.6919519 0.066216591
 pigauto_multitrait    0.5         y4 0.5992370 0.052715996
    rphylopars_blup    0.5         y4 0.5474226 0.036610105
        column_mean    1.0         y4 1.0516603 0.027206585
      lm_crosstrait    1.0         y4 0.8920458 0.054988986
 pigauto_multitrait    1.0         y4 0.7572924 0.041701386
    rphylopars_blup    1.0         y4 0.7099185 0.028082433
```
