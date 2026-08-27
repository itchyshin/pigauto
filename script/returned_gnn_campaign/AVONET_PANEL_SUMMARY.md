# AVONET300 multi-seed external validation panel

Generated: 2026-08-27 14:08:54 UTC
Candidate SHA: `6fddd79`
Driver SHA: `HEAD`
Dataset: avonet300 (300 species, 7 traits), miss_frac = 30% MCAR (paired per seed).
Seeds: 2026, 2027, 2028, 2029, 2030
Jobs: 15 total, 0 failures

**Framing:** Real-data corroboration only — not pre-registered in Phase A.
BACE may win on continuous morphometrics; pigauto's claim is unified mixed-type
imputation with calibrated gates. `baseline_pigauto_fixed1` is pigauto's internal
phylogenetic baseline (Rphylopars joint path when installed).

## Per-trait metrics (mean ± MCSE across seeds)

| Method | Trait | Type | Metric | Mean | MCSE | n |
|---|---|---|---|---:|---:|---:|
| pigauto_fixed1 | Mass | continuous | rmse | 0.3616 | 0.0420 | 5 |
| pigauto_fixed1 | Mass | continuous | pearson_r | 0.9517 | 0.0060 | 5 |
| pigauto_fixed1 | Mass | continuous | mae | 0.2781 | 0.0188 | 5 |
| pigauto_fixed1 | Mass | continuous | conformal_coverage_95 | 0.9670 | 0.0204 | 5 |
| pigauto_fixed1 | Beak.Length_Culmen | continuous | rmse | 0.5985 | 0.0445 | 5 |
| pigauto_fixed1 | Beak.Length_Culmen | continuous | pearson_r | 0.7760 | 0.0793 | 5 |
| pigauto_fixed1 | Beak.Length_Culmen | continuous | mae | 0.4343 | 0.0109 | 5 |
| pigauto_fixed1 | Beak.Length_Culmen | continuous | conformal_coverage_95 | 0.9299 | 0.0347 | 5 |
| pigauto_fixed1 | Tarsus.Length | continuous | rmse | 0.5049 | 0.0650 | 5 |
| pigauto_fixed1 | Tarsus.Length | continuous | pearson_r | 0.8561 | 0.0260 | 5 |
| pigauto_fixed1 | Tarsus.Length | continuous | mae | 0.3497 | 0.0360 | 5 |
| pigauto_fixed1 | Tarsus.Length | continuous | conformal_coverage_95 | 0.9608 | 0.0187 | 5 |
| pigauto_fixed1 | Wing.Length | continuous | rmse | 0.5224 | 0.0318 | 5 |
| pigauto_fixed1 | Wing.Length | continuous | pearson_r | 0.8443 | 0.0307 | 5 |
| pigauto_fixed1 | Wing.Length | continuous | mae | 0.4135 | 0.0250 | 5 |
| pigauto_fixed1 | Wing.Length | continuous | conformal_coverage_95 | 0.9326 | 0.0395 | 5 |
| pigauto_fixed1 | Trophic.Level | categorical | accuracy | 0.7326 | 0.0916 | 5 |
| pigauto_fixed1 | Trophic.Level | categorical | accuracy_Carnivore | 0.9800 | 0.0200 | 5 |
| pigauto_fixed1 | Trophic.Level | categorical | accuracy_Herbivore | 0.6143 | 0.1857 | 5 |
| pigauto_fixed1 | Trophic.Level | categorical | accuracy_Omnivore | 0.3700 | 0.1212 | 5 |
| pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy | 0.7766 | 0.0564 | 5 |
| pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Aerial | 0.8500 | 0.1000 | 5 |
| pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Generalist | 0.5000 | 0.2887 | 4 |
| pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Insessorial | 0.9618 | 0.0234 | 5 |
| pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Terrestrial | 0.5467 | 0.0871 | 5 |
| pigauto_fixed1 | Migration | ordinal | rmse | 1.0392 | 0.1142 | 5 |
| pigauto_fixed1 | Migration | ordinal | spearman_rho | 0.4703 | 0.0652 | 5 |
| pigauto_fixed1 | Migration | ordinal | conformal_coverage_95 | 0.8973 | 0.0401 | 5 |
| baseline_pigauto_fixed1 | Mass | continuous | rmse | 0.4678 | 0.0702 | 5 |
| baseline_pigauto_fixed1 | Mass | continuous | pearson_r | 0.8867 | 0.0324 | 5 |
| baseline_pigauto_fixed1 | Mass | continuous | mae | 0.3453 | 0.0412 | 5 |
| baseline_pigauto_fixed1 | Beak.Length_Culmen | continuous | rmse | 0.5984 | 0.0502 | 5 |
| baseline_pigauto_fixed1 | Beak.Length_Culmen | continuous | pearson_r | 0.7667 | 0.0788 | 5 |
| baseline_pigauto_fixed1 | Beak.Length_Culmen | continuous | mae | 0.4318 | 0.0150 | 5 |
| baseline_pigauto_fixed1 | Tarsus.Length | continuous | rmse | 0.5505 | 0.1045 | 5 |
| baseline_pigauto_fixed1 | Tarsus.Length | continuous | pearson_r | 0.8078 | 0.0477 | 5 |
| baseline_pigauto_fixed1 | Tarsus.Length | continuous | mae | 0.3672 | 0.0486 | 5 |
| baseline_pigauto_fixed1 | Wing.Length | continuous | rmse | 0.5900 | 0.0374 | 5 |
| baseline_pigauto_fixed1 | Wing.Length | continuous | pearson_r | 0.7956 | 0.0274 | 5 |
| baseline_pigauto_fixed1 | Wing.Length | continuous | mae | 0.4478 | 0.0240 | 5 |
| baseline_pigauto_fixed1 | Trophic.Level | categorical | accuracy | 0.8063 | 0.0483 | 5 |
| baseline_pigauto_fixed1 | Trophic.Level | categorical | accuracy_Carnivore | 0.9550 | 0.0278 | 5 |
| baseline_pigauto_fixed1 | Trophic.Level | categorical | accuracy_Herbivore | 0.7921 | 0.1072 | 5 |
| baseline_pigauto_fixed1 | Trophic.Level | categorical | accuracy_Omnivore | 0.3700 | 0.1212 | 5 |
| baseline_pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy | 0.7766 | 0.0564 | 5 |
| baseline_pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Aerial | 0.8500 | 0.1000 | 5 |
| baseline_pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Generalist | 0.5000 | 0.2887 | 4 |
| baseline_pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Insessorial | 0.9618 | 0.0234 | 5 |
| baseline_pigauto_fixed1 | Primary.Lifestyle | categorical | accuracy_Terrestrial | 0.5467 | 0.0871 | 5 |
| baseline_pigauto_fixed1 | Migration | ordinal | rmse | 1.0457 | 0.1203 | 5 |
| baseline_pigauto_fixed1 | Migration | ordinal | spearman_rho | 0.4600 | 0.0557 | 5 |
| pigauto_bayes | Mass | continuous | rmse | 0.3550 | 0.0382 | 5 |
| pigauto_bayes | Mass | continuous | pearson_r | 0.9501 | 0.0061 | 5 |
| pigauto_bayes | Mass | continuous | mae | 0.2723 | 0.0172 | 5 |
| pigauto_bayes | Mass | continuous | conformal_coverage_95 | 0.9761 | 0.0153 | 5 |
| pigauto_bayes | Beak.Length_Culmen | continuous | rmse | 0.6040 | 0.0438 | 5 |
| pigauto_bayes | Beak.Length_Culmen | continuous | pearson_r | 0.7675 | 0.0714 | 5 |
| pigauto_bayes | Beak.Length_Culmen | continuous | mae | 0.4359 | 0.0157 | 5 |
| pigauto_bayes | Beak.Length_Culmen | continuous | conformal_coverage_95 | 0.9299 | 0.0347 | 5 |
| pigauto_bayes | Tarsus.Length | continuous | rmse | 0.5110 | 0.0609 | 5 |
| pigauto_bayes | Tarsus.Length | continuous | pearson_r | 0.8540 | 0.0235 | 5 |
| pigauto_bayes | Tarsus.Length | continuous | mae | 0.3512 | 0.0309 | 5 |
| pigauto_bayes | Tarsus.Length | continuous | conformal_coverage_95 | 0.9608 | 0.0187 | 5 |
| pigauto_bayes | Wing.Length | continuous | rmse | 0.5383 | 0.0296 | 5 |
| pigauto_bayes | Wing.Length | continuous | pearson_r | 0.8378 | 0.0274 | 5 |
| pigauto_bayes | Wing.Length | continuous | mae | 0.4229 | 0.0185 | 5 |
| pigauto_bayes | Wing.Length | continuous | conformal_coverage_95 | 0.9150 | 0.0346 | 5 |
| pigauto_bayes | Trophic.Level | categorical | accuracy | 0.7326 | 0.0916 | 5 |
| pigauto_bayes | Trophic.Level | categorical | accuracy_Carnivore | 0.9800 | 0.0200 | 5 |
| pigauto_bayes | Trophic.Level | categorical | accuracy_Herbivore | 0.6143 | 0.1857 | 5 |
| pigauto_bayes | Trophic.Level | categorical | accuracy_Omnivore | 0.3700 | 0.1212 | 5 |
| pigauto_bayes | Primary.Lifestyle | categorical | accuracy | 0.7766 | 0.0564 | 5 |
| pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Aerial | 0.8500 | 0.1000 | 5 |
| pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Generalist | 0.5000 | 0.2887 | 4 |
| pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Insessorial | 0.9618 | 0.0234 | 5 |
| pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Terrestrial | 0.5467 | 0.0871 | 5 |
| pigauto_bayes | Migration | ordinal | rmse | 1.0849 | 0.0817 | 5 |
| pigauto_bayes | Migration | ordinal | spearman_rho | 0.4552 | 0.0279 | 5 |
| pigauto_bayes | Migration | ordinal | conformal_coverage_95 | 0.9183 | 0.0428 | 5 |
| baseline_pigauto_bayes | Mass | continuous | rmse | 0.4663 | 0.0644 | 5 |
| baseline_pigauto_bayes | Mass | continuous | pearson_r | 0.8907 | 0.0290 | 5 |
| baseline_pigauto_bayes | Mass | continuous | mae | 0.3545 | 0.0388 | 5 |
| baseline_pigauto_bayes | Beak.Length_Culmen | continuous | rmse | 0.5951 | 0.0477 | 5 |
| baseline_pigauto_bayes | Beak.Length_Culmen | continuous | pearson_r | 0.7630 | 0.0775 | 5 |
| baseline_pigauto_bayes | Beak.Length_Culmen | continuous | mae | 0.4293 | 0.0163 | 5 |
| baseline_pigauto_bayes | Tarsus.Length | continuous | rmse | 0.5314 | 0.0900 | 5 |
| baseline_pigauto_bayes | Tarsus.Length | continuous | pearson_r | 0.8277 | 0.0333 | 5 |
| baseline_pigauto_bayes | Tarsus.Length | continuous | mae | 0.3605 | 0.0452 | 5 |
| baseline_pigauto_bayes | Wing.Length | continuous | rmse | 0.5932 | 0.0320 | 5 |
| baseline_pigauto_bayes | Wing.Length | continuous | pearson_r | 0.7974 | 0.0285 | 5 |
| baseline_pigauto_bayes | Wing.Length | continuous | mae | 0.4559 | 0.0203 | 5 |
| baseline_pigauto_bayes | Trophic.Level | categorical | accuracy | 0.8063 | 0.0483 | 5 |
| baseline_pigauto_bayes | Trophic.Level | categorical | accuracy_Carnivore | 0.9550 | 0.0278 | 5 |
| baseline_pigauto_bayes | Trophic.Level | categorical | accuracy_Herbivore | 0.7921 | 0.1072 | 5 |
| baseline_pigauto_bayes | Trophic.Level | categorical | accuracy_Omnivore | 0.3700 | 0.1212 | 5 |
| baseline_pigauto_bayes | Primary.Lifestyle | categorical | accuracy | 0.7766 | 0.0564 | 5 |
| baseline_pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Aerial | 0.8500 | 0.1000 | 5 |
| baseline_pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Generalist | 0.5000 | 0.2887 | 4 |
| baseline_pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Insessorial | 0.9618 | 0.0234 | 5 |
| baseline_pigauto_bayes | Primary.Lifestyle | categorical | accuracy_Terrestrial | 0.5467 | 0.0871 | 5 |
| baseline_pigauto_bayes | Migration | ordinal | rmse | 1.0905 | 0.0806 | 5 |
| baseline_pigauto_bayes | Migration | ordinal | spearman_rho | 0.4502 | 0.0280 | 5 |
| bace | Mass | continuous | rmse | 1759.8047 | 217.1571 | 5 |
| bace | Mass | continuous | pearson_r | 0.6359 | 0.1135 | 5 |
| bace | Beak.Length_Culmen | continuous | rmse | 20.0210 | 2.1808 | 5 |
| bace | Beak.Length_Culmen | continuous | pearson_r | 0.8586 | 0.0170 | 5 |
| bace | Tarsus.Length | continuous | rmse | 22.2066 | 3.7346 | 5 |
| bace | Tarsus.Length | continuous | pearson_r | 0.8798 | 0.0265 | 5 |
| bace | Wing.Length | continuous | rmse | 52.3590 | 3.4050 | 5 |
| bace | Wing.Length | continuous | pearson_r | 0.9128 | 0.0071 | 5 |
| bace | Trophic.Level | categorical | accuracy | 0.5089 | 0.0184 | 5 |
| bace | Primary.Lifestyle | categorical | accuracy | 0.3378 | 0.0171 | 5 |
| bace | Migration | ordinal | accuracy | 0.6267 | 0.0215 | 5 |

## Bayes λ sensitivity (continuous traits: pigauto_bayes vs pigauto_fixed1)

| Trait | RMSE fixed_1 | RMSE bayes | Δ (bayes − fixed) |
|---|---:|---:|---:|
| Mass | 0.3616 | 0.3550 | -0.0066 |
| Beak.Length_Culmen | 0.5985 | 0.6040 | 0.0055 |
| Tarsus.Length | 0.5049 | 0.5110 | 0.0061 |
| Wing.Length | 0.5224 | 0.5383 | 0.0159 |

## pigauto blend vs internal baseline (fixed_1)

| Trait | Type | Metric | Blend | Baseline | Δ (blend − base) |
|---|---|---|---:|---:|---:|
| Mass | continuous | rmse | 0.3616 | 0.4678 | 0.1062 |
| Beak.Length_Culmen | continuous | rmse | 0.5985 | 0.5984 | -1e-04 |
| Tarsus.Length | continuous | rmse | 0.5049 | 0.5505 | 0.0456 |
| Wing.Length | continuous | rmse | 0.5224 | 0.5900 | 0.0676 |
| Trophic.Level | categorical | accuracy | 0.7326 | 0.8063 | 0.0737 |
| Primary.Lifestyle | categorical | accuracy | 0.7766 | 0.7766 | 0.0000 |

## BACE vs pigauto_fixed1 (held-out cells)

**Scale note:** pigauto RMSE is on log/z-transformed latent scale; BACE RMSE is on
raw trait units — do not compare RMSE across methods. Pearson *r* and discrete
accuracy are on comparable scales.

### Continuous traits (Pearson *r*; higher is better)

| Trait | pigauto | BACE | Δ (BACE − pigauto) | Winner |
|---|---:|---:|---:|---|
| Mass | 0.9517 | 0.6359 | -0.3158 | pigauto |
| Beak.Length_Culmen | 0.7760 | 0.8586 | 0.0826 | BACE |
| Tarsus.Length | 0.8561 | 0.8798 | 0.0237 | BACE |
| Wing.Length | 0.8443 | 0.9128 | 0.0685 | BACE |

### Discrete traits (accuracy; higher is better)

| Trait | Type | pigauto | BACE | Δ (BACE − pigauto) |
|---|---|---|---:|---:|---:|
| Trophic.Level | categorical | 0.7326 | 0.5089 | -0.2237 |
| Primary.Lifestyle | categorical | 0.7766 | 0.3378 | -0.4388 |

## Wall time by method (mean seconds per job)

```
          method  fit_sec
1           bace  42.0698
2  pigauto_bayes 122.6188
3 pigauto_fixed1 116.3962
[1] FALSE
```

Artifacts: `script/returned_gnn_campaign/results_avonet_panel`
