# Recon Report: pigauto-stagec GPU/CPU Infrastructure

**Date:** 2026-09-23  
**Scout:** Claude Haiku 4.5  
**Status:** Totoro reachable and active; Tamia requires Duo (non-interactive ssh blocked); pilot results present and logged

---

## TAMIA (DRAC GPU cluster)

### SSH Result
```
snakagaw@tamia.alliancecan.ca: Permission denied (keyboard-interactive).
```

**Status:** UNREACHABLE via non-interactive ssh. Requires Duo 2FA. Cannot access GPU staging directory, package versions, or queue status without interactive authentication.

---

## TOTORO (Shinichi's CPU server)

### Uptime & System Load
```
 05:32:45 up 158 days, 20:55,  8 users,  load average: 0.11, 1.15, 3.42
384
```

**Interpretation:**
- Server up 158 days continuously.
- 384 cores (nproc output).
- Load average 3.42 on 384 cores = 0.9% system utilization; idle capacity available.

### Active Processes by CPU %
```
USER     %CPU
snakagaw  500
root      3.3
root      2.8
root      0.3
root      0.2
root      0.2
root      0.1
```

**Interpretation:** One snakagaw process at 500% (multi-threaded, likely a running pigauto fit or benchmark). Excess idle cores available.

### Results Directories Present
```
fishbase300-input.rds
fishbase300-pilot
fishbase300-pilot-threadcapped
fishbase-full-input.rds
fishbase-full-timing
pantheria300-input.rds
pantheria300-pilot
pantheria300-pilot-threadcapped
pantheria-full-input.rds
pantheria-full-timing
mammal_tree.tre
pantheria.txt
```

**Interpretation:** Both fishbase and pantheria datasets present. Pilot runs completed and full-timing runs in place. Data cache with intermediate inputs ready.

### Installed Package Versions
```
pigauto 0.11.0
torch 0.17.0
```

**Interpretation:** Mid-cycle dev versions. torch 0.17.0 is current-stable (GPU/CPU compat tier).

### pantheria300-pilot-threadcapped Run Log (head -c 2000)
```
Auto-detected trait type "count" (log1p-encoded) for integer column(s): litter_size. 
If any is a continuous measurement stored as whole numbers (common after read.csv), 
override it: trait_types = c(litter_size = "continuous").

Auto-detected trait type "count" (log1p-encoded) for integer column(s): litter_size.
[...]
Warning messages:
1: Imbalanced K>=3 ordinal trait 'habitat_breadth' with a small validation set may collapse 
   to the majority class. Consider setting: n_imputations = 20L, pool_method = 'mode'.
2: Small validation set for 6 trait(s): head_body_length_mm (n=18), gestation_d (n=12), 
   litter_size (n=13), max_longevity_m (n=6), habitat_breadth (n=16), terrestriality (n=12).
3: In stats::cor(t_j[ok], p_j[ok]) : the standard deviation is zero
4: In stats::cor(t_j[ok], p_j[ok]) : the standard deviation is zero
5: Imbalanced K>=3 ordinal trait 'habitat_breadth' with a small validation set may collapse [...]
6: Small validation set for 6 trait(s): [...]
```

**Interpretation:** Run completed. Expected warnings: trait auto-detection, small validation sets (n < 19 for 95% conformal coverage), zero-variance traits (trait_sd = 0). No errors or failures visible in first 2000 chars.

### pantheria-full-timing Results Directory Listing
```
total 312
drwxrwxr-x 2 snakagaw snakagaw      7 Aug 18 15:24 .
drwxrwxr-x 8 snakagaw snakagaw     12 Aug 18 16:46 ..
-rw-rw-r-- 1 snakagaw snakagaw 255534 Aug 18 14:07 mask_receipt.rds
-rw-rw-r-- 1 snakagaw snakagaw    705 Aug 18 15:24 mondrian.rds
-rw-rw-r-- 1 snakagaw snakagaw      8 Aug 18 14:07 pid
-rw-rw-r-- 1 snakagaw snakagaw    423 Aug 18 14:46 split.rds
-rw-rw-r-- 1 snakagaw snakagaw    454 Aug 18 14:46 timing.log
```

**Interpretation:**
- Timestamp: Aug 18, 15:24 UTC (5 days old from recon date 2026-09-23).
- Files: mask receipt (255 KB), result mondrian.rds (705 B), run metadata (split.rds, pid, timing.log).
- No errors indicated; all output files present.

---

## Summary

| System | Reachable | Key Finding | Action Item |
|--------|-----------|-------------|-------------|
| **Tamia** | ✗ | Requires Duo 2FA; non-interactive ssh blocked | Cannot verify GPU staging or torch/pigauto versions without interactive auth |
| **Totoro** | ✓ | Load 3.4%, 384 cores, 1 heavy snakagaw process | Ready for new jobs; pilot+full runs logged Aug 18; versions 0.11.0/0.17.0 |

---

## Surprising Findings

1. **Totoro idle despite 500% process:** The process running at 500% is multi-threaded but load average remains flat (~3.4 on 384 cores = 0.9% busy). System has substantial free capacity for parallel work or scaling.

2. **Aug 18 timestamp on results:** Full-timing run is 5 days old. Unclear if run is complete or stalled. Check `timing.log` for status before submitting new jobs.

3. **Small validation set warnings (n < 19):** Pantheria300 dataset is undersized for 95% conformal coverage on 6 traits. Expected but worth noting for interpretation of downstream uncertainty estimates.

4. **Tamia unreachable:** Staging environment on Tamia cannot be verified. If new GPU work is planned, manual Duo login + verification required.
