#!/usr/bin/env bash
# One-line-per-host status of the rubin-freq-bace campaign (run from the Mac; uses the ~/.ssh/cm-* sockets).
# Targets (option B): BACE n100 1200, n300 1200, n1000 600 fits; freq 1200 per n (3600). Each fit = one rds.
# Frequentist results of record come from fir ~/pigauto_rubin_f3 (freq-v3); the earlier freq results are superseded.
#   bash script/rubin_status.sh
set -uo pipefail
q() { ssh -o BatchMode=yes -o ConnectTimeout=15 "$1" "$2" 2>/dev/null || echo "$1: unreachable"; }
cnt='for a in bace freq; do for n in 100 300 1000; do c=$(find "$R/results/$a" -name "*_n${n}_M20_*.rds" 2>/dev/null | wc -l); [ "$c" -gt 0 ] && printf "%s_n%s=%s " $a $n $c; done; done; printf "queued=%s failed=%s\n" "$(squeue -u $USER -h -r 2>/dev/null | wc -l)" "$(sacct -u $USER -S 2026-09-24T13:30 -n -o JobName%30,State%20 2>/dev/null | grep rubin_ | grep -cE "FAILED|TIMEOUT|OUT_OF_ME|CANCELLED")"'
echo "nibi    $(q nibi "R=~/projects/def-snakagaw/snakagaw/pigauto_rubin; $cnt") | freq-cmp $(q nibi 'find ~/projects/def-snakagaw/snakagaw/pigauto_rubin_freq/results/freq -name "*.rds" | wc -l')"
echo "fir     $(q fir "R=~/pigauto_rubin; $cnt") | freq-v3 $(q fir "R=~/pigauto_rubin_f3; $cnt")"
echo "rorqual $(q rorqual "R=~/projects/def-snakagaw/snakagaw/pigauto_rubin; $cnt")"
echo "totoro  $(q totoro 'c=$(find ~/pigauto_rubin/results/bace -name "*_n1000_M20_*.rds" 2>/dev/null | wc -l); r=$(pgrep -u $USER -fc "rubin_cell.R"); echo "bace_n1000=$c running=$r"')"
