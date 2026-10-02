#!/usr/bin/env bash
# Round-2 campaign driver (design.md 5e): twins 25-40 x 200 reps, then the 44
# previously non-converged cells of regimes 1-24, all at code 9597e18.
set -u
B=/home/snakagaw/pigauto_mi_posterior/9597e18b79
C=$B/code; O=$B/campaign/sim
export R_LIBS_USER=/home/snakagaw/R/lib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
echo "PHASE1_START $(date)"
bash $C/script/mi_gls/12_totoro_campaign.sh "$C" "$O" 140 25-40 200
echo "PHASE1_END $(date)"
echo "PHASE2_START $(date)"
declare -A R=( [1]="118" [3]="154,155" [5]="59,65,120,182,183" [7]="124,176,193" [9]="110,140" [11]="37" [13]="127,132,184" [21]="15,20,28,68,73,146,157,175,182,188,190" [23]="5,6,37,59,76,79,87,88,111,144,153,170,174,186,196,198" )
for reg in "${!R[@]}"; do
  n=$(echo "${R[$reg]}" | tr "," "\n" | wc -l)
  bash $C/script/mi_gls/12_totoro_campaign.sh "$C" "$O" "$n" "$reg" "${R[$reg]}" > $B/campaign/phase2_regime_${reg}.log 2>&1 &
  sleep 2
done
wait
echo "PHASE2_END $(date)"
echo "DRIVER_DONE $(ls $O/*.rds | wc -l) rep files"
