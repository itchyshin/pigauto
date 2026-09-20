#!/usr/bin/env bash
# Pool per-cell results from every machine into one directory, then aggregate and verify.
#
#   bash script/campaign_sim_pool.sh core     ~/pigauto_sim/pooled
#   bash script/campaign_sim_pool.sh factorial ~/pigauto_sim/pooled
#
# Run from the Mac. The fast arms land on Totoro and the Bayesian arm on the DRAC clusters, so no
# single machine holds a complete cell: the completeness gates only mean anything against the pool.
# rsync is additive and idempotent, so this is safe to re-run while the campaign is still going.
#
# Results are kept in PER-HOST subdirectories, never flattened. The same (cell, seed) exists on two
# machines carrying DIFFERENT arms under the SAME filename: Totoro holds the fast arms, the clusters
# hold BACE. Flattening with --ignore-existing silently discarded every Bayesian result (measured
# 2026-09-20: 1034 nibi cells pooled to 0). The aggregator and the completeness gates recurse.
set -euo pipefail
STAGE="${1:?stage: core|factorial|avonet}"
DEST="${2:?destination directory on this machine}"
mkdir -p "$DEST/$STAGE"
SSH="ssh -o BatchMode=yes -o ConnectTimeout=20"

pull () {   # host root label
  local host="$1" root="$2" tag="${1#*@}"; tag="${tag%%.*}"
  if $SSH "$host" "[ -d $root/results/$STAGE ]" 2>/dev/null; then
    local n
    n=$($SSH "$host" "ls $root/results/$STAGE/*.rds 2>/dev/null | wc -l" || echo 0)
    printf '  %-34s %6s cells\n' "$host" "$n"
    mkdir -p "$DEST/$STAGE/$tag"
    rsync -az -e "$SSH" "$host:$root/results/$STAGE/" "$DEST/$STAGE/$tag/" 2>/dev/null || true
  else
    printf '  %-34s %6s\n' "$host" "(no $STAGE dir)"
  fi
}

echo "pooling stage=$STAGE into $DEST/$STAGE"
pull snakagaw@totoro.biology.ualberta.ca '~/pigauto_sim'
pull snakagaw@nibi.alliancecan.ca        '~/projects/def-snakagaw/snakagaw/pigauto_sim'
pull snakagaw@fir.alliancecan.ca         '/home/snakagaw/pigauto_sim'
# rorqual held 9 salvaged cells before it was retired for quota; harmless to re-check.
pull snakagaw@rorqual.alliancecan.ca     '/project/def-snakagaw/snakagaw/pigauto_sim'
echo "pooled total: $(find "$DEST/$STAGE" -name '*.rds' 2>/dev/null | wc -l) cell-files across $(ls -1 "$DEST/$STAGE" 2>/dev/null | wc -l) host dirs"
