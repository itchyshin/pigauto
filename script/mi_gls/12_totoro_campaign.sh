#!/usr/bin/env bash
# script/mi_gls/12_totoro_campaign.sh
#
# PRIMARY route for the posterior-MI simulation campaign (arc/mi-posterior):
# runs script/mi_gls/01_cell_v2.R for every (regime, rep) cell on Totoro,
# N cells at a time, from a FROZEN copy of the code (a git archive of one
# commit, with its SHA in a file). 11_fir_array.sbatch is the DRAC fallback.
#
# Usage:
#   script/mi_gls/12_totoro_campaign.sh <code_dir> <out_dir> <n_jobs> [regimes] [reps]
#     code_dir  frozen code: `git archive <sha>` unpacked, plus a file SHA holding <sha>
#     out_dir   results (regime_<id>_rep_<k>.rds), logs/ (one log per cell), CODE_SHA
#     n_jobs    cells run in parallel, 1..150 (the Totoro cap for this lane, D-143);
#               each cell is single-threaded (BLAS pinned to 1 thread, chains serial)
#     regimes   ids or ranges from script/mi_gls/regimes.R, e.g. 1-40 (default: the planned
#               grid, 1-24 plus the in-model twins 25-40), 25-40, 22,24, 9-16,21
#     reps      N for reps 1..N (default 200), or a range a-b
#
# Cells whose .rds already exists are skipped, so rerunning the same
# command resumes. An out_dir is tied to one code SHA (out_dir/CODE_SHA);
# a different SHA is refused. The n = 1000 regimes go first (longest
# cells first; both-missing regimes, which also run posterior_none, first
# within each size). A failing cell prints
#   CELL_FAILED regime=<id> rep=<k> status=<exit code>
# and the others carry on. The last line is a CAMPAIGN_DONE count line;
# the exit status is 1 if any cell failed.
#
# Bootstrap (commit first: git archive holds committed files only):
#   Mac:    SHA=$(git rev-parse HEAD)
#           git archive --format=tar.gz -o mi_post_$SHA.tar.gz "$SHA"
#           scp mi_post_$SHA.tar.gz totoro:~/hsq_work/mi_post/
#   Totoro: C=~/hsq_work/mi_post/code_$SHA
#           mkdir -p "$C" && tar -xzf ~/hsq_work/mi_post/mi_post_$SHA.tar.gz -C "$C"
#           echo "$SHA" > "$C/SHA"
#           nohup "$C/script/mi_gls/12_totoro_campaign.sh" "$C" ~/hsq_work/mi_post/out_$SHA 150 \
#             > ~/hsq_work/mi_post/campaign_$SHA.log 2>&1 &
# code_dir must be an unpacked archive: a git checkout (a .git directory or
# worktree file) is refused, because load_all() would also source untracked
# or ignored files that no SHA records.
#
# Then summarise on Totoro into the paths the GATES.md G6/G7 CHECK lines
# read, copy those three files back into the worktree, and gate there:
#   Totoro: D=docs/dev-log/mi-posterior
#           cd "$C" && Rscript script/mi_gls/03_summarise_v2.R ~/hsq_work/mi_post/out_$SHA \
#             $D/sim_summary.csv $D/sim_summary.md $D/cell_coverage.csv
#   Mac (worktree root):
#           D=docs/dev-log/mi-posterior; R=totoro:hsq_work/mi_post/code_$SHA/$D
#           rsync -av "$R/sim_summary.csv" "$R/sim_summary.md" "$R/cell_coverage.csv" "$D/"
#           Rscript script/mi_gls/04_acceptance.R $D/sim_summary.csv     # G6
#           Rscript script/mi_gls/05_cell_coverage.R $D/cell_coverage.csv # G7
#
# Before the full run, time one cell of the slowest regime (n = 1000, both
# missing) and estimate: wall = cells x per-cell time / n_jobs.
#   script/mi_gls/12_totoro_campaign.sh "$C" ~/hsq_work/mi_post/probe_$SHA 1 24 1

set -euo pipefail

die() { echo "12_totoro_campaign: $*" >&2; exit 2; }

[[ $# -ge 3 ]] || die "usage: $0 <code_dir> <out_dir> <n_jobs> [regimes (default 1-40)] [reps (default 200)]"
CODE_DIR=$1
OUT_DIR=$2
N_JOBS=$3
REGIMES_SPEC=${4:-1-40}
REPS_SPEC=${5:-200}

# Decimal only, no leading zero, at most 3 digits: bash arithmetic reads a
# leading 0 as octal (0200 -> 128 would pass the cap) while xargs -P reads
# base 10, and a long digit string can overflow the comparison.
[[ "$N_JOBS" =~ ^[1-9][0-9]{0,2}$ ]] ||
  die "n_jobs must be a decimal integer 1..150 without a leading zero; got '$N_JOBS'"
N_JOBS=$(( 10#$N_JOBS ))
(( N_JOBS <= 150 )) || die "n_jobs=$N_JOBS refused: Totoro cap for this lane is 150 cores"

[[ -d "$CODE_DIR" ]] || die "code_dir '$CODE_DIR' does not exist"
CODE_DIR=$(cd "$CODE_DIR" && pwd)
[[ -f "$CODE_DIR/SHA" ]] || die "missing $CODE_DIR/SHA (write the archived commit's SHA there)"
[[ -f "$CODE_DIR/DESCRIPTION" && -f "$CODE_DIR/script/mi_gls/01_cell_v2.R" ]] ||
  die "$CODE_DIR is not a pigauto code tree with script/mi_gls/01_cell_v2.R"
SHA=$(tr -d '[:space:]' < "$CODE_DIR/SHA")
[[ "$SHA" =~ ^[0-9a-f]{7,40}$ ]] || die "SHA file does not hold a git SHA: '$SHA'"
# A checkout is refused outright (-e also catches a worktree's .git file):
# untracked and ignored files in it would be sourced by load_all() but are
# recorded by no SHA. Unpack a git archive instead (see the header).
[[ ! -e "$CODE_DIR/.git" ]] ||
  die "code_dir is a git checkout ($CODE_DIR/.git exists); unpack 'git archive $SHA' into a fresh directory"

mkdir -p "$OUT_DIR/logs"
OUT_DIR=$(cd "$OUT_DIR" && pwd)
if [[ -f "$OUT_DIR/CODE_SHA" ]]; then
  old=$(tr -d '[:space:]' < "$OUT_DIR/CODE_SHA")
  [[ "$old" == "$SHA" ]] || die "out_dir holds results from code $old, not $SHA; use a fresh out_dir"
else
  echo "$SHA" > "$OUT_DIR/CODE_SHA"
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MI_POST_SHA="$SHA"

# ---- spec parsing ------------------------------------------------------------
expand_ids() {   # "1-3,7" -> one id per line
  local part
  local -a parts
  IFS=',' read -ra parts <<< "$1"
  for part in "${parts[@]}"; do
    # 10# keeps a leading 0 decimal (07 -> 7), so file names match the cell's.
    if [[ "$part" =~ ^([0-9]{1,6})-([0-9]{1,6})$ ]]; then
      seq $(( 10#${BASH_REMATCH[1]} )) $(( 10#${BASH_REMATCH[2]} ))
    elif [[ "$part" =~ ^[0-9]{1,6}$ ]]; then
      echo $(( 10#$part ))
    else
      die "cannot parse '$1' (use e.g. 1-24 or 1,21)"
    fi
  done
}
REGIME_IDS=$(expand_ids "$REGIMES_SPEC" | sort -n -u | tr '\n' ' ')
if [[ "$REPS_SPEC" =~ ^[0-9]{1,6}$ ]]; then
  REPS_SPEC=$(( 10#$REPS_SPEC ))
  (( REPS_SPEC >= 1 )) || die "reps must be >= 1"
  REP_IDS=$(seq 1 "$REPS_SPEC" | tr '\n' ' ')
else
  REP_IDS=$(expand_ids "$REPS_SPEC" | sort -n -u | tr '\n' ' ')
fi
[[ -n "${REP_IDS// /}" ]] || die "no reps in '$REPS_SPEC'"

cd "$CODE_DIR"

# ---- preflight: deps (never installed here), the package loads, the API is live --
# (No backslashes in these R snippets: Rscript -e rewrites them.)
Rscript -e '
  imp <- read.dcf("DESCRIPTION", fields = "Imports")[1, 1]
  imp <- trimws(sub("[(].*$", "", strsplit(gsub("[[:space:]]+", " ", imp), ",")[[1]]))
  need <- unique(c("phylolm", "ape", "nlme", "devtools", imp[nzchar(imp)]))
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("MISSING R packages: ", paste(miss, collapse = " "))
    message("libPaths: ", paste(.libPaths(), collapse = " "))
    quit(status = 1)
  }
  suppressMessages(devtools::load_all(quiet = TRUE))
  if (!("posterior" %in% eval(formals(multi_impute)$draws_method))) {
    message("multi_impute() in this code tree has no draws_method = posterior")
    quit(status = 1)
  }' || die "preflight failed (see above); fix the library, not this script"

# Longest cells first: n = 1000 before n = 300, both-missing first within
# each size. Unknown regime ids are refused.
ORDER=$(Rscript -e '
  source(file.path("script", "mi_gls", "regimes.R"))
  ids <- as.integer(commandArgs(trailingOnly = TRUE))
  bad <- setdiff(ids, regimes$regime_id)
  if (length(bad)) { message("unknown regime id(s): ", paste(bad, collapse = " ")); quit(status = 1) }
  r <- regimes[match(ids, regimes$regime_id), ]
  cat(r$regime_id[order(-r$n, r$missing != "both", r$regime_id)])' $REGIME_IDS) ||
  die "regime list rejected"

STAMP=$(date +%Y%m%d-%H%M%S)
JOBLIST="$OUT_DIR/logs/joblist_$STAMP.txt"
FAILED_FILE="$OUT_DIR/logs/failed_$STAMP.txt"
: > "$JOBLIST"
: > "$FAILED_FILE"
N_PLANNED=0
N_SKIPPED=0
for r in $ORDER; do
  for p in $REP_IDS; do
    N_PLANNED=$(( N_PLANNED + 1 ))
    if [[ -f "$OUT_DIR/regime_${r}_rep_${p}.rds" ]]; then
      N_SKIPPED=$(( N_SKIPPED + 1 ))
    else
      echo "$r $p" >> "$JOBLIST"
    fi
  done
done
N_RUN=$(wc -l < "$JOBLIST" | tr -d ' ')

echo "CAMPAIGN_START sha=$SHA code_dir=$CODE_DIR out_dir=$OUT_DIR n_jobs=$N_JOBS"
echo "  regimes (run order)=$ORDER"
echo "  reps=$REPS_SPEC planned=$N_PLANNED skipped_existing=$N_SKIPPED to_run=$N_RUN"
echo "  host=$(hostname) nproc=$(nproc 2>/dev/null || echo NA) load=$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null || echo NA)"
echo "  started $(date)"

run_cell() {
  local r=$1 p=$2 st=0
  nice -n 10 Rscript script/mi_gls/01_cell_v2.R "$r" "$p" "$OUT_DIR" \
    > "$OUT_DIR/logs/regime_${r}_rep_${p}.log" 2>&1 || st=$?
  if (( st != 0 )); then
    echo "CELL_FAILED regime=$r rep=$p status=$st"
    echo "$r $p $st" >> "$FAILED_FILE"
  fi
  return 0   # never stop the other cells
}
export -f run_cell
export OUT_DIR FAILED_FILE

if (( N_RUN > 0 )); then
  xargs -P "$N_JOBS" -n 2 bash -c 'run_cell "$1" "$2"' _ < "$JOBLIST"
fi

N_OK=0
while read -r r p; do
  if [[ -f "$OUT_DIR/regime_${r}_rep_${p}.rds" ]]; then N_OK=$(( N_OK + 1 )); fi
done < "$JOBLIST"
N_FAILED=$(wc -l < "$FAILED_FILE" | tr -d ' ')
echo "CAMPAIGN_DONE sha=$SHA planned=$N_PLANNED skipped_existing=$N_SKIPPED ran=$N_RUN ok=$N_OK failed=$N_FAILED finished=$(date)"
(( N_FAILED == 0 )) || exit 1
