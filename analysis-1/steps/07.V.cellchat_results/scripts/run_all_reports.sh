#!/usr/bin/env bash
# Run all three CellChat Rmds × 4 pairs × 2 filter modes = 24 jobs in parallel.
#
# Output HTMLs: rundir/cellchat_{results,circle,hierarchy}.{pair}.{sfx}.linked.html

set -euo pipefail

THISSCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$(git -C "$THISSCRIPTDIR" rev-parse --show-toplevel)"
KNIT2HTML=$(realpath scripts/knit2html)

RMDS=(cellchat_results cellchat_circle cellchat_hierarchy)
PAIRS=(A_C_vs_NC 6m_C_vs_NC NC_A_vs_6m C_A_vs_6m)
SFXS=(all pfilt)

SCRIPTS=analysis-1/steps/07.V.cellchat_results/scripts
OUTDIR=$(realpath analysis-1/steps/07.V.cellchat_results/rundir)
mkdir -p "$OUTDIR"

# Copy Rmds to rundir with per-combination names so knit2html produces distinct outputs.
# Copies (not symlinks) avoid broken-symlink failures in rmarkdown's post-knit processing.
for rmd in "${RMDS[@]}"; do
  for pair in "${PAIRS[@]}"; do
    for sfx in "${SFXS[@]}"; do
      cp "${SCRIPTS}/${rmd}.Rmd" "${OUTDIR}/${rmd}.${pair}.${sfx}.Rmd"
    done
  done
done

# Remove copies on exit (success or failure)
cleanup() {
  for rmd in "${RMDS[@]}"; do
    for pair in "${PAIRS[@]}"; do
      for sfx in "${SFXS[@]}"; do
        rm -f "${OUTDIR}/${rmd}.${pair}.${sfx}.Rmd"
      done
    done
  done
}
trap cleanup EXIT

pids=()
labels=()

for rmd in "${RMDS[@]}"; do
  for pair in "${PAIRS[@]}"; do
    for sfx in "${SFXS[@]}"; do
      LOG="${OUTDIR}/${rmd}.${pair}.${sfx}.log"
      echo "$(date '+%F %T')  START  ${rmd}  pair=${pair}  sfx=${sfx}"
      "$KNIT2HTML" "${OUTDIR}/${rmd}.${pair}.${sfx}.Rmd" "pair=${pair}" "sfx=${sfx}" "filter_mode=${sfx}" \
        > "$LOG" 2>&1 &
      pids+=($!)
      labels+=("${rmd}.${pair}.${sfx}")
    done
  done
done

echo "Launched ${#pids[@]} jobs."
echo "Waiting for all to complete..."

failed=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "$(date '+%F %T')  OK     ${labels[$i]}"
  else
    echo "$(date '+%F %T')  FAIL   ${labels[$i]} — see ${OUTDIR}/${labels[$i]}.log"
    failed=$((failed+1))
  fi
done

if [ "$failed" -gt 0 ]; then
  echo "$failed job(s) failed."
  exit 1
fi

echo "All ${#pids[@]} reports done."
echo "HTML files in: ${OUTDIR}"
