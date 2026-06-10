#!/usr/bin/env bash
# Run deseq_tables_and_vulcanos.Rmd for all 4 pairs in parallel.
# Output: analysis-1/steps/05.V.de_results/rundir/deseq_tables_and_vulcanos.{pair}.linked.html

set -euo pipefail

THISSCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$(git -C "$THISSCRIPTDIR" rev-parse --show-toplevel)"
KNIT2HTML=$(realpath scripts/knit2html)

PAIRS=(A M6 NC C)
SCRIPTS=analysis-1/steps/05.V.de_results/scripts
OUTDIR=$(realpath analysis-1/steps/05.V.de_results/rundir)
mkdir -p "$OUTDIR"

for pair in "${PAIRS[@]}"; do
  ln -sf deseq_tables_and_vulcanos.Rmd "${SCRIPTS}/deseq_tables_and_vulcanos.${pair}.Rmd"
done

cleanup() {
  for pair in "${PAIRS[@]}"; do
    rm -f "${SCRIPTS}/deseq_tables_and_vulcanos.${pair}.Rmd"
  done
}
trap cleanup EXIT

pids=()
for pair in "${PAIRS[@]}"; do
  LOG="${OUTDIR}/deseq_tables_and_vulcanos.${pair}.log"
  echo "$(date '+%F %T')  START  pair=${pair}"
  "$KNIT2HTML" "${SCRIPTS}/deseq_tables_and_vulcanos.${pair}.Rmd" "pair=${pair}" \
    > "$LOG" 2>&1 &
  pids+=($!)
done

echo "Launched ${#pids[@]} jobs: ${pids[*]}"
echo "Waiting for all to complete..."

failed=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "$(date '+%F %T')  OK     pair=${PAIRS[$i]}"
  else
    echo "$(date '+%F %T')  FAIL   pair=${PAIRS[$i]} - see ${OUTDIR}/deseq_tables_and_vulcanos.${PAIRS[$i]}.log"
    failed=$((failed + 1))
  fi
done

if [ "$failed" -gt 0 ]; then
  echo "$failed job(s) failed."
  exit 1
fi

echo "All ${#PAIRS[@]} reports done."
echo "HTML files in: ${OUTDIR}"
