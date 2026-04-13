#!/usr/bin/env bash
# Run deseq_tables_and_vulcanos.Rmd for all 4 pairs in parallel = 4 HTML reports.
# Output: analysis-1/steps/05.V.de_results/rundir/deseq_tables_and_vulcanos.{pair}.html

set -euo pipefail

cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"

RMD=analysis-1/steps/05.V.de_results/scripts/deseq_tables_and_vulcanos.Rmd
OUTDIR=$(realpath analysis-1/steps/05.V.de_results/rundir)
ROOT=$(pwd)

PAIRS=(A M6 NC C)

pids=()

for pair in "${PAIRS[@]}"; do
  OUT_HTML="${OUTDIR}/deseq_tables_and_vulcanos.${pair}.html"
  LOG_FILE="${OUTDIR}/deseq_tables_and_vulcanos.${pair}.log"
  echo "$(date '+%F %T')  START  pair=${pair}"
  Rscript -e "
    withr::with_dir('${ROOT}', {
      rmarkdown::render(
        '${RMD}',
        output_format = 'html_document',
        output_file   = '${OUT_HTML}',
        params        = list(pair='${pair}')
      )
    })
  " > "${LOG_FILE}" 2>&1 &
  pids+=($!)
done

echo "Launched ${#pids[@]} jobs: ${pids[*]}"
echo "Waiting for all to complete..."

failed=0
for i in "${!pids[@]}"; do
  pid=${pids[$i]}
  if wait "$pid"; then
    echo "$(date '+%F %T')  OK     pair=${PAIRS[$i]}"
  else
    echo "$(date '+%F %T')  FAIL   pair=${PAIRS[$i]} — see ${OUTDIR}/deseq_tables_and_vulcanos.${PAIRS[$i]}.log"
    failed=$((failed+1))
  fi
done

if [ "$failed" -gt 0 ]; then
  echo "$failed job(s) failed."
  exit 1
fi

echo "All 4 reports done."
echo "HTML files in: ${OUTDIR}"
