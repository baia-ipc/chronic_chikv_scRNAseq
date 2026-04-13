#!/usr/bin/env bash
# Run all three CellChat Rmds × 4 pairs × 2 filter modes = 24 jobs in parallel.
#
# Output HTMLs: rundir/cellchat_{results,circle,hierarchy}.{pair}.{sfx}.html
# PDFs:         results/{pair}/{sfx}/plots/

set -euo pipefail

cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"

SCRIPTS=analysis-1/steps/07.V.cellchat_results/scripts
OUTDIR=$(realpath analysis-1/steps/07.V.cellchat_results/rundir)
ROOT=$(pwd)

RMDS=(cellchat_results cellchat_circle cellchat_hierarchy)
PAIRS=(A_NC_vs_C 6m_NC_vs_C NC_A_vs_6m C_A_vs_6m)
SFXS=(all pfilt)

pids=()
labels=()

for rmd in "${RMDS[@]}"; do
  for pair in "${PAIRS[@]}"; do
    for sfx in "${SFXS[@]}"; do
      OUT_HTML="${OUTDIR}/${rmd}.${pair}.${sfx}.html"
      LOG_FILE="${OUTDIR}/${rmd}.${pair}.${sfx}.log"
      echo "$(date '+%F %T')  START  ${rmd}  pair=${pair}  sfx=${sfx}"
      INTER_DIR="/tmp/rmd_intermediates/${rmd}.${pair}.${sfx}"
      mkdir -p "${INTER_DIR}"
      Rscript -e "
        withr::with_dir('${ROOT}', {
          rmarkdown::render(
            '${SCRIPTS}/${rmd}.Rmd',
            output_format    = 'html_document',
            output_file      = '${OUT_HTML}',
            intermediates_dir = '${INTER_DIR}',
            params           = list(pair='${pair}', sfx='${sfx}')
          )
        })
      " > "${LOG_FILE}" 2>&1 &
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

echo "All 24 reports done."
echo "HTML files in: ${OUTDIR}"
