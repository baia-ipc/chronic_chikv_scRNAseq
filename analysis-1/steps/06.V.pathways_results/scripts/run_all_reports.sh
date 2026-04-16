#!/usr/bin/env bash
# Run all 06.V pathway visualization reports in parallel.
# Each database × analysis combination gets its own output HTML file in rundir/.
# Output: analysis-1/steps/06.V.pathways_results/rundir/<db>.<analysis>.linked.html

set -euo pipefail

THISSCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$(git -C "$THISSCRIPTDIR" rev-parse --show-toplevel)"
KNIT2HTML=$(realpath scripts/knit2html)

SCRIPTS=$(realpath analysis-1/steps/06.V.pathways_results/scripts)
OUTDIR=$(realpath analysis-1/steps/06.V.pathways_results/rundir)
mkdir -p "$OUTDIR"

DATABASES=(GO KEGG MSigDb Reactome)
ANALYSES=(ORA GSEA)

pids=()
report_names=()

for analysis in "${ANALYSES[@]}"; do
  for db in "${DATABASES[@]}"; do
    # Create a database-specific symlink so each run gets a distinct output filename
    DEST_RMD="${OUTDIR}/viz_pathway_analysis.${db}.${analysis}.Rmd"
    ln -sf "${SCRIPTS}/viz_pathway_analysis.${analysis}.Rmd" "$DEST_RMD"
    LOG="${OUTDIR}/viz_pathway_analysis.${db}.${analysis}.log"
    echo "$(date '+%F %T')  START  db=${db} analysis=${analysis}"
    "$KNIT2HTML" "$DEST_RMD" "database=${db}" > "$LOG" 2>&1 &
    pids+=($!)
    report_names+=("${db}.${analysis}")
  done
done

# KEGGMaps runs from scripts/ directly (no database param needed)
LOG="${OUTDIR}/viz_pathway_analysis.KEGGMaps.log"
echo "$(date '+%F %T')  START  KEGGMaps"
"$KNIT2HTML" "${SCRIPTS}/viz_pathway_analysis.KEGGMaps.Rmd" > "$LOG" 2>&1 &
pids+=($!)
report_names+=("KEGGMaps")

echo "Launched ${#pids[@]} jobs: ${pids[*]}"
echo "Waiting for all to complete..."

failed=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "$(date '+%F %T')  OK     ${report_names[$i]}"
  else
    echo "$(date '+%F %T')  FAIL   ${report_names[$i]} — see ${OUTDIR}/viz_pathway_analysis.${report_names[$i]}.log"
    failed=$((failed+1))
  fi
done

if [ "$failed" -gt 0 ]; then
  echo "$failed job(s) failed."
  exit 1
fi

echo "All ${#report_names[@]} reports done."
echo "HTML files in: ${OUTDIR}"
