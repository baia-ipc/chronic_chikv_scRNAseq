#!/bin/bash
#
# Expose the cellranger outs/ directories under cellranger/results/, using the
# final SampleID as name, so that the downstream analysis can read each sample
# at cellranger/results/<SampleID>/ (this directory IS the cellranger outs/).
#
# The links point into the local run directory produced by run_count.sh
# (cellranger/rundir/<SampleID>/outs).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELLRANGER_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${CELLRANGER_DIR}/results"
RUNDIR="${CELLRANGER_DIR}/rundir"

mkdir -p "${RESULTS_DIR}"

for outs in "${RUNDIR}"/*/outs; do
  [ -d "${outs}" ] || continue
  sample_id="$(basename "$(dirname "${outs}")")"
  link="${RESULTS_DIR}/${sample_id}"
  if [ -e "${link}" ]; then
    echo "Skipping ${sample_id}: ${link} already exists"
    continue
  fi
  ln -s "../rundir/${sample_id}/outs" "${link}"
  echo "Linked ${sample_id} -> ../rundir/${sample_id}/outs"
done
