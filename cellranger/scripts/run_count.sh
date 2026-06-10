#!/bin/bash
#
# Run cellranger count for every sample in the project metadata table.
#
# Each sample produces a run directory cellranger/rundir/<SampleID> whose outs/
# subdirectory holds the results. Already-computed samples are skipped, so this
# is safe to re-run. Afterwards run ./link_results.sh to expose the outs/
# directories under cellranger/results/, then ./run_qc.sh.
#
# Configure the external inputs via environment variables (see config.sh).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

CELLRANGER_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PRJ_DIR="$(cd "${CELLRANGER_DIR}/.." && pwd)"

METADATA="${PRJ_DIR}/metadata/samples.tsv"
RUNDIR="${CELLRANGER_DIR}/rundir"

mkdir -p "${RUNDIR}"

"${SCRIPT_DIR}/run_cellranger_count.py" \
  "${RUNDIR}" \
  "${METADATA}" \
  "${TRANSCRIPTOME}" \
  "${READS_DIR}" \
  --localcores "${LOCALCORES}" \
  --localvmem "${LOCALVMEM}" \
  --verbose "$@"
