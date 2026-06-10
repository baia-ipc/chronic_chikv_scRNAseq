#!/bin/bash
#
# Collect the cellranger metrics_summary.csv of every sample into a single QC
# table (cellranger/qc_output/metrics_summaries.tsv), annotated with the sample
# metadata and with the metric columns renamed to short identifiers.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELLRANGER_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PRJ_DIR="$(cd "${CELLRANGER_DIR}/.." && pwd)"

RUNDIR="${CELLRANGER_DIR}/rundir"
QC_DIR="${CELLRANGER_DIR}/qc_output"
METADATA="${PRJ_DIR}/metadata/samples.tsv"
COLNAMES="${CELLRANGER_DIR}/config/metrics_summary_fields_renaming.tsv"

mkdir -p "${QC_DIR}"

"${SCRIPT_DIR}/collect_metric_summaries.py" \
  "${RUNDIR}" \
  --metadata "${METADATA}" \
  --colnames "${COLNAMES}" \
  > "${QC_DIR}/metrics_summaries.tsv"
