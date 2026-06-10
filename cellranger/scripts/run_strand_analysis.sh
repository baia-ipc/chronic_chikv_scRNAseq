#!/bin/bash
#
# For every sample, count UMIs mapping to the CHIKV genome split by strand,
# using strand_analysis.py on the reads aligned to the CHIKV accession. Results
# are written per sample to cellranger/strand_analysis_out/, plus a combined
# summary.txt.

set -euo pipefail

CHIKV_ACCESSION="NC_004162.2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELLRANGER_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUNDIR="${CELLRANGER_DIR}/rundir"
OUT_DIR="${CELLRANGER_DIR}/strand_analysis_out"

mkdir -p "${OUT_DIR}"

for sample_dir in "${RUNDIR}"/*; do
  sample="$(basename "${sample_dir}")"
  bam="${sample_dir}/outs/possorted_genome_bam.bam"
  [ -f "${bam}" ] || continue
  echo "==== strand analysis: ${sample} ====="
  samtools view -h -b "${bam}" "${CHIKV_ACCESSION}" | \
    "${SCRIPT_DIR}/strand_analysis.py" /dev/stdin --consolidate | \
    tee "${OUT_DIR}/${sample}.txt"
done

# Consolidate results into a single report:
rm -f "${OUT_DIR}/summary.txt"
for sample_dir in "${RUNDIR}"/*; do
  sample="$(basename "${sample_dir}")"
  [ -f "${RUNDIR}/${sample}/outs/possorted_genome_bam.bam" ] || continue
  echo "Sample: ${sample}" >> "${OUT_DIR}/summary.txt"
  if [ ! -s "${OUT_DIR}/${sample}.txt" ]; then
    echo "  None" >> "${OUT_DIR}/summary.txt"
  else
    sed 's/^/  /' "${OUT_DIR}/${sample}.txt" >> "${OUT_DIR}/summary.txt"
  fi
done
