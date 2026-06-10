#!/bin/bash
#
# Build the Cell Ranger reference transcriptome with the CHIKV genome added to
# the human GRCh38-2020-A reference. The resulting "hsa_chikv" directory is the
# transcriptome to pass to cellranger count (CHIKV_TRANSCRIPTOME in config.sh).
#
# Inputs (external, machine specific -> provide via environment variables):
#
#   CHIKV_HSA_REFDATA   10x prebuilt human reference directory
#                       "refdata-gex-GRCh38-2020-A" (download from 10x Genomics),
#                       containing fasta/genome.fa and genes/genes.gtf.
#
#   CHIKV_GENOME_DIR    Directory with the CHIKV RefSeq genome files
#                         GCF_000854045.1_ViralProj14998_genomic.fna.gz
#                         GCF_000854045.1_ViralProj14998_genomic.gtf.gz
#                       (RefSeq assembly GCF_000854045.1, accession NC_004162.2;
#                        ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/854/045/
#                        GCF_000854045.1_ViralProj14998/).
#
# Compute resources are taken from config.sh (LOCALCORES / LOCALVMEM).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELLRANGER_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

HSAREF_DIR="${CHIKV_HSA_REFDATA:?Set CHIKV_HSA_REFDATA to the refdata-gex-GRCh38-2020-A directory}"
CHIKVREF_DIR="${CHIKV_GENOME_DIR:?Set CHIKV_GENOME_DIR to the CHIKV RefSeq genome directory}"
LOCALCORES="${CHIKV_LOCALCORES:-60}"
LOCALVMEM="${CHIKV_LOCALVMEM:-200}"

OUTDIR="${CELLRANGER_DIR}/reference"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Combined genome FASTA: human + CHIKV
cp "${CHIKVREF_DIR}/GCF_000854045.1_ViralProj14998_genomic.fna.gz" CHIKV.fa.gz
gunzip -f CHIKV.fa.gz
cat "${HSAREF_DIR}/fasta/genome.fa" CHIKV.fa > genome.fa

# Combined gene annotation: human + CHIKV. The CHIKV GTF only annotates "gene"
# features, so edit_gtf_for_cellranger.py adds the transcript/exon features that
# cellranger mkref (STAR) requires.
cp "${CHIKVREF_DIR}/GCF_000854045.1_ViralProj14998_genomic.gtf.gz" CHIKV.gtf.gz
gunzip -f CHIKV.gtf.gz
"${SCRIPT_DIR}/edit_gtf_for_cellranger.py" CHIKV.gtf > CHIKV.exon.gtf
cat "${HSAREF_DIR}/genes/genes.gtf" CHIKV.exon.gtf > genes.gtf

cellranger mkref --genome=hsa_chikv \
  --fasta=genome.fa \
  --genes=genes.gtf \
  --ref-version=GRCh38-2020-A \
  --memgb="${LOCALVMEM}" \
  --nthreads="${LOCALCORES}"
