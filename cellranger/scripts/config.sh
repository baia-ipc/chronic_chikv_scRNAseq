#!/bin/bash
#
# Shared configuration for the Cell Ranger reproduction pipeline.
#
# This file is *sourced* by the run_*.sh scripts. It only declares the inputs
# that live outside of this repository and that are therefore specific to the
# machine the pipeline is run on. Nothing here is hard-coded to a particular
# host: set the variables below in your environment before running, e.g.
#
#   export CHIKV_READS_DIR=/path/to/10x/reads
#   export CHIKV_TRANSCRIPTOME=/path/to/cellranger_ref_gex_hsa_chikv
#
# (or export them once in your shell profile).

set -euo pipefail

# Root directory of the 10x raw reads. It must contain one subdirectory per
# sequencing run (RunID), each holding a "<RunID>_10X_RawData_Outs" folder with
# the per-sample FASTQ directories. See the project metadata table for the
# RunID / DemultiplexedID / FlowcellID of each sample.
READS_DIR="${CHIKV_READS_DIR:?Set CHIKV_READS_DIR to the root directory of the 10x reads}"

# Cell Ranger reference transcriptome with the CHIKV genome added. This is the
# output of build_reference.sh (the "hsa_chikv" directory it produces).
TRANSCRIPTOME="${CHIKV_TRANSCRIPTOME:?Set CHIKV_TRANSCRIPTOME to the cellranger_ref_gex_hsa_chikv directory}"

# Local compute resources handed to cellranger. These reproduce the values used
# for the original run; override them to fit the machine you run on.
LOCALCORES="${CHIKV_LOCALCORES:-60}"
LOCALVMEM="${CHIKV_LOCALVMEM:-200}"
