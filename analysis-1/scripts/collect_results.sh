#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${ANALYSIS_DIR}/results"

mkdir -p "${RESULTS_DIR}"
cd "${RESULTS_DIR}"

ln -s ../steps/02.V.filtered/rundir/after_filtering_merged.html 1_preprocessing.html
ln -s ../steps/03.V.integrated/rundir/post_harmony_analysis.html 2_after_integration.html

mkdir -p 3_proportion_tests
cd 3_proportion_tests
ln -s ../../steps/04.V.prop_results/results/plots
ln -s ../../steps/04.V.prop_results/results/tables
ln -s ../../steps/04.V.prop_results/rundir/proportion_analysis.html report.html
cd ..

mkdir -p 4_deseq
cd 4_deseq
ln -s ../../steps/05.V.de_results/results/tables/p_adj_filt p_adj_filtered_tables
ln -s ../../steps/05.V.de_results/results/plots/vulcano vulcano_plots
ln -s ../../steps/05.V.de_results/rundir/deseq_tables_and_vulcanos.html report.html
cd ..

mkdir -p 5_pathways
cd 5_pathways
ln -s ../../steps/06.V.pathways_results/results/* .
