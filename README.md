# CHIKV scRNA-seq analysis

This repository contains the single-cell RNA-seq analysis for our Chikungunya
longitudinal study. The main comparisons are:

- chronic vs non-chronic participants
- acute vs 6-month samples

The repository currently has two main directories:

- `cellranger/`: Cell Ranger output links and scripts.
- `analysis-1/`: a Seurat-based scRNA-seq analysis pipeline using the Cell
  Ranger outputs as input.

## Repository layout

```text
.
|-- cellranger/          # Cell Ranger analysis/output links
|-- metadata/            # sample metadata used by analysis-1
|-- analysis-1/          # Seurat-based scRNA-seq pipeline
|   |-- steps/           # numbered processing and visualization steps
|   |-- scripts/         # analysis-level helper scripts
|   `-- results/         # collected report/result links
|-- analysis-2/          # Scripts of Lien de Caluwe'
|-- scripts/             # project-level rendering helpers
|-- renv/                # renv activation
|-- DESCRIPTION          # R package dependencies
`-- renv.lock            # locked R package state
```

Each `analysis-1/steps/<step>/` directory contains:

- `scripts/`: source R Markdown scripts and optional batch shell wrappers
- `rundir/`: rendered HTML reports, logs, and knit intermediates
- `results/`: step outputs such as RDS objects, plots, and tables

## Cell Ranger

The `cellranger/` directory contains the reproducible Cell Ranger workflow that
produces the input used by `analysis-1`. It uses the project sample sheet at
`metadata/samples.tsv` and creates one Cell Ranger `outs/` directory per sample.
The downstream Seurat pipeline reads files exposed as
`cellranger/results/<SampleID>/filtered_feature_bc_matrix.h5`.

Run Cell Ranger commands from the `cellranger/` directory:

```bash
cd cellranger
```

Configure external input locations and run resources:

```bash
export CHIKV_READS_DIR=/path/to/10x/reads
export CHIKV_TRANSCRIPTOME=/path/to/hsa_chikv

# Required only when rebuilding the human+CHIKV reference:
export CHIKV_HSA_REFDATA=/path/to/refdata-gex-GRCh38-2020-A
export CHIKV_GENOME_DIR=/path/to/chikv_refseq_GCF_000854045

export CHIKV_LOCALCORES=60
export CHIKV_LOCALVMEM=200
```

Requirements are Cell Ranger, `samtools`, and Python 3 with `docopt`, `sh`, and
`loguru`. Cell Ranger v9.0.1 was used for this project.

Build the human+CHIKV reference once, if it is not already available:

```bash
./scripts/build_reference.sh
```

Run `cellranger count` for all samples. Already-computed samples are skipped:

```bash
./scripts/run_count.sh
```

Expose the generated `outs/` directories under `cellranger/results/`:

```bash
./scripts/link_results.sh
```

Collect per-sample Cell Ranger metrics:

```bash
./scripts/run_qc.sh
```

Optionally run per-sample CHIKV strand analysis:

```bash
./scripts/run_strand_analysis.sh
```

`run_count.sh` accepts selection and dry-run options forwarded to
`scripts/run_cellranger_count.py`, for example:

```bash
./scripts/run_count.sh --samples 037-A,037-6m
./scripts/run_count.sh --years 2025 --dry
```

Generated Cell Ranger outputs are written under `cellranger/reference/`,
`cellranger/rundir/`, `cellranger/results/`, `cellranger/qc_output/`, and
`cellranger/strand_analysis_out/`. See `cellranger/README.md` for the detailed
script reference.

## analysis-1 Seurat pipeline

Run commands from the repository root unless noted otherwise.

The pipeline uses `scripts/knit2html` to render R Markdown steps. When an Rmd is
inside a step `scripts/` directory, `knit2html` runs it from the sibling
`rundir/` directory and writes the linked HTML report there.

Restore the R environment before running analysis steps:

```bash
Rscript -e 'source("renv/activate.R"); renv::restore()'
```

### Step order

1. Create per-sample Seurat objects.
2. Render pre-filtering QC reports.
3. Filter per-sample Seurat objects.
4. Render merged post-filtering QC reports.
5. Integrate filtered objects with Harmony and post-process the integrated object.
6. Render integrated-object reports.
7. Run cell proportion tests and reports.
8. Run pseudobulk DESeq2 and reports.
9. Run pathway analysis and reports.
10. Run CellChat analysis and reports.

### Run the pipeline

Create Seurat objects for all configured samples:

```bash
cd analysis-1/steps/01.P.create_so/scripts
bash run_all_samples.sh
cd -
```

Render pre-filtering QC reports:

```bash
cd analysis-1/steps/01.V.unfiltered/scripts
bash run_all_samples.sh
../../../../scripts/knit2html before_filtering_merged.Rmd
cd -
```

Filter Seurat objects for all configured samples:

```bash
cd analysis-1/steps/02.P.filter/scripts
bash run_all_samples.sh
cd -
```

Render the merged post-filtering report:

```bash
scripts/knit2html analysis-1/steps/02.V.filtered/scripts/after_filtering_merged.Rmd
```

Run Harmony integration and post-processing:

```bash
scripts/knit2html analysis-1/steps/03.P.integration/scripts/harmony.Rmd
scripts/knit2html analysis-1/steps/03.P.integration/scripts/harmony_post_processing.Rmd
```

Render integrated-object reports:

```bash
scripts/knit2html analysis-1/steps/03.V.integrated/scripts/post_harmony_analysis.Rmd
scripts/knit2html analysis-1/steps/03.V.integrated/scripts/UMAP_single_cell_types.Rmd
```

Run cell proportion tests and report:

```bash
scripts/knit2html analysis-1/steps/04.P.prop_test/scripts/cell_proportions.Rmd
scripts/knit2html analysis-1/steps/04.V.prop_results/scripts/proportion_analysis.Rmd
```

Run pseudobulk DESeq2 and all DE reports:

```bash
scripts/knit2html analysis-1/steps/05.P.compute_de/scripts/compute_de.Rmd
bash analysis-1/steps/05.V.de_results/scripts/run_all_reports.sh
```

Run pathway computation and all pathway reports:

```bash
scripts/knit2html analysis-1/steps/06.P.compute_pathways/scripts/compute_pathway_analysis_GO.Rmd
scripts/knit2html analysis-1/steps/06.P.compute_pathways/scripts/compute_pathway_analysis_KEGG.Rmd
scripts/knit2html analysis-1/steps/06.P.compute_pathways/scripts/compute_pathway_analysis_MSigDB.Rmd
scripts/knit2html analysis-1/steps/06.P.compute_pathways/scripts/compute_pathway_analysis_Reactome.Rmd
bash analysis-1/steps/06.V.pathways_results/scripts/run_all_reports.sh
```

Run CellChat computation and all CellChat reports:

```bash
scripts/knit2html analysis-1/steps/07.P.compute_cellchat/scripts/compute_cellchat.Rmd
bash analysis-1/steps/07.V.cellchat_results/scripts/run_all_reports.sh
```

Collect selected result links under `analysis-1/results/`:

```bash
bash analysis-1/scripts/collect_results.sh
```

## Lien de Caluwe' analysis scripts

Lien de Caluwe' (Immunology Unit, Institut Pasteur du Cambodge)
first author of the submitted manuscript, performed some
furhter analyses. Her scripts are contained in analysis-2.
Please contact her for more information.

## Notes

- The per-sample wrappers define the current sample list internally.
- If GNU `parallel` is available, the per-sample wrappers use it; otherwise they
  run samples serially.
- Intermediate objects are written under each step's `results/vars/`.
- Linked HTML reports are written under each step's `rundir/`.

