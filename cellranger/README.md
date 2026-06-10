# Cell Ranger

This directory holds the Cell Ranger inputs for the downstream analysis: one
directory per sample under `results/<SampleID>/`, each being a Cell Ranger
`outs/` directory (the analysis reads `results/<SampleID>/filtered_feature_bc_matrix.h5`
etc. directly).

The samples and their sequencing metadata are defined in the project sample
sheet [`../metadata/samples.tsv`](../metadata/samples.tsv).

## Run

The scripts in `scripts/` allow to run the analysis.
External inputs are supplied through environment variables (see
[`scripts/config.sh`](scripts/config.sh)).

Requirements: `cellranger` (we used v.9.0.1), `samtools`, and Python 3 with
`docopt`, `sh` and `loguru`.

Configure the external inputs:

```bash
export CHIKV_READS_DIR=/path/to/10x/reads      # one <RunID>/ subdir per run
export CHIKV_TRANSCRIPTOME=/path/to/hsa_chikv  # reference to use

# to rebuild the reference (GRCh38 + CHIKV):
export CHIKV_HSA_REFDATA=/path/to/refdata-gex-GRCh38-2020-A
export CHIKV_GENOME_DIR=/path/to/chikv_refseq_GCF_000854045

# run parameters:
export CHIKV_LOCALCORES=60 CHIKV_LOCALVMEM=200
```

Then, from this directory:

```bash
# 1. (necessary once only) build the human+CHIKV reference transcriptome -> reference/hsa_chikv
./scripts/build_reference.sh

# 2. run cellranger count for every sample -> rundir/<SampleID>/outs
#    (already-computed samples are skipped automatically)
./scripts/run_count.sh

# 3. expose the outs/ directories as results/<SampleID>
./scripts/link_results.sh

# 4. QC: collect the per-sample metrics into qc_output/metrics_summaries.tsv
./scripts/run_qc.sh

# 5. (optional) per-sample CHIKV strand analysis -> strand_analysis_out/
./scripts/run_strand_analysis.sh
```

`run_count.sh` accepts extra selection options forwarded to
`run_cellranger_count.py` (e.g. `--samples 037-A,037-6m`, `--years 2025`,
`--dry`); see `./scripts/run_cellranger_count.py --help`.

### Scripts

| Script | Purpose |
| --- | --- |
| `config.sh` | Shared configuration; declares the external inputs (sourced by the run scripts). |
| `build_reference.sh` | `cellranger mkref` for the human+CHIKV reference. |
| `edit_gtf_for_cellranger.py` | Adds transcript/exon features to the CHIKV GTF for STAR. |
| `run_cellranger_count.py` | Runs `cellranger count` for the samples in the metadata table. |
| `run_count.sh` | Wrapper running the count step over `../metadata/samples.tsv`. |
| `link_results.sh` | Links each `rundir/<SampleID>/outs` to `results/<SampleID>`. |
| `collect_metric_summaries.py` | Merges the per-sample `metrics_summary.csv` files. |
| `run_qc.sh` | Wrapper producing `qc_output/metrics_summaries.tsv`. |
| `run_strand_analysis.sh` | Per-cell CHIKV UMI strand analysis over all samples. |
| `strand_analysis.py` | SAM/BAM strand counting used by `run_strand_analysis.sh`. |

The generated directories (`reference/`, `rundir/`, `results/`, `qc_output/`,
`strand_analysis_out/`) are git-ignored.
