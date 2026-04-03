For the moment the output directory of the chikv analysis-2 and analysis-4 steps running cellranger
are linked here as ``prev_analysis_1`` and ``prev_analysis_2`` (output directories).

```
prev_analysis_1 -> /srv/immunology/chikv/analysis-3/steps/002.0.cellranger_count.chikv_hsa_ref/output/
prev_analysis_2 -> /srv/immunology/chikv/analysis-4/steps/002.0.cellranger_count.chikv_hsa_ref/output/
```

The single samples outputs of cellranger are then linked here as well.
The commands in ``make_links.sh`` were used to create the links
