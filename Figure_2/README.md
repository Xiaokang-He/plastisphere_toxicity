# Figure 2 — Transcriptomic functional enrichment

`Fig2.R` generates the panel-a bubble plot from `Fig2a_DEG_Functional_Enrichment.csv`. Each point represents an enriched term in one exposure–control comparison. Point size maps to gene count and fill maps to the supplied BH-adjusted P value.

Run from the repository root with `ggplot2`, `dplyr`, `readr`, `ggh4x` and Cairo graphics support.

```sh
Rscript Figure_2/Fig2.R
```

The script writes `Fig2a.pdf` and `Fig2a.png` to `Figure_2/results/`. The GSEA plots in panel b were exported from Majorbio Cloud.
