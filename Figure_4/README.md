# Figure 4 — Transcriptome–proteome integration

## Panels a–b

`Fig4ab.R` reads `Fig4a_Nine_Quadrant_Microplastic.csv` and `Fig4b_Nine_Quadrant_Plastisphere.csv`. The points, differential-expression category counts and two-sided Pearson tests use every matched gene–protein pair in each comparison. Both plots share the same axes.

Run from the repository root with `ggplot2` and `patchwork`.

```sh
Rscript Figure_4/Fig4ab.R
```

The script writes `Fig4ab.pdf` and `Fig4ab.png` to `Figure_4/results/`.

## Panel c

`Fig4c.R` reads `Fig4c_Transcriptome_Proteome_Correlation.csv`. Each tile shows a term-level Pearson correlation between transcript and protein log₂ fold changes. Asterisks mark BH-adjusted P values below 0.05. Term labels include matched gene–protein pair counts.

Run from the repository root with `ggplot2` and Cairo graphics support.

```sh
Rscript Figure_4/Fig4c.R
```

The script writes `Fig4c.pdf` and `Fig4c.png` to `Figure_4/results/`.
