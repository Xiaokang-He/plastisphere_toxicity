# Figure 6 — Gut microbial KEGG ReporterScore

`Fig6.R` reads one row per exposure–control comparison and KEGG pathway from `Fig6_KEGG_ReporterScore.csv`. Tile color maps to the supplied ReporterScore. The tile label gives the total number of KOs followed in parentheses by the count with higher mean abundance under exposure. An asterisk marks |ReporterScore| ≥ 1.65.

Run from the repository root with `ggplot2` and Cairo graphics support.

```sh
Rscript Figure_6/Fig6.R
```

The script writes `Fig6_KEGG_ReporterScore.pdf`, `Fig6_KEGG_ReporterScore.png` and `Fig6_KEGG_ReporterScore.tiff` to `Figure_6/results/Fig6/`. The raster files use 600 dpi.
