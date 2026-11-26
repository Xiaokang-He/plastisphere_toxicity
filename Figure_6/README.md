# Figure 6 — Gut Microbial KEGG ReporterScore

`Fig6_KEGG_ReporterScore.csv` contains the ReporterScore values and pathway metadata exported from Majorbio Cloud. The R script uses these values to render the heatmap; ReporterScores are not recalculated here.

## Reproducing the heatmap

Run the script from the repository root:

```sh
Rscript Figure_6/Fig6.R
```

Required package: `ggplot2`. Cairo graphics support is required for the PDF export.

The script saves the following files to `Figure_6/results/Fig6/`:

- `Fig6_KEGG_ReporterScore.pdf`
- `Fig6_KEGG_ReporterScore.png` (600 dpi)
- `Fig6_KEGG_ReporterScore.tiff` (600 dpi, LZW compression)
