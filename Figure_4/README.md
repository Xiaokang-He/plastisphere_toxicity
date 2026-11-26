# Figure 4 — Multi-omics Integration

- `Fig4a_Nine_Quadrant_Microplastic.csv`: panel a nine-quadrant input.
- `Fig4b_Nine_Quadrant_Plastisphere.csv`: panel b nine-quadrant input.
- `Fig4c_Transcriptome_Proteome_Correlation.csv`: Majorbio-exported correlation coefficients and adjusted P values used to render the panel-c heatmap.

The numeric values in panel c are imported results; the R workflow is for heatmap visualization.

## Reproduce Fig. 4c

Run the following command from the repository root with ggplot2 installed (tested with version 3.5.2) and Cairo graphics support.

```sh
Rscript Figure_4/Fig4c.R
```

The script reads the supplied correlation table and marks adjusted P values below 0.05 with an asterisk. It saves `Fig4c.pdf` and a 600 dpi `Fig4c.png` in `Figure_4/results/`.
