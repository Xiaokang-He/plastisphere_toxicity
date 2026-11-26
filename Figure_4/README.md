# Main Fig. 4 — Multi-omics Integration

## Figure panels and provenance

- **Panels a–b:** Nine-quadrant plots showing gene–protein fold-change relationships for the microplastic-versus-control and plastisphere-versus-control comparisons. R scripts will generate these plots.
- **Panel c:** Pearson correlation coefficients were calculated and exported from Majorbio Cloud. The R script uses the exported correlation results to render the heatmap; it does not recalculate the platform-derived correlation coefficients.
- **Panels d–e:** Schematic illustrations of the proposed subcellular toxicity pathways. These panels are not included in the planned R scripts.

## Files

- `Fig4.R`: R script for panels a–c (to be added).
- Figure-specific gene–protein data and the Majorbio-exported correlation table required by the script will be stored directly in this folder.

## Run

After the script and input files have been added, run from the repository root in RStudio:

```r
source("Figure_4/Fig4.R")
```
