# Main Fig. 6 — Gut Microbial KEGG Pathway Enrichment

## Data provenance

ReporterScore values and the associated pathway-level data were calculated and exported from Majorbio Cloud. The R script uses these exported values to generate the heatmap. It does not recalculate the ReporterScores.

## Files

- `Fig6.R`: R script for heatmap visualization (to be added).
- The Majorbio-exported ReporterScore table and any other figure-specific input files will be stored directly in this folder.

## Run

After the script and input files have been added, run from the repository root in RStudio:

```r
source("Figure_6/Fig6.R")
```
