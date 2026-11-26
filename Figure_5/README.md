# Main Fig. 5 — Gut 16S Alpha Diversity

## Files

- `Fig5.R`: R script for the main-text Fig. 5.
- `Gut_16S_Alpha_Diversity.csv`: input dataset required by the script.

Keep the R script and the input CSV file directly in this folder. Do not place them in separate subfolders.

## Run

From the repository root in RStudio:

```r
source("Figure_5/Fig5.R")
```

The script reads the input CSV from this folder and writes generated figures and analysis outputs to `Figure_5/results/`.

## Requirements

The script checks for its required R packages and installs missing packages when needed. Times New Roman is used for figure text; availability of this font may depend on the operating system.
