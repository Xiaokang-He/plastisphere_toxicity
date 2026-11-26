# Main Fig. 1 — Physiological and Biochemical Endpoints

## Files

- `Fig1.R`: R script for the main-text Fig. 1.
- `Clam_Physiological_and_Biochemical_Endpoints.csv`: input dataset required by the script.

Keep the R script and the input CSV file directly in this folder. Do not place them in separate subfolders.

## Run

From the repository root in RStudio:

```r
source("Figure_1/Fig1.R")
```

The script reads the input CSV from this folder and writes generated figures and analysis outputs to `Figure_1/results/Fig1/`.

## Requirements

The script checks for its required R packages and installs missing packages when needed. Times New Roman is used for figure text; availability of this font may depend on the operating system.
