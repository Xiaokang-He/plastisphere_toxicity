# Figure 1 — Physiological and Biochemical Endpoints

`Fig1a_Clam_Physiological_and_Biochemical_Endpoints.csv` contains 15 sample-level records for Fig. 1a–j: siphoning efficiency, SOD, CAT, GSH, GSH-Px, MDA, AChE, lysozyme, D-Lac and DAO.

The group is derived from the `sample` identifier when the CSV does not contain a separate `group` column. The plotting script renames `lysozyme` to `LZM` and `D_Lac` to its internal plotting name.

## Reproducing Fig. 1

Run from the repository root:

```r
source("Figure_1/Fig1.R")
```

The script also works when the working directory is already `Figure_1`:

```r
source("Fig1.R")
```

Required packages: `ggplot2`, `ggpattern` and `WRS2`. The script uses `lincon()` for the ten pairwise comparisons among the five groups for each endpoint, with `tr = 0.2` and no multiple-comparison adjustment. Significant raw P values are calculated from the supplied CSV data and added to the figure.

The script exports only:

- `Fig1.pdf`: vector PDF.
- `Fig1.jpg`: 300-dpi JPEG.

Both files are saved alongside the input CSV in `Figure_1/`. No machine-specific absolute paths or `setwd()` commands are used.
