# R Scripts and Figure-Specific Data

This repository contains the R scripts and cleaned, figure-specific input datasets used to generate selected panels of the main-text figures for the study:

*Plastisphere formation reshapes the ecological toxicity of microplastics in* *Corbicula fluminea*.

The repository documents the code-based visualization workflows and the provenance of panels generated using external analysis platforms. It does not include raw sequencing data or code for Supplementary Figures.

## Repository organization

The repository is organized by main-text figure. Each `Figure_1/`–`Figure_7/` folder contains the R script(s), the corresponding input CSV file(s), and a folder-level README. Figure-specific code and data are kept together; there are no shared `data/` or `scripts/` directories.

| Figure folder | R scripts and panel coverage | Data and processing provenance |
|---|---|---|
| `Figure_1/` | `Fig1.R` — panels a–j | Uses the sample-level physiological and biochemical endpoint data supplied in this folder. |
| `Figure_2/` | `Fig2.R` — panel a | Panel b GSEA plots were exported from Majorbio Cloud; no R script for panel b is included. |
| `Figure_3/` | `Fig3.R` — panels a–b | Panels c–d were constructed from Majorbio-exported DEP data using STRING and visualized/exported with Cytoscape; no R script for these panels is included. |
| `Figure_4/` | `Fig4ab.R` — panels a–b; `Fig4c.R` — panel c | Panels a–b use the supplied nine-quadrant input tables. Panel c uses Pearson correlation coefficients and adjusted P values exported from Majorbio Cloud; the R script renders the heatmap. |
| `Figure_5/` | `Fig5.R` — panels a–c; `Fig5f.R` — panel f; `Fig5h-k.R` — panels h–k; `Fig5l.R` — panel l | Panels d–e use PCoA/NMDS results generated in Majorbio Cloud. Panel g is a Majorbio-generated co-occurrence network laid out and rendered in Gephi. The exported node and edge tables are supplied as inputs for the R calculations in panels h–l. |
| `Figure_6/` | `Fig6.R` — KEGG ReporterScore heatmap | ReporterScore values and pathway metadata were exported from Majorbio Cloud. The R script visualizes these values; it does not recalculate ReporterScores. |
| `Figure_7/` | `Fig7.R` — panels a–c | Panels a–c use Majorbio-exported association results for R-based visualization. Pearson, Mantel, and sample-level input tables for panels d–e are also supplied, but no R plotting script for panels d–e is included. |

## Running the scripts

Run scripts from the repository root unless the corresponding folder-level README specifies otherwise. The scripts use repository-relative paths and do not require machine-specific file locations.

For Figure 1, run the script in a fresh R session from the repository root:

```r
source("Figure_1/Fig1.R")
```

For the other R-generated panels, run:

```sh
Rscript Figure_2/Fig2.R
Rscript Figure_3/Fig3.R
Rscript Figure_4/Fig4ab.R
Rscript Figure_4/Fig4c.R
Rscript Figure_5/Fig5.R
Rscript Figure_5/Fig5f.R
Rscript Figure_5/Fig5h-k.R
Rscript Figure_5/Fig5l.R
Rscript Figure_6/Fig6.R
Rscript Figure_7/Fig7.R
```

Each figure folder README provides panel-specific input, package, and output details. Some scripts require Cairo graphics support; see the relevant folder README before running them.

## Input data and provenance

- Figure-specific CSV inputs are stored alongside their corresponding scripts.
- `data_manifest.csv` provides a file-level inventory, including source worksheet, panel coverage, dimensions, and the current committed file size in bytes.
- `DATA_PROVENANCE_AND_QA.md` describes the source data, data-preparation steps, and quality-control notes, including any item that remains pending source verification.
- CSV files use UTF-8 encoding (including UTF-8 BOM where applicable); scripts specify the relevant encoding when reading input files.

The supplied datasets are prepared for figure reproduction. They do not contain raw sequencing reads or replace the complete underlying transcriptomic, proteomic, 16S, or metagenomic datasets.
