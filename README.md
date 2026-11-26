# R Scripts and Figure-Specific Data

This repository contains figure-specific R scripts and cleaned input datasets for the study:

*Plastisphere formation reshapes the ecological toxicity of microplastics in Corbicula fluminea.*

## Organization

Each main-text figure has its own folder. The R script and its figure-specific input data are stored directly in that folder. CSV files are UTF-8 encoded and can be read in R with:

```r
read.csv("Figure_1/Clam_Physiological_and_Biochemical_Endpoints.csv",
         fileEncoding = "UTF-8-BOM", check.names = FALSE)
```

## Figure coverage and data provenance

| Figure | R-code coverage planned/current | Data and external processing |
|---|---|---|
| Fig. 1 | a–j; `Figure_1/Fig1.R` is present | Sample-level siphoning and biochemical endpoints. |
| Fig. 2 | a; `Figure_2/Fig2.R` is present | Panel b GSEA plots were exported from Majorbio Cloud. |
| Fig. 3 | a–b; script to be added | Panels c–d were constructed from Majorbio-exported DEP data using STRING and visualized/exported with Cytoscape; no R script is supplied for those panels. |
| Fig. 4 | a–c; script to be added | Panel c correlation coefficients and adjusted P values were calculated/exported by Majorbio Cloud; code renders the heatmap. |
| Fig. 5 | a–l data are available; current `Figure_5/Fig5.R` covers alpha-diversity panels a–c | Panel g node and edge data were exported from Majorbio Cloud. Topology summaries and robustness data are supplied separately. |
| Fig. 6 | heatmap; script to be added | ReporterScore values were calculated/exported by Majorbio Cloud; code renders the heatmap. |
| Fig. 7 | a–e; script to be added | Panels a–c contain Majorbio-exported broad association results for code-based visualization. Panels d–e include separate Pearson and Mantel result tables. |

## Existing scripts

- `Figure_1/Fig1.R`: physiological and biochemical endpoints, panels a–j.
- `Figure_2/Fig2.R`: GO and KEGG enrichment bubble plot, panel a.
- `Figure_5/Fig5.R`: gut 16S alpha diversity, panels a–c. The additional Fig. 5 input tables are provided for subsequent code completion.

## Input data

The cleaned datasets are stored in the corresponding `Figure_1/` to `Figure_7/` folders. `data_manifest.csv` lists each data file, source worksheet, panel coverage and row count. `DATA_PROVENANCE_AND_QA.md` records the cleaning rules and data-quality notes.

For Fig. 1 and Fig. 5 sample-level datasets, a `group` field was added from the sample ID for compatibility with the existing scripts: Control = CTR, Microplastic_Low = MP1, Microplastic_High = MP10, Plastisphere_Low = MPP1, and Plastisphere_High = MPP10. Original sample identifiers and measured values are retained.

## Source-data note

The cleaned files were prepared against the latest matching source workbook available in the project Library during this task: `Source Data to Main Figures-DESKTOP-JPDUFHL(4).xlsx`. Numeric values were not recalculated. Where the separately uploaded CSVs had reduced precision, blank rows/columns, or corrupted characters, the source-workbook values were retained.

One taxonomy value requires manual verification: in Fig. 5g, the `Order` field for `Plastisphere-Low / OTU1656` is stored as numeric `46350` in the source workbook (the uploaded CSV displayed `11月24日`). It has not been guessed or silently corrected.
