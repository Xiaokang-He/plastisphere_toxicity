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
| Fig. 3 | a–b; `Figure_3/Fig3.R` is present | Panels c–d were constructed from Majorbio-exported DEP data using STRING and visualized/exported with Cytoscape; no R script is supplied for those panels. |
| Fig. 4 | c; `Figure_4/Fig4c.R` is present | The heatmap uses Pearson correlation coefficients and adjusted P values exported from Majorbio Cloud. |
| Fig. 5 | a–c, f and h–k; `Figure_5/Fig5.R`, `Figure_5/Fig5f.R` and `Figure_5/Fig5h-k.R` are present | Majorbio Cloud supplied the PCoA and NMDS coordinates for panels d–e and the co-occurrence network for panel g. Gephi was used for the network layout and rendering. Panels h–k calculate network topology metrics from the exported node and edge tables; network correlations are not recalculated. |
| Fig. 6 | heatmap; script to be added | ReporterScore values were calculated/exported by Majorbio Cloud; code renders the heatmap. |
| Fig. 7 | a–e; sample-level inputs and exported association-result tables are present; script to be added | Panels a–c contain Majorbio-exported broad association results for code-based visualization. Panels d–e use sample-level gut microbiota metrics and host biochemical biomarkers, alongside separate Pearson and Mantel result tables. |

## Existing scripts

- `Figure_1/Fig1.R`: physiological and biochemical endpoints, panels a–j.
- `Figure_2/Fig2.R`: GO and KEGG enrichment bubble plot, panel a.
- `Figure_3/Fig3.R`: proteomic functional enrichment, panels a–b.
- `Figure_4/Fig4c.R`: transcriptome-proteome Pearson correlation heatmap, panel c.
- `Figure_5/Fig5.R`: gut 16S alpha diversity, panels a–c.
- `Figure_5/Fig5f.R` plots genus-level relative abundance for panel f.\n- `Figure_5/Fig5h-k.R` calculates and plots network topology metrics for panels h–k from the supplied node and edge tables.

## Input data

The cleaned datasets are stored in the corresponding `Figure_1/` to `Figure_7/` folders. `data_manifest.csv` lists each data file, source worksheet, panel coverage and row count. `DATA_PROVENANCE_AND_QA.md` records the cleaning rules and data-quality notes.

For Fig. 1 and Fig. 5 sample-level datasets, a `group` field was added from the sample ID for compatibility with the existing scripts: Control = CTR, Microplastic_Low = MP1, Microplastic_High = MP10, Plastisphere_Low = MPP1, and Plastisphere_High = MPP10. Original sample identifiers and measured values are retained.

## Source-data note

The cleaned files were prepared against the latest matching source workbook available in the project Library during this task: `Source Data to Main Figures-DESKTOP-JPDUFHL(4).xlsx`. Numeric values were not recalculated. Where the separately uploaded CSVs had reduced precision, blank rows/columns, or corrupted characters, the source-workbook values were retained.

The `Order` field for `Plastisphere-Low / OTU1656` in `Figure_5/Network_Topology_Raw_Nodes.csv` requires manual verification. The table retains the numeric value `46350` from the source workbook (the uploaded CSV displayed `11月24日`).
