# Figure 1 — Physiological and Biochemical Endpoints

`Clam_Physiological_and_Biochemical_Endpoints.csv` contains the sample-level data for Fig. 1a–j: siphoning efficiency, SOD, CAT, GSH, GSH-Px, MDA, AChE, lysozyme, D-Lac and DAO.

The file includes `sample` and an R-compatible `group` field (CTR, MP1, MP10, MPP1, MPP10). The group code was derived from the original sample identifier; measured values are unchanged.

Run the existing script from the repository root:

```r
source("Figure_1/Fig1.R")
```
