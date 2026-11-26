# Figure 7 — Gut microbiota–host associations

## Panels a–c

`Fig7.R` draws the association sets described in the three tables below. Symbol size maps to the number of microbial KOs or host genes or proteins in the set. The association tables use Mantel r > 0.70 and BH-adjusted P < 0.05 for both exposures (a) or for one exposure with BH-adjusted P ≥ 0.05 in the other (b–c).

| Panel | Association table |
| --- | --- |
| a | `Fig7a_Shared_Associations.csv` |
| b | `Fig7b_Microplastic_Specific_Associations.csv` |
| c | `Fig7c_Plastisphere_Specific_Associations.csv` |

Run `Rscript Figure_7/Fig7.R` from the repository root with R's `grid` package and Cairo. The script writes `Fig7_Association_Sets_Compact.pdf` and `Fig7_Association_Sets_Compact.png` to `Figure_7/results/Fig7/`.

## Panels d–e

Panel d uses the plastisphere comparison and panel e uses the microplastic comparison. Each analysis includes the same three control tanks and three tanks at each concentration of the given exposure. The scripts align the two input tables by tank ID. Diversity and biomarker values are tank measurements, and each average-degree value comes from the network for that treatment.

| Panel | Microbiota input | Biomarker input |
| --- | --- | --- |
| d | `Fig7d_Gut_Microbiota_Metrics_Plastisphere.csv` | `Fig7d_Host_Biochemical_Biomarkers_Plastisphere.csv` |
| e | `Fig7e_Gut_Microbiota_Metrics_Microplastic.csv` | `Fig7e_Host_Biochemical_Biomarkers_Microplastic.csv` |

Ellipses show pairwise Pearson correlations among the nine biomarkers, with stars based on their P values. Line widths and colors show Mantel r and unadjusted permutation P values, respectively. Each Mantel test compares Bray–Curtis dissimilarities for one microbiota metric and Euclidean distances for one biomarker. Each pair uses 999 permutations with a fixed permutation seed.

Run from the repository root with `ggplot2`, `dplyr` and `vegan`.

```sh
Rscript Figure_7/Fig7d_Plastisphere.R
Rscript Figure_7/Fig7e_Microplastic.R
```

The scripts write `Fig7d_Plastisphere_ellipse_final.pdf` and `Fig7e_Microplastic_ellipse_final.pdf` to `Figure_7/`.
