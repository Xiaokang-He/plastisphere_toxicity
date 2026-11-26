# Figure 7 — Gut Microbiota–Host Associations

## Panels a–c: association networks

- `Fig7a_Shared_Associations.csv`: panel a shared associations.
- `Fig7b_Microplastic_Specific_Associations.csv`: panel b microplastic-specific associations.
- `Fig7c_Plastisphere_Specific_Associations.csv`: panel c plastisphere-specific associations.

Run `Fig7.R` from the repository root or the `Figure_7` folder. The script reads the three association CSV files in this folder and saves the PDF and PNG outputs to `Figure_7/results/Fig7/` (or `results/Fig7/` when run from this folder).

## Panels d–e: microbiota–biomarker associations

Panel d represents the plastisphere exposure comparison, and panel e represents the microplastic exposure comparison. Each panel uses matched sample-level inputs for the three groups (Control, Low and High; three samples per group).

### Panel d — Plastisphere

- `Fig7d_Gut_Microbiota_Metrics_Plastisphere.csv`: Shannon, Pielou, ACE and average-degree metrics.
- `Fig7d_Host_Biochemical_Biomarkers_Plastisphere.csv`: host biochemical biomarkers.
- `Fig7d_Plastisphere.R`: plotting and statistical-analysis script.

### Panel e — Microplastic

- `Fig7e_Gut_Microbiota_Metrics_Microplastic.csv`: Shannon, Pielou, ACE and average-degree metrics.
- `Fig7e_Host_Biochemical_Biomarkers_Microplastic.csv`: host biochemical biomarkers.
- `Fig7e_Microplastic.R`: plotting and statistical-analysis script.

The scripts calculate Pearson correlations among the nine host biochemical biomarkers and pairwise Mantel tests between each gut microbiota metric and each host biomarker. Mantel tests use 999 permutations. The scripts retain the plotting settings in the supplied source code and export separate PDF files for panels d and e.

Run from the repository root:

```sh
Rscript Figure_7/Fig7d_Plastisphere.R
Rscript Figure_7/Fig7e_Microplastic.R
```

Alternatively, change the working directory to `Figure_7/` and run the corresponding `.R` file directly. Required packages are `ggplot2`, `dplyr` and `vegan`.