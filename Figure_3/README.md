# Figure 3 — Proteomic Enrichment and PPI Networks

`Fig3.R` reproduces the enrichment plots for panels a and b. From the repository root, run `Rscript Figure_3/Fig3.R` with `ggplot2` and Cairo support in R. The script writes `Fig3a.pdf`, `Fig3b.pdf`, and `Fig3ab.pdf` to `Figure_3/results/`.

- `Fig3a_DEP_Functional_Enrichment.csv`: Fig. 3a DEP enrichment results for the microplastic comparison.
- `Fig3b_DEP_Functional_Enrichment.csv`: Fig. 3b DEP enrichment results for the plastisphere comparison.

Panels c–d were constructed by submitting Majorbio-exported DEP data to STRING and visualizing/exporting the PPI networks with Cytoscape. No R script is supplied for those panels.
