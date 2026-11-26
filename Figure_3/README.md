# Figure 3 — Proteomic enrichment and interaction networks

`Fig3.R` generates the GO and KEGG enrichment plots for microplastic versus control (a) and plastisphere versus control (b). It reads `Fig3a_DEP_Functional_Enrichment.csv` and `Fig3b_DEP_Functional_Enrichment.csv`. The bars show −log₁₀ of the supplied adjusted P values, and the circles show the protein counts.

Run from the repository root with `ggplot2` and Cairo graphics support.

```sh
Rscript Figure_3/Fig3.R
```

The script writes `Fig3a.pdf`, `Fig3b.pdf` and `Fig3ab.pdf` to `Figure_3/results/`. Panels c–d show STRING protein–protein interaction networks drawn in Cytoscape using DEP results exported from Majorbio Cloud.
