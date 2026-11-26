# Figure 5 — Gut Microbiota

- `Gut_16S_Alpha_Diversity.csv`: panels a–c; input filename and columns match the current `Fig5.R`.
- `Fig5d-e_Ordination_Coordinates.csv`: PCoA and NMDS sample coordinates.
- `Fig5d_PCoA_Explained_Variance.csv`: PCoA explained variance.
- `Fig5d_ANOSIM_Results.csv`: ANOSIM summary.
- `Fig5e_NMDS_Stress.csv`: NMDS stress.
- `Fig5f_Genus_Relative_Abundance.csv`: genus-level relative abundance.
- `Fig5g_Network_Nodes.csv` and `Fig5g_Network_Edges.csv`: separate network node taxonomy and edge tables.
- `Fig5h-k_Network_Topology.csv`: network topology summaries.
- `Fig5l_Network_Robustness.csv`: natural-connectivity robustness data.

Panel-g nodes and edges were exported from Majorbio Cloud. Panels h–k and l were extracted from the source workbook.

## Reproduce Fig. 5f

Run the following command from the repository root with `ggplot2`, `ggalluvial`, and Cairo graphics support.

```sh
Rscript Figure_5/Fig5f.R
```

The script reads the supplied relative abundances as percentages and saves `Fig5f.pdf` and a 600 dpi `Fig5f.png` in `Figure_5/results/`.
