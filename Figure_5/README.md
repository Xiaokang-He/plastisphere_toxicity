# Figure 5 — Gut microbiota

Run all commands from the repository root.

| Panels | Command | Input CSVs | Output |
| --- | --- | --- | --- |
| a–c | `Rscript Figure_5/Fig5.R` | `Fig5a-c_Gut_16S_Alpha_Diversity.csv` | `Figure_5/results/Fig5/Fig5abc_Alpha_Diversity_test.pdf` |
| f | `Rscript Figure_5/Fig5f.R` | `Fig5f_Genus_Relative_Abundance.csv` | `Fig5f.pdf` and `Fig5f.png` in `Figure_5/results/` |
| h–k | `Rscript Figure_5/Fig5h-k.R` | `Fig5g_Network_Topology_Raw_Nodes.csv` and `Fig5g_Network_Topology_Raw_Edges.csv` | `Figure_5/results/Fig5h-k_Network_Topology.pdf` |
| l | `Rscript Figure_5/Fig5l.R` | `Fig5l_Network_Robustness.csv` | `Robustness.pdf` and `Robustness_data.csv` in `Figure_5/` |

**Panels a–c.** Shannon diversity, Pielou's evenness and ACE richness are read by tank. The bars show means ± SD, and the brackets show pairwise P values below 0.05. With three tanks per group and 20% trimming, the unadjusted pairwise tests reduce to Welch comparisons. The script uses R's `grid` package and Cairo.

**Panel f.** The five relative-abundance profiles use the percentages in the input table. The script uses `ggplot2`, `ggalluvial` and Cairo.

**Panels h–k.** The script derives the percentage of positive edges, number of distinct nodes, number of unique undirected edges and average degree (2E/N), in that order. It uses `dplyr`, `ggplot2`, `ggpattern`, `cowplot` and Cairo.

**Panel l.** The input table gives natural connectivity for each of the five networks across 0–99% random node removal. The script plots those values, adds a smooth curve for each treatment, and exports the plotted points in `Robustness_data.csv`. It uses `ggplot2` and R's `grid` package.

Panels d–e show Bray–Curtis PCoA and NMDS from Majorbio Cloud. Panel g shows the OTU co-occurrence networks laid out in Gephi, with nodes colored by phylum. The exported node and edge tables also supply the calculations in panels h–k.
