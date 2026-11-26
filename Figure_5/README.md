# Figure 5 — Gut Microbiota

- `Gut_16S_Alpha_Diversity.csv`: panels a–c; input filename and columns match the current `Fig5.R`.
- `Fig5f_Genus_Relative_Abundance.csv`: genus-level relative abundance.
- `Network_Topology_Raw_Nodes.csv` and `Network_Topology_Raw_Edges.csv` contain the exported node taxonomy and edge lists used as inputs for network topology calculations.

## Ordination and network provenance

The PCoA and NMDS analyses for panels d–e were performed in Majorbio Cloud, and the exported coordinates were used directly. Majorbio Cloud computed the co-occurrence network for panel g and exported its nodes and edges. Gephi was used to lay out and draw the exported network.

## Reproduce Fig. 5f

Run the following command from the repository root with `ggplot2`, `ggalluvial`, and Cairo graphics support.

```sh
Rscript Figure_5/Fig5f.R
```

The script reads the supplied relative abundances as percentages and saves `Fig5f.pdf` and a 600 dpi `Fig5f.png` in `Figure_5/results/`.


## Reproduce Fig. 5h-k

Run the following command from the repository root with `dplyr`, `ggplot2`, `ggpattern`, `cowplot`, and Cairo graphics support.

```sh
Rscript Figure_5/Fig5h-k.R
```

The script uses the supplied Majorbio-exported node and edge tables to calculate the number of distinct nodes, the number of unique undirected edges, the positive-correlation ratio, and the average degree (2E/N) for each group. It does not recalculate network correlations or reconstruct the co-occurrence network. The combined figure is saved as `Figure_5/results/Fig5h-k_Network_Topology.pdf`.
