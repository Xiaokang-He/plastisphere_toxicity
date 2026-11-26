# Figure 5 — Gut Microbiota

This folder contains the R scripts and input tables required to reproduce the panels generated with R.

## File inventory

| File | Purpose | Used by |
|---|---|---|
| `Fig5.R` | Gut 16S alpha-diversity panels (a–c) | `Fig5a-c_Gut_16S_Alpha_Diversity.csv` |
| `Fig5f.R` | Genus-level relative-abundance panel (f) | `Fig5f_Genus_Relative_Abundance.csv` |
| `Fig5h-k.R` | Network topology panels (h–k) | `Fig5g_Network_Topology_Raw_Nodes.csv`, `Fig5g_Network_Topology_Raw_Edges.csv` |
| `Fig5l.R` | Network robustness panel (l), based on natural connectivity after random node removal | `Fig5g_Network_Topology_Raw_Nodes.csv`, `Fig5g_Network_Topology_Raw_Edges.csv` |

### Input data

- `Fig5a-c_Gut_16S_Alpha_Diversity.csv`: alpha-diversity measurements for panels a–c.
- `Fig5f_Genus_Relative_Abundance.csv`: genus-level relative abundances for panel f.
- `Fig5g_Network_Topology_Raw_Nodes.csv`: exported node table used for panels h–l.
- `Fig5g_Network_Topology_Raw_Edges.csv`: exported edge table used for panels h–l.

All four CSV files are used by at least one script and should be retained.

## Reproducing the R-generated panels

Run the commands below from the repository root.

### Panels a–c: Alpha diversity

Run with base R graphics support.

```sh
Rscript Figure_5/Fig5.R
```

Output: `Figure_5/results/Fig5/Fig5abc_Alpha_Diversity_test.pdf`.

### Panel f: Genus-level relative abundance

Required packages: `ggplot2`, `ggalluvial`, and Cairo graphics support.

```sh
Rscript Figure_5/Fig5f.R
```

Outputs: `Figure_5/results/Fig5f.pdf` and `Figure_5/results/Fig5f.png`.

### Panels h–k: Network topology

Required packages: `dplyr`, `ggplot2`, `ggpattern`, `cowplot`, and Cairo graphics support.

The script uses the supplied node and edge tables to calculate the number of distinct nodes, the number of unique undirected edges, the positive-correlation ratio, and the average degree (2E/N). It does not recalculate network correlations or reconstruct the co-occurrence network.

```sh
Rscript Figure_5/Fig5h-k.R
```

Output: `Figure_5/results/Fig5h-k_Network_Topology.pdf`.

**Package compatibility:** `ggpattern` requires a compatible `ggplot2` version. If it fails to load with a namespace/version error, update `ggplot2` to meet the installed `ggpattern` dependency and reinstall `ggpattern`.

### Panel l: Network robustness

Required package: `ggplot2` (the script also uses R's built-in `grid` package).

The script builds an unweighted, undirected adjacency matrix from the supplied node and edge tables, randomly removes nodes, and calculates natural connectivity and average degree. The plotted panel shows natural connectivity as a function of the node-removal ratio.

```sh
Rscript Figure_5/Fig5l.R
```

Outputs: `Figure_5/Robustness.pdf` and `Figure_5/Robustness_data.csv`.

## Analyses performed outside R

The PCoA and NMDS analyses for panels d–e were performed in Majorbio Cloud using the exported coordinates. Majorbio Cloud generated the co-occurrence network for panel g and exported its nodes and edges; Gephi was used to lay out and draw that network.
