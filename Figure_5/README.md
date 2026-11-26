# Main Fig. 5 — Gut Microbiota Diversity, Composition and Co-occurrence Networks

## Figure panels

The planned R script covers panels **a–l**:

- a–c: Alpha-diversity indices (Shannon, Pielou's evenness and ACE)
- d–e: Beta-diversity analyses (PCoA and NMDS), with ANOSIM
- f: Genus-level relative abundance
- g: Gut microbial co-occurrence networks
- h–k: Network topological properties (positive correlation ratio, number of nodes, number of edges and average degree)
- l: Network robustness evaluated by natural connectivity under random node removal

## Data provenance

The node and edge data used for panel g were calculated and exported from Majorbio Cloud. Figure-specific R code is used for the topology-related outputs, including the panel-i node-count plot, according to the current workflow. The exact input files and the division between imported network data, calculated metrics and plotted summaries will be documented after the complete script is checked.

## Files

- `Fig5.R`: R script for panels a–l.
- `Gut_16S_Alpha_Diversity.csv`: input dataset currently specified in the existing script.
- Additional figure-specific input files will be stored directly in this folder.

## Run

From the repository root in RStudio:

```r
source("Figure_5/Fig5.R")
```

Generated figures and analysis outputs are written to `Figure_5/results/Fig5/`.

## Requirements

The script checks for its required R packages and installs missing packages when needed. Times New Roman is used for figure text; font availability may depend on the operating system.
