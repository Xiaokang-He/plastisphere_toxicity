# R Scripts and Figure-Specific Data

This repository contains the R scripts and input datasets used to generate the main-text figures for the study:

*Plastisphere formation reshapes the ecological toxicity of microplastics in Corbicula fluminea.*

## Organization

Each main-text figure has its own folder. The R script and all input data required for that figure are stored together directly inside the same folder. No separate `scripts/` or `data/` directories are used, and datasets are not shared between figure folders.

```text
Figure_1/
  Fig1.R
  <figure-specific input data files>
  README.md
Figure_2/
  Fig2.R
  <figure-specific input data files>
  README.md
Figure_3/
  Fig3.R
  <figure-specific input data files>
  README.md
Figure_4/
  Fig4.R
  <figure-specific input data files>
  README.md
Figure_5/
  Fig5.R
  <figure-specific input data files>
  README.md
Figure_6/
  Fig6.R
  <figure-specific input data files>
  README.md
Figure_7/
  Fig7.R
  <figure-specific input data files>
  README.md
Plastisphere_Microplastic_Toxicity.Rproj
README.md
.gitignore
```

The folder names follow the main-text figure numbers. Each folder is self-contained and includes only the code and input data needed for its corresponding figure. The `results/` subfolder is created locally when a script is run; generated outputs are excluded from version control.

## Current status

| Figure | Analysis | R script | Input data |
|---|---|---|---|
| Fig. 1 | Physiological and biochemical endpoints | `Figure_1/Fig1.R` | `Clam_Physiological_and_Biochemical_Endpoints.csv` |
| Fig. 2 | To be added | `Figure_2/Fig2.R` | To be added |
| Fig. 3 | To be added | `Figure_3/Fig3.R` | To be added |
| Fig. 4 | To be added | `Figure_4/Fig4.R` | To be added |
| Fig. 5 | Gut 16S alpha-diversity indices | `Figure_5/Fig5.R` | `Gut_16S_Alpha_Diversity.csv` |
| Fig. 6 | To be added | `Figure_6/Fig6.R` | To be added |
| Fig. 7 | To be added | `Figure_7/Fig7.R` | To be added |

## Requirements

- R and RStudio.
- The R packages required by each figure-specific script.
- The input data files listed in the corresponding figure folder's README.

Each script should be run from the repository root so that its relative paths resolve correctly. For example:

```r
source("Figure_1/Fig1.R")
```

The script reads its input data from `Figure_1/` and writes generated figures and analysis outputs to `Figure_1/results/Fig1/`. Other figure scripts follow the same pattern.

## Data and code

Input datasets are stored alongside their corresponding R scripts. Each figure folder is independent; no shared input-data directory is used.

Repository: https://github.com/Xiaokang-He/plastisphere_toxicity
