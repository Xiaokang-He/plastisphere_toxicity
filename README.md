# Analysis Scripts

This repository contains R scripts used for data analysis and figure generation for the study:

*Plastisphere formation reshapes the ecological toxicity of microplastics in Corbicula fluminea.*

## Included scripts

| Figure | Analysis | R script | Proposed input data name |
|---|---|---|---|
| Fig. 1 | Physiological and biochemical endpoints | `scripts/main_figures/fig1_physiology_biochemistry.R` | `Clam_Physiological_and_Biochemical_Endpoints` |
| Fig. 5 | Gut 16S alpha-diversity indices | `scripts/main_figures/fig5_gut_16s_alpha_diversity.R` | `Gut_16S_Alpha_Diversity` |

Additional analysis and figure-generation scripts will be added as they are organized.

## Input data

Input data files are not included in this repository. Place the required data files in `data/input/`.

The English names in the table above are standardized names proposed for the next input-data reorganization. The current R scripts still refer to the original working-file names. The script paths will be synchronized with the reorganized input files when that process is completed. Retain the appropriate file extension (for example, `.csv` or `.xlsx`) when saving the actual data files.

## Running the scripts

1. Install R and RStudio.
2. Open `Plastisphere_Microplastic_Toxicity.Rproj` in RStudio.
3. Place the corresponding input data file in `data/input/`.
4. Run `scripts/install_packages.R` to install the required R packages.
5. Run the required script from the repository root. Figures and other outputs are saved under `results/`.

The scripts assume that the repository root is the working directory. They do not modify the input data files.

## Repository structure

```text
scripts/
  install_packages.R
  main_figures/
    fig1_physiology_biochemistry.R
    fig5_gut_16s_alpha_diversity.R
data/
  input/
    README.md
results/
Plastisphere_Microplastic_Toxicity.Rproj
README.md
README_zh-CN.md
CODE_INVENTORY.md
CODE_AVAILABILITY_DRAFT.md
```

## Code availability statement (draft)

R scripts used for data analysis and figure generation are being organized in this GitHub repository: https://github.com/Xiaokang-He/plastisphere_toxicity. Additional scripts will be added as they are prepared.
