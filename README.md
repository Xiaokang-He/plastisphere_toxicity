# Plastisphere formation reshapes the ecological toxicity of microplastics in *Corbicula fluminea*

This repository contains selected R scripts used for statistical analysis, ordination/source-data preparation, and figure generation for the manuscript submitted to *Nature Communications*.

**Manuscript:** “Plastisphere formation reshapes the ecological toxicity of microplastics in *Corbicula fluminea*.”

## Scope

This repository is being organized to accompany the manuscript. It currently contains selected R scripts for statistical analyses, data processing, source-data preparation, and figure generation.

The scripts currently cover:
- Physiological and biochemical endpoint analyses for Main Fig. 1.
- Gut 16S alpha-diversity analyses for Main Fig. 5.
- RT–qPCR validation.
- EPS quantification and PARAFAC fluorescence-component analysis.
- The preliminary siphoning-efficiency time course.
- Gut microbial composition analyses involving *Acinetobacter* and *Burkholderia*.
- VFDB-related abundance summaries and heatmaps.
- RDA fitting, envfit analysis, and source-data export.
- Mantel-association coverage visualization.

The R scripts and related documentation are still being organized, and additional materials will be added as they are prepared. Transcriptomic, gut 16S rRNA gene, and shotgun-metagenomic bioinformatics analyses were conducted using the Majorbio Cloud platform; the corresponding tools and version information are provided in the Methods and Supporting Information.

## Repository structure

```text
scripts/
  install_packages.R
  main_figures/
    fig1_physiology_biochemistry.R
    fig5_gut_16s_alpha_diversity.R
  validation/
    rt_qpcr_validation.R
  ordination/
    supp_fig24_rda_source_data.R
  supplementary/
    supp_fig2_eps_quantification.R
    supp_fig3_parafac_fmax.R
    supp_siphoning_timecourse.R
    supp_fig16_acinetobacter_burkholderia.R
    supp_fig19_vfg_cumulative_abundance.R
    supp_fig20_vfg_heatmap.R
    supp_fig23_mantel_coverage.R
data/
  input/                 # author-supplied input files are not bundled
results/                 # generated outputs
```

## Included scripts and required input files

Place the required input files in `data/input/` using the filenames below. Input data are intentionally not bundled in this code-only repository. The filenames below are the original analysis-input filenames; the published Source Data workbooks contain figure-level data but are not guaranteed to be drop-in replacements for these working input tables. Check the required columns and formats in each script before running it. Do not upload restricted or unpublished data to a public repository without checking the applicable data-sharing permissions.

| Script | Manuscript item / analysis | Required input file(s) |
|---|---|---|
| `fig1_physiology_biochemistry.R` | Main Fig. 1: physiological and biochemical endpoints; WRS2 tests | `酶活数据_最终 - 8.30.csv` |
| `fig5_gut_16s_alpha_diversity.R` | Main Fig. 5a–c: gut 16S alpha diversity | `α多样性数据.csv` |
| `rt_qpcr_validation.R` | Supplementary RT–qPCR validation plots and statistics | `rt-qPCR.csv` |
| `supp_fig2_eps_quantification.R` | Supplementary Fig. 2: EPS components | `EPS含量.csv` |
| `supp_fig3_parafac_fmax.R` | Supplementary Fig. 3e: PARAFAC component Fmax | `3-component_result - 副本.xlsx` (sheet `FMax`) |
| `supp_siphoning_timecourse.R` | Preliminary time-course siphoning assay | `虹吸效率.csv` |
| `supp_fig16_acinetobacter_burkholderia.R` | Supplementary Fig. 16: genus/species composition | `Acinetobacter丰度表.csv` |
| `supp_fig19_vfg_cumulative_abundance.R` | Supplementary Fig. 19: cumulative Core VFG abundance | `vfdb_core_VFDB_ID_profile.xls` (tab-delimited text despite the `.xls` extension) |
| `supp_fig20_vfg_heatmap.R` | Supplementary Fig. 20: Core VFG heatmaps and WRS2 comparisons | `vfdb_core_VFDB_ID_profile.xls`; `VFDB毒力因子丰度表_VFDB_Origin.csv` |
| `supp_fig24_rda_source_data.R` | Supplementary Fig. 24: Hellinger-transformed OTU RDA and envfit/source-data export | `OTU - 副本.csv`; `环境因子 - 副本.csv`; `分组依据 - 副本.csv` |
| `supp_fig23_mantel_coverage.R` | Supplementary Fig. 23: Mantel-association coverage visualization | No external input; reported counts are defined in the script |

The RDA script includes Hellinger transformation, RDA fitting, envfit analysis with 999 permutations, and source-data export. The Mantel-related script prepares the association-coverage visualization.

## Running the scripts

1. Install R (record the exact version used for the final analysis before release).
2. Open `Plastisphere_Microplastic_Toxicity.Rproj` in RStudio. This sets the working directory to the repository root.
3. Place the required input tables in `data/input/`.
4. Run `scripts/install_packages.R` once.
5. Run the relevant script individually. Outputs are written to the corresponding subfolders under `results/`.

Scripts assume that they are run with the repository root as the working directory. They do not modify the input tables.

## Statistical methods represented in the included scripts

Several scripts use `WRS2::t1way()` for the omnibus robust one-way comparison and `WRS2::lincon()` for pairwise comparisons, with a trimming proportion of `tr = 0.20`. Where specified in the source scripts, `lincon(method = "none")` is used and the pairwise P values are unadjusted. Figure-specific annotation gates are retained as implemented in the individual scripts; consult the comments and exported analysis notes in each script rather than assuming one annotation rule applies to every figure.

With three biological replicates per group, `floor(0.20 × 3) = 0`; therefore, the 20% trimming setting does not remove observations from these groups.

## Ongoing organization

The repository is being progressively organized. The scripts currently included are listed above, and additional analysis scripts and documentation will be added as they are prepared. Input datasets are not bundled in this repository; the required input filenames and formats are listed in the table above.

## Code availability statement (draft)

Selected R scripts for statistical analyses and figure generation are being organized in this GitHub repository: https://github.com/Xiaokang-He/plastisphere_toxicity.

Transcriptomic, gut 16S rRNA gene, and shotgun-metagenomic bioinformatics analyses were conducted using the Majorbio Cloud platform. The corresponding software tools and version information are provided in the Methods and Supporting Information.
