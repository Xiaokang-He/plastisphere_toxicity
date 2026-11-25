# Plastisphere formation reshapes the ecological toxicity of microplastics in *Corbicula fluminea*

This repository contains selected R scripts used for statistical analysis, ordination/source-data preparation, and figure generation for the manuscript submitted to *Nature Communications*.

**Manuscript:** “Plastisphere formation reshapes the ecological toxicity of microplastics in *Corbicula fluminea*.”

## Scope

The repository contains selected downstream R scripts. It is **not yet a complete reimplementation of every custom analysis in the study**. Transcriptomic, 16S rRNA gene, and shotgun-metagenomic read processing and primary bioinformatics analyses were conducted using the Majorbio Cloud platform and the software/parameters described in the Methods and Supporting Information. Platform-internal workflows are not reproduced here.

The full custom R workflow that generated the cross-omics Pearson Mantel tests and the associated Benjamini–Hochberg correction has not yet been identified in the authors’ archived scripts. The repository therefore does not currently reproduce the Mantel association statistics underlying Supplementary Fig. 22 and related analyses. The included Supplementary Fig. 23 script reproduces the coverage plot from the reported counts; it does **not** rerun Mantel tests or FDR correction. Do not describe this repository as a complete computational reproduction until the original Mantel analysis code is added or its omission is explicitly resolved.

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
| `supp_fig23_mantel_coverage.R` | Supplementary Fig. 23: coverage plot from reported Mantel association counts | No external input; reported counts are defined in the script |

The RDA script performs the Hellinger transformation, fits `vegan::rda()` with the environmental variables, calculates `vegan::envfit()` with 999 permutations, and exports figure/source-data tables. It does not reproduce the original figure-drawing code if that was performed separately. The Supplementary Fig. 23 script uses the counts already reported for the selected maps; it does not recompute the underlying Mantel tests.

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

## Important limitations of this draft

- This package has undergone a **static path and content review only**. It has not been executed end-to-end in this environment because an R runtime and the required input datasets were not available. Run and verify each script in the authors’ R environment before public release or submission.
- The package covers selected downstream analyses and figures; it does not contain the Majorbio Cloud sequencing pipelines.
- The full custom R code for the Pearson Mantel tests and BH correction used for the microbial–host cross-omics associations is still missing. The topology and coverage scripts are downstream of precomputed association results and must not be presented as substitutes for the original test-generation workflow. Add the original script and its required inputs before claiming that the repository fully reproduces the Mantel analyses.
- The RDA fit/envfit and source-data exporter is included, but the RDA plot-generation script itself was not identified.
- The particle-property script containing the carbonyl-index panel is not included because the final carbonyl-index input/version needs to be checked against the latest Supporting Information before release.
- The saved KEGG Level 3 PCA script is not included because its saved version also calculates additional PERMANOVA/PERMDISP/ANOSIM outputs that are not described in the current Supporting Information text. Reconcile the script with the final reported workflow before adding it.
- Package versions and the exact R session information have not yet been frozen. Record `sessionInfo()` from the verified analysis environment before final release.
- No reuse license has been selected. Confirm the authors’ preferred license before adding a `LICENSE` file.

## Code availability statement (draft)

Selected custom R scripts for statistical analyses and figure-generation workflows are available at https://github.com/hexiaokang1214-creator/plastisphere-microplastic-toxicity. Primary sequencing-data processing and bioinformatics analyses were performed using the Majorbio Cloud platform and the software described in the Methods and Supporting Information.

**This is not yet a final Code Availability statement.** Before using it in the manuscript, add the missing custom Mantel test/BH workflow or obtain an explicit decision on its omission, verify every included script against the final manuscript and Source Data, replace the placeholder with the public repository URL, and record the verified R/package versions.