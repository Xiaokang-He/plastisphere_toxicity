# Code inventory and release decisions

## Included scripts

| Repository script | Source working file | Scope / release note |
|---|---|---|
| `scripts/main_figures/fig1_physiology_biochemistry.R` | `Fig1_生理生化指标_WRS2_修复SVG输出路径版.R` | WRS2-based Fig. 1 statistical/plotting workflow. |
| `scripts/main_figures/fig5_gut_16s_alpha_diversity.R` | `肠道微生物16S_alpha多样性_WRS2统一绘图版_240x80mm.R` | WRS2-based Shannon, Pielou and ACE panels. |
| `scripts/validation/rt_qpcr_validation.R` | `RTqPCR_WRS2_t1way_lincon_完整修订代码.R` | WRS2-based RT-qPCR analysis and plotting. |
| `scripts/supplementary/supp_fig2_eps_quantification.R` | `EPS定量_WRS2_t1way_lincon_一行三面板_570x200mm_PDF版.R` | WRS2-based EPS component analysis. |
| `scripts/supplementary/supp_fig3_parafac_fmax.R` | `C1_FMax_WRS2_t1way_lincon_未校正P值_修复print报错版.R` | WRS2-based PARAFAC Fmax analysis. |
| `scripts/supplementary/supp_siphoning_timecourse.R` | `虹吸效率_WRS2_t1way_lincon_未校正P值_完整修正版.R` | WRS2-based preliminary time-course siphoning analysis. |
| `scripts/supplementary/supp_fig16_acinetobacter_burkholderia.R` | `Acinetobacter_Burkholderia_WRS2_robust_SD_errorbars_corrected.R` | WRS2-based genus/species composition plots. |
| `scripts/supplementary/supp_fig19_vfg_cumulative_abundance.R` | `Core105_VFDB_WRS2_t1way_lincon_no_adjustment.R` | WRS2-based cumulative Core VFG abundance. |
| `scripts/supplementary/supp_fig20_vfg_heatmap.R` | `Core105_VFDB_WRS2_t1way_lincon_rawP(1).R` | WRS2-based per-ID Core VFG heatmaps and comparisons. |
| `scripts/ordination/supp_fig24_rda_source_data.R` | `Export_RDA_Figure_Source_Data.R` | Standalone Hellinger/RDA/envfit fit and figure/source-data export; plot-generation code not found. |
| `scripts/supplementary/supp_fig23_mantel_coverage.R` | `Supplementary_Fig23_updated.R` | Recreates the coverage plot from reported counts; does not rerun Mantel tests or BH correction. |

For the nine selected WRS2 scripts, the repository copies adapt local absolute paths to the project-relative `data/input/` and `results/` structure. The statistical and plotting operations were retained. The two additional scripts adapt only their file paths; the RDA model/envfit and reported coverage counts are retained.

## Not included yet — needs author verification

| Candidate workflow | Reason to hold |
|---|---|
| Full microbial–host Pearson Mantel tests and global BH correction | The saved R files reviewed include source-data exporters, network topology calculations and downstream visualizations, but no original script that generates the full Mantel r, permutation P and BH-adjusted P results. This is a material gap because the current SI describes these analyses. |
| `Mantel_network_topology.R` and `Export_Mantel_Network_Source_Data.R` | These read existing association results; they do not recalculate Mantel tests or BH correction. Include only if the corresponding derived inputs and their role in the final paper are documented. |
| `seven_sets_half_A4_size_linear.R` and `Supplementary_Fig23_updated.R` alternatives | The selected Fig. 23 coverage script is included. Other network visualization variants need confirmation against the final SI figures before release. |
| `Core_VFDB_UpSet_from_local_raw.R`, `Core_VFDB_4group_scatter_from_local_raw.R` | Saved versions state that they were not executed in the generation environment; verify that they generated the final Supplementary Fig. 21 before release. |
| `塑料表面理化性质_六面板_WRS2_robust_精修版_导出修复(1).R` | Includes the carbonyl-index panel. Confirm that its input file and values match the latest final SI after the FTIR update. |
| `KEGG_Level3_Hellinger_PCA_Final_Corrected.R` | Saved version also calculates PERMANOVA, PERMDISP and ANOSIM; current SI text describes PCA but does not describe those additional tests. Reconcile before adding. |

## Required final checks

- Add the original Mantel test/BH script and its required inputs, or obtain an explicit author/editorial decision on its omission.
- Verify each included script against the final figure and final Source Data.
- Run each script from the repository root with the corresponding input file.
- Save the actual R `sessionInfo()` and package versions.
- Confirm the published Source Data can be mapped to each required input table, or document the necessary extraction/formatting steps.
- Select an appropriate reuse license only after author confirmation.