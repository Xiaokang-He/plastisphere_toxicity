# Figure 1 — Siphoning and biochemical endpoints

`Fig1a_Clam_Physiological_and_Biochemical_Endpoints.csv` contains fifteen tank-level records, with three independent tanks per treatment. Panels a–j show siphoning efficiency, SOD, CAT, GSH, GSH-Px, MDA, AChE, lysozyme, D-Lac and DAO, in that order.

`Fig1.R` reads the group from the `sample` ID. Points show tank measurements and bars show means ± SD. Brackets show unadjusted pairwise P values below 0.05 from `WRS2::lincon()` at 20% trimming.

Run from the repository root with `ggplot2`, `ggpattern`, `WRS2` and Cairo graphics support.

```sh
Rscript Figure_1/Fig1.R
```

The script writes `Fig1.pdf` and `Fig1.jpg` (300 dpi) alongside the input table. It can also be run as `Rscript Fig1.R` from `Figure_1/`.
