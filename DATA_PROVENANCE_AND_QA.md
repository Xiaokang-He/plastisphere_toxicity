# Data provenance and quality-control notes

## Source workbook
`Source Data to Main Figures-DESKTOP-JPDUFHL(4).xlsx` (18 worksheets; latest matching workbook available in the project Library during this task).

## Processing
- Removed blank-only rows/columns and merged-table title rows.
- Retained one valid column-header row per CSV.
- Split Fig. 3ab into panels a and b.
- Split Fig. 5g into node and edge tables.
- Split Fig. 5a–e into alpha-diversity, ordination coordinates, PCoA variance, ANOSIM and NMDS stress tables.
- Split Fig. 7a–c by panel and Fig. 7de into Pearson/Mantel tables, then panels d/e.
- Extracted Fig. 5h–k and Fig. 5l directly from the workbook because separate CSVs were not supplied.
- Added R-compatible group codes to Fig. 1 and Fig. 5 sample-level tables. No measured values or statistical results were recalculated.

## Cross-check notes
- Fig. 1 uploaded CSV headers contained encoding damage in units; source-workbook headers were used.
- Fig. 2a, Fig. 4a/b, Fig. 4c, Fig. 5 ordination, Fig. 5f and Fig. 7 uploaded CSVs showed reduced numeric precision in some columns; full-precision source-workbook values were retained.
- Fig. 4b had a blank trailing row/column.
- Fig. 5a–e contained several vertically stacked data blocks; these were separated into independent CSVs.
- Fig. 5g contained side-by-side node and edge tables; these were separated.
- Fig. 7a–c contained blank rows and an encoding-damaged arrow; the source-workbook text was retained.
- Fig. 7de contained side-by-side Pearson and Mantel tables; these were separated by table and panel.
- Fig. 5g: the `Order` field for `Plastisphere-Low / OTU1656` is numeric `46350` in the source workbook, while the uploaded CSV displayed `11月24日`. The value was preserved exactly as stored in the workbook and requires verification against the original taxonomy export.

See `data_manifest.csv` for file-level source sheets, panel coverage and row counts.
