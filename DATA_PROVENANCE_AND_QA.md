# Data provenance and quality-control notes

## Source workbook

The cleaned CSV datasets were prepared from an 18-worksheet source workbook titled `Source Data to Main Figures`. The workbook is not included in this repository; the working-copy filename suffix has been omitted.

## Processing

- Removed blank-only rows/columns and merged-table title rows.
- Retained one valid column-header row per CSV.
- Split Fig. 3ab into panels a and b.
- Split the Fig. 5g worksheet into node and edge tables for network topology calculations.
- Extracted the alpha-diversity data block from the Fig. 5a–e worksheet.
- Split Fig. 7a–c by panel and Fig. 7de into Pearson/Mantel tables, then panels d/e.
- Added R-compatible group codes to Fig. 1 and Fig. 5 sample-level tables. No measured values or statistical results were recalculated.

## Cross-check notes

- Fig. 1 uploaded CSV headers contained encoding damage in units; source-workbook headers were used.
- Fig. 2a, Fig. 4a/b, Fig. 4c, Fig. 5f and Fig. 7 uploaded CSVs showed reduced numeric precision in some columns; full-precision source-workbook values were retained.
- Fig. 4b had a blank trailing row/column.
- Fig. 7a–c contained blank rows and an encoding-damaged arrow; the source-workbook text was retained.
- Fig. 7de contained side-by-side Pearson and Mantel tables; these were separated by table and panel.
- **Pending verification:** the `Order` field for `Plastisphere-Low / OTU1656` in `Figure_5/Fig5g_Network_Topology_Raw_Nodes.csv` retains the workbook value `46350`, whereas an earlier uploaded CSV displayed `11月24日`. This entry has not been independently verified against the original taxonomy export. The value is retained as recorded in the source workbook and should not be treated as a confirmed taxonomic assignment until checked against that export.

See `data_manifest.csv` for file-level source sheets, panel coverage, row counts, column counts and the current committed file sizes in bytes.
