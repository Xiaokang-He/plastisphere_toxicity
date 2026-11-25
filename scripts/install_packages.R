# Install the R packages used by the scripts in this repository.
# Run from the repository root in an R session.

packages <- c(
  "dplyr", "tidyr", "ggplot2", "ggpattern", "scales", "WRS2",
  "patchwork", "cowplot", "ggalluvial", "ggsignif", "readxl",
  "stringr", "svglite", "showtext", "vegan", "writexl",
  "pdftools", "rstudioapi", "vegan", "writexl"
)

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

cat("Package installation check complete.\n")