# Install the R packages used by the main-text figure scripts.
# Run from the repository root in an R session.

packages <- c(
  "dplyr",
  "ggplot2",
  "ggpattern",
  "WRS2",
  "svglite",
  "showtext",
  "rstudioapi",
  "pdftools"
)

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

cat("Package installation check complete.\n")
