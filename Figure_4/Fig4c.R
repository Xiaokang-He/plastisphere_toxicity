# Fig. 4c | Transcriptome-proteome Pearson correlation heatmap
# Run from the repository root or from the Figure_4 folder.

library(ggplot2)

input_name <- "Fig4c_Transcriptome_Proteome_Correlation.csv"
figure_dir <- if (file.exists(file.path("Figure_4", input_name))) "Figure_4" else "."
input_file <- file.path(figure_dir, input_name)
output_dir <- file.path(figure_dir, "results")
plot_width_mm <- 130
plot_height_mm <- 197
plot_family <- "Times New Roman"

if (!file.exists(input_file)) {
  stop(
    "Could not find ", input_name, ". Run from the repository root or the Figure_4 folder.",
    call. = FALSE
  )
}
dat <- read.csv(input_file, fileEncoding = "UTF-8-BOM", check.names = FALSE)

database_order <- c("GO", "KEGG")
comparison_labels <- c("Microplastic vs Control", "Plastisphere vs Control")
r_cols <- c("Microplastic_Pearson_r", "Plastisphere_Pearson_r")
q_cols <- c("Microplastic_Adjusted_P", "Plastisphere_Adjusted_P")
required_cols <- c(
  "Database", "Term_ID", "Term_Name", "Annotated_Gene_Protein_Pair_Count",
  r_cols, q_cols
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols)) {
  stop("Missing CSV columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (!setequal(unique(dat$Database), database_order)) {
  stop("The Database column must contain GO and KEGG rows.", call. = FALSE)
}

# Stack GO above KEGG, keeping the source order within each section.
dat <- dat[order(match(dat$Database, database_order), seq_len(nrow(dat))), ]
dat$y <- rev(seq_len(nrow(dat))) + 0.35 * (dat$Database == "GO")
dat$term_label <- sprintf("%s (%d)", dat$Term_Name, dat$Annotated_Gene_Protein_Pair_Count)
long_term <- dat$Term_ID == "GO:0019829"
dat$term_label[long_term] <- sub(
  "monoatomic ", "monoatomic\n", dat$term_label[long_term], fixed = TRUE
)

heatmap <- data.frame(
  x = rep(seq_along(comparison_labels), each = nrow(dat)),
  y = rep(dat$y, times = length(comparison_labels)),
  r = unlist(dat[r_cols], use.names = FALSE),
  q = unlist(dat[q_cols], use.names = FALSE)
)
heatmap$label <- sub("^-", "\u2212", sprintf("%+.2f", heatmap$r))
significant <- !is.na(heatmap$q) & heatmap$q < 0.05
heatmap$label[significant] <- paste0(heatmap$label[significant], "*")

section_rows <- split(dat$y, factor(dat$Database, levels = database_order))
strips <- data.frame(
  title = c("GO term", "KEGG pathway"),
  ymin = vapply(section_rows, min, numeric(1)) - 0.5,
  ymax = vapply(section_rows, max, numeric(1)) + 0.5
)
strips$y <- (strips$ymin + strips$ymax) / 2
x_labels <- data.frame(x = seq_along(comparison_labels), label = comparison_labels)

# Space for the rotated comparison labels and the shared colour bar.
x_limits <- c(-1.26, 7.74)
y_limits <- c(-3.10, max(dat$y) + 0.80)
legend_anchor <- (c(2.44, -1.95) - c(x_limits[1], y_limits[1])) /
  c(diff(x_limits), diff(y_limits))

p <- ggplot() +
  geom_tile(
    data = heatmap, aes(x, y, fill = r),
    width = 1, height = 1, colour = "white", linewidth = 0.45
  ) +
  geom_text(
    data = heatmap, aes(x, y, label = label),
    family = plot_family, size = 3.30
  ) +
  geom_rect(
    data = strips, aes(ymin = ymin, ymax = ymax),
    xmin = 0.16, xmax = 0.45, fill = c("#C3E2E7", "#DEC5DB")
  ) +
  geom_text(
    data = strips, aes(y = y, label = title), x = 0.305,
    angle = 90, family = plot_family, size = 4.15
  ) +
  geom_text(
    data = dat, aes(y = y, label = term_label), x = 2.54,
    hjust = 0, lineheight = 0.86, family = plot_family, size = 3.48
  ) +
  geom_text(
    data = x_labels, aes(x, label = label), y = 0.25,
    angle = 30, hjust = 1, family = plot_family, size = 3.45
  ) +
  annotate(
    "text", x = -0.03, y = max(dat$y) + 0.37, label = "(c)",
    hjust = 0, family = plot_family, size = 5.80
  ) +
  scale_fill_gradient2(
    low = "#EC6D79", mid = "white", high = "#5EA04C", midpoint = 0,
    limits = c(-0.81, 0.81), breaks = seq(-0.8, 0.8, by = 0.4),
    labels = function(x) sub("^-", "\u2212", sprintf("%.1f", x)),
    name = "Pearson r",
    guide = guide_colourbar(
      direction = "horizontal", display = "gradient",
      theme = theme(
        legend.key.width = grid::unit(39, "mm"),
        legend.key.height = grid::unit(2.5, "mm")
      )
    )
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  theme_void(base_family = plot_family) +
  theme(
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    legend.position = "inside",
    legend.position.inside = legend_anchor,
    legend.justification.inside = c(0, 0),
    legend.title.position = "top",
    legend.title = element_text(
      face = "bold", size = 11, hjust = 0, margin = margin(b = 0.5, unit = "mm")
    ),
    legend.text = element_text(size = 9.5, margin = margin(t = 0.5, unit = "mm")),
    legend.frame = element_rect(colour = "black", linewidth = 0.25),
    legend.ticks = element_line(colour = "black", linewidth = 0.20),
    legend.ticks.length = grid::unit(0.65, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "Fig4c.pdf"), p, device = grDevices::cairo_pdf,
  width = plot_width_mm, height = plot_height_mm, units = "mm", bg = "white"
)
ggsave(
  file.path(output_dir, "Fig4c.png"), p, device = grDevices::png, type = "cairo",
  width = plot_width_mm, height = plot_height_mm, units = "mm", dpi = 600, bg = "white"
)

if (interactive()) {
  print(p)
}
