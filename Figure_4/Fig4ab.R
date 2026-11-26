# Fig. 4a,b | Transcriptome-proteome nine-quadrant plots
# Run from the repository root or from the Figure_4 folder.
# Install once if needed: install.packages(c("ggplot2", "patchwork"))

library(ggplot2)
library(patchwork)

data_dir <- if (file.exists("Figure_4/Fig4a_Nine_Quadrant_Microplastic.csv")) {
  "Figure_4"
} else {
  "."
}
output_dir <- file.path(data_dir, "results")
output_name <- "Fig4ab"

x_limits <- c(-10.5, 10.5)
y_limits <- c(-15.5, 15.5)

font_family <- "serif"
if (.Platform$OS.type == "windows" && capabilities("cairo")) {
  font_family <- "Times New Roman"
  grDevices::windowsFonts(
    "Times New Roman" = grDevices::windowsFont("TT Times New Roman")
  )
}

plot_colours <- c(
  "DEGs(Up) and DEPs(Up)"     = "#C44E2B",
  "DEGs(Up) and DEPs(Down)"   = "#7186B8",
  "DEGs(Down) and DEPs(Down)" = "#59AD62",
  "DEGs(Down) and DEPs(Up)"   = "#74C8DF",
  "Others"                   = "#D9D9D9"
)

make_panel <- function(csv_name, comparison, panel_tag) {
  csv_path <- file.path(data_dir, csv_name)
  if (!file.exists(csv_path)) stop("CSV file not found: ", csv_path, call. = FALSE)
  points <- read.csv(csv_path, fileEncoding = "UTF-8-BOM")

  # Use the same data for the plot, corner counts and Pearson test.
  points <- subset(
    points,
    Protein_log2FC >= x_limits[1] & Protein_log2FC <= x_limits[2] &
      Gene_log2FC >= y_limits[1] & Gene_log2FC <= y_limits[2]
  )
  highlighted <- points$Plot_Category != "Others"
  category_counts <- table(factor(points$Plot_Category, levels = names(plot_colours)))
  corner_categories <- names(plot_colours)[c(2, 1, 3, 4)]
  count_df <- data.frame(
    x = rep(c(x_limits[1] + 1, x_limits[2] - 1), times = 2),
    y = rep(c(y_limits[2] - 2, y_limits[1] + 2), each = 2),
    label = sprintf("(%d)", as.integer(category_counts[corner_categories])),
    hjust = c(0, 1, 0, 1)
  )

  pearson <- cor.test(points$Protein_log2FC, points$Gene_log2FC, method = "pearson")
  r_label <- sprintf("Pearson correlation: %.3f", unname(pearson$estimate))
  p_label <- if (pearson$p.value < 0.001) {
    "italic(P) < 0.001"
  } else {
    sprintf("italic(P) == %.3f", pearson$p.value)
  }

  ggplot(points, aes(x = Protein_log2FC, y = Gene_log2FC)) +
    geom_point(
      data = points[!highlighted, ],
      colour = plot_colours[["Others"]], size = 1.8, alpha = 0.75
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.55) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.55) +
    annotate("rect", xmin = -6.2, xmax = 6.2, ymin = 11.75, ymax = 15.4, fill = "white") +
    geom_point(
      data = points[highlighted, ], aes(fill = Plot_Category),
      shape = 21, colour = "#282828", size = 2.4, stroke = 0.32,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = plot_colours, limits = names(plot_colours), drop = FALSE, name = NULL,
      guide = guide_legend(ncol = 1, override.aes = list(size = 3.5))
    ) +
    annotate(
      "text", x = 0, y = 14.45, label = comparison,
      fontface = "bold", family = font_family, size = 6
    ) +
    annotate("text", x = 0, y = 13.1, label = r_label, family = font_family, size = 4.7) +
    annotate(
      "text", x = 0, y = 12.15, label = p_label,
      parse = TRUE, family = font_family, size = 4.5
    ) +
    geom_text(
      data = count_df, aes(x = x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE, family = font_family, size = 5
    ) +
    scale_x_continuous(breaks = seq(-9, 9, by = 3)) +
    scale_y_continuous(breaks = seq(-15, 15, by = 5)) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +
    labs(
      x = "Log2(Fold Change) of Protein",
      y = "Log2(Fold Change) of Gene",
      tag = panel_tag
    ) +
    theme_bw(base_size = 14, base_family = font_family) +
    theme(
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(colour = "black", size = 12),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(colour = "#D0D0D0", linewidth = 0.4),
      panel.grid.minor = element_line(colour = "#E7E7E7", linewidth = 0.3),
      legend.text = element_text(size = 11),
      legend.key.size = grid::unit(6.5, "mm"),
      plot.tag = element_text(size = 21, hjust = 0, vjust = 1),
      plot.tag.position = "topleft",
      plot.margin = margin(t = 10, r = 7, b = 6, l = 14)
    )
}

p_a <- make_panel("Fig4a_Nine_Quadrant_Microplastic.csv", "Microplastic vs Control", "(a)")
p_b <- make_panel("Fig4b_Nine_Quadrant_Plastisphere.csv", "Plastisphere vs Control", "(b)")

note <- paste(
  "Each point represents a matched gene\u2013protein pair. Points within the plot limits were used",
  "for the counts and two-sided Pearson tests. Dashed lines indicate log2FC cutoffs of \u00b11.",
  sep = "\n"
)

figure_4ab <- (
  p_a / p_b +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = note,
      theme = theme(
        plot.caption = element_text(family = font_family, size = 10, hjust = 0),
        plot.margin = margin(5, 7, 6, 7)
      )
    )
) & theme(legend.position = "bottom")

dir.create(output_dir, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(
  file.path(output_dir, paste0(output_name, ".pdf")),
  plot = figure_4ab, device = pdf_device, width = 8.3, height = 16.5, units = "in", bg = "white"
)
ggsave(
  file.path(output_dir, paste0(output_name, ".png")),
  plot = figure_4ab, width = 8.3, height = 16.5, units = "in", dpi = 600, bg = "white"
)

if (interactive()) print(figure_4ab)
