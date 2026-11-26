# Fig. 3a-b | Proteomic functional enrichment
# Run from the repository root or the Figure_3 folder.
# PDF output goes to Figure_3/results/.

library(ggplot2)

input_name <- "Fig3a_DEP_Functional_Enrichment.csv"
figure_dir <- if (file.exists(file.path("Figure_3", input_name))) "Figure_3" else "."
output_dir <- file.path(figure_dir, "results")
font_family <- if (.Platform$OS.type == "windows") "Times New Roman" else "Liberation Serif"

category_colors <- c(BP = "#c3e1e6", CC = "#dcc6dc", MF = "#f3dfb7", KEGG = "#96c38e")
protein_breaks <- c(4, 8, 12, 16)

read_enrichment <- function(panel) {
  filename <- paste0("Fig3", panel, "_DEP_Functional_Enrichment.csv")
  read.csv(file.path(figure_dir, filename), fileEncoding = "UTF-8-BOM", check.names = FALSE)
}

enrichment <- list(a = read_enrichment("a"), b = read_enrichment("b"))
x_max <- ceiling(max(enrichment$a$Minus_Log10_Adjusted_P_Value,
                      enrichment$b$Minus_Log10_Adjusted_P_Value))
max_protein_count <- max(protein_breaks, enrichment$a$Protein_Count, enrichment$b$Protein_Count)

plot_enrichment <- function(terms, title, legend_top) {
  terms$Category <- factor(terms$Category, levels = names(category_colors))
  terms <- terms[order(terms$Category, -terms$Minus_Log10_Adjusted_P_Value), ]
  terms$Row <- rev(seq_len(nrow(terms)))
  terms$Label <- sub(
    " (SSU-rRNA, 5.8S rRNA, LSU-rRNA)", "\n(SSU-rRNA, 5.8S rRNA, LSU-rRNA)",
    terms$Term_Description, fixed = TRUE
  )

  bar_width <- 0.82
  category_rows <- split(terms$Row, terms$Category, drop = TRUE)
  strips <- data.frame(
    Category = names(category_rows),
    ymin = vapply(category_rows, min, numeric(1)) - bar_width / 2,
    ymax = vapply(category_rows, max, numeric(1)) + bar_width / 2
  )

  ggplot(terms, aes(x = Minus_Log10_Adjusted_P_Value, y = Row)) +
    geom_col(aes(fill = Category), orientation = "y", width = bar_width, alpha = 0.8,
             show.legend = TRUE) +
    geom_rect(
      aes(xmin = -0.70, xmax = -0.40, ymin = ymin, ymax = ymax, fill = Category),
      data = strips, inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_text(
      aes(x = -0.55, y = (ymin + ymax) / 2, label = Category),
      data = strips, inherit.aes = FALSE,
      family = font_family, size = 9.8, size.unit = "pt"
    ) +
    geom_point(aes(x = -0.20, size = Protein_Count, fill = Category),
               shape = 21, colour = "#292929", stroke = 0.24,
               show.legend = c(fill = FALSE, size = TRUE)) +
    geom_text(aes(x = -0.20, label = Protein_Count),
              family = font_family, size = 10.4, size.unit = "pt", show.legend = FALSE) +
    geom_text(aes(x = 0.046, label = Label), hjust = 0, lineheight = 0.82,
              family = font_family, size = 11.0, size.unit = "pt", show.legend = FALSE) +
    scale_fill_manual(values = category_colors, limits = names(category_colors), guide = "none") +
    scale_size_area(
      name = "Protein Number", max_size = 8.2,
      limits = c(0, max_protein_count), breaks = protein_breaks,
      guide = guide_legend(override.aes = list(fill = "white"))
    ) +
    scale_x_continuous(
      limits = c(-0.76, x_max), breaks = seq.int(0, x_max),
      expand = expansion(mult = 0), guide = guide_axis(cap = "both")
    ) +
    scale_y_continuous(limits = c(0.39, nrow(terms) + 0.50),
                       breaks = NULL, expand = expansion(mult = 0)) +
    labs(title = title, x = expression(-log[10] * "(adjusted " * italic(P) * ")"), y = NULL) +
    theme_classic(base_family = font_family, base_size = 11) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 13.2, margin = margin(b = 3)),
      plot.margin = margin(2, 22, 2, 15),
      axis.line.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.34, colour = "#222222"),
      axis.ticks.x = element_line(linewidth = 0.34, colour = "#222222"),
      axis.ticks.length.x = grid::unit(1.4, "mm"),
      axis.text.x = element_text(size = 10.4, colour = "black", margin = margin(t = 1)),
      axis.title.x = element_text(size = 11.6, hjust = 0.60, margin = margin(t = 1)),
      legend.position = "inside",
      legend.position.inside = c(1, legend_top),
      legend.justification.inside = c(0.84, 1),
      legend.title = element_text(size = 12, margin = margin(b = 0.5)),
      legend.text = element_text(size = 10.4),
      legend.key.height = grid::unit(6.4, "mm"),
      legend.key.width = grid::unit(6.4, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_blank()
    )
}

plots <- list(
  a = plot_enrichment(enrichment$a, "(a) Microplastic vs. Control: Enrichment analysis", 0.56),
  b = plot_enrichment(enrichment$b, "(b) Plastisphere vs. Control: Enrichment analysis", 0.78)
)

panel_heights <- 0.74 + 0.333 * vapply(enrichment, nrow, integer(1))

dir.create(output_dir, showWarnings = FALSE)

save_pdf <- function(filename, height, draw) {
  grDevices::cairo_pdf(file.path(output_dir, filename), width = 6.2, height = height, bg = "white")
  on.exit(grDevices::dev.off())
  draw()
}

for (panel in names(plots)) {
  plot <- plots[[panel]]
  if (interactive()) print(plot)
  save_pdf(paste0("Fig3", panel, ".pdf"), panel_heights[[panel]], function() print(plot))
}

combined_a <- plots$a +
  guides(fill = guide_legend(title = "Category", override.aes = list(alpha = 1)), size = "none") +
  theme(
    legend.position.inside = c(0.84, 0.56),
    legend.justification.inside = c(0, 1),
    legend.key.height = grid::unit(4.8, "mm"),
    legend.key.width = grid::unit(4.8, "mm"),
    legend.key.spacing.y = grid::unit(0.65, "mm")
  )

combined_heights <- c(panel_heights[["a"]], 0.12, panel_heights[["b"]])
save_pdf("Fig3ab.pdf", sum(combined_heights), function() {
  layout <- grid::grid.layout(3, heights = grid::unit(combined_heights, "null"))
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = layout))
  print(combined_a, newpage = FALSE, vp = grid::viewport(layout.pos.row = 1))
  print(plots$b, newpage = FALSE, vp = grid::viewport(layout.pos.row = 3))
  grid::popViewport()
})
