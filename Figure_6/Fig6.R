# Fig. 6 | KEGG ReporterScore heatmap
# Run from the repository root or the Figure_6 folder.
# Input data: Fig6_KEGG_ReporterScore.csv
# Output: Figure_6/results/Fig6/ (when run from the repository root)

library(ggplot2)

input_name <- "Fig6_KEGG_ReporterScore.csv"
figure_dir <- if (file.exists(file.path("Figure_6", input_name))) "Figure_6" else "."
output_dir <- file.path(figure_dir, "results", "Fig6")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(file.path(figure_dir, input_name),
                fileEncoding = "UTF-8-BOM", check.names = FALSE,
                stringsAsFactors = FALSE)

# Pathway order and functional domains (top to bottom)
pathways <- data.frame(
  KEGG_Pathway_ID = c(
    "ko05100", "ko05130", "ko05132",
    "ko03070", "ko02010", "ko02020",
    "ko00500", "ko01200", "ko00630", "ko00640", "ko00071", "ko04973",
    "ko00970", "ko00290", "ko00270", "ko00564",
    "ko00910", "ko00720"
  ),
  Group = c(
    rep("I. Host–interaction–associated annotations", 3),
    rep("II. Environmental sensing and transport", 3),
    rep("III. Carbon and nutrient metabolism", 6),
    rep("IV. Amino–acid and biosynthetic metabolism", 4),
    rep("V. Energy metabolism", 2)
  ),
  stringsAsFactors = FALSE
)
pathways$Pathway <- dat$KEGG_Pathway_Name[match(pathways$KEGG_Pathway_ID, dat$KEGG_Pathway_ID)]

# Row positions with small gaps between functional domains
group_id <- match(pathways$Group, unique(pathways$Group))
pathways$y <- nrow(pathways) - seq_len(nrow(pathways)) + 1 -
  (group_id - 1) * 0.22

# Italicize species names in pathway labels
pathways$plot_label <- paste0('"', pathways$Pathway, '"')
pathways$plot_label[pathways$KEGG_Pathway_ID == "ko05130"] <-
  'paste("Pathogenic", italic("Escherichia coli"), "infection")'
pathways$plot_label[pathways$KEGG_Pathway_ID == "ko05132"] <-
  'paste(italic("Salmonella"), "infection")'

# Prepare the two comparison columns: MP first, plastisphere second
comparison_x <- c(Microplastic = 8.55, Plastisphere = 9.65)
dat$x <- unname(comparison_x[dat$Exposure_Type])
dat$y <- pathways$y[match(dat$KEGG_Pathway_ID, pathways$KEGG_Pathway_ID)]
if (anyNA(dat$x) || anyNA(dat$y) ||
    anyDuplicated(dat[c("Exposure_Type", "KEGG_Pathway_ID")]) ||
    nrow(dat) != length(comparison_x) * nrow(pathways)) {
  stop("Expected one row per exposure and pathway.")
}
dat$Label <- paste0(
  dat$Total_KO_Count, " (", dat$Upregulated_KO_Count, ")",
  ifelse(abs(dat$Reporter_Score) >= 1.65, "*", "")
)

# Functional-domain blocks
group_names <- unique(pathways$Group)
group_blocks <- data.frame(
  Group = group_names,
  ymin = sapply(group_names, function(g)
    min(pathways$y[pathways$Group == g]) - 0.45),
  ymax = sapply(group_names, function(g)
    max(pathways$y[pathways$Group == g]) + 0.45)
)
group_blocks$ymid <- (group_blocks$ymin + group_blocks$ymax) / 2

# Custom horizontal Reporter Score colour scale
bar_scores <- seq(-3, 3, length.out = 100)
bar_data <- data.frame(
  x = seq(0.55, 3.75, length.out = 100),
  y = -1.48,
  score = bar_scores
)
bar_tick_scores <- seq(-3, 3, by = 1)
bar_ticks <- data.frame(
  score = bar_tick_scores,
  x = 0.55 + (bar_tick_scores + 3) / 6 * 3.20
)

# Comparison labels below the heatmap
comparison_labels <- data.frame(
  x = unname(comparison_x) + 0.4,
  y = -0.62,
  label = paste(names(comparison_x), "vs Control")
)

font_family <- "Times New Roman"

p <- ggplot() +
  geom_rect(
    data = group_blocks,
    aes(xmin = 0.10, xmax = 4.20, ymin = ymin, ymax = ymax),
    fill = "#E6E6E6", colour = NA
  ) +
  geom_text(
    data = group_blocks,
    aes(x = 2.15, y = ymid, label = Group),
    family = font_family, size = 3.5
  ) +
  geom_text(
    data = pathways,
    aes(x = 7.95, y = y, label = plot_label),
    parse = TRUE, hjust = 1,
    family = font_family, size = 3.5
  ) +
  geom_tile(
    data = dat,
    aes(x = x, y = y, fill = Reporter_Score),
    width = 0.98, height = 0.86,
    colour = "white", linewidth = 0.35
  ) +
  geom_text(
    data = dat,
    aes(x = x, y = y, label = Label),
    family = font_family, size = 3.35
  ) +
  geom_text(
    data = comparison_labels,
    aes(x = x, y = y, label = label),
    angle = 35, hjust = 1, vjust = 1,
    family = font_family, size = 3.5
  ) +
  geom_tile(
    data = bar_data,
    aes(x = x, y = y, fill = score),
    width = 0.033, height = 0.20
  ) +
  annotate(
    "text", x = 2.15, y = -1.12, label = "Reporter score",
    family = font_family, fontface = "bold", size = 3.5
  ) +
  geom_segment(
    data = bar_ticks,
    aes(x = x, xend = x, y = -1.61, yend = -1.68),
    linewidth = 0.25
  ) +
  geom_text(
    data = bar_ticks,
    aes(x = x, y = -1.86, label = score),
    family = font_family, size = 3.2
  ) +
  scale_fill_gradient2(
    low = "#EC6D79", mid = "white", high = "#5EA04C",
    midpoint = 0, limits = c(-3, 3), guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0, 10.40), ylim = c(-2.20, 18.65),
    expand = FALSE, clip = "off"
  ) +
  theme_void(base_family = font_family) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    text = element_text(colour = "black")
  )

if (interactive()) print(p)

ggsave(
  file.path(output_dir, "Fig6_KEGG_ReporterScore.pdf"),
  p, width = 9, height = 6, units = "in", device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Fig6_KEGG_ReporterScore.png"),
  p, width = 9, height = 6, units = "in", dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Fig6_KEGG_ReporterScore.tiff"),
  p, width = 9, height = 6, units = "in", dpi = 600,
  compression = "lzw", bg = "white"
)
