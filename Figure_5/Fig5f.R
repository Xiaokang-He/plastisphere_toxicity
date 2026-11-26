# Fig. 5f | Genus-level relative abundance
# Run from the repository root or from the Figure_5 folder.

library(ggplot2)
library(ggalluvial)

input_name <- "Fig5f_Genus_Relative_Abundance.csv"
figure_dir <- if (file.exists(file.path("Figure_5", input_name))) "Figure_5" else "."
input_file <- file.path(figure_dir, input_name)
output_dir <- file.path(figure_dir, "results")
plot_width_mm <- 105
plot_height_mm <- 130
plot_family <- "Times New Roman"

group_columns <- c(
  Control = "Control_percent",
  "Microplastic-Low" = "Microplastic_Low_percent",
  "Microplastic-High" = "Microplastic_High_percent",
  "Plastisphere-Low" = "Plastisphere_Low_percent",
  "Plastisphere-High" = "Plastisphere_High_percent"
)

# Legend order, from the top of the stack to the bottom.
genus_colors <- c(
  others = "#ADC3D9",
  unclassified_f__Dermacoccaceae = "#EC0000",
  Caulobacter = "#41B43F",
  norank_f__Caulobacteraceae = "#0099B3",
  unclassified_f__Peptostreptococcaceae = "#925E9F",
  Romboutsia = "#FDAE91",
  norank_o__Chloroplast = "#AC002A",
  Enhydrobacter = "#ACB6B6",
  Bacillus = "#9E3B4F",
  Enterobacter = "#B2811E",
  Methyloversatilis = "#40A67D",
  Paenibacillus = "#6B7EA9",
  Stenotrophomonas = "#D7675A",
  Aerococcus = "#6E4D6D",
  Sphingobium = "#D1793E",
  Bosea = "#5EAD73",
  Aquabacterium = "#6694AE",
  Epulopiscium = "#AE80A1",
  Staphylococcus = "#EAA392",
  Corynebacterium = "#949CAB",
  Clostridium = "#C65356",
  Blastomonas = "#8D9D4B",
  Sphingomonas = "#4D9F9B",
  Brevundimonas = "#8A7CA4",
  Methylobacterium = "#E3ABA9",
  Bradyrhizobium = "#C26161",
  Pseudomonas = "#B69C96",
  Roseateles = "#EC7676",
  "Burkholderia-Caballeronia-Paraburkholderia" = "#5AA6B3",
  Agrobacterium = "#7BB47A",
  Acinetobacter = "#00468B"
)
genus_order <- names(genus_colors)
legend_labels <- as.expression(lapply(genus_order, function(genus) {
  if (genus == "others") return("Other genera")
  if (grepl("^(norank|unclassified)_", genus)) return(genus)
  bquote(italic(.(genus)))
}))

if (!file.exists(input_file)) {
  stop(
    "Could not find ", input_name, ". Run from the repository root or the Figure_5 folder.",
    call. = FALSE
  )
}
dat <- read.csv(input_file, fileEncoding = "UTF-8-BOM", check.names = FALSE)
missing_cols <- setdiff(c("Genus", unname(group_columns)), names(dat))
if (length(missing_cols)) {
  stop("Missing CSV columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(dat$Genus) || !setequal(dat$Genus, genus_order)) {
  stop("The CSV must contain one row for each taxon in the colour palette.", call. = FALSE)
}

# The supplied abundances are already percentages.
abundance <- unlist(dat[unname(group_columns)], use.names = FALSE)
if (!is.numeric(abundance) || any(!is.finite(abundance)) || any(abundance < 0)) {
  stop("Relative abundances must be finite, non-negative numbers.", call. = FALSE)
}
group_order <- names(group_columns)
long <- data.frame(
  group = factor(rep(group_order, each = nrow(dat)), levels = group_order),
  Genus = factor(rep(dat$Genus, times = length(group_order)), levels = genus_order),
  abundance = abundance
)

p <- ggplot(long, aes(x = group, y = abundance, fill = Genus,
                     stratum = Genus, alluvium = Genus)) +
  geom_flow(alpha = 0.2, width = 0.6) +
  geom_stratum(alpha = 0.8, width = 0.6, colour = "black", linewidth = 0.16) +
  scale_fill_manual(
    values = genus_colors, breaks = genus_order, labels = legend_labels,
    name = "Genus", guide = guide_legend(ncol = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = expansion(add = 2)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Relative abundance (%)", tag = "(f)") +
  theme_grey(base_family = plot_family) +
  theme(
    panel.grid.major = element_line(colour = "white", linewidth = 0.28),
    panel.grid.minor = element_line(colour = "white", linewidth = 0.18),
    axis.title.y = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(colour = "black", size = 8),
    axis.text.x = element_text(colour = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_line(colour = "grey25", linewidth = 0.18),
    axis.ticks.length = grid::unit(0.75, "mm"),
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(face = "bold", size = 10, margin = margin(b = 0.2, unit = "mm")),
    legend.text = element_text(size = 6.5, margin = margin(l = 0.7, unit = "mm")),
    legend.key = element_rect(fill = "#EBEBEB", colour = NA),
    legend.key.size = grid::unit(3.5, "mm"),
    legend.key.spacing.y = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(0.9, "mm"),
    legend.box.margin = margin(t = -1.1, unit = "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.tag = element_text(size = 10.5, hjust = 0, vjust = 1),
    plot.tag.position = c(0, 0.997),
    plot.margin = margin(4, 3.5, 2.4, 2, unit = "mm")
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "Fig5f.pdf"), p, device = grDevices::cairo_pdf,
  width = plot_width_mm, height = plot_height_mm, units = "mm", bg = "white"
)
ggsave(
  file.path(output_dir, "Fig5f.png"), p, device = grDevices::png, type = "cairo",
  width = plot_width_mm, height = plot_height_mm, units = "mm", dpi = 600, bg = "white"
)

if (interactive()) {
  print(p)
}
