# ============================================================
# Fig. 2a | GO and KEGG enrichment bubble plot
# Repository script for the manuscript figure
# Run from the repository root or the Figure_2 folder.
# Input:  Figure_2/Fig2a_DEG_Functional_Enrichment.csv
# Output: Figure_2/results/Fig2a.pdf and Fig2a.png
# ============================================================

library(ggplot2)
library(dplyr)
library(readr)
library(ggh4x)

# 1. Input and output ------------------------------------------------------

# Run from the repository root or directly from the Figure_2 folder.
working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

figure_dir <- if (basename(working_dir) == "Figure_2") {
  working_dir
} else {
  file.path(working_dir, "Figure_2")
}

data_file <- file.path(figure_dir, "Fig2a_DEG_Functional_Enrichment.csv")
output_dir <- file.path(figure_dir, "results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop(
    "CSV file not found: ", data_file,
    "\nRun this script from the repository root or the Figure_2 folder."
  )
}

# Figure dimensions in millimetres.
fig_width  <- 130
fig_height <- 197

# 2. Display order ---------------------------------------------------------

group_order <- c(
  "Microplastic_Low_vs_Control",
  "Microplastic_High_vs_Control",
  "Plastisphere_Low_vs_Control",
  "Plastisphere_High_vs_Control"
)

group_labels <- c(
  "Microplastic_Low_vs_Control"  = "Microplastic-Low vs Control",
  "Microplastic_High_vs_Control" = "Microplastic-High vs Control",
  "Plastisphere_Low_vs_Control"  = "Plastisphere-Low vs Control",
  "Plastisphere_High_vs_Control" = "Plastisphere-High vs Control"
)

category_order <- c("BP", "CC", "MF", "KEGG")

category_cols <- c(
  BP   = "#B8DCE2",
  CC   = "#DCC8D9",
  MF   = "#F2DDA8",
  KEGG = "#91C18C"
)

# Top-to-bottom order within each category, following the manuscript figure.
term_order <- c(
  # BP
  "cell cycle process",
  "reproductive process",
  "mitotic cell cycle process",
  "meiotic cell cycle process",
  "developmental process involved in reproduction",
  "chromosome organization",
  "chromosome organization involved in meiotic cell cycle",
  "multicellular organismal reproductive process",
  "glucan metabolic process",
  "defense response to bacterium",
  "disruption of anatomical structure in another organism",
  "disruption of cell in another organism",
  "killing of cells of another organism",
  "phosphorylation",
  "intracellular signal transduction",
  "meiosis I cell cycle process",
  "synaptonemal complex organization",
  "synaptonemal complex assembly",
  "microtubule cytoskeleton organization",
  "organelle organization",
  "multicellular organismal process",
  "regulation of neuron projection development",

  # CC
  "extracellular region",
  "external encapsulating structure",
  "extracellular matrix",
  "collagen-containing extracellular matrix",
  "collagen trimer",
  "lysosome",
  "lytic vacuole",
  "other organism part",
  "other organism cell membrane",
  "phagophore assembly site membrane",
  "cullin-RING ubiquitin ligase complex",
  "ubiquitin ligase complex",
  "SCF ubiquitin ligase complex",
  "Cul4-RING E3 ubiquitin ligase complex",
  "serine/threonine protein kinase complex",
  "protein kinase complex",
  "protein-DNA complex",
  "kinetochore",
  "central element",

  # MF
  "chitin binding",
  "cysteine-type endopeptidase regulator activity involved in apoptotic process",
  "peptidase activator activity involved in apoptotic process",
  "cysteine-type endopeptidase activator activity involved in apoptotic process",
  "cellulase activity",
  "structural constituent of chromatin",

  # KEGG
  "Cell cycle",
  "Purine metabolism",
  "Starch and sucrose metabolism",
  "Cholesterol metabolism",
  "Pyrimidine metabolism",
  "Nucleotide metabolism",
  "Glycolysis / Gluconeogenesis",
  "Protein digestion and absorption",
  "Primary immunodeficiency",
  "Antigen processing and presentation",
  "Apoptosis",
  "Necroptosis",
  "Lysosome",
  "Hedgehog signaling pathway",
  "FoxO signaling pathway",
  "AMPK signaling pathway",
  "Small cell lung cancer",
  "Acute myeloid leukemia",
  "Adherens junction",
  "Sulfur relay system",
  "Tyrosine metabolism"
)

# 3. Read and check the source data ---------------------------------------

data <- readr::read_csv(data_file, show_col_types = FALSE)

required_cols <- c(
  "Comparison", "Term_Description", "Category", "Gene_Count", "Adjusted_P"
)

if (!all(required_cols %in% names(data))) {
  stop("The input CSV does not contain all required columns.")
}

unknown_terms <- setdiff(unique(data$Term_Description), term_order)
if (length(unknown_terms) > 0) {
  stop(
    "Terms missing from term_order: ",
    paste(unknown_terms, collapse = "; ")
  )
}

data <- data %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = group_order,
      labels = unname(group_labels[group_order])
    ),
    Category = factor(Category, levels = category_order),
    Term_Description = factor(Term_Description, levels = rev(term_order)),
    Gene_Count = as.numeric(Gene_Count),
    Adjusted_P = as.numeric(Adjusted_P)
  )

# Adjusted_P is read from the supplied table as-is.
# BH adjustment is not recalculated in this plotting script.
p_legend_title <- "Adjusted P"

# Invisible points create the coloured square Category legend.
category_key <- tibble::tibble(
  Category = factor(category_order, levels = category_order),
  Comparison = factor(
    rep(unname(group_labels[group_order[1]]), length(category_order)),
    levels = unname(group_labels[group_order])
  ),
  Term_Description = factor(
    c(
      "cell cycle process",
      "extracellular region",
      "chitin binding",
      "Cell cycle"
    ),
    levels = rev(term_order)
  )
)

# 4. Plot ------------------------------------------------------------------

category_strips <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = unname(category_cols),
    colour = "white"
  ),
  text_y = ggh4x::elem_list_text(
    colour = "black",
    face = "bold",
    size = 7
  )
)

p <- ggplot(data, aes(x = Comparison, y = Term_Description)) +
  geom_point(
    aes(size = Gene_Count^1.5, fill = Adjusted_P),
    shape = 21,
    colour = "black",
    stroke = 0.12
  ) +
  geom_point(
    data = category_key,
    aes(x = Comparison, y = Term_Description, shape = Category),
    inherit.aes = FALSE,
    alpha = 0,
    size = 3,
    show.legend = TRUE
  ) +
  scale_fill_gradientn(
    colours = c("#009A63", "#E4AD77", "#EB5A56"),
    values = scales::rescale(c(0, 0.023, 0.05)),
    limits = c(0, 0.05),
    breaks = c(0.01, 0.02, 0.03, 0.04),
    name = p_legend_title
  ) +
  scale_size_continuous(
    range = c(1.8, 6),
    limits = c(4^1.5, 218^1.5),
    breaks = c(4, 30, 100, 218)^1.5,
    labels = c("4", "30", "100", "218"),
    name = "Gene number"
  ) +
  scale_shape_manual(
    name = "Category",
    values = setNames(rep(22, length(category_order)), category_order),
    breaks = category_order
  ) +
  scale_y_discrete(position = "right", drop = TRUE) +
  facet_grid2(
    rows = vars(Category),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    strip = category_strips
  ) +
  guides(
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      barheight = grid::unit(22, "mm"),
      barwidth = grid::unit(5, "mm")
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        alpha = 1,
        fill = unname(category_cols),
        colour = unname(category_cols),
        size = 4
      )
    ),
    size = guide_legend(
      order = 3,
      override.aes = list(
        fill = "#202020",
        colour = "#202020"
      )
    )
  ) +
  labs(x = NULL, y = NULL, tag = "(a)") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    panel.spacing.y = grid::unit(0, "mm"),

    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 7
    ),

    axis.text.x = element_text(
      colour = "black",
      size = 6.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y.right = element_text(
      colour = "black",
      size = 6.2,
      hjust = 0,
      margin = margin(l = 1, unit = "mm")
    ),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.ticks.length.y = grid::unit(0.7, "mm"),
    axis.line = element_blank(),

    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "left",
    legend.title = element_text(size = 7, face = "plain"),
    legend.text = element_text(size = 6.5),
    legend.key = element_blank(),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.spacing.y = grid::unit(3, "mm"),

    plot.tag = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 10
    ),
    plot.tag.position = c(0, 1),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  ) +
  coord_cartesian(clip = "off")

if (interactive()) {
  print(p)
}

# 5. Export figure files ---------------------------------------------------

ggsave(
  filename = file.path(output_dir, "Fig2a.pdf"),
  device = grDevices::cairo_pdf,
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "mm",
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Fig2a.png"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "mm",
  dpi = 600,
  bg = "white"
)
