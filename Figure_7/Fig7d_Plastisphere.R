# Figure 7d | Plastisphere
# Reproduces the Pearson-correlation and Mantel-association plot.


rm(list = ls())

library(ggplot2)
library(dplyr)
library(vegan)

set.seed(123)

# Run from the repository root or the Figure_7 folder.
figure_dir <- if (file.exists(file.path("Figure_7", "Fig7d_Gut_Microbiota_Metrics_Plastisphere.csv"))) {
  "Figure_7"
} else {
  "."
}


# 1. Read data ------------------------------------------------

biochem_data <- read.csv(
  file.path(figure_dir, "Fig7d_Host_Biochemical_Biomarkers_Plastisphere.csv"),
  row.names = 1,
  check.names = FALSE
)

alpha_data <- read.csv(
  file.path(figure_dir, "Fig7d_Gut_Microbiota_Metrics_Plastisphere.csv"),
  row.names = 1,
  check.names = FALSE
)

sample_id <- intersect(
  rownames(biochem_data),
  rownames(alpha_data)
)

biochem_data <- biochem_data[sample_id, , drop = FALSE]
alpha_data <- alpha_data[sample_id, 1:4, drop = FALSE]

biomarker_names <- c(
  "SOD", "CAT", "GSH", "GSH-Px", "MDA",
  "AChE", "Lysozyme", "D-Lac", "DAO"
)

alpha_names <- c(
  "Shannon", "Pielou", "ACE", "Average Degree"
)

colnames(biochem_data) <- biomarker_names
colnames(alpha_data) <- alpha_names

biochem_data[] <- lapply(biochem_data, as.numeric)
alpha_data[] <- lapply(alpha_data, as.numeric)

# 2. Pearson correlation and P values -------------------------

nvar <- length(biomarker_names)

cor_mat <- cor(
  biochem_data,
  method = "pearson",
  use = "pairwise.complete.obs"
)

p_mat <- matrix(
  NA_real_,
  nrow = nvar,
  ncol = nvar,
  dimnames = list(biomarker_names, biomarker_names)
)

for (i in seq_len(nvar - 1)) {
  for (j in (i + 1):nvar) {
    
    ok <- complete.cases(
      biochem_data[[i]],
      biochem_data[[j]]
    )
    
    p_value <- cor.test(
      biochem_data[[i]][ok],
      biochem_data[[j]][ok],
      method = "pearson"
    )$p.value
    
    p_mat[i, j] <- p_value
    p_mat[j, i] <- p_value
  }
}

# 3. Pairwise Mantel tests ------------------------------------

pairs <- expand.grid(
  alpha = alpha_names,
  biomarker = biomarker_names,
  stringsAsFactors = FALSE
)

mantel_df <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  
  a <- pairs$alpha[i]
  b <- pairs$biomarker[i]
  
  ok <- complete.cases(
    alpha_data[[a]],
    biochem_data[[b]]
  )
  
  fit <- vegan::mantel(
    dist(alpha_data[ok, a, drop = FALSE]),
    dist(biochem_data[ok, b, drop = FALSE]),
    method = "pearson",
    permutations = 999
  )
  
  data.frame(
    alpha = a,
    biomarker = b,
    r = as.numeric(fit$statistic),
    p = fit$signif
  )
})) %>%
  mutate(
    r_group = cut(
      r,
      breaks = c(-Inf, 0.25, 0.5, Inf),
      labels = c("< 0.25", "0.25 - 0.5", ">= 0.5"),
      right = FALSE
    ),
    p_group = cut(
      p,
      breaks = c(-Inf, 0.01, 0.05, Inf),
      labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"),
      right = FALSE
    )
  )

# 4. Correlation matrix ---------------------------------------

cell_df <- expand.grid(
  row = seq_len(nvar),
  col = seq_len(nvar)
) %>%
  filter(col > row) %>%
  mutate(
    x = col,
    y = nvar - row + 1,
    r = cor_mat[cbind(row, col)],
    p = p_mat[cbind(row, col)],
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 5. Generate correlation ellipses ----------------------------

theta <- seq(0, 2 * pi, length.out = 100)

ellipse_df <- bind_rows(lapply(seq_len(nrow(cell_df)), function(i) {
  
  r <- cell_df$r[i]
  
  eig <- eigen(
    matrix(c(1, r, r, 1), nrow = 2),
    symmetric = TRUE
  )
  
  xy <- eig$vectors %*%
    diag(sqrt(pmax(eig$values, 0))) %*%
    rbind(cos(theta), sin(theta))
  
  data.frame(
    x = cell_df$x[i] + 0.45 * xy[1, ],
    y = cell_df$y[i] + 0.45 * xy[2, ],
    group = i,
    r = r
  )
}))

# 6. Mantel network nodes -------------------------------------

bio_nodes <- data.frame(
  biomarker = biomarker_names,
  x_to = seq_len(nvar),
  y_to = nvar:1
)

alpha_nodes <- data.frame(
  alpha = alpha_names,
  x_from = c(-3.2, -1.2, 0.8, 2.8),
  y_from = c(6.8, 4.8, 2.8, 0.8)
)

mantel_df <- mantel_df %>%
  left_join(alpha_nodes, by = "alpha") %>%
  left_join(bio_nodes, by = "biomarker")

# 7. Original correlation color palette -----------------------

color_palette <- c(
  colorRampPalette(c("#5EA04C", "white"))(100),
  colorRampPalette(c("white", "#EC6D79"))(100)
)

# 8. Plot ------------------------------------------------------

p <- ggplot() +
  
  # Mantel connections
  geom_curve(
    data = mantel_df,
    aes(
      x = x_from,
      y = y_from,
      xend = x_to,
      yend = y_to,
      colour = p_group,
      linewidth = r_group
    ),
    curvature = 0.18,
    lineend = "round",
    alpha = 0.95
  ) +
  
  # Matrix grid
  geom_tile(
    data = cell_df,
    aes(x = x, y = y),
    width = 1,
    height = 1,
    fill = "white",
    colour = "grey78",
    linewidth = 0.22
  ) +
  
  # Correlation ellipses
  geom_polygon(
    data = ellipse_df,
    aes(
      x = x,
      y = y,
      group = group,
      fill = r
    ),
    colour = "grey45",
    linewidth = 0.12
  ) +
  
  # Significance stars
  geom_text(
    data = cell_df,
    aes(
      x = x,
      y = y,
      label = sig
    ),
    family = "serif",
    size = 4.2,
    fontface = "bold",
    colour = "black"
  ) +
  
  # Alpha diversity nodes
  geom_point(
    data = alpha_nodes,
    aes(x = x_from, y = y_from),
    shape = 21,
    size = 2.2,
    fill = "#FF0066",
    colour = "#3049B8",
    stroke = 0.5
  ) +
  
  # Biomarker nodes
  geom_point(
    data = bio_nodes,
    aes(x = x_to, y = y_to),
    shape = 21,
    size = 2.2,
    fill = "#FF0066",
    colour = "#3049B8",
    stroke = 0.5
  ) +
  
  # Alpha diversity labels
  geom_text(
    data = alpha_nodes,
    aes(
      x = x_from - 0.18,
      y = y_from,
      label = alpha
    ),
    hjust = 1,
    family = "serif",
    size = 3.5
  ) +
  
  # Top biomarker labels
  geom_text(
    data = bio_nodes,
    aes(
      x = x_to,
      y = 10.05,
      label = biomarker
    ),
    angle = 90,
    family = "serif",
    size = 3.5,
    hjust = 0.5,
    vjust = 0.5
  ) +
  
  # Right biomarker labels
  geom_text(
    data = bio_nodes,
    aes(
      x = 9.65,
      y = y_to,
      label = biomarker
    ),
    hjust = 0,
    family = "serif",
    size = 3.5
  ) +
  
  # Panel label
  annotate(
    "text",
    x = -4.65,
    y = 10.25,
    label = "(d)",
    family = "serif",
    size = 5
  ) +
  
  # Pearson correlation color scale
  scale_fill_gradientn(
    colours = color_palette,
    limits = c(-1, 1),
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("-0.8", "-0.4", "0.0", "0.4", "0.8"),
    oob = scales::squish,
    name = "Pearson's r"
  ) +
  
  # Mantel P color scale
  scale_colour_manual(
    values = c(
      "< 0.01" = "#4CAF61",
      "0.01 - 0.05" = "#C9C45A",
      ">= 0.05" = "#D9D9D9"
    ),
    name = "Mantel's p",
    drop = FALSE
  ) +
  
  # Mantel r line width
  scale_linewidth_manual(
    values = c(
      "< 0.25" = 0.45,
      "0.25 - 0.5" = 0.9,
      ">= 0.5" = 1.5
    ),
    name = "Mantel's r",
    drop = FALSE
  ) +
  
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(linewidth = 1)
    ),
    linewidth = guide_legend(
      order = 2,
      override.aes = list(colour = "grey35")
    ),
    fill = guide_colorbar(
      order = 3
    )
  ) +
  
  coord_fixed(
    ratio = 1,
    xlim = c(-5.1, 10.8),
    ylim = c(0.2, 10.6),
    expand = FALSE,
    clip = "off"
  ) +
  
  theme_void(base_family = "serif") +
  
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(
      family = "serif",
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      family = "serif",
      size = 9
    ),
    plot.margin = margin(8, 14, 8, 8),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    )
  )

# 9. Display --------------------------------------------------

print(p)

# 10. Export PDF ----------------------------------------------

ggsave(
  file.path(figure_dir, "Fig7d_Plastisphere_ellipse_final.pdf"),
  plot = p,
  width = 220,
  height = 150,
  units = "mm",
  device = grDevices::pdf,
  family = "serif",
  useDingbats = FALSE,
  bg = "white"
)