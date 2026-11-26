library(ggplot2)
library(grid)

# Read the natural connectivity series for the manuscript figure.
dat <- read.csv(
  "Figure_5/Fig5l_Network_Robustness.csv",
  fileEncoding = "UTF-8-BOM", check.names = FALSE
)

# Define the group order and colour palette
groups <- c(
  "Control", "Microplastic-Low", "Microplastic-High",
  "Plastisphere-Low", "Plastisphere-High"
)

colors <- c(
  "Control" = "#b6dcd5",
  "Microplastic-Low" = "#ffe0e2",
  "Microplastic-High" = "#ffaab9",
  "Plastisphere-Low" = "#baa5d5",
  "Plastisphere-High" = "#8c69b9"
)

dat$network <- factor(dat$network, levels = groups)

# Draw the robustness plot
p <- ggplot(dat, aes(remove_node, values)) +
  geom_smooth(
    aes(color = network),
    se = FALSE,
    show.legend = TRUE
  ) +
  geom_line(
    aes(group = network, color = network),
    linewidth = 0.5,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = network),
    color = "black",
    alpha = 0.6,
    size = 3,
    shape = 21,
    stroke = 0.1,
    show.legend = TRUE
  ) +
  scale_color_manual(
    name = "Network",
    values = colors,
    breaks = groups
  ) +
  scale_fill_manual(
    name = "Network",
    values = colors,
    breaks = groups
  ) +
  scale_x_continuous(
    name = "Removal node ratio (%)",
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Natural connectivity of the remaining network",
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, 100),
    ylim = c(0, 100),
    expand = FALSE
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.width = unit(14, "mm"),
    legend.key.height = unit(7, "mm"),
    legend.spacing.y = unit(2, "mm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5
    ),
    plot.title = element_blank(),
    plot.margin = margin(5, 8, 5, 5)
  )

if (interactive()) print(p)

# Export the figure and calculated data
ggsave(
  "Figure_5/Robustness.pdf",
  plot = p,
  width = 160,
  height = 105,
  units = "mm"
)

write.csv(
  dat,
  "Figure_5/Robustness_data.csv",
  row.names = FALSE
)
