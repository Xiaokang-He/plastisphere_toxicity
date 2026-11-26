library(dplyr)
library(ggplot2)
library(ggpattern)
library(cowplot)

# Read the exported network node and edge tables
nodes <- read.csv(
  "Figure_5/Fig5g_Network_Topology_Raw_Nodes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

edges <- read.csv(
  "Figure_5/Fig5g_Network_Topology_Raw_Edges.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Group order, colours and stripe angles
group_order <- c(
  "Control", "Microplastic-Low", "Microplastic-High",
  "Plastisphere-Low", "Plastisphere-High"
)

group_colors <- c(
  "Control" = "#b6dcd5",
  "Microplastic-Low" = "#ffe0e2",
  "Microplastic-High" = "#ffaab9",
  "Plastisphere-Low" = "#baa5d5",
  "Plastisphere-High" = "#8c69b9"
)

group_angles <- c(
  "Control" = 0,
  "Microplastic-Low" = 45,
  "Microplastic-High" = 45,
  "Plastisphere-Low" = 135,
  "Plastisphere-High" = 135
)

# Calculate the number of distinct nodes in each group
node_data <- nodes %>%
  filter(Subplot %in% group_order) %>%
  group_by(Subplot) %>%
  summarise(
    Number_of_nodes = n_distinct(Node_Name),
    .groups = "drop"
  )

# Count unique undirected edges and calculate the positive-edge ratio
edge_data <- edges %>%
  filter(Subplot %in% group_order) %>%
  mutate(
    Node_1 = pmin(Source_Node_Name, Target_Node_Name),
    Node_2 = pmax(Source_Node_Name, Target_Node_Name)
  ) %>%
  distinct(Subplot, Node_1, Node_2, .keep_all = TRUE) %>%
  group_by(Subplot) %>%
  summarise(
    Number_of_edges = n(),
    Positive_correlation_ratio =
      mean(Correlation_Type == "Positive") * 100,
    .groups = "drop"
  )

# Combine the topology metrics
topology_data <- node_data %>%
  left_join(edge_data, by = "Subplot") %>%
  mutate(
    group = factor(Subplot, levels = group_order),
    Average_degree = 2 * Number_of_edges / Number_of_nodes
  ) %>%
  arrange(group)

print(topology_data)

# Create a consistent bar-plot style for all four metrics
make_plot <- function(metric, y_label, y_max, y_breaks) {

  ggplot(
    topology_data,
    aes(x = group, y = .data[[metric]], fill = group)
  ) +
    geom_col_pattern(
      aes(
        pattern = group,
        pattern_fill = group,
        pattern_angle = group
      ),
      pattern_spacing = 0.16,
      pattern_size = 0.2,
      pattern_density = 0.5,
      pattern_colour = "black",
      colour = "black",
      linewidth = 0.45,
      width = 0.72
    ) +
    scale_pattern_manual(
      values = setNames(rep("stripe", 5), group_order)
    ) +
    scale_pattern_angle_manual(values = group_angles) +
    scale_pattern_fill_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(limits = group_order, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = y_label) +
    theme_bw(
      base_size = 10,
      base_family = "Times New Roman"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        size = 9.5,
        colour = "black",
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        size = 9.5,
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 11,
        face = "bold",
        colour = "black",
        margin = margin(r = 5)
      ),
      axis.title.x = element_blank(),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.4
      ),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.grid.major = element_line(
        colour = "#D0D0D0",
        linewidth = 0.35
      ),
      panel.grid.minor = element_line(
        colour = "#E8E8E8",
        linewidth = 0.25
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.45
      ),
      plot.margin = margin(
        t = 20,
        r = 4,
        b = 28,
        l = 4,
        unit = "pt"
      )
    )
}

# Plot the four network topology metrics
p_h <- make_plot(
  "Positive_correlation_ratio",
  "Positive correlation ratio (%)",
  100,
  seq(0, 100, 25)
)

p_i <- make_plot(
  "Number_of_nodes",
  "Number of nodes",
  400,
  seq(0, 400, 100)
)

p_j <- make_plot(
  "Number_of_edges",
  "Number of edges",
  8000,
  seq(0, 8000, 2000)
)

p_k <- make_plot(
  "Average_degree",
  "Average degree",
  80,
  seq(0, 80, 20)
)

# Combine panels h-k
combined_plot <- plot_grid(
  p_h, p_i, p_j, p_k,
  ncol = 4,
  labels = c("(h)", "(i)", "(j)", "(k)"),
  label_size = 12,
  label_fontface = "plain",
  label_fontfamily = "Times New Roman",
  label_x = 0.12,
  label_y = 0.98,
  hjust = 0,
  vjust = 1,
  align = "hv",
  axis = "tblr"
)

if (interactive()) print(combined_plot)

# Export the combined figure
dir.create("Figure_5/results", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "Figure_5/results/Fig5h-k_Network_Topology.pdf",
  plot = combined_plot,
  device = cairo_pdf,
  width = 390,
  height = 130,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)
