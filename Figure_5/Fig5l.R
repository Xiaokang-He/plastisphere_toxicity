rm(list = ls())

library(ggplot2)
library(grid)

# Read the network node and edge tables from the repository
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

# Build an unweighted, undirected adjacency matrix for each group
build_adj <- function(group) {
  node <- nodes$Node_Name[nodes$Subplot == group]
  edge <- edges[edges$Subplot == group, ]

  adj <- matrix(
    0, length(node), length(node),
    dimnames = list(node, node)
  )

  i <- match(edge$Source_Node_Name, node)
  j <- match(edge$Target_Node_Name, node)

  adj[cbind(i, j)] <- 1
  adj[cbind(j, i)] <- 1

  adj
}

# Calculate natural connectivity from the eigenvalues of the adjacency matrix
nc <- function(adj) {
  lambda <- eigen(adj, symmetric = TRUE, only.values = TRUE)$values
  m <- max(lambda)
  m + log(mean(exp(lambda - m)))
}

# Randomly remove nodes and calculate network robustness metrics
calc_network <- function(group) {
  adj <- build_adj(group)
  n <- nrow(adj)

  result <- vector("list", 100)

  for (i in 0:99) {
    if (i == 0) {
      sub <- adj
    } else {
      remove_n <- floor(n * i / 100)
      keep <- sample.int(n, n - remove_n)
      sub <- adj[keep, keep, drop = FALSE]
    }

    result[[i + 1]] <- data.frame(
      remove_node = i,
      variable = c("natural_connectivity", "average_degree"),
      values = c(
        nc(sub),
        mean(rowSums(sub))
      ),
      network = group
    )
  }

  do.call(rbind, result)
}

# Calculate robustness metrics for all five groups
set.seed(123)
dat <- do.call(rbind, lapply(groups, calc_network))

dat$network <- factor(dat$network, levels = groups)

# Retain natural connectivity for the final plot
plot_data <- subset(dat, variable == "natural_connectivity")

# Draw the robustness plot
p <- ggplot(plot_data, aes(remove_node, values)) +
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

print(p)

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
