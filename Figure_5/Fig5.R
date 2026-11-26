# Fig. 5a-c | Gut 16S alpha diversity
# Run from the repository root or the Figure_5 folder.
# The PDF is written to Figure_5/results/Fig5/.

library(grid)

input_name <- "Gut_16S_Alpha_Diversity.csv"
figure_dir <- if (file.exists(file.path("Figure_5", input_name))) "Figure_5" else "."
output_dir <- file.path(figure_dir, "results", "Fig5")
dat <- read.csv(file.path(figure_dir, input_name),
                fileEncoding = "UTF-8-BOM")
groups <- c("CTR", "MP1", "MP10", "MPP1", "MPP10")
labels <- c("Control", "Microplastic-Low", "Microplastic-High",
            "Plastisphere-Low", "Plastisphere-High")
titles <- c("Shannon", "Pielou", "ACE")
axis_titles <- paste(c("Shannon", "Pielou-e", "ACE"), "index of OTU level")
limits <- c(3.5, 0.8, 210)
steps <- c(0.5, 0.1, 30)

# The source CSV has three independent replicates per group.
group_id <- rep(seq_along(groups), each = 3L)
values <- as.matrix(dat[order(match(trimws(dat$group), groups)), c("shannon", "Pielou_e", "ace")])
means <- rowsum(values, group_id) / 3
sds <- sqrt(rowsum((values - means[group_id, ])^2, group_id) / 2)

# With n = 3 and tr = 0.2, these Welch P values match the unadjusted lincon values.
pairs <- combn(seq_along(groups), 2L)
pairs <- pairs[, order(pairs[1, ], -pairs[2, ])]
variance_a <- sds[pairs[1, ], ]^2 / 3
variance_b <- sds[pairs[2, ], ]^2 / 3
se_squared <- variance_a + variance_b
degrees_freedom <- 2 * se_squared^2 / (variance_a^2 + variance_b^2)
t_value <- (means[pairs[1, ], ] - means[pairs[2, ], ]) / sqrt(se_squared)
p_values <- 2 * pt(-abs(t_value), degrees_freedom)

# Panel layout in mm.
mm <- function(x) unit(x, "mm")
family <- "Times New Roman"
font <- function(size = 25.5, face = "plain")
  gpar(fontfamily = family, fontsize = size, fontface = face)
outline <- gpar(col = "black", fill = NA, lwd = 2.1)
segment <- function(x0, y0, x1, y1, gp = outline)
  grid.segments(mm(x0), mm(y0), mm(x1), mm(y1), gp = gp)
major_grid <- gpar(col = "#EBEBEB", lwd = 2.1)
minor_grid <- gpar(col = "#EBEBEB", lwd = 1.05)
panel_width <- 135
panel_height <- 128.5
title_height <- 17
group_step <- panel_width / 5.2
centers <- (1:5 - 0.4) * group_step
minor_x <- (seq(0.5, 5.5) - 0.4) * group_step
bar_width <- 0.8 * group_step
bar_left <- centers - bar_width / 2
bar_fill <- gpar(fill = c("#B6DCD5", "#FFE0E2", "#FFAAB9", "#BAA5D5", "#8C69B9"), col = NA)
bottom_left <- c("left", "bottom")
cap <- bar_width / 8

# Set the spacing and phase of the horizontal and diagonal hatches.
control_y <- seq(6.4, panel_height, by = 7.7)
pitch <- 7.7 * sqrt(2)
hatch <- expand.grid(row = seq.int(-3L, ceiling((panel_height + bar_width) / pitch)),
                     bar = 2:5)
slope <- c(1, 1, -1, -1)[hatch$bar - 1L]
hatch_x <- bar_left[hatch$bar]
hatch_y <- (ifelse(slope > 0, 2, 6.2) + slope * hatch_x) %% pitch + hatch$row * pitch
hatch_pen <- gpar(col = "black", lwd = 1.26, lineend = "butt")

p_italic <- font(16, "italic")
p_regular <- font(16)
p_width <- grobWidth(textGrob("P", gp = p_italic))
number_widths <- unit.c(grobWidth(textGrob(" = 0.000", gp = p_regular)),
                       grobWidth(textGrob(" < 0.001", gp = p_regular)))

if ("showtext" %in% loadedNamespaces()) showtext::showtext_auto(FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(output_dir, "Fig5abc_Alpha_Diversity_test.pdf"),
          width = 570 / 25.4, height = 200 / 25.4, family = family, bg = "white")
grid.newpage()

tryCatch(for (i in seq_along(titles)) {
  scale <- panel_height / limits[i]
  tick_values <- seq(0, limits[i], by = steps[i])
  ticks <- tick_values * scale
  minor_y <- head(ticks, -1) + diff(ticks) / 2
  pushViewport(viewport(mm(62 + 167.5 * (i - 1)), mm(44), mm(panel_width),
                        mm(panel_height), just = bottom_left))
  segment(0, minor_y, panel_width, minor_y, minor_grid)
  segment(minor_x, 0, minor_x, panel_height, minor_grid)
  segment(0, ticks, panel_width, ticks, major_grid)
  segment(centers, 0, centers, panel_height, major_grid)

  bar_top <- means[, i] * scale
  grid.rect(mm(bar_left), mm(0), mm(bar_width), mm(bar_top), just = bottom_left, gp = bar_fill)
  control_lines <- control_y[control_y < bar_top[1]]
  segment(bar_left[1], control_lines, bar_left[1] + bar_width, control_lines, hatch_pen)
  hatch_top <- bar_top[hatch$bar]
  start <- pmax(0, ifelse(slope > 0, -hatch_y, hatch_y - hatch_top))
  end <- pmin(bar_width, ifelse(slope > 0, hatch_top - hatch_y, hatch_y))
  keep <- start < end
  segment(hatch_x[keep] + start[keep], hatch_y[keep] + slope[keep] * start[keep],
          hatch_x[keep] + end[keep], hatch_y[keep] + slope[keep] * end[keep], hatch_pen)
  grid.rect(mm(bar_left), mm(0), mm(bar_width), mm(bar_top), just = bottom_left, gp = outline)

  lower <- (means[, i] - sds[, i]) * scale
  upper <- (means[, i] + sds[, i]) * scale
  segment(centers, lower, centers, upper)
  segment(centers - cap, lower, centers + cap, lower)
  segment(centers - cap, upper, centers + cap, upper)
  point_y <- values[, i] * scale
  grid.circle(mm(centers[group_id] + c(-cap, 0, cap)), mm(point_y),
              r = mm(0.72), gp = gpar(fill = "white", col = "black", lwd = 0.75))

  significant <- which(p_values[, i] < 0.05)
  if (length(significant)) {
    left <- centers[pairs[1, significant]]
    right <- centers[pairs[2, significant]]
    y <- max(c(94.3, 85.3, 0)[i], upper, point_y + 2.5) + 6.9 * (seq_along(significant) - 1)
    grid.polyline(mm(c(rbind(left, left, right, right))), mm(c(rbind(y - 1.5, y, y, y - 1.5))),
                  id.lengths = rep(4L, length(y)), gp = outline)
    p <- p_values[significant, i]
    small <- p < 0.001
    suffix <- ifelse(small, " < 0.001", sprintf(" = %.3f", p))
    text_x <- mm((left + right) / 2) + (p_width - number_widths[1L + small]) / 2
    text_rows <- rep(seq_along(y), each = 2L)
    grid.text(c(rbind("P", suffix)), text_x[text_rows], mm(y[text_rows] + 1.15),
              hjust = c(1, 0), vjust = 0, gp = font(16, c("italic", "plain")))
  }

  grid.rect(mm(0), mm(0), mm(panel_width), mm(panel_height + title_height),
            just = bottom_left, gp = outline)
  segment(0, panel_height, panel_width, panel_height)
  title_y <- mm(panel_height + title_height / 2)
  grid.text(titles[i], mm(panel_width / 2), title_y, gp = font(32.5, "bold"))
  grid.text(paste0("(", letters[i], ")"), mm(2), title_y, just = "left", gp = font(36))

  segment(-1.5, ticks, 0, ticks)
  tick_format <- c("%.1f", "%.1f", "%.0f")[i]
  grid.text(sprintf(tick_format, tick_values), mm(-2.4), mm(ticks), just = "right", gp = font())
  grid.text(axis_titles[i], mm(c(-19.5, -19.5, -21.7)[i]), mm(panel_height / 2),
            rot = 90, gp = font(29, "bold"))
  segment(centers, 0, centers, -1.5)
  grid.text(labels, mm(centers), mm(-0.8), just = c("right", "top"), rot = 30, gp = font())
  popViewport()
}, finally = invisible(dev.off()))
