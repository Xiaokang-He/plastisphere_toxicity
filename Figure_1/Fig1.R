library(ggplot2)
library(grid)

data_file <- "Fig1a_Clam_Physiological_and_Biochemical_Endpoints.csv"
input_file <- if (file.exists(data_file)) data_file else file.path("Figure_1", data_file)
if (!file.exists(input_file)) stop("Cannot find input data: ", data_file)
output_dir <- dirname(input_file)

dat <- read.csv(input_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(dat) <- sub("^lysozyme$", "LZM", sub("^D_Lac.*", "D_Lac", names(dat)))

groups <- c(CTR = "Control", MP1 = "Microplastic_Low", MP10 = "Microplastic_High",
            MPP1 = "Plastisphere_Low", MPP10 = "Plastisphere_High")
sample_group <- if (is.null(dat$group)) sub("_[0-9]+$", "", dat$sample) else groups[trimws(dat$group)]
dat$group <- factor(sample_group, levels = groups)
colors <- setNames(c("#b6dcd5", "#ffe0e2", "#ffaab9", "#baa5d5", "#8c69b9"), groups)

variables <- c("Siphonage_efficiency", "SOD", "CAT", "GSH", "GSH_PX", "MDA", "AchE", "LZM", "D_Lac", "DAO")
titles <- c("Siphoning", "SOD", "CAT", "GSH", "GSH_PX", "MDA", "AchE", "Lysozyme", "D_Lac", "DAO")
y_labels <- c("Siphonage efficiency (ml/clam/h)", "SOD activity (U/mgprot)", "CAT activity (U/mgprot)",
              "GSH content (\u03bcmol/gprot)", "GSH-Px activity (U/mgprot)", "MDA content (nmol/mgprot)",
              "AchE activity (U/mgprot)", "Lysozyme activity (U/mgprot)", "D-Lac content (\u03bcmol/mgprot)",
              "DAO activity (U/g prot)")
y_steps <- c(15, 150, 1.5, 5, 3, 0.6, 0.7, 7, 0.08, 7)
digits <- c(0, 0, 1, 0, 0, 1, 1, 0, 2, 0)

column_pitch <- 135
row_pitch <- 198
panel_width <- 89.277
panel_height <- 136.515
title_height <- 18.238
page_width <- 44 + 4 * column_pitch + panel_width + 8
page_height <- 46 + row_pitch + panel_height + title_height + 4
font <- "Times New Roman"
axis_text <- gpar(fontfamily = font, fontsize = 28.4)
axis_title <- gpar(fontfamily = font, fontsize = 28.4, fontface = "bold")
title_text <- gpar(fontfamily = font, fontsize = 36.59, fontface = "bold")
line <- gpar(lwd = 1.788 / 0.75)
line_width <- line$lwd / ggplot2::.pt
mm <- function(x) unit(x, "mm")

draw_panel <- function(i, left, bottom) {
  obs <- na.omit(data.frame(group = dat$group, value = as.numeric(dat[[variables[i]]])))
  obs$x <- as.integer(obs$group)
  obs$jitter_x <- obs$x + ave(obs$x, obs$group, FUN = function(x) 0.1 * (seq_along(x) - (length(x) + 1) / 2))
  values <- split(obs$value, obs$group)
  bars <- data.frame(group = factor(groups, levels = groups), x = 1:5, angle = c(0, 45, 45, 135, 135),
                     mean = vapply(values, mean, 0.0), sd = vapply(values, sd, 0.0))

  pairs <- WRS2::lincon(value ~ group, data = obs, tr = 0.2, method = "none")
  idx <- pairs$comp[, 1:2, drop = FALSE]
  sig <- data.frame(x1 = match(pairs$fnames[idx[, 1]], groups),
                    x2 = match(pairs$fnames[idx[, 2]], groups), p = pairs$comp[, 6])
  sig <- subset(sig, p < 0.05)
  sig <- sig[order(sig$x1, -sig$x2), ]
  right_edge <- rep(-Inf, nrow(sig))
  sig$lane <- integer(nrow(sig))
  for (j in seq_len(nrow(sig))) {
    lane <- which(right_edge < sig$x1[j])[1]
    right_edge[lane] <- sig$x2[j]
    sig$lane[j] <- lane
  }

  step <- y_steps[i]
  data_top <- max(obs$value, bars$mean + bars$sd)
  ymax <- max(4 * step, ceiling(data_top / step) * step)
  lanes <- max(0, sig$lane)
  font_sizes <- seq(10, 18, 0.5)
  spacing <- font_sizes * 25.4 / 72 * 0.9 + 0.60
  room <- (1 - data_top / ymax) * panel_height
  p_size <- max(10, font_sizes[0.75 + lanes * spacing <= room])
  bracket_step <- p_size * 25.4 / 72 * 0.9 + 0.60
  if (lanes > 0) {
    needed <- data_top / (1 - (0.75 + lanes * bracket_step) / panel_height)
    ymax <- max(ymax, ceiling(needed / (step / 4)) * (step / 4))
  }
  sig$y <- data_top + (sig$lane - 1) * bracket_step * ymax / panel_height + 0.50 * ymax / panel_height
  sig$hook <- sig$y - 0.012 * ymax
  ticks <- seq(0, ymax, by = step)

  p <- ggplot() +
    ggpattern::geom_col_pattern(aes(x, mean, fill = group, pattern_fill = group, pattern_angle = angle),
      data = bars, width = if (i == 1) 0.9 else 0.8, colour = "black", linewidth = line_width,
      pattern = "stripe", pattern_spacing = 0.12, pattern_density = 0.5, pattern_colour = "black") +
    geom_errorbar(aes(x, ymin = mean - sd, ymax = mean + sd), data = bars, width = 0.2, linewidth = line_width) +
    geom_point(aes(jitter_x, value), data = obs, shape = 21, fill = "white", size = 2, stroke = 0.45) +
    geom_segment(aes(x1, y, xend = x2, yend = y), data = sig, linewidth = line_width) +
    geom_segment(aes(x1, y, xend = x1, yend = hook), data = sig, linewidth = line_width) +
    geom_segment(aes(x2, y, xend = x2, yend = hook), data = sig, linewidth = line_width) +
    scale_fill_manual(values = colors) + ggpattern::scale_pattern_fill_manual(values = colors) +
    ggpattern::scale_pattern_angle_identity() +
    scale_x_continuous(breaks = 1:5, minor_breaks = seq(0.5, 5.5, 0.5)) +
    scale_y_continuous(breaks = ticks, minor_breaks = ticks + step / 2) +
    coord_cartesian(xlim = c(0.4, 5.6), ylim = c(0, ymax), expand = FALSE, clip = "off") +
    theme_void() + theme(legend.position = "none",
      panel.grid.major = element_line(colour = "#EBEBEB", linewidth = line_width),
      panel.grid.minor = element_line(colour = "#EBEBEB", linewidth = line_width * 0.886 / 1.788),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = line_width))

  chart <- ggplotGrob(p)
  pushViewport(viewport(mm(left), mm(bottom), mm(panel_width), mm(panel_height), just = c("left", "bottom")))
  grid.draw(chart$grobs[[which(chart$layout$name == "panel")]])
  popViewport()

  x_pos <- function(x) mm(left + (x - 0.4) / 5.2 * panel_width)
  y_pos <- function(y) mm(bottom + y / ymax * panel_height)
  grid.segments(mm(left - 1.8), y_pos(ticks), mm(left), y_pos(ticks), gp = line)
  grid.text(formatC(ticks, format = "f", digits = digits[i]), mm(left - 4), y_pos(ticks), just = "right", gp = axis_text)
  grid.text(y_labels[i], mm(left - 31), mm(bottom + panel_height / 2), rot = 90, gp = axis_title)
  grid.segments(x_pos(1:5), mm(bottom), x_pos(1:5), mm(bottom - 1.8), gp = line)
  grid.text(gsub("_", "-", groups), x_pos(1:5), mm(bottom - 0.1), just = c("right", "top"), rot = 30, gp = axis_text)

  if (nrow(sig) > 0) {
    numbers <- ifelse(round(sig$p, 3) >= 0.05, sprintf("%.4f", sig$p), sprintf("%.3f", sig$p))
    labels <- ifelse(sig$p < 0.001, 'italic(P)<"0.001"', sprintf('italic(P)=="%s"', numbers))
    grid.text(parse(text = labels), x_pos((sig$x1 + sig$x2) / 2), y_pos(sig$y) + mm(0.35),
              just = "bottom", gp = gpar(fontfamily = font, fontsize = p_size))
  }

  title_x <- mm(left + panel_width / 2)
  title_y <- mm(bottom + panel_height + title_height / 2)
  grid.rect(title_x, title_y, mm(panel_width), mm(title_height), gp = gpar(fill = "white", lwd = line$lwd))
  grid.text(sprintf("(%s) %s", letters[i], titles[i]), title_x, title_y, gp = title_text)
}

draw_figure <- function() {
  grid.newpage()
  for (i in seq_along(variables)) {
    left <- 44 + ((i - 1) %% 5) * column_pitch
    bottom <- 46 + (i <= 5) * row_pitch
    draw_panel(i, left, bottom)
  }
}

grDevices::cairo_pdf(file.path(output_dir, "Fig1.pdf"),
  width = page_width / 25.4, height = page_height / 25.4, family = font)
draw_figure()
invisible(grDevices::dev.off())

grDevices::jpeg(file.path(output_dir, "Fig1.jpg"),
  width = page_width, height = page_height, units = "mm",
  res = 300, quality = 100, type = "cairo")
draw_figure()
invisible(grDevices::dev.off())
