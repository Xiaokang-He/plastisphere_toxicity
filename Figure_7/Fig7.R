# Fig. 7a-c | Gut microbial-host functional association sets
# Run from the repository root or the Figure_7 folder.
# Coordinates were measured from the reference figure (792 x 600 px).
# Node counts and labels are read from the three panel-specific CSV files.

library(grid)

figure_dir <- if (file.exists(file.path("Figure_7", "Fig7a_Shared_Associations.csv"))) {
  "Figure_7"
} else {
  "."
}

shared <- read.csv(file.path(figure_dir, "Fig7a_Shared_Associations.csv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE)
mp <- read.csv(file.path(figure_dir, "Fig7b_Microplastic_Specific_Associations.csv"),
               fileEncoding = "UTF-8-BOM", check.names = FALSE)
mpp <- read.csv(file.path(figure_dir, "Fig7c_Plastisphere_Specific_Associations.csv"),
                fileEncoding = "UTF-8-BOM", check.names = FALSE)

make_nodes <- function(d, scope) {
  is_mp <- scope == "Microplastic"
  scope_label <- if (is_mp) "MP" else "MPP"
  fill <- if (is_mp) "#FFB3C0" else "#9878C0"

  microbial <- data.frame(
    subpanel = d$Figure_Subpanel,
    side = "Microbe",
    term_id = d$Gut_Microbiota_KEGG_Pathway_ID,
    name = d$Gut_Microbiota_KEGG_Pathway,
    n = d[[paste0("N_KOs_in_Gut_Microbiota_Pathway_", scope)]],
    shape = 24, scope = scope_label, fill = fill,
    stringsAsFactors = FALSE
  )
  host <- data.frame(
    subpanel = d$Figure_Subpanel,
    side = "Host",
    term_id = d$Host_Function_Term_ID,
    name = d$Host_Function_Term_Name,
    n = d[[paste0("N_Host_DEGs_or_DEPs_", scope)]],
    shape = ifelse(grepl("DEPs", d$Host_Omics_Set_Type), 22, 21),
    scope = scope_label, fill = fill,
    stringsAsFactors = FALSE
  )
  rbind(microbial, host)
}

nodes <- unique(rbind(
  make_nodes(shared, "Microplastic"),
  make_nodes(shared, "Plastisphere"),
  make_nodes(mp, "Microplastic"),
  make_nodes(mpp, "Plastisphere")
))

# Node centres measured in pixels from the reference figure.
positions <- read.table(text = "
subpanel side term_id x_px y_px
a1 Microbe map02020 122.5 79
a1 Microbe map02010 122.5 109
a1 Host GO:0007155 250.5 80
a1 Host GO:0031012 325.5 80
a1 Host GO:0009986 250.5 110
a1 Host GO:0097060 325.5 110
b1 Microbe map00515 88 211
b1 Microbe map02020 88 257
b1 Host map04814 227.5 220
b1 Host map04820 301 220
b1 Host map04022 269.5 257
b2 Microbe map03013 455 211
b2 Microbe map04710 455 257
b2 Host GO:0005911 603 220
b2 Host map00590 678 220
b2 Host map00770 639 257
c1 Microbe map00190 90 359
c1 Microbe map02010 90 389
c1 Microbe map02020 90 419
c1 Host map04620 239 359
c1 Host map04625 315 359
c1 Host map04530 239 389
c1 Host GO:0098631 315 389
c1 Host GO:0061919 239 419
c1 Host map04810 315 419
c2 Microbe map04120 452 359
c2 Microbe map04137 452 419
c2 Host map04144 603 359
c2 Host map04137 700 359
c2 Host GO:0009056 603 389
c2 Host GO:0016788 700 389
c2 Host map04141 603 419
c2 Host map04120 700 419
c3 Microbe map00190 90 505
c3 Microbe map02010 90 534
c3 Microbe map04120 90 563
c3 Host map04145 232 505
c3 Host map04210 299 505
c3 Host map04148 232 534
c3 Host map04140 299 534
c3 Host GO:0008219 232 563
c3 Host map04217 299 563
c4 Microbe map00515 452 505
c4 Microbe map00190 452 534
c4 Microbe map02010 452 563
c4 Host map00564 603 505
c4 Host map04151 700 505
c4 Host GO:1990204 650 534
c4 Host map00620 603 563
c4 Host map00071 700 563
", header = TRUE, stringsAsFactors = FALSE)

key <- function(d) paste(d$subpanel, d$side, d$term_id, sep = "|")
idx <- match(key(nodes), key(positions))
nodes$x_px <- positions$x_px[idx]
nodes$y_px <- positions$y_px[idx]

label_map <- c(
  "ABC transporters" = "ABC transporters",
  "Two-component system" = "Two-component\nsystem",
  "Mannose type O-glycan biosynthesis" = "Mannose-type\nO-glycan biosynth.",
  "Nucleocytoplasmic transport" = "Nucleocytoplasmic\ntransport",
  "Circadian rhythm" = "Circadian rhythm",
  "Motor proteins" = "Motor proteins",
  "Cytoskeleton in muscle cells" = "Cytoskeleton in\nmuscle cells",
  "cGMP-PKG signaling pathway" = "cGMP-PKG\nsignaling",
  "cell-cell junction" = "Cell-cell junction",
  "Arachidonic acid metabolism" = "Arachidonic acid\nmetabolism",
  "Pantothenate and CoA biosynthesis" = "Pantothenate /\nCoA biosynth.",
  "Oxidative phosphorylation" = "Oxidative\nphosphorylation",
  "Toll-like receptor signaling pathway" = "Toll-like receptor\nsignaling",
  "C-type lectin receptor signaling pathway" = "C-type lectin\nreceptor signaling",
  "Tight junction" = "Tight junction",
  "cell adhesion mediator activity" = "Cell adhesion\nmediator activity",
  "process utilizing autophagic mechanism" = "Autophagic\nmechanism",
  "Regulation of actin cytoskeleton" = "Actin cytoskeleton\nregulation",
  "Ubiquitin mediated proteolysis" = "Ubiquitin-mediated\nproteolysis",
  "Mitophagy - animal" = "Mitophagy - animal",
  "Endocytosis" = "Endocytosis",
  "catabolic process" = "Catabolic process",
  "hydrolase activity, acting on ester bonds" = "Ester-bond\nhydrolase",
  "Protein processing in endoplasmic reticulum" = "ER protein\nprocessing",
  "Phagosome" = "Phagosome",
  "Apoptosis" = "Apoptosis",
  "Efferocytosis" = "Efferocytosis",
  "Autophagy - animal" = "Autophagy - animal",
  "cell death" = "Cell death",
  "Necroptosis" = "Necroptosis",
  "Glycerophospholipid metabolism" = "Glycerophospholipid\nmetabolism",
  "PI3K-Akt signaling pathway" = "PI3K-Akt\nsignaling",
  "oxidoreductase complex" = "Oxidoreductase\ncomplex",
  "Pyruvate metabolism" = "Pyruvate\nmetabolism",
  "Fatty acid degradation" = "Fatty acid\ndegradation",
  "cell surface" = "Cell surface",
  "cell adhesion" = "Cell adhesion",
  "extracellular matrix" = "Extracellular\nmatrix",
  "synaptic membrane" = "Synaptic membrane"
)

nodes$label <- unname(label_map[nodes$name])
nodes$label[is.na(nodes$label)] <- nodes$name[is.na(nodes$label)]

# The reference image is 792:600 (1.32). The new canvas preserves this ratio.
page_width_mm <- 180
page_height_mm <- page_width_mm * 600 / 792
px_mm <- page_width_mm / 792
font_family <- "Times New Roman"
label_font_pt <- 5.0
panel_font_pt <- 6.7
subtitle_font_pt <- 4.8
legend_font_pt <- 3.8
outline_padding_mm <- 0.55
shared_gap_mm <- 12.5 * px_mm

x_mm <- function(x) x * px_mm
y_mm <- function(y) y * px_mm
grid_y <- function(y) unit(page_height_mm - y, "mm")

size_from_n <- function(n) {
  z <- (log1p(n) - log1p(3)) / (log1p(122) - log1p(3))
  1.4 + 2.8 * z
}

nodes$x_mm <- x_mm(nodes$x_px)
nodes$y_mm <- y_mm(nodes$y_px)
nodes$size_mm <- size_from_n(nodes$n)

text_bounds <- function(x, y, label, fontsize = label_font_pt) {
  g <- textGrob(label, gp = gpar(fontfamily = font_family,
                                  fontsize = fontsize, lineheight = 0.92))
  w <- convertWidth(grobWidth(g), "mm", valueOnly = TRUE)
  h <- convertHeight(grobHeight(g), "mm", valueOnly = TRUE)
  c(left = x - w / 2, right = x + w / 2,
    top = y - h / 2, bottom = y + h / 2)
}

symbol_bounds <- function(x, y, shape, diameter) {
  if (shape == 21) {
    r <- diameter / 2
    c(left = x - r, right = x + r, top = y - r, bottom = y + r)
  } else if (shape == 22) {
    h <- diameter / 2
    c(left = x - h, right = x + h, top = y - h, bottom = y + h)
  } else {
    w <- diameter * 0.866
    c(left = x - w / 2, right = x + w / 2,
      top = y - diameter * 2 / 3, bottom = y + diameter / 3)
  }
}

fit_ellipse <- function(boxes, padding = outline_padding_mm) {
  boxes <- as.matrix(boxes)
  cx <- mean(range(boxes[, 1:2]))
  cy <- mean(range(boxes[, 3:4]))
  corners <- do.call(rbind, lapply(seq_len(nrow(boxes)), function(i) {
    expand.grid(
      x = boxes[i, 1:2] + c(-padding, padding),
      y = boxes[i, 3:4] + c(-padding, padding)
    )
  }))
  dx <- corners$x - cx
  dy <- corners$y - cy
  objective <- function(log_ratio) {
    ratio <- exp(log_ratio)
    pi * ratio * max((dx / ratio)^2 + dy^2)
  }
  ratio <- exp(optimize(objective, interval = log(c(0.25, 30)))$minimum)
  ry <- sqrt(max((dx / ratio)^2 + dy^2))
  list(cx = cx, cy = cy, rx = ratio * ry, ry = ry)
}

group_boxes <- function(d, shared = FALSE) {
  labs <- unique(d[c("term_id", "x_mm", "y_mm", "label")])
  do.call(rbind, lapply(seq_len(nrow(labs)), function(i) {
    z <- labs[i, ]
    symbols <- d[d$term_id == z$term_id, ]
    b <- list(text_bounds(z$x_mm, z$y_mm, z$label))
    for (j in seq_len(nrow(symbols))) {
      shift <- if (shared) {
        if (symbols$scope[j] == "MP") -shared_gap_mm else shared_gap_mm
      } else 0
      b[[length(b) + 1]] <- symbol_bounds(
        z$x_mm + shift, z$y_mm, symbols$shape[j], symbols$size_mm[j]
      )
    }
    do.call(rbind, b)
  }))
}

draw_text <- function(x, y, label, fontsize, face = "plain", just = "centre") {
  grid.text(label, x = unit(x, "mm"), y = grid_y(y), just = just,
            gp = gpar(fontfamily = font_family, fontsize = fontsize,
                      fontface = face, lineheight = 0.92))
}

draw_symbol <- function(x, y, shape, diameter, fill) {
  gp <- gpar(fill = fill, col = "#555555", lwd = 0.45)
  if (shape == 21) {
    grid.circle(unit(x, "mm"), grid_y(y), r = unit(diameter / 2, "mm"), gp = gp)
  } else if (shape == 22) {
    grid.rect(unit(x, "mm"), grid_y(y),
              width = unit(diameter, "mm"), height = unit(diameter, "mm"), gp = gp)
  } else {
    h <- diameter
    w <- diameter * 0.866
    grid.polygon(
      x = unit(x + c(0, -w / 2, w / 2), "mm"),
      y = grid_y(y + c(-2 * h / 3, h / 3, h / 3)),
      gp = gp
    )
  }
}

draw_ellipse <- function(e) {
  t <- seq(0, 2 * pi, length.out = 181)
  grid.lines(
    unit(e$cx + e$rx * cos(t), "mm"),
    grid_y(e$cy + e$ry * sin(t)),
    gp = gpar(col = "#AAB6BD", lty = "dashed", lwd = 0.55)
  )
}

draw_connector <- function(a, b) {
  # Draw the two common external tangents between the ellipses.
  # side = -1: upper tangent; side = 1: lower tangent.
  tangent <- function(side) {
    support <- function(e, m) sqrt((e$rx * m)^2 + e$ry^2)
    intercept <- function(e, m) e$cy - m * e$cx + side * support(e, m)
    f <- function(m) intercept(a, m) - intercept(b, m)

    slope <- uniroot(f, interval = c(-1000, 1000))$root
    sa <- support(a, slope)
    sb <- support(b, slope)

    pa <- c(
      x = a$cx - side * a$rx^2 * slope / sa,
      y = a$cy + side * a$ry^2 / sa
    )
    pb <- c(
      x = b$cx - side * b$rx^2 * slope / sb,
      y = b$cy + side * b$ry^2 / sb
    )

    grid.lines(
      unit(c(pa["x"], pb["x"]), "mm"),
      grid_y(c(pa["y"], pb["y"])),
      gp = gpar(col = "#AAB6BD", lty = "dashed", lwd = 0.55)
    )
  }

  tangent(-1)
  tangent(1)
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  panels <- c("a1", "b1", "b2", "c1", "c2", "c3", "c4")
  ellipses <- vector("list", length(panels))
  names(ellipses) <- panels

  for (sp in panels) {
    d <- nodes[nodes$subpanel == sp, ]
    shared_panel <- sp == "a1"
    microbe <- d[d$side == "Microbe", ]
    host <- d[d$side == "Host", ]
    ellipses[[sp]] <- list(
      microbe = fit_ellipse(group_boxes(microbe, shared_panel)),
      host = fit_ellipse(group_boxes(host, shared_panel))
    )
    draw_connector(ellipses[[sp]]$microbe, ellipses[[sp]]$host)
    draw_ellipse(ellipses[[sp]]$microbe)
    draw_ellipse(ellipses[[sp]]$host)
  }

  # Draw symbols first, then labels over their centres.
  for (i in seq_len(nrow(nodes))) {
    z <- nodes[i, ]
    shift <- if (z$subpanel == "a1") {
      if (z$scope == "MP") -shared_gap_mm else shared_gap_mm
    } else 0
    draw_symbol(z$x_mm + shift, z$y_mm, z$shape, z$size_mm, z$fill)
  }

  labels <- unique(nodes[c("subpanel", "side", "term_id", "x_mm", "y_mm", "label")])
  for (i in seq_len(nrow(labels))) {
    z <- labels[i, ]
    draw_text(z$x_mm, z$y_mm, z$label, label_font_pt)
  }

  # Main section headings and subpanel headings.
  panel_titles <- data.frame(
    x = c(26, 26, 26), y = c(22, 150, 286),
    label = c("(a) Shared associations",
              "(b) Microplastic-only associations",
              "(c) Plastisphere-only associations")
  )
  for (i in seq_len(nrow(panel_titles))) {
    z <- panel_titles[i, ]
    draw_text(z$x * px_mm, z$y * px_mm, z$label, panel_font_pt, "bold", "left")
  }

  subtitle_positions <- data.frame(
    subpanel = panels,
    x = c(120, 50, 434, 50, 434, 50, 434),
    y = c(39, 167, 167, 313, 313, 458, 458)
  )
  subtitle_source <- unique(rbind(
    shared[c("Figure_Subpanel", "Figure_Subpanel_Title")],
    mp[c("Figure_Subpanel", "Figure_Subpanel_Title")],
    mpp[c("Figure_Subpanel", "Figure_Subpanel_Title")]
  ))
  names(subtitle_source) <- c("subpanel", "title")
  subtitles <- subtitle_positions
  subtitles$title <- subtitle_source$title[match(
    subtitles$subpanel, subtitle_source$subpanel
  )]
  for (i in seq_len(nrow(subtitles))) {
    z <- subtitles[i, ]
    draw_text(z$x * px_mm, z$y * px_mm,
              paste0(z$subpanel, ". ", z$title),
              subtitle_font_pt, just = "left")
  }

  # Compact legend, positioned to match the reference figure.
  grid.rect(
    x = unit(581 * px_mm, "mm"), y = grid_y(84.5 * px_mm),
    width = unit(308 * px_mm, "mm"), height = unit(87 * px_mm, "mm"),
    gp = gpar(fill = "#FCEDEE", col = NA)
  )

  draw_text(452 * px_mm, 53 * px_mm, "Exposure type", legend_font_pt, "bold", "left")
  draw_text(539 * px_mm, 53 * px_mm, "Functional set type", legend_font_pt, "bold", "left")
  draw_text(685 * px_mm, 53 * px_mm, "Set size (n)", legend_font_pt, "bold", "left")

  legend_items <- data.frame(
    x = c(451, 451, 541, 541, 541),
    y = c(73, 95, 73, 95, 117),
    shape = c(22, 22, 24, 21, 22),
    fill = c("#FFB3C0", "#9878C0", rep("#B5BEC4", 3)),
    label = c("Microplastic", "Plastisphere",
              "Gut microbial KEGG map sets",
              "Host transcript sets (DEGs)",
              "Host protein sets (DEPs)")
  )
  for (i in seq_len(nrow(legend_items))) {
    z <- legend_items[i, ]
    draw_symbol(z$x * px_mm, z$y * px_mm, z$shape, 2.1, z$fill)
    draw_text((z$x + 18) * px_mm, z$y * px_mm, z$label,
              legend_font_pt, just = "left")
  }

  counts <- c(3, 18, 122)
  for (i in seq_along(counts)) {
    y <- c(73, 95, 117)[i]
    draw_symbol(697 * px_mm, y * px_mm, 21,
                size_from_n(counts[i]), "#B5BEC4")
    draw_text(716 * px_mm, y * px_mm, as.character(counts[i]),
              legend_font_pt, just = "left")
  }
}

output_dir <- file.path(figure_dir, "results", "Fig7")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cairo_pdf(file.path(output_dir, "Fig7_Association_Sets_Compact.pdf"),
          width = page_width_mm / 25.4,
          height = page_height_mm / 25.4,
          family = font_family, pointsize = 8, bg = "white")
draw_figure()
dev.off()

png(file.path(output_dir, "Fig7_Association_Sets_Compact.png"),
    width = page_width_mm, height = page_height_mm, units = "mm",
    res = 600, type = "cairo", pointsize = 8, bg = "white")
draw_figure()
dev.off()