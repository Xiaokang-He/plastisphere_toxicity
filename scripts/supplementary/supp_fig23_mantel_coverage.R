# Repository copy of the Supplementary Fig. 23 coverage plot.
# Repository-only adaptation: output path is relative to the project root.
# The reported association counts and plotting logic are retained.

# 11-map association coverage: Times New Roman, compact rows,
# black bar outlines and black numeric labels. Base R only.
# Qualifying links: Mantel r >= 0.70 and BH q_free < 0.05.
# This script plots existing results; it does not rerun Mantel tests.

# 1. Output directory and appearance settings.
repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(repo_root, "results")
out_dir <- file.path(results_dir, "supp_fig23_mantel_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
row_gap <- 0.78             # Original within-group row spacing: 1.00.
group_gap <- 0.62           # Extra spacing between row groups.
bar_half_height <- 0.14
bar_offset <- 0.155
bar_border_lwd <- 0.70
label_cex <- 8.8 / 12

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(
    "Times New Roman" = grDevices::windowsFont("Times New Roman")
  )
}
if (!isTRUE(capabilities("cairo"))) {
  stop("Cairo graphics support is required for this font export. Check capabilities('cairo').")
}

# 2. Data, in the same order as the supplied PDF.
dat <- data.frame(
  Map = c("map03040", "map03013", "map01240", "map04151", "map04710",
          "map04120", "map04137", "map00190", "map02010", "map02020", "map00515"),
  Term = c("Spliceosome", "Nucleocytoplasmic transport", "Biosynthesis of cofactors",
           "PI3K-Akt signaling pathway", "Circadian rhythm", "Ubiquitin-mediated proteolysis",
           "Mitophagy - animal", "Oxidative phosphorylation", "ABC transporters",
           "Two-component system", "Mannose-type O-glycan biosynthesis"),
  MP_RNA = c(195, 75, 31, 131, 13, 16, 57, 2, 17, 47, 9),
  MPP_RNA = c(165, 0, 5, 4, 7, 206, 342, 351, 296, 352, 46),
  MP_Protein = c(109, 33, 72, 0, 106, 4, 0, 1, 1, 20, 42),
  MPP_Protein = c(19, 0, 0, 0, 0, 33, 8, 1, 0, 1, 2),
  stringsAsFactors = FALSE
)

# Denominators are all tested host functional sets in each scope/layer.
tested <- c(MP_RNA = 511, MPP_RNA = 527, MP_Protein = 111, MPP_Protein = 35)
for (nm in names(tested)) {
  stopifnot(all(dat[[nm]] >= 0), all(dat[[nm]] <= tested[[nm]]))
  dat[[paste0(nm, "_pct")]] <- 100 * dat[[nm]] / tested[[nm]]
}

# Derive row positions, headings and separators from the spacing settings.
group_id <- c(1, rep(2, 4), rep(3, 4), rep(4, 2))
group_start <- which(!duplicated(group_id))
step <- c(0, rep(row_gap, nrow(dat) - 1))
step[group_start[-1]] <- step[group_start[-1]] + group_gap
dat$y <- cumsum(step)
stopifnot(row_gap > 2 * bar_offset + 2 * bar_half_height + 0.06)

group_headers <- data.frame(
  y = dat$y[group_start] - c(0.77, 0.65, 0.65, 0.65),
  label = c("Coverage in both groups", "Higher mainly in MP",
            "Higher mainly in MPP", "Opposite layer preferences"),
  stringsAsFactors = FALSE
)
# Optional clearer wording for the last group:
# group_headers$label[4] <- "RNA: MPP higher; Protein: MP higher"

cols <- c(MP = "#FFB3C0", MPP = "#9878C0")
ink <- "#173244"
muted <- "#627681"

y_top <- -1.15
y_bottom <- max(dat$y) + 0.73
separator_y <- group_headers$y[-1] - 0.15
stripe_half_height <- row_gap / 2 - 0.02
fig_width <- 13.2
# Shrink the physical figure too, so compressed rows are not stretched back.
fig_height <- 10.4 * (y_bottom - y_top) / (12.68 + 1.15)

# 3. Draw the figure on the current graphics device.
draw_coverage <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mar = rep(0, 4), oma = rep(0, 4), family = font_family,
      xpd = NA, xaxs = "i", yaxs = "i", lheight = 1.08, cex = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  text(0.035, 0.966,
       "C  Gut microbial-host association coverage in both groups",
       adj = c(0, 1), cex = 18.2 / 12, font = 2, col = "#102B3B")
  text(0.035, 0.922,
       "Each map retains the full coverage of MP and MPP; both shared and exposure-biased patterns remain visible.",
       adj = c(0, 0.5), cex = 10.4 / 12, col = muted)
  legend(0.64, 0.905, xjust = 0.5, yjust = 1,
         legend = c("MP - microplastic exposure", "MPP - plastisphere exposure"),
         fill = unname(cols), border = "black", bty = "n", ncol = 2,
         cex = 10.4 / 12, text.col = "black", x.intersp = 0.7)

  # Figure coordinates match the original two-panel layout.
  panels <- list(c(0.305, 0.609545, 0.165, 0.810),
                 c(0.670455, 0.975, 0.165, 0.810))
  layers <- c("RNA", "Protein")
  titles <- c("Host transcriptome", "Host proteome")
  subtitles <- c("Tested functions: MP 511 | MPP 527",
                 "Tested functions: MP 111 | MPP 35")

  for (j in seq_along(panels)) {
    xmid <- mean(panels[[j]][1:2])
    text(xmid, 0.857, titles[j], cex = 11.8 / 12, font = 2, col = ink)
    text(xmid, 0.831, subtitles[j], cex = 11.8 / 12, font = 2, col = ink)
  }
  text(0.64, 0.115, "Qualifying associations / all tested host functions",
       cex = 11.3 / 12, col = ink)
  notes <- c(
    "Labels: coverage percentage (number of qualifying associations). Two long bars indicate broad coverage in both exposures.",
    "Qualifying: Mantel r >= 0.70 and BH q_free < 0.05 (unrestricted permutations). Row groups describe coverage patterns.",
    "Scope-specific DEG/DEP sets; no formal between-scope comparison. Within-treatment FDR support was not obtained."
  )
  text(0.035, c(0.075, 0.051, 0.027), notes,
       adj = c(0, 0.5), cex = 8.7 / 12, col = muted)

  for (j in seq_along(layers)) {
    par(fig = panels[[j]], mar = rep(0, 4), new = TRUE, family = font_family,
        xpd = NA, mgp = c(2, 0.7, 0))
    plot.new()
    plot.window(xlim = c(0, 102), ylim = c(y_bottom, y_top),
                xaxs = "i", yaxs = "i")

    # Alternating row backgrounds.
    for (i in seq_len(nrow(dat))) {
      if (i %% 2 == 1) {
        rect(0, dat$y[i] - stripe_half_height, 102, dat$y[i] + stripe_half_height,
             col = "#F3F6F8", border = NA)
      }
    }
    abline(v = seq(0, 100, 25), col = "#DAE3E8", lwd = 0.65)
    abline(v = 0, col = "#A1B3BE", lwd = 0.8)
    abline(h = separator_y, col = "#BCCCD5", lwd = 0.7)

    for (scope in c("MP", "MPP")) {
      nm <- paste(scope, layers[j], sep = "_")
      pct <- dat[[paste0(nm, "_pct")]]
      n <- dat[[nm]]
      yy <- dat$y + if (scope == "MP") -bar_offset else bar_offset
      rect(0, yy - bar_half_height, pct, yy + bar_half_height,
           col = cols[[scope]], border = "black", lwd = bar_border_lwd)
      labels <- sprintf("%.1f%% (%d)", pct, as.integer(n))
      # Fit labels to the actual Times New Roman text width.
      label_pad <- 1.4
      text_width <- strwidth(labels, cex = label_cex, font = 2)
      inside <- pct > 0 & pct >= text_width + 2 * label_pad
      if (any(inside)) {
        text(pct[inside] - label_pad, yy[inside], labels[inside],
             adj = c(1, 0.5), col = "black", font = 2, cex = label_cex)
      }
      if (any(!inside)) {
        text(pct[!inside] + label_pad, yy[!inside], labels[!inside],
             adj = c(0, 0.5), col = "black", cex = label_cex)
      }
    }

    axis(1, pos = y_bottom, at = seq(0, 100, 25),
         labels = paste0(seq(0, 100, 25), "%"),
         tcl = -0.25, cex.axis = 10 / 12, las = 1,
         col = "#AABBC5", col.axis = muted, lwd = 0.8)

    if (j == 1) {
      text(-4, dat$y, paste(dat$Term, dat$Map, sep = "\n"),
           adj = c(1, 0.5), cex = 9.8 / 12, col = ink)
      text(-2.55, group_headers$y, group_headers$label,
           adj = c(1, 0.5), cex = 9.5 / 12, font = 2, col = "#3D5665")
    }
  }
  invisible(NULL)
}

# 4. Cairo devices use the requested system font in both outputs.
grDevices::cairo_pdf(filename = file.path(out_dir, "map_coverage_absolute_comparison_R.pdf"),
                     width = fig_width, height = fig_height, pointsize = 12,
                     family = font_family, bg = "white")
tryCatch(draw_coverage(), finally = grDevices::dev.off())

png_args <- list(filename = file.path(out_dir, "map_coverage_absolute_comparison_R.png"),
                 width = fig_width, height = fig_height, units = "in", res = 600,
                 pointsize = 12, bg = "white", type = "cairo", family = font_family)
do.call(grDevices::png, png_args)
tryCatch(draw_coverage(), finally = grDevices::dev.off())

message("Saved figures to: ", normalizePath(out_dir, winslash = "/"))