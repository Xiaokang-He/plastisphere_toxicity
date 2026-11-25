# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# Core 105 VFDB VFG cumulative abundance
# Figure style follows the siphoning-efficiency reference plot:
#   - bars = mean +/- SD
#   - overlay the 3 individual sample totals per group
#   - show exact pairwise WRS2 lincon P values on brackets (not stars)
#   - retain the original five-group colors and hatch directions
#
# DATA:
#   Each sample total = sum of the 105 ORIGINAL core VFG abundances.
#   The 105 VFG entries are NOT treated as biological replicates.
#
# STATISTICS:
#   Overall: WRS2 t1way() robust one-way test
#   Pairwise: WRS2 lincon() for all 10 comparisons
#   Trim level: tr = 0.20
#   P-value adjustment: method = "none"
#   No overall t1way significance gate
#   No multiple-comparison P-value adjustment
#   Only pairwise comparisons with raw P < 0.05 are annotated.
#
# IMPORTANT:
#   The original sample mapping and abundance calculation are retained.
#   The overall t1way result is reported separately and is NOT a
#   prerequisite for displaying pairwise lincon comparisons.
#   lincon confidence intervals remain simultaneous intervals; the
#   method = "none" setting applies to P-value adjustment.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ============================================================
# 1. Settings
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- input_dir
input_file <- file.path(input_dir, "vfdb_core_VFDB_ID_profile.xls")
out_dir <- file.path(results_dir, "supp_fig19_vfg_cumulative_abundance")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

groups <- c("CTR", "MP1", "MP10", "MPP1", "MPP10")

prefix <- c(
  CTR = "T11",
  MP1 = "T6",
  MP10 = "T10",
  MPP1 = "T4",
  MPP10 = "T3"
)

group_labels <- c(
  CTR = "Control",
  MP1 = "Microplastic-Low",
  MP10 = "Microplastic-High",
  MPP1 = "Plastisphere-Low",
  MPP10 = "Plastisphere-High"
)

# Original five-group colors
colors <- c(
  CTR  = "#b6dcd5",
  MP1  = "#ffe0e2",
  MP10 = "#ffaab9",
  MPP1 = "#baa5d5",
  MPP10 = "#8c69b9"
)

# Original hatch directions
# Control: horizontal; MP: diagonal /; MPP: diagonal \
angles <- c(
  CTR = 0,
  MP1 = 45,
  MP10 = 45,
  MPP1 = 135,
  MPP10 = 135
)

# Plot settings
show_points <- TRUE
show_only_significant <- TRUE

figure_width <- 6.4
figure_height <- 6.5
figure_dpi <- 600L

# Raw sample points: horizontally separated within each bar
point_offsets <- c(-0.13, 0, 0.13)

# ============================================================
# 2. Packages
# ============================================================

packages <- c(
  "ggplot2",
  "ggpattern",
  "WRS2"
)

missing <- packages[
  !vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpattern)
  library(WRS2)
})

# WRS2 >= 1.1-3 is required for lincon(method = "none").
if (utils::packageVersion("WRS2") < "1.1.3") {
  stop(
    "Please update WRS2 to version 1.1-3 or later. ",
    "The method argument in lincon() is required for unadjusted P values."
  )
}

# ============================================================
# 3. Read the tab-delimited .xls file
# ============================================================

if (!file.exists(input_file)) {
  stop(
    "Input file not found:\n",
    input_file
  )
}

dat <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

names(dat) <- trimws(
  gsub("\ufeff", "", names(dat), fixed = TRUE)
)

if (anyDuplicated(names(dat))) {
  stop("Duplicated column names.")
}

if (!"VFDB_ID" %in% names(dat)) {
  stop("Missing VFDB_ID column.")
}

# Retain the VFDB core subset if a Database column is present
if ("Database" %in% names(dat)) {
  dat <- dat[
    tolower(trimws(dat$Database)) == "core",
    ,
    drop = FALSE
  ]
}

ids <- trimws(as.character(dat$VFDB_ID))

if (
  nrow(dat) != 105L ||
  anyNA(ids) ||
  any(!nzchar(ids)) ||
  anyDuplicated(ids)
) {
  stop(
    "Expected exactly 105 nonmissing, unique core VFG entries."
  )
}

# ============================================================
# 4. Sample columns and group information
# ============================================================

sample_cols <- unlist(
  lapply(
    groups,
    function(g) {
      paste0(prefix[g], "_M_", 1:3)
    }
  ),
  use.names = FALSE
)

if (!all(sample_cols %in% names(dat))) {
  stop(
    "Missing sample columns: ",
    paste(
      setdiff(sample_cols, names(dat)),
      collapse = ", "
    )
  )
}

X <- as.matrix(
  data.frame(
    lapply(
      dat[, sample_cols, drop = FALSE],
      function(v) {
        suppressWarnings(
          as.numeric(as.character(v))
        )
      }
    ),
    check.names = FALSE
  )
)

colnames(X) <- sample_cols
rownames(X) <- ids

if (any(!is.finite(X)) || any(X < 0)) {
  stop("Invalid, missing or negative abundance.")
}

# Verify the supplied Total column, if available
if ("Total" %in% names(dat)) {
  supplied_total <- suppressWarnings(
    as.numeric(as.character(dat$Total))
  )
  
  if (
    !isTRUE(
      all.equal(
        unname(rowSums(X)),
        supplied_total,
        tolerance = 1e-8
      )
    )
  ) {
    stop(
      "The Total column does not match the sum of the 15 sample columns."
    )
  }
}

# ============================================================
# 5. Calculate cumulative abundance for each biological sample
# ============================================================

# One total per biological sample: sum of 105 original VFG entries.
sample_total <- data.frame(
  Sample = sample_cols,
  Group = factor(
    rep(groups, each = 3),
    levels = groups
  ),
  Replicate = rep(1:3, times = 5),
  Cumulative_VFG_abundance = unname(colSums(X))
)

# Group mean, SD and median; n = 3 biological samples per group
group_summary <- do.call(
  rbind,
  lapply(
    groups,
    function(g) {
      v <- sample_total$Cumulative_VFG_abundance[
        sample_total$Group == g
      ]
      
      data.frame(
        Group = g,
        n = length(v),
        Mean = mean(v),
        SD = sd(v),
        Median = median(v)
      )
    }
  )
)

group_summary$Group <- factor(
  group_summary$Group,
  levels = groups
)

group_summary$Figure_label <- unname(
  group_labels[as.character(group_summary$Group)]
)

group_summary$Plot_x <- seq_along(groups)

# Slight horizontal separation of the three raw sample points
sample_total$Plot_x <- rep(
  seq_along(groups),
  each = 3
) + rep(
  point_offsets,
  times = length(groups)
)

# ============================================================
# 6. WRS2 t1way() + lincon()
#    All 10 pairwise comparisons; no P-value adjustment
# ============================================================

TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05

# The analysis unit is one cumulative VFG abundance total
# per independent biological sample (n = 3 per group).
#
# t1way() is reported as the overall robust one-way test.
# Its significance is NOT used as a gate for pairwise tests.
#
# lincon(method = "none") calculates all pairwise comparisons
# without multiple-comparison adjustment of P values.

run_wrs2 <- function() {

  d <- sample_total
  d$Group <- factor(
    as.character(d$Group),
    levels = groups
  )

  # Confirm that all five groups are represented.
  group_n <- table(
    factor(d$Group, levels = groups)
  )

  if (any(group_n < 2)) {
    stop(
      "At least one group has fewer than two biological samples. ",
      "Please check the sample mapping."
    )
  }

  # ----------------------------------------------------------
  # 6.1 Overall robust one-way test: WRS2 t1way()
  # ----------------------------------------------------------

  overall_fit <- WRS2::t1way(
    Cumulative_VFG_abundance ~ Group,
    data = d,
    tr = TRIM_LEVEL
  )

  overall_result <- data.frame(
    Test = "WRS2_t1way",
    Trim_level = TRIM_LEVEL,
    Statistic_F = as.numeric(overall_fit$test),
    df1 = as.numeric(overall_fit$df1),
    df2 = as.numeric(overall_fit$df2),
    P_t1way = as.numeric(overall_fit$p.value),
    Effect_size = if (!is.null(overall_fit$effsize)) {
      as.numeric(overall_fit$effsize)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # 6.2 All pairwise robust comparisons: WRS2 lincon()
  #
  # method = "none": no P-value adjustment.
  # All 10 pairwise comparisons are retained.
  # No overall t1way significance gate is applied.
  # ----------------------------------------------------------

  lincon_fit <- WRS2::lincon(
    Cumulative_VFG_abundance ~ Group,
    data = d,
    tr = TRIM_LEVEL,
    alpha = ALPHA_LEVEL,
    method = "none"
  )

  comp <- as.matrix(lincon_fit$comp)

  if (ncol(comp) != 6) {
    stop(
      "Unexpected lincon result structure. ",
      "Please check your installed WRS2 version."
    )
  }

  # WRS2 lincon output:
  # Column 1: first group index
  # Column 2: second group index
  # Column 3: trimmed mean difference (psihat)
  # Column 4: simultaneous CI lower limit
  # Column 5: simultaneous CI upper limit
  # Column 6: P value under the specified method
  #
  # With method = "none", column 6 contains unadjusted P values.
  # The confidence intervals remain simultaneous confidence intervals.

  group_names <- as.character(lincon_fit$fnames)

  if (
    length(group_names) != length(groups) ||
    !setequal(group_names, groups)
  ) {
    stop(
      "The group names returned by WRS2 lincon() do not match ",
      "the five expected groups."
    )
  }

  group1_index <- as.integer(comp[, 1])
  group2_index <- as.integer(comp[, 2])

  if (
    anyNA(group1_index) ||
    anyNA(group2_index) ||
    any(group1_index < 1 | group1_index > length(group_names)) ||
    any(group2_index < 1 | group2_index > length(group_names))
  ) {
    stop("Invalid group indices in WRS2 lincon() results.")
  }

  lincon_result <- data.frame(
    Group1 = group_names[group1_index],
    Group2 = group_names[group2_index],
    psihat = as.numeric(comp[, 3]),
    CI_lower_simultaneous = as.numeric(comp[, 4]),
    CI_upper_simultaneous = as.numeric(comp[, 5]),
    P_raw = as.numeric(comp[, 6]),
    P_adjustment = "none",
    Trim_level = TRIM_LEVEL,
    P_t1way = as.numeric(overall_fit$p.value),
    stringsAsFactors = FALSE
  )

  # Five groups must produce exactly 10 pairwise comparisons.
  if (nrow(lincon_result) != choose(length(groups), 2)) {
    stop(
      "WRS2 lincon() did not return all 10 pairwise comparisons."
    )
  }

  if (
    anyNA(lincon_result$P_raw) ||
    any(
      lincon_result$P_raw < 0 |
      lincon_result$P_raw > 1
    )
  ) {
    stop("Invalid P values in WRS2 lincon() results.")
  }

  list(
    overall = overall_result,
    pairwise = lincon_result
  )
}

# Run WRS2 analysis.
wrs2_results <- run_wrs2()

t1way_results <- wrs2_results$overall
lincon_results <- wrs2_results$pairwise

# IMPORTANT:
# Do not filter pairwise comparisons by the overall t1way P value.
# Do not adjust pairwise P values.
# Only raw lincon P < 0.05 determines figure annotations.

lincon_significant <- lincon_results[
  lincon_results$P_raw < ALPHA_LEVEL,
  ,
  drop = FALSE
]

cat("\n========== WRS2 t1way overall result ==========\n")
print(
  t1way_results,
  row.names = FALSE,
  digits = 6
)

cat("\n========== WRS2 lincon: all 10 pairwise comparisons ==========\n")
print(
  lincon_results,
  row.names = FALSE,
  digits = 6
)

cat("\n========== Pairwise comparisons with raw P < 0.05 ==========\n")
if (nrow(lincon_significant) > 0) {
  print(
    lincon_significant,
    row.names = FALSE,
    digits = 6
  )
} else {
  cat("No pairwise comparison has raw P < 0.05.\n")
}


# ============================================================
# 7. Select comparisons to annotate and format exact P labels
# ============================================================

if (show_only_significant) {
  brackets <- lincon_significant
} else {
  brackets <- lincon_results
}

if (nrow(brackets) > 0) {

  # Map the returned group names to the plotting x coordinates.
  brackets$X1 <- match(
    as.character(brackets$Group1),
    groups
  )

  brackets$X2 <- match(
    as.character(brackets$Group2),
    groups
  )

  if (anyNA(brackets$X1) || anyNA(brackets$X2)) {
    stop("Could not map WRS2 lincon group names to plotting positions.")
  }

  # Ensure the left and right bracket coordinates are in ascending order.
  x_left <- pmin(brackets$X1, brackets$X2)
  x_right <- pmax(brackets$X1, brackets$X2)
  brackets$X1 <- x_left
  brackets$X2 <- x_right

  brackets$Label1 <- unname(
    group_labels[groups[brackets$X1]]
  )

  brackets$Label2 <- unname(
    group_labels[groups[brackets$X2]]
  )

  # Shorter brackets are placed lower; longer brackets higher.
  brackets <- brackets[
    order(
      brackets$X2 - brackets$X1,
      brackets$X1,
      brackets$X2
    ),
    ,
    drop = FALSE
  ]
}

# Display raw P values directly, rounded to three decimal places.
# Plotmath is used so that P appears in italics.
format_p_plotmath <- function(p) {

  if (is.na(p)) {
    return('"P = NA"')
  }

  if (p < 0.001) {
    return("italic(P) < 0.001")
  }

  paste0(
    "italic(P) == ",
    sprintf("%.3f", p)
  )
}

format_p_text <- function(p) {

  if (is.na(p)) {
    return("P = NA")
  }

  if (p < 0.001) {
    return("P < 0.001")
  }

  sprintf("P = %.3f", p)
}


# ============================================================
# 8. Calculate bracket heights and y-axis upper limit
# ============================================================

data_top <- max(
  group_summary$Mean + group_summary$SD,
  sample_total$Cumulative_VFG_abundance,
  na.rm = TRUE
)

plot_span <- max(data_top, 1)

if (nrow(brackets) > 0) {
  
  # Vertical spacing between stacked comparison brackets
  bracket_step <- plot_span * 0.075
  
  brackets$Y <- (
    data_top +
      plot_span * 0.055 +
      (seq_len(nrow(brackets)) - 1) * bracket_step
  )
  
  # Text position above each bracket
  brackets$Label_y <- brackets$Y +
    plot_span * 0.025
  
  brackets$P_label <- vapply(
    brackets$P_raw,
    format_p_plotmath,
    character(1)
  )
  
  # Leave enough space above the highest P-value label
  y_required <- max(
    brackets$Label_y + plot_span * 0.065,
    data_top * 1.10
  )
  
} else {
  
  brackets$Y <- numeric(0)
  brackets$Label_y <- numeric(0)
  brackets$P_label <- character(0)
  
  y_required <- data_top * 1.12
}

# Round the upper limit upward to a convenient 50-unit interval
y_top <- ceiling(
  y_required / 50
) * 50

if (y_top <= data_top) {
  y_top <- ceiling(
    data_top * 1.15 / 50
  ) * 50
}

tip_height <- plot_span * 0.018

# ============================================================
# 9. Statistical summary and caption
# ============================================================

t1way_label <- paste0(
  "WRS2 t1way: F = ",
  sprintf("%.2f", t1way_results$Statistic_F[1]),
  ", df1 = ",
  sprintf("%.0f", t1way_results$df1[1]),
  ", df2 = ",
  sprintf("%.2f", t1way_results$df2[1]),
  ", ",
  format_p_text(t1way_results$P_t1way[1])
)

caption_text <- paste0(
  "Bars represent mean +/- SD (n = 3); points represent individual sample totals. ",
  t1way_label,
  ".\n",
  "Pairwise comparisons: WRS2 lincon (tr = 0.20; method = none), ",
  "with no P-value adjustment and no overall-test significance gate. ",
  "Only comparisons with raw P < 0.05 are annotated."
)

cat("\n============================================\n")
cat("WRS2 t1way overall result\n")
cat("============================================\n")
print(
  t1way_results,
  row.names = FALSE,
  digits = 6
)

cat("\nAll 10 WRS2 lincon pairwise comparisons (unadjusted raw P):\n")
print(
  lincon_results[
    ,
    c(
      "Group1",
      "Group2",
      "psihat",
      "CI_lower_simultaneous",
      "CI_upper_simultaneous",
      "P_raw",
      "P_adjustment"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat("\nComparisons displayed on the figure:\n")
if (nrow(brackets) > 0) {
  print(
    brackets[
      ,
      c(
        "Label1",
        "Label2",
        "P_raw"
      )
    ],
    row.names = FALSE,
    digits = 6
  )
} else {
  cat("No pairwise comparison has raw P < 0.05.\n")
}


# ============================================================
# 10. Font settings
# ============================================================

font_export <- "serif"

if (
  .Platform$OS.type == "windows" &&
  isTRUE(capabilities("cairo"))
) {
  grDevices::windowsFonts(
    "Times New Roman" = grDevices::windowsFont("Times New Roman")
  )
  
  font_export <- "Times New Roman"
}

# ============================================================
# 11. Plot function
# ============================================================

make_plot <- function(family = "serif") {
  
  p <- ggplot() +
    
    # --------------------------------------------------------
  # Bars: mean abundance
  # --------------------------------------------------------
  ggpattern::geom_col_pattern(
    data = group_summary,
    
    aes(
      x = Plot_x,
      y = Mean,
      fill = Group,
      pattern_fill = Group,
      pattern_angle = Group
    ),
    
    width = 0.82,
    colour = "black",
    linewidth = 0.65,
    
    pattern = "stripe",
    pattern_colour = "black",
    pattern_size = 0.22,
    pattern_density = 0.010,
    pattern_spacing = 0.085,
    
    show.legend = FALSE
  ) +
    
    # --------------------------------------------------------
  # Individual sample points: three per group
  # White-filled circles with black outlines, as in Fig. 2
  # --------------------------------------------------------
  {
    if (show_points) {
      geom_point(
        data = sample_total,
        
        aes(
          x = Plot_x,
          y = Cumulative_VFG_abundance
        ),
        
        shape = 21,
        size = 2.7,
        stroke = 0.65,
        fill = "white",
        colour = "black",
        alpha = 1,
        
        show.legend = FALSE
      )
    }
  } +
    
    # --------------------------------------------------------
  # Mean +/- SD error bars
  # --------------------------------------------------------
  geom_errorbar(
    data = group_summary,
    
    aes(
      x = Plot_x,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    
    width = 0.18,
    colour = "black",
    linewidth = 0.65
  ) +
    
    # --------------------------------------------------------
  # Manual color and hatch settings
  # --------------------------------------------------------
  scale_fill_manual(
    values = colors,
    drop = FALSE,
    guide = "none"
  ) +
    
    scale_pattern_fill_manual(
      values = colors,
      drop = FALSE,
      guide = "none"
    ) +
    
    scale_pattern_angle_manual(
      values = angles,
      drop = FALSE,
      guide = "none"
    ) +
    
    # --------------------------------------------------------
  # X axis
  # --------------------------------------------------------
  scale_x_continuous(
    breaks = seq_along(groups),
    labels = unname(group_labels[groups]),
    limits = c(0.4, 5.6),
    expand = expansion(mult = 0)
  ) +
    
    # --------------------------------------------------------
  # Y axis
  # --------------------------------------------------------
  scale_y_continuous(
    limits = c(0, y_top),
    
    breaks = seq(
      0,
      y_top,
      by = 100
    ),
    
    minor_breaks = seq(
      0,
      y_top,
      by = 50
    ),
    
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
    
    labs(
      title = "Cumulative core VFG abundance",
      
      x = NULL,
      
      y = "Cumulative core VFG abundance\n(input profile units)",
      
      caption = caption_text
    ) +
    
    theme_bw(
      base_family = family,
      base_size = 13
    ) +
    
    theme(
      
      # Title: centered and bold, similar to the reference figure
      plot.title = element_text(
        family = family,
        size = 17,
        face = "bold",
        hjust = 0.5,
        colour = "black",
        margin = margin(b = 10)
      ),
      
      # No legend: groups are identified by the x-axis
      legend.position = "none",
      
      # Panel border
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.8
      ),
      
      # Grid lines
      panel.grid.major = element_line(
        colour = "#E6E6E6",
        linewidth = 0.4
      ),
      
      panel.grid.minor.x = element_blank(),
      
      panel.grid.minor.y = element_line(
        colour = "#F0F0F0",
        linewidth = 0.3
      ),
      
      # Axis text
      axis.text = element_text(
        colour = "black",
        family = family
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 30,
        hjust = 1,
        vjust = 1
      ),
      
      axis.text.y = element_text(
        size = 12
      ),
      
      # Axis titles
      axis.title.y = element_text(
        size = 15,
        face = "bold",
        colour = "black",
        margin = margin(r = 9)
      ),
      
      # Axis ticks
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks.length = grid::unit(
        0.12,
        "cm"
      ),
      
      # Caption: statistical method and overall WRS2 t1way result
      plot.caption = element_text(
        family = family,
        size = 8.5,
        hjust = 0,
        lineheight = 1.15,
        colour = "black",
        margin = margin(t = 9)
      ),
      
      plot.margin = margin(
        8,
        12,
        8,
        8
      )
    )
  
  # ----------------------------------------------------------
  # Pairwise WRS2 lincon brackets with exact raw P values (not asterisks)
  # ----------------------------------------------------------
  if (nrow(brackets) > 0) {
    
    p <- p +
      
      # Horizontal bracket lines
      geom_segment(
        data = brackets,
        
        aes(
          x = X1,
          xend = X2,
          y = Y,
          yend = Y
        ),
        
        inherit.aes = FALSE,
        linewidth = 0.65,
        colour = "black"
      ) +
      
      # Left vertical bracket ends
      geom_segment(
        data = brackets,
        
        aes(
          x = X1,
          xend = X1,
          y = Y,
          yend = Y - tip_height
        ),
        
        inherit.aes = FALSE,
        linewidth = 0.65,
        colour = "black"
      ) +
      
      # Right vertical bracket ends
      geom_segment(
        data = brackets,
        
        aes(
          x = X2,
          xend = X2,
          y = Y,
          yend = Y - tip_height
        ),
        
        inherit.aes = FALSE,
        linewidth = 0.65,
        colour = "black"
      ) +
      
      # Actual P values; italic P, three decimal places
      geom_text(
        data = brackets,
        
        aes(
          x = (X1 + X2) / 2,
          y = Label_y,
          label = P_label
        ),
        
        inherit.aes = FALSE,
        parse = TRUE,
        family = family,
        size = 4.0,
        colour = "black",
        vjust = 0
      )
  }
  
  p
}

# ============================================================
# 12. Generate and preview plot
# ============================================================

p_final <- make_plot(font_export)

print(p_final)

# ============================================================
# 13. Export PDF and PNG
# ============================================================

stem <- "Core105_cumulative_abundance_meanSD_points_exactP"

if (isTRUE(capabilities("cairo"))) {
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(stem, ".pdf")
    ),
    
    plot = p_final,
    
    device = grDevices::cairo_pdf,
    family = font_export,
    
    width = figure_width,
    height = figure_height,
    units = "in",
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(stem, ".png")
    ),
    
    plot = p_final,
    
    width = figure_width,
    height = figure_height,
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    type = "cairo"
  )
  
} else {
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(stem, ".pdf")
    ),
    
    plot = p_final,
    
    device = "pdf",
    family = "Times",
    
    width = figure_width,
    height = figure_height,
    units = "in",
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(stem, ".png")
    ),
    
    plot = p_final,
    
    width = figure_width,
    height = figure_height,
    units = "in",
    dpi = figure_dpi,
    bg = "white"
  )
}

# ============================================================
# 14. Export all numerical results
# ============================================================

# 14.1 Individual sample totals used in WRS2 tests.
write.csv(
  sample_total[
    ,
    c(
      "Sample",
      "Group",
      "Replicate",
      "Cumulative_VFG_abundance"
    )
  ],
  file.path(
    out_dir,
    "Core105_sample_totals.csv"
  ),
  row.names = FALSE
)

# 14.2 Descriptive statistics for each group.
write.csv(
  group_summary[
    ,
    c(
      "Group",
      "Figure_label",
      "n",
      "Mean",
      "SD",
      "Median"
    )
  ],
  file.path(
    out_dir,
    "Core105_group_summary.csv"
  ),
  row.names = FALSE
)

# 14.3 Overall WRS2 t1way result.
write.csv(
  t1way_results,
  file.path(
    out_dir,
    "WRS2_t1way_overall_results.csv"
  ),
  row.names = FALSE
)

# 14.4 All 10 WRS2 lincon pairwise comparisons (unadjusted raw P).
write.csv(
  lincon_results,
  file.path(
    out_dir,
    "WRS2_lincon_all10_pairwise_rawP_no_adjustment.csv"
  ),
  row.names = FALSE
)

# 14.5 Pairwise comparisons displayed on the figure.
write.csv(
  brackets,
  file.path(
    out_dir,
    "WRS2_lincon_displayed_pairwise_rawP.csv"
  ),
  row.names = FALSE
)

# 14.6 Group sample sizes used for the robust analysis.
WRS2_group_n <- data.frame(
  Group = factor(groups, levels = groups),
  n = as.integer(table(factor(sample_total$Group, levels = groups)))
)

write.csv(
  WRS2_group_n,
  file.path(
    out_dir,
    "WRS2_group_sample_sizes.csv"
  ),
  row.names = FALSE
)

# 14.7 Input information.
write.csv(
  data.frame(
    Input_file = input_file,
    N_core_entries = nrow(X),
    Number_of_biological_samples = ncol(X),
    Sample_columns = paste(
      sample_cols,
      collapse = ";"
    )
  ),
  file.path(
    out_dir,
    "Input_used.csv"
  ),
  row.names = FALSE
)

# Analysis notes.
analysis_notes <- c(
  "Sample total = sum of the 105 ORIGINAL VFDB core VFG abundances.",
  "No log transformation or row-Z transformation was applied before summation.",
  "The 105 VFG entries are not biological replicates.",
  "Each biological sample contributes one cumulative VFG abundance value.",
  "Bars represent arithmetic mean +/- SD across three biological sample totals per group.",
  "Individual sample totals are shown as white-filled circles with black outlines.",
  "Overall test: WRS2 t1way().",
  "Pairwise test: WRS2 lincon().",
  "Trimmed mean proportion: tr = 0.20.",
  "All five experimental groups were included in both WRS2 analyses.",
  "All 10 pairwise comparisons were calculated.",
  "P-value adjustment: method = none; no multiple-comparison P-value adjustment was applied.",
  "No overall t1way significance gate was applied to pairwise comparisons.",
  "Only pairwise comparisons with raw lincon P < 0.05 are annotated by default.",
  "The overall WRS2 t1way result is reported separately from pairwise lincon P values.",
  "P values on brackets are shown directly to three decimal places; P < 0.001 is shown as P < 0.001.",
  "WRS2 lincon confidence intervals remain simultaneous confidence intervals; method = none applies to P-value adjustment."
)

writeLines(
  analysis_notes,
  con = file.path(
    out_dir,
    "WRS2_analysis_notes.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "R_session_info.txt"
  )
)


# ============================================================
# 15. Final console summary
# ============================================================

cat("\n============================================\n")
cat("Core 105 VFG cumulative abundance plot completed.\n")
cat("Bars: mean +/- SD; points: individual sample totals.\n")
cat("Overall test: WRS2 t1way().\n")
cat("Pairwise test: WRS2 lincon(), method = none.\n")
cat("No overall-test gate; no P-value adjustment.\n")

cat("\nWRS2 t1way overall result:\n")
print(
  t1way_results,
  row.names = FALSE,
  digits = 6
)

cat("\nAll 10 WRS2 lincon pairwise comparisons:\n")
print(
  lincon_results[
    ,
    c(
      "Group1",
      "Group2",
      "P_raw",
      "P_adjustment"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat("\nDisplayed pairwise raw P values:\n")
if (nrow(brackets) > 0) {
  for (i in seq_len(nrow(brackets))) {
    cat(
      brackets$Label1[i],
      " vs ",
      brackets$Label2[i],
      ": ",
      format_p_text(brackets$P_raw[i]),
      "\n",
      sep = ""
    )
  }
} else {
  cat("None (all raw pairwise P >= 0.05).\n")
}

cat("\nOutput directory:\n")
cat(
  normalizePath(
    out_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)
cat("============================================\n")