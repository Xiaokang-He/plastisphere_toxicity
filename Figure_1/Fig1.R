# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# Figure 1: physiological and biochemical endpoints; compact WRS2 annotations with the original Y-axis limits preserved
#
# Page: 950 × 400 mm
# Layout: 5 columns × 2 rows; each cell measures 190 × 200 mm
#
# Exact dimensions for each panel:
#   XY plotting area: 89.277 × 136.515 mm
#   Title box: 89.277 × 18.238 mm
#
# Typography (Times New Roman):
#   Title: 36.59 pt, bold
#   Y-axis tick labels: 28.4 pt
#   Y-axis title: 28.4 pt, bold
#   X-axis tick labels: 28.4 pt, rotated 30°
#   Significance P values: 18 pt by default; shrink adaptively to a minimum of 10 pt when space is limited; italicize only the P character
#
# Line widths:
#   Bar borders, error bars, XY panel border, major grid lines, and significance brackets: 1.788 pt
#   Minor grid lines: 0.886 pt
#
# Bar plots:
#   Mean ± SD
#   Show one white-filled, black-outlined point for each independent biological replicate
#   Apply slight horizontal jitter only; preserve the original y values
#
# Statistics:
#   WRS2 t1way() + lincon()
#   t1way: heteroscedasticity-robust one-way omnibus test across five groups (tr = 0.2)
#   lincon: all 10 pairwise comparisons among the five groups; no multiple-comparison adjustment of P values
#   Display the exact unadjusted P value when P < 0.05; do not use asterisks
#   Note: with n = 3 per group, tr = 0.2 does not trim any observations.
#
# Font output:
#   Illustrator-editable version: svglite exports native SVG <text> elements
#   Auxiliary PDF version: Cairo PDF preserves text objects
#   Disable automatic showtext conversion to prevent text from becoming vector paths
#   Do not export any font-outlined files
#
# Y-axis limits: prioritize the ymax values specified in the original plot_config.
# Fit significance annotations into the existing upper Y-axis space whenever possible; reduce the P-value font size adaptively if needed,
# and expand the Y-axis upper limit only by the minimum amount when the minimum font size is still insufficient.
# Do not substantially stretch the Y-axis merely to accommodate significance brackets.
#
# This version does not use ggplot's default axis-text layout or cowplot scaling.
# Extract each ggplot panel directly and position it manually on a 950 × 400 mm page:
#   title box, Y-axis title, Y-axis tick labels, X-axis tick labels, and XY plotting area.
# ============================================================


# ============================================================
# 0. Clear the workspace and load dependencies
# ============================================================

rm(list = ls())

required_packages <- c(
  "dplyr",
  "ggplot2",
  "ggpattern",
  "WRS2",
  "svglite"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpattern)
  library(WRS2)
  library(grid)
})


# ============================================================
# 1. Input data and output directories
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
figure_dir <- file.path(repo_root, "Figure_1")
input_dir <- figure_dir
results_dir <- file.path(figure_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- figure_dir

input_file <- file.path(
  figure_dir,
  "Clam_Physiological_and_Biochemical_Endpoints.csv"
)

if (!file.exists(input_file)) {
  stop(
    "Raw data file not found:\n",
    input_file,
    "\nPlease check the file name and path."
  )
}

# On some Windows systems with Chinese locale settings, the SVG device may not reliably open output paths containing non-ASCII characters.
# Export results to an ASCII-only path; continue reading the source CSV from base_dir.
output_dir <- file.path(results_dir, "Fig1")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop("Cannot create output directory:\n", output_dir)
}

# Check whether the output directory is writable to avoid errors when opening the SVG/PDF device.
test_file <- file.path(output_dir, "__write_test.tmp")
can_write <- tryCatch({
  writeLines("test", test_file, useBytes = TRUE)
  unlink(test_file)
  TRUE
}, error = function(e) FALSE)

if (!can_write) {
  stop("Output directory is not writable. Check path permissions:\n", output_dir)
}


# ============================================================
# 2. Register Times New Roman fonts
#
# Use svglite to export native <text> elements for Illustrator-editable SVG output.
# Do not use showtext or sysfonts; keep text as editable text rather than vector outlines.
#
# Important: if showtext_auto(TRUE) was previously enabled in the same R session,
# it may continue to affect subsequent devices. Explicitly disable automatic showtext conversion here.
# ============================================================

font_family <- "Times New Roman"

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(
    "Times New Roman" = grDevices::windowsFont("Times New Roman")
  )
}

# Prevent a showtext automatic mode left over from a previous session from converting text to paths.
if (requireNamespace("showtext", quietly = TRUE)) {
  try(
    showtext::showtext_auto(FALSE),
    silent = TRUE
  )
}

if (!isTRUE(capabilities("cairo"))) {
  stop(
    "The current R environment does not support the Cairo PDF device.",
    "\nPlease use a Windows build of R with Cairo support."
  )
}

# Check that the Windows font files exist.
font_dir <- file.path(Sys.getenv("WINDIR"), "Fonts")
font_files <- c(
  regular = file.path(font_dir, "times.ttf"),
  bold = file.path(font_dir, "timesbd.ttf"),
  italic = file.path(font_dir, "timesi.ttf"),
  bolditalic = file.path(font_dir, "timesbi.ttf")
)

if (!all(file.exists(font_files))) {
  stop(
    "The complete Times New Roman font files were not found:\n",
    paste(font_files[!file.exists(font_files)], collapse = "\n")
  )
}

# ============================================================
# 3. Page and physical figure dimensions
# ============================================================

PAGE_WIDTH_MM  <- 950
PAGE_HEIGHT_MM <- 400

N_COLS <- 5
N_ROWS <- 2

CELL_WIDTH_MM  <- PAGE_WIDTH_MM / N_COLS
CELL_HEIGHT_MM <- PAGE_HEIGHT_MM / N_ROWS

# User-specified XY plotting-area dimensions
PANEL_WIDTH_MM  <- 89.277
PANEL_HEIGHT_MM <- 136.515

# User-specified title-box dimensions
TITLE_WIDTH_MM  <- 89.277
TITLE_HEIGHT_MM <- 18.238

# Font sizes in points
TITLE_FONT_PT <- 36.59
AXIS_TEXT_PT  <- 28.4
AXIS_TITLE_PT <- 28.4

# Significance P-value font size: 18 pt by default; if the original Y-axis range of a panel does not provide enough space,
# shrink the font size for that panel only, to no smaller than 10 pt.
P_VALUE_FONT_PT <- 18
P_VALUE_MIN_FONT_PT <- 10
P_VALUE_FONT_STEP_PT <- 0.5

# Physical spacing for compact annotations (mm).
# Use compact text clearance, first-bracket spacing, and top margin to preserve the original Y-axis limits whenever possible.
P_VALUE_CLEARANCE_MM <- 0.35
P_VALUE_TEXT_HEIGHT_MM <- P_VALUE_FONT_PT * 25.4 / 72 * 0.90
P_BRACKET_EXTRA_GAP_MM <- 0.25
P_BRACKET_STEP_MM <- P_VALUE_TEXT_HEIGHT_MM +
  P_VALUE_CLEARANCE_MM + P_BRACKET_EXTRA_GAP_MM
P_FIRST_BRACKET_GAP_MM <- 0.50
P_VALUE_TOP_MARGIN_MM <- 0.50

# Estimate physical text height at different P-value font sizes; actual drawing uses grid font metrics.
estimate_p_text_height_mm <- function(font_pt) {
  font_pt * 25.4 / 72 * 0.90
}

# Line widths are specified by the user in points
MAIN_LINE_PT  <- 1.788
MINOR_LINE_PT <- 0.886

# ggplot2's linewidth argument is specified in mm and is internally converted to grid's lwd;
# one grid lwd unit corresponds to 1/96 inch (0.75 pt).
# Therefore, using pt * 25.4 / 72 directly would reduce the final PDF line width to about 75% of the requested value.
# Reverse-calculate using ggplot2's .pt conversion factor (72.27/25.4) and grid's 0.75 pt/lwd
# to ensure that the physical line width in the PDF matches the specified point value.

pt_to_ggplot_linewidth_mm <- function(pt) {
  pt * 25.4 / (72.27 * 0.75)
}

# geom_text(size = ) uses mm, so the standard pt-to-mm conversion is retained
pt_to_mm <- function(pt) {
  pt * 25.4 / 72
}

MAIN_LINE_MM  <- pt_to_ggplot_linewidth_mm(MAIN_LINE_PT)
MINOR_LINE_MM <- pt_to_ggplot_linewidth_mm(MINOR_LINE_PT)

# grid::gpar(lwd) directly uses units of 1/96 inch
pt_to_grid_lwd <- function(pt) {
  pt * 96 / 72
}

MAIN_LINE_GRID <- pt_to_grid_lwd(MAIN_LINE_PT)
MINOR_LINE_GRID <- pt_to_grid_lwd(MINOR_LINE_PT)

# Use the specified light gray consistently for major and minor grid lines
GRID_MAJOR_COLOR <- "#EBEBEB"
GRID_MINOR_COLOR <- "#EBEBEB"

# Height of the X-axis label area within each cell
# 200 - 136.515 - 18.238 = 45.247 mm
X_LABEL_ZONE_MM <- CELL_HEIGHT_MM -
  PANEL_HEIGHT_MM -
  TITLE_HEIGHT_MM

if (X_LABEL_ZONE_MM <= 0) {
  stop("The page height is insufficient for the specified title box and XY plotting area.")
}

# Center the plotting area horizontally within each 190 mm cell
PANEL_LEFT_IN_CELL_MM <- (
  CELL_WIDTH_MM - PANEL_WIDTH_MM
) / 2

# Align the left and right edges of the title box exactly with the plotting area
TITLE_LEFT_IN_CELL_MM <- PANEL_LEFT_IN_CELL_MM

# Manually position Y-axis tick labels and the Y-axis title
Y_TICK_LABEL_X_OFFSET_MM <- 4
Y_AXIS_TITLE_X_OFFSET_MM <- 31

# Place X-axis labels immediately below the plotting area
X_LABEL_ANCHOR_OFFSET_MM <- 0.10


# ============================================================
# 4. Groups, colors, and pattern settings
# ============================================================

group_levels <- c(
  "CTR",
  "MP1",
  "MP10",
  "MPP1",
  "MPP10"
)

# WRS2 robust one-way ANOVA settings
TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05
LINCON_P_ADJUST_METHOD <- "none"

# Retain the previously selected colors
fill_colors <- c(
  CTR  = "#b6dcd5",
  MP1  = "#ffe0e2",
  MP10 = "#ffaab9",
  MPP1 = "#baa5d5",
  MPP10 = "#8c69b9"
)

# Retain the previously selected pattern angles
angle_values <- c(
  CTR  = 0,
  MP1  = 45,
  MP10 = 45,
  MPP1 = 135,
  MPP10 = 135
)

pattern_values <- c(
  CTR  = "stripe",
  MP1  = "stripe",
  MP10 = "stripe",
  MPP1 = "stripe",
  MPP10 = "stripe"
)

# Use full group names on the X axis
group_labels <- c(
  CTR  = "Control",
  MP1  = "Microplastic-Low",
  MP10 = "Microplastic-High",
  MPP1 = "Plastisphere-Low",
  MPP10 = "Plastisphere-High"
)


# ============================================================
# 5. Endpoint parameters and final order
#
# Order follows the continuous panel labels (a)–(j):
# Top row: (a)–(e)
# Bottom row: (f)–(j)
#
# Shorten the first panel title to Siphoning:
# The 36.59 pt font cannot fit the full title within the 89.277 mm title box
# "Siphoning_efficiency"; retain the full endpoint name and unit on the Y axis.
# ============================================================

plot_config <- list(
  
  Siphonage_efficiency = list(
    tag = "(a)",
    title = "Siphoning",
    y_label = "Siphonage efficiency (ml/clam/h)",
    ymax = 60,
    break_by = 15,
    digits = 0,
    bar_width = 0.9
  ),
  
  SOD = list(
    tag = "(b)",
    title = "SOD",
    y_label = "SOD activity (U/mgprot)",
    ymax = 600,
    break_by = 150,
    digits = 0,
    bar_width = 0.8
  ),
  
  CAT = list(
    tag = "(c)",
    title = "CAT",
    y_label = "CAT activity (U/mgprot)",
    ymax = 6,
    break_by = 1.5,
    digits = 1,
    bar_width = 0.8
  ),
  
  GSH = list(
    tag = "(d)",
    title = "GSH",
    y_label = "GSH content (μmol/gprot)",
    ymax = 20,
    break_by = 5,
    digits = 0,
    bar_width = 0.8
  ),
  
  GSH_PX = list(
    tag = "(e)",
    title = "GSH_PX",
    y_label = "GSH-Px activity (U/mgprot)",
    ymax = 12,
    break_by = 3,
    digits = 0,
    bar_width = 0.8
  ),
  
  MDA = list(
    tag = "(f)",
    title = "MDA",
    y_label = "MDA content (nmol/mgprot)",
    ymax = 2.4,
    break_by = 0.6,
    digits = 1,
    bar_width = 0.8
  ),
  
  AchE = list(
    tag = "(g)",
    title = "AchE",
    y_label = "AchE activity (U/mgprot)",
    ymax = 2.8,
    break_by = 0.7,
    digits = 1,
    bar_width = 0.8
  ),
  
  LZM = list(
    tag = "(h)",
    title = "Lysozyme",
    y_label = "Lysozyme activity (U/mgprot)",
    ymax = 28,
    break_by = 7,
    digits = 0,
    bar_width = 0.8
  ),
  
  D_Lac = list(
    tag = "(i)",
    title = "D_Lac",
    y_label = "D-Lac content (μmol/mgprot)",
    ymax = 0.32,
    break_by = 0.08,
    digits = 2,
    bar_width = 0.8
  ),
  
  DAO = list(
    tag = "(j)",
    title = "DAO",
    y_label = "DAO activity (U/g prot)",
    ymax = 28,
    break_by = 7,
    digits = 0,
    bar_width = 0.8
  )
)

plot_order <- names(plot_config)


# ============================================================
# 6. Read data and validate groups
# ============================================================

dat_raw <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

required_columns <- c(
  "group",
  plot_order
)

missing_columns <- setdiff(
  required_columns,
  names(dat_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "The raw data are missing the following columns:",
    paste(missing_columns, collapse = ", ")
  )
}

dat_raw$group <- factor(
  trimws(as.character(dat_raw$group)),
  levels = group_levels
)

if (anyNA(dat_raw$group)) {
  stop(
    "The group column contains unrecognized group labels.",
    "\nExpected: CTR, MP1, MP10, MPP1, and MPP10."
  )
}

cat("\n========== Raw data groups ==========\n")
print(table(dat_raw$group))


# ============================================================
# 7. Format P values and axis tick labels
# ============================================================

format_p_value_suffix <- function(p) {
  
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p < 0.001) {
    return("< 0.001")
  }
  
# Prevent significant P values from being rounded to P = 0.050
  if (round(p, 3) >= 0.05) {
    return(sprintf("= %.4f", p))
  }
  
  sprintf("= %.3f", p)
}

# Full text used in the CSV statistical results.
format_p_value <- function(p) {
  paste("P", format_p_value_suffix(p))
}


format_axis_value <- function(x, digits = 0) {
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


# ============================================================
# 8. Build a panel grob and statistical results for one endpoint
#
# This function generates only the XY plotting area:
#   it does not generate the default ggplot X/Y-axis text;
#   it does not generate the default ggplot title;
#   it does not generate the default ggplot outer margins.
#
# Axis titles, tick labels, and title boxes are drawn later using a manual grid layout.
# ============================================================

make_panel <- function(variable, config) {
  
  # ----------------------------------------------------------
# 8.1 Raw data
  # ----------------------------------------------------------
  
  dat <- dat_raw %>%
    transmute(
      group = factor(
        as.character(group),
        levels = group_levels
      ),
      concentration = as.numeric(.data[[variable]])
    ) %>%
    filter(
      !is.na(group),
      !is.na(concentration)
    ) %>%
    group_by(group) %>%
    mutate(
      Replicate_order = row_number(),
      
# Apply slight horizontal jitter only
# Preserve the original observed y values exactly
      Plot_x = as.numeric(first(group)) +
        c(-0.10, 0, 0.10)[Replicate_order]
    ) %>%
    ungroup()
  
# Check the number of independent biological replicates in each group
  sample_check <- dat %>%
    count(group, name = "n")
  
  if (nrow(sample_check) != 5 ||
      any(sample_check$n != 3)) {
    
    stop(
      variable,
      ": requires three valid independent biological replicates in each of the five groups.\n",
      paste(
        capture.output(print(sample_check)),
        collapse = "\n"
      )
    )
  }
  
  # ----------------------------------------------------------
  # 8.2 Mean ± SD
  # ----------------------------------------------------------
  
  summary_data <- dat %>%
    group_by(group) %>%
    summarise(
      Mean = mean(concentration),
      SD = sd(concentration),
      n = n(),
      x = as.numeric(first(group)),
      .groups = "drop"
    )
  
  # ----------------------------------------------------------
# 8.3 WRS2 t1way omnibus test
  # ----------------------------------------------------------
# t1way tests the overall difference among the five groups and allows unequal variances.
# With n = 3 per group and tr = 0.20, floor(n * tr) = 0,
# so no observations are actually trimmed; at this sample size, the test is equivalent to
# a Welch-type heteroscedastic one-way test based on the raw means.
  
  t1way_test <- WRS2::t1way(
    concentration ~ group,
    data = dat,
    tr = TRIM_LEVEL,
    alpha = ALPHA_LEVEL
  )
  
  t1way_result <- data.frame(
    Indicator = variable,
    F = unname(as.numeric(t1way_test$test)),
    df1 = unname(as.numeric(t1way_test$df1)),
    df2 = unname(as.numeric(t1way_test$df2)),
    P = unname(as.numeric(t1way_test$p.value)),
    Effect_size = unname(as.numeric(t1way_test$effsize)),
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
# 8.4 WRS2 lincon pairwise comparisons
  # ----------------------------------------------------------
# Ten pairwise comparisons are performed across the five groups.
# method = "none": lincon returns raw P values without multiple-comparison adjustment.
# Note: lincon confidence intervals remain simultaneous; only raw P values are used for annotations in this figure.
# Also export trimmed-mean differences and their simultaneous confidence intervals.
  
  lincon_test <- WRS2::lincon(
    concentration ~ group,
    data = dat,
    tr = TRIM_LEVEL,
    alpha = ALPHA_LEVEL,
    method = LINCON_P_ADJUST_METHOD
  )
  
  lincon_matrix <- as.matrix(lincon_test$comp)
  
  if (ncol(lincon_matrix) < 6 || nrow(lincon_matrix) != 10) {
    stop(
      variable,
      ": WRS2::lincon() did not return the expected 10 pairwise comparisons."
    )
  }
  
  lincon_group_names <- as.character(lincon_test$fnames)
  
  lincon_result <- data.frame(
    Indicator = variable,
    group1 = lincon_group_names[as.integer(lincon_matrix[, 1])],
    group2 = lincon_group_names[as.integer(lincon_matrix[, 2])],
    Estimate = as.numeric(lincon_matrix[, 3]),
    CI_lower = as.numeric(lincon_matrix[, 4]),
    CI_upper = as.numeric(lincon_matrix[, 5]),
    P_raw = as.numeric(lincon_matrix[, 6]),
    P_adjustment = "None",
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
# 8.5 Significant comparisons and bracket ordering
  # ----------------------------------------------------------
# Select and annotate comparisons using unadjusted raw P values only.
# Comparisons with non-overlapping horizontal spans may share the same bracket height,
# avoiding unnecessary Y-axis expansion caused by stacking every significant comparison.
  
  sig <- lincon_result %>%
    filter(P_raw < ALPHA_LEVEL) %>%
    mutate(
      x1 = match(group1, group_levels),
      x2 = match(group2, group_levels),
      span = abs(x2 - x1),
      P_suffix = vapply(
        P_raw,
        format_p_value_suffix,
        character(1)
      ),
      P_label = paste0(
        "P ",
        vapply(P_raw, format_p_value_suffix, character(1))
      )
    ) %>%
    arrange(
      x1,
      desc(x2),
      span
    )
  
# Bracket levels: only horizontally non-overlapping intervals may share the same height.
# If the end of one bracket coincides with the start of another, assign them to different levels
# to prevent the two bracket endpoints from touching.
  if (nrow(sig) > 0) {
    lane_ends <- numeric(0)
    sig$lane <- integer(nrow(sig))
    
    for (i in seq_len(nrow(sig))) {
      available_lane <- which(lane_ends < sig$x1[i])
      
      if (length(available_lane) > 0) {
        lane_id <- available_lane[1]
        lane_ends[lane_id] <- sig$x2[i]
      } else {
        lane_ends <- c(lane_ends, sig$x2[i])
        lane_id <- length(lane_ends)
      }
      
      sig$lane[i] <- lane_id
    }
    
    sig <- sig %>% arrange(lane, x1, x2)
  }
  
  # ----------------------------------------------------------
# 8.6 Determine the data top and prioritize the original Y-axis limits
  # ----------------------------------------------------------
# data_top must be defined before calculating base_y_upper and significance-bracket positions.
# Consider both the upper ends of Mean + SD error bars and all raw observations to prevent clipping.
  data_top <- max(
    c(
      summary_data$Mean + summary_data$SD,
      dat$concentration
    ),
    na.rm = TRUE
  )

  if (!is.finite(data_top)) {
    stop(variable, ": cannot calculate a valid data_top. Please check the raw data.")
  }

# First use the ymax originally specified in plot_config as the Y-axis upper limit.
# Increase the limit for the data itself only if the raw data (Mean + SD or individual observations) exceed it.
# Use the blank space above data_top for significance brackets whenever possible; if space is insufficient,
# progressively reduce the P-value font size for that panel (18 pt to 10 pt) instead of expanding the Y axis first.
# Expand the Y-axis upper limit by the minimum required amount only if the 10 pt font still does not fit.

  base_y_upper <- max(
    config$ymax,
    ceiling(data_top / config$break_by) * config$break_by
  )

  p_font_pt <- P_VALUE_FONT_PT
  p_clearance_mm <- P_VALUE_CLEARANCE_MM
  p_text_height_mm <- estimate_p_text_height_mm(p_font_pt)
  p_bracket_step_mm <- p_text_height_mm +
    p_clearance_mm + P_BRACKET_EXTRA_GAP_MM

  if (nrow(sig) > 0) {

    n_lanes <- max(sig$lane)

# Calculate the available physical height above data_top without changing the Y-axis upper limit.
    available_height_mm <- (
      (base_y_upper - data_top) / base_y_upper
    ) * PANEL_HEIGHT_MM

# Start at the default 18 pt and test progressively smaller sizes to select the largest font that fits within the original Y-axis range.
    candidate_fonts <- seq(
      P_VALUE_FONT_PT,
      P_VALUE_MIN_FONT_PT,
      by = -P_VALUE_FONT_STEP_PT
    )

    selected_font <- NA_real_
    selected_text_height <- NA_real_
    selected_step <- NA_real_

    for (font_try in candidate_fonts) {
      text_height_try <- estimate_p_text_height_mm(font_try)
      step_try <- text_height_try +
        p_clearance_mm + P_BRACKET_EXTRA_GAP_MM

      required_height_try <-
        P_FIRST_BRACKET_GAP_MM +
        (n_lanes - 1) * step_try +
        p_clearance_mm +
        text_height_try +
        P_VALUE_TOP_MARGIN_MM

      if (required_height_try <= available_height_mm) {
        selected_font <- font_try
        selected_text_height <- text_height_try
        selected_step <- step_try
        break
      }
    }

    if (is.finite(selected_font)) {
# The original ymax accommodates all annotations: keep the Y-axis limits unchanged.
      y_upper <- base_y_upper
      p_font_pt <- selected_font
      p_text_height_mm <- selected_text_height
      p_bracket_step_mm <- selected_step

    } else {
# Calculate the minimum required Y-axis expansion only if the annotations still do not fit at the minimum font size.
      p_font_pt <- P_VALUE_MIN_FONT_PT
      p_text_height_mm <- estimate_p_text_height_mm(p_font_pt)
      p_bracket_step_mm <- p_text_height_mm +
        p_clearance_mm + P_BRACKET_EXTRA_GAP_MM

      required_height_mm <-
        P_FIRST_BRACKET_GAP_MM +
        (n_lanes - 1) * p_bracket_step_mm +
        p_clearance_mm +
        p_text_height_mm +
        P_VALUE_TOP_MARGIN_MM

      reserved_fraction <- required_height_mm / PANEL_HEIGHT_MM

      if (reserved_fraction >= 0.80) {
        stop(
          variable,
          ": Too many significance annotations; even after reducing the P-value font size to",
          P_VALUE_MIN_FONT_PT,
          " pt, they still cannot be accommodated appropriately."
        )
      }

      y_upper_needed <- data_top / (1 - reserved_fraction)

# Increase the limit in smaller break_by/4 increments to avoid excessive expansion from a full break_by step.
      y_rounding_step <- config$break_by / 4
      y_upper <- max(
        base_y_upper,
        ceiling(y_upper_needed / y_rounding_step) * y_rounding_step
      )
    }

# Convert the physical bracket spacing to data coordinates using the final Y-axis range.
    first_gap_data <- (
      P_FIRST_BRACKET_GAP_MM / PANEL_HEIGHT_MM
    ) * y_upper

    bracket_step_data <- (
      p_bracket_step_mm / PANEL_HEIGHT_MM
    ) * y_upper

    label_center_offset_data <- (
      p_clearance_mm + p_text_height_mm / 2
    ) / PANEL_HEIGHT_MM * y_upper

    sig$y_position <- data_top +
      first_gap_data +
      (sig$lane - 1) * bracket_step_data

    sig$label_y <- sig$y_position + label_center_offset_data

  } else {

    sig$y_position <- numeric(0)
    sig$label_y <- numeric(0)

# When there are no significance annotations, strictly retain the original ymax unless the data themselves exceed it.
    y_upper <- base_y_upper
  }

  y_breaks <- seq(
    0,
    y_upper,
    by = config$break_by
  )
  
  y_minor_all <- seq(
    0,
    y_upper,
    by = config$break_by / 2
  )
  
# Remove minor breaks that overlap major ticks
  y_minor <- y_minor_all[
    !vapply(
      y_minor_all,
      function(v) {
        any(abs(y_breaks - v) < 1e-8)
      },
      logical(1)
    )
  ]
  
  # ----------------------------------------------------------
# 8.7 Draw the contents inside the panel
  # ----------------------------------------------------------
  
  p <- ggplot() +
    
# Bars: Mean
# Bar border width: 1.788 pt
    
    ggpattern::geom_col_pattern(
      data = summary_data,
      
      aes(
        x = x,
        y = Mean,
        fill = group,
        pattern = group,
        pattern_fill = group,
        pattern_angle = group
      ),
      
      pattern = "stripe",
      pattern_spacing = 0.12,
      pattern_density = 0.5,
      pattern_colour = "black",
      
      colour = "black",
      linewidth = MAIN_LINE_MM,
      
      width = config$bar_width
    ) +
    
# Error bars: Mean ± SD
    
    geom_errorbar(
      data = summary_data,
      
      aes(
        x = x,
        ymin = Mean - SD,
        ymax = Mean + SD
      ),
      
      width = 0.20,
      linewidth = MAIN_LINE_MM,
      colour = "black"
    ) +
    
# Raw points: white fill with black outlines
    
    geom_point(
      data = dat,
      
      aes(
        x = Plot_x,
        y = concentration
      ),
      
      inherit.aes = FALSE,
      
      shape = 21,
      fill = "white",
      colour = "black",
      size = 2.0,
      stroke = 0.45
    ) +
    
# Significance brackets: 1.788 pt, matching the main axis lines
    
    {
      if (nrow(sig) > 0) {
        list(
          geom_segment(
            data = sig,
            aes(
              x = x1,
              xend = x2,
              y = y_position,
              yend = y_position
            ),
            inherit.aes = FALSE,
            linewidth = MAIN_LINE_MM,
            colour = "black"
          ),
          
          geom_segment(
            data = sig,
            aes(
              x = x1,
              xend = x1,
              y = y_position,
              yend = y_position - y_upper * 0.012
            ),
            inherit.aes = FALSE,
            linewidth = MAIN_LINE_MM,
            colour = "black"
          ),
          
          geom_segment(
            data = sig,
            aes(
              x = x2,
              xend = x2,
              y = y_position,
              yend = y_position - y_upper * 0.012
            ),
            inherit.aes = FALSE,
            linewidth = MAIN_LINE_MM,
            colour = "black"
          )
        )
      }
    } +
    
    scale_fill_manual(
      values = fill_colors,
      drop = FALSE
    ) +
    
    scale_pattern_manual(
      values = pattern_values,
      drop = FALSE
    ) +
    
    scale_pattern_fill_manual(
      values = fill_colors,
      drop = FALSE
    ) +
    
    scale_pattern_angle_manual(
      values = angle_values,
      drop = FALSE
    ) +
    
    scale_x_continuous(
      limits = c(0.4, 5.6),
      breaks = 1:5,
      minor_breaks = seq(0.5, 5.5, by = 0.5),
      expand = expansion(mult = 0)
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
      minor_breaks = y_minor,
      expand = expansion(mult = c(0, 0))
    ) +
    
    coord_cartesian(
      xlim = c(0.4, 5.6),
      ylim = c(0, y_upper),
      expand = FALSE,
      clip = "off"
    ) +
    
    theme_void(
      base_family = font_family
    ) +
    
    theme(
      
      legend.position = "none",
      
# Major grid lines: 1.788 pt
      panel.grid.major.x = element_line(
        colour = GRID_MAJOR_COLOR,
        linewidth = MAIN_LINE_MM
      ),
      
      panel.grid.major.y = element_line(
        colour = GRID_MAJOR_COLOR,
        linewidth = MAIN_LINE_MM
      ),
      
# Minor grid lines: 0.886 pt
      panel.grid.minor.x = element_line(
        colour = GRID_MINOR_COLOR,
        linewidth = MINOR_LINE_MM
      ),
      
      panel.grid.minor.y = element_line(
        colour = GRID_MINOR_COLOR,
        linewidth = MINOR_LINE_MM
      ),
      
# XY plotting-area border: 1.788 pt
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = MAIN_LINE_MM
      ),
      
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
# All default axis text is drawn manually using grid
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      
      plot.margin = margin(
        0, 0, 0, 0,
        unit = "mm"
      )
    )
  
  # ----------------------------------------------------------
# 8.8 Extract the unique panel grob
  # ----------------------------------------------------------
  
  gt <- ggplot2::ggplotGrob(p)
  
  panel_index <- which(
    gt$layout$name == "panel"
  )
  
  if (length(panel_index) != 1) {
    stop(
      variable,
      ": Could not uniquely extract the XY plotting-area panel grob."
    )
  }
  
# Disable clipping at the panel-level gtable.
# Bars, error bars, points, grid lines, and significance brackets share the same panel grob;
# do not create separate clipping viewports for individual geom layers. P-value text is drawn separately
# in the outer page viewport and remains a native text object rather than being converted to paths by showtext.
  gt$layout$clip[panel_index] <- "off"
  
  panel_grob <- gt$grobs[[panel_index]]
  
  list(
    variable = variable,
    config = config,
    panel_grob = panel_grob,
    summary = summary_data,
    t1way = t1way_result,
    lincon = lincon_result,
    significant = sig,
    y_upper = y_upper,
    y_breaks = y_breaks,
    p_font_pt = p_font_pt,
    p_clearance_mm = p_clearance_mm,
    p_bracket_step_mm = p_bracket_step_mm
  )
}


# ============================================================
# 9. Generate all 10 endpoints
# ============================================================

all_results <- lapply(
  plot_order,
  function(nm) {
    
    message(
      "Generating panel: ",
      nm
    )
    
    make_panel(
      variable = nm,
      config = plot_config[[nm]]
    )
  }
)

names(all_results) <- plot_order


# ============================================================
# 10. Draw an individual panel manually
# ============================================================

# Significance annotations use two separate native text objects:
#   "P"：Times New Roman Italic
#   " = 0.012" / " < 0.001"：Times New Roman Regular
#
# this italicizes only the P character while keeping the remaining symbols and numbers in regular type;
# the SVG therefore contains editable text elements rather than outlined paths.

draw_editable_p_value <- function(
    center_x_mm,
    bracket_y_mm,
    suffix,
    fontsize = P_VALUE_FONT_PT,
    clearance_mm = P_VALUE_CLEARANCE_MM
) {
  
  gp_p <- gpar(
    fontfamily = font_family,
    fontface = "italic",
    fontsize = fontsize,
    col = "black"
  )
  
  gp_suffix <- gpar(
    fontfamily = font_family,
    fontface = "plain",
    fontsize = fontsize,
    col = "black"
  )
  
  grob_p <- textGrob("P", gp = gp_p)
  grob_suffix <- textGrob(suffix, gp = gp_suffix)
  
  width_p <- convertWidth(
    grobWidth(grob_p),
    "mm",
    valueOnly = TRUE
  )
  
  width_suffix <- convertWidth(
    grobWidth(grob_suffix),
    "mm",
    valueOnly = TRUE
  )
  
# Measure text height using the actual font metrics of the current output device.
# Vertically center P and its suffix; use the greater of the two heights as the total text height.
  height_p <- convertHeight(
    grobHeight(grob_p),
    "mm",
    valueOnly = TRUE
  )
  
  height_suffix <- convertHeight(
    grobHeight(grob_suffix),
    "mm",
    valueOnly = TRUE
  )
  
  text_height_mm <- max(
    height_p,
    height_suffix
  )
  
# Keep the specified clear gap between the bottom of the text and the horizontal bracket line.
  center_y_mm <- bracket_y_mm +
    clearance_mm +
    text_height_mm / 2
  
  total_width <- width_p + width_suffix
  left_x <- center_x_mm - total_width / 2
  
# P: italic, separate native text object
  grid.text(
    label = "P",
    x = unit(left_x, "mm"),
    y = unit(center_y_mm, "mm"),
    just = c("left", "centre"),
    gp = gp_p
  )
  
# Equals/less-than signs and numbers: regular type, separate native text object
  grid.text(
    label = suffix,
    x = unit(left_x + width_p, "mm"),
    y = unit(center_y_mm, "mm"),
    just = c("left", "centre"),
    gp = gp_suffix
  )
  
  invisible(
    list(
      text_height_mm = text_height_mm,
      clearance_mm = clearance_mm,
      center_y_mm = center_y_mm
    )
  )
}



#
# Each cell measures 190 × 200 mm.
#
# Vertical allocation:
#   Title box: 18.238 mm
#   XY panel: 136.515 mm
#   X-axis label area: 45.247 mm
#
# Horizontal allocation:
#   Center the panel; panel width: 89.277 mm
#   Place the Y-axis title and tick labels to the left of the panel
#
# Do not scale the entire figure to fit the cell.
# ============================================================

draw_one_chart <- function(
    result,
    cell_x,
    cell_y,
    draw_title = TRUE
) {
  
  cfg <- result$config
  
# Left edge of the panel within the cell
  panel_left <- cell_x +
    PANEL_LEFT_IN_CELL_MM
  
# Bottom edge of the panel within the cell
  panel_bottom <- cell_y +
    X_LABEL_ZONE_MM
  
  panel_top <- panel_bottom +
    PANEL_HEIGHT_MM
  
  # ----------------------------------------------------------
# 10.1 Draw the XY panel
  # ----------------------------------------------------------
  
  pushViewport(
    viewport(
      x = unit(panel_left, "mm"),
      y = unit(panel_bottom, "mm"),
      width = unit(PANEL_WIDTH_MM, "mm"),
      height = unit(PANEL_HEIGHT_MM, "mm"),
      just = c("left", "bottom"),
# Use only one positioning viewport per panel and disable clipping
      clip = "off"
    )
  )
  
  grid.draw(result$panel_grob)
  
  popViewport()
  
  
  # ----------------------------------------------------------
# 10.2 Draw Y-axis tick marks and tick labels
  # ----------------------------------------------------------
  
  y_values <- result$y_breaks
  
  for (yv in y_values) {
    
    y_fraction <- yv / result$y_upper
    
    y_position <- panel_bottom +
      y_fraction * PANEL_HEIGHT_MM
    
# Short Y-axis tick marks
    grid.segments(
      x0 = unit(
        panel_left - 1.8,
        "mm"
      ),
      x1 = unit(
        panel_left,
        "mm"
      ),
      y0 = unit(y_position, "mm"),
      y1 = unit(y_position, "mm"),
      
      gp = gpar(
        col = "black",
        lwd = MAIN_LINE_GRID
      )
    )
    
# Y-axis tick labels
    grid.text(
      label = format_axis_value(
        yv,
        cfg$digits
      ),
      
      x = unit(
        panel_left - Y_TICK_LABEL_X_OFFSET_MM,
        "mm"
      ),
      
      y = unit(
        y_position,
        "mm"
      ),
      
      just = c("right", "centre"),
      
      gp = gpar(
        fontfamily = font_family,
        fontsize = AXIS_TEXT_PT,
        col = "black"
      )
    )
  }
  
  
  # ----------------------------------------------------------
# 10.3 Y-axis title
  #
# Rotate by 90 degrees and position manually.
# Do not use the space reserved by ggplot's default axis-title layout.
  # ----------------------------------------------------------
  
  grid.text(
    label = cfg$y_label,
    
    x = unit(
      panel_left - Y_AXIS_TITLE_X_OFFSET_MM,
      "mm"
    ),
    
    y = unit(
      panel_bottom + PANEL_HEIGHT_MM / 2,
      "mm"
    ),
    
    rot = 90,
    
    just = "centre",
    
    gp = gpar(
      fontfamily = font_family,
      fontsize = AXIS_TITLE_PT,
      fontface = "bold",
      col = "black"
    )
  )
  
  
  # ----------------------------------------------------------
# 10.4 Significance P-value text
  #
# Bracket lines are drawn inside the panel and use the same width as the main axis.
# Draw P-value text in the page coordinate system of the same panel:
#   italicize P;
#   keep numbers and comparison symbols in regular type;
#   adapt the font size to the panel and maintain a compact gap below the bracket line.
  # ----------------------------------------------------------
  
  if (nrow(result$significant) > 0) {
    
    for (j in seq_len(nrow(result$significant))) {
      
      one_sig <- result$significant[j, ]
      
      x_mid <- (one_sig$x1 + one_sig$x2) / 2
      
      x_fraction <- (x_mid - 0.4) / 5.2
      
      p_x_mm <- panel_left +
        x_fraction * PANEL_WIDTH_MM
      
# Actual position of the horizontal line in page coordinates.
# The P-value text function measures the actual font height and positions the text using the compact gap specified for the panel.
      bracket_y_mm <- panel_bottom +
        (one_sig$y_position / result$y_upper) *
        PANEL_HEIGHT_MM
      
      draw_editable_p_value(
        center_x_mm = p_x_mm,
        bracket_y_mm = bracket_y_mm,
        suffix = one_sig$P_suffix,
        fontsize = result$p_font_pt,
        clearance_mm = result$p_clearance_mm
      )
    }
  }
  
  
  # ----------------------------------------------------------
# 10.5 Draw X-axis tick marks and group labels rotated by 30°
  # ----------------------------------------------------------
  
  x_values <- 1:5
  
  x_positions <- panel_left +
    (
      (x_values - 0.4) / 5.2
    ) * PANEL_WIDTH_MM
  
  for (i in seq_along(x_values)) {
    
# Short X-axis tick marks
    grid.segments(
      x0 = unit(x_positions[i], "mm"),
      x1 = unit(x_positions[i], "mm"),
      
      y0 = unit(panel_bottom, "mm"),
      y1 = unit(panel_bottom - 1.8, "mm"),
      
      gp = gpar(
        col = "black",
        lwd = MAIN_LINE_GRID
      )
    )
    
# Group labels
# 28.4 pt, rotated by 30 degrees
# Align the top-right corner with the corresponding tick position
    grid.text(
      label = unname(
        group_labels[group_levels[i]]
      ),
      
      x = unit(
        x_positions[i],
        "mm"
      ),
      
      y = unit(
        panel_bottom - X_LABEL_ANCHOR_OFFSET_MM,
        "mm"
      ),
      
      just = c("right", "top"),
      
      rot = 30,
      
      gp = gpar(
        fontfamily = font_family,
        fontsize = AXIS_TEXT_PT,
        col = "black"
      )
    )
  }
  
  
  # ----------------------------------------------------------
# 10.6 Title box
  #
# Only the title has an outer border.
# The title box has exactly the same width as the XY panel.
  # ----------------------------------------------------------
  
  if (draw_title) {
    
    title_center_x <- panel_left +
      TITLE_WIDTH_MM / 2
    
    title_center_y <- panel_top +
      TITLE_HEIGHT_MM / 2
    
    grid.rect(
      x = unit(title_center_x, "mm"),
      y = unit(title_center_y, "mm"),
      
      width = unit(TITLE_WIDTH_MM, "mm"),
      height = unit(TITLE_HEIGHT_MM, "mm"),
      
      gp = gpar(
        fill = "white",
        col = "black",
        lwd = MAIN_LINE_GRID
      )
    )
    
    grid.text(
      label = paste(
        cfg$tag,
        cfg$title
      ),
      
      x = unit(title_center_x, "mm"),
      y = unit(title_center_y, "mm"),
      
      just = "centre",
      
      gp = gpar(
        fontfamily = font_family,
        fontsize = TITLE_FONT_PT,
        fontface = "bold",
        col = "black"
      )
    )
  }
  
  invisible(NULL)
}


# ============================================================
# 11. Manual full-page layout function
# ============================================================

draw_full_page <- function(newpage = TRUE) {
  
  if (isTRUE(newpage)) {
    grid.newpage()
  }
  
  for (i in seq_along(plot_order)) {
    
    row_id <- ceiling(i / N_COLS)
    
    col_id <- ((i - 1) %% N_COLS) + 1
    
# The top row is placed in the upper half of the page
# The bottom row is placed in the lower half of the page
    cell_x <- (col_id - 1) * CELL_WIDTH_MM
    
    cell_y <- PAGE_HEIGHT_MM -
      row_id * CELL_HEIGHT_MM
    
    draw_one_chart(
      result = all_results[[plot_order[i]]],
      cell_x = cell_x,
      cell_y = cell_y
    )
  }
  
  invisible(NULL)
}


# ============================================================
# 12. Check physical panel dimensions and text dimensions
#
# Do not measure the natural dimensions of the complete ggplot gtable here.
# Check the manually specified physical dimensions and the actual text widths directly.
# ============================================================

measure_text_mm <- function(
    label,
    fontsize,
    fontface = "plain"
) {
  
  tg <- textGrob(
    label,
    gp = gpar(
      fontfamily = font_family,
      fontsize = fontsize,
      fontface = fontface
    )
  )
  
  c(
    width_mm = convertWidth(
      grobWidth(tg),
      "mm",
      valueOnly = TRUE
    ),
    
    height_mm = convertHeight(
      grobHeight(tg),
      "mm",
      valueOnly = TRUE
    )
  )
}


check_text_layout <- function() {
  
  cat("\n========== Fixed physical dimensions ==========\n")
  
  cat(
    "Page: ",
    PAGE_WIDTH_MM,
    " × ",
    PAGE_HEIGHT_MM,
    " mm\n"
  )
  
  cat(
    "Cell: ",
    CELL_WIDTH_MM,
    " × ",
    CELL_HEIGHT_MM,
    " mm\n"
  )
  
  cat(
    "XY plotting area: ",
    PANEL_WIDTH_MM,
    " × ",
    PANEL_HEIGHT_MM,
    " mm\n"
  )
  
  cat(
    "Title box: ",
    TITLE_WIDTH_MM,
    " × ",
    TITLE_HEIGHT_MM,
    " mm\n"
  )
  
  cat(
    "X-axis label area: ",
    round(X_LABEL_ZONE_MM, 3),
    " mm\n"
  )
  
  cat("\n========== Title text width check ==========\n")
  
  title_metrics <- lapply(
    all_results,
    function(res) {
      
      label <- paste(
        res$config$tag,
        res$config$title
      )
      
      measure_text_mm(
        label,
        TITLE_FONT_PT,
        "bold"
      )
    }
  )
  
  title_table <- data.frame(
    Indicator = names(all_results),
    Title_width_mm = vapply(
      title_metrics,
      function(x) x["width_mm"],
      numeric(1)
    )
  )
  
  print(title_table)
  
  if (any(title_table$Title_width_mm > TITLE_WIDTH_MM)) {
    warning(
      "Some title text exceeds the title-box width. Please check the title labels."
    )
  }
  
  cat("\n========== Y-axis title length check ==========\n")
  
  y_title_metrics <- lapply(
    all_results,
    function(res) {
      measure_text_mm(
        res$config$y_label,
        AXIS_TITLE_PT,
        "bold"
      )
    }
  )
  
  y_title_table <- data.frame(
    Indicator = names(all_results),
    Rotated_title_length_mm = vapply(
      y_title_metrics,
      function(x) x["width_mm"],
      numeric(1)
    )
  )
  
  print(y_title_table)
  
  if (any(
    y_title_table$Rotated_title_length_mm >
    PANEL_HEIGHT_MM
  )) {
    warning(
      "Some rotated Y-axis titles exceed the height of the XY plotting area."
    )
  }
  
  cat("\n========== X-axis label area check ==========\n")
  
  x_metrics <- lapply(
    unname(group_labels[group_levels]),
    function(label) {
      measure_text_mm(
        label,
        AXIS_TEXT_PT,
        "plain"
      )
    }
  )
  
# Vertical projection after 30° rotation:
  # width*sin(30°) + height*cos(30°)
  
  x_projected_height <- vapply(
    x_metrics,
    function(x) {
      x["width_mm"] * sin(pi / 6) +
        x["height_mm"] * cos(pi / 6)
    },
    numeric(1)
  )
  
  x_table <- data.frame(
    Group = unname(group_labels[group_levels]),
    Rotated_height_mm = x_projected_height,
    Available_height_mm = X_LABEL_ZONE_MM
  )
  
  print(x_table)
  
  if (any(x_projected_height > X_LABEL_ZONE_MM)) {
    warning(
      "Some rotated X-axis labels exceed the reserved vertical space.",
      "\nPlease check that the bottom labels are fully visible in the final PDF."
    )
  }
  
  invisible(
    list(
      title = title_table,
      y_axis = y_title_table,
      x_axis = x_table
    )
  )
}


# ============================================================
# 13. Automatic backup function for the current R script
#
# When this script is open in RStudio and run using Source,
# automatically copy the current R file to the output directory.
# ============================================================

get_current_script_path <- function() {
  
  candidates <- character(0)
  
# Currently edited file in RStudio
  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()
  ) {
    
    p <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    
    if (length(p) == 1 && nzchar(p)) {
      candidates <- c(candidates, p)
    }
  }
  
# source() call stack
  frames <- sys.frames()
  
  for (fr in rev(frames)) {
    if (!is.null(fr$ofile)) {
      candidates <- c(
        candidates,
        as.character(fr$ofile)
      )
    }
  }
  
  # Rscript --file
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  
  if (length(file_arg) > 0) {
    candidates <- c(
      candidates,
      sub("^--file=", "", file_arg[1])
    )
  }
  
  candidates <- unique(
    candidates[nzchar(candidates)]
  )
  
  for (p in candidates) {
    if (file.exists(p)) {
      return(
        normalizePath(
          p,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }
  
  NULL
}


script_backup_path <- file.path(
  output_dir,
  "Fig1_WRS2_script_backup.R"
)

current_script_path <- get_current_script_path()

if (!is.null(current_script_path)) {
  
  copied <- file.copy(
    from = current_script_path,
    to = script_backup_path,
    overwrite = TRUE
  )
  
  if (copied) {
    cat(
      "\nR script backup completed successfully:\n",
      normalizePath(
        script_backup_path,
        winslash = "/"
      ),
      "\n"
    )
  } else {
    warning("Automatic R script backup failed. Please copy the current R file manually.")
  }
  
} else {
  
  warning(
    "Could not automatically identify the path of the current R script.\n",
    "Open this R file in RStudio and click Source to run it; ",
    "the script will then be backed up to the output directory."
  )
}


# ============================================================
# 14. Output file paths
# ============================================================

# Main Illustrator-editable file: native SVG text objects
svg_file <- file.path(
  output_dir,
  "Fig1_editable.svg"
)

# Auxiliary PDF: Cairo output with preserved font text objects
pdf_file <- file.path(
  output_dir,
  "Fig1_editable.pdf"
)

png_file <- file.path(
  output_dir,
  "Fig1_preview.png"
)


# ============================================================
# 15. Export PDF
#
# Cairo PDF：
#   Preserve vector graphics and font text objects;
#   do not call showtext;
#   do not convert text to outlines in advance.
# ============================================================

grDevices::cairo_pdf(
  filename = pdf_file,
  
  width = PAGE_WIDTH_MM / 25.4,
  height = PAGE_HEIGHT_MM / 25.4,
  
  family = font_family,
  symbolfamily = font_family,
  pointsize = 12,
  bg = "white"
)

tryCatch(
  {
    grid.newpage()
    check_text_layout()
    draw_full_page(newpage = FALSE)
  },
  finally = {
    grDevices::dev.off()
  }
)

# Optionally check whether fonts are embedded in the PDF.
# pdftools is used only to read PDF font information and does not modify the figure.
if (requireNamespace("pdftools", quietly = TRUE)) {
  
  pdf_font_info <- pdftools::pdf_fonts(pdf_file)
  
  cat("\n========== PDF font embedding check ==========\n")
  print(pdf_font_info)
  
  write.csv(
    pdf_font_info,
    file.path(output_dir, "PDF_Font_Embedding_Check.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  if ("embedded" %in% names(pdf_font_info) &&
      any(!pdf_font_info$embedded)) {
    warning(
      "The PDF contains fonts that are not embedded. Check that Times New Roman is installed correctly ",
      "and that Cairo is calling this font successfully."
    )
  }
  
} else {
  message(
    "To check PDF font embedding, install pdftools and run the script again:",
    "install.packages('pdftools')"
  )
}


# ============================================================
# 16. Export Illustrator-editable SVG (primary deliverable)
#
# Use svglite to export native SVG <text> elements:
#   * Do not use showtext;
#   * Do not convert Times New Roman text to paths;
#   * Titles, axes, and P values remain SVG text objects;
#   * Italicize only the P character; keep the numeric part in regular type.
#
# In Illustrator, use File > Open to open this SVG,
# rather than selecting the earlier font-outline/vector-path file.
# Times New Roman must be installed on the computer.
# ============================================================

svglite::svglite(
  file = svg_file,
  width = PAGE_WIDTH_MM / 25.4,
  height = PAGE_HEIGHT_MM / 25.4,
  bg = "white",
  pointsize = 12,
  standalone = TRUE
)

tryCatch(
  {
    draw_full_page()
  },
  finally = {
    grDevices::dev.off()
  }
)

# Check whether the SVG contains native text elements.
svg_lines <- readLines(
  svg_file,
  warn = FALSE,
  encoding = "UTF-8"
)

svg_text_count <- sum(
  grepl("<text([[:space:]>])", svg_lines)
)

svg_font_family_present <- any(
  grepl("Times New Roman", svg_lines, fixed = TRUE)
)

cat("\n========== Illustrator SVG text-object check ==========\n")
cat("Number of native SVG text elements: ", svg_text_count, "\n")
cat("Contains a Times New Roman font declaration: ",
    svg_font_family_present, "\n")

if (svg_text_count == 0) {
  stop(
    "No native <text> elements were detected in the SVG.",
    "\nDo not use this file as the Illustrator-editable text version."
  )
}

if (!svg_font_family_present) {
  warning(
    "No Times New Roman font declaration was detected in the SVG. ",
    "Check the font mapping of the current R graphics device."
  )
}

write.csv(
  data.frame(
    SVG_text_elements_detected = svg_text_count,
    Times_New_Roman_declared = svg_font_family_present,
    SVG_file = svg_file
  ),
  file.path(output_dir, "Illustrator_SVG_Text_Object_Check.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 17. Export PNG preview
# ============================================================

grDevices::png(
  filename = png_file,
  
  width = PAGE_WIDTH_MM,
  height = PAGE_HEIGHT_MM,
  
  units = "mm",
  res = 300,
  
  type = "cairo",
  family = font_family,
  pointsize = 12,
  bg = "white"
)

tryCatch(
  {
    draw_full_page()
  },
  finally = {
    grDevices::dev.off()
  }
)


# ============================================================
# 18. Export 10 individual PDF files
#
# Page size for each individual panel: 190 × 200 mm
# The panel and title box retain their specified physical dimensions.
# ============================================================

single_dir <- file.path(
  output_dir,
  "Individual_Endpoint_PDFs"
)

dir.create(
  single_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (nm in plot_order) {
  
  single_file <- file.path(
    single_dir,
    paste0(nm, ".pdf")
  )
  
  grDevices::cairo_pdf(
    filename = single_file,
    
    width = CELL_WIDTH_MM / 25.4,
    height = CELL_HEIGHT_MM / 25.4,
    
    family = font_family,
    symbolfamily = font_family,
    pointsize = 12,
    bg = "white"
  )
  
  tryCatch(
    {
      grid.newpage()
      
      draw_one_chart(
        result = all_results[[nm]],
        cell_x = 0,
        cell_y = 0
      )
    },
    finally = {
      grDevices::dev.off()
    }
  )
}


# ============================================================
# 19. Summarize and export statistical results
# ============================================================

t1way_results <- bind_rows(
  lapply(
    all_results,
    function(x) x$t1way
  )
)

lincon_results <- bind_rows(
  lapply(
    all_results,
    function(x) x$lincon
  )
)

significant_lincon_results <- lincon_results %>%
  filter(P_raw < ALPHA_LEVEL) %>%
  mutate(
    P_label = vapply(
      P_raw,
      format_p_value,
      character(1)
    )
  )

summary_results <- bind_rows(
  lapply(
    names(all_results),
    function(nm) {
      all_results[[nm]]$summary %>%
        mutate(Indicator = nm)
    }
  )
)

write.csv(
  t1way_results,
  file.path(
    output_dir,
    "WRS2_t1way_Overall_Test_Results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  lincon_results,
  file.path(
    output_dir,
    "WRS2_lincon_All_10_Pairwise_Comparisons_Unadjusted_P.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  significant_lincon_results,
  file.path(
    output_dir,
    "WRS2_lincon_Significant_Comparisons_Unadjusted_P.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  summary_results,
  file.path(
    output_dir,
    "Group_Mean_SD_Summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 20. Save plotting parameter records
# ============================================================

# Record the final Y-axis upper limit, bracket levels, and P-value font size for each panel separately,
# to identify panels that retained the original ymax and those adjusted slightly because of data or annotation-space constraints.
y_axis_annotation_table <- bind_rows(
  lapply(names(all_results), function(nm) {
    x <- all_results[[nm]]
    data.frame(
      Indicator = nm,
      Configured_Ymax = x$config$ymax,
      Final_Ymax = x$y_upper,
      Significant_comparisons = nrow(x$significant),
      Bracket_lanes = if (nrow(x$significant) > 0) max(x$significant$lane) else 0,
      P_value_font_pt = if (nrow(x$significant) > 0) x$p_font_pt else NA_real_,
      Bracket_step_mm = if (nrow(x$significant) > 0) x$p_bracket_step_mm else NA_real_,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  y_axis_annotation_table,
  file.path(output_dir, "Endpoint_Y_Axis_Limits_and_Significance_Parameters.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

parameter_table <- data.frame(
  Parameter = c(
    "Page width",
    "Page height",
    "Layout",
    "Cell width",
    "Cell height",
    "Panel width",
    "Panel height",
    "Title frame width",
    "Title frame height",
    "Title font",
    "Y-axis tick font",
    "Y-axis title font",
    "X-axis tick font",
    "X-axis angle",
    "Main line width",
    "Minor line width",
    "Major grid color",
    "Minor grid color",
    "Significance bracket width",
    "Panel clipping",
    "Font family",
    "PDF device",
    "Illustrator editable SVG",
    "Cairo PDF text",
    "P-value font",
    "P character style",
    "P-value to bracket clearance",
    "Adjacent bracket line spacing",
    "Bracket lane allocation",
    "Scatter",
    "Bar summary",
    "Statistical test",
    "Multiple comparison correction",
    "Significance display"
  ),
  
  Value = c(
    "950 mm",
    "400 mm",
    "5 columns × 2 rows",
    "190 mm",
    "200 mm",
    "89.277 mm",
    "136.515 mm",
    "89.277 mm",
    "18.238 mm",
    "36.59 pt, bold",
    "28.4 pt",
    "28.4 pt, bold",
    "28.4 pt",
    "30 degrees",
    "1.788 pt",
    "0.886 pt",
    "#EBEBEB",
    "#EBEBEB",
    "1.788 pt",
    "off",
    "Times New Roman",
    "Cairo PDF",
    "svglite native SVG <text> elements",
    "Cairo PDF embedded font text",
    "18 pt default; adaptive down to 10 pt only when needed to preserve configured Y-axis range",
    "Italic P; regular operator and values",
    paste0(P_VALUE_CLEARANCE_MM, " mm; panel-specific"),
    "Panel-specific; based on selected P-value font size and available Y-axis headroom",
    "Non-overlapping x-intervals share a lane; overlapping intervals use separate lanes",
    "White fill, black outline",
    "Mean ± SD",
    "WRS2 t1way (tr = 0.20) + lincon",
    "No multiple-comparison adjustment for 10 pairwise comparisons per indicator",
    "Unadjusted raw P values, no stars"
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  parameter_table,
  file.path(
    output_dir,
    "Final_Plot_Parameter_Record.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 21. Final output summary
# ============================================================

cat("\n")
cat("====================================================\n")
cat("Physiological and biochemical endpoint figure: manual physical-layout version completed\n")
cat("====================================================\n")

cat(
  "Page: ",
  PAGE_WIDTH_MM,
  " × ",
  PAGE_HEIGHT_MM,
  " mm\n"
)

cat("Layout: 5 columns × 2 rows\n")

cat(
  "XY plotting area: ",
  PANEL_WIDTH_MM,
  " × ",
  PANEL_HEIGHT_MM,
  " mm\n"
)

cat(
  "Title box: ",
  TITLE_WIDTH_MM,
  " × ",
  TITLE_HEIGHT_MM,
  " mm\n"
)

cat("Font: Times New Roman\n")
cat("Primary Illustrator-editable file: native SVG text objects from svglite\n")
cat("Auxiliary PDF: Cairo PDF with preserved font text objects\n")
cat("showtext is disabled; no font-outline paths are exported\n")
cat("P-value font size: ", P_VALUE_FONT_PT, " pt by default; adaptively reduced based on available space in each panel, down to ", P_VALUE_MIN_FONT_PT, " pt; only the P character is italicized\n", sep = "")
cat("Y-axis limits: retain the original ymax in plot_config whenever possible; increase only minimally when required by the data or annotations\n")
cat("Major line width: ", MAIN_LINE_PT, " pt; minor line width: ", MINOR_LINE_PT, " pt\n")
cat("Major/minor grid line color: #EBEBEB\n")
cat("Significance bracket line width: ", MAIN_LINE_PT, " pt (same as the main axis lines)\n")
cat("Panel clipping: disabled; no separate clipping viewports are created for individual geom layers\n")
cat("Bar height: Mean; error bars: SD\n")
cat("Points: white fill with black outlines\n")
cat("Statistics: WRS2 t1way() + lincon(); tr = ", TRIM_LEVEL, "\n", sep = "")
cat("Pairwise comparisons: all 10 pairs; no multiple-comparison adjustment of P values\n")
cat("Significance annotations: raw unadjusted P values; no asterisks\n")
cat("Final Y-axis limits and P-value font sizes for each endpoint: see Endpoint_Y_Axis_Limits_and_Significance_Parameters.csv\n")
cat("Note: with n = 3 per group, tr = 0.20 does not trim any observations\n\n")

cat("Output directory:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/"
  ),
  "\n"
)

cat("\nPDF：\n")
cat(pdf_file, "\n")

cat("\nIllustrator-editable SVG with native text objects:\n")
cat(svg_file, "\n")

cat("\nPNG：\n")
cat(png_file, "\n")

cat("\nR script backup:\n")
cat(script_backup_path, "\n")

cat("====================================================\n")