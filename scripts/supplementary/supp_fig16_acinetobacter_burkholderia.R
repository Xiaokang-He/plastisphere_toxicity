# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.


# ============================================================
# Acinetobacter and Burkholderia species composition
#
# Figure 1: g__Acinetobacter
# Figure 2: g__Burkholderia-Caballeronia-Paraburkholderia
#
# Plot:
#   Stacked bar: Species-level mean relative abundance
#   Error bars: Total genus abundance, Mean +/- SD
#   Flow: Species-level alluvial flow
#
# Statistics:
#   Overall: WRS2 t1way()
#   Pairwise: WRS2 lincon()
#   Trim level: tr = 0.20
#   P adjustment: method = "none"
#   No overall significance gate
#   No multiple-comparison adjustment
#   Annotation: raw pairwise P < 0.05
#   P >= 0.001: italic P, 3 decimal places
#   P < 0.001: P < 0.001
#
# Analysis unit:
#   Total genus abundance of each independent sample
#
# Output:
#   Two individual figures
#   One combined figure (2 columns x 1 row)
#   WRS2 statistical result tables
# ============================================================


# -------------------- 0. Clean environment -------------------

rm(list = ls())


# -------------------- 1. Load packages -----------------------

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "ggalluvial",
  "ggsignif",
  "cowplot",
  "WRS2"
)

for (pkg in packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(ggsignif)
library(cowplot)
library(WRS2)


if (utils::packageVersion("WRS2") < "1.1.3") {
  stop(
    "Please update WRS2 to version 1.1-3 or later. ",
    "The method argument in lincon() is required for unadjusted P values."
  )
}


# -------------------- 2. Paths and data ----------------------

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- file.path(input_dir, "Acinetobacter丰度表.csv")
output_dir <- file.path(results_dir, "supp_fig16_acinetobacter_burkholderia")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_columns <- c(
  "Genus",
  "Species",
  "Abundance",
  "group",
  "Sample",
  "Sample2"
)

missing_columns <- setdiff(
  required_columns,
  names(data)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# Clean character columns.

data$Genus <- trimws(
  as.character(data$Genus)
)

data$Species <- trimws(
  as.character(data$Species)
)

data$group <- trimws(
  as.character(data$group)
)


# Convert abundance to numeric.

abundance_original <- data$Abundance

data$Abundance <- suppressWarnings(
  as.numeric(as.character(data$Abundance))
)

if (
  any(
    is.na(data$Abundance) &
    !is.na(abundance_original)
  )
) {
  stop(
    "Abundance contains non-numeric values. ",
    "Please check the original abundance table."
  )
}


data <- data %>%
  filter(
    !is.na(Abundance),
    is.finite(Abundance)
  )


# Experimental group order.

groups <- c(
  "CTR",
  "MP1",
  "MP10",
  "MPP1",
  "MPP10"
)

if (
  any(
    !data$group %in% groups
  )
) {
  stop(
    "Unexpected group names: ",
    paste(
      unique(data$group[!data$group %in% groups]),
      collapse = ", "
    )
  )
}


data$group <- factor(
  data$group,
  levels = groups
)


# Target genera.

genus_acinetobacter <- "g__Acinetobacter"

genus_burkholderia <- "g__Burkholderia-Caballeronia-Paraburkholderia"

genus_targets <- c(
  genus_acinetobacter,
  genus_burkholderia
)


if (
  !all(
    genus_targets %in% unique(data$Genus)
  )
) {
  stop(
    "One or both target genera are missing from the input data."
  )
}


# ============================================================
# 3. Prepare sample-level total genus abundance
# ============================================================

# IMPORTANT:
# Statistical tests are performed using total abundance
# of each genus in each independent biological sample.
#
# Species-level abundance is summed within each sample.
#
# This is different from testing species-level group means.
#
# Samples in which the target genus is absent are assigned
# an abundance of zero rather than being omitted.


sample_info <- data %>%
  distinct(
    Sample,
    Sample2,
    group
  )


genus_sample_totals <- data %>%
  filter(
    Genus %in% genus_targets
  ) %>%
  group_by(
    Genus,
    Sample,
    Sample2,
    group
  ) %>%
  summarise(
    Abundance = sum(
      Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Complete all sample x genus combinations.

genus_sample_data <- tidyr::crossing(
  sample_info,
  Genus = genus_targets
) %>%
  left_join(
    genus_sample_totals,
    by = c(
      "Genus",
      "Sample",
      "Sample2",
      "group"
    )
  ) %>%
  mutate(
    Abundance = coalesce(
      Abundance,
      0
    ),
    group = factor(
      group,
      levels = groups
    )
  ) %>%
  arrange(
    Genus,
    group,
    Sample,
    Sample2
  )


# Check sample sizes.

WRS2_group_n <- genus_sample_data %>%
  group_by(
    Genus,
    group
  ) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )


cat(
  "\n========== Sample sizes ==========\n"
)

print(
  as.data.frame(WRS2_group_n),
  row.names = FALSE
)


if (any(WRS2_group_n$n < 2)) {
  
  stop(
    "At least one genus x group combination has fewer than 2 samples. ",
    "Please check the sample identifiers and original abundance table."
  )
}


# ============================================================
# 4. WRS2 t1way() + lincon()
# ============================================================

TRIM_LEVEL <- 0.20

ALPHA_LEVEL <- 0.05


run_wrs2 <- function(genus_name) {
  
  d <- genus_sample_data %>%
    filter(
      Genus == genus_name
    ) %>%
    mutate(
      group = factor(
        as.character(group),
        levels = groups
      )
    )
  
  
  # Check all five groups.
  
  group_n <- table(
    factor(
      d$group,
      levels = groups
    )
  )
  
  
  if (any(group_n < 2)) {
    
    stop(
      "Insufficient samples for genus: ",
      genus_name
    )
  }
  
  
  # ----------------------------------------------------------
  # 4.1 Overall robust one-way ANOVA
  #
  # WRS2 t1way()
  # tr = 0.20
  #
  # This overall test is calculated and reported.
  # It is NOT used as a prerequisite for pairwise comparisons.
  # ----------------------------------------------------------
  
  overall_fit <- WRS2::t1way(
    Abundance ~ group,
    data = d,
    tr = TRIM_LEVEL
  )
  
  
  overall_result <- data.frame(
    
    Genus = genus_name,
    
    Test = "WRS2_t1way",
    
    Trim_level = TRIM_LEVEL,
    
    Statistic_F = as.numeric(
      overall_fit$test
    ),
    
    df1 = as.numeric(
      overall_fit$df1
    ),
    
    df2 = as.numeric(
      overall_fit$df2
    ),
    
    P_t1way = as.numeric(
      overall_fit$p.value
    ),
    
    Effect_size = as.numeric(
      overall_fit$effsize
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 4.2 All pairwise robust comparisons
  #
  # WRS2 lincon()
  # tr = 0.20
  # method = "none"
  #
  # No multiple-comparison adjustment of P values.
  # All 10 pairwise comparisons are retained.
  # ----------------------------------------------------------
  
  lincon_fit <- WRS2::lincon(
    Abundance ~ group,
    data = d,
    tr = TRIM_LEVEL,
    alpha = ALPHA_LEVEL,
    method = "none"
  )
  
  
  comp <- as.matrix(
    lincon_fit$comp
  )
  
  
  if (ncol(comp) != 6) {
    
    stop(
      "Unexpected lincon result structure for ",
      genus_name,
      ". Please check your WRS2 version."
    )
  }
  
  
  # WRS2 lincon output:
  #
  # Column 1: First group index
  # Column 2: Second group index
  # Column 3: Trimmed mean difference
  # Column 4: Simultaneous CI lower limit
  # Column 5: Simultaneous CI upper limit
  # Column 6: P value after the specified method
  #
  # method = "none" means column 6 contains unadjusted P values.
  
  
  group_names <- as.character(
    lincon_fit$fnames
  )
  
  
  group1_index <- as.integer(
    comp[, 1]
  )
  
  group2_index <- as.integer(
    comp[, 2]
  )
  
  
  if (
    anyNA(group1_index) ||
    anyNA(group2_index)
  ) {
    stop(
      "Invalid group indices in lincon results for ",
      genus_name
    )
  }
  
  
  lincon_result <- data.frame(
    
    Genus = genus_name,
    
    Group1 = group_names[group1_index],
    
    Group2 = group_names[group2_index],
    
    psihat = as.numeric(
      comp[, 3]
    ),
    
    CI_lower_simultaneous = as.numeric(
      comp[, 4]
    ),
    
    CI_upper_simultaneous = as.numeric(
      comp[, 5]
    ),
    
    P_raw = as.numeric(
      comp[, 6]
    ),
    
    P_adjustment = "none",
    
    Trim_level = TRIM_LEVEL,
    
    P_t1way = overall_result$P_t1way[1],
    
    stringsAsFactors = FALSE
  )
  
  
  # Check all 10 comparisons.
  
  if (
    nrow(lincon_result) != choose(
      length(groups),
      2
    )
  ) {
    
    stop(
      "lincon did not return all 10 pairwise comparisons for ",
      genus_name
    )
  }
  
  
  # Check P values.
  
  if (
    anyNA(lincon_result$P_raw) ||
    any(
      lincon_result$P_raw < 0 |
      lincon_result$P_raw > 1
    )
  ) {
    
    stop(
      "Invalid P values in lincon results for ",
      genus_name
    )
  }
  
  
  list(
    overall = overall_result,
    pairwise = lincon_result
  )
}


# Run the analyses separately for both genera.

wrs_acinetobacter <- run_wrs2(
  genus_acinetobacter
)

wrs_burkholderia <- run_wrs2(
  genus_burkholderia
)


# Combine overall results.

t1way_results <- bind_rows(
  wrs_acinetobacter$overall,
  wrs_burkholderia$overall
)


# Combine all pairwise results.

lincon_results <- bind_rows(
  wrs_acinetobacter$pairwise,
  wrs_burkholderia$pairwise
)


# ------------------------------------------------------------
# IMPORTANT:
#
# No overall t1way significance gate.
# No P-value adjustment.
#
# Only raw pairwise P < 0.05 determines whether
# a comparison is displayed on the figure.
# ------------------------------------------------------------

lincon_significant <- lincon_results %>%
  filter(
    P_raw < ALPHA_LEVEL
  )


cat(
  "\n========== WRS2 t1way overall results ==========\n"
)

print(
  as.data.frame(t1way_results),
  row.names = FALSE,
  digits = 6
)


cat(
  "\n========== WRS2 lincon all pairwise results ==========\n"
)

print(
  as.data.frame(lincon_results),
  row.names = FALSE,
  digits = 6
)


cat(
  "\n========== Significant pairwise results ==========\n"
)

print(
  as.data.frame(lincon_significant),
  row.names = FALSE,
  digits = 6
)


# ============================================================
# 5. P-value label formatting
# ============================================================

format_p_label <- function(p) {
  
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p < 0.001) {
    
    return(
      "italic(P) < 0.001"
    )
  }
  
  sprintf(
    "italic(P) == %.3f",
    p
  )
}


# ============================================================
# 6. Species-level plotting function
# ============================================================

make_genus_plot <- function(
    
  genus_name,
  
  species_levels,
  
  fill_values,
  
  flow_alpha,
  
  y_axis_max,
  
  y_axis_breaks,
  
  bracket_offset,
  
  bracket_step,
  
  output_prefix
  
) {
  
  
  # ----------------------------------------------------------
  # 6.1 Filter genus-specific data
  # ----------------------------------------------------------
  
  df <- data %>%
    filter(
      Genus == genus_name
    ) %>%
    group_by(
      group,
      Species
    ) %>%
    filter(
      sum(
        Abundance,
        na.rm = TRUE
      ) != 0
    ) %>%
    ungroup()
  
  
  if (nrow(df) == 0) {
    
    stop(
      "No plotting data for ",
      genus_name
    )
  }
  
  
  unknown_species <- setdiff(
    unique(as.character(df$Species)),
    species_levels
  )
  
  
  if (length(unknown_species) > 0) {
    
    stop(
      "Unexpected species in ",
      genus_name,
      ": ",
      paste(unknown_species, collapse = ", ")
    )
  }
  
  
  # Group order.
  
  df$group <- factor(
    as.character(df$group),
    levels = groups
  )
  
  
  # Original species order.
  
  df$Species <- factor(
    df$Species,
    levels = species_levels
  )
  
  
  # ----------------------------------------------------------
  # 6.2 Species-level mean and SD
  #
  # The stacked bars show the arithmetic mean abundance
  # of each species within each group.
  # Species-level SD is retained in the exported summary,
  # but is NOT used to draw the genus-total error bars.
  # ----------------------------------------------------------
  
  dt <- df %>%
    group_by(
      group,
      Species,
      Genus
    ) %>%
    summarise(
      n = sum(!is.na(Abundance)),
      m = mean(Abundance, na.rm = TRUE),
      sd = sd(Abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  dt$group <- factor(
    dt$group,
    levels = groups
  )
  
  
  # Reverse species order for stacking,
  # consistent with the original code.
  
  dt$Species <- factor(
    dt$Species,
    levels = rev(species_levels)
  )
  
  
  # Calculate cumulative positions.
  
  dt <- dt %>%
    arrange(
      group,
      Species
    ) %>%
    group_by(
      group,
      Genus
    ) %>%
    mutate(
      position = cumsum(m)
    ) %>%
    ungroup()
  
  
  # ----------------------------------------------------------
  # 6.2.1 Total genus abundance: Mean +/- SD
  #
  # IMPORTANT:
  # Calculate the error bars from the total abundance of
  # each independent sample, i.e. the sum of all species
  # belonging to this genus within that sample.
  #
  # This is the same sample-level quantity used by WRS2.
  # One error bar is drawn per experimental group.
  # ----------------------------------------------------------
  
  genus_total_summary <- genus_sample_data %>%
    filter(
      Genus == genus_name
    ) %>%
    group_by(
      group
    ) %>%
    summarise(
      n = n(),
      total_mean = mean(Abundance, na.rm = TRUE),
      total_sd = sd(Abundance, na.rm = TRUE),
      mean_plus_sd = total_mean + total_sd,
      .groups = "drop"
    ) %>%
    mutate(
      group = factor(
        as.character(group),
        levels = groups
      )
    ) %>%
    arrange(group)
  
  
  if (
    nrow(genus_total_summary) != length(groups) ||
    anyNA(genus_total_summary$total_mean) ||
    anyNA(genus_total_summary$total_sd)
  ) {
    stop(
      "Incomplete total-genus summary for ",
      genus_name,
      ". Check sample-level abundance data."
    )
  }
  
  
  # Export total-genus Mean +/- SD used for error bars.
  
  write.csv(
    genus_total_summary,
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_total_genus_mean_SD.csv"
      )
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  
  # Export species-level descriptive statistics.
  
  write.csv(
    dt,
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_species_mean_SD.csv"
      )
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  
  # ----------------------------------------------------------
  # 6.3 Select significant pairwise comparisons
  #
  # No overall t1way gate.
  # No P adjustment.
  # ----------------------------------------------------------
  
  ann <- lincon_significant %>%
    filter(
      Genus == genus_name
    ) %>%
    mutate(
      comparison_order = match(
        paste(
          Group1,
          Group2,
          sep = "__"
        ),
        
        vapply(
          combn(
            groups,
            2,
            simplify = FALSE
          ),
          
          function(z) {
            paste(z, collapse = "__")
          },
          
          character(1)
        )
      )
    ) %>%
    arrange(
      comparison_order
    )
  
  
  # ----------------------------------------------------------
  # 6.4 Bracket positions
  #
  # Preserve the original genus-specific spacing.
  # ----------------------------------------------------------
  
  # Use the same total-genus Mean + SD values that define
  # the tops of the error bars when positioning brackets.
  
  y_start <- max(
    genus_total_summary$mean_plus_sd,
    na.rm = TRUE
  )
  
  
  if (nrow(ann) > 0) {
    
    ann$y_position <- pmax(
      0,
      y_start -
        bracket_offset +
        (seq_len(nrow(ann)) - 1) *
        bracket_step
    )
    
    
    ann$p_label <- vapply(
      ann$P_raw,
      format_p_label,
      character(1)
    )
    
  }
  
  
  # Extend the plotting range when required by either
  # the total-genus Mean + SD error bars or P-value brackets.
  
  plot_ymax <- max(
    y_axis_max,
    max(genus_total_summary$mean_plus_sd, na.rm = TRUE)
  )
  
  if (nrow(ann) > 0) {
    
    plot_ymax <- max(
      plot_ymax,
      max(ann$y_position) + y_axis_max * 0.04
    )
    
  }
  
  
  # ----------------------------------------------------------
  # 6.5 Main alluvial stacked bar plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    
    data = dt,
    
    aes(
      x = group,
      y = m,
      fill = Species,
      stratum = Species,
      alluvium = Species
    )
    
  ) +
    
    
    # Stacked bars.
    
    geom_stratum(
      alpha = 1,
      lwd = 0.45,
      width = 0.6
    ) +
    
    
    # Species flow.
    
    geom_flow(
      alpha = flow_alpha,
      width = 0.6
    ) +
    
    
    # Original species colors.
    
    scale_fill_manual(
      values = fill_values,
      drop = FALSE
    ) +
    
    
    # One total-genus Mean +/- SD error bar per group.
    # These values are calculated from sample-level total
    # genus abundance, not from individual stacked species.
    
    geom_errorbar(
      data = genus_total_summary,
      mapping = aes(
        x = group,
        ymin = total_mean - total_sd,
        ymax = total_mean + total_sd
      ),
      inherit.aes = FALSE,
      width = 0.2,
      linewidth = 0.45
    ) +
    
    
    theme_test() +
    
    
    labs(
      y = "Relative abundance (%)",
      x = ""
    ) +
    
    
    scale_y_continuous(
      
      expand = expansion(
        mult = c(0, 0)
      ),
      
      breaks = y_axis_breaks
      
    ) +
    
    
    coord_cartesian(
      ylim = c(0, plot_ymax),
      clip = "off"
    ) +
    
    
    theme_bw() +
    
    
    theme(
      
      legend.position = "none",
      
      axis.text.x = element_text(
        size = 20,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = 20,
        colour = "black"
      ),
      
      axis.title.y = element_text(
        size = 26,
        face = "bold"
      ),
      
      axis.title.x = element_text(
        size = 26,
        face = "bold"
      ),
      
      plot.title = element_text(
        size = 26,
        face = "bold"
      )
      
    )
  
  
  # ----------------------------------------------------------
  # 6.6 Add WRS2 raw P-value brackets
  #
  # Only raw P < 0.05 is displayed.
  #
  # The annotation labels are not significance stars.
  # ----------------------------------------------------------
  
  if (nrow(ann) > 0) {
    
    p <- p +
      
      geom_signif(
        
        data = ann,
        
        mapping = aes(
          
          xmin = Group1,
          
          xmax = Group2,
          
          annotations = p_label,
          
          y_position = y_position
          
        ),
        
        manual = TRUE,
        
        inherit.aes = FALSE,
        
        parse = TRUE,
        
        map_signif_level = FALSE,
        
        tip_length = 0.015,
        
        color = "black",
        
        vjust = 0.65,
        
        size = 0.35,
        
        textsize = 6,
        
        margin_top = 0,
        
        step_increase = 0
        
      )
    
  }
  
  
  return(
    list(
      plot = p,
      summary = dt,
      annotations = ann
    )
  )
  
}


# ============================================================
# 7. Figure 1: Acinetobacter
# ============================================================

species_acinetobacter <- c(
  
  "s__Acinetobacter_johnsonii",
  
  "s__unclassified_g__Acinetobacter",
  
  "s__Acinetobacter_lwoffii",
  
  "s__Acinetobacter_pittii"
  
)


# Original palette.

mycol_acinetobacter <- rev(
  c(
    "#f2f6f9",
    "#ccdae8",
    "#92afcc",
    "#266099"
  )
)


# Original factor order is reversed for plotting.

fill_acinetobacter <- setNames(
  rev(mycol_acinetobacter),
  rev(species_acinetobacter)
)


result_acinetobacter <- make_genus_plot(
  
  genus_name = genus_acinetobacter,
  
  species_levels = species_acinetobacter,
  
  fill_values = fill_acinetobacter,
  
  flow_alpha = 0.30,
  
  y_axis_max = 1.05,
  
  y_axis_breaks = seq(
    0,
    1,
    0.25
  ),
  
  bracket_offset = 0.03,
  
  bracket_step = 0.03,
  
  output_prefix = "Acinetobacter"
  
)


p3_1 <- result_acinetobacter$plot

print(p3_1)


# ============================================================
# 8. Figure 2: Burkholderia-Caballeronia-Paraburkholderia
# ============================================================

species_burkholderia <- c(
  
  "s__Burkholderia_cenocepacia",
  
  "s__unclassified_g__Burkholderia-Caballeronia-Paraburkholderia"
  
)


# Original palette.

mycol_burkholderia <- rev(
  c(
    "#85d4e3",
    "#5aa6b3"
  )
)


fill_burkholderia <- setNames(
  rev(mycol_burkholderia),
  rev(species_burkholderia)
)


result_burkholderia <- make_genus_plot(
  
  genus_name = genus_burkholderia,
  
  species_levels = species_burkholderia,
  
  fill_values = fill_burkholderia,
  
  flow_alpha = 0.20,
  
  y_axis_max = 0.063,
  
  y_axis_breaks = seq(
    0,
    0.06,
    0.015
  ),
  
  bracket_offset = 0.008,
  
  bracket_step = 0.002,
  
  output_prefix = "Burkholderia"
  
)


p3_2 <- result_burkholderia$plot

print(p3_2)


# ============================================================
# 9. Combine two panels
# ============================================================

plots <- list(
  p3_1,
  p3_2
)


combined_plot <- cowplot::plot_grid(
  
  plots[[1]],
  
  plots[[2]],
  
  ncol = 2,
  
  nrow = 1,
  
  rel_widths = c(1, 1),
  
  rel_heights = c(1, 1)
  
)


print(combined_plot)


# ============================================================
# 10. Export figures
# ============================================================

# Individual figures.

ggsave(
  
  filename = file.path(
    output_dir,
    "Acinetobacter_WRS2_robust.pdf"
  ),
  
  plot = p3_1,
  
  units = "mm",
  
  height = 170,
  
  width = 180,
  
  bg = "white"
  
)


ggsave(
  
  filename = file.path(
    output_dir,
    "Burkholderia_WRS2_robust.pdf"
  ),
  
  plot = p3_2,
  
  units = "mm",
  
  height = 170,
  
  width = 180,
  
  bg = "white"
  
)


# Combined figure: 2 columns x 1 row.

ggsave(
  
  filename = file.path(
    output_dir,
    "Acinetobacter_Burkholderia_WRS2_robust_combined.pdf"
  ),
  
  plot = combined_plot,
  
  units = "mm",
  
  height = 170,
  
  width = 360,
  
  bg = "white"
  
)


ggsave(
  
  filename = file.path(
    output_dir,
    "Acinetobacter_Burkholderia_WRS2_robust_combined.png"
  ),
  
  plot = combined_plot,
  
  units = "mm",
  
  height = 170,
  
  width = 360,
  
  dpi = 600,
  
  bg = "white"
  
)


# ============================================================
# 11. Export statistical results
# ============================================================


# 11.1 Overall WRS2 t1way results.

write.csv(
  
  t1way_results,
  
  file.path(
    output_dir,
    "WRS2_t1way_overall_results.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# 11.2 All 10 pairwise comparisons for each genus.

write.csv(
  
  lincon_results,
  
  file.path(
    output_dir,
    "WRS2_lincon_all_pairwise_rawP_no_adjustment.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# 11.3 Significant pairwise comparisons.
#
# Selection:
#   P_raw < 0.05 ONLY.
#
# No overall t1way gate.
# No P-value adjustment.

write.csv(
  
  lincon_significant,
  
  file.path(
    output_dir,
    "WRS2_lincon_significant_rawP.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# 11.4 Sample sizes.

write.csv(
  
  WRS2_group_n,
  
  file.path(
    output_dir,
    "WRS2_group_sample_sizes.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# 11.5 Sample-level genus abundance used for statistical tests.

write.csv(
  
  genus_sample_data,
  
  file.path(
    output_dir,
    "WRS2_sample_level_genus_abundance.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# 11.6 Statistical annotation summary.

annotation_summary <- bind_rows(
  
  result_acinetobacter$annotations,
  
  result_burkholderia$annotations
  
)


write.csv(
  
  annotation_summary,
  
  file.path(
    output_dir,
    "WRS2_significant_comparisons_displayed.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# ============================================================
# 12. Statistical analysis notes
# ============================================================

analysis_notes <- c(
  
  "Overall test: WRS2 t1way().",
  
  "Pairwise test: WRS2 lincon().",
  
  "Trimmed mean proportion: tr = 0.20.",
  
  "All five experimental groups were included.",
  
  "All ten pairwise comparisons were calculated separately for each genus.",
  
  "P-value adjustment: method = none.",
  
  "No overall t1way significance gate was applied to pairwise comparisons.",
  
  "Pairwise annotations were based exclusively on raw lincon P < 0.05.",
  
  "P values >= 0.001 were displayed to three decimal places.",
  
  "P values < 0.001 were displayed as P < 0.001.",
  
  "The statistical analysis used sample-level total genus abundance.",
  
  "Species-level stacked bars show arithmetic mean abundance.",
  
  "Error bars show standard deviation (SD) of sample-level total genus abundance.",
  
  "WRS2 lincon confidence intervals remain simultaneous confidence intervals; method = none applies to P-value adjustment."
)


writeLines(
  
  analysis_notes,
  
  file.path(
    output_dir,
    "WRS2_statistical_analysis_notes.txt"
  ),
  
  useBytes = TRUE
  
)


# ============================================================
# 13. Final output summary
# ============================================================

cat(
  "\n============================================\n",
  "WRS2 analysis completed.\n",
  "============================================\n"
)


cat(
  "\nOverall t1way results:\n"
)

print(
  as.data.frame(t1way_results),
  row.names = FALSE,
  digits = 6
)


cat(
  "\nAll pairwise lincon results:\n"
)

print(
  as.data.frame(lincon_results),
  row.names = FALSE,
  digits = 6
)


cat(
  "\nSignificant comparisons (raw P < 0.05):\n"
)

print(
  as.data.frame(lincon_significant),
  row.names = FALSE,
  digits = 6
)


cat(
  "\nOutput directory:\n",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)