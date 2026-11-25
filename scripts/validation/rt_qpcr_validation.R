# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# RT-qPCR grouped horizontal bar plots
# Descriptive statistics: Mean +/- SD
# Statistics: WRS2 t1way() + WRS2 lincon()
# 20% trimmed-mean robust one-way ANOVA and all pairwise contrasts
# Raw pairwise P values (method = "none"; no multiplicity adjustment)
# Significant pairwise comparisons are annotated with brackets
# and actual P values rather than letters or asterisks.
#
# Input columns:
#   Gene | group1 | relative_CT
# Groups:
#   Ctr, MP1, MP10, MPP1, MPP10
# ============================================================

rm(list = ls())

# -------------------- 1. Paths and packages ------------------
repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- file.path(input_dir, "rt-qPCR.csv")
output_dir <- file.path(results_dir, "rt_qpcr_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("WRS2", quietly = TRUE)) {
  install.packages("WRS2", repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(patchwork)

# WRS2 robust-analysis settings used consistently in the other figures.
TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

data2 <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c("Gene", "group1", "relative_CT")
missing_cols <- setdiff(required_cols, names(data2))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -------------------- 2. Clean and validate data --------------
data2$Gene <- trimws(as.character(data2$Gene))
data2$group1 <- trimws(as.character(data2$group1))
data2$relative_CT <- suppressWarnings(
  as.numeric(as.character(data2$relative_CT))
)

if (any(!is.na(data2$relative_CT) & !is.finite(data2$relative_CT))) {
  stop(
    "The relative_CT column contains Inf or -Inf values. ",
    "Please check the input CSV before calculating summary statistics."
  )
}

data2 <- data2[
  !is.na(data2$Gene) &
    !is.na(data2$group1) &
    !is.na(data2$relative_CT),
  ,
  drop = FALSE
]

groups_stat <- c("Ctr", "MP1", "MP10", "MPP1", "MPP10")

unexpected_groups <- setdiff(unique(data2$group1), groups_stat)
if (length(unexpected_groups) > 0) {
  stop(
    "Unexpected group names in group1: ",
    paste(unexpected_groups, collapse = ", "),
    ". Expected: ",
    paste(groups_stat, collapse = ", ")
  )
}

data2$group1 <- factor(data2$group1, levels = groups_stat)

genes_left <- c("CxI", "GST-2", "CAM", "NPC", "Profiln", "F-ACTIN")
genes_right <- c(
  "Femtin", "CASPASE-9", "CASPASE-7",
  "Cathepsins", "NAGA", "FADD"
)
genes_all <- c(genes_left, genes_right)

missing_genes <- setdiff(genes_all, unique(data2$Gene))
if (length(missing_genes) > 0) {
  stop(
    "The following target genes are missing from the input file: ",
    paste(missing_genes, collapse = ", ")
  )
}

if (any(data2$relative_CT < 0, na.rm = TRUE)) {
  warning(
    "Some relative_CT values are negative. ",
    "Check whether this column contains relative expression values."
  )
}

# Check biological replicate number in each gene x group.
replicate_check <- aggregate(
  relative_CT ~ Gene + group1,
  data = data2,
  FUN = length
)
names(replicate_check)[3] <- "n"
if (any(replicate_check$n != 3)) {
  warning(
    "Some Gene x group combinations do not contain exactly 3 observations. ",
    "The analysis will use the available non-missing observations."
  )
}

# -------------------- 3. Mean and SD --------------------------
summary_list <- lapply(genes_all, function(g) {
  do.call(rbind, lapply(groups_stat, function(gr) {
    v <- data2$relative_CT[
      data2$Gene == g & data2$group1 == gr
    ]
    
    if (length(v) == 0) {
      stop("No observations for Gene = ", g, ", group = ", gr)
    }
    
    data.frame(
      Gene = g,
      group1 = gr,
      n = length(v),
      mean = mean(v),
      SD = if (length(v) > 1) sd(v) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
})

result_summary <- do.call(rbind, summary_list)

# -------------------- 4. WRS2 robust test function ----------------
# Each target gene is analyzed separately.
# Overall test: WRS2::t1way(), tr = 0.20.
# Pairwise tests: WRS2::lincon(), tr = 0.20, method = "none".
# All 10 pairwise comparisons are performed for every gene.
# Plot annotation requires BOTH:
#   overall t1way P < 0.05 AND pairwise raw lincon P < 0.05.

run_wrs2_one_gene <- function(df_gene, gene_name) {

  df_gene <- df_gene[
    !is.na(df_gene$relative_CT) &
      is.finite(df_gene$relative_CT) &
      !is.na(df_gene$group1),
    ,
    drop = FALSE
  ]

  df_gene$group1 <- factor(
    as.character(df_gene$group1),
    levels = groups_stat
  )

  group_n <- table(
    factor(df_gene$group1, levels = groups_stat)
  )

  if (any(group_n < 2)) {
    stop(
      "Gene = ", gene_name,
      " does not have at least 2 valid observations in every group. ",
      "Current sample sizes: ",
      paste(
        names(group_n),
        as.integer(group_n),
        sep = " = ",
        collapse = "; "
      )
    )
  }

  # ----------------------------------------------------------
  # 4.1 Overall robust test: WRS2::t1way()
  # ----------------------------------------------------------
  t1way_fit <- tryCatch(
    WRS2::t1way(
      relative_CT ~ group1,
      data = df_gene,
      tr = TRIM_LEVEL
    ),
    error = function(e) {
      stop(
        "WRS2::t1way() failed for Gene = ", gene_name, ":\n",
        conditionMessage(e)
      )
    }
  )

  t1way_result <- data.frame(
    Gene = gene_name,
    Test = "WRS2 t1way",
    Trim_level = TRIM_LEVEL,
    Statistic_F = as.numeric(t1way_fit$test),
    df1 = as.numeric(t1way_fit$df1),
    df2 = as.numeric(t1way_fit$df2),
    p_t1way = as.numeric(t1way_fit$p.value),
    Effect_size = as.numeric(t1way_fit$effsize),
    stringsAsFactors = FALSE
  )

  if (
    length(t1way_result$p_t1way) != 1 ||
      !is.finite(t1way_result$p_t1way) ||
      t1way_result$p_t1way < 0 ||
      t1way_result$p_t1way > 1
  ) {
    stop(
      "WRS2::t1way() returned an invalid overall P value for Gene = ",
      gene_name, "."
    )
  }

  # ----------------------------------------------------------
  # 4.2 All pairwise robust contrasts: WRS2::lincon()
  # method = "none": no multiple-comparison adjustment.
  # The first two columns of $comp are group indices;
  # column 6 contains the unadjusted P value.
  # ----------------------------------------------------------
  lincon_fit <- tryCatch(
    WRS2::lincon(
      relative_CT ~ group1,
      data = df_gene,
      tr = TRIM_LEVEL,
      alpha = ALPHA_LEVEL,
      method = "none"
    ),
    error = function(e) {
      stop(
        "WRS2::lincon() failed for Gene = ", gene_name, ":\n",
        conditionMessage(e)
      )
    }
  )

  comp <- as.data.frame(lincon_fit$comp)

  if (ncol(comp) < 6) {
    stop(
      "WRS2::lincon() returned fewer than 6 columns for Gene = ",
      gene_name, "; cannot extract all pairwise contrasts."
    )
  }

  group_index_1 <- as.integer(comp[[1]])
  group_index_2 <- as.integer(comp[[2]])

  group_names <- lincon_fit$fnames

  if (
    is.null(group_names) ||
      length(group_names) != length(groups_stat)
  ) {
    group_names <- levels(
      droplevels(df_gene$group1)
    )
  }

  if (
    anyNA(group_index_1) ||
      anyNA(group_index_2) ||
      any(
        group_index_1 < 1 |
          group_index_1 > length(group_names)
      ) ||
      any(
        group_index_2 < 1 |
          group_index_2 > length(group_names)
      )
  ) {
    stop(
      "Could not map WRS2::lincon() group indices back to the ",
      "five experimental groups for Gene = ", gene_name, "."
    )
  }

  lincon_result <- data.frame(
    Gene = gene_name,
    Group1 = as.character(group_names[group_index_1]),
    Group2 = as.character(group_names[group_index_2]),
    psihat = as.numeric(comp[[3]]),
    ci_lower = as.numeric(comp[[4]]),
    ci_upper = as.numeric(comp[[5]]),
    P_raw = as.numeric(comp[[6]]),
    p_t1way = t1way_result$p_t1way[1],
    Test = "WRS2 lincon",
    Trim_level = TRIM_LEVEL,
    P_adjustment = "none",
    stringsAsFactors = FALSE
  )

  expected_comparisons <- choose(length(groups_stat), 2)

  if (nrow(lincon_result) != expected_comparisons) {
    stop(
      "WRS2::lincon() did not return all ",
      expected_comparisons,
      " pairwise comparisons for Gene = ", gene_name,
      "; returned ", nrow(lincon_result), "."
    )
  }

  lincon_result$x1 <- match(
    lincon_result$Group1,
    groups_stat
  )
  lincon_result$x2 <- match(
    lincon_result$Group2,
    groups_stat
  )

  if (
    anyNA(lincon_result$x1) ||
      anyNA(lincon_result$x2) ||
      any(!is.finite(lincon_result$P_raw)) ||
      any(
        lincon_result$P_raw < 0 |
          lincon_result$P_raw > 1
      )
  ) {
    stop(
      "WRS2::lincon() returned invalid group labels or P values for Gene = ",
      gene_name, "."
    )
  }

  # Confirm that all unique pairwise combinations are present.
  pair_key <- function(a, b) {
    paste(sort(c(a, b)), collapse = "__")
  }

  expected_pairs <- combn(
    groups_stat,
    2,
    FUN = function(z) pair_key(z[1], z[2])
  )

  observed_pairs <- vapply(
    seq_len(nrow(lincon_result)),
    function(i) {
      pair_key(
        lincon_result$Group1[i],
        lincon_result$Group2[i]
      )
    },
    character(1)
  )

  if (!setequal(expected_pairs, observed_pairs)) {
    stop(
      "WRS2::lincon() did not return the complete set of unique ",
      "pairwise comparisons for Gene = ", gene_name, "."
    )
  }

  lincon_result <- lincon_result[
    order(
      lincon_result$x2 - lincon_result$x1,
      lincon_result$x1
    ),
    ,
    drop = FALSE
  ]

  list(
    t1way = t1way_result,
    lincon = lincon_result,
    group_n = data.frame(
      Gene = gene_name,
      Group = names(group_n),
      n = as.integer(group_n),
      stringsAsFactors = FALSE
    )
  )
}


# -------------------- 5. Run WRS2 tests for each gene ------------
test_results <- lapply(genes_all, function(g) {
  df_gene <- data2[
    data2$Gene == g,
    ,
    drop = FALSE
  ]
  run_wrs2_one_gene(df_gene, g)
})

names(test_results) <- genes_all

t1way_results <- do.call(
  rbind,
  lapply(genes_all, function(g) {
    test_results[[g]]$t1way
  })
)

lincon_results <- do.call(
  rbind,
  lapply(genes_all, function(g) {
    test_results[[g]]$lincon
  })
)

WRS2_group_n <- do.call(
  rbind,
  lapply(genes_all, function(g) {
    test_results[[g]]$group_n
  })
)

# All pairwise tests are retained in lincon_results.
# Only comparisons meeting BOTH significance conditions are plotted.
lincon_significant <- lincon_results[
  !is.na(lincon_results$P_raw) &
    lincon_results$P_raw < ALPHA_LEVEL &
    !is.na(lincon_results$p_t1way) &
    lincon_results$p_t1way < ALPHA_LEVEL,
  ,
  drop = FALSE
]

cat("\n========== WRS2 t1way overall tests for all genes ==========\n")
print(
  as.data.frame(t1way_results),
  row.names = FALSE
)

cat("\n========== WRS2 lincon all pairwise comparisons ==========\n")
print(
  as.data.frame(lincon_results),
  row.names = FALSE
)

cat("\n========== Comparisons displayed on plots ==========\n")
print(
  as.data.frame(lincon_significant),
  row.names = FALSE
)


# -------------------- 6. Plotting positions -------------------
# The plotting order is reversed so that, after coord_flip(),
# the visible top-to-bottom order is:
# Control, MP-Low, MP-High, Plastisphere-Low, Plastisphere-High.

groups_plot <- rev(groups_stat)

fill_cols <- c(
  "Ctr"  = "#b6dcd5",
  "MP1"  = "#ffe0e2",
  "MP10" = "#ffaab9",
  "MPP1" = "#baa5d5",
  "MPP10"= "#8c69b9"
)

group_labels <- c(
  "Ctr" = "Control",
  "MP1" = "Microplastic-Low",
  "MP10" = "Microplastic-High",
  "MPP1" = "Plastisphere-Low",
  "MPP10" = "Plastisphere-High"
)

# Bar centers are manually positioned to keep the exact group order
# and to make the lincon brackets align with the corresponding bars.
dodge_span <- 0.72
group_offsets <- seq(
  -dodge_span / 2,
  dodge_span / 2,
  length.out = length(groups_plot)
)
names(group_offsets) <- groups_plot

bar_width <- 0.13

# -------------------- 7. P-value label formatting --------------
format_p_label <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) {
    return("italic(P) < 0.001")
  }
  paste0(
    "italic(P) == ",
    formatC(p, format = "f", digits = 3)
  )
}

# -------------------- 8. Construct bracket coordinates --------
make_brackets <- function(gene_set, gene_index_map, summary_data,
                          lincon_data) {

  # Annotate only when BOTH the overall t1way test and the
  # corresponding unadjusted lincon pairwise test are significant.
  sig <- lincon_data[
    !is.na(lincon_data$P_raw) &
      lincon_data$P_raw < ALPHA_LEVEL &
      !is.na(lincon_data$p_t1way) &
      lincon_data$p_t1way < ALPHA_LEVEL,
    ,
    drop = FALSE
  ]
  
  sig <- sig[sig$Gene %in% gene_set, , drop = FALSE]
  
  if (nrow(sig) == 0) {
    return(data.frame(
      Gene = character(),
      x1 = numeric(),
      x2 = numeric(),
      xmin = numeric(),
      xmax = numeric(),
      xmid = numeric(),
      y_position = numeric(),
      y_label = numeric(),
      tick = numeric(),
      p_label = character(),
      P_raw = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  global_top <- max(
    summary_data$mean + summary_data$SD,
    na.rm = TRUE
  )
  
  # Dynamic spacing is based on the scale of the expression values.
  start_gap <- max(0.20, global_top * 0.045)
  bracket_step <- max(0.28, global_top * 0.075)
  tick_size <- max(0.07, global_top * 0.018)
  label_gap <- max(0.06, global_top * 0.018)
  
  bracket_list <- list()
  
  for (g in gene_set) {
    
    sg <- sig[sig$Gene == g, , drop = FALSE]
    if (nrow(sg) == 0) next
    
    # Longer comparisons are placed closer to the bars;
    # shorter nested comparisons are placed progressively farther out.
    sg$span <- abs(
      group_offsets[sg$Group1] -
        group_offsets[sg$Group2]
    )
    
    sg <- sg[
      order(-sg$span, sg$P_raw),
      ,
      drop = FALSE
    ]
    
    gene_top <- max(
      summary_data$mean[summary_data$Gene == g] +
        summary_data$SD[summary_data$Gene == g],
      na.rm = TRUE
    )
    
    gene_x <- unname(gene_index_map[g])
    
    sg$x1 <- gene_x + unname(group_offsets[sg$Group1])
    sg$x2 <- gene_x + unname(group_offsets[sg$Group2])
    sg$xmin <- pmin(sg$x1, sg$x2)
    sg$xmax <- pmax(sg$x1, sg$x2)
    sg$xmid <- (sg$xmin + sg$xmax) / 2
    
    sg$y_position <- gene_top + start_gap +
      (seq_len(nrow(sg)) - 1) * bracket_step
    
    sg$y_label <- sg$y_position + label_gap
    sg$tick <- tick_size
    sg$p_label <- vapply(
      sg$P_raw,
      format_p_label,
      character(1)
    )
    
    bracket_list[[g]] <- sg[, c(
      "Gene", "x1", "x2", "xmin", "xmax", "xmid",
      "y_position", "y_label", "tick", "p_label", "P_raw"
    )]
  }
  
  if (length(bracket_list) == 0) {
    return(data.frame(
      Gene = character(),
      x1 = numeric(),
      x2 = numeric(),
      xmin = numeric(),
      xmax = numeric(),
      xmid = numeric(),
      y_position = numeric(),
      y_label = numeric(),
      tick = numeric(),
      p_label = character(),
      P_raw = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  do.call(rbind, bracket_list)
}

# -------------------- 9. Plot function ------------------------
make_gene_plot <- function(gene_set, plot_title = NULL) {
  
  # The first visible gene is at the top after coord_flip().
  gene_levels <- rev(gene_set)
  
  gene_index_map <- seq_along(gene_levels)
  names(gene_index_map) <- gene_levels
  
  plot_data <- result_summary[
    result_summary$Gene %in% gene_set,
    ,
    drop = FALSE
  ]
  
  plot_data$group1 <- factor(
    plot_data$group1,
    levels = groups_plot
  )
  
  plot_data$Gene <- factor(
    plot_data$Gene,
    levels = gene_levels
  )
  
  plot_data$gene_index <- unname(
    gene_index_map[as.character(plot_data$Gene)]
  )
  
  plot_data$plot_x <- plot_data$gene_index +
    unname(group_offsets[as.character(plot_data$group1)])
  
  brackets <- make_brackets(
    gene_set = gene_set,
    gene_index_map = gene_index_map,
    summary_data = plot_data,
    lincon_data = lincon_results
  )
  
  observed_top <- max(
    plot_data$mean + plot_data$SD,
    na.rm = TRUE
  )
  
  if (nrow(brackets) > 0) {
    observed_top <- max(
      observed_top,
      brackets$y_label,
      na.rm = TRUE
    )
  }
  
  # Preserve the full mean +/- SD interval. The lower axis remains at zero
  # unless an SD interval genuinely extends below zero.
  observed_bottom <- min(
    plot_data$mean - plot_data$SD,
    na.rm = TRUE
  )
  y_lower <- min(0, observed_bottom)
  
  # Extra room to the right of the longest P-value label.
  y_upper <- observed_top + max(0.35, observed_top * 0.10)
  
  p <- ggplot() +
    
    # Mean bars.
    geom_col(
      data = plot_data,
      aes(
        x = plot_x,
        y = mean,
        fill = group1
      ),
      width = bar_width,
      colour = "black",
      linewidth = 0.30
    ) +
    
    # Mean +/- SD error bars.
    geom_errorbar(
      data = plot_data,
      aes(
        x = plot_x,
        ymin = mean - SD,
        ymax = mean + SD
      ),
      width = 0.075,
      linewidth = 0.30,
      colour = "black",
      na.rm = TRUE
    ) +
    
    # Main horizontal bracket after coord_flip:
    # connects the two compared treatment bars.
    geom_segment(
      data = brackets,
      aes(
        x = xmin,
        xend = xmax,
        y = y_position,
        yend = y_position
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    
    # Left bracket cap.
    geom_segment(
      data = brackets,
      aes(
        x = xmin,
        xend = xmin,
        y = y_position,
        yend = y_position - tick
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    
    # Right bracket cap.
    geom_segment(
      data = brackets,
      aes(
        x = xmax,
        xend = xmax,
        y = y_position,
        yend = y_position - tick
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    
    # Actual, unadjusted WRS2 lincon P values; labels are vertical.
    geom_text(
      data = brackets,
      aes(
        x = xmid,
        y = y_label,
        label = p_label
      ),
      inherit.aes = FALSE,
      parse = TRUE,
      angle = 270,         # Rotate vertical P-value labels by 180 degrees
      size = 2.0,
      family = "Times New Roman",
      colour = "black",
      hjust = 0.5,
      vjust = 0.5
    ) +
    
    scale_fill_manual(
      name = "Experimental group",
      values = fill_cols,
      breaks = groups_stat,
      labels = group_labels[groups_stat],
      drop = FALSE
    ) +
    
    scale_x_continuous(
      breaks = seq_along(gene_levels),
      labels = gene_levels,
      limits = c(0.5, length(gene_levels) + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = pretty(c(y_lower, y_upper), n = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    
    coord_flip(clip = "off") +
    
    labs(
      x = "Target gene",
      y = expression("Relative mRNA expression (" * 2^{-Delta*Delta*Ct} * ")"),
      title = plot_title
    ) +
    
    theme_bw(
      base_size = 12,
      base_family = "Times New Roman"
    ) +
    
    theme(
      panel.grid.major = element_line(
        colour = "grey88",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.65
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 13,
        colour = "black"
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 13,
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 10.5,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 11,
        colour = "black"
      ),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold",
        size = 10.5
      ),
      legend.text = element_text(
        size = 9.5,
        colour = "black"
      ),
      legend.key = element_rect(
        colour = "black",
        linewidth = 0.25
      ),
      legend.key.size = grid::unit(0.42, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      plot.margin = margin(
        t = 8, r = 14, b = 8, l = 8
      )
    )
  
  return(p)
}

# -------------------- 10. Generate the two panels --------------
p1 <- make_gene_plot(
  genes_left,
  plot_title = NULL
)

p2 <- make_gene_plot(
  genes_right,
  plot_title = NULL
)

print(p1)
print(p2)

# Combined figure with one shared legend.
p_combined <- patchwork::wrap_plots(
  p1, p2,
  ncol = 2,
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

print(p_combined)

# -------------------- 11. Save figures -------------------------
# Use Cairo PDF when available; otherwise fall back to the standard PDF device.
pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

# Individual panels.
ggsave(
  file.path(output_dir, "RT-qPCR_panel_1.pdf"),
  p1,
  width = 85,
  height = 170,
  units = "mm",
  device = pdf_device,
  bg = "white"
)

ggsave(
  file.path(output_dir, "RT-qPCR_panel_2.pdf"),
  p2,
  width = 85,
  height = 170,
  units = "mm",
  device = pdf_device,
  bg = "white"
)

# Combined two-panel figure.
# Cairo successfully exported the two individual panels but failed on the
# patchwork-combined object. Therefore, try Cairo first and automatically
# fall back to the standard vector PDF device if Cairo cannot finish the file.

combined_pdf <- file.path(output_dir, "RT-qPCR_combined.pdf")

# Remove an incomplete PDF left by a previous failed export.
if (file.exists(combined_pdf)) {
  try(unlink(combined_pdf), silent = TRUE)
}

combined_pdf_saved <- FALSE

if (capabilities("cairo")) {
  combined_pdf_saved <- tryCatch(
    {
      ggplot2::ggsave(
        filename = combined_pdf,
        plot = p_combined,
        width = 175,
        height = 175,
        units = "mm",
        device = grDevices::cairo_pdf,
        bg = "white"
      )
      file.exists(combined_pdf) &&
        !is.na(file.info(combined_pdf)$size) &&
        file.info(combined_pdf)$size > 0
    },
    error = function(e) {
      message(
        "Cairo PDF export failed for the combined figure: ",
        conditionMessage(e),
        "\\nTrying the standard vector PDF device instead."
      )
      FALSE
    }
  )
}

# Fallback: standard PDF remains vector-based; useDingbats = FALSE
# improves compatibility when opening/editing the PDF in Illustrator.
if (!combined_pdf_saved) {
  if (file.exists(combined_pdf)) {
    try(unlink(combined_pdf), silent = TRUE)
  }
  
  ggplot2::ggsave(
    filename = combined_pdf,
    plot = p_combined,
    width = 175,
    height = 175,
    units = "mm",
    device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "white"
  )
  
  if (!file.exists(combined_pdf) || file.info(combined_pdf)$size <= 0) {
    stop("The combined PDF was not generated successfully. Check whether the file is open or the output folder is writable.")
  }
  
  message("Combined PDF saved using the standard vector PDF device: ", combined_pdf)
} else {
  message("Combined PDF saved using Cairo PDF: ", combined_pdf)
}

# High-resolution PNG export.
ggplot2::ggsave(
  filename = file.path(output_dir, "RT-qPCR_combined.png"),
  plot = p_combined,
  width = 175,
  height = 175,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# -------------------- 12. Export statistical results ------------
write.csv(
  result_summary,
  file.path(output_dir, "RTqPCR_mean_SD.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  WRS2_group_n,
  file.path(output_dir, "RTqPCR_WRS2_group_sample_sizes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  t1way_results,
  file.path(output_dir, "RTqPCR_WRS2_t1way_overall_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  lincon_results,
  file.path(output_dir, "RTqPCR_WRS2_lincon_all_pairwise_rawP_no_adjustment.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Export only comparisons actually displayed on the plots:
# overall t1way P < 0.05 AND pairwise raw lincon P < 0.05.
write.csv(
  lincon_significant,
  file.path(output_dir, "RTqPCR_WRS2_significant_comparisons_displayed.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Save the analysis rules for reproducibility.
analysis_notes <- c(
  "Descriptive statistics: arithmetic mean +/- standard deviation (SD).",
  "Each target gene was analyzed separately.",
  "Overall test: WRS2::t1way() robust one-way ANOVA on 20% trimmed means.",
  "Trimming proportion: tr = 0.20.",
  "Pairwise test: WRS2::lincon() for all 10 pairwise group contrasts per gene.",
  "Pairwise P-value adjustment: method = 'none'; no multiplicity adjustment was applied.",
  "All pairwise lincon contrasts were calculated regardless of the corresponding overall t1way result.",
  "A pairwise comparison is displayed only when BOTH overall t1way P < 0.05 and pairwise raw lincon P < 0.05.",
  "Displayed labels are actual unadjusted pairwise P values: P = 0.xxx or P < 0.001.",
  "No significance letters or asterisks are used.",
  "Bars show arithmetic group means; error bars show mean +/- SD; points show individual observations.",
  "WRS2 inference is based on 20% trimmed means and is distinct from the arithmetic mean +/- SD used for visualization."
)

writeLines(
  analysis_notes,
  file.path(output_dir, "Statistical_analysis_notes.txt"),
  useBytes = TRUE
)

cat("\n========== WRS2 t1way overall results ==========\n")
print(
  as.data.frame(t1way_results),
  row.names = FALSE
)

cat("\n========== WRS2 lincon all pairwise results (raw P) ==========\n")
print(
  as.data.frame(lincon_results),
  row.names = FALSE
)

cat("\n========== Significant comparisons displayed (both P < 0.05) ==========\n")
print(
  as.data.frame(lincon_significant),
  row.names = FALSE
)

cat(
  "\nAnalysis and figures saved to:\n",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE),
  "\n"
)