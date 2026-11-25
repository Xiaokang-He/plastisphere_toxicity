# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# Core105 VFDB heatmaps: log2FC vs Control and log2(mean + 1)
# Statistics: WRS2 t1way + lincon
# lincon P values are unadjusted (method = "none").
# No t1way omnibus significance gate is applied to lincon comparisons.
# All 105 VFDB IDs are retained; input abundances are never changed.
# ============================================================

# -------------------- 1. User settings -----------------------
repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- input_dir
out_dir <- file.path(results_dir, "supp_fig20_vfg_heatmap")

# Star rules based only on unadjusted WRS2 lincon P values:
# ** raw lincon P < 0.01
# *  raw lincon P < 0.05
# No overall t1way screening gate and no P-value adjustment.
alpha_star <- 0.05
alpha_double_star <- 0.01

# WRS2 trimming level (default = 20%).
# Note: with n = 3 per group, floor(0.2 * 3) = 0 observations
# are trimmed; the trimmed means are therefore equal to sample means.
wrs2_tr <- 0.20

fc_limit <- 2
row_order <- "category"   # "category", "input", or "Control_abundance"
make_pseudocount_sensitivity <- TRUE
pseudocount <- NULL       # NULL = 1/2 global minimum positive input abundance
font_family <- if (.Platform$OS.type == "windows" && capabilities("cairo")) {
  "Times New Roman"
} else {
  "serif"
}

group_levels <- c("Control", "MP1", "MP10", "MPP1", "MPP10")
group_prefix <- c(Control="T11", MP1="T6", MP10="T10", MPP1="T4", MPP10="T3")
core_file <- file.path(base_dir, "vfdb_core_VFDB_ID_profile.xls")
annotation_file <- file.path(base_dir, "VFDB毒力因子丰度表_VFDB_Origin.csv")

if (!requireNamespace("ggplot2", quietly=TRUE)) {
  install.packages("ggplot2", repos="https://cloud.r-project.org")
}
if (!requireNamespace("WRS2", quietly=TRUE)) {
  install.packages("WRS2", repos="https://cloud.r-project.org")
}
library(ggplot2)
library(WRS2)

# The current WRS2 API exposes method in lincon(); method="none"
# is required here to prevent its default Hochberg P-value adjustment.
if (!"method" %in% names(formals(WRS2::lincon))) {
  stop(
    "The installed WRS2 version does not expose lincon(method=...). ",
    "Please update WRS2 before running this script."
  )
}

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(core_file)) stop("Core profile file not found: ", core_file)
stopifnot(is.finite(fc_limit), fc_limit > 0)
stopifnot(row_order %in% c("category", "input", "Control_abundance"))

# -------------------- 2. Read and validate -------------------
# The .xls profile is tab-delimited text, not a binary Excel workbook.
dat <- read.delim(
  core_file, header=TRUE, sep="\t", check.names=FALSE,
  quote="", comment.char="", fileEncoding="UTF-8-BOM",
  stringsAsFactors=FALSE
)

if (!"VFDB_ID" %in% names(dat)) stop("Missing VFDB_ID column.")
if (nrow(dat) != 105L) stop("Expected 105 core IDs; found ", nrow(dat))
if (anyDuplicated(dat$VFDB_ID)) stop("Duplicated VFDB_IDs detected.")

sample_cols <- unlist(lapply(group_levels, function(g) {
  paste0(unname(group_prefix[g]), "_M_", 1:3)
}), use.names=FALSE)

missing_cols <- setdiff(sample_cols, names(dat))
if (length(missing_cols)) {
  stop("Missing samples: ", paste(missing_cols, collapse=", "))
}

X <- as.matrix(dat[, sample_cols, drop=FALSE])
storage.mode(X) <- "double"
rownames(X) <- dat$VFDB_ID

if (any(!is.finite(X)) || any(X < 0)) {
  stop("Missing, invalid, or negative abundance values detected.")
}

if ("Total" %in% names(dat)) {
  if (!isTRUE(all.equal(
    unname(rowSums(X)), as.numeric(dat$Total), tolerance=1e-8
  ))) {
    stop("Total column does not match the 15 sample columns.")
  }
}

gene_ids <- rownames(X)
n_gene <- nrow(X)
sample_group <- factor(rep(group_levels, each=3), levels=group_levels)

group_mean <- sapply(group_levels, function(g) {
  rowMeans(X[, sample_group==g, drop=FALSE])
})
group_sd <- sapply(group_levels, function(g) {
  apply(X[, sample_group==g, drop=FALSE], 1, sd)
})
group_n_detected <- sapply(group_levels, function(g) {
  rowSums(X[, sample_group==g, drop=FALSE] > 0)
})

rownames(group_mean) <- rownames(group_sd) <-
  rownames(group_n_detected) <- gene_ids
colnames(group_mean) <- colnames(group_sd) <-
  colnames(group_n_detected) <- group_levels

# -------------------- 3. Optional annotations ----------------
gene_name <- rep("", n_gene)
category <- rep("Unannotated", n_gene)

if (file.exists(annotation_file)) {
  ann <- read.csv(
    annotation_file, check.names=FALSE, stringsAsFactors=FALSE,
    fileEncoding="UTF-8-BOM"
  )
  id_col <- intersect(c("VFDB ID", "VFDB_ID"), names(ann))
  if (length(id_col)) {
    if (anyDuplicated(ann[[id_col[1]]])) {
      stop("Duplicated VFDB IDs in annotation file.")
    }
    mi <- match(gene_ids, ann[[id_col[1]]])
    if ("Related genes" %in% names(ann)) {
      gene_name <- as.character(ann[["Related genes"]][mi])
    }
    if ("VFcategory" %in% names(ann)) {
      category <- as.character(ann[["VFcategory"]][mi])
    }
    gene_name[is.na(gene_name)] <- ""
    category[is.na(category) | category==""] <- "Unannotated"
  }
}

# -------------------- 4. WRS2 t1way + lincon ------------------
# For each VFDB_ID:
# 1) t1way tests the overall difference among the five groups.
# 2) lincon computes pairwise comparisons using trimmed means.
# 3) Only Control vs MP1/MP10/MPP1/MPP10 are retained below.
# 4) method = "none" explicitly disables lincon P-value adjustment.
# 5) The t1way result is exported but is NOT used as a screening gate.
#
# Important:
# WRS2::lincon() computes all pairwise contrasts internally. The
# comparisons retained for the heatmaps are the four contrasts with
# Control. With method="none", their P values are not multiplicity-
# adjusted. The returned confidence intervals are simultaneous CIs.
#
# t1way() may be non-estimable when any group's Winsorized variance is
# zero. Such failures are recorded, and lincon() is still attempted
# independently because no omnibus gate is applied.

t1way_list <- vector("list", n_gene)
lincon_list <- vector("list", n_gene)

for (i in seq_len(n_gene)) {
  values <- as.numeric(X[i, ])

  analysis_dat <- data.frame(
    Abundance=values,
    Group=factor(sample_group, levels=group_levels)
  )

  # Overall robust one-way test. Its P value is descriptive/exported
  # only; it does not determine whether lincon comparisons are shown.
  t1 <- tryCatch(
    WRS2::t1way(
      Abundance ~ Group,
      data=analysis_dat,
      tr=wrs2_tr,
      nboot=100
    ),
    error=function(e) e
  )

  if (inherits(t1, "error")) {
    t1way_list[[i]] <- data.frame(
      VFDB_ID=gene_ids[i],
      T1way_F=NA_real_,
      T1way_df1=NA_real_,
      T1way_df2=NA_real_,
      T1way_P_raw=NA_real_,
      T1way_effect_size=NA_real_,
      T1way_error=conditionMessage(t1),
      stringsAsFactors=FALSE
    )
  } else {
    t1way_list[[i]] <- data.frame(
      VFDB_ID=gene_ids[i],
      T1way_F=as.numeric(t1$test),
      T1way_df1=as.numeric(t1$df1),
      T1way_df2=as.numeric(t1$df2),
      T1way_P_raw=as.numeric(t1$p.value),
      T1way_effect_size=if (!is.null(t1$effsize)) {
        as.numeric(t1$effsize)
      } else {
        NA_real_
      },
      T1way_error="",
      stringsAsFactors=FALSE
    )
  }

  # Pairwise robust contrasts. method="none" is essential:
  # the current WRS2 default is method="hochberg".
  lc <- tryCatch(
    WRS2::lincon(
      Abundance ~ Group,
      data=analysis_dat,
      tr=wrs2_tr,
      alpha=0.05,
      method="none"
    ),
    error=function(e) e
  )

  if (inherits(lc, "error")) {
    lincon_list[[i]] <- data.frame(
      VFDB_ID=gene_ids[i],
      Group=group_levels[-1],
      Contrast=paste("Control vs", group_levels[-1]),
      lincon_psihat_Control_minus_treatment=NA_real_,
      lincon_CI_lower_simultaneous=NA_real_,
      lincon_CI_upper_simultaneous=NA_real_,
      lincon_P_raw=NA_real_,
      lincon_error=conditionMessage(lc),
      stringsAsFactors=FALSE
    )
  } else {
    comp <- as.matrix(lc$comp)

    if (ncol(comp) < 6L) {
      stop(
        "Unexpected WRS2::lincon output structure for ",
        gene_ids[i],
        ". Expected six columns in lc$comp."
      )
    }

    control_idx <- match("Control", lc$fnames)
    if (is.na(control_idx)) {
      stop("Control group was not found in WRS2::lincon output.")
    }

    group1_idx <- as.integer(comp[, 1])
    group2_idx <- as.integer(comp[, 2])

    one_gene_lincon <- lapply(group_levels[-1], function(g) {
      treatment_idx <- match(g, lc$fnames)
      hit <- which(
        group1_idx == control_idx &
          group2_idx == treatment_idx
      )

      if (length(hit) != 1L) {
        return(data.frame(
          VFDB_ID=gene_ids[i],
          Group=g,
          Contrast=paste("Control vs", g),
          lincon_psihat_Control_minus_treatment=NA_real_,
          lincon_CI_lower_simultaneous=NA_real_,
          lincon_CI_upper_simultaneous=NA_real_,
          lincon_P_raw=NA_real_,
          lincon_error="Expected Control-treatment contrast not found exactly once.",
          stringsAsFactors=FALSE
        ))
      }

      z <- comp[hit, , drop=FALSE]

      data.frame(
        VFDB_ID=gene_ids[i],
        Group=g,
        Contrast=paste("Control vs", g),
        lincon_psihat_Control_minus_treatment=as.numeric(z[1, 3]),
        lincon_CI_lower_simultaneous=as.numeric(z[1, 4]),
        lincon_CI_upper_simultaneous=as.numeric(z[1, 5]),
        lincon_P_raw=as.numeric(z[1, 6]),
        lincon_error="",
        stringsAsFactors=FALSE
      )
    })

    lincon_list[[i]] <- do.call(rbind, one_gene_lincon)
  }
}

t1way_results <- do.call(rbind, t1way_list)
tests <- do.call(rbind, lincon_list)

# Add the overall t1way P value for reference only.
# It is deliberately NOT used to filter or gate lincon significance.
tests$T1way_P_raw <- t1way_results$T1way_P_raw[
  match(tests$VFDB_ID, t1way_results$VFDB_ID)
]

# Significance stars use raw, unadjusted lincon P values only.
tests$stars <- ""
tests$stars[
  is.finite(tests$lincon_P_raw) &
    tests$lincon_P_raw < alpha_star
] <- "*"
tests$stars[
  is.finite(tests$lincon_P_raw) &
    tests$lincon_P_raw < alpha_double_star
] <- "**"

tests$mark_significant <- tests$stars != ""

# -------------------- 5. Effect sizes and heatmap data --------
# Main matrix: true log2(mean_treatment / mean_Control), without
# pseudocount. Zero denominators/numerators are explicitly flagged.
plot_dat <- expand.grid(
  VFDB_ID=gene_ids,
  Group=group_levels,
  KEEP.OUT.ATTRS=FALSE,
  stringsAsFactors=FALSE
)

plot_dat$Mean <- as.vector(group_mean)
plot_dat$SD <- as.vector(group_sd)
plot_dat$n_detected <- as.vector(group_n_detected)
plot_dat$Control_mean <- group_mean[
  match(plot_dat$VFDB_ID, gene_ids), "Control"
]

plot_dat$status <- "finite"
plot_dat$status[plot_dat$Mean==0 & plot_dat$Control_mean==0] <-
  "both_not_detected"
plot_dat$status[plot_dat$Mean>0 & plot_dat$Control_mean==0] <-
  "Control_zero_treatment_detected"
plot_dat$status[plot_dat$Mean==0 & plot_dat$Control_mean>0] <-
  "treatment_not_detected"
plot_dat$status[
  plot_dat$Group=="Control" & plot_dat$Control_mean>0
] <- "reference"

plot_dat$log2FC <- NA_real_
finite_idx <- plot_dat$Mean>0 & plot_dat$Control_mean>0
plot_dat$log2FC[finite_idx] <- log2(
  plot_dat$Mean[finite_idx] / plot_dat$Control_mean[finite_idx]
)
plot_dat$log2_mean_plus1 <- log2(plot_dat$Mean+1)

if (is.null(pseudocount)) {
  positive_values <- X[X>0]
  if (!length(positive_values)) stop("All input abundances are zero.")
  pseudocount <- min(positive_values)/2
}
stopifnot(
  length(pseudocount)==1L,
  is.finite(pseudocount),
  pseudocount>0
)

plot_dat$log2FC_shared_pseudocount <- log2(
  (plot_dat$Mean+pseudocount) /
    (plot_dat$Control_mean+pseudocount)
)
plot_dat$pseudocount <- pseudocount
plot_dat$zero_affected <- plot_dat$Mean==0 |
  plot_dat$Control_mean==0

# Map each treatment-vs-Control lincon result to its heatmap cell.
test_key <- paste(tests$VFDB_ID, tests$Group, sep="::")
plot_key <- paste(plot_dat$VFDB_ID, plot_dat$Group, sep="::")
pi <- match(plot_key, test_key)

plot_dat$lincon_P_raw <- tests$lincon_P_raw[pi]
plot_dat$T1way_P_raw <- tests$T1way_P_raw[pi]
plot_dat$stars <- tests$stars[pi]
plot_dat$stars[is.na(plot_dat$stars)] <- ""

# The Control column is a reference and is not assigned a pairwise star.
plot_dat$stars[plot_dat$Group=="Control"] <- ""

# -------------------- 6. Stable row order and labels ----------
ord <- switch(
  row_order,
  category=order(category, -group_mean[, "Control"], gene_ids),
  input=seq_len(n_gene),
  Control_abundance=order(-group_mean[, "Control"], gene_ids)
)

# Row labels: show gene names only (e.g., tufA, groEL, htpB).
# Do not display VFDB IDs or Control abundance values.
# If a gene name is unavailable, leave that row label blank.
label <- gene_name
label[is.na(label)] <- ""
names(label) <- gene_ids

plot_dat$Row <- factor(
  plot_dat$VFDB_ID,
  levels=rev(gene_ids[ord])
)
plot_dat$Group <- factor(plot_dat$Group, levels=group_levels)
plot_dat$Gene_name <- gene_name[match(plot_dat$VFDB_ID, gene_ids)]
plot_dat$Category <- category[match(plot_dat$VFDB_ID, gene_ids)]

# -------------------- 7. Heatmap functions --------------------
star_note <- paste0(
  "Stars are based on raw two-sided WRS2 lincon P values for Control vs each treatment; ",
  "method='none' (no P-value adjustment): ** P < 0.01; * P < 0.05. ",
  "No t1way omnibus significance gate was applied."
)

base_theme <- theme_minimal(
  base_family=font_family,
  base_size=11
) +
  theme(
    panel.grid=element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_text(size=6.6, color="black"),
    axis.text.x=element_text(size=11, color="black"),
    plot.title=element_text(face="bold", size=14),
    plot.subtitle=element_text(size=9),
    plot.caption=element_text(size=8, hjust=0),
    legend.position="top",
    plot.margin=margin(10, 15, 10, 10)
  )

make_fc_plot <- function(df, value_column, sensitivity=FALSE) {
  df$FC_display <- pmax(
    -fc_limit,
    pmin(fc_limit, df[[value_column]])
  )

  df$cell_label <- ""

  if (!sensitivity) {
    df$cell_label[df$status=="both_not_detected"] <- "ND"
    df$cell_label[
      df$status=="Control_zero_treatment_detected"
    ] <- "+"
    df$cell_label[df$status=="treatment_not_detected"] <- "-"
  } else {
    df$cell_label[df$zero_affected] <- "z"
  }

  good <- is.finite(df[[value_column]])
  df$cell_label[
    good & df[[value_column]]>fc_limit
  ] <- paste0(
    df$cell_label[good & df[[value_column]]>fc_limit], ">"
  )
  df$cell_label[
    good & df[[value_column]]< -fc_limit
  ] <- paste0(
    df$cell_label[good & df[[value_column]]< -fc_limit], "<"
  )

  df$cell_label <- paste0(df$cell_label, df$stars)
  df$text_color <- ifelse(
    is.finite(df$FC_display) & abs(df$FC_display)>1.5,
    "white", "black"
  )

  zero_note <- if (!sensitivity) {
    paste(
      "Gray: non-finite ratio.",
      "+: Control zero, treatment detected;",
      "-: treatment zero; ND: both zero."
    )
  } else {
    paste0(
      "Shared pseudocount = ",
      format(pseudocount, digits=6),
      "; z: a zero input mean affects this ratio."
    )
  }

  ggplot(df, aes(x=Group, y=Row)) +
    geom_tile(
      aes(fill=FC_display),
      color="white",
      linewidth=.18
    ) +
    geom_text(
      aes(label=cell_label, color=text_color),
      size=2.4,
      family=font_family
    ) +
    scale_color_identity(guide="none") +
    scale_y_discrete(labels=label, drop=FALSE) +
    scale_x_discrete(
      drop=FALSE,
      position="top",
      expand=expansion(add=.05)
    ) +
    scale_fill_gradient2(
      low="#2C6E9B",
      mid="white",
      high="#B84640",
      midpoint=0,
      limits=c(-fc_limit, fc_limit),
      breaks=seq(-fc_limit, fc_limit, length.out=5),
      na.value="#D5D5D5",
      name=if (sensitivity) {
        "Regularized log2 ratio"
      } else {
        "log2FC vs Control"
      }
    ) +
    labs(
      title=if (sensitivity) {
        "Core105: shared-pseudocount sensitivity"
      } else {
        "Core105: log2 fold change relative to Control"
      },
      subtitle="Group means; no row Z standardization.",
      caption=paste(
        zero_note,
        paste0(
          "Fixed scale +/-", fc_limit,
          "; < or > marks values beyond the color limits.",
          " Full values are exported."
        ),
        star_note,
        sep="\n"
      )
    ) +
    base_theme
}

# Heatmap 1: fold change relative to Control.
p_fc <- make_fc_plot(plot_dat, "log2FC")

# Heatmap 2: abundance display, log2(mean + 1).
# Stars are still treatment-vs-Control tests, not abundance transformations.
abund_dat <- plot_dat
abund_dat$text_color <- ifelse(
  abund_dat$log2_mean_plus1 >
    .7*max(abund_dat$log2_mean_plus1),
  "white", "black"
)

p_abundance <- ggplot(abund_dat, aes(x=Group, y=Row)) +
  geom_tile(
    aes(fill=log2_mean_plus1),
    color="white",
    linewidth=.18
  ) +
  geom_text(
    aes(label=stars, color=text_color),
    size=2.4,
    family=font_family
  ) +
  scale_color_identity(guide="none") +
  scale_y_discrete(labels=label, drop=FALSE) +
  scale_x_discrete(
    drop=FALSE,
    position="top",
    expand=expansion(add=.05)
  ) +
  scale_fill_gradient(
    low="white",
    high="#356D96",
    name="log2(mean + 1)"
  ) +
  labs(
    title="Core105: abundance on one common scale",
    subtitle="Group means; log2(mean + 1); no row Z standardization.",
    caption=paste(
      "Stars indicate significant treatment-vs-Control differences in the raw abundance data.",
      star_note,
      sep="\n"
    )
  ) +
  base_theme

# -------------------- 8. Save heatmaps -----------------------
fig_width <- 12
fig_height <- max(14, .18*n_gene+3.5)

save_plot <- function(p, stem) {
  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    p,
    width=fig_width,
    height=fig_height,
    units="in",
    dpi=300,
    bg="white",
    limitsize=FALSE
  )

  if (capabilities("cairo")) {
    ggsave(
      file.path(out_dir, paste0(stem, ".pdf")),
      p,
      device=grDevices::cairo_pdf,
      width=fig_width,
      height=fig_height,
      units="in",
      family=font_family,
      limitsize=FALSE
    )
  } else {
    p_pdf <- p + theme(text=element_text(family="Times"))
    ggsave(
      file.path(out_dir, paste0(stem, ".pdf")),
      p_pdf,
      device="pdf",
      family="Times",
      width=fig_width,
      height=fig_height,
      units="in",
      limitsize=FALSE
    )
  }
}

save_plot(p_fc, "01_Core105_true_log2FC_vs_Control_WRS2_lincon_rawP")
save_plot(p_abundance, "02_Core105_log2_mean_plus1_WRS2_lincon_rawP")

if (make_pseudocount_sensitivity) {
  p_pc <- make_fc_plot(
    plot_dat,
    "log2FC_shared_pseudocount",
    sensitivity=TRUE
  )
  save_plot(
    p_pc,
    "03_Core105_shared_pseudocount_sensitivity_WRS2_lincon_rawP"
  )
}

# -------------------- 9. Export numerical results -------------
write_csv <- function(x, name) {
  write.csv(
    x,
    file.path(out_dir, name),
    row.names=FALSE,
    na="NA",
    fileEncoding="UTF-8"
  )
}

write_csv(
  data.frame(VFDB_ID=gene_ids, X, check.names=FALSE),
  "Core105_original_sample_values.csv"
)
write_csv(
  data.frame(VFDB_ID=gene_ids, group_mean, check.names=FALSE),
  "Core105_group_mean.csv"
)
write_csv(t1way_results, "Core105_WRS2_t1way_rawP_and_errors.csv")
write_csv(tests, "Core105_WRS2_lincon_Control_rawP_noAdjustment.csv")
write_csv(plot_dat, "Core105_heatmap_values_and_WRS2_lincon_rawP_stars.csv")

summary_list <- lapply(group_levels[-1], function(g) {
  z <- plot_dat[plot_dat$Group==g, ]
  f <- is.finite(z$log2FC)
  values <- z$log2FC[f]

  data.frame(
    Group=g,
    all_IDs=n_gene,
    finite_log2FC_IDs=sum(f),
    median_abs_log2FC=if (length(values)) median(abs(values)) else NA_real_,
    n_within_1_5_fold=sum(abs(values)<=log2(1.5)),
    fraction_within_1_5_fold=mean(abs(values)<=log2(1.5)),
    n_within_2_fold=sum(abs(values)<=1),
    fraction_within_2_fold=mean(abs(values)<=1),
    Control_zero_treatment_detected=sum(
      z$status=="Control_zero_treatment_detected"
    ),
    treatment_not_detected=sum(
      z$status=="treatment_not_detected"
    ),
    both_not_detected=sum(z$status=="both_not_detected"),
    n_single_star=sum(z$stars=="*"),
    n_double_star=sum(z$stars=="**"),
    stringsAsFactors=FALSE
  )
})

effect_summary <- do.call(rbind, summary_list)
write_csv(effect_summary, "Core105_effect_size_summary.csv")

notes <- c(
  "All 105 Core VFDB IDs are retained.",
  "Main matrix: log2(mean_treatment / mean_Control), without pseudocount or row Z standardization.",
  "The abundance heatmap displays log2(mean + 1); this is not a fold change.",
  "Means are calculated from the three biological replicates in each group.",
  "For each VFDB_ID, WRS2 t1way tests the overall difference among the five groups.",
  "WRS2 lincon computes pairwise trimmed-mean contrasts; only Control vs MP1, MP10, MPP1, and MPP10 are retained for reporting and heatmap annotation.",
  "WRS2 lincon was called with method='none'; the reported lincon P values are not adjusted.",
  "No t1way omnibus significance gate is applied to the lincon comparisons.",
  "No multiple-testing correction is applied across the four Control contrasts or across the 105 VFDB IDs.",
  paste0("WRS2 trimming level: tr = ", wrs2_tr, "."),
  "With n=3 per group and tr=0.20, floor(n*tr)=0; therefore no observations are trimmed and the trimmed means equal the sample means.",
  "The confidence intervals returned by lincon are simultaneous confidence intervals; the P values used for stars are unadjusted because method='none'.",
  "If t1way cannot be computed because a group's Winsorized variance is zero, the error is recorded in the t1way results file; lincon is still attempted independently.",
  paste0("Star rule: ** raw lincon P < ", alpha_double_star,
         "; * raw lincon P < ", alpha_star, "."),
  "Control cells are the reference and do not receive pairwise stars.",
  "Zero values are retained; non-finite fold changes are flagged rather than replaced.",
  paste0("Pseudocount sensitivity uses one shared constant: ",
         format(pseudocount, digits=12)),
  "Because 105 IDs are tested across four treatment-vs-Control contrasts without P-value correction, star annotations should be interpreted as exploratory.",
  "Input abundance units and cross-sample normalization are not documented in the supplied files.",
  "n=3 per group; estimates and inferential results should be interpreted cautiously given the small sample size."
)

writeLines(
  notes,
  file.path(out_dir, "Analysis_notes.txt"),
  useBytes=TRUE
)

capture.output(
  sessionInfo(),
  file=file.path(out_dir, "R_session_info.txt")
)

print(effect_summary)
message(
  "Saved heatmaps and numerical results to: ",
  normalizePath(out_dir, winslash="/", mustWork=FALSE)
)