# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# PARAFAC C1 FMax 柱状图
#
# 柱高：Mean；误差棒：Mean ± SD
# 散点：每个独立重复的原始值，仅沿水平方向轻微错开
#
# 统计：
#   WRS2::t1way() + WRS2::lincon()
#   t1way：四组总体异方差稳健单因素截尾均值检验
#   lincon：四组全部 6 次两两比较
#   tr = 0.20：20% 截尾均值
#   method = "none"：两两比较不进行多重比较校正
#   标注条件：总体 t1way P < 0.05 且 lincon 原始 P < 0.05
#   标注内容：实际 P 值；保留三位小数；P < 0.001 显示为 <0.001
#
# 样本映射：
#   N_S3 = Day7；N_S4 = Day14
#   N_S2 = Day30；N_S1 = Day60
#
# 输出：PNG/PDF，108 × 115 mm
# 同时导出 t1way、全部 lincon 两两比较及描述性统计 CSV。
# ============================================================


# ============================================================
# 0. 清理环境
# ============================================================

rm(list = ls())


# ============================================================
# 1. 加载依赖
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
  "ggplot2",
  "ggpattern",
  "WRS2",
  "scales"
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
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggpattern)
  library(WRS2)
  library(scales)
  library(grid)
})


# ============================================================
# 2. 设置路径并读取数据
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- input_dir
output_dir <- file.path(results_dir, "supp_fig3_parafac_fmax")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

file_path <- file.path(
  input_dir,
  "3-component_result - 副本.xlsx"
)

if (!dir.exists(base_dir)) {
  stop("数据文件夹不存在：\n", base_dir)
}

if (!file.exists(file_path)) {
  stop("找不到原始 Excel 文件：\n", file_path)
}

raw_data <- readxl::read_excel(
  file_path,
  sheet = "FMax"
)

if (!all(c("EEM", "FMax1") %in% names(raw_data))) {
  stop("FMax 工作表必须包含 EEM 和 FMax1 两列，请检查列名。")
}

if (!is.numeric(raw_data$FMax1)) {
  stop("FMax1 列必须为数值型，请检查原始数据。")
}

cat("\n========== 原始数据 ==========\n")
print(raw_data)


# ============================================================
# 3. 样本映射与数据整理
# ============================================================

groups <- c(
  "Day7",
  "Day14",
  "Day30",
  "Day60"
)

plot_data <- raw_data %>%
  mutate(
    EEM = as.character(EEM),
    Time = case_when(
      str_detect(EEM, fixed("N_S3")) ~ "Day7",
      str_detect(EEM, fixed("N_S4")) ~ "Day14",
      str_detect(EEM, fixed("N_S2")) ~ "Day30",
      str_detect(EEM, fixed("N_S1")) ~ "Day60",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Time),
    !is.na(FMax1),
    is.finite(FMax1)
  ) %>%
  mutate(
    Time = factor(Time, levels = groups)
  )

if (nrow(plot_data) == 0) {
  stop(
    "没有识别到有效样本，请检查 FMax 工作表 EEM 列中的样本编号。"
  )
}

group_n <- table(
  factor(plot_data$Time, levels = groups)
)

if (any(group_n < 2)) {
  stop(
    "每个时间组至少需要 2 个有效独立观测值。各组当前样本量：\n",
    paste(
      names(group_n),
      as.integer(group_n),
      sep = " = ",
      collapse = "; "
    )
  )
}

cat("\n========== 各组样本量 ==========\n")
print(group_n)

cat("\n========== 整理后的分析数据 ==========\n")
print(plot_data)


# ============================================================
# 4. 描述性统计：Mean ± SD
# ============================================================

summary_data <- plot_data %>%
  group_by(Time) %>%
  summarise(
    n = sum(!is.na(FMax1)),
    Mean = mean(FMax1, na.rm = TRUE),
    SD = sd(FMax1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Time)

cat("\n========== Mean ± SD ==========\n")
print(summary_data, n = Inf)


# ============================================================
# 5. WRS2 差异检验
#
# 与前述 EPS 定量分析保持一致：
#   总体检验：WRS2::t1way()
#   两两比较：WRS2::lincon()
#   截尾比例：tr = 0.20
#   两两比较：method = "none"，不进行多重比较校正
# ============================================================

TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05


# ------------------------------------------------------------
# 5.1 总体检验：t1way
# ------------------------------------------------------------

t1way_fit <- tryCatch(
  WRS2::t1way(
    FMax1 ~ Time,
    data = plot_data,
    tr = TRIM_LEVEL
  ),
  error = function(e) {
    stop(
      "WRS2::t1way() 运行失败：\n",
      conditionMessage(e),
      "\n请检查各组样本量、组内重复值及 Winsorized variance。"
    )
  }
)

t1way_result <- data.frame(
  Test = "WRS2 t1way",
  Trim_level = TRIM_LEVEL,
  Statistic_F = as.numeric(t1way_fit$test),
  df1 = as.numeric(t1way_fit$df1),
  df2 = as.numeric(t1way_fit$df2),
  p_t1way = as.numeric(t1way_fit$p.value),
  Effect_size = as.numeric(t1way_fit$effsize),
  stringsAsFactors = FALSE
)

if (length(t1way_result$p_t1way) != 1 ||
    !is.finite(t1way_result$p_t1way) ||
    t1way_result$p_t1way < 0 ||
    t1way_result$p_t1way > 1) {
  stop("t1way 总体检验返回了无效 P 值，请检查数据及 WRS2 输出。")
}

cat("\n========== WRS2 t1way 总体检验 ==========\n")
print(t1way_result)


# ------------------------------------------------------------
# 5.2 两两比较：lincon
#
# 对四组进行全部 6 次两两比较。
# method = "none" 表示不进行多重比较校正。
#
# WRS2::lincon() 的 comp 矩阵：
#   第1列、第2列：比较组别在因子水平中的序号
#   第3列：psihat（截尾均值差）
#   第4列：ci.lower
#   第5列：ci.upper
#   第6列：p.value（未校正 P 值）
# ------------------------------------------------------------

lincon_fit <- tryCatch(
  WRS2::lincon(
    FMax1 ~ Time,
    data = plot_data,
    tr = TRIM_LEVEL,
    method = "none"
  ),
  error = function(e) {
    stop(
      "WRS2::lincon() 运行失败：\n",
      conditionMessage(e),
      "\n请检查各组样本量、组内重复值及 Winsorized variance。"
    )
  }
)

comp <- as.data.frame(lincon_fit$comp)

if (ncol(comp) < 6) {
  stop(
    "lincon 结果列数不足 6 列，无法正确提取组别、置信区间和 P 值。"
  )
}

group_index_1 <- as.integer(comp[[1]])
group_index_2 <- as.integer(comp[[2]])

group_names <- lincon_fit$fnames

if (is.null(group_names) || length(group_names) != length(groups)) {
  # 兼容部分 WRS2 版本：按输入因子水平顺序恢复组名
  group_names <- levels(droplevels(plot_data$Time))
}

if (anyNA(group_index_1) ||
    anyNA(group_index_2) ||
    any(group_index_1 < 1 | group_index_1 > length(group_names)) ||
    any(group_index_2 < 1 | group_index_2 > length(group_names))) {
  stop(
    "lincon 结果无法映射回 Day7、Day14、Day30、Day60，",
    "请检查 WRS2 输出结构。"
  )
}

lincon_result <- data.frame(
  group1 = as.character(group_names[group_index_1]),
  group2 = as.character(group_names[group_index_2]),
  psihat = as.numeric(comp[[3]]),
  ci_lower = as.numeric(comp[[4]]),
  ci_upper = as.numeric(comp[[5]]),
  p_raw = as.numeric(comp[[6]]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    x1 = match(group1, groups),
    x2 = match(group2, groups)
  ) %>%
  arrange(
    x2 - x1,
    x1
  )

if (nrow(lincon_result) != choose(length(groups), 2)) {
  stop(
    "lincon 未返回预期的全部 6 次两两比较；实际返回 ",
    nrow(lincon_result),
    " 次。"
  )
}

if (anyNA(lincon_result$x1) ||
    anyNA(lincon_result$x2) ||
    any(!is.finite(lincon_result$p_raw)) ||
    any(lincon_result$p_raw < 0 | lincon_result$p_raw > 1)) {
  stop("lincon 结果存在无法识别的组别或无效 P 值。")
}

# 将总体 t1way P 值并入每一条两两比较结果
lincon_result <- lincon_result %>%
  mutate(
    p_t1way = t1way_result$p_t1way[1],
    Test = "WRS2 lincon",
    Trim_level = TRIM_LEVEL,
    P_adjustment = "none"
  ) %>%
  select(
    Test,
    Trim_level,
    group1,
    group2,
    psihat,
    ci_lower,
    ci_upper,
    p_raw,
    p_t1way,
    P_adjustment,
    x1,
    x2
  )

cat("\n========== WRS2 lincon 全部 6 次两两比较 ==========\n")
print(as.data.frame(lincon_result), row.names = FALSE)


# ------------------------------------------------------------
# 5.3 P 值格式
#
# P >= 0.001：保留三位小数，例如 P = 0.023
# P < 0.001：统一标注为 P < 0.001，不使用科学计数法
# P 字符斜体；数值保留三位小数。
# ------------------------------------------------------------

format_p_label <- function(p) {

  vapply(
    p,
    function(z) {

      if (is.na(z)) {
        return("italic(P) == 'NA'")
      }

      if (z < 0.001) {
        return("italic(P) < 0.001")
      }

      sprintf(
        "italic(P) == '%.3f'",
        z
      )
    },
    character(1)
  )
}


# ------------------------------------------------------------
# 5.4 显著性筛选
#
# 与前述 EPS 分析保持一致：
#   总体 t1way P < 0.05
#   且对应 lincon 原始 P < 0.05
#
# 不使用校正后的 P 值，不进行多重比较校正。
# ------------------------------------------------------------

ann <- lincon_result %>%
  filter(
    p_t1way < ALPHA_LEVEL,
    p_raw < ALPHA_LEVEL
  ) %>%
  mutate(
    p_label = format_p_label(p_raw)
  ) %>%
  arrange(
    x2 - x1,
    x1
  )

cat("\n========== 满足标注条件的两两比较 ==========\n")
print(as.data.frame(ann), row.names = FALSE)


# ============================================================
# 6. 配色和柱形纹理
# ============================================================

colors <- c(
  Day7 = "#92CDE3",
  Day14 = "#95C3AA",
  Day30 = "#82A2CD",
  Day60 = "#CAB9D5"
)

patterns <- c(
  Day7 = "stripe",
  Day14 = "stripe",
  Day30 = "crosshatch",
  Day60 = "centered_grid"
)

angles <- c(
  Day7 = 45,
  Day14 = 135,
  Day30 = 45,
  Day60 = 0
)


# ============================================================
# 7. Day60 居中网格纹理
# ============================================================

centered_grid <- function(
    params,
    boundary_df,
    aspect_ratio,
    legend = FALSE
) {

  xlim <- grid::convertX(
    grid::unit(range(boundary_df$x), "npc"),
    "mm",
    valueOnly = TRUE
  )

  ylim <- grid::convertY(
    grid::unit(range(boundary_df$y), "npc"),
    "mm",
    valueOnly = TRUE
  )

  step <- grid::convertWidth(
    grid::unit(
      params$pattern_spacing / 2,
      params$pattern_units
    ),
    "mm",
    valueOnly = TRUE
  )

  centered_positions <- function(limits) {

    n <- floor(diff(limits) / (2 * step))

    if (!is.finite(n) || n < 0) {
      return(mean(limits))
    }

    mean(limits) + seq.int(-n, n) * step
  }

  line_width <- params$pattern_linewidth

  if (is.null(line_width)) {
    line_width <- params$pattern_size
  }

  gp <- grid::gpar(
    col = params$pattern_colour,
    lwd = line_width * 72.27 / 25.4
  )

  grid::grobTree(
    grid::segmentsGrob(
      x0 = centered_positions(xlim),
      x1 = centered_positions(xlim),
      y0 = ylim[1],
      y1 = ylim[2],
      default.units = "mm",
      gp = gp
    ),
    grid::segmentsGrob(
      x0 = xlim[1],
      x1 = xlim[2],
      y0 = centered_positions(ylim),
      y1 = centered_positions(ylim),
      default.units = "mm",
      gp = gp
    )
  )
}

pattern_funcs <- getOption(
  "ggpattern_geometry_funcs",
  list()
)

pattern_funcs$centered_grid <- centered_grid

options(
  ggpattern_geometry_funcs = pattern_funcs
)


# ============================================================
# 8. 实际 P 值标注位置
# ============================================================

base_y <- max(
  c(
    summary_data$Mean + summary_data$SD,
    plot_data$FMax1
  ),
  na.rm = TRUE
)

ann <- ann %>%
  mutate(
    bracket_y = base_y * (
      1.06 + (row_number() - 1) * 0.100
    ),
    tip_y = bracket_y - base_y * 0.018,
    text_y = bracket_y + base_y * 0.025,
    text_x = (x1 + x2) / 2
  )

upper <- max(
  c(
    base_y * 1.108,
    ann$text_y + base_y * 0.070
  ),
  na.rm = TRUE
)


# ============================================================
# 9. 绘图：Mean ± SD + 原始散点 + WRS2 实际 P 值
#
# 保留原始图形参数，不因更换统计方法改变绘图样式。
# ============================================================

p <- ggplot(
  summary_data,
  aes(
    x = Time,
    y = Mean,
    fill = Time
  )
) +

  # 均值柱与纹理
  geom_col_pattern(
    aes(
      pattern = Time,
      pattern_angle = Time
    ),
    pattern_spacing = 0.08,
    pattern_units = "snpc",
    pattern_density = 0.50,
    pattern_fill = NA,
    pattern_fill2 = NA,
    pattern_colour = "black",
    pattern_size = 0.22,
    pattern_grid = "square",
    width = 0.78,
    colour = "black",
    linewidth = 0.55
  ) +

  # Mean ± SD
  geom_errorbar(
    aes(
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.17,
    linewidth = 0.55
  ) +

  # 原始散点：白色填充、黑色边框
  # height = 0：保持原始纵坐标
  geom_point(
    data = plot_data,
    aes(
      x = Time,
      y = FMax1
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 2026
    ),
    shape = 21,
    size = 2.4,
    stroke = 0.55,
    colour = "black",
    fill = "white"
  ) +

  # 比较括号横线
  geom_segment(
    data = ann,
    aes(
      x = x1,
      xend = x2,
      y = bracket_y,
      yend = bracket_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +

  # 比较括号左侧短线
  geom_segment(
    data = ann,
    aes(
      x = x1,
      xend = x1,
      y = tip_y,
      yend = bracket_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +

  # 比较括号右侧短线
  geom_segment(
    data = ann,
    aes(
      x = x2,
      xend = x2,
      y = tip_y,
      yend = bracket_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +

  # 实际 P 值
  geom_text(
    data = ann,
    aes(
      x = text_x,
      y = text_y,
      label = p_label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    size = 3.5,
    fontface = "plain"
  ) +

  scale_fill_manual(
    values = colors
  ) +

  scale_pattern_manual(
    values = patterns
  ) +

  scale_pattern_angle_manual(
    values = angles
  ) +

  scale_x_discrete(
    limits = groups,
    labels = groups
  ) +

  scale_y_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = function(limits, breaks) {
      head(breaks, -1) + diff(breaks) / 2
    },
    labels = scales::label_number()
  ) +

  coord_cartesian(
    ylim = c(0, upper),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = expression(C[1] ~ F[max] ~ "(a.u.)"),
    title = "C1 FMax"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    legend.position = "none",

    panel.border = element_rect(
      colour = "black",
      linewidth = 0.55
    ),

    panel.grid.major.x = element_line(
      colour = "grey89",
      linewidth = 0.35
    ),

    panel.grid.major.y = element_line(
      colour = "grey89",
      linewidth = 0.35
    ),

    panel.grid.minor.x = element_blank(),

    panel.grid.minor.y = element_line(
      colour = "grey89",
      linewidth = 0.18
    ),

    axis.text.x = element_text(
      size = 11,
      angle = 28,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),

    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),

    axis.title.y = element_text(
      size = 13,
      face = "bold",
      colour = "black"
    ),

    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),

    plot.margin = margin(7, 9, 9, 7)
  )

print(p)


# ============================================================
# 10. 导出统计结果
# ============================================================

write.csv(
  summary_data,
  file.path(
    output_dir,
    "C1_FMax_Mean_SD_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  t1way_result,
  file.path(
    output_dir,
    "C1_FMax_WRS2_t1way_overall.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  lincon_result,
  file.path(
    output_dir,
    "C1_FMax_WRS2_lincon_all_6_pairs_rawP_no_adjustment.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  ann,
  file.path(
    output_dir,
    "C1_FMax_WRS2_significant_pairs_for_plot.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 11. 保存 PNG/PDF，导出尺寸 108 × 115 mm
# ============================================================

png_file <- file.path(
  output_dir,
  "C1_FMax_散点_WRS2_t1way_lincon_未校正P值_SD.png"
)

pdf_file <- file.path(
  output_dir,
  "C1_FMax_散点_WRS2_t1way_lincon_未校正P值_SD.pdf"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 108,
  height = 115,
  units = "mm",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 108,
  height = 115,
  units = "mm",
  bg = "white"
)


# ============================================================
# 12. 控制台输出与导出检查
# ============================================================

cat("\n============================================================\n")
cat("C1 FMax WRS2 差异检验与绘图完成\n")
cat("============================================================\n")

cat(
  "总体检验：WRS2::t1way()；tr = ",
  TRIM_LEVEL,
  "\n",
  sep = ""
)

cat(
  "两两比较：WRS2::lincon()；method = none；全部 6 对；不进行多重比较校正\n"
)

cat(
  "总体 t1way P 值：",
  format.pval(t1way_result$p_t1way, digits = 6),
  "\n",
  sep = ""
)

cat(
  "标注规则：总体 P < 0.05 且两两比较原始 P < 0.05\n"
)

cat(
  "P 值格式：保留三位小数；P < 0.001 显示为 <0.001\n"
)

cat(
  "PNG：",
  png_file,
  "\n",
  sep = ""
)

cat(
  "PDF：",
  pdf_file,
  "\n",
  sep = ""
)

cat("============================================================\n")