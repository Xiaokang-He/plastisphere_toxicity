# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# EPS 定量：Protein、Polysaccharide、Humic substances
#
# 柱高：Mean；误差棒：Mean ± SD
# 散点：每个独立重复的原始值，仅沿水平方向轻微错开
#
# 统计：
#   WRS2::t1way() + WRS2::lincon()
#   t1way：四组总体异方差稳健单因素检验（tr = 0.20）
#   lincon：四组全部 6 次两两比较
#   method = "none"：两两比较P值不进行多重比较校正
#   标注条件：该指标总体 t1way P < 0.05 且 lincon 原始 P < 0.05
#   标注内容：实际 P 值；保留三位小数；小于 0.001 时标注为 <0.001
#
# 填充：
#   Day7 45°斜线；Day14 135°斜线；
#   Day30 45°交叉网格；Day60 居中的横竖网格
#
# 输出：3列 × 1行；页面570 × 200 mm；Cairo矢量PDF、300 dpi PNG及单图PDF
# ============================================================


# ============================================================
# 0. 加载依赖
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggpattern",
  "WRS2",
  "cowplot",
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpattern)
  library(WRS2)
  library(cowplot)
})

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(
    "Times New Roman" = grDevices::windowsFont("Times New Roman")
  )
}

# 防止此前R会话遗留的showtext自动模式将文字转换成路径。
if (requireNamespace("showtext", quietly = TRUE)) {
  try(showtext::showtext_auto(FALSE), silent = TRUE)
}


# ============================================================
# 1. 参数、路径与数据读取
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- input_dir

input_file <- file.path(
  base_dir,
  "EPS含量.csv"
)

if (!file.exists(input_file)) {
  stop(
    "找不到原始数据文件：\n",
    input_file,
    "\n请确认文件名和路径。"
  )
}

if (!dir.exists(base_dir)) {
  stop("数据文件夹不存在：\n", base_dir)
}

# WRS2 统计参数
TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05

# 与肠道微生物16S α多样性代码一致的横向一排导出规格
PAGE_WIDTH_MM <- 570
PAGE_HEIGHT_MM <- 200
N_COLS <- 3
N_ROWS <- 1
CELL_WIDTH_MM <- PAGE_WIDTH_MM / N_COLS
CELL_HEIGHT_MM <- PAGE_HEIGHT_MM / N_ROWS
PNG_DPI <- 300
FONT_FAMILY <- "Times New Roman"

# 固定随机种子，保证图中散点水平抖动位置可复现。
set.seed(2026)

dat <- read.csv(
  input_file,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

# 输入表的五列依次为：时间、重复、蛋白质、多糖、腐殖质
if (ncol(dat) != 5) {
  stop(
    "EPS含量.csv 应包含 5 列：时间、重复、蛋白质、多糖、腐殖质。"
  )
}

names(dat) <- c(
  "Time", "Replicate", "Protein", "Polysaccharide", "Humic"
)

groups <- c("Day7", "Day14", "Day30", "Day60")
metrics <- c("Protein", "Polysaccharide", "Humic")

dat$Time <- factor(
  gsub("\\s+", "", as.character(dat$Time)),
  levels = groups
)

# 检查时间组别是否均能正确识别
if (any(is.na(dat$Time))) {
  bad_times <- unique(
    as.character(dat$Time[is.na(dat$Time)])
  )
  stop(
    "存在无法识别的时间组别。请检查 EPS含量.csv 中的 Time 列，",
    "应为 Day7、Day14、Day30、Day60。"
  )
}

# 检查三个指标是否为数值型
for (m in metrics) {
  if (!is.numeric(dat[[m]])) {
    stop(
      "指标列 ", m, " 不是数值型。请检查 EPS含量.csv 中是否存在",
      "非数值字符、单位或异常空格。"
    )
  }
}

long <- dat %>%
  pivot_longer(
    all_of(metrics),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = factor(Index, levels = metrics),
    Time = factor(Time, levels = groups)
  )

if (any(!is.finite(long$Value))) {
  stop(
    "Protein、Polysaccharide 或 Humic 中存在 NA、NaN、Inf 或 -Inf。",
    "\n请先检查原始数据，确保三个指标均为有限数值。"
  )
}

# 确保每个指标均包含四个时间组
group_check <- long %>%
  count(Index, Time, .drop = FALSE)

if (any(group_check$n == 0)) {
  print(group_check)
  stop("至少一个指标缺少某个时间组的数据，无法完成四组 WRS2 分析。")
}


# ============================================================
# 2. 计算各组 Mean 和 SD；误差棒使用 Mean ± SD
# ============================================================

sm <- long %>%
  group_by(Index, Time) %>%
  summarise(
    n = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

if (any(!is.finite(sm$sd))) {
  stop(
    "至少一个指标的某个时间组无法计算 SD。",
    "\n请检查该组是否只有一个有效重复或存在异常数据。"
  )
}


# ============================================================
# 3. WRS2：t1way 总体检验 + lincon 全部两两比较
#
# 重要：
#   1) 每个指标独立进行一次 t1way 总体检验。
#   2) 每个指标独立进行一次 lincon，获得四组全部 6 次比较。
#   3) method = \"none\" 明确关闭 lincon 的 P 值多重比较校正。
#   4) 图中标注仍保留原代码的双重条件：
#        总体 t1way P < 0.05
#        且该组间 lincon 原始 P < 0.05
#   5) 不把 Mean ± SD 改为 trimmed mean；柱高和误差棒仍按原代码计算。
# ============================================================

run_wrs2_one <- function(data_one, metric_name, tr = TRIM_LEVEL) {

  d <- data_one %>%
    filter(Index == metric_name) %>%
    mutate(
      Time = factor(as.character(Time), levels = groups)
    ) %>%
    droplevels()

  # 统计前检查四组及每组样本量
  observed_groups <- levels(d$Time)
  n_by_group <- table(
    factor(d$Time, levels = groups)
  )

  if (length(observed_groups) != length(groups) ||
      any(n_by_group == 0)) {
    stop(
      "指标 ", metric_name,
      " 未包含完整的 Day7、Day14、Day30、Day60 四组。"
    )
  }

  if (any(n_by_group < 2)) {
    stop(
      "指标 ", metric_name,
      " 至少有一个时间组少于 2 个有效观测值，",
      "无法可靠执行当前 WRS2 分析。各组 n：",
      paste(names(n_by_group), n_by_group, sep = "=", collapse = "; ")
    )
  }

  # ----------------------------------------------------------
  # 3.1 总体检验：异方差稳健单因素截尾均值检验
  # ----------------------------------------------------------

  t1_fit <- tryCatch(
    WRS2::t1way(
      Value ~ Time,
      data = d,
      tr = tr
    ),
    error = function(e) {
      stop(
        "指标 ", metric_name, " 的 WRS2::t1way() 运行失败：\n",
        conditionMessage(e),
        "\n请检查各组样本量、组内重复值及 Winsorized variance。"
      )
    }
  )

  t1_table <- data.frame(
    Index = metric_name,
    Test = "WRS2 t1way",
    Trim_level = tr,
    Statistic_F = as.numeric(t1_fit$test),
    df1 = as.numeric(t1_fit$df1),
    df2 = as.numeric(t1_fit$df2),
    p_t1way = as.numeric(t1_fit$p.value),
    Effect_size = as.numeric(t1_fit$effsize),
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # 3.2 两两比较：四组全部 6 次 lincon
  #
  # WRS2::lincon() 的 comp 矩阵包含：
  #   第1列、第2列：组别在因子水平中的序号
  #   第3列：psihat（截尾均值差）
  #   第4列：ci.lower
  #   第5列：ci.upper
  #   第6列：p.value
  #
  # method = "none"：第6列保留未进行多重比较校正的原始 P 值。
  # ----------------------------------------------------------

  lin_fit <- tryCatch(
    WRS2::lincon(
      Value ~ Time,
      data = d,
      tr = tr,
      method = "none"
    ),
    error = function(e) {
      stop(
        "指标 ", metric_name, " 的 WRS2::lincon() 运行失败：\n",
        conditionMessage(e),
        "\n请检查各组样本量、组内重复值及 Winsorized variance。"
      )
    }
  )

  comp <- as.data.frame(lin_fit$comp)

  if (ncol(comp) < 6) {
    stop(
      "指标 ", metric_name,
      " 的 lincon 结果列数不足 6 列，无法正确提取组别和 P 值。"
    )
  }

  group_index_1 <- as.integer(comp[[1]])
  group_index_2 <- as.integer(comp[[2]])

  group_names <- lin_fit$fnames

  if (is.null(group_names) || length(group_names) != length(groups)) {
    # 兼容部分 WRS2 版本：按输入因子水平顺序恢复组名
    group_names <- levels(d$Time)
  }

  if (anyNA(group_index_1) ||
      anyNA(group_index_2) ||
      any(group_index_1 < 1 | group_index_1 > length(group_names)) ||
      any(group_index_2 < 1 | group_index_2 > length(group_names))) {
    stop(
      "指标 ", metric_name,
      " 的 lincon 结果无法映射回原始组别，请检查 WRS2 输出结构。"
    )
  }

  lin_table <- data.frame(
    Index = metric_name,
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
    )

  if (nrow(lin_table) != choose(length(groups), 2)) {
    stop(
      "指标 ", metric_name,
      " 的 lincon 未返回预期的全部 6 次两两比较。"
    )
  }

  if (any(!is.finite(lin_table$p_raw)) ||
      any(lin_table$p_raw < 0 | lin_table$p_raw > 1)) {
    stop(
      "指标 ", metric_name,
      " 的 lincon 结果中存在无效 P 值。"
    )
  }

  list(
    t1way = t1_table,
    lincon = lin_table,
    t1way_object = t1_fit,
    lincon_object = lin_fit
  )
}


# 对三个指标分别执行 WRS2 分析
wrs2_results <- setNames(
  lapply(
    metrics,
    function(m) {
      run_wrs2_one(
        data_one = long,
        metric_name = m,
        tr = TRIM_LEVEL
      )
    }
  ),
  metrics
)

# 总体检验结果
t1way_result <- bind_rows(
  lapply(wrs2_results, `[[`, "t1way")
)

# 所有两两比较结果：每个指标 6 次，共 18 次
lincon_result <- bind_rows(
  lapply(wrs2_results, `[[`, "lincon")
) %>%
  left_join(
    t1way_result %>%
      select(Index, p_t1way),
    by = "Index"
  ) %>%
  arrange(
    factor(Index, levels = metrics),
    x2 - x1,
    x1
  )

# ------------------------------------------------------------
# 3.3 P 值格式
#
# P >= 0.001：保留三位小数，例如 P = 0.023
# P < 0.001：统一标注为 P < 0.001，不使用科学计数法
# P 字符斜体，等号和数值为常规体。
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

      sprintf("italic(P) == '%.3f'", z)
    },
    character(1)
  )
}


# ------------------------------------------------------------
# 3.4 显著性筛选
#
# 保留原代码的双重筛选规则：
#   总体 t1way P < 0.05
#   且对应 lincon 未校正 P < 0.05
# ------------------------------------------------------------

sig <- lincon_result %>%
  filter(
    p_t1way < ALPHA_LEVEL,
    p_raw < ALPHA_LEVEL
  ) %>%
  mutate(
    p_label = format_p_label(p_raw)
  ) %>%
  arrange(
    factor(Index, levels = metrics),
    x2 - x1,
    x1
  )


# ============================================================
# 4. 配色与填充：保持原代码不变
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


# Day60：按柱中心对称排列，横纵采用相同实际间距。
# 此自定义图案用于本代码中的矩形柱。
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

  # 与模板中透明条带、density = 0.50 对应的细线间距一致
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
# 5. 单图函数
#
# 除统计结果对象及P值字段外，保留原代码的绘图设置。
# ============================================================

make_plot <- function(metric, title, y_label) {

  s <- filter(sm, Index == metric)
  d <- filter(long, Index == metric)

  # 括号放在误差棒和全部原始散点上方
  base_y <- max(
    c(
      s$mean + s$sd,
      d$Value
    ),
    na.rm = TRUE
  )

  ann <- sig %>%
    filter(Index == metric) %>%
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
      base_y * 1.140,
      ann$text_y + base_y * 0.070
    ),
    na.rm = TRUE
  )

  ggplot(
    s,
    aes(
      x = Time,
      y = mean,
      fill = Time
    )
  ) +

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

    geom_errorbar(
      aes(
        ymin = pmax(0, mean - sd),
        ymax = mean + sd
      ),
      width = 0.17,
      linewidth = 0.55
    ) +

    # 白色填充、黑色边框，便于在条纹柱上辨认每个重复
    geom_point(
      data = d,
      aes(
        x = Time,
        y = Value
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
      fontface = "plain",
      family = FONT_FAMILY
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
      y = y_label,
      title = title
    ) +

    theme_bw(
      base_size = 12,
      base_family = FONT_FAMILY
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
}


# ============================================================
# 6. 三面板图
# ============================================================

plots <- list(
  Protein = make_plot(
    "Protein",
    "Protein",
    "Protein (mg/g PS)"
  ),

  Polysaccharide = make_plot(
    "Polysaccharide",
    "Polysaccharide",
    "Polysaccharide (mg/g PS)"
  ),

  Humic = make_plot(
    "Humic",
    "Humic substances",
    "Humic substances\n(mg HA-equivalent/g PS)"
  )
)

# 三个面板固定按 Protein → Polysaccharide → Humic 从左至右排列
fig <- cowplot::plot_grid(
  plotlist = plots,
  ncol = N_COLS,
  nrow = N_ROWS,
  labels = c("(a)", "(b)", "(c)"),
  label_size = 14
)

print(fig)


# ============================================================
# 7. 按照α多样性代码的风格导出合并图
#
# 页面：570 × 200 mm
# 排列：3列 × 1行；每个单元格190 × 200 mm
# 输出：
#   1) Cairo矢量PDF（主交付文件）
#   2) PNG预览（300 dpi）
#   3) 三张单独PDF（每张190 × 200 mm）
#
# 本版本仅导出PDF和PNG，不生成SVG文件。
# ============================================================

output_dir <- file.path(
  results_dir,
  "supp_fig2_eps_quantification"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(output_dir)) {
  stop(
    "无法创建输出文件夹，请检查 输出目录及文件夹权限：\n",
    output_dir
  )
}

output_stem <- paste0(
  "EPS定量_Protein_Polysaccharide_Humic_",
  "WRS2_t1way_lincon_三位小数_570x200mm"
)

pdf_file <- file.path(
  output_dir,
  paste0(output_stem, "_Cairo矢量PDF.pdf")
)

png_file <- file.path(
  output_dir,
  paste0(output_stem, "_预览.png")
)


# ------------------------------------------------------------
# 7.1 导出主图：Cairo矢量PDF
# ------------------------------------------------------------

if (!isTRUE(capabilities("cairo"))) {
  stop(
    "当前 R 环境不支持 Cairo PDF 设备，无法按当前设置导出 Cairo 矢量 PDF。",
    "\n请检查 R 的 Cairo 支持后重新运行。"
  )
}

grDevices::cairo_pdf(
  filename = pdf_file,
  width = PAGE_WIDTH_MM / 25.4,
  height = PAGE_HEIGHT_MM / 25.4,
  family = FONT_FAMILY,
  symbolfamily = FONT_FAMILY,
  pointsize = 12,
  bg = "white"
)

tryCatch(
  {
    print(fig)
  },
  finally = {
    grDevices::dev.off()
  }
)

# 检查PDF是否成功生成且文件非空
if (!file.exists(pdf_file)) {
  stop(
    "PDF文件未成功生成，请检查输出目录权限及磁盘可用空间：\n",
    pdf_file
  )
}

if (file.info(pdf_file)$size <= 0) {
  stop("PDF文件已生成，但文件大小为0：\n", pdf_file)
}

cat("\n========== Cairo PDF导出检查 ==========\n")
cat("PDF文件：", pdf_file, "\n")
cat("PDF文件大小：", round(file.info(pdf_file)$size / 1024, 1), " KB\n")
cat("页面尺寸：", PAGE_WIDTH_MM, " × ", PAGE_HEIGHT_MM, " mm\n")
cat("排列方式：", N_COLS, "列 × ", N_ROWS, "行\n")
cat("字体设置：", FONT_FAMILY, "\n")


# ------------------------------------------------------------
# 7.2 导出PNG预览（300 dpi）
# ------------------------------------------------------------

grDevices::png(
  filename = png_file,
  width = PAGE_WIDTH_MM,
  height = PAGE_HEIGHT_MM,
  units = "mm",
  res = PNG_DPI,
  type = "cairo",
  family = FONT_FAMILY,
  pointsize = 12,
  bg = "white"
)

tryCatch(
  {
    print(fig)
  },
  finally = {
    grDevices::dev.off()
  }
)

if (!file.exists(png_file)) {
  warning(
    "PNG预览未成功生成，但主PDF已经导出：\n",
    png_file
  )
}


# ============================================================
# 8. 输出3张单独PDF
#    每张单图页面：190 × 200 mm
# ============================================================

single_dir <- file.path(
  output_dir,
  "EPS定量_单张指标PDF_190x200mm"
)

dir.create(
  single_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(single_dir)) {
  stop("无法创建单图PDF输出文件夹：\n", single_dir)
}

for (nm in names(plots)) {

  single_file <- file.path(
    single_dir,
    paste0(nm, "_WRS2.pdf")
  )

  grDevices::cairo_pdf(
    filename = single_file,
    width = CELL_WIDTH_MM / 25.4,
    height = CELL_HEIGHT_MM / 25.4,
    family = FONT_FAMILY,
    symbolfamily = FONT_FAMILY,
    pointsize = 12,
    bg = "white"
  )

  tryCatch(
    {
      print(plots[[nm]])
    },
    finally = {
      grDevices::dev.off()
    }
  )

  if (!file.exists(single_file) || file.info(single_file)$size <= 0) {
    stop("单张PDF未成功生成：\n", single_file)
  }
}


# 8. 控制台输出统计结果
# ============================================================

cat("\n")
cat("============================================================\n")
cat("EPS 定量：WRS2 统计分析完成\n")
cat("============================================================\n")

cat(
  "总体检验：WRS2::t1way()；tr = ",
  TRIM_LEVEL,
  "\n",
  sep = ""
)

cat(
  "两两比较：WRS2::lincon()；method = \"none\"\n"
)

cat(
  "每个指标：4组全部 6 次两两比较；3个指标共 18 次比较\n"
)

cat(
  "图中标注条件：总体 t1way P < ",
  ALPHA_LEVEL,
  " 且对应 lincon 原始 P < ",
  ALPHA_LEVEL,
  "\n",
  sep = ""
)

cat(
  "柱高：Mean；误差棒：Mean ± SD；散点：原始独立重复值\n"
)

cat(
  "P值显示：原始未校正P值；保留三位小数；P < 0.001 显示为 <0.001；不使用星号\n"
)

cat("页面：", PAGE_WIDTH_MM, " × ", PAGE_HEIGHT_MM, " mm\n", sep = "")
cat("排列：", N_COLS, "列 × ", N_ROWS, "行\n", sep = "")
cat("每个单元格：", CELL_WIDTH_MM, " × ", CELL_HEIGHT_MM, " mm\n", sep = "")
cat("PNG：", PNG_DPI, " dpi\n", sep = "")
cat("主图：Cairo矢量PDF\n")
cat("PNG预览：", PNG_DPI, " dpi\n", sep = "")

cat("\n========== WRS2 t1way 总体检验结果 ==========\n")
print(t1way_result, row.names = FALSE)

cat("\n========== WRS2 lincon 全部两两比较结果 ==========\n")
print(lincon_result, row.names = FALSE)

cat("\n========== 满足双重条件的显著比较 ==========\n")
if (nrow(sig) == 0) {
  cat("没有同时满足总体 t1way P < 0.05 和两两比较原始 P < 0.05 的比较。\n")
} else {
  print(
    sig %>%
      select(
        Index,
        group1,
        group2,
        p_t1way,
        p_raw,
        p_label
      ),
    row.names = FALSE
  )
}

cat("\n主图Cairo矢量PDF：\n", pdf_file, "\n", sep = "")
cat("\nPNG预览：\n", png_file, "\n", sep = "")
cat("\n单张PDF文件夹：\n", single_dir, "\n", sep = "")
cat("\n全部输出文件夹：\n", output_dir, "\n", sep = "")

cat("============================================================\n")