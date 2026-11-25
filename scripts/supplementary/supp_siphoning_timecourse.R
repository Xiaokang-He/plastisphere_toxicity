# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.


# ============================================================
# 虹吸效率柱状图（CTR + Day0 + Day7 + Day14 + Day30 + Day60）
#
# 柱高：Mean；误差棒：Mean ± SD
# 散点：每个独立重复的原始值，仅沿水平方向轻微错开
#
# 统计：WRS2 t1way() + lincon()，不进行多重比较校正
# 截尾比例：tr = 0.20（20% trimmed means）
# 比较范围：六组共15次两两比较
# 标注条件：总体 t1way P < 0.05 且 lincon 原始 P < 0.05
# 标注内容：P 保留三位小数；P < 0.001 显示为 <0.001
#
# Day0：
# 浅红色 #FFB3C0
# 竖直条纹，柱中心对称排列
# 竖线间距 = Day60 网格线间距 × 1.5
#
# 其他组：保持原来的颜色和填充图案
# 输出：PNG/PDF，120 × 115 mm（按原有 ggsave 参数）
# ============================================================

library(dplyr)
library(ggplot2)
library(ggpattern)
library(WRS2)

# ============================================================
# 1. 读取数据
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- input_dir
output_dir <- file.path(results_dir, "supp_siphoning_timecourse")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c(
  "CTR",
  "Day0",
  "Day7",
  "Day14",
  "Day30",
  "Day60"
)

dat <- read.csv(
  file.path(input_dir, "虹吸效率.csv"),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) %>%
  transmute(
    Group = factor(
      gsub("\\s+", "", as.character(Group)),
      levels = groups
    ),
    Rep = Rep,
    Value = .data[["m (mL/clam/h)"]]
  )

# 检查分组
if (any(is.na(dat$Group))) {
  stop("存在无法识别的分组名称，请检查 CSV 文件。")
}

# 计算均值和标准差
sm <- dat %>%
  group_by(Group) %>%
  summarise(
    n = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 2. WRS2 t1way + lincon 差异检验
#
# 总体检验：WRS2::t1way()
# 两两比较：WRS2::lincon()
# 截尾比例：tr = 0.20
# method = "none"：不进行多重比较校正
# 六组共 15 次两两比较
# ============================================================

TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05

# 检查数据类型及有效观测值
if (!is.numeric(dat$Value)) {
  stop("虹吸效率 Value 列必须为数值型，请检查 CSV 文件。")
}

if (any(!is.na(dat$Value) & !is.finite(dat$Value))) {
  stop("虹吸效率数据中存在 Inf 或 -Inf，请检查原始数据。")
}

analysis_dat <- dat %>%
  filter(
    !is.na(Value),
    is.finite(Value)
  )

# WRS2 检验要求每组均有有效独立观测值
group_n <- table(
  factor(analysis_dat$Group, levels = groups)
)

if (any(group_n < 2)) {
  stop(
    "WRS2 检验要求每组至少有 2 个有效观测值。当前各组样本量：\n",
    paste(
      names(group_n),
      as.integer(group_n),
      sep = " = ",
      collapse = "; "
    )
  )
}

cat("\n========== WRS2 分组有效样本量 ==========\n")
print(group_n)


# ------------------------------------------------------------
# 2.1 总体检验：WRS2 t1way
# ------------------------------------------------------------

t1way_fit <- tryCatch(
  WRS2::t1way(
    Value ~ Group,
    data = analysis_dat,
    tr = TRIM_LEVEL
  ),
  error = function(e) {
    stop(
      "WRS2::t1way() 运行失败：\n",
      conditionMessage(e)
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
  stop("WRS2 t1way 返回了无效的总体 P 值。")
}

cat("\n========== WRS2 t1way 总体检验 ==========\n")
print(t1way_result)


# ------------------------------------------------------------
# 2.2 两两比较：WRS2 lincon
#
# method = "none"：P 值不进行多重比较校正。
# comp 的列顺序为：
# Group、Group、psihat、ci.lower、ci.upper、p.value
# 第1、2列是组别序号；第6列是未校正 P 值。
# ------------------------------------------------------------

lincon_fit <- tryCatch(
  WRS2::lincon(
    Value ~ Group,
    data = analysis_dat,
    tr = TRIM_LEVEL,
    method = "none"
  ),
  error = function(e) {
    stop(
      "WRS2::lincon() 运行失败：\n",
      conditionMessage(e)
    )
  }
)

comp <- as.data.frame(lincon_fit$comp)

if (ncol(comp) < 6) {
  stop(
    "WRS2 lincon 返回的结果列数不足 6 列，无法提取两两比较结果。"
  )
}

group_index_1 <- as.integer(comp[[1]])
group_index_2 <- as.integer(comp[[2]])

group_names <- lincon_fit$fnames

if (is.null(group_names) ||
    length(group_names) != length(groups)) {
  group_names <- levels(droplevels(analysis_dat$Group))
}

if (anyNA(group_index_1) ||
    anyNA(group_index_2) ||
    any(group_index_1 < 1 | group_index_1 > length(group_names)) ||
    any(group_index_2 < 1 | group_index_2 > length(group_names))) {
  stop(
    "无法将 lincon 的组别序号映射回 CTR、Day0、Day7、Day14、Day30、Day60。"
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
    x2 = match(group2, groups),
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
  ) %>%
  arrange(
    x2 - x1,
    x1
  )

expected_comparisons <- choose(length(groups), 2)

if (nrow(lincon_result) != expected_comparisons) {
  stop(
    "lincon 未返回预期的全部 ",
    expected_comparisons,
    " 次两两比较；实际返回 ",
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

cat("\n========== WRS2 lincon 全部 15 次两两比较 ==========\n")
print(
  as.data.frame(lincon_result),
  row.names = FALSE
)


# ------------------------------------------------------------
# 2.3 P 值格式化
#
# P >= 0.001：保留三位小数，例如 P = 0.023
# P < 0.001：显示为 P < 0.001
# 不使用科学计数法。
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

# ============================================================
# 3. 显著性标注位置
# ============================================================

base_y <- max(
  c(sm$mean + sm$sd, dat$Value),
  na.rm = TRUE
)

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
  ) %>%
  mutate(
    bracket_y = base_y * (
      1.06 + (row_number() - 1) * 0.100
    ),
    tip_y = bracket_y - base_y * 0.018,
    text_y = bracket_y + base_y * 0.025,
    text_x = (x1 + x2) / 2
  )

upper <- max(c(
  base_y * 1.105,
  ann$text_y + base_y * 0.070
))

# ============================================================
# 4. 配色和填充图案
# ============================================================

colors <- c(
  CTR   = "#b6dcd5",
  Day0  = "#FFB3C0",
  Day7  = "#92CDE3",
  Day14 = "#95C3AA",
  Day30 = "#82A2CD",
  Day60 = "#CAB9D5"
)

patterns <- c(
  CTR   = "stripe",
  Day0  = "centered_vertical",
  Day7  = "stripe",
  Day14 = "stripe",
  Day30 = "crosshatch",
  Day60 = "centered_grid"
)

angles <- c(
  CTR   = 0,
  Day0  = 90,
  Day7  = 45,
  Day14 = 135,
  Day30 = 45,
  Day60 = 0
)

# ============================================================
# 5. Day60：原有居中横竖网格函数
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
  
  # 保持原有 Day60 网格间距计算方式
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
    
    # 竖线
    grid::segmentsGrob(
      x0 = centered_positions(xlim),
      x1 = centered_positions(xlim),
      y0 = ylim[1],
      y1 = ylim[2],
      default.units = "mm",
      gp = gp
    ),
    
    # 横线
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

# ============================================================
# 6. Day0：居中竖线函数
# 以柱子中心为基准，左右对称排列
# ============================================================

centered_vertical <- function(
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
  
  # 直接使用 Day0 的实际间距
  # 注意：这里不再除以2
  step <- grid::convertWidth(
    grid::unit(
      params$pattern_spacing,
      params$pattern_units
    ),
    "mm",
    valueOnly = TRUE
  )
  
  # 以柱子中心为对称轴
  n <- floor(diff(xlim) / (2 * step))
  
  x_positions <- mean(xlim) +
    seq.int(-n, n) * step
  
  line_width <- params$pattern_linewidth
  
  if (is.null(line_width)) {
    line_width <- params$pattern_size
  }
  
  gp <- grid::gpar(
    col = params$pattern_colour,
    lwd = line_width * 72.27 / 25.4
  )
  
  # 只绘制竖线
  grid::segmentsGrob(
    x0 = x_positions,
    x1 = x_positions,
    y0 = ylim[1],
    y1 = ylim[2],
    default.units = "mm",
    gp = gp
  )
}

# ============================================================
# 7. 注册两种自定义图案
# ============================================================

pattern_funcs <- getOption(
  "ggpattern_geometry_funcs",
  list()
)

pattern_funcs$centered_grid <- centered_grid

pattern_funcs$centered_vertical <- centered_vertical

options(
  ggpattern_geometry_funcs = pattern_funcs
)

# ============================================================
# 8. 填充间距设置
# ============================================================

# Day60 原始 pattern_spacing = 0.08
# centered_grid 内部除以2
# 因此实际网格间距 = 0.04

day60_spacing <- 0.08 / 2

# 按照要求：Day0 间距为 Day60 的1.5倍
day0_spacing <- day60_spacing * 1

cat(
  "\nDay60 网格线间距：", day60_spacing,
  "\nDay0 竖线间距：", day0_spacing,
  "\nDay0 / Day60 间距比例：",
  day0_spacing / day60_spacing,
  "\n"
)

# ============================================================
# 9. 绘图
# ============================================================

p <- ggplot(
  sm,
  aes(
    x = Group,
    y = mean,
    fill = Group
  )
) +
  
  # ----------------------------------------------------------
# 9.1 原有五组：保持原始填充方式
# ----------------------------------------------------------

geom_col_pattern(
  data = sm %>% filter(Group != "Day0"),
  
  aes(
    pattern = Group,
    pattern_angle = Group
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
  
  # ----------------------------------------------------------
# 9.2 Day0：浅红色 + 居中竖线
# ----------------------------------------------------------

geom_col_pattern(
  data = sm %>% filter(Group == "Day0"),
  
  aes(
    pattern = Group,
    pattern_angle = Group
  ),
  
  pattern_spacing = day0_spacing,
  pattern_units = "snpc",
  pattern_density = 0.35,
  pattern_fill = NA,
  pattern_fill2 = NA,
  pattern_colour = "black",
  pattern_size = 0.18,
  pattern_grid = "square",
  
  width = 0.78,
  colour = "black",
  linewidth = 0.55
) +
  
  # ----------------------------------------------------------
# 9.3 Mean ± SD
# ----------------------------------------------------------

geom_errorbar(
  aes(
    ymin = pmax(0, mean - sd),
    ymax = mean + sd
  ),
  width = 0.17,
  linewidth = 0.55
) +
  
  # ----------------------------------------------------------
# 9.4 原始散点
# ----------------------------------------------------------

geom_point(
  data = dat,
  
  aes(
    x = Group,
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
  
  # ----------------------------------------------------------
# 9.5 显著性括号
# ----------------------------------------------------------

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
  
  # ----------------------------------------------------------
# 9.6 实际 P 值
# ----------------------------------------------------------

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
  
  # ----------------------------------------------------------
# 9.7 颜色和图案映射
# ----------------------------------------------------------

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
    
    labels = c(
      CTR = "Control",
      Day0 = "Day0",
      Day7 = "Day7",
      Day14 = "Day14",
      Day30 = "Day30",
      Day60 = "Day60"
    )
  ) +
  
  # ----------------------------------------------------------
# 9.8 Y轴
# ----------------------------------------------------------

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
    
    y = expression(atop(
      "Siphoning efficiency",
      "(mL " * plain(clam)^{-1} * " " * plain(h)^{-1} * ")"
    )),
    
    title = "Siphoning efficiency"
  ) +
  
  # ----------------------------------------------------------
# 9.9 主题
# ----------------------------------------------------------

theme_bw(base_size = 12) +
  
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

# ============================================================
# 10. 显示图形
# ============================================================

print(p)

# ============================================================
# 11. 导出 PNG 和 PDF
# 文件夹、文件名和尺寸保持不变
# ============================================================

ggsave(
  file.path(
    output_dir,
    "虹吸效率_柱状图_散点_实际P值_未校正_SD.png"
  ),
  
  p,
  
  width = 120,
  height = 115,
  units = "mm",
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "虹吸效率_柱状图_散点_实际P值_未校正_SD.pdf"
  ),
  
  p,
  
  width = 120,
  height = 115,
  units = "mm",
  bg = "white"
)

# ============================================================
# 12. 输出统计结果
# ============================================================

cat("\n========== WRS2 t1way 总体检验 ==========\n")
print(t1way_result)

cat("\n========== WRS2 lincon 全部 15 次两两比较（未校正 P 值） ==========\n")
print(
  as.data.frame(lincon_result),
  row.names = FALSE
)

cat("\n========== 满足显著性标注条件的两两比较 ==========\n")

if (nrow(ann) > 0) {
  print(
    as.data.frame(
      ann %>%
        select(
          group1,
          group2,
          p_raw,
          p_label
        )
    ),
    row.names = FALSE
  )
} else {
  cat(
    "没有比较同时满足：总体 t1way P < 0.05 且 lincon 原始 P < 0.05。\n"
  )
}

cat(
  "\n图形已保存至：\n",
  output_dir,
  "\n"
)