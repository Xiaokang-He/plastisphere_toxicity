# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# Fig.1 生理生化指标柱状图：WRS2紧凑标注及初始Y轴范围保留版
#
# 页面：950 × 400 mm
# 排列：5列 × 2行；每个单元格190 × 200 mm
#
# 每张图的严格尺寸：
#   XY绘图区：89.277 × 136.515 mm
#   标题框：  89.277 × 18.238 mm
#
# 字体（Times New Roman）：
#   标题：36.59 pt，加粗
#   Y轴刻度：28.4 pt
#   Y轴标题：28.4 pt，加粗
#   X轴刻度：28.4 pt，倾斜30°
#   显著性P值：默认18 pt；空间不足时按panel自适应缩小至最低10 pt；仅P字符斜体
#
# 线宽：
#   柱边框、误差棒、XY绘图区边框、主网格线、显著性括号：1.788 pt
#   次网格线：0.886 pt
#
# 柱状图：
#   Mean ± SD
#   每个独立生物学重复显示一个白色填充、黑色边框散点
#   散点仅水平方向轻微错开，纵坐标保留原始值
#
# 统计：
#   WRS2 t1way() + lincon()
#   t1way：五组总体异方差稳健单因素检验（tr = 0.2）
#   lincon：五组全部10组两两比较；P值不进行多重比较校正
#   原始 P < 0.05 显示实际P值，不使用星号
#   注意：每组 n = 3 时，tr = 0.2 不会实际截尾任何观测值。
#
# 字体输出：
#   Illustrator编辑版：svglite输出原生SVG <text>对象
#   辅助PDF版：Cairo PDF保留字体文字对象
#   禁用showtext自动转换，避免文字变成矢量路径
#   不再输出任何字体轮廓化文件
#
# Y轴范围：优先保留最初版本plot_config中设置的ymax。
# 显著性标注优先压缩在既有Y轴上部空间内；必要时自适应缩小P值字号，
# 仅当最小字号仍无法容纳时才对Y轴上限作最小幅度扩展。
# 不再为了容纳显著性括号而整体大幅拉伸Y轴。
#
# 本版不使用ggplot默认轴文字布局，也不使用cowplot缩放。
# 直接抽取ggplot面板，在950 × 400 mm页面上手工放置：
#   标题框、Y轴标题、Y轴刻度、X轴刻度、XY绘图区。
# ============================================================


# ============================================================
# 0. 清理环境、加载依赖
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
# 1. 原始数据位置及输出位置
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- input_dir

input_file <- file.path(
  input_dir,
  "酶活数据_最终 - 8.30.csv"
)

if (!file.exists(input_file)) {
  stop(
    "找不到原始数据文件：\n",
    input_file,
    "\n请确认文件名和路径。"
  )
}

# SVG设备在部分Windows中文区域设置下可能无法稳定打开含中文的输出路径。
# 将结果统一输出到纯英文路径；原始CSV仍从base_dir读取。
output_dir <- file.path(results_dir, "fig1_physiology_biochemistry")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop("无法创建输出文件夹：\n", output_dir)
}

# 检查输出目录是否可写，避免到打开SVG/PDF设备时才报错。
test_file <- file.path(output_dir, "__write_test.tmp")
can_write <- tryCatch({
  writeLines("test", test_file, useBytes = TRUE)
  unlink(test_file)
  TRUE
}, error = function(e) FALSE)

if (!can_write) {
  stop("输出文件夹不可写，请检查路径权限：\n", output_dir)
}


# ============================================================
# 2. Times New Roman字体注册
#
# Illustrator可编辑SVG使用svglite输出原生<text>元素。
# 不使用showtext/sysfonts，不把文字转换为矢量轮廓。
#
# 重要：如果此前在同一个R会话中运行过showtext_auto(TRUE)，
# 它可能继续影响后续设备。这里显式关闭showtext自动转换。
# ============================================================

font_family <- "Times New Roman"

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(
    "Times New Roman" = grDevices::windowsFont("Times New Roman")
  )
}

# 防止此前会话遗留的showtext自动模式将文字转换成路径。
if (requireNamespace("showtext", quietly = TRUE)) {
  try(
    showtext::showtext_auto(FALSE),
    silent = TRUE
  )
}

if (!isTRUE(capabilities("cairo"))) {
  stop(
    "当前R环境不支持Cairo PDF设备。",
    "\n请使用支持Cairo的Windows版R。"
  )
}

# 检查Windows字体文件是否存在。
font_dir <- file.path(Sys.getenv("WINDIR"), "Fonts")
font_files <- c(
  regular = file.path(font_dir, "times.ttf"),
  bold = file.path(font_dir, "timesbd.ttf"),
  italic = file.path(font_dir, "timesi.ttf"),
  bolditalic = file.path(font_dir, "timesbi.ttf")
)

if (!all(file.exists(font_files))) {
  stop(
    "未找到完整的Times New Roman字体文件：\n",
    paste(font_files[!file.exists(font_files)], collapse = "\n")
  )
}

# ============================================================
# 3. 页面及图形物理尺寸
# ============================================================

PAGE_WIDTH_MM  <- 950
PAGE_HEIGHT_MM <- 400

N_COLS <- 5
N_ROWS <- 2

CELL_WIDTH_MM  <- PAGE_WIDTH_MM / N_COLS
CELL_HEIGHT_MM <- PAGE_HEIGHT_MM / N_ROWS

# 用户指定的XY绘图区尺寸
PANEL_WIDTH_MM  <- 89.277
PANEL_HEIGHT_MM <- 136.515

# 用户指定的标题框尺寸
TITLE_WIDTH_MM  <- 89.277
TITLE_HEIGHT_MM <- 18.238

# 文字字号，单位pt
TITLE_FONT_PT <- 36.59
AXIS_TEXT_PT  <- 28.4
AXIS_TITLE_PT <- 28.4

# 显著性P值字号：默认18 pt；若某一panel的原始Y轴范围内空间不足，
# 仅该panel自动缩小字号，最小不低于10 pt。
P_VALUE_FONT_PT <- 18
P_VALUE_MIN_FONT_PT <- 10
P_VALUE_FONT_STEP_PT <- 0.5

# 以下为紧凑标注的物理间距（mm）。
# 采用较小的文字净间距、首条括号间距和顶部余量，尽量保留初始Y轴范围。
P_VALUE_CLEARANCE_MM <- 0.35
P_VALUE_TEXT_HEIGHT_MM <- P_VALUE_FONT_PT * 25.4 / 72 * 0.90
P_BRACKET_EXTRA_GAP_MM <- 0.25
P_BRACKET_STEP_MM <- P_VALUE_TEXT_HEIGHT_MM +
  P_VALUE_CLEARANCE_MM + P_BRACKET_EXTRA_GAP_MM
P_FIRST_BRACKET_GAP_MM <- 0.50
P_VALUE_TOP_MARGIN_MM <- 0.50

# 用于估算不同P值字号的物理文字高度；绘制时仍以grid实际字体度量为准。
estimate_p_text_height_mm <- function(font_pt) {
  font_pt * 25.4 / 72 * 0.90
}

# 线宽，用户提供单位为pt
MAIN_LINE_PT  <- 1.788
MINOR_LINE_PT <- 0.886

# ggplot2的linewidth参数以mm输入，但内部会先转换为grid的lwd；
# grid的1个lwd单位实际对应1/96英寸（0.75 pt）。
# 因此，直接使用pt * 25.4 / 72会使最终PDF线宽缩小到约75%。
# 这里按ggplot2的.pt换算系数(72.27/25.4)和grid的0.75 pt/lwd
# 反算，确保PDF中的实际物理线宽等于指定的pt值。

pt_to_ggplot_linewidth_mm <- function(pt) {
  pt * 25.4 / (72.27 * 0.75)
}

# geom_text(size=)使用mm，仍采用常规pt到mm换算
pt_to_mm <- function(pt) {
  pt * 25.4 / 72
}

MAIN_LINE_MM  <- pt_to_ggplot_linewidth_mm(MAIN_LINE_PT)
MINOR_LINE_MM <- pt_to_ggplot_linewidth_mm(MINOR_LINE_PT)

# grid::gpar(lwd)直接使用1/96英寸单位
pt_to_grid_lwd <- function(pt) {
  pt * 96 / 72
}

MAIN_LINE_GRID <- pt_to_grid_lwd(MAIN_LINE_PT)
MINOR_LINE_GRID <- pt_to_grid_lwd(MINOR_LINE_PT)

# 主网格线和次网格线统一使用指定的浅灰色
GRID_MAJOR_COLOR <- "#EBEBEB"
GRID_MINOR_COLOR <- "#EBEBEB"

# 每个单元格中，X轴标签区域高度
# 200 - 136.515 - 18.238 = 45.247 mm
X_LABEL_ZONE_MM <- CELL_HEIGHT_MM -
  PANEL_HEIGHT_MM -
  TITLE_HEIGHT_MM

if (X_LABEL_ZONE_MM <= 0) {
  stop("页面高度不足以容纳指定的标题框和XY绘图区。")
}

# 绘图区在每个190 mm单元格中水平居中
PANEL_LEFT_IN_CELL_MM <- (
  CELL_WIDTH_MM - PANEL_WIDTH_MM
) / 2

# 标题框与绘图区左右边缘完全对齐
TITLE_LEFT_IN_CELL_MM <- PANEL_LEFT_IN_CELL_MM

# Y轴刻度数字和Y轴标题手工定位
Y_TICK_LABEL_X_OFFSET_MM <- 4
Y_AXIS_TITLE_X_OFFSET_MM <- 31

# X轴标签紧接绘图区底边
X_LABEL_ANCHOR_OFFSET_MM <- 0.10


# ============================================================
# 4. 组别、颜色、条纹设置
# ============================================================

group_levels <- c(
  "CTR",
  "MP1",
  "MP10",
  "MPP1",
  "MPP10"
)

# WRS2稳健单因素方差分析设置
TRIM_LEVEL <- 0.20
ALPHA_LEVEL <- 0.05
LINCON_P_ADJUST_METHOD <- "none"

# 保留此前确定的原始颜色
fill_colors <- c(
  CTR  = "#b6dcd5",
  MP1  = "#ffe0e2",
  MP10 = "#ffaab9",
  MPP1 = "#baa5d5",
  MPP10 = "#8c69b9"
)

# 保留此前确定的条纹角度
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

# X轴使用完整组别名称
group_labels <- c(
  CTR  = "Control",
  MP1  = "Microplastic-Low",
  MP10 = "Microplastic-High",
  MPP1 = "Plastisphere-Low",
  MPP10 = "Plastisphere-High"
)


# ============================================================
# 5. 指标参数及最终顺序
#
# 顺序依据(a)-(j)连续编号：
# 第一排：(a)-(e)
# 第二排：(f)-(j)
#
# 第一项标题缩短为Siphoning：
# 36.59 pt字体无法在89.277 mm标题框中完整容纳
# "Siphoning_efficiency"；完整指标名称和单位保留在Y轴。
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
# 6. 读取数据并检查分组
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
    "原始数据缺少以下列：",
    paste(missing_columns, collapse = ", ")
  )
}

dat_raw$group <- factor(
  trimws(as.character(dat_raw$group)),
  levels = group_levels
)

if (anyNA(dat_raw$group)) {
  stop(
    "group列中存在无法识别的组别。",
    "\n应为：CTR、MP1、MP10、MPP1、MPP10。"
  )
}

cat("\n========== 原始数据分组 ==========\n")
print(table(dat_raw$group))


# ============================================================
# 7. P值和坐标刻度格式化
# ============================================================

format_p_value_suffix <- function(p) {
  
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p < 0.001) {
    return("< 0.001")
  }
  
  # 避免显著P值四舍五入后显示为P = 0.050
  if (round(p, 3) >= 0.05) {
    return(sprintf("= %.4f", p))
  }
  
  sprintf("= %.3f", p)
}

# 用于CSV统计结果的完整文本。
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
# 8. 创建单个指标的panel grob及统计结果
#
# 这里只生成XY绘图区本身：
#   不生成ggplot默认X/Y轴文字；
#   不生成ggplot默认标题；
#   不生成ggplot默认外部边距。
#
# 轴标题、轴刻度、标题框在后续grid手工布局中绘制。
# ============================================================

make_panel <- function(variable, config) {
  
  # ----------------------------------------------------------
  # 8.1 原始数据
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
      
      # 只沿水平方向轻微错开
      # 纵坐标完全保留原始观测值
      Plot_x = as.numeric(first(group)) +
        c(-0.10, 0, 0.10)[Replicate_order]
    ) %>%
    ungroup()
  
  # 检查每组独立生物学重复数
  sample_check <- dat %>%
    count(group, name = "n")
  
  if (nrow(sample_check) != 5 ||
      any(sample_check$n != 3)) {
    
    stop(
      variable,
      "：要求五组各3个有效独立生物学重复。\n",
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
  # 8.3 WRS2 t1way总体检验
  # ----------------------------------------------------------
  # t1way检验五组总体差异，允许组间方差不齐。
  # 每组n = 3且tr = 0.20时，floor(n * tr) = 0，
  # 因此本数据不会实际截尾，检验在当前样本量下相当于
  # 基于原始均值的Welch型异方差单因素检验。
  
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
  # 8.4 WRS2 lincon全部两两比较
  # ----------------------------------------------------------
  # 五组共10组两两比较。
  # method = "none"：lincon返回未进行多重比较校正的原始P值。
  # 注意：lincon的置信区间仍为同时置信区间；本图仅依据原始P值标注。
  # 同时输出截尾均值差及其同时置信区间。
  
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
      "：WRS2::lincon()未返回预期的10组两两比较结果。"
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
  # 8.5 显著比较及括号排序
  # ----------------------------------------------------------
  # 仅根据未校正的原始P值筛选和标注。
  # 对水平方向不重叠的比较，允许共用同一括号高度，
  # 避免所有显著比较都逐条堆叠而不必要地拉高Y轴。
  
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
  
  # 区间分层：只有横向区间互不重叠的括号才共用同一高度。
  # 若前一括号终点与后一括号起点相同，也分配到不同层，
  # 避免两个括号端点在同一位置相接。
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
  # 8.6 计算数据顶部，并优先保留初始Y轴范围
  # ----------------------------------------------------------
  # data_top必须在计算base_y_upper和显著性括号位置之前定义。
  # 同时考虑Mean + SD误差棒上端及全部原始观测值，避免任何数据元素被裁切。
  data_top <- max(
    c(
      summary_data$Mean + summary_data$SD,
      dat$concentration
    ),
    na.rm = TRUE
  )

  if (!is.finite(data_top)) {
    stop(variable, "：无法计算有效的data_top，请检查原始数据。")
  }

  # 先以plot_config中最初设定的ymax作为Y轴上限。
  # 若原始数据（Mean + SD或单个原始值）本身超过该范围，才因数据本身上调。
  # 显著性括号优先使用data_top以上的空白区域；空间不足时逐步缩小该panel的
  # P值字号（18 pt至10 pt），而不是先扩大Y轴。
  # 只有10 pt仍无法容纳时，才按所需空间对Y轴上限作最小幅度扩展。

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

    # 在不改变Y轴上限的情况下，计算data_top以上可用的物理高度。
    available_height_mm <- (
      (base_y_upper - data_top) / base_y_upper
    ) * PANEL_HEIGHT_MM

    # 从默认18 pt开始逐级尝试，选择在初始Y轴范围内能够容纳的最大字号。
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
      # 原始ymax足以容纳全部标注：保持Y轴范围不变。
      y_upper <- base_y_upper
      p_font_pt <- selected_font
      p_text_height_mm <- selected_text_height
      p_bracket_step_mm <- selected_step

    } else {
      # 即使字号缩小到下限仍放不下时，才计算最小必要Y轴扩展量。
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
          "：显著性标注过多，即使P值字号缩小至",
          P_VALUE_MIN_FONT_PT,
          " pt仍无法合理容纳。"
        )
      }

      y_upper_needed <- data_top / (1 - reserved_fraction)

      # 只按较细的break_by/4步长上调，避免按整格break_by上调造成过度扩展。
      y_rounding_step <- config$break_by / 4
      y_upper <- max(
        base_y_upper,
        ceiling(y_upper_needed / y_rounding_step) * y_rounding_step
      )
    }

    # 根据最终Y轴范围，将括号的物理间距换算为数据坐标。
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

    # 没有显著性标注时，严格保留初始ymax（除非数据本身超出）。
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
  
  # 去掉与主刻度重合的次刻度
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
  # 8.7 绘制panel内部内容
  # ----------------------------------------------------------
  
  p <- ggplot() +
    
    # 柱子：Mean
    # 柱边框线宽1.788 pt
    
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
    
    # 误差棒：Mean ± SD
    
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
    
    # 原始散点：白色填充、黑色边框
    
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
    
    # 显著性括号：与主轴线一致，1.788 pt
    
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