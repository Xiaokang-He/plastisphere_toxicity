# Repository copy of the analysis script from the authors' working files.
# Repository-only adaptation: local absolute paths were replaced with paths
# relative to the project root. The statistical and plotting operations were
# otherwise retained. Review the associated input-data requirements in README.md.

# ============================================================
# 肠道微生物16S α多样性柱状图：手工物理布局紧凑标注版
#
# 页面：240 × 80 mm
# 排列：3列 × 1行；每个单元格80 × 80 mm
#
# 每张图的严格尺寸：
#   XY绘图区：35.7108 × 54.606 mm（按第二段代码布局等比例缩小）
#   标题框：  35.7108 × 7.2952 mm
#
# 字体（Times New Roman）：
#   标题：14.636 pt，加粗
#   Y轴刻度：11.36 pt
#   Y轴标题：11.36 pt，加粗
#   X轴刻度：11.36 pt，倾斜30°
#   显著性P值：7.2 pt；仅P字符斜体，其余数字和符号为常规体
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
#   不输出字体轮廓化文件
#
# 本版与生理生化指标图使用相同的物理布局、字体、线宽、
# 网格、柱形纹理、散点、误差棒、显著性标注和统计流程。
# 仅将展示指标替换为Shannon、Pielou_e和ACE三个α多样性指标。
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
input_file <- file.path(input_dir, "α多样性数据.csv")
output_dir <- file.path(results_dir, "fig5_gut_16s_alpha_diversity")

if (!file.exists(input_file)) {
  stop(
    "找不到原始数据文件：\n",
    input_file,
    "\n请确认文件名和路径。"
  )
}

# Results are written to results/fig5_gut_16s_alpha_diversity.

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop(
    "无法创建输出文件夹：\n", output_dir,
    "\n请检查输出目录是否存在、是否有写入权限。"
  )
}

output_dir <- normalizePath(
  output_dir,
  winslash = "/",
  mustWork = TRUE
)


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

# 第一段原始代码的合并图导出尺寸为 width = 240 mm、height = 80 mm；
# 三张图各占一个80 × 80 mm单元格。单张图也按80 × 80 mm导出。
PAGE_WIDTH_MM  <- 240
PAGE_HEIGHT_MM <- 80

N_COLS <- 3
N_ROWS <- 1

CELL_WIDTH_MM  <- PAGE_WIDTH_MM / N_COLS
CELL_HEIGHT_MM <- PAGE_HEIGHT_MM / N_ROWS

# 第二段代码的单元格为190 × 200 mm，而第一段合并图每格为80 × 80 mm。
# 为保持第二段手工布局的相对比例，将面板、标题框、文字和间距按0.4等比例缩小。
# 线宽仍严格保留第二段代码的1.788 pt和0.886 pt，不随页面缩小。
LAYOUT_SCALE <- 0.4

PANEL_WIDTH_MM  <- 89.277 * LAYOUT_SCALE
PANEL_HEIGHT_MM <- 136.515 * LAYOUT_SCALE

TITLE_WIDTH_MM  <- 89.277 * LAYOUT_SCALE
TITLE_HEIGHT_MM <- 18.238 * LAYOUT_SCALE

# 文字字号按页面比例缩小，保证在80 × 80 mm单图中不溢出。
TITLE_FONT_PT <- 36.59 * LAYOUT_SCALE
AXIS_TEXT_PT  <- 28.4 * LAYOUT_SCALE
AXIS_TITLE_PT <- 28.4 * LAYOUT_SCALE
P_VALUE_FONT_PT <- 18 * LAYOUT_SCALE

# P值与括号的间距也按0.4等比例缩小。
P_VALUE_CLEARANCE_MM <- 1.5 * LAYOUT_SCALE

# 按缩放后的P值字号计算Times New Roman文字的预留高度。
# 原设置采用1.20倍字号高度，过度放大了括号层间距和Y轴上限；
# 这里改为0.90倍字号高度。绘制时仍通过grid实际测量文字高度，
# 因此该值主要用于预估Y轴顶部空间。
P_VALUE_TEXT_HEIGHT_MM <- P_VALUE_FONT_PT * 25.4 / 72 * 0.90

# 相邻显著性括号横线之间的物理间距：
# 预留文字高度 + 文字/横线间距 + 等比例缩放后的额外空隙。
P_BRACKET_EXTRA_GAP_MM <- 0.5 * LAYOUT_SCALE
P_BRACKET_STEP_MM <- P_VALUE_TEXT_HEIGHT_MM +
  P_VALUE_CLEARANCE_MM + P_BRACKET_EXTRA_GAP_MM

# 最低一条显著性括号与数据最高点之间的物理间距
P_FIRST_BRACKET_GAP_MM <- 2.5 * LAYOUT_SCALE

# 最上方P值文字与绘图区顶部之间的预留空间
P_VALUE_TOP_MARGIN_MM <- 2.5 * LAYOUT_SCALE

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

# 每个80 mm高单元格中，X轴标签区域高度自动计算。
X_LABEL_ZONE_MM <- CELL_HEIGHT_MM -
  PANEL_HEIGHT_MM -
  TITLE_HEIGHT_MM

if (X_LABEL_ZONE_MM <= 0) {
  stop("页面高度不足以容纳指定的标题框和XY绘图区。")
}

# 绘图区在每个80 mm单元格中水平居中
PANEL_LEFT_IN_CELL_MM <- (
  CELL_WIDTH_MM - PANEL_WIDTH_MM
) / 2

# 标题框与绘图区左右边缘完全对齐
TITLE_LEFT_IN_CELL_MM <- PANEL_LEFT_IN_CELL_MM

# Y轴刻度数字和Y轴标题手工定位
Y_TICK_LABEL_X_OFFSET_MM <- 4 * LAYOUT_SCALE
Y_AXIS_TITLE_X_OFFSET_MM <- 31 * LAYOUT_SCALE

# X轴标签紧接绘图区底边
X_LABEL_ANCHOR_OFFSET_MM <- 0.10 * LAYOUT_SCALE


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
# 5. α多样性指标参数及最终顺序
#
# 顺序与原始合并图一致：
#   (a) Shannon
#   (b) Pielou_e
#   (c) ACE
#
# Y轴基础范围沿用第一段α多样性代码中的设置；
# 若显著性标注需要额外空间，程序会像第二段代码一样自动扩展Y轴上限。
# ============================================================

plot_config <- list(

  shannon = list(
    tag = "(a)",
    title = "Shannon",
    y_label = "Shannon index of OTU level",
    ymax = 3.25,
    break_by = 0.5,
    digits = 1,
    bar_width = 0.8
  ),

  Pielou_e = list(
    tag = "(b)",
    title = "Pielou_e",
    y_label = "Pielou_e index of OTU level",
    ymax = 0.65,
    break_by = 0.1,
    digits = 1,
    bar_width = 0.8
  ),

  ace = list(
    tag = "(c)",
    title = "ACE",
    y_label = "ACE index of OTU level",
    ymax = 195,
    break_by = 30,
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
  "sample",
  "group",
  plot_order
)

missing_columns <- setdiff(
  required_columns,
  names(dat_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "原始α多样性数据缺少以下列：",
    paste(missing_columns, collapse = ", ")
  )
}

# 保留样本编号，并检查每行是否对应唯一的独立样本。
dat_raw$sample <- trimws(as.character(dat_raw$sample))

if (anyNA(dat_raw$sample) || any(!nzchar(dat_raw$sample))) {
  stop("sample列存在空值，请检查原始α多样性数据。")
}

if (anyDuplicated(dat_raw$sample)) {
  stop(
    "sample列存在重复样本编号：",
    paste(unique(dat_raw$sample[duplicated(dat_raw$sample)]),
          collapse = ", ")
  )
}

# 第一段代码中的原始分组编码与第二段代码使用的组别名称对应关系：
# G1 = MPP10；G2 = MP10；G3 = MPP1；G4 = CTR；G5 = MP1。
# 同时兼容数据文件中已经写为CTR、MP1、MP10、MPP1、MPP10的情况。
group_map <- c(
  G1 = "MPP10",
  G2 = "MP10",
  G3 = "MPP1",
  G4 = "CTR",
  G5 = "MP1"
)

group_original <- trimws(as.character(dat_raw$group))

dat_raw$group <- ifelse(
  group_original %in% names(group_map),
  unname(group_map[group_original]),
  group_original
)

dat_raw$group <- factor(
  trimws(as.character(dat_raw$group)),
  levels = group_levels
)

if (anyNA(dat_raw$group)) {
  stop(
    "group列中存在无法识别的组别。",
    "\\n原始组别应为G1、G2、G3、G4、G5，",
    "或已经转换为CTR、MP1、MP10、MPP1、MPP10。"
  )
}

# 三个α多样性指标均要求五组各有3个独立生物学重复；
# 每个指标均执行相同的t1way总体检验和lincon全部10组两两比较，
# 并按未校正原始P < 0.05标注实际P值。
# make_panel()中还会对每个指标的有效数据再次逐组检查。
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
  # 8.6 自动计算Y轴上限和括号高度
  # ----------------------------------------------------------
  
  data_top <- max(
    c(
      summary_data$Mean + summary_data$SD,
      dat$concentration
    ),
    na.rm = TRUE
  )
  
  if (nrow(sig) > 0) {
    
    # 根据实际需要的括号层数，而不是显著比较总数，
    # 计算从data_top到最上方P值文字所需的物理高度。
    n_lanes <- max(sig$lane)
    
    reserved_height_mm <-
      P_FIRST_BRACKET_GAP_MM +
      (n_lanes - 1) * P_BRACKET_STEP_MM +
      P_VALUE_CLEARANCE_MM +
      P_VALUE_TEXT_HEIGHT_MM +
      P_VALUE_TOP_MARGIN_MM
    
    reserved_fraction <- reserved_height_mm /
      PANEL_HEIGHT_MM
    
    if (reserved_fraction >= 0.80) {
      stop(
        variable,
        "：显著性标注数量过多，预留高度超过绘图区可用高度。"
      )
    }
    
    y_upper_needed <- data_top /
      (1 - reserved_fraction)
    
    y_upper <- max(
      config$ymax,
      ceiling(
        y_upper_needed / config$break_by
      ) * config$break_by
    )
    
    # 将物理mm间距换算成当前Y轴的数据坐标。
    first_gap_data <- (
      P_FIRST_BRACKET_GAP_MM /
        PANEL_HEIGHT_MM
    ) * y_upper
    
    bracket_step_data <- (
      P_BRACKET_STEP_MM /
        PANEL_HEIGHT_MM
    ) * y_upper
    
    label_center_offset_data <- (
      P_VALUE_CLEARANCE_MM +
        P_VALUE_TEXT_HEIGHT_MM / 2
    ) / PANEL_HEIGHT_MM * y_upper
    
    sig$y_position <- data_top +
      first_gap_data +
      (sig$lane - 1) * bracket_step_data
    
    # label_y为文字中心位置；绘制时再按实际字体高度定位。
    sig$label_y <- sig$y_position +
      label_center_offset_data
    
  } else {
    
    sig$y_position <- numeric(0)
    sig$label_y <- numeric(0)
    
    y_upper <- max(
      config$ymax,
      ceiling(
        data_top * 1.05 / config$break_by
      ) * config$break_by
    )
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
      
      # 主网格线：1.788 pt
      panel.grid.major.x = element_line(
        colour = GRID_MAJOR_COLOR,
        linewidth = MAIN_LINE_MM
      ),
      
      panel.grid.major.y = element_line(
        colour = GRID_MAJOR_COLOR,
        linewidth = MAIN_LINE_MM
      ),
      
      # 次网格线：0.886 pt
      panel.grid.minor.x = element_line(
        colour = GRID_MINOR_COLOR,
        linewidth = MINOR_LINE_MM
      ),
      
      panel.grid.minor.y = element_line(
        colour = GRID_MINOR_COLOR,
        linewidth = MINOR_LINE_MM
      ),
      
      # XY绘图区边框：1.788 pt
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