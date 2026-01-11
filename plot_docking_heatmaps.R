args <- commandArgs(trailingOnly = TRUE)

## ==========================
## 用户可自定义参数（优先修改这里）
## ==========================

cfg <- list(
  ## 输入输出
  ## - input_csv: 对接结果 CSV 路径
  ## - output_root: 输出根目录（为空则默认 CSV 所在目录）
  input_csv_default = "/Users/dingtunan/生信分析/批量对接脚本/0109对接结果/results/docking_results.csv",
  output_root_default = NULL,

  ## 输出文件夹名称（会自动创建）
  output_dirs = list(
    style1 = "heatmap_style1",
    style2 = "heatmap_style2",
    style3 = "heatmap_style3"
  ),

  ## 输出文件名（不带扩展名）
  ## - base_name: 正常方向
  ## - swapped_name: xy 轴对调版本
  output_files = list(
    base_name = "docking_heatmap",
    swapped_name = "docking_heatmap_xy_swapped"
  ),

  ## 图像保存参数
  save = list(
    dpi = 300,
    bg = "white"
  ),

  ## 画布尺寸（单位 inch）
  ## - width = max(width_min, n_x * width_per_x)
  ## - height = max(height_min, n_y * height_per_y)
  size = list(
    width_min = 8,
    width_per_x = 1.0,
    height_min = 8,
    height_per_y = 0.55
  ),

  ## 标题与坐标轴（按你要求：纯英文）
  labels = list(
    title = "Docking Heatmap",
    x = "Protein",
    y = "Ligand",
    x_swapped = "Ligand",
    y_swapped = "Protein"
  ),

  ## 文字显示
  text = list(
    value_digits = 2,
    value_size = 3,
    value_color = "black"
  ),

  ## 主题与坐标轴样式
  theme = list(
    base_size = 11,
    title_hjust = 0.5,
    x_text_angle = 45,
    x_text_hjust = 1,
    x_text_vjust = 1,
    y_text_size = 9,
    remove_grid = TRUE
  ),

  ## 颜色（顶刊常用红-白-蓝发散）
  ## - low/mid/high: 对应负值/0/正值
  colors = list(
    low = "red",
    mid = "white",
    high = "blue",
    positive_na = "grey75",
    tile_border_12 = "grey80",
    tile_border_3 = "grey70",
    point_border = "grey20"
  ),

  ## 风格2/3 的“颜色范围”处理
  ## - 你的需求：范围取 min/max，然后上下取整；但颜色分配以 0 为中心
  scale = list(
    round_min = TRUE,
    round_max = TRUE
  ),

  ## 风格1（仅负值区间上色，正值灰）
  style1 = list(
    tile_linewidth = 0.4,
    use_negative_midpoint = TRUE
  ),

  ## 风格2（全矩阵上色 + 数值）
  style2 = list(
    tile_linewidth = 0.4
  ),

  ## 风格3（方格边框 + 圆点填充 + 数值，圆大小分层）
  style3 = list(
    tile_linewidth = 0.45,
    point_stroke = 0.3,
    ## 分层阈值（可自行改为更细/更粗）
    tier_breaks = c(-Inf, -10, -8, -6, -4, Inf),
    tier_labels = c("<= -10", "(-10, -8]", "(-8, -6]", "(-6, -4]", "> -4"),
    ## 圆点大小：最大圆直径≈格子边长（建议<=1，避免压边）
    max_diameter_to_tile = 1.00,
    ## 各层级相对大小（从最强负值到最弱/正值）
    tier_size_ratios = c(1.00, 0.85, 0.70, 0.55, 0.45)
  ),

  ## 是否生成 xy 轴对调版本
  output_swapped = TRUE
)

## 命令行参数优先生效：
## - 第1个参数：input_csv
## - 第2个参数：output_root
input_csv <- if (length(args) >= 1) args[[1]] else cfg$input_csv_default
output_root <- if (length(args) >= 2) args[[2]] else cfg$output_root_default
if (is.null(output_root) || !nzchar(output_root)) output_root <- dirname(input_csv)

required_pkgs <- c("ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "Missing R packages: %s. Please install them first.",
      paste(missing_pkgs, collapse = ", ")
    ),
    call. = FALSE
  )
}

## 创建文件夹（若不存在则递归创建），并返回规范化路径
safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

## 读取对接结果 CSV，并校验关键列与数值格式
read_results <- function(path) {
  if (!file.exists(path)) stop(sprintf("Input file not found: %s", path), call. = FALSE)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  required_cols <- c("Receptor", "Ligand", "Affinity")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in CSV: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  df$Affinity <- suppressWarnings(as.numeric(df$Affinity))
  if (anyNA(df$Affinity)) {
    bad_n <- sum(is.na(df$Affinity))
    stop(sprintf("Affinity has %d NA values after numeric conversion.", bad_n), call. = FALSE)
  }
  df
}

## 根据行列数量估算画布尺寸
compute_plot_size <- function(n_x, n_y, cfg) {
  width <- max(cfg$size$width_min, n_x * cfg$size$width_per_x)
  height <- max(cfg$size$height_min, n_y * cfg$size$height_per_y)
  list(width = width, height = height)
}

## 依据画布宽度与列数，估算单个格子边长对应的圆点半径（mm）
estimate_max_point_radius_mm <- function(plot_width_in, n_x, max_diameter_to_tile = 1.0) {
  if (!is.finite(plot_width_in) || !is.finite(n_x) || n_x <= 0) return(4)
  tile_side_mm <- (plot_width_in * 25.4) / n_x
  (tile_side_mm / 2) * max_diameter_to_tile
}

## 同时保存 PNG 与 PDF
save_plot_dual <- function(p, out_base, width, height, cfg) {
  ggplot2::ggsave(
    filename = paste0(out_base, ".png"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = cfg$save$dpi,
    bg = cfg$save$bg
  )
  ggplot2::ggsave(
    filename = paste0(out_base, ".pdf"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = cfg$save$bg
  )
}

## 将数值按“负值/正值分别归一化”为对称区间 [-1, 1]（0 永远在中心）
to_symmetric_fill <- function(x, lim) {
  if (length(lim) != 2L || any(!is.finite(lim))) return(rep(NA_real_, length(x)))
  minv <- lim[[1]]
  maxv <- lim[[2]]
  if (!is.finite(minv) || !is.finite(maxv) || minv == maxv) return(rep(0, length(x)))

  x_clamped <- pmin(pmax(x, minv), maxv)

  if (minv >= 0) {
    if (maxv == 0) return(rep(0, length(x_clamped)))
    return((x_clamped / maxv))
  }
  if (maxv <= 0) {
    if (minv == 0) return(rep(0, length(x_clamped)))
    return((x_clamped / abs(minv)))
  }

  out <- rep(NA_real_, length(x_clamped))
  neg_idx <- x_clamped <= 0
  out[neg_idx] <- x_clamped[neg_idx] / abs(minv)
  pos_idx <- x_clamped > 0
  out[pos_idx] <- x_clamped[pos_idx] / maxv
  out
}

## 风格1：只对负值上色（正值灰），色带范围仅覆盖负值区间
make_style1 <- function(df, swap_xy = FALSE, lim = NULL, neg_midpoint = NULL, cfg) {
  df$FillValue <- ifelse(df$Affinity > 0, NA_real_, df$Affinity)
  df$Label <- sprintf(paste0("%.", cfg$text$value_digits, "f"), df$Affinity)
  x_var <- if (swap_xy) "Ligand" else "Receptor"
  y_var <- if (swap_xy) "Receptor" else "Ligand"

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = FillValue),
      color = cfg$colors$tile_border_12,
      linewidth = cfg$style1$tile_linewidth
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label),
      size = cfg$text$value_size,
      color = cfg$text$value_color
    ) +
    ggplot2::scale_fill_gradient2(
      low = cfg$colors$low,
      mid = cfg$colors$mid,
      high = cfg$colors$high,
      midpoint = neg_midpoint,
      limits = lim,
      na.value = cfg$colors$positive_na,
      name = "Affinity"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = cfg$labels$title,
      x = if (swap_xy) cfg$labels$x_swapped else cfg$labels$x,
      y = if (swap_xy) cfg$labels$y_swapped else cfg$labels$y
    ) +
    ggplot2::theme_minimal(base_size = cfg$theme$base_size) +
    ggplot2::theme(
      panel.grid = if (cfg$theme$remove_grid) ggplot2::element_blank() else ggplot2::element_line(),
      axis.text.x = ggplot2::element_text(
        angle = cfg$theme$x_text_angle,
        hjust = cfg$theme$x_text_hjust,
        vjust = cfg$theme$x_text_vjust
      ),
      axis.text.y = ggplot2::element_text(size = cfg$theme$y_text_size),
      plot.title = ggplot2::element_text(hjust = cfg$theme$title_hjust)
    )
}

## 风格2：全数据红-白-蓝，范围为 [最小值, 最大值]，0 固定为白色中心
make_style2 <- function(df, swap_xy = FALSE, lim = NULL, cfg) {
  df$Label <- sprintf(paste0("%.", cfg$text$value_digits, "f"), df$Affinity)
  x_var <- if (swap_xy) "Ligand" else "Receptor"
  y_var <- if (swap_xy) "Receptor" else "Ligand"

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = AffinityFillSym),
      color = cfg$colors$tile_border_12,
      linewidth = cfg$style2$tile_linewidth
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label),
      size = cfg$text$value_size,
      color = cfg$text$value_color
    ) +
    ggplot2::scale_fill_gradient2(
      low = cfg$colors$low,
      mid = cfg$colors$mid,
      high = cfg$colors$high,
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      labels = c(sprintf("%g", lim[[1]]), "0", sprintf("%g", lim[[2]])),
      name = "Affinity"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = cfg$labels$title,
      x = if (swap_xy) cfg$labels$x_swapped else cfg$labels$x,
      y = if (swap_xy) cfg$labels$y_swapped else cfg$labels$y
    ) +
    ggplot2::theme_minimal(base_size = cfg$theme$base_size) +
    ggplot2::theme(
      panel.grid = if (cfg$theme$remove_grid) ggplot2::element_blank() else ggplot2::element_line(),
      axis.text.x = ggplot2::element_text(
        angle = cfg$theme$x_text_angle,
        hjust = cfg$theme$x_text_hjust,
        vjust = cfg$theme$x_text_vjust
      ),
      axis.text.y = ggplot2::element_text(size = cfg$theme$y_text_size),
      plot.title = ggplot2::element_text(hjust = cfg$theme$title_hjust)
    )
}

## 将亲和度分层，用于风格3的圆点大小
make_size_tier <- function(x, cfg) {
  cut(
    x,
    breaks = cfg$style3$tier_breaks,
    labels = cfg$style3$tier_labels,
    right = TRUE,
    include.lowest = TRUE
  )
}

## 风格3：方格边框 + 圆点填充颜色 + 数值标注，圆大小按分层变化
make_style3 <- function(df, swap_xy = FALSE, lim = NULL, max_radius_mm = NULL, cfg) {
  df$SizeTier <- make_size_tier(df$Affinity, cfg)
  size_levels <- levels(df$SizeTier)
  if (is.null(max_radius_mm) || !is.finite(max_radius_mm)) max_radius_mm <- 6
  size_values <- max_radius_mm * cfg$style3$tier_size_ratios
  names(size_values) <- size_levels
  df$Label <- sprintf(paste0("%.", cfg$text$value_digits, "f"), df$Affinity)
  x_var <- if (swap_xy) "Ligand" else "Receptor"
  y_var <- if (swap_xy) "Receptor" else "Ligand"

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_tile(fill = NA, color = cfg$colors$tile_border_3, linewidth = cfg$style3$tile_linewidth) +
    ggplot2::geom_point(
      ggplot2::aes(fill = AffinityFillSym, size = SizeTier),
      shape = 21,
      color = cfg$colors$point_border,
      stroke = cfg$style3$point_stroke
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label),
      size = cfg$text$value_size,
      color = cfg$text$value_color
    ) +
    ggplot2::scale_fill_gradient2(
      low = cfg$colors$low,
      mid = cfg$colors$mid,
      high = cfg$colors$high,
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      labels = c(sprintf("%g", lim[[1]]), "0", sprintf("%g", lim[[2]])),
      name = "Affinity"
    ) +
    ggplot2::scale_size_manual(values = size_values, name = "Affinity Tier") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = cfg$labels$title,
      x = if (swap_xy) cfg$labels$x_swapped else cfg$labels$x,
      y = if (swap_xy) cfg$labels$y_swapped else cfg$labels$y
    ) +
    ggplot2::theme_minimal(base_size = cfg$theme$base_size) +
    ggplot2::theme(
      panel.grid = if (cfg$theme$remove_grid) ggplot2::element_blank() else ggplot2::element_line(),
      axis.text.x = ggplot2::element_text(
        angle = cfg$theme$x_text_angle,
        hjust = cfg$theme$x_text_hjust,
        vjust = cfg$theme$x_text_vjust
      ),
      axis.text.y = ggplot2::element_text(size = cfg$theme$y_text_size),
      plot.title = ggplot2::element_text(hjust = cfg$theme$title_hjust)
    )
}

df <- read_results(input_csv)

receptor_levels <- unique(df$Receptor)
ligand_levels <- unique(df$Ligand)
df$Receptor <- factor(df$Receptor, levels = receptor_levels)
df$Ligand <- factor(df$Ligand, levels = rev(ligand_levels))

lim_all <- c(min(df$Affinity), max(df$Affinity))
lim_all_rounded <- lim_all
if (isTRUE(cfg$scale$round_min)) lim_all_rounded[[1]] <- floor(lim_all[[1]])
if (isTRUE(cfg$scale$round_max)) lim_all_rounded[[2]] <- ceiling(lim_all[[2]])
df$AffinityFillSym <- to_symmetric_fill(df$Affinity, lim = lim_all_rounded)
neg_vals <- df$Affinity[df$Affinity < 0]
if (length(neg_vals) == 0) stop("No negative Affinity values found; style1 needs negatives.", call. = FALSE)
lim_neg <- c(min(neg_vals), max(neg_vals))
neg_midpoint <- if (isTRUE(cfg$style1$use_negative_midpoint)) mean(lim_neg) else 0

style1_dir <- safe_mkdir(file.path(output_root, cfg$output_dirs$style1))
style2_dir <- safe_mkdir(file.path(output_root, cfg$output_dirs$style2))
style3_dir <- safe_mkdir(file.path(output_root, cfg$output_dirs$style3))

plot_size <- compute_plot_size(length(receptor_levels), length(ligand_levels), cfg)
max_radius_normal <- estimate_max_point_radius_mm(
  plot_width_in = plot_size$width,
  n_x = length(receptor_levels),
  max_diameter_to_tile = cfg$style3$max_diameter_to_tile
)
max_radius_swapped <- estimate_max_point_radius_mm(
  plot_width_in = plot_size$width,
  n_x = length(ligand_levels),
  max_diameter_to_tile = cfg$style3$max_diameter_to_tile
)

p1 <- make_style1(df, swap_xy = FALSE, lim = lim_neg, neg_midpoint = neg_midpoint, cfg = cfg)
save_plot_dual(
  p1,
  file.path(style1_dir, cfg$output_files$base_name),
  plot_size$width,
  plot_size$height,
  cfg
)
if (isTRUE(cfg$output_swapped)) {
  p1_t <- make_style1(df, swap_xy = TRUE, lim = lim_neg, neg_midpoint = neg_midpoint, cfg = cfg)
  save_plot_dual(
    p1_t,
    file.path(style1_dir, cfg$output_files$swapped_name),
    plot_size$width,
    plot_size$height,
    cfg
  )
}

p2 <- make_style2(df, swap_xy = FALSE, lim = lim_all_rounded, cfg = cfg)
save_plot_dual(
  p2,
  file.path(style2_dir, cfg$output_files$base_name),
  plot_size$width,
  plot_size$height,
  cfg
)
if (isTRUE(cfg$output_swapped)) {
  p2_t <- make_style2(df, swap_xy = TRUE, lim = lim_all_rounded, cfg = cfg)
  save_plot_dual(
    p2_t,
    file.path(style2_dir, cfg$output_files$swapped_name),
    plot_size$width,
    plot_size$height,
    cfg
  )
}

p3 <- make_style3(
  df,
  swap_xy = FALSE,
  lim = lim_all_rounded,
  max_radius_mm = max_radius_normal,
  cfg = cfg
)
save_plot_dual(
  p3,
  file.path(style3_dir, cfg$output_files$base_name),
  plot_size$width,
  plot_size$height,
  cfg
)
if (isTRUE(cfg$output_swapped)) {
  p3_t <- make_style3(
    df,
    swap_xy = TRUE,
    lim = lim_all_rounded,
    max_radius_mm = max_radius_swapped,
    cfg = cfg
  )
  save_plot_dual(
    p3_t,
    file.path(style3_dir, cfg$output_files$swapped_name),
    plot_size$width,
    plot_size$height,
    cfg
  )
}

message("Done. Output folders:")
message("- ", style1_dir)
message("- ", style2_dir)
message("- ", style3_dir)
