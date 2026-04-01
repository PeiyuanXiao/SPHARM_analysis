# ==============================================================================
# viz3d_utils.R
# Plotly 三维可视化公共函数
# 被以下脚本引用：
#   - 01_alignment/align_svd.R
#   - 01_alignment/align_lin2024.R
#
# 依赖包：plotly, jsonlite
# ==============================================================================


# ------------------------------------------------------------------------------
# add_scars_3d()
# 在 Plotly 3D 图中绘制所有刮痕（线段 + 箭头锥体）
#
# 参数：
#   fig          : plotly 图形对象
#   sx,sy,sz     : 刮痕起点坐标向量
#   ex,ey,ez     : 刮痕终点坐标向量
#   highlight_idx: （可选）需要高亮的刮痕索引，默认 NULL 即全部统一样式
#                  Lin 2024 方法中用于高亮最长刮痕
#
# 说明：
#   - 不传 highlight_idx → 全部刮痕统一显示为 steelblue（SVD 方法）
#   - 传入 highlight_idx → 该刮痕显示为粉色加粗（Lin 2024 方法）
# ------------------------------------------------------------------------------
add_scars_3d <- function(fig, sx, sy, sz, ex, ey, ez,
                         highlight_idx = NULL) {
  for (i in seq_along(sx)) {
    is_hl <- !is.null(highlight_idx) && i == highlight_idx
    clr   <- if (is_hl) "pink"  else "steelblue"
    lwd   <- if (is_hl) 6      else 2
    csz   <- if (is_hl) 0.08   else 0.05
    
    dx <- ex[i] - sx[i]
    dy <- ey[i] - sy[i]
    dz <- ez[i] - sz[i]
    
    fig <- fig %>% add_trace(
      type = "scatter3d", mode = "lines",
      x = c(sx[i], ex[i]), y = c(sy[i], ey[i]), z = c(sz[i], ez[i]),
      line = list(color = clr, width = lwd),
      showlegend = FALSE
    )
    fig <- fig %>% add_trace(
      type       = "cone",
      x = list(ex[i]), y = list(ey[i]), z = list(ez[i]),
      u = list(dx),    v = list(dy),    w = list(dz),
      sizemode   = "scaled", sizeref = csz,
      colorscale = list(list(0, clr), list(1, clr)),
      cmin = 0, cmax = 1,
      showscale = FALSE, showlegend = FALSE,
      anchor    = "tip",
      lighting  = list(ambient = 0.9, diffuse = 0.5)
    )
  }
  fig
}


# ------------------------------------------------------------------------------
# add_arrow_3d()
# 在 Plotly 3D 图中绘制一个带箭头的方向指示器
#
# 参数：
#   fig       : plotly 图形对象
#   origin    : 箭头起点，长度为3的数值向量
#   direction : 箭头方向向量（会被归一化后乘以 scale）
#   scale     : 箭头长度缩放系数
#   color     : 颜色字符串，默认 "red"
# ------------------------------------------------------------------------------
add_arrow_3d <- function(fig, origin, direction, scale, color = "red") {
  d   <- direction / sqrt(sum(direction^2)) * scale
  tip <- origin + d
  
  fig <- fig %>% add_trace(
    type = "scatter3d", mode = "lines",
    x = c(origin[1], tip[1]),
    y = c(origin[2], tip[2]),
    z = c(origin[3], tip[3]),
    line = list(color = color, width = 6),
    showlegend = FALSE
  )
  fig %>% add_trace(
    type       = "cone",
    x = list(tip[1]), y = list(tip[2]), z = list(tip[3]),
    u = list(d[1]),   v = list(d[2]),   w = list(d[3]),
    sizemode   = "scaled", sizeref = 0.2,
    colorscale = list(list(0, color), list(1, color)),
    cmin = 0, cmax = 1,
    showscale = FALSE, showlegend = FALSE,
    anchor    = "tail",
    lighting  = list(ambient = 0.9, diffuse = 0.5)
  )
}


# ------------------------------------------------------------------------------
# add_tilted_plane_3d()
# 在 Plotly 3D 图中绘制一个任意法线方向的平面（四边形网格）
#
# 参数：
#   fig       : plotly 图形对象
#   center    : 平面中心点，长度为3的数值向量
#   normal    : 平面法向量（会被归一化）
#   half_size : 平面半边长（控制显示大小）
# ------------------------------------------------------------------------------
add_tilted_plane_3d <- function(fig, center, normal, half_size) {
  n   <- normal / sqrt(sum(normal^2))
  ref <- if (abs(n[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
  u   <- ref - sum(ref * n) * n
  u   <- u / sqrt(sum(u^2))
  v   <- c(
    n[2] * u[3] - n[3] * u[2],
    n[3] * u[1] - n[1] * u[3],
    n[1] * u[2] - n[2] * u[1]
  )
  corners <- rbind(
    center + half_size * u + half_size * v,
    center - half_size * u + half_size * v,
    center - half_size * u - half_size * v,
    center + half_size * u - half_size * v
  )
  fig %>% add_trace(
    type = "mesh3d",
    x = corners[, 1], y = corners[, 2], z = corners[, 3],
    i = c(0, 0), j = c(1, 2), k = c(2, 3),
    intensity  = c(0, 0, 0, 0),
    colorscale = list(list(0, "lightgray"), list(1, "lightgray")),
    cmin = 0, cmax = 1,
    opacity = 0.35, flatshading = TRUE,
    showscale = FALSE, showlegend = FALSE
  )
}


# ------------------------------------------------------------------------------
# add_plane_3d()
# 在 Plotly 3D 图中绘制一个水平平面（平行于 XY 平面）
#
# 参数：
#   fig       : plotly 图形对象
#   cx, cy    : 平面中心的 X、Y 坐标
#   z0        : 平面所在的 Z 高度
#   half_size : 平面半边长
# ------------------------------------------------------------------------------
add_plane_3d <- function(fig, cx, cy, z0, half_size) {
  h  <- half_size
  xs <- c(cx - h, cx + h, cx + h, cx - h)
  ys <- c(cy - h, cy - h, cy + h, cy + h)
  zs <- rep(z0, 4)
  fig %>% add_trace(
    type = "mesh3d",
    x = xs, y = ys, z = zs,
    i = c(0, 0), j = c(1, 2), k = c(2, 3),
    intensity  = c(0, 0, 0, 0),
    colorscale = list(list(0, "lightgray"), list(1, "lightgray")),
    cmin = 0, cmax = 1,
    opacity = 0.35, flatshading = TRUE,
    showscale = FALSE, showlegend = FALSE
  )
}


# ------------------------------------------------------------------------------
# make_scene()
# 返回 Plotly 3D 场景的通用布局参数（相机角度、坐标轴设置等）
# ------------------------------------------------------------------------------
make_scene <- function() {
  list(
    camera     = list(eye = list(x = 1.6, y = 1.6, z = 1.1)),
    xaxis      = list(title = "X", showgrid = TRUE, zeroline = TRUE),
    yaxis      = list(title = "Y", showgrid = TRUE, zeroline = TRUE),
    zaxis      = list(title = "Z", showgrid = TRUE, zeroline = TRUE),
    aspectmode = "data"
  )
}


# ------------------------------------------------------------------------------
# panel_layout()
# 返回单个 Plotly 面板的 layout 参数（标题、边距、背景色）
#
# 参数：
#   title_text : 面板标题字符串（支持 HTML 加粗标签）
# ------------------------------------------------------------------------------
panel_layout <- function(title_text) {
  list(
    title         = list(text = title_text, font = list(size = 13),
                         x = 0.5, xanchor = "center"),
    margin        = list(t = 50, b = 5, l = 5, r = 5),
    paper_bgcolor = "#f5f7fa",
    scene         = make_scene()
  )
}


# ------------------------------------------------------------------------------
# get_panel_json()
# 将 Plotly 图形对象序列化为 JSON 字符串，用于嵌入 HTML 导出
#
# 参数：
#   p : plotly 图形对象
#
# 返回：
#   包含 data 和 layout 的 JSON 字符串
# ------------------------------------------------------------------------------
get_panel_json <- function(p) {
  built       <- plotly_build(p)
  data_json   <- toJSON(built$x$data,   auto_unbox = TRUE, null = "null", force = TRUE)
  layout_json <- toJSON(built$x$layout, auto_unbox = TRUE, null = "null", force = TRUE)
  paste0('{"data":', data_json, ',"layout":', layout_json, '}')
}