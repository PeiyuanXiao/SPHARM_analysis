# viz3d_utils.R
# Plotly 3D helpers shared by align_svd.R and align_lin2024.R. Requires: plotly, jsonlite.

# Draw all scars as line segments with cone arrowheads.
# highlight_idx (optional): index of one scar to emphasise (used by Lin 2024 for the longest scar).
add_scars_3d <- function(fig, sx, sy, sz, ex, ey, ez,
                         highlight_idx = NULL) {
  for (i in seq_along(sx)) {
    is_hl <- !is.null(highlight_idx) && i == highlight_idx
    clr   <- if (is_hl) "pink" else "steelblue"
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

# Draw a single direction arrow: `direction` is normalised then scaled by `scale`.
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

# Draw a quad plane with arbitrary normal, centred at `center`.
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

# Draw a horizontal plane (parallel to XY) at height z0.
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

# Shared 3D scene layout (camera, axes, aspect).
make_scene <- function() {
  list(
    camera     = list(eye = list(x = 1.6, y = 1.6, z = 1.1)),
    xaxis      = list(title = "X", showgrid = TRUE, zeroline = TRUE),
    yaxis      = list(title = "Y", showgrid = TRUE, zeroline = TRUE),
    zaxis      = list(title = "Z", showgrid = TRUE, zeroline = TRUE),
    aspectmode = "data"
  )
}

# Per-panel layout; title_text accepts HTML bold tags.
panel_layout <- function(title_text) {
  list(
    title         = list(text = title_text, font = list(size = 13),
                         x = 0.5, xanchor = "center"),
    margin        = list(t = 50, b = 5, l = 5, r = 5),
    paper_bgcolor = "#f5f7fa",
    scene         = make_scene()
  )
}

# Serialise a plotly object to a JSON string for HTML embedding.
get_panel_json <- function(p) {
  built       <- plotly_build(p)
  data_json   <- toJSON(built$x$data,   auto_unbox = TRUE, null = "null", force = TRUE)
  layout_json <- toJSON(built$x$layout, auto_unbox = TRUE, null = "null", force = TRUE)
  paste0('{"data":', data_json, ',"layout":', layout_json, '}')
}
