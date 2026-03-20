library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Rfast2)
library(Directional)
library(plotly)
library(htmltools)
library(jsonlite)
library(listviewer)
library(ks) 
library(rgl) 

# === Load Data ===
raw_data <- read_excel("analysis/data/raw_data/Scar_orientation_data.xlsx")
# ==============================================================================
# === Alignment Function (SVD normal vector) ===
align_svd <- function(df_group) {
  
  # === Step 1: Rotation — Align SVD normal to the Z-axis ===
  # Use the 3rd PC from the SVD of the scar vectors as the normal vector
  
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  
  if (sum(valid) >= 3) {
    U      <- cbind(dx[valid]/len[valid],
                    dy[valid]/len[valid],
                    dz[valid]/len[valid])
    normal <- svd(U)$v[, 3]
  }
  
  # Ensure the Z component of the normal vector is positive
  normal <- normal / sqrt(sum(normal^2))
  if (normal[3] < 0) normal <- -normal
  
  # Rotation matrix: align the normal vector with the Z-axis
  get_rot_matrix <- function(a, b) {
    a <- a / sqrt(sum(a^2))
    b <- b / sqrt(sum(b^2))
    cos_theta <- sum(a * b)
    if (cos_theta < -1 + 1e-10) {
      perp <- if (abs(a[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
      v <- perp - sum(perp * a) * a
      v <- v / sqrt(sum(v^2))
      return(2 * outer(v, v) - diag(3))
    }
    if (cos_theta > 1 - 1e-10) return(diag(3))
    v <- c(a[2]*b[3] - a[3]*b[2],
           a[3]*b[1] - a[1]*b[3],
           a[1]*b[2] - a[2]*b[1])
    v_skew <- matrix(c( 0,    -v[3],  v[2],
                        v[3],  0,    -v[1],
                        -v[2],  v[1],  0   ), 3, 3, byrow = TRUE)
    diag(3) + v_skew + v_skew %*% v_skew * ((1 - cos_theta) / sum(v^2))
  }
  
  R1 <- get_rot_matrix(normal, c(0, 0, 1))
  
  # Rotate all vectors
  S <- as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  
  df_group$s_x <- S[,1]; df_group$s_y <- S[,2]; df_group$s_z <- S[,3]
  df_group$e_x <- E[,1]; df_group$e_y <- E[,2]; df_group$e_z <- E[,3]
  df_group$d_x <- D[,1]; df_group$d_y <- D[,2]; df_group$d_z <- D[,3]
  
  # === Step 2: Translation — Center the rotated plane at the origin ===
  center         <- as.numeric(df_group[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  center_rotated <- as.numeric(R1 %*% center)
  
  df_group <- df_group %>%
    mutate(
      s_x = s_x - center_rotated[1],
      s_y = s_y - center_rotated[2],
      s_z = s_z - center_rotated[3],
      e_x = e_x - center_rotated[1],
      e_y = e_y - center_rotated[2],
      e_z = e_z - center_rotated[3]
    )
  
  return(df_group)
}

# --- Batch Execution ---
aligned_data <- raw_data %>%
  group_by(ID) %>%
  group_modify(~ align_svd(.x)) %>%
  ungroup()

# --- Verify the coordinates of the fitted plane normal ---
aligned_data %>%
  mutate(
    dx  = e_x - s_x, dy = e_y - s_y, dz = e_z - s_z,
    len = sqrt(dx^2 + dy^2 + dz^2),
    uz  = dz / len
  ) %>%
  filter(len > 1e-10) %>%
  group_by(ID, Typology) %>%
  summarise(
    N        = n(),
    uz_mean  = round(mean(uz),  3),
    uz_sd    = round(sd(uz),    3),
    .groups  = "drop"
  ) %>%
  print()

# ==============================================================================
# --- Model alignment visualization ---
# === Step 1: Plot features ===
get_rot_matrix <- function(a, b) {
  a <- a / sqrt(sum(a^2)); b <- b / sqrt(sum(b^2))
  cos_theta <- sum(a * b)
  if (cos_theta < -1 + 1e-10) {
    perp <- if (abs(a[1]) < 0.9) c(1,0,0) else c(0,1,0)
    v <- perp - sum(perp*a)*a; v <- v / sqrt(sum(v^2))
    return(2*outer(v,v) - diag(3))
  }
  if (cos_theta > 1 - 1e-10) return(diag(3))
  v <- c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
  v_skew <- matrix(c(0,-v[3],v[2],v[3],0,-v[1],-v[2],v[1],0), 3,3,byrow=TRUE)
  diag(3) + v_skew + v_skew %*% v_skew * ((1-cos_theta)/sum(v^2))
}

add_scars_3d <- function(fig, sx, sy, sz, ex, ey, ez,
                         lengths, highlight_idx = NULL) {
  for (i in seq_along(sx)) {
    is_hl <- !is.null(highlight_idx) && i == highlight_idx
    clr   <- if (is_hl) "pink" else "steelblue"
    lwd   <- if (is_hl) 6 else 2
    csz   <- if (is_hl) 0.08 else 0.05
    dx <- ex[i] - sx[i]; dy <- ey[i] - sy[i]; dz <- ez[i] - sz[i]
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

add_arrow_3d <- function(fig, origin, direction, scale, color = "red") {
  d   <- direction / sqrt(sum(direction^2)) * scale
  tip <- origin + d
  fig <- fig %>% add_trace(
    type = "scatter3d", mode = "lines",
    x = c(origin[1], tip[1]), y = c(origin[2], tip[2]), z = c(origin[3], tip[3]),
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

add_tilted_plane_3d <- function(fig, center, normal, half_size) {
  n   <- normal / sqrt(sum(normal^2))
  ref <- if (abs(n[1]) < 0.9) c(1,0,0) else c(0,1,0)
  u   <- ref - sum(ref*n)*n; u <- u / sqrt(sum(u^2))
  v   <- c(n[2]*u[3]-n[3]*u[2], n[3]*u[1]-n[1]*u[3], n[1]*u[2]-n[2]*u[1])
  corners <- rbind(
    center + half_size*u + half_size*v,
    center - half_size*u + half_size*v,
    center - half_size*u - half_size*v,
    center + half_size*u - half_size*v
  )
  fig %>% add_trace(
    type = "mesh3d",
    x = corners[,1], y = corners[,2], z = corners[,3],
    i = c(0,0), j = c(1,2), k = c(2,3),
    intensity  = c(0,0,0,0),
    colorscale = list(list(0,"lightgray"), list(1,"lightgray")),
    cmin = 0, cmax = 1,
    opacity = 0.35, flatshading = TRUE,
    showscale = FALSE, showlegend = FALSE
  )
}

add_plane_3d <- function(fig, cx, cy, z0, half_size) {
  h  <- half_size
  xs <- c(cx-h, cx+h, cx+h, cx-h)
  ys <- c(cy-h, cy-h, cy+h, cy+h)
  zs <- rep(z0, 4)
  fig %>% add_trace(
    type = "mesh3d",
    x = xs, y = ys, z = zs,
    i = c(0,0), j = c(1,2), k = c(2,3),
    intensity  = c(0,0,0,0),
    colorscale = list(list(0,"lightgray"), list(1,"lightgray")),
    cmin = 0, cmax = 1,
    opacity = 0.35, flatshading = TRUE,
    showscale = FALSE, showlegend = FALSE
  )
}

# === Step 2: Scene layout ===
make_scene <- function() {
  list(
    camera     = list(eye = list(x=1.6, y=1.6, z=1.1)),
    xaxis      = list(title="X", showgrid=TRUE, zeroline=TRUE),
    yaxis      = list(title="Y", showgrid=TRUE, zeroline=TRUE),
    zaxis      = list(title="Z", showgrid=TRUE, zeroline=TRUE),
    aspectmode = "data"
  )
}

panel_layout <- function(title_text) {
  list(
    title         = list(text=title_text, font=list(size=13), x=0.5, xanchor="center"),
    margin        = list(t=50, b=5, l=5, r=5),
    paper_bgcolor = "#f5f7fa",
    scene         = make_scene()
  )
}

# === Step 3: Create interactive graphics ===
build_panel <- function(demo_id) {
  df <- raw_data %>% filter(ID == demo_id)
  
  s0 <- as.matrix(df[, c("Start_X","Start_Y","Start_Z")])
  e0 <- as.matrix(df[, c("End_X",  "End_Y",  "End_Z"  )])
  center_raw <- as.numeric(df[1, c("Pos_X","Pos_Y","Pos_Z")])
  
  longest_idx <- which.max(df$Length)
  arr_scale   <- max(dist(s0)) * 0.25
  half_sz     <- max(dist(s0)) * 0.55
  
  # --- SVD normal vector ---
  dx  <- e0[,1] - s0[,1]
  dy  <- e0[,2] - s0[,2]
  dz  <- e0[,3] - s0[,3]
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  
  if (sum(valid) >= 3) {
    U          <- cbind(dx[valid]/len[valid],
                        dy[valid]/len[valid],
                        dz[valid]/len[valid])
    normal_svd <- svd(U)$v[, 3]
  } else {
    normal_svd <- as.numeric(df[1, c("Norm_X","Norm_Y","Norm_Z")])
  }
  normal_svd <- normal_svd / sqrt(sum(normal_svd^2))
  if (normal_svd[3] < 0) normal_svd <- -normal_svd
  
  # Step 1: Rotate the SVD normal vector to the Z-axis
  R1        <- get_rot_matrix(normal_svd, c(0,0,1))
  s1        <- s0 %*% t(R1)
  e1        <- e0 %*% t(R1)
  center_r1 <- as.numeric(R1 %*% center_raw)
  
  # Step 2: Moving the SVD plane center to the global origin
  s2 <- sweep(s1, 2, center_r1, "-")
  e2 <- sweep(e1, 2, center_r1, "-")
  
  # --- Panel 0: Raw data with SVD plane and normal ---
  p0 <- plot_ly() %>%
    add_scars_3d(s0[,1],s0[,2],s0[,3],
                 e0[,1],e0[,2],e0[,3], df$Length, longest_idx) %>%
    add_arrow_3d(center_raw, normal_svd, arr_scale) %>%
    add_tilted_plane_3d(center_raw, normal_svd, half_sz) %>%
    layout(panel_layout("<b>Step 0</b>: Raw data — SVD normal shown"))
  
  # --- Panel 1: Align SVD normal to Z-axis ---
  p1 <- plot_ly() %>%
    add_scars_3d(s1[,1],s1[,2],s1[,3],
                 e1[,1],e1[,2],e1[,3], df$Length, longest_idx) %>%
    add_arrow_3d(center_r1, c(0,0,1), arr_scale) %>%
    add_plane_3d(center_r1[1], center_r1[2], center_r1[3], half_sz) %>%
    layout(panel_layout("<b>Step 1</b>: Rotate — SVD normal aligned to Z-axis"))
  
  # --- Panel 2: Moving the SVD plane center to the global origin (0,0,0) ---
  p2 <- plot_ly() %>%
    add_scars_3d(s2[,1],s2[,2],s2[,3],
                 e2[,1],e2[,2],e2[,3], df$Length, longest_idx) %>%
    add_arrow_3d(c(0,0,0), c(0,0,1), arr_scale) %>%
    add_plane_3d(0, 0, 0, half_sz) %>%
    layout(panel_layout("<b>Step 2</b>: Translate — center moved to origin"))
  
  list(p0=p0, p1=p1, p2=p2)
}

# === Step 4: Combine panels ===
get_panel_json <- function(p) {
  built       <- plotly_build(p)
  data_json   <- toJSON(built$x$data,   auto_unbox=TRUE, null="null", force=TRUE)
  layout_json <- toJSON(built$x$layout, auto_unbox=TRUE, null="null", force=TRUE)
  paste0('{"data":', data_json, ',"layout":', layout_json, '}')
}

all_ids     <- unique(raw_data$ID)
panels_list <- lapply(as.character(all_ids), build_panel)
names(panels_list) <- as.character(all_ids)

js_data_lines <- sapply(as.character(all_ids), function(id) {
  ps <- panels_list[[id]]
  paste0(
    'allPanels["', id, '"] = {',
    '"p0":', get_panel_json(ps$p0), ',',
    '"p1":', get_panel_json(ps$p1), ',',
    '"p2":', get_panel_json(ps$p2),
    '};'
  )
})
js_data_block <- paste(js_data_lines, collapse="\n")
ids_json      <- toJSON(as.character(all_ids), auto_unbox=FALSE)

# === Step 5: Export figure (HTML format) ===
grid <- browsable(
  tagList(
    tags$script(src="https://cdn.plot.ly/plotly-2.27.0.min.js"),
    
    tags$h3(
      style = "font-family:sans-serif; text-align:center; margin:16px 0 4px;",
      "Core Alignment Pipeline  (SVD normal)"
    ),
    
    tags$div(
      style = "text-align:center; margin-bottom:10px;",
      tags$label("Select specimen: ",
                 style = "font-family:sans-serif; font-size:13px;"),
      tags$select(
        id    = "specimenSelect",
        style = "font-size:13px; padding:3px 8px;",
        lapply(as.character(all_ids), function(id) tags$option(value=id, id))
      )
    ),
    
    tags$p(
      style = "font-family:sans-serif; text-align:center; color:#666;
               margin:0 0 12px; font-size:13px;",
      HTML("&#9642; <b style='color:steelblue'>Blue</b> = Flaking scars &nbsp;|&nbsp;
            <b style='color:pink'>Pink</b> = Longest scar &nbsp;|&nbsp;
            <b style='color:red'>Red arrow</b> = SVD normal &nbsp;|&nbsp;
            <b style='color:lightgray'>Gray plane</b> = SVD best-fit plane")
    ),
    
    tags$div(
      style = "display:grid; grid-template-columns:1fr 1fr 1fr;
               gap:8px; padding:0 12px 12px;",
      tags$div(id="plot0", style="height:480px;"),
      tags$div(id="plot1", style="height:480px;"),
      tags$div(id="plot2", style="height:480px;")
    ),
    
    tags$script(HTML(paste0(
      "var allPanels = {};\n",
      "var allIds = ", ids_json, ";\n",
      js_data_block, "\n",
      "
      function renderSpecimen(id) {
        var ps = allPanels[id];
        ['p0','p1','p2'].forEach(function(k, i) {
          Plotly.react('plot' + i, ps[k].data, ps[k].layout);
        });
      }
      renderSpecimen(allIds[0]);
      document.getElementById('specimenSelect').addEventListener('change', function() {
        renderSpecimen(this.value);
      });
      "
    )))
  )
)

htmltools::save_html(grid, "scar_alignment_interactive.html")
browseURL("scar_alignment_interactive.html")

# ==============================================================================
# Compute descriptive statistics of scar orientations for each core
directional_stats <- aligned_data %>%
  group_by(ID) %>%
  summarise(
    n = n(),
    
    # --- Mean direction vector ---
    mean_dx = mean(d_x),
    mean_dy = mean(d_y),
    mean_dz = mean(d_z),
    
    # --- Mean resultant length R̄ ---
    R_bar = sqrt(mean_dx^2 + mean_dy^2 + mean_dz^2),
    
    # --- Mean azimuth ---
    mean_azimuth = (atan2(mean_dy, mean_dx) * 180/pi + 360) %% 360,
    
    # --- Mean inclination ---
    mean_plunge = asin(pmax(-1, pmin(1,
                                     mean_dz / sqrt(mean_dx^2 + mean_dy^2 + mean_dz^2)
    ))) * 180/pi,
    
    # --- Spherical variance (1 − R̄) ---
    spherical_var = 1 - R_bar,
    
    # --- Angular dispersion ---
    angular_sd_deg = sqrt(-2 * log(pmax(R_bar, 1e-10))) * 180/pi,
  
    azimuth_list = list((atan2(d_y, d_x) * 180/pi + 360) %% 360),
    plunge_list  = list(asin(pmax(-1, pmin(1, d_z))) * 180/pi),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(c(R_bar, mean_azimuth, mean_plunge,
             spherical_var, angular_sd_deg), ~ round(.x, 3))
  )

cat("===== Spherical directional descriptive statistics =====\n\n")
directional_stats %>%
  select(ID, n, mean_azimuth, mean_plunge, R_bar,
         spherical_var, angular_sd_deg) %>%
  print(n = Inf)

# Scar orientation visualization
# Lambert equal-area projection
lambert_project <- function(dx, dy, dz) {
  flip <- dz > 0
  dx[flip] <- -dx[flip]
  dy[flip] <- -dy[flip]
  dz[flip] <- -dz[flip]
  
  r  <- sqrt(2 / pmax(1e-10, 1 - dz))
  px <- dx * r / 2
  py <- dy * r / 2
  
  data.frame(px = px, py = py, flipped = flip)
}

# Prepare projection data
proj_data <- aligned_data %>%
  rowwise() %>%
  mutate(
    proj = list(lambert_project(d_x, d_y, d_z)),
    px   = lambert_project(d_x, d_y, d_z)$px,
    py   = lambert_project(d_x, d_y, d_z)$py
  ) %>%
  ungroup()

# Mean direction projection
mean_proj <- directional_stats %>%
  rowwise() %>%
  mutate(
    mpx = lambert_project(
      cos(mean_azimuth * pi/180) * cos(mean_plunge * pi/180),
      sin(mean_azimuth * pi/180) * cos(mean_plunge * pi/180),
      sin(mean_plunge  * pi/180))$px,
    mpy = lambert_project(
      cos(mean_azimuth * pi/180) * cos(mean_plunge * pi/180),
      sin(mean_azimuth * pi/180) * cos(mean_plunge * pi/180),
      sin(mean_plunge  * pi/180))$py
  ) %>%
  ungroup()

# Reference circle
circle_df <- data.frame(
  x = cos(seq(0, 2*pi, length.out = 200)),
  y = sin(seq(0, 2*pi, length.out = 200))
)

# Compute convex hull
hull_data <- proj_data %>%
  group_by(ID) %>%
  slice(chull(px, py)) %>% 
  ungroup()

# 2D Schmidt net plot
ggplot() +
  geom_path(data = circle_df, aes(x = x, y = y),
            color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_polygon(data = hull_data,
               aes(x = px, y = py, group = factor(ID), fill = factor(ID)),
               alpha = 0.12, color = NA) +
  geom_polygon(data = hull_data,
               aes(x = px, y = py, group = factor(ID), color = factor(ID)),
               fill = NA, linewidth = 0.4, linetype = "dashed") +
  geom_point(data = proj_data,
             aes(x = px, y = py, color = factor(ID)),
             size = 2.5, alpha = 0.8) +
  geom_point(data = mean_proj,
             aes(x = mpx, y = mpy, fill = factor(ID)),
             shape = 23, size = 4, color = "white", stroke = 1) +
  annotate("text", x = 0,     y =  1.08, label = "0°",   size = 3, color = "grey40") +
  annotate("text", x =  1.12, y =  0,    label = "90°",  size = 3, color = "grey40") +
  annotate("text", x =  0,    y = -1.08, label = "180°", size = 3, color = "grey40") +
  annotate("text", x = -1.12, y =  0,    label = "270°", size = 3, color = "grey40") +
  scale_color_viridis_d(option = "viridis", name = "Core ID") +
  scale_fill_viridis_d(option = "viridis", name = "Core ID") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text     = element_blank(),
    panel.grid    = element_blank()
  )

# ============================================================
# === Spherical KDE-PCoA analysis ===
stopifnot("Typology" %in% colnames(raw_data))
cat("Data successfully loaded,", nrow(raw_data), "rows in total,",
    length(unique(raw_data$ID)), "cores,",
    length(unique(raw_data$Typology)), "types\n")
cat("Type distribution:\n")
print(table(raw_data$Typology))

# --- Step 1: SVD Alignment ---
get_rot_matrix <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  cos_theta <- sum(a * b)
  if (cos_theta < -1 + 1e-10) {
    perp <- if (abs(a[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
    v <- perp - sum(perp * a) * a
    v <- v / sqrt(sum(v^2))
    return(2 * outer(v, v) - diag(3))
  }
  if (cos_theta > 1 - 1e-10) return(diag(3))
  v <- c(a[2]*b[3] - a[3]*b[2],
         a[3]*b[1] - a[1]*b[3],
         a[1]*b[2] - a[2]*b[1])
  v_skew <- matrix(c( 0,    -v[3],  v[2],
                      v[3],  0,    -v[1],
                      -v[2],  v[1],  0   ), 3, 3, byrow = TRUE)
  diag(3) + v_skew + v_skew %*% v_skew * ((1 - cos_theta) / sum(v^2))
}

align_svd <- function(df_group) {
  
  # --- SVD Rotation---
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  
  if (sum(valid) >= 3) {
    U      <- cbind(dx[valid]/len[valid],
                    dy[valid]/len[valid],
                    dz[valid]/len[valid])
    normal <- svd(U)$v[, 3]
  } else
  
  normal <- normal / sqrt(sum(normal^2))
  if (normal[3] < 0) normal <- -normal
  
  R1 <- get_rot_matrix(normal, c(0, 0, 1))
  
  S <- as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  
  df_group$s_x <- S[,1]; df_group$s_y <- S[,2]; df_group$s_z <- S[,3]
  df_group$e_x <- E[,1]; df_group$e_y <- E[,2]; df_group$e_z <- E[,3]
  df_group$d_x <- D[,1]; df_group$d_y <- D[,2]; df_group$d_z <- D[,3]

  center         <- as.numeric(df_group[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  center_rotated <- as.numeric(R1 %*% center)
  
  df_group <- df_group %>%
    mutate(
      s_x = s_x - center_rotated[1],
      s_y = s_y - center_rotated[2],
      s_z = s_z - center_rotated[3],
      e_x = e_x - center_rotated[1],
      e_y = e_y - center_rotated[2],
      e_z = e_z - center_rotated[3]
    )
  
  return(df_group)
}

aligned_data <- raw_data %>%
  group_by(ID) %>%
  group_modify(~ align_svd(.x)) %>%
  ungroup()

# --- Step 2: Calculate unit direction vectors ---
directions <- aligned_data %>%
  mutate(
    dx  = e_x - s_x, dy = e_y - s_y, dz = e_z - s_z,
    len = sqrt(dx^2 + dy^2 + dz^2),
    ux  = dx/len, uy = dy/len, uz = dz/len
  ) %>%
  filter(len > 1e-10)

# --- Step 3: Spherical KDE ---
make_sphere_grid <- function(n_bearing = 72, n_plunge = 36) {
  bearing_seq <- seq(0, 2*pi, length.out = n_bearing + 1)[-(n_bearing+1)]
  plunge_seq  <- seq(-pi/2 * 0.95, pi/2 * 0.95, length.out = n_plunge)
  grid <- expand.grid(bearing = bearing_seq, plunge = plunge_seq)
  grid$x <- cos(grid$plunge) * cos(grid$bearing)
  grid$y <- cos(grid$plunge) * sin(grid$bearing)
  grid$z <- sin(grid$plunge)
  grid
}

sphere_grid <- make_sphere_grid(72, 36)
n_grid <- nrow(sphere_grid)
G <- as.matrix(sphere_grid[, c("x","y","z")])

fit_vmf_kde <- function(ux, uy, uz, G, bandwidth = 0.35) {
  X       <- cbind(ux, uy, uz)
  kappa   <- 1 / bandwidth^2
  dot_mat <- G %*% t(X)
  density <- rowMeans(exp(kappa * dot_mat))
  density / sum(density)
}

all_ids      <- unique(directions$ID)
n_cores      <- length(all_ids)
kde_matrix   <- matrix(NA, nrow = n_cores, ncol = n_grid,
                       dimnames = list(as.character(all_ids), NULL))
typology_vec <- character(n_cores)
n_scars_vec  <- integer(n_cores)

cat("\nFitting spherical KDE for", n_cores, "cores...\n")
for (i in seq_along(all_ids)) {
  id_i   <- all_ids[i]
  df_i   <- directions %>% filter(ID == id_i)
  kde_matrix[i, ] <- fit_vmf_kde(df_i$ux, df_i$uy, df_i$uz, G, bandwidth = 0.35)
  typology_vec[i] <- df_i$Typology[1]
  n_scars_vec[i]  <- nrow(df_i)
  cat(sprintf("  [%d/%d] ID=%-6s  Type=%-20s  Number of scars=%d\n",
              i, n_cores, id_i, typology_vec[i], nrow(df_i)))
}

# --- Step 4: Hellinger distance matrix ---
cat("\nCalculate Hellinger distance matrix...\n")

sqrt_kde <- sqrt(kde_matrix)

hellinger_dist <- matrix(0, nrow = n_cores, ncol = n_cores,
                         dimnames = list(as.character(all_ids),
                                         as.character(all_ids)))
for (i in 1:(n_cores - 1)) {
  for (j in (i+1):n_cores) {
    d <- sqrt(sum((sqrt_kde[i,] - sqrt_kde[j,])^2))
    hellinger_dist[i, j] <- d
    hellinger_dist[j, i] <- d
  }
}

cat("Hellinger distance matrix (range)：",
    round(min(hellinger_dist[hellinger_dist > 0]), 4), "~",
    round(max(hellinger_dist), 4), "\n")

# --- Step 5: PCoA ---
pcoa_result <- cmdscale(hellinger_dist, k = min(n_cores - 1, 10), eig = TRUE)

eigenvalues <- pcoa_result$eig
pos_eig     <- pmax(eigenvalues, 0)
var_exp     <- pos_eig / sum(pos_eig)
cumvar      <- cumsum(var_exp)

cat("\n===== PCoA Variance Explained =====\n")
n_show <- min(8, length(var_exp))
for (k in 1:n_show) {
  cat(sprintf("  PCo%d: %5.1f%%  (Cumulative %5.1f%%)\n",
              k, var_exp[k]*100, cumvar[k]*100))
}

neg_eig_ratio <- sum(pmin(eigenvalues, 0)) / sum(abs(eigenvalues))
if (abs(neg_eig_ratio) > 0.05) {
  cat(sprintf("\nNote: Negative eigenvalue ratio = %.1f%%, the distance matrix shows some non-Euclidean properties\n",
              abs(neg_eig_ratio) * 100))
} else {
  cat(sprintf("\nNegative eigenvalue ratio = %.1f%%, the distance matrix is approximately Euclidean, and the PCoA results are reliable\n",
              abs(neg_eig_ratio) * 100))
}

pcoa_coords <- pcoa_result$points
colnames(pcoa_coords) <- paste0("PCo", 1:ncol(pcoa_coords))

# --- Step 6: PERMANOVA ---
permanova_test <- function(dist_mat, groups, n_perm = 999) {
  n <- nrow(dist_mat)
  calc_F <- function(d, g) {
    ug <- unique(g); k <- length(ug)
    ss_w <- sum(sapply(ug, function(gi) {
      idx <- which(g == gi)
      if (length(idx) < 2) return(0)
      sum(d[idx, idx]^2) / (2 * length(idx))
    }))
    ss_t <- sum(d^2) / (2 * n)
    ss_b <- ss_t - ss_w
    (ss_b / (k - 1)) / (ss_w / (n - k))
  }
  F_obs  <- calc_F(dist_mat, groups)
  F_perm <- replicate(n_perm, calc_F(dist_mat, sample(groups)))
  p_val  <- (sum(F_perm >= F_obs) + 1) / (n_perm + 1)
  list(F = F_obs, p = p_val, n_perm = n_perm, F_dist = F_perm)
}

cat("\nRunning PERMANOVA (999 permutations)...\n")
perm_result <- permanova_test(hellinger_dist, typology_vec, n_perm = 999)

cat("\n===== PERMANOVA Results =====\n")
cat(sprintf("F statistic : %.4f\n", perm_result$F))
cat(sprintf("p-value     : %.4f  (%d permutations)\n", perm_result$p, perm_result$n_perm))
cat(sprintf("Conclusion  : %s\n",
            ifelse(perm_result$p < 0.05,
                   "Significant differences in scar distribution patterns among core types (p < 0.05)",
                   "No significant differences among types detected (p >= 0.05)")))

# --- Step 7: Visualization ---
theme_arch <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

all_types   <- unique(typology_vec)
n_types     <- length(all_types)
type_colors <- c("#2196F3","#E53935","#43A047","#FB8C00","#8E24AA","#00ACC1")
names(type_colors) <- all_types

# --- Panel A：scree plot ---
scree_df <- data.frame(
  Axis = factor(paste0("PCo", 1:n_show), levels = paste0("PCo", 1:n_show)),
  Var  = var_exp[1:n_show] * 100,
  Cum  = cumvar[1:n_show] * 100
)

scale_factor <- 100 / max(scree_df$Cum)

p_scree <- ggplot(scree_df, aes(x = Axis)) +
  geom_col(aes(y = Var), fill = "steelblue", alpha = 0.8, width = 0.6) +
  geom_line(aes(y = Cum * scale_factor, group = 1),
            color = "#E53935", linewidth = 1) +
  geom_point(aes(y = Cum * scale_factor),
             color = "#E53935", size = 2.5) +
  geom_hline(yintercept = 80 * scale_factor,
             linetype = "dotted", color = "grey50") +
  scale_y_continuous(
    name     = "Variance Explained (%)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cumulative (%)")
  ) +
  labs(title = "PCoA Scree Plot", x = "Principal Coordinate") +
  theme_arch()

# --- Panel B：PCoA satter plot ---
pcoa_df <- data.frame(
  ID       = as.character(all_ids),
  PCo1     = pcoa_coords[, 1],
  PCo2     = pcoa_coords[, 2],
  Typology = typology_vec,
  N_scars  = n_scars_vec
)

p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Typology)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(size = N_scars), alpha = 0.85) +
  geom_text(aes(label = ID), size = 2.8, vjust = -1, show.legend = FALSE) +
  scale_color_manual(values = type_colors) +
  scale_size_continuous(name = "N scars", range = c(3, 8)) +
  labs(
    title    = "PCoA",
    x        = paste0("PCo1 (", round(var_exp[1]*100, 1), "%)"),
    y        = paste0("PCo2 (", round(var_exp[2]*100, 1), "%)"),
    color    = "Typology"
  ) +
  theme_arch()

# --- Panel C：heat plot ---
order_idx    <- order(typology_vec)
ordered_ids  <- as.character(all_ids)[order_idx]
ordered_type <- typology_vec[order_idx]

heat_df <- expand.grid(Row = ordered_ids, Col = ordered_ids,
                       stringsAsFactors = FALSE)
heat_df$dist <- as.vector(hellinger_dist[order_idx, order_idx])
heat_df$Row  <- factor(heat_df$Row, levels = ordered_ids)
heat_df$Col  <- factor(heat_df$Col, levels = ordered_ids)

p_heatmap <- ggplot(heat_df, aes(x = Col, y = Row, fill = dist)) +
  geom_tile() +
  scale_fill_gradient2(
    low  = "white", mid = "#90CAF9", high = "#0D47A1",
    midpoint = max(heat_df$dist) / 2,
    name = "Hellinger\ndistance"
  ) +
  labs(
    title    = "Hellinger Distance Matrix",
    subtitle = "sorted by typology",
    x = NULL, y = NULL
  ) +
  theme_arch() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 7),
    legend.position = "right"
  )

# --- Panel D：PERM distribution ---
perm_df <- data.frame(F_perm = perm_result$F_dist)

p_perm <- ggplot(perm_df, aes(x = F_perm)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.75, color = "white") +
  geom_vline(xintercept = perm_result$F,
             color = "#E53935", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
           x     = perm_result$F,
           y     = Inf,
           label = sprintf("Observed F = %.3f", perm_result$F),
           hjust = -0.08, vjust = 1.5,
           color = "#E53935", size = 3.5) +
  labs(
    title    = "PERMANOVA Permutation Distribution",
    subtitle = sprintf("F = %.3f  |  p = %.3f  |  %d permutations",
                       perm_result$F, perm_result$p, perm_result$n_perm),
    x = "F statistic", y = "Count"
  ) +
  theme_arch()

# --- Step 8: patchwork ---
fig <- (p_scree | p_pcoa) / (p_heatmap | p_perm) +
  plot_annotation(
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

print(fig)






output_dir <- "H:/SDG_Lithic_Analysis/analysis/data/drived_data"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(sphere_grid, 
          file.path(output_dir, "sphere_grid.csv"), 
          row.names = FALSE)

kde_df <- as.data.frame(kde_matrix)
kde_df$ID <- rownames(kde_matrix)
kde_df$Typology <- typology_vec
write.csv(kde_df, 
          file.path(output_dir, "kde_matrix.csv"), 
          row.names = FALSE)

