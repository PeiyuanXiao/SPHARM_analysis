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
# ==============================================================================
# --- Load Data ---
raw_data <- read_excel("analysis/data/raw_data/Scar_orientation_data.xlsx")

# --- Alignment Function (based on Lin et al. 2024) ---
align_lin2024 <- function(df_group) {
  
  # === Step 1: Rotation — Align the normal to the Z-axis ===
  normal <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal <- normal / sqrt(sum(normal^2))
  
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
    return(diag(3) + v_skew + v_skew %*% v_skew * ((1 - cos_theta) / sum(v^2)))
  }
  
  R1 <- get_rot_matrix(normal, c(0, 0, 1))
  
  # Rotate all points and direction vectors
  S <- as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  
  df_group$s_x <- S[,1]; df_group$s_y <- S[,2]; df_group$s_z <- S[,3]
  df_group$e_x <- E[,1]; df_group$e_y <- E[,2]; df_group$e_z <- E[,3]
  df_group$d_x <- D[,1]; df_group$d_y <- D[,2]; df_group$d_z <- D[,3]
  
  # === Step 2: Translation — Move longest scar start point to (0, 0) ===
  # Find the longest scar and use the rotated starting point XY coordinates as the translation offset
  longest_idx   <- which.max(df_group$Length)
  shift_x       <- df_group$s_x[longest_idx]
  shift_y       <- df_group$s_y[longest_idx]
  
  df_group <- df_group %>%
    mutate(
      s_x = s_x - shift_x,
      s_y = s_y - shift_y,
      e_x = e_x - shift_x,
      e_y = e_y - shift_y
    )
  
  return(df_group)
}

# --- Batch Execution ---
aligned_data_lin2024 <- raw_data %>%
  group_by(ID) %>%
  group_modify(~ align_lin2024(.x)) %>%
  ungroup()

# --- Verification ---
# --- check whether the starting point of the longest scar is at (0, 0) ---
aligned_data_lin2024 %>%
  group_by(ID) %>%
  slice(which.max(Length)) %>%
  select(ID, s_x, s_y, s_z) %>%
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
build_panels <- function(demo_id) {
  df <- raw_data %>% filter(ID == demo_id)
  
  s0 <- as.matrix(df[, c("Start_X","Start_Y","Start_Z")])
  e0 <- as.matrix(df[, c("End_X",  "End_Y",  "End_Z"  )])
  normal_raw <- as.numeric(df[1, c("Norm_X","Norm_Y","Norm_Z")])
  normal_raw <- normal_raw / sqrt(sum(normal_raw^2))
  center_raw <- as.numeric(df[1, c("Pos_X","Pos_Y","Pos_Z")])
  
  longest_idx <- which.max(df$Length)
  arr_scale   <- max(dist(s0)) * 0.25
  half_sz     <- max(dist(s0)) * 0.55
  
  # Run alignment pipeline
  R1        <- get_rot_matrix(normal_raw, c(0,0,1))
  s1        <- s0 %*% t(R1); e1 <- e0 %*% t(R1)
  center_r1 <- as.numeric(R1 %*% center_raw)
  shift_x <- s1[longest_idx, 1]
  shift_y <- s1[longest_idx, 2]
  s2 <- s1; e2 <- e1
  s2[,1] <- s1[,1] - shift_x; s2[,2] <- s1[,2] - shift_y
  e2[,1] <- e1[,1] - shift_x; e2[,2] <- e1[,2] - shift_y
  z_longest <- s2[longest_idx, 3]
  
  # --- Panel 0: Raw data ---
  p0 <- plot_ly() %>%
    add_scars_3d(s0[,1],s0[,2],s0[,3], e0[,1],e0[,2],e0[,3], df$Length, longest_idx) %>%
    add_arrow_3d(center_raw, normal_raw, arr_scale) %>%
    add_tilted_plane_3d(center_raw, normal_raw, half_sz) %>%
    layout(panel_layout("<b>Step 0</b>: Raw data — arbitrary orientation"))
  
  # --- Panel 1: Rotate the normal vector to the Z-axis ---
  p1 <- plot_ly() %>%
    add_scars_3d(s1[,1],s1[,2],s1[,3], e1[,1],e1[,2],e1[,3], df$Length, longest_idx) %>%
    add_arrow_3d(center_r1, c(0,0,1), arr_scale) %>%
    add_plane_3d(center_r1[1], center_r1[2], center_r1[3], half_sz) %>%
    layout(panel_layout("<b>Step 1</b>: Rotate — normal aligned to Z-axis"))
  
  # --- Panel 2: Translate the starting point of the longest scar to (0,0) ---
  p2 <- plot_ly() %>%
    add_scars_3d(s2[,1],s2[,2],s2[,3], e2[,1],e2[,2],e2[,3], df$Length, longest_idx) %>%
    add_arrow_3d(c(0, 0, z_longest), c(0,0,1), arr_scale) %>%
    add_plane_3d(0, 0, z_longest, half_sz) %>%
    layout(panel_layout("<b>Step 2 (Lin 2024)</b>: Translate — longest scar start to (0,0)"))
  
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
panels_list <- lapply(as.character(all_ids), build_panels)
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
      "Core Alignment Pipeline"
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
      style = "font-family:sans-serif; text-align:center; color:#666; margin:0 0 12px; font-size:13px;",
      HTML("&#9642; <b style='color:steelblue'>Blue</b> = Flaking scars &nbsp;|&nbsp;
      <b style='color:pink'>Pink</b> = Longest scar &nbsp;|&nbsp;
      <b style='color:red'>Red arrow</b> = Plane normal &nbsp;|&nbsp;
      <b style='color:lightgray'>Gray plane</b> = Best fitting plane &nbsp;|&nbsp;
      Panel 4 = Top-down view of final alignment")
    ),
    
    tags$div(
      style = "display:grid; grid-template-columns:1fr 1fr 1fr; gap:8px; padding:0 12px 12px;",
      tags$div(id="plot0", style="height:500px;"),
      tags$div(id="plot1", style="height:500px;"),
      tags$div(id="plot2", style="height:500px;")
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

htmltools::save_html(grid, "scar_alignment_interactive_Lin_method.html")
browseURL("scar_alignment_interactive_Lin_method.html")

# ============================================================
# === Spherical KDE-PCoA analysis ===
raw_data <- read_excel("analysis/data/Scar_orientation_data.xlsx")

stopifnot("Typology" %in% colnames(raw_data))
cat("Data successfully loaded,", nrow(raw_data), "rows in total,",
    length(unique(raw_data$ID)), "cores,",
    length(unique(raw_data$Typology)), "types\n")
cat("Type distribution:\n")
print(table(raw_data$Typology))

# --- Step 1: Geomagic Alignment ---
get_rot_matrix <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  cos_theta <- sum(a * b)
  if (cos_theta < -1 + 1e-10) {
    perp <- if (abs(a[1]) < 0.9) c(1,0,0) else c(0,1,0)
    v <- perp - sum(perp*a)*a; v <- v / sqrt(sum(v^2))
    return(2*outer(v,v) - diag(3))
  }
  if (cos_theta > 1 - 1e-10) return(diag(3))
  v <- c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
  v_skew <- matrix(c(0,-v[3],v[2], v[3],0,-v[1], -v[2],v[1],0), 3,3,byrow=TRUE)
  diag(3) + v_skew + v_skew %*% v_skew * ((1-cos_theta)/sum(v^2))
}

align_lin2024 <- function(df_group) {
  normal <- as.numeric(df_group[1, c("Norm_X","Norm_Y","Norm_Z")])
  normal <- normal / sqrt(sum(normal^2))
  R1 <- get_rot_matrix(normal, c(0,0,1))
  S <- as.matrix(df_group[, c("Start_X","Start_Y","Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",  "End_Y",  "End_Z" )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  df_group$s_x <- S[,1]; df_group$s_y <- S[,2]; df_group$s_z <- S[,3]
  df_group$e_x <- E[,1]; df_group$e_y <- E[,2]; df_group$e_z <- E[,3]
  df_group$d_x <- D[,1]; df_group$d_y <- D[,2]; df_group$d_z <- D[,3]
  longest_idx <- which.max(df_group$Length)
  shift_x <- df_group$s_x[longest_idx]
  shift_y <- df_group$s_y[longest_idx]
  df_group %>% mutate(
    s_x = s_x - shift_x, s_y = s_y - shift_y,
    e_x = e_x - shift_x, e_y = e_y - shift_y
  )
}

aligned_data <- raw_data %>%
  group_by(ID) %>%
  group_modify(~ align_lin2024(.x)) %>%
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

heat_df$Row <- factor(heat_df$Row, levels = ordered_ids)
heat_df$Col <- factor(heat_df$Col, levels = ordered_ids)

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

# ==============================================================================
raw_data <- read_excel("analysis/data/Scar_orientation_data.xlsx")

stopifnot("Typology" %in% colnames(raw_data))
cat("Data successfully loaded,", nrow(raw_data), "rows in total,",
    length(unique(raw_data$ID)), "cores,",
    length(unique(raw_data$Typology)), "types\n")

# Step 1: Alignment
get_rot_matrix <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  cos_theta <- sum(a * b)
  if (cos_theta < -1 + 1e-10) {
    perp <- if (abs(a[1]) < 0.9) c(1,0,0) else c(0,1,0)
    v <- perp - sum(perp*a)*a; v <- v / sqrt(sum(v^2))
    return(2*outer(v,v) - diag(3))
  }
  if (cos_theta > 1 - 1e-10) return(diag(3))
  v <- c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
  v_skew <- matrix(c(0,-v[3],v[2], v[3],0,-v[1],-v[2],v[1],0), 3,3,byrow=TRUE)
  diag(3) + v_skew + v_skew %*% v_skew * ((1-cos_theta)/sum(v^2))
}

align_lin2024 <- function(df_group) {
  normal <- as.numeric(df_group[1, c("Norm_X","Norm_Y","Norm_Z")])
  normal <- normal / sqrt(sum(normal^2))
  R1 <- get_rot_matrix(normal, c(0,0,1))
  S <- as.matrix(df_group[, c("Start_X","Start_Y","Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",  "End_Y",  "End_Z" )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  df_group$s_x <- S[,1]; df_group$s_y <- S[,2]; df_group$s_z <- S[,3]
  df_group$e_x <- E[,1]; df_group$e_y <- E[,2]; df_group$e_z <- E[,3]
  df_group$d_x <- D[,1]; df_group$d_y <- D[,2]; df_group$d_z <- D[,3]
  longest_idx <- which.max(df_group$Length)
  shift_x <- df_group$s_x[longest_idx]
  shift_y <- df_group$s_y[longest_idx]
  df_group %>% mutate(
    s_x = s_x - shift_x, s_y = s_y - shift_y,
    e_x = e_x - shift_x, e_y = e_y - shift_y
  )
}

aligned_data <- raw_data %>%
  group_by(ID) %>%
  group_modify(~ align_lin2024(.x)) %>%
  ungroup()

# Step 2: Calculate unit direction vectors
directions <- aligned_data %>%
  mutate(
    dx  = e_x - s_x, dy = e_y - s_y, dz = e_z - s_z,
    len = sqrt(dx^2 + dy^2 + dz^2),
    ux  = dx/len, uy = dy/len, uz = dz/len
  ) %>%
  filter(len > 1e-10)

# Step 3: Calculate eigenvalues and fabric indices
calc_fabric_indices <- function(ux, uy, uz) {
  # Construct matrix of unit direction vectors (n×3)
  U <- cbind(ux, uy, uz)
  n <- nrow(U)
  
  # Orientation matrix
  T_mat <- t(U) %*% U / n
  
  # Eigen decomposition, eigenvalues sorted in descending order
  eig    <- eigen(T_mat, symmetric = TRUE)
  evals  <- sort(eig$values, decreasing = TRUE) 
  evals  <- pmax(evals, 0)                       
  
  # Normalize
  evals <- evals / sum(evals)
  e1 <- evals[1]; e2 <- evals[2]; e3 <- evals[3]
  
  # Shape indices (Lin et al. 2024, p.5; after Benn 1994)
  elongation <- 1 - (e2 / e1) 
  isotropy   <- e3 / e1 
  
  list(
    e1 = e1, e2 = e2, e3 = e3,
    elongation = elongation,
    isotropy   = isotropy,
    n_scars    = n
  )
}

all_ids      <- unique(directions$ID)
results_list <- vector("list", length(all_ids))

cat("\nCalculating Fabric Indices for", length(all_ids), "cores...\n")

for (i in seq_along(all_ids)) {
  id_i  <- all_ids[i]
  df_i  <- directions %>% filter(ID == id_i)
  idx   <- calc_fabric_indices(df_i$ux, df_i$uy, df_i$uz)
  
  results_list[[i]] <- data.frame(
    ID         = id_i,
    Typology   = df_i$Typology[1],
    N_scars    = idx$n_scars,
    e1         = idx$e1,
    e2         = idx$e2,
    e3         = idx$e3,
    Elongation = idx$elongation,
    Isotropy   = idx$isotropy
  )
  
  cat(sprintf("  [%d/%d] ID=%-6s  Type=%-20s  N=%2d  Elong=%.3f  Isot=%.3f\n",
              i, length(all_ids), id_i, df_i$Typology[1],
              idx$n_scars, idx$elongation, idx$isotropy))
}

fabric_df <- do.call(rbind, results_list)

cat("\n===== Fabric Indices Summary =====\n")
print(fabric_df %>%
        group_by(Typology) %>%
        summarise(
          n          = n(),
          Elong_mean = round(mean(Elongation), 3),
          Elong_sd   = round(sd(Elongation), 3),
          Isot_mean  = round(mean(Isotropy), 3),
          Isot_sd    = round(sd(Isotropy), 3),
          .groups = "drop"
        ))

# Step 5: Visualization
# Color scheme (consistent with KDE-PCoA code)
all_types   <- unique(fabric_df$Typology)
n_types     <- length(all_types)
type_colors <- c("#2196F3","#E53935","#43A047","#FB8C00","#8E24AA","#00ACC1")
names(type_colors) <- all_types

theme_arch <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

# Benn triangle plot
h <- sqrt(3) / 2

# Triangle vertices and edges
triangle_vertices <- data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, h,   0)
)

# Grid lines (draw contour every 20%)
make_grid_lines <- function(step = 0.2) {
  lines <- list()
  vals  <- seq(step, 1 - step, by = step)
  
  for (v in vals) {
    # Contour line of Elongation (E = v): I from 0 to 1-v
    i_seq <- seq(0, 1 - v, length.out = 50)
    lines[[length(lines)+1]] <- data.frame(
      x    = v * 1 + i_seq * 0.5,
      y    = v * 0 + i_seq * h,
      type = "E", val = v
    )
    # Contour line of Isotropy (I = v): E from 0 to 1-v
    e_seq <- seq(0, 1 - v, length.out = 50)
    lines[[length(lines)+1]] <- data.frame(
      x    = e_seq * 1 + v * 0.5,
      y    = e_seq * 0 + v * h,
      type = "I", val = v
    )
    # Contour line of Planar (P = v, i.e., 1-E-I = v): E from 0 to 1-v
    e_seq2 <- seq(0, 1 - v, length.out = 50)
    i_seq2 <- 1 - v - e_seq2
    lines[[length(lines)+1]] <- data.frame(
      x    = e_seq2 * 1 + i_seq2 * 0.5,
      y    = e_seq2 * 0 + i_seq2 * h,
      type = "P", val = v
    )
  }
  do.call(rbind, lines)
}

grid_lines <- make_grid_lines(0.2)

# Cartesian coordinates for cores
fabric_df <- fabric_df %>%
  mutate(
    benn_x = Elongation * 1 + Isotropy * 0.5,
    benn_y = Isotropy * h
  )

# Vertex labels
vertex_labels <- data.frame(
  x     = c(-0.05, 1.05, 0.5),
  y     = c(-0.05, -0.05, h + 0.05),
  label = c("Planar", "Linear", "Isotropic")
)

# Tick labels
tick_labels <- data.frame(
  x     = c(seq(0.2, 0.8, 0.2) * 0.5,       
            seq(0.2, 0.8, 0.2),                    
            1 - seq(0.2, 0.8, 0.2) * 0.5),         
  y     = c(seq(0.2, 0.8, 0.2) * h,
            rep(-0.06, 4),
            seq(0.2, 0.8, 0.2) * h),
  label = c(paste0(seq(20,80,20), "%"),
            paste0(seq(20,80,20), "%"),
            paste0(seq(20,80,20), "%")),
  axis  = rep(c("I","E","P"), each = 4)
)

p_benn <- ggplot() +
  geom_path(data = grid_lines,
            aes(x = x, y = y, group = interaction(type, val)),
            color = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_path(data = triangle_vertices,
            aes(x = x, y = y),
            color = "grey30", linewidth = 0.7) +
  geom_point(data = fabric_df,
             aes(x = benn_x, y = benn_y,
                 color = Typology, size = N_scars),
             alpha = 0.85) +
  geom_text(data = fabric_df,
            aes(x = benn_x, y = benn_y, label = ID,
                color = Typology),
            size = 2.5, vjust = -1, show.legend = FALSE) +
  geom_text(data = vertex_labels,
            aes(x = x, y = y, label = label),
            fontface = "bold", size = 4, color = "grey20") +
  scale_color_manual(values = type_colors) +
  scale_size_continuous(name = "N scars", range = c(3, 8)) +
  coord_fixed(ratio = 1) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
    legend.position = "right",
    plot.margin   = margin(20, 20, 20, 20)
  )

p_benn




















