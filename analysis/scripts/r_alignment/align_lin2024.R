# align_lin2024.R
# Scar-direction alignment (Lin et al. 2024 method). Sourced by the `align_lin2024_csv` target.
# Used only by the rotational-invariance validation pipeline.
#
# Algorithm (two steps):
#   Step 1  Rotate    — align the measured morphological normal (Norm_X/Y/Z) to Z.
#   Step 2  Translate — move the longest scar's start point to (0, 0).
#
# Differs from align_svd.R: rotation uses the measured normal (not an SVD-fitted one),
# translation anchors the longest scar (not the centroid), and there is no XY rotation step.
#
# Input : analysis/data/raw_data/Scar_orientation_data.xlsx
# Output: analysis/data/derived_data/directions_aligned_lin2024.csv
#         analysis/output/html/scar_alignment_lin2024.html (interactive 3-panel diagnostic)

library(here)
library(readxl)
library(dplyr)
library(plotly)
library(htmltools)
library(jsonlite)
library(readr)

conflicted::conflicts_prefer(plotly::layout)

source(here("analysis/scripts/r_utils/geometry_utils.R"))
source(here("analysis/scripts/r_utils/viz3d_utils.R"))


# Scar length: use the Length column if present, else compute from endpoints.
get_scar_length <- function(df) {
  if ("Length" %in% names(df)) return(df$Length)
  dx <- df$End_X - df$Start_X
  dy <- df$End_Y - df$Start_Y
  dz <- df$End_Z - df$Start_Z
  sqrt(dx^2 + dy^2 + dz^2)
}


# 1. Load data -----------------------------------------------------------------

xlsx_path <- here("analysis/data/raw_data/Scar_orientation_data.xlsx")
data_1   <- read_excel(xlsx_path, sheet = 1)
data_2   <- read_excel(xlsx_path, sheet = 2)
raw_data <- bind_rows(data_1, data_2)


# 2. Alignment (Lin 2024) ------------------------------------------------------

align_lin2024 <- function(df_group) {
  
  # Unit direction vectors (same as align_svd).
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  
  df_group$Direct_X <- ifelse(valid, dx / len, 0)
  df_group$Direct_Y <- ifelse(valid, dy / len, 0)
  df_group$Direct_Z <- ifelse(valid, dz / len, 0)
  
  # Step 1: rotate morphological normal onto Z.
  normal <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal <- normal / sqrt(sum(normal^2))
  
  R1 <- get_rot_matrix(normal, c(0, 0, 1))
  
  S <- as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  
  df_group$s_x <- S[, 1]; df_group$s_y <- S[, 2]; df_group$s_z <- S[, 3]
  df_group$e_x <- E[, 1]; df_group$e_y <- E[, 2]; df_group$e_z <- E[, 3]
  df_group$d_x <- D[, 1]; df_group$d_y <- D[, 2]; df_group$d_z <- D[, 3]
  
  # Step 2: translate longest scar's start to (0, 0).
  longest_idx <- which.max(get_scar_length(df_group))
  shift_x     <- df_group$s_x[longest_idx]
  shift_y     <- df_group$s_y[longest_idx]
  
  df_group <- df_group %>%
    mutate(
      s_x = s_x - shift_x,
      s_y = s_y - shift_y,
      e_x = e_x - shift_x,
      e_y = e_y - shift_y
    )
  
  df_group
}


# 3. Run alignment -------------------------------------------------------------

aligned_data_lin2024 <- raw_data %>%
  group_by(ID) %>%
  group_modify(~ align_lin2024(.x)) %>%
  ungroup()


# 4. Check: longest scar's start should sit at (0, 0) --------------------------

aligned_data_lin2024 %>%
  group_by(ID) %>%
  slice(which.max(get_scar_length(cur_data()))) %>%
  select(ID, s_x, s_y, s_z) %>%
  print()


# 5. Interactive 3-panel diagnostic --------------------------------------------

build_panels <- function(demo_id) {
  df <- raw_data %>% filter(ID == demo_id)
  
  s0         <- as.matrix(df[, c("Start_X", "Start_Y", "Start_Z")])
  e0         <- as.matrix(df[, c("End_X",   "End_Y",   "End_Z"  )])
  normal_raw <- as.numeric(df[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_raw <- normal_raw / sqrt(sum(normal_raw^2))
  center_raw <- as.numeric(df[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  
  longest_idx <- which.max(get_scar_length(df))
  arr_scale   <- max(dist(s0)) * 0.25
  half_sz     <- max(dist(s0)) * 0.55
  
  # Replay the alignment to visualise each intermediate state.
  R1        <- get_rot_matrix(normal_raw, c(0, 0, 1))
  s1        <- s0 %*% t(R1)
  e1        <- e0 %*% t(R1)
  center_r1 <- as.numeric(R1 %*% center_raw)
  
  shift_x   <- s1[longest_idx, 1]
  shift_y   <- s1[longest_idx, 2]
  s2        <- s1; e2 <- e1
  s2[, 1]   <- s1[, 1] - shift_x; s2[, 2] <- s1[, 2] - shift_y
  e2[, 1]   <- e1[, 1] - shift_x; e2[, 2] <- e1[, 2] - shift_y
  z_longest <- s2[longest_idx, 3]
  
  # highlight_idx marks the longest scar in pink.
  p0 <- plot_ly() %>%
    add_scars_3d(s0[,1], s0[,2], s0[,3],
                 e0[,1], e0[,2], e0[,3],
                 highlight_idx = longest_idx) %>%
    add_arrow_3d(center_raw, normal_raw, arr_scale) %>%
    add_tilted_plane_3d(center_raw, normal_raw, half_sz) %>%
    layout(panel_layout("<b>Step 0</b>: Raw data — arbitrary orientation"))
  
  p1 <- plot_ly() %>%
    add_scars_3d(s1[,1], s1[,2], s1[,3],
                 e1[,1], e1[,2], e1[,3],
                 highlight_idx = longest_idx) %>%
    add_arrow_3d(center_r1, c(0, 0, 1), arr_scale) %>%
    add_plane_3d(center_r1[1], center_r1[2], center_r1[3], half_sz) %>%
    layout(panel_layout("<b>Step 1</b>: Rotate — normal aligned to Z-axis"))
  
  p2 <- plot_ly() %>%
    add_scars_3d(s2[,1], s2[,2], s2[,3],
                 e2[,1], e2[,2], e2[,3],
                 highlight_idx = longest_idx) %>%
    add_arrow_3d(c(0, 0, z_longest), c(0, 0, 1), arr_scale) %>%
    add_plane_3d(0, 0, z_longest, half_sz) %>%
    layout(panel_layout("<b>Step 2 (Lin 2024)</b>: Translate — longest scar start to (0,0)"))
  
  list(p0 = p0, p1 = p1, p2 = p2)
}


# 6. Build panels and export HTML ----------------------------------------------

all_ids            <- unique(raw_data$ID)
panels_list        <- lapply(as.character(all_ids), build_panels)
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
js_data_block <- paste(js_data_lines, collapse = "\n")
ids_json      <- toJSON(as.character(all_ids), auto_unbox = FALSE)

grid <- browsable(
  tagList(
    tags$script(src = "https://cdn.plot.ly/plotly-2.27.0.min.js"),
    
    tags$h3(
      style = "font-family:sans-serif; text-align:center; margin:16px 0 4px;",
      "Core Alignment Pipeline (Lin 2024)"
    ),
    
    tags$div(
      style = "text-align:center; margin-bottom:10px;",
      tags$label("Select specimen: ",
                 style = "font-family:sans-serif; font-size:13px;"),
      tags$select(
        id    = "specimenSelect",
        style = "font-size:13px; padding:3px 8px;",
        lapply(as.character(all_ids), function(id) tags$option(value = id, id))
      )
    ),
    
    tags$p(
      style = "font-family:sans-serif; text-align:center; color:#666;
               margin:0 0 12px; font-size:13px;",
      HTML("&#9642; <b style='color:steelblue'>Blue</b> = Flaking scars &nbsp;|&nbsp;
            <b style='color:pink'>Pink</b> = Longest scar &nbsp;|&nbsp;
            <b style='color:red'>Red arrow</b> = Plane normal &nbsp;|&nbsp;
            <b style='color:lightgray'>Gray plane</b> = Best fitting plane")
    ),
    
    tags$div(
      style = "display:grid; grid-template-columns:1fr 1fr 1fr;
               gap:8px; padding:0 12px 12px;",
      tags$div(id = "plot0", style = "height:500px;"),
      tags$div(id = "plot1", style = "height:500px;"),
      tags$div(id = "plot2", style = "height:500px;")
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

out_path <- here("analysis/output/html/scar_alignment_lin2024.html")
htmltools::save_html(grid, out_path)
# browseURL(out_path)  # disabled: avoids opening a browser during tar_make()


# 7. Export aligned directions for the Python SPHARM step ----------------------

aligned_lin2024_export <- aligned_data_lin2024 %>%
  select(ID, any_of("Typology"), ux = d_x, uy = d_y, uz = d_z)

write_csv(aligned_lin2024_export,
          here("analysis/data/derived_data/directions_aligned_lin2024.csv"))
cat("Saved: directions_aligned_lin2024.csv\n")
