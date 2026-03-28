library(readxl)
library(dplyr)
library(ggplot2)

# === Input data ===
# ideal model
IM_data <- read_excel("analysis/data/raw_data/Scar_orientation_data.xlsx", sheet = 1)

# archeological model
arch_data <- read_excel("analysis/data/raw_data/Scar_orientation_data.xlsx", sheet = 2)
metric_data <- read_excel("analysis/data/raw_data/SDG_core_metric.xlsx")

merged_arch_data <- arch_data %>%
  left_join(
    metric_data %>% select(ID, Layer, Raw_mat, Core_type_Li_merged), 
    by = "ID" 
  )

head(merged_arch_data)


# === Computing plane angle and distance ===
compute_plane_deviation <- function(df_group) {
  
  # --- Morphological best-fitting plane ---
  normal_geo <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_geo <- normal_geo / sqrt(sum(normal_geo^2))
  p0 <- as.numeric(df_group[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  
  # --- Technological best-fitting plane ---
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  if (sum(valid) < 3) return(NULL)
  
  U          <- cbind(dx[valid]/len[valid],
                      dy[valid]/len[valid],
                      dz[valid]/len[valid])
  normal_svd <- svd(U)$v[, 3]
  normal_svd <- normal_svd / sqrt(sum(normal_svd^2))
  
  # 统一法线方向
  if (sum(normal_geo * normal_svd) < 0) normal_svd <- -normal_svd
  
  # --- 夹角 ---
  cos_angle <- min(1, max(-1, sum(normal_geo * normal_svd)))
  angle_deg <- acos(cos_angle) * 180 / pi
  
  # --- 片疤端点到 Geomagic 形态平面的平均垂直距离 ---      # ← 修改
  endpoints <- rbind(                                        # ← 修改
    as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]),  # ← 修改
    as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )])   # ← 修改
  )                                                          # ← 修改
  # 点到平面距离 = |(point - p0) · normal_geo|              # ← 修改
  dist_to_plane <- abs(                                      # ← 修改
    (endpoints - matrix(p0, nrow(endpoints), 3, byrow = TRUE)) %*% normal_geo  # ← 修改
  )                                                          # ← 修改
  mean_dist <- mean(dist_to_plane)                           # ← 修改
  
  data.frame(
    angle_deg = round(angle_deg, 2),
    mean_dist = round(mean_dist, 4)                          # ← 修改
  )
}

deviation_results <- IM_data %>%
  group_by(ID) %>%
  group_map(~ compute_plane_deviation(.x), .keep = TRUE) %>%
  bind_rows()

print(deviation_results)


deviation_results_arch <- merged_arch_data %>%
  group_by(Layer) %>% 
  group_modify(~ compute_plane_deviation(.x))

print(deviation_results_arch)

deviation_results <- merged_arch_data %>%
  group_by(Raw_mat) %>%
  group_modify(~ compute_plane_deviation(.x))

print(deviation_results)


table(merged_arch_data$Layer, merged_arch_data$Raw_mat)

deviation_results_controlled <- merged_arch_data %>%
  filter(Layer %in% c(3, 4)) %>%
  group_by(Layer, Raw_mat) %>%
  group_modify(~ compute_plane_deviation(.x))

print(deviation_results_controlled)
