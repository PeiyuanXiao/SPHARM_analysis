# ==============================================================================
# plane_angle.R
# 形态平面与技术平面的偏差分析
#
# 计算每个标本两个最优拟合平面之间的夹角和距离：
#   - 形态平面：Geomagic 软件拟合的形态法线（Norm_X/Y/Z）
#   - 技术平面：从片疤向量 SVD 分解得到的最优拟合平面
#
# 分组分析：
#   - 理想模型（IM）：按标本逐个计算
#   - 考古标本：按 Layer、Raw_mat、Layer × Raw_mat 分组计算
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 1: IM，sheet 2: 考古）
#   - analysis/data/raw_data/SDG_core_metric.xlsx
# 输出：
#   控制台打印各分组的偏差结果（angle_deg, mean_dist）
# ==============================================================================

library(here)
library(readxl)
library(dplyr)


# ==============================================================================
# 1. 读取数据
# ==============================================================================

IM_data   <- read_excel(here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 1)
arch_data <- read_excel(here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 2)
metric_data <- read_excel(here("analysis/data/raw_data/SDG_core_metric.xlsx"))

merged_arch_data <- arch_data %>%
  left_join(
    metric_data %>% select(ID, Layer, Raw_mat, Core_type_Li_merged),
    by = "ID"
  )


# ==============================================================================
# 2. 偏差计算函数
# ==============================================================================

# compute_plane_deviation()
# 计算形态平面与技术平面之间的夹角和平均垂直距离
#
# 参数：
#   df_group : 单个分组的数据框，需包含
#              Norm_X/Y/Z（形态法线）、Pos_X/Y/Z（平面参考点）
#              Start_X/Y/Z、End_X/Y/Z（刮痕端点）
#
# 返回：
#   data.frame，含：
#     angle_deg : 两平面法线夹角（度）
#     mean_dist : 刮痕端点到形态平面的平均垂直距离

compute_plane_deviation <- function(df_group) {
  
  # --- 形态最优拟合平面（Geomagic 法线）---
  normal_geo <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_geo <- normal_geo / sqrt(sum(normal_geo^2))
  p0 <- as.numeric(df_group[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  
  # --- 技术最优拟合平面（刮痕向量 SVD）---
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  if (sum(valid) < 3) return(NULL)
  
  U          <- cbind(dx[valid] / len[valid],
                      dy[valid] / len[valid],
                      dz[valid] / len[valid])
  normal_svd <- svd(U)$v[, 3]
  normal_svd <- normal_svd / sqrt(sum(normal_svd^2))
  
  # 统一法线方向（确保两法线在同侧）
  if (sum(normal_geo * normal_svd) < 0) normal_svd <- -normal_svd
  
  # --- 夹角（度）---
  cos_angle <- min(1, max(-1, sum(normal_geo * normal_svd)))
  angle_deg <- acos(cos_angle) * 180 / pi
  
  # --- 刮痕端点到形态平面的平均垂直距离 ---
  endpoints <- rbind(
    as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]),
    as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )])
  )
  dist_to_plane <- abs(
    (endpoints - matrix(p0, nrow(endpoints), 3, byrow = TRUE)) %*% normal_geo
  )
  mean_dist <- mean(dist_to_plane)
  
  data.frame(
    angle_deg = round(angle_deg, 2),
    mean_dist = round(mean_dist, 4)
  )
}


# ==============================================================================
# 3. 理想模型（IM）：按标本逐个计算
# ==============================================================================

cat("===== 理想模型：各标本平面偏差 =====\n")
deviation_IM <- IM_data %>%
  group_by(ID) %>%
  group_modify(~ compute_plane_deviation(.x))
print(deviation_IM)


# ==============================================================================
# 4. 考古标本：按 Layer 分组
# ==============================================================================

cat("\n===== 考古标本：按 Layer 分组 =====\n")
deviation_arch_layer <- merged_arch_data %>%
  group_by(Layer) %>%
  group_modify(~ compute_plane_deviation(.x))
print(deviation_arch_layer)


# ==============================================================================
# 5. 考古标本：按 Raw_mat 分组
# ==============================================================================

cat("\n===== 考古标本：按 Raw_mat 分组 =====\n")
deviation_arch_rawmat <- merged_arch_data %>%
  group_by(Raw_mat) %>%
  group_modify(~ compute_plane_deviation(.x))
print(deviation_arch_rawmat)


# ==============================================================================
# 6. 考古标本：Layer × Raw_mat 交叉分组（控制变量）
# ==============================================================================

cat("\n===== 考古标本：Layer × Raw_mat 交叉分组（Layer 3 & 4）=====\n")
cat("各组样本量：\n")
print(table(merged_arch_data$Layer, merged_arch_data$Raw_mat))

deviation_arch_controlled <- merged_arch_data %>%
  filter(Layer %in% c(3, 4)) %>%
  group_by(Layer, Raw_mat) %>%
  group_modify(~ compute_plane_deviation(.x))

print(deviation_arch_controlled)

