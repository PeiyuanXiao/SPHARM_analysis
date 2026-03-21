# === plane_deviation_analysis.R ===
library(readxl)
library(dplyr)
library(ggplot2)

raw_data <- read_excel("analysis/data/raw_data/Scar_orientation_data.xlsx")

compute_plane_deviation <- function(df_group) {
  
  # Geomagic 法线
  normal_geo <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_geo <- normal_geo / sqrt(sum(normal_geo^2))
  
  # SVD 法线
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
  if (sum(normal_geo * normal_svd) < 0) normal_svd <- -normal_svd
  
  # 计算夹角
  cos_angle <- min(1, max(-1, sum(normal_geo * normal_svd)))
  angle_deg <- acos(cos_angle) * 180 / pi
  
  data.frame(
    ID        = df_group$ID[1],
    angle_deg = round(angle_deg, 2)
  )
}

deviation_results <- raw_data %>%
  group_by(ID) %>%
  group_map(~ compute_plane_deviation(.x), .keep = TRUE) %>%
  bind_rows()

print(deviation_results)








# 添加分组标签
deviation_results <- deviation_results %>%
  mutate(
    group = if_else(startsWith(as.character(ID), "IM"), "理想模型", "真实标本"),
    ID    = factor(ID, levels = deviation_results$ID[order(deviation_results$angle_deg)])
  )

ggplot(deviation_results, aes(y = ID)) +
  
  # 哑铃横线：从 0 到角度值
  geom_segment(
    aes(x = 0, xend = angle_deg, yend = ID, color = group),
    linewidth = 0.8, alpha = 0.6
  ) +
  
  # 左端点（0°参考点）
  geom_point(
    aes(x = 0, color = group),
    size = 2.5, shape = 21, fill = "white", stroke = 1.2
  ) +
  
  # 右端点（实际角度）
  geom_point(
    aes(x = angle_deg, color = group),
    size = 3
  ) +
  
  # 角度数值标签
  geom_text(
    aes(x = angle_deg, label = paste0(angle_deg, "°")),
    hjust = -0.3, size = 3, color = "grey40"
  ) +
  
  # 参考线：90° 虚线
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  
  scale_color_manual(
    values = c("理想模型" = "#1D9E75", "真实标本" = "#D85A30"),
    name   = NULL
  ) +
  
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 30, 60, 90),
    labels = c("0°", "30°", "60°", "90°")
  ) +
  
  labs(
    x = "两平面法线夹角",
    y = NULL,
  ) +
  
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 11, face = "bold"),
    axis.text.y        = element_text(size = 10),
    plot.caption       = element_text(color = "grey50", size = 9)
  )




