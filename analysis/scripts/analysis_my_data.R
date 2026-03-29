library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)
library(htmlwidgets)
library(vegan)
library(ggrepel)
# ============================================================
# 读取数据
# ============================================================
SPHARM_direction <- 
  read_csv("analysis/data/derived_data/SPHARM_direction.csv")
SPHARM_morphology <- 
  read_csv("analysis/data/derived_data/SPHARM_morphology.csv")
variance_direction <- 
  read_csv("analysis/data/derived_data/SPHARM_direction_variance_per_degree.csv")
variance_morphology <- 
  read_csv("analysis/data/derived_data/variance_per_degree.csv")

# ==============================================================================
# 方差贡献折线图-筛选阶数
# ==============================================================================
# 绘图函数----------------------------------------------------------------------
plot_spharm_variance <- function(df_dir, df_mor, use_log_y = TRUE) {
  
  df_dir <- df_dir %>% mutate(Feature = "Direction")
  df_mor <- df_mor %>% mutate(Feature = "Morphology")
  
  df_combined <- bind_rows(df_dir, df_mor)
  
  p <- ggplot(df_combined, aes(x = degree, y = variance, 
                               color = Feature, shape = Feature)) +
    geom_line(size = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    scale_x_continuous(breaks = seq(min(df_combined$degree), 
                                    max(df_combined$degree), by = 1)) +
    scale_color_manual(values = c("Direction" = "#FFBAE0", 
                                  "Morphology" = "#A1C2E6")) +
    theme_bw() +
    labs(
      title = "SPHARM Variance per Degree",
      x = "Spherical Harmonic Degree (l)",
      y = "Variance",
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.minor.x = element_blank() 
    )
  
}

# 出图函数----------------------------------------------------------------------
variance_plot <- plot_spharm_variance(variance_direction, variance_morphology)
print(variance_plot)

ggsave("analysis/data/derived_data/Variance_Comparison.png", 
       plot = variance_plot, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 筛选1-5阶
# ==============================================================================
SPHARM_direction_filter <- SPHARM_direction %>%
  select(
    ID,          
    SHE,                 
    spectral_entropy,
    power_l1:power_l5  
  )

SPHARM_morphology_filter <- SPHARM_morphology %>%
  select(
    specimen_id,   
    SHE,                  
    spectral_entropy,
    power_degree_1:power_degree_5  
  )

print(SPHARM_direction_filter)
print(SPHARM_morphology_filter)

# ============================================================
# UMAP
# ============================================================
# 1. 准备数据：提取数值矩阵并进行简单预处理
umap_input <- df %>% 
  select(all_of(power_cols)) %>% 
  as.matrix()

# 2. 运行 UMAP (设置随机种子保证结果可重复)
set.seed(42) 
umap_config <- umap.defaults
umap_results <- umap(umap_input, config = umap_config)

# 3. 将结果合并回原始数据框
df_umap <- df %>%
  mutate(
    UMAP1 = umap_results$layout[, 1],
    UMAP2 = umap_results$layout[, 2]
  )

# 4. 绘制 UMAP 散点图
df_hull <- df_umap %>%
  group_by(Typology) %>%
  slice(chull(UMAP1, UMAP2)) %>%
  ungroup()

  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = transform(df_umap, Typology = NULL), 
             color = "grey90", size = 1) +
  geom_polygon(data = df_hull, aes(fill = Typology), alpha = 0.3) +
  geom_point(aes(color = Typology), size = 1.5) +
  facet_wrap(~Typology, ncol = 4) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.background = element_blank())




  # 【极度重要】：UMAP 是随机算法，设置随机种子保证每次跑出来的图一模一样
  set.seed(42) 
  
  # ============================================================
  # Step 1: 形态学数据 (Morphology) UMAP 降维
  # ============================================================
  # 1. 提取用于降维的纯数字特征矩阵 (只用 1-5 阶)
  morph_features <- SPHARM_morphology_filter %>%
    select(power_degree_1:power_degree_5) %>%
    as.matrix()
  
  # 2. 运行 UMAP 算法
  # 注意：因为你的理想模型样本量很小(约10个)，我们需要调小 n_neighbors 参数，否则 UMAP 会报错
  umap_morph <- umap(morph_features, n_neighbors = 3, random_state = 42)
  
  # 3. 构建绘图数据框
  df_umap_morph <- data.frame(
    ID    = SPHARM_morphology_filter$specimen_id,
    UMAP1 = umap_morph$layout[, 1],
    UMAP2 = umap_morph$layout[, 2]
  )
  
  # 4. 绘制形态学 UMAP 图
  p_morph <- ggplot(df_umap_morph, aes(x = UMAP1, y = UMAP2, color = ID)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE, max.overlaps = 20) +
    scale_color_viridis_d(option = "turbo") + # 使用高对比度离散配色
    theme_minimal() +
    labs(
      title = "UMAP: Overall Morphology",
      subtitle = "Based on SPHARM Degrees 1-5",
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "none", # 名字已经在图上了，关掉图例让图面更干净
      panel.border  = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # ============================================================
  # Step 2: 片疤方向数据 (Scar Direction) UMAP 降维
  # ============================================================
  # 1. 提取特征
  scar_features <- SPHARM_direction_filter %>%
    select(power_l1:power_l5) %>%
    as.matrix()
  
  # 2. 运行 UMAP (同样调整 n_neighbors)
  umap_scar <- umap(scar_features, n_neighbors = 3, random_state = 42)
  
  # 3. 构建绘图数据框
  df_umap_scar <- data.frame(
    ID    = SPHARM_direction_filter$ID,
    UMAP1 = umap_scar$layout[, 1],
    UMAP2 = umap_scar$layout[, 2]
  )
  
  # 4. 绘制片疤方向 UMAP 图
  p_scar <- ggplot(df_umap_scar, aes(x = UMAP1, y = UMAP2, color = ID)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE, max.overlaps = 20) +
    scale_color_viridis_d(option = "turbo") +
    theme_minimal() +
    labs(
      title = "UMAP: Scar Direction Pattern",
      subtitle = "Based on SPHARM Degrees 1-5",
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "none",
      panel.border  = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # ============================================================
  # Step 3: 合并图像并输出
  # ============================================================
  # 使用 patchwork 包极其简单的加号语法拼接图像
  p_combined <- p_morph + p_scar
  
  # 打印显示
  print(p_combined)
  