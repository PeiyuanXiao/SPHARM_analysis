# ==============================================================================
# spharm_analysis.R
# SPHARM 特征分析：方差贡献筛选阶数 + UMAP 降维可视化
#
# 分析流程：
#   1. 读取 SPHARM 功率谱数据（方向 + 形态）及样本元数据
#   2. 绘制各阶方差贡献折线图，辅助判断保留阶数
#   3. 筛选 1–5 阶功率谱特征
#   4. UMAP 降维，按 Layer / Core type / Raw material 分别着色出图
#   5. 保存筛选后的数据供 mantel_cia.R 使用
#
# 输入：
#   - analysis/data/derived_data/SPHARM_direction.csv
#   - analysis/data/derived_data/SPHARM_morphology.csv
#   - analysis/data/derived_data/SPHARM_direction_variance_per_degree.csv
#   - analysis/data/derived_data/variance_per_degree.csv
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#
# 输出：
#   - analysis/output/figures/Variance_Comparison.png
#   - analysis/output/figures/UMAP_by_Layer.png
#   - analysis/output/figures/UMAP_by_CoreType.png
#   - analysis/output/figures/UMAP_by_RawMat.png
#   - analysis/data/derived_data/SPHARM_direction_filter.rds  ← 供 mantel_cia.R 读取
#   - analysis/data/derived_data/SPHARM_morphology_filter.rds ← 供 mantel_cia.R 读取
# ==============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)
library(vegan)
library(ggrepel)
library(readxl)

set.seed(42)


# ==============================================================================
# 1. 读取数据
# ==============================================================================

SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"))
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"))

variance_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction_variance_per_degree.csv"))
variance_morphology <- read_csv(here("analysis/data/derived_data/variance_per_degree.csv"))

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat)


# ==============================================================================
# 2. 方差贡献折线图 — 辅助判断保留阶数
# ==============================================================================

plot_spharm_variance <- function(df_dir, df_mor) {
  
  df_dir <- df_dir %>% mutate(Feature = "Direction")
  df_mor <- df_mor %>% mutate(Feature = "Morphology")
  df_combined <- bind_rows(df_dir, df_mor)
  
  ggplot(df_combined, aes(x = degree, y = variance,
                          color = Feature, shape = Feature)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    scale_x_continuous(
      breaks = seq(min(df_combined$degree), max(df_combined$degree), by = 1)
    ) +
    scale_color_manual(values = c("Direction" = "#FFBAE0",
                                  "Morphology" = "#A1C2E6")) +
    theme_bw() +
    labs(
      title = "SPHARM Variance per Degree",
      x     = "Spherical Harmonic Degree (l)",
      y     = "Variance"
    ) +
    theme(
      plot.title         = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position    = "right",
      legend.title       = element_text(face = "bold"),
      panel.grid.minor.x = element_blank()
    )
}

variance_plot <- plot_spharm_variance(variance_direction, variance_morphology)
print(variance_plot)

ggsave(
  here("analysis/output/figures/Variance_Comparison.png"),
  plot = variance_plot, width = 8, height = 6, dpi = 300
)


# ==============================================================================
# 3. 筛选 1–5 阶功率谱特征
# ==============================================================================

SPHARM_direction_filter <- SPHARM_direction %>%
  select(ID, SHE, spectral_entropy, power_l1:power_l5) %>%
  left_join(metric_data, by = "ID")

SPHARM_morphology_filter <- SPHARM_morphology %>%
  select(ID, SHE, spectral_entropy, power_l1:power_l5) %>%
  left_join(metric_data, by = "ID")

print(SPHARM_direction_filter)
print(SPHARM_morphology_filter)


# ==============================================================================
# 4. UMAP 降维可视化
# ==============================================================================

run_umap_pair <- function(morph_filter, scar_filter, color_var, color_title) {
  
  # ── Morphology ──────────────────────────────────────────────────────────────
  morph_features <- morph_filter %>% select(power_l1:power_l5) %>% as.matrix()
  umap_morph     <- umap(morph_features, n_neighbors = 5, random_state = 42)
  
  df_umap_morph <- data.frame(
    ID    = morph_filter$ID,
    color = as.factor(morph_filter[[color_var]]),
    UMAP1 = umap_morph$layout[, 1],
    UMAP2 = umap_morph$layout[, 2]
  )
  
  # ── Scar Direction ───────────────────────────────────────────────────────────
  scar_features <- scar_filter %>% select(power_l1:power_l5) %>% as.matrix()
  umap_scar     <- umap(scar_features, n_neighbors = 5, random_state = 42)
  
  df_umap_scar <- data.frame(
    ID    = scar_filter$ID,
    color = as.factor(scar_filter[[color_var]]),
    UMAP1 = umap_scar$layout[, 1],
    UMAP2 = umap_scar$layout[, 2]
  )
  
  # ── 统一颜色（两图用同一套离散色板）───────────────────────────────────────
  all_levels <- union(
    as.character(unique(df_umap_morph$color)),
    as.character(unique(df_umap_scar$color))
  ) %>% sort()
  
  color_palette <- setNames(
    scales::hue_pal()(length(all_levels)),
    all_levels
  )
  
  # ── 单图绘制函数 ────────────────────────────────────────────────────────────
  make_umap_plot <- function(df, subtitle) {
    ggplot(df, aes(x = UMAP1, y = UMAP2, color = color)) +
      geom_point(size = 3.5, alpha = 0.85) +
      geom_text_repel(
        data = df %>% filter(startsWith(as.character(ID), "IM_")),
        aes(label = ID),
        size          = 2.8,
        show.legend   = FALSE,
        max.overlaps  = 20,
        segment.color = "grey60",
        fontface      = "italic"
      ) +
      scale_color_manual(values = color_palette, name = color_title,
                         drop = FALSE) +
      theme_minimal(base_size = 11) +
      labs(subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
      theme(
        plot.subtitle    = element_text(hjust = 0.5, color = "grey30", size = 10),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 8),
        panel.border     = element_rect(color = "black", fill = NA,
                                        linewidth = 0.5),
        panel.grid.minor = element_blank()
      )
  }
  
  p_morph <- make_umap_plot(df_umap_morph, "Overall Morphology (Degrees 1–5)")
  p_scar  <- make_umap_plot(df_umap_scar,  "Scar Direction Pattern (Degrees 1–5)")
  
  # ── 拼图 ────────────────────────────────────────────────────────────────────
  (p_morph + p_scar) +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)
      )
    )
}

# ── 生成三张图 ────────────────────────────────────────────────────────────────
p_layer    <- run_umap_pair(SPHARM_morphology_filter, SPHARM_direction_filter,
                            color_var = "Layer",              color_title = "Layer")
p_coretype <- run_umap_pair(SPHARM_morphology_filter, SPHARM_direction_filter,
                            color_var = "Core_type_Li_merged", color_title = "Core Type")
p_rawmat   <- run_umap_pair(SPHARM_morphology_filter, SPHARM_direction_filter,
                            color_var = "Raw_mat",            color_title = "Raw Material")

print(p_layer)
print(p_coretype)
print(p_rawmat)

ggsave(here("analysis/output/figures/UMAP_by_Layer.png"),
       plot = p_layer,    width = 12, height = 5.5, dpi = 300)
ggsave(here("analysis/output/figures/UMAP_by_CoreType.png"),
       plot = p_coretype, width = 12, height = 5.5, dpi = 300)
ggsave(here("analysis/output/figures/UMAP_by_RawMat.png"),
       plot = p_rawmat,   width = 12, height = 5.5, dpi = 300)


# ==============================================================================
# 5. 保存筛选后数据供下游脚本使用
# mantel_cia.R 会通过 readRDS() 读取，无需重新运行本脚本
# ==============================================================================

saveRDS(SPHARM_direction_filter,
        here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
saveRDS(SPHARM_morphology_filter,
        here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

cat("已保存：SPHARM_direction_filter.rds\n")
cat("已保存：SPHARM_morphology_filter.rds\n")

