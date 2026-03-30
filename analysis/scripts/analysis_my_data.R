library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)
library(htmlwidgets)
library(vegan)
library(ggrepel)
library(readxl)
library(glue)
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

metric_data <- 
  read_xlsx("analysis/data/raw_data/SDG_core_metric.xlsx") %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat)


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
  ) %>%
  left_join(metric_data, by = "ID")

SPHARM_morphology_filter <- SPHARM_morphology %>%
  select(
    ID,   
    SHE,                  
    spectral_entropy,
    power_l1:power_l5  
  ) %>%
  left_join(metric_data, by = "ID")

print(SPHARM_direction_filter)
print(SPHARM_morphology_filter)

set.seed(42)

# ============================================================
# UMAP
# ============================================================
run_umap_pair <- function(morph_filter, scar_filter, color_var, color_title) {
  
  # ── Morphology ──────────────────────────────────────────
  morph_features <- morph_filter %>%
    select(power_l1:power_l5) %>%
    as.matrix()
  
  umap_morph <- umap(morph_features, n_neighbors = 5, random_state = 42)
  
  df_umap_morph <- data.frame(
    ID    = morph_filter$ID,
    color = as.factor(morph_filter[[color_var]]),
    UMAP1 = umap_morph$layout[, 1],
    UMAP2 = umap_morph$layout[, 2]
  )
  
  # ── Scar Direction ───────────────────────────────────────
  scar_features <- scar_filter %>%
    select(power_l1:power_l5) %>%
    as.matrix()
  
  umap_scar <- umap(scar_features, n_neighbors = 5, random_state = 42)
  
  df_umap_scar <- data.frame(
    ID    = scar_filter$ID,
    color = as.factor(scar_filter[[color_var]]),
    UMAP1 = umap_scar$layout[, 1],
    UMAP2 = umap_scar$layout[, 2]
  )
  
  # ── 统一颜色 ─────────────────────────────────────────────
  # 两图用同一套离散色板，保证颜色含义一致
  all_levels <- union(
    as.character(unique(df_umap_morph$color)),
    as.character(unique(df_umap_scar$color))
  ) %>% sort()
  
  color_palette <- setNames(
    scales::hue_pal()(length(all_levels)),
    all_levels
  )
  
  # ── 绘图函数（内部复用）──────────────────────────────────
  make_umap_plot <- function(df, subtitle) {
    ggplot(df, aes(x = UMAP1, y = UMAP2, color = color)) +
      geom_point(size = 3.5, alpha = 0.85) +
      geom_text_repel(
        data     = df %>% filter(startsWith(as.character(ID), "IM_")),
        aes(label = ID),
        size          = 2.8,
        show.legend   = FALSE,
        max.overlaps  = 20,
        segment.color = "grey60",
        fontface      = "italic"
      ) +
      scale_color_manual(
        values = color_palette,
        name   = color_title,
        drop   = FALSE
      ) +
      theme_minimal(base_size = 11) +
      labs(
        subtitle = subtitle,
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      theme(
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 10),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 8),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank()
      )
  }
  
  p_morph <- make_umap_plot(df_umap_morph, "Overall Morphology (Degrees 1–5)")
  p_scar  <- make_umap_plot(df_umap_scar,  "Scar Direction Pattern (Degrees 1–5)")
  
  # ── 拼图 ─────────────────────────────────────────────────
  p_combined <- (p_morph + p_scar) +
    plot_layout(guides = "collect") +         
    plot_annotation(
      theme   = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)
      )
    )
  
  return(p_combined)
}

# ============================================================
# 生成三张图
# ============================================================

p_layer <- run_umap_pair(
  SPHARM_morphology_filter,
  SPHARM_direction_filter,
  color_var   = "Layer",
  color_title = "Layer"
)

p_coretype <- run_umap_pair(
  SPHARM_morphology_filter,
  SPHARM_direction_filter,
  color_var   = "Core_type_Li_merged",
  color_title = "Core Type"
)

p_rawmat <- run_umap_pair(
  SPHARM_morphology_filter,
  SPHARM_direction_filter,
  color_var   = "Raw_mat",
  color_title = "Raw Material"
)

# ============================================================
# 打印 & 保存
# ============================================================
print(p_layer)
print(p_coretype)
print(p_rawmat)

ggsave(
  "analysis/data/derived_data/UMAP_by_Layer.png",
  plot = p_layer, width = 12, height = 5.5, dpi = 300
)
ggsave(
  "analysis/data/derived_data/UMAP_by_CoreType.png",
  plot = p_coretype, width = 12, height = 5.5, dpi = 300
)
ggsave(
  "analysis/data/derived_data/UMAP_by_RawMat.png",
  plot = p_rawmat, width = 12, height = 5.5, dpi = 300
)









