# ==============================================================================
# spharm_analysis.R
# SPHARM 特征分析：方差贡献筛选阶数 + 降维可视化 + 统计检验
#
# 分析流程：
#   1. 读取 SPHARM 功率谱数据（方向 + 形态）及样本元数据
#   2. 绘制各阶方差贡献折线图，辅助判断保留阶数
#   3. 筛选 1–4 阶功率谱特征，划分 EXP+IM / SDG+IM 子集
#   4. UMAP 降维可视化
#   5. EXP 标本：z-score 标准化 + LDA 可视化
#   6. 方向统计量（SPI、Fabric E+I）：KW + Dunn + PERMANOVA
#   7. 两两事后检验汇总热图
#   8. 保存筛选后数据供下游脚本使用
#
# 输入：
#   - analysis/data/derived_data/SPHARM_direction.csv
#   - analysis/data/derived_data/SPHARM_morphology.csv
#   - analysis/data/derived_data/SPHARM_direction_variance_per_degree.csv
#   - analysis/data/derived_data/variance_per_degree.csv
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#   - analysis/data/raw_data/Scar_orientation_data.xlsx
#
# 输出：
#   - analysis/output/figures/Variance_Comparison.png
#   - analysis/output/figures/UMAP_EXP_by_Typology.png
#   - analysis/output/figures/UMAP_SDG_by_Layer.png
#   - analysis/output/figures/UMAP_SDG_by_CoreType.png
#   - analysis/output/figures/UMAP_SDG_by_RawMat.png
#   - analysis/output/figures/LDA_concat_by_Typology.png
#   - analysis/output/figures/LDA_concat_loadings.png
#   - analysis/output/figures/LDA_morph_by_Typology.png
#   - analysis/output/figures/LDA_dir_by_Typology.png
#   - analysis/output/figures/LDA_EI_by_Typology.png
#   - analysis/output/figures/Boxplot_SPI_by_Typology.png
#   - analysis/output/figures/Heatmap_posthoc_summary.png
#   - analysis/data/derived_data/SPHARM_direction_filter.rds
#   - analysis/data/derived_data/SPHARM_morphology_filter.rds
# ==============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)
library(ggrepel)
library(readxl)
library(MASS)
library(vegan)
library(FSA)
library(RVAideMemoire)

conflicted::conflicts_prefer(ggplot2::margin)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)

set.seed(42)

# ==============================================================================
# 全局参数
# ==============================================================================

POWER_COLS      <- paste0("power_l", 1:5)
N_NEIGHBORS     <- 10
EXCLUDE_TYPES   <- c("Biface")
LEVALLOIS_MERGE <- c("Levallois convergent", "Levallois laminar",
                     "Levallois preferential", "Levallois recurrent")

# ==============================================================================
# 1. 读取数据
# ==============================================================================

SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"))
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"))

variance_direction  <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction_variance_per_degree.csv"))
variance_morphology <- read_csv(
  here("analysis/data/derived_data/variance_per_degree.csv"))

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

# ==============================================================================
# 2. 方差贡献折线图
# ==============================================================================

plot_spharm_variance <- function(df_dir, df_mor) {
  bind_rows(
    df_dir %>% mutate(Feature = "Direction"),
    df_mor %>% mutate(Feature = "Morphology")
  ) %>%
    ggplot(aes(x = degree, y = variance, color = Feature, shape = Feature)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    scale_x_continuous(breaks = function(x) seq(min(x), max(x), by = 1)) +
    scale_color_manual(values = c("Direction"  = "#FFBAE0",
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
ggsave(here("analysis/output/figures/Variance_Comparison.png"),
       plot = variance_plot, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 3. 筛选特征 + 划分子集
# ==============================================================================

filter_spharm <- function(df, meta = NULL) {
  result <- df %>%
    select(ID, Typology, SHE, spectral_entropy, all_of(POWER_COLS))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}

SPHARM_direction_filter  <- filter_spharm(SPHARM_direction,  metric_data)
SPHARM_morphology_filter <- filter_spharm(SPHARM_morphology, metric_data)

split_by_group <- function(df) {
  list(
    exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
    sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_"))
  )
}

dir_splits   <- split_by_group(SPHARM_direction_filter)
morph_splits <- split_by_group(SPHARM_morphology_filter)

# ==============================================================================
# 4. UMAP 可视化
# ==============================================================================

compute_umap <- function(df, n_neighbors = N_NEIGHBORS, seed = 42) {
  features <- df %>% select(all_of(POWER_COLS)) %>% as.matrix()
  result   <- umap(features, n_neighbors = n_neighbors, random_state = seed)
  df %>% mutate(
    UMAP1 = result$layout[, 1],
    UMAP2 = result$layout[, 2],
    is_IM = str_starts(ID, "IM_")
  )
}

make_umap_plot <- function(df_umap, color_var, color_title, subtitle) {
  df_im  <- df_umap %>% filter(is_IM)
  df_reg <- df_umap %>% filter(!is_IM)
  
  color_levels  <- sort(unique(as.character(df_reg[[color_var]])))
  color_palette <- setNames(scales::hue_pal()(length(color_levels)), color_levels)
  
  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
    geom_point(
      data  = df_reg,
      aes(color = .data[[color_var]]),
      size = 3, alpha = 0.85, shape = 16
    ) +
    geom_point(
      data  = df_im,
      shape = 17, size = 3.5, color = "grey30", alpha = 0.9
    ) +
    geom_text_repel(
      data         = df_im,
      aes(label    = ID %>%
            str_remove("^IM_") %>%
            str_replace_all("_", " ") %>%
            str_to_sentence()),
      size          = 2.5,
      color         = "grey30",
      fontface      = "italic",
      max.overlaps  = 20,
      segment.color = "grey60"
    ) +
    scale_color_manual(values = color_palette, name = color_title) +
    theme_minimal(base_size = 11) +
    labs(subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    theme(
      plot.subtitle    = element_text(hjust = 0.5, color = "grey30", size = 10),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 8),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank()
    )
}

make_umap_pair <- function(df_morph_umap, df_dir_umap, color_var, color_title) {
  p_morph <- make_umap_plot(df_morph_umap, color_var, color_title,
                            "Overall Morphology (Degrees 1–4)")
  p_dir   <- make_umap_plot(df_dir_umap,   color_var, color_title,
                            "Scar Direction Pattern (Degrees 1–4)")
  (p_morph + p_dir) +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
      )
    )
}

# EXP + IM
umap_exp_dir   <- compute_umap(dir_splits$exp_im)
umap_exp_morph <- compute_umap(morph_splits$exp_im)

p_exp <- make_umap_pair(umap_exp_morph, umap_exp_dir,
                        color_var = "Typology", color_title = "Core Type")
print(p_exp)
ggsave(here("analysis/output/figures/UMAP_EXP_by_Typology.png"),
       plot = p_exp, width = 12, height = 5.5, dpi = 300)

# SDG + IM
umap_sdg_dir   <- compute_umap(dir_splits$sdg_im)
umap_sdg_morph <- compute_umap(morph_splits$sdg_im)

sdg_plots <- list(
  Layer    = make_umap_pair(umap_sdg_morph, umap_sdg_dir,
                            "Layer",               "Layer"),
  CoreType = make_umap_pair(umap_sdg_morph, umap_sdg_dir,
                            "Core_type_Li_merged", "Core Type"),
  RawMat   = make_umap_pair(umap_sdg_morph, umap_sdg_dir,
                            "Raw_mat",             "Raw Material")
)

walk(sdg_plots, print)
ggsave(here("analysis/output/figures/UMAP_SDG_by_Layer.png"),
       plot = sdg_plots$Layer,    width = 12, height = 5.5, dpi = 300)
ggsave(here("analysis/output/figures/UMAP_SDG_by_CoreType.png"),
       plot = sdg_plots$CoreType, width = 12, height = 5.5, dpi = 300)
ggsave(here("analysis/output/figures/UMAP_SDG_by_RawMat.png"),
       plot = sdg_plots$RawMat,   width = 12, height = 5.5, dpi = 300)

# ==============================================================================
# 5. EXP 标本：z-score 标准化 + LDA 可视化
# ==============================================================================

df_exp_dir   <- dir_splits$exp_im
df_exp_morph <- morph_splits$exp_im

scale_features <- function(df_target, cols = POWER_COLS) {
  ref_mat  <- df_target %>%
    filter(!str_starts(ID, "IM_")) %>%
    select(all_of(cols)) %>% as.matrix()
  col_mean <- colMeans(ref_mat)
  col_sd   <- apply(ref_mat, 2, sd)
  mat      <- df_target %>% select(all_of(cols)) %>% as.matrix()
  scale(mat, center = col_mean, scale = col_sd)
}

z_dir   <- scale_features(df_exp_dir)
z_morph <- scale_features(df_exp_morph)
colnames(z_dir)   <- paste0("dir_",   POWER_COLS)
colnames(z_morph) <- paste0("morph_", POWER_COLS)
z_concat <- cbind(z_morph, z_dir)

df_exp_only <- df_exp_dir %>%
  filter(!str_starts(ID, "IM_"),
         !Typology %in% EXCLUDE_TYPES) %>%
  mutate(
    Typology = case_when(
      Typology %in% LEVALLOIS_MERGE ~ "Levallois",
      TRUE ~ Typology
    ),
    Typology = droplevels(as.factor(Typology))
  )

non_im_idx <- !str_starts(df_exp_dir$ID, "IM_") &
  !df_exp_dir$Typology %in% EXCLUDE_TYPES

X_concat   <- z_concat[non_im_idx, ]
y_typology <- df_exp_only$Typology

cat(sprintf("EXP 保留标本数: %d，类别数: %d\n",
            nrow(df_exp_only), nlevels(y_typology)))
print(table(y_typology))

run_lda_plot <- function(X, y, ids, title_str, subtitle_str, filename) {
  fit      <- lda(X, grouping = y)
  prop_var <- fit$svd^2 / sum(fit$svd^2)
  
  scores <- predict(fit)$x %>%
    as.data.frame() %>%
    mutate(Typology = y, ID = ids)
  
  hull_df <- scores %>%
    group_by(Typology) %>%
    slice(chull(LD1, LD2)) %>%
    ungroup()
  
  p <- ggplot(scores, aes(x = LD1, y = LD2, color = Typology)) +
    geom_polygon(data = hull_df,
                 aes(fill = Typology, group = Typology),
                 alpha = 0.15, color = NA) +
    geom_polygon(data = hull_df,
                 aes(color = Typology, group = Typology),
                 fill = NA, linewidth = 0.6) +
    geom_point(size = 3, alpha = 0.85) +
    theme_minimal(base_size = 11) +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x = sprintf("LD1 (%.1f%%)", prop_var[1] * 100),
      y = sprintf("LD2 (%.1f%%)", prop_var[2] * 100)
    ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  print(p)
  ggsave(here("analysis/output/figures", filename),
         plot = p, width = 7, height = 5.5, dpi = 300)
  
  invisible(list(fit = fit, scores = scores, prop_var = prop_var))
}

# LDA 1：拼接谱
res_concat <- run_lda_plot(
  X = X_concat, y = y_typology, ids = df_exp_only$ID,
  title_str    = "LDA: EXP Specimens by Typology",
  subtitle_str = "Concatenated SPHARM (Morphology + Direction, Degrees 1–4)",
  filename     = "LDA_concat_by_Typology.png"
)

# 拼接谱载荷图
as.data.frame(res_concat$fit$scaling) %>%
  rownames_to_column("Feature") %>%
  mutate(Type = ifelse(str_starts(Feature, "morph"), "Morphology", "Direction")) %>%
  ggplot(aes(x = LD1, y = LD2, label = Feature, color = Type)) +
  geom_segment(aes(xend = 0, yend = 0),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "LDA Feature Loadings — Concatenated SPHARM")
ggsave(here("analysis/output/figures/LDA_concat_loadings.png"),
       width = 7, height = 6, dpi = 300)

# LDA 2：形态谱
res_morph <- run_lda_plot(
  X = z_morph[non_im_idx, ], y = y_typology, ids = df_exp_only$ID,
  title_str    = "LDA: EXP Specimens by Typology",
  subtitle_str = "SPHARM Morphology only (Degrees 1–4)",
  filename     = "LDA_morph_by_Typology.png"
)

# LDA 3：方向谱
res_dir <- run_lda_plot(
  X = z_dir[non_im_idx, ], y = y_typology, ids = df_exp_only$ID,
  title_str    = "LDA: EXP Specimens by Typology",
  subtitle_str = "SPHARM Direction only (Degrees 1–4)",
  filename     = "LDA_dir_by_Typology.png"
)

# ==============================================================================
# 6. 方向统计量：SPI + Fabric (E+I)
# ==============================================================================

compute_SPI <- function(ux, uy, uz) {
  sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2) / length(ux)
}

compute_EI <- function(ux, uy, uz) {
  mat <- cbind(ux, uy, uz)
  S   <- t(mat) %*% mat / nrow(mat)
  eig <- sort(eigen(S, symmetric = TRUE)$values, decreasing = TRUE)
  list(
    E       = 1 - eig[2] / eig[1],
    I       = 1 - eig[3] / eig[2],
    lambda1 = eig[1], lambda2 = eig[2], lambda3 = eig[3]
  )
}

raw_dirs <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 3
) %>%
  mutate(
    dx     = End_X - Start_X,
    dy     = End_Y - Start_Y,
    dz     = End_Z - Start_Z,
    length = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(length > 1e-10) %>%
  mutate(
    ux = dx / length,
    uy = dy / length,
    uz = dz / length
  )

results <- raw_dirs %>%
  group_by(ID) %>%
  summarise(
    n_scars = n(),
    SPI     = compute_SPI(ux, uy, uz),
    EI      = list(compute_EI(ux, uy, uz)),
    .groups = "drop"
  ) %>%
  mutate(
    E       = map_dbl(EI, "E"),
    I       = map_dbl(EI, "I"),
    lambda1 = map_dbl(EI, "lambda1"),
    lambda2 = map_dbl(EI, "lambda2"),
    lambda3 = map_dbl(EI, "lambda3")
  ) %>%
  select(-EI) %>%
  arrange(ID)

cat("SPI/E/I 计算完成，标本数：", nrow(results), "\n")
print(results %>%
        select(ID, n_scars, SPI, E, I) %>%
        mutate(across(c(SPI, E, I), \(x) round(x, 4))),
      n = Inf)

results_typed <- results %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID") %>%
  filter(str_starts(ID, "EXP"),
         !Typology %in% EXCLUDE_TYPES) %>%
  mutate(
    Typology = case_when(
      Typology %in% LEVALLOIS_MERGE ~ "Levallois",
      TRUE ~ Typology
    ),
    Typology = droplevels(as.factor(Typology))
  ) %>%
  filter(complete.cases(SPI, E, I))

y_rei <- results_typed$Typology

cat(sprintf("方向统计量可用标本: %d，类别数: %d\n",
            nrow(results_typed), nlevels(y_rei)))
print(table(y_rei))

# --- SPI：Kruskal-Wallis + Dunn（Holm）+ 箱线图 ---
kw_SPI <- kruskal.test(SPI ~ Typology, data = results_typed)
cat("\n--- Kruskal-Wallis 检验（SPI）---\n")
print(kw_SPI)

cat("\n--- Dunn 事后检验（Holm 校正）：SPI ---\n")
dunn_SPI <- FSA::dunnTest(SPI ~ Typology, data = results_typed, method = "holm")
print(dunn_SPI)

p_SPI_box <- ggplot(results_typed,
                    aes(x = Typology, y = SPI, fill = Typology)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6) +
  annotate("text",
           x     = Inf, y = Inf,
           label = sprintf("Kruskal-Wallis\nH = %.2f, df = %d, p < 0.001",
                           kw_SPI$statistic, kw_SPI$parameter),
           hjust = 1.05, vjust = 1.2, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(title = "SPI (Scar Pattern Index) by Typology", x = NULL, y = "SPI")

print(p_SPI_box)
ggsave(here("analysis/output/figures/Boxplot_SPI_by_Typology.png"),
       plot = p_SPI_box, width = 6, height = 5, dpi = 300)

# --- PERMANOVA 辅助函数 ---
run_permanova <- function(X, group_vec, label) {
  df_grp <- data.frame(Typology = group_vec)
  
  set.seed(42)
  perm_global <- adonis2(X ~ Typology, data = df_grp,
                         method = "euclidean", permutations = 999)
  cat(sprintf("\n--- PERMANOVA：%s ---\n", label))
  print(perm_global)
  
  set.seed(42)
  pairwise_res <- pairwise.perm.manova(
    dist(X, method = "euclidean"),
    group_vec,
    nperm    = 999,
    p.method = "holm"
  )
  cat(sprintf("\n--- 两两 PERMANOVA 事后检验（Holm）：%s ---\n", label))
  print(pairwise_res$p.value)
  
  invisible(list(global = perm_global, pairwise = pairwise_res))
}

# --- Fabric (E+I)：PERMANOVA + LDA ---
X_EI    <- results_typed %>% select(E, I) %>% as.matrix()
perm_EI <- run_permanova(X_EI, y_rei, "Fabric (E+I)")

lda_fit_EI  <- lda(X_EI, grouping = y_rei)
prop_var_EI <- lda_fit_EI$svd^2 / sum(lda_fit_EI$svd^2)

scores_EI <- predict(lda_fit_EI)$x %>%
  as.data.frame() %>%
  mutate(Typology = y_rei, ID = results_typed$ID)

n_ld <- ncol(scores_EI %>% select(starts_with("LD")))

if (n_ld == 1) {
  p_ei <- ggplot(scores_EI,
                 aes(x = LD1, y = Typology, fill = Typology)) +
    ggridges::geom_density_ridges(alpha = 0.6, scale = 0.9) +
    theme_minimal(base_size = 11) +
    labs(
      title    = "LDA: EXP Specimens by Typology",
      subtitle = "Fabric (E + I)",
      x        = sprintf("LD1 (%.1f%%)", prop_var_EI[1] * 100),
      y        = NULL
    ) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")
} else {
  hull_EI <- scores_EI %>%
    group_by(Typology) %>%
    slice(chull(LD1, LD2)) %>%
    ungroup()
  
  p_ei <- ggplot(scores_EI, aes(x = LD1, y = LD2, color = Typology)) +
    geom_polygon(data = hull_EI,
                 aes(fill = Typology, group = Typology),
                 alpha = 0.15, color = NA) +
    geom_polygon(data = hull_EI,
                 aes(color = Typology, group = Typology),
                 fill = NA, linewidth = 0.6) +
    geom_point(size = 3, alpha = 0.85) +
    theme_minimal(base_size = 11) +
    labs(
      title    = "LDA: EXP Specimens by Typology",
      subtitle = "Fabric (E + I)",
      x = sprintf("LD1 (%.1f%%)", prop_var_EI[1] * 100),
      y = sprintf("LD2 (%.1f%%)", prop_var_EI[2] * 100)
    ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

print(p_ei)
ggsave(here("analysis/output/figures/LDA_EI_by_Typology.png"),
       plot = p_ei, width = 7, height = 5.5, dpi = 300)

# --- SPHARM PERMANOVA ---
perm_concat <- run_permanova(X_concat,             y_typology, "SPHARM 拼接")
perm_morph  <- run_permanova(z_morph[non_im_idx,], y_typology, "SPHARM 形态谱")
perm_dir    <- run_permanova(z_dir[non_im_idx,],   y_typology, "SPHARM 方向谱")

cat("\n========== PERMANOVA 汇总 ==========\n")
cat(sprintf("SPHARM 拼接 : R² = %.3f, p = %.3f\n",
            perm_concat$global$R2[1], perm_concat$global$`Pr(>F)`[1]))
cat(sprintf("SPHARM 形态 : R² = %.3f, p = %.3f\n",
            perm_morph$global$R2[1],  perm_morph$global$`Pr(>F)`[1]))
cat(sprintf("SPHARM 方向 : R² = %.3f, p = %.3f\n",
            perm_dir$global$R2[1],    perm_dir$global$`Pr(>F)`[1]))
cat(sprintf("Fabric (E+I): R² = %.3f, p = %.3f\n",
            perm_EI$global$R2[1],     perm_EI$global$`Pr(>F)`[1]))
cat(sprintf("SPI KW      : H  = %.2f, df = %d, p < 0.001\n",
            kw_SPI$statistic, kw_SPI$parameter))

# ==============================================================================
# 7. 热图：两两事后检验汇总
# ==============================================================================

pairs_full <- c(
  "Bidirectional – Discoid",
  "Bidirectional – Levallois",
  "Bidirectional – Multiplatform",
  "Bidirectional – Unidirectional",
  "Discoid – Levallois",
  "Discoid – Multiplatform",
  "Discoid – Unidirectional",
  "Levallois – Multiplatform",
  "Levallois – Unidirectional",
  "Multiplatform – Unidirectional"
)

method_levels <- c(
  "Morphology\nSPHARM",
  "Direction\nSPHARM",
  "Concatenated\nSPHARM",
  "Fabric\n(E+I)",
  "SPI\n(Dunn)"
)

p_adj_vec <- c(
  # Morphology SPHARM — 全部 ns
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  
  # Direction SPHARM（Holm）
  0.788, 0.010, 0.024, 0.010,
  0.024, 0.316, 0.027,
  0.010, 0.024,
  0.014,
  
  # Concatenated SPHARM（Holm）
  0.884, 0.010, 0.140, 0.064,
  0.063, 0.838, 0.140,
  0.144, 0.441,
  0.140,
  
  # Fabric E+I（Holm）
  0.030, 0.010, 0.224, 0.080,
  0.243, 0.586, 0.010,
  0.761, 0.010,
  0.010,
  
  # SPI Dunn（Holm）
  # Bidir–Dis, Bidir–Lev, Bidir–Multi, Bidir–Uni
  1.0000, 0.6534, 0.9665, 0.0058,
  # Dis–Lev, Dis–Multi, Dis–Uni
  0.1455, 0.1572, 0.0001,
  # Lev–Multi, Lev–Uni
  0.8816, 0.0040,
  # Multi–Uni
  0.0421
)

plot_data <- tibble(
  Pair   = factor(rep(pairs_full, 5), levels = rev(pairs_full)),
  Method = factor(rep(method_levels, each = length(pairs_full)),
                  levels = method_levels),
  p_adj  = p_adj_vec
) %>%
  mutate(
    Sig = case_when(
      p_adj <= 0.001 ~ "p ≤ 0.001 ***",
      p_adj <= 0.01  ~ "p ≤ 0.01 **",
      p_adj <= 0.05  ~ "p ≤ 0.05 *",
      TRUE           ~ "ns"
    ),
    Sig = factor(Sig, levels = c(
      "p ≤ 0.001 ***",
      "p ≤ 0.01 **",
      "p ≤ 0.05 *",
      "ns"
    )),
    Label = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

sig_colors <- c(
  "p ≤ 0.001 ***" = "#0F6E56",
  "p ≤ 0.01 **"   = "#1D9E75",
  "p ≤ 0.05 *"    = "#97C459",
  "ns"            = "#E8E6DF"
)

label_colors <- c(
  "p ≤ 0.001 ***" = "#ffffff",
  "p ≤ 0.01 **"   = "#ffffff",
  "p ≤ 0.05 *"    = "#27500A",
  "ns"            = "#888780"
)

p_heatmap <- ggplot(plot_data, aes(x = Method, y = Pair, fill = Sig)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Label, color = Sig),
            size = 3.8, fontface = "bold") +
  scale_fill_manual(values = sig_colors,
                    name = NULL) +
  scale_color_manual(values = label_colors, guide = "none") +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0, size = 12),
    plot.subtitle   = element_text(color = "grey50", size = 8.5,
                                   margin = margin(b = 10)),
    plot.caption    = element_text(color = "grey60", size = 8,
                                   hjust = 0, margin = margin(t = 8)),
    axis.text.x     = element_text(face = "bold", size = 9.5,
                                   lineheight = 1.1),
    axis.text.y     = element_text(size = 9.5),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(0.55, "cm"),
    plot.margin     = margin(10, 10, 10, 10)
  ) +
  labs(
    title    = "Pairwise post-hoc significance across methods (Holm–Bonferroni)"
  )

p_heatmap
ggsave(here("analysis/output/figures/Heatmap_posthoc_summary.png"),
       plot = p_heatmap, width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# 8. 保存筛选后数据供下游脚本使用
# ==============================================================================

saveRDS(SPHARM_direction_filter,
        here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
saveRDS(SPHARM_morphology_filter,
        here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

cat("已保存：SPHARM_direction_filter.rds\n")
cat("已保存：SPHARM_morphology_filter.rds\n")