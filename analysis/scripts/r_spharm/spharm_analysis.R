# ==============================================================================
# spharm_analysis.R
# SPHARM 特征分析：方差贡献筛选阶数 + UMAP 降维可视化
#
# 分析流程：
#   1. 读取 SPHARM 功率谱数据（方向 + 形态）及样本元数据
#   2. 绘制各阶方差贡献折线图，辅助判断保留阶数
#   3. 筛选 1–4 阶功率谱特征，划分 EXP+IM / SDG+IM 子集
#   4. UMAP 降维：
#        EXP+IM → 按 Typology 着色
#        SDG+IM → 按 Layer / Core type / Raw material 分别着色
#      两套图中 IM 标本均以参考点形式呈现
#   5. 保存筛选后数据供 mantel_cia.R 使用
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
#   - analysis/output/figures/UMAP_EXP_by_Typology.png
#   - analysis/output/figures/UMAP_SDG_by_Layer.png
#   - analysis/output/figures/UMAP_SDG_by_CoreType.png
#   - analysis/output/figures/UMAP_SDG_by_RawMat.png
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

set.seed(42)

# ==============================================================================
# 全局参数
# ==============================================================================

POWER_COLS  <- paste0("power_l", 1:4)   
N_NEIGHBORS <- 10                        


# ==============================================================================
# 1. 读取数据
# ==============================================================================

SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"))
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"))

variance_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction_variance_per_degree.csv"))
variance_morphology <- read_csv(here("analysis/data/derived_data/variance_per_degree.csv"))

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

# ==============================================================================
# 2. 方差贡献折线图
# ==============================================================================

plot_spharm_variance <- function(df_dir, df_mor) {
  df_combined <- bind_rows(
    df_dir %>% mutate(Feature = "Direction"),
    df_mor %>% mutate(Feature = "Morphology")
  )
  
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
ggsave(here("analysis/output/figures/Variance_Comparison.png"),
       plot = variance_plot, width = 8, height = 6, dpi = 300)


# ==============================================================================
# 3. 筛选特征 + 划分子集
# ==============================================================================

# 基础筛选：保留 ID、Typology、熵值、功率谱列，连接元数据
filter_spharm <- function(df, meta = NULL) {
  result <- df %>%
    select(ID, Typology, SHE, spectral_entropy, all_of(POWER_COLS))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}

SPHARM_direction_filter  <- filter_spharm(SPHARM_direction,  metric_data)
SPHARM_morphology_filter <- filter_spharm(SPHARM_morphology, metric_data)

# 划分子集
#   EXP+IM：实验标本 + 理想模型（EXP 前缀不含下划线，如 EXP01_）
#   SDG+IM：考古标本 + 理想模型
split_by_group <- function(df) {
  list(
    exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
    sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_"))
  )
}

dir_splits  <- split_by_group(SPHARM_direction_filter)
morph_splits <- split_by_group(SPHARM_morphology_filter)


# ==============================================================================
# 4. UMAP 工具函数
# ==============================================================================

# --- 4-1: 计算 UMAP，附加坐标和分组标记 ---
compute_umap <- function(df, n_neighbors = N_NEIGHBORS, seed = 42) {
  features <- df %>% select(all_of(POWER_COLS)) %>% as.matrix()
  result   <- umap(features, n_neighbors = n_neighbors, random_state = seed)
  
  df %>% mutate(
    UMAP1 = result$layout[, 1],
    UMAP2 = result$layout[, 2],
    is_IM = str_starts(ID, "IM_")
  )
}

# --- 4-2: 单张 UMAP 图 ---
# 普通标本按 color_var 着色（圆点），IM 标本灰色三角形并标注标签
make_umap_plot <- function(df_umap, color_var, color_title, subtitle) {
  
  df_im  <- df_umap %>% filter(is_IM)
  df_reg <- df_umap %>% filter(!is_IM)
  
  # 从非 IM 标本构建颜色映射
  color_levels  <- sort(unique(as.character(df_reg[[color_var]])))
  color_palette <- setNames(scales::hue_pal()(length(color_levels)), color_levels)
  
  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
    # 普通标本：圆点，按变量着色
    geom_point(
      data  = df_reg,
      aes(color = .data[[color_var]]),
      size = 3, alpha = 0.85, shape = 16
    ) +
    # IM 标本：灰色三角形
    geom_point(
      data  = df_im,
      shape = 17, size = 3.5, color = "grey30", alpha = 0.9
    ) +
    # IM 标本标签
    geom_text_repel(
      data          = df_im,
      aes(label     = ID %>%
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

# --- 4-3: 形态 + 方向并排拼图 ---
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


# ==============================================================================
# 5. 生成 UMAP 图
# ==============================================================================

# --- EXP + IM：按 Typology 着色（59 EXP + 12 IM = 71 标本）---
umap_exp_dir   <- compute_umap(dir_splits$exp_im)
umap_exp_morph <- compute_umap(morph_splits$exp_im)

p_exp <- make_umap_pair(umap_exp_morph, umap_exp_dir,
                        color_var   = "Typology",
                        color_title = "Core Type")
print(p_exp)
ggsave(here("analysis/output/figures/UMAP_EXP_by_Typology.png"),
       plot = p_exp, width = 12, height = 5.5, dpi = 300)


# --- SDG + IM：按三种元数据着色（55 SDG + 12 IM = 67 标本）---
umap_sdg_dir   <- compute_umap(dir_splits$sdg_im)
umap_sdg_morph <- compute_umap(morph_splits$sdg_im)

p_sdg_layer    <- make_umap_pair(umap_sdg_morph, umap_sdg_dir,
                                 color_var   = "Layer",
                                 color_title = "Layer")
p_sdg_coretype <- make_umap_pair(umap_sdg_morph, umap_sdg_dir,
                                 color_var   = "Core_type_Li_merged",
                                 color_title = "Core Type")
p_sdg_rawmat   <- make_umap_pair(umap_sdg_morph, umap_sdg_dir,
                                 color_var   = "Raw_mat",
                                 color_title = "Raw Material")

print(p_sdg_layer)
print(p_sdg_coretype)
print(p_sdg_rawmat)

ggsave(here("analysis/output/figures/UMAP_SDG_by_Layer.png"),
       plot = p_sdg_layer,    width = 12, height = 5.5, dpi = 300)
ggsave(here("analysis/output/figures/UMAP_SDG_by_CoreType.png"),
       plot = p_sdg_coretype, width = 12, height = 5.5, dpi = 300)
ggsave(here("analysis/output/figures/UMAP_SDG_by_RawMat.png"),
       plot = p_sdg_rawmat,   width = 12, height = 5.5, dpi = 300)


# ==============================================================================
# 5b. EXP 标本：z-score 标准化 + 拼接LDA
# ==============================================================================

library(MASS)

# --- 从 splits 里取出 EXP+IM 子集（这两个对象才是脚本里实际存在的）---
df_exp_dir   <- dir_splits$exp_im     # 补上这行
df_exp_morph <- morph_splits$exp_im   # 补上这行

# --- z-score 标准化（以 EXP 标本的均值/标准差为基准）---
scale_features <- function(df_target, cols = POWER_COLS) {
  ref_mat  <- df_target %>%
    dplyr::filter(!str_starts(ID, "IM_")) %>%
    dplyr::select(all_of(cols)) %>% as.matrix()
  col_mean <- colMeans(ref_mat)
  col_sd   <- apply(ref_mat, 2, sd)
  mat      <- df_target %>% 
    dplyr::select(dplyr::all_of(cols)) %>% as.matrix()
  scale(mat, center = col_mean, scale = col_sd)
}

z_dir   <- scale_features(df_exp_dir)
z_morph <- scale_features(df_exp_morph)

colnames(z_dir)   <- paste0("dir_",   POWER_COLS)
colnames(z_morph) <- paste0("morph_", POWER_COLS)
z_concat <- cbind(z_morph, z_dir)    # 补上这行

# --- 之后接你现有的过滤/合并/LDA 代码 ---
EXCLUDE_TYPES   <- c("Multiplatform", "Biface","Unidirectional", "Bidirectional", "Discoid")
# LEVALLOIS_MERGE <- c("Levallois convergent", "Levallois laminar",
#                      "Levallois preferential", "Levallois recurrent")

df_exp_only <- df_exp_dir %>%
  filter(!str_starts(ID, "IM_"),
         !Typology %in% EXCLUDE_TYPES) 
# %>%
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

cat(sprintf("保留标本数: %d，类别数: %d\n",
            nrow(df_exp_only), nlevels(y_typology)))
print(table(y_typology))

# --- LDA ---
lda_fit <- lda(X_concat, grouping = y_typology)

# 查看各判别轴解释的方差比例
prop_var <- lda_fit$svd^2 / sum(lda_fit$svd^2)
cat("各LD轴方差贡献：\n")
print(round(prop_var, 3))

# 投影到 LD 空间
lda_scores <- predict(lda_fit)$x %>%
  as.data.frame() %>%
  mutate(Typology = y_typology,
         ID = df_exp_only$ID)

hull_df <- lda_scores %>%
  group_by(Typology) %>%
  slice(chull(LD1, LD2)) %>%   # 关键：计算凸包点
  ungroup()

# --- 可视化 LD1 vs LD2 ---
ggplot(lda_scores, aes(x = LD1, y = LD2, color = Typology)) +
  geom_polygon(data = hull_df,
               aes(x = LD1, y = LD2, fill = Typology, group = Typology),
               alpha = 0.2, color = NA) +
  geom_polygon(data = hull_df,
               aes(x = LD1, y = LD2, color = Typology, group = Typology),
               fill = NA, linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.85) +
  # 数据点标签
  geom_text_repel(aes(label = ID),
                  size = 2.5,
                  max.overlaps = 50,
                  segment.color = "grey60") +
  theme_minimal(base_size = 11) +
  labs(title = "LDA: EXP Specimens by Typology",
       subtitle = "Concatenated SPHARM features (Morphology + Direction, Degrees 1–4)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- 留一法交叉验证（LOO-CV）：最诚实的小样本评估 ---
lda_cv  <- lda(X_concat, grouping = y_typology, CV = TRUE)
conf_mat <- table(Predicted = lda_cv$class, Actual = y_typology)
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
cat(sprintf("\nLOO-CV 准确率: %.1f%%\n", accuracy * 100))
print(conf_mat)



# 载荷矩阵：哪些特征对 LD1/LD2 贡献最大
loadings <- lda_fit$scaling
print(round(loadings[, 1:2], 3))

# 可视化载荷
as.data.frame(loadings) %>%
  rownames_to_column("Feature") %>%
  mutate(Type = ifelse(str_starts(Feature, "morph"), "Morphology", "Direction")) %>%
  ggplot(aes(x = LD1, y = LD2, label = Feature, color = Type)) +
  geom_segment(aes(xend = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")),
               alpha = 0.7) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "LDA Feature Loadings")

# 只用形态谱
lda_cv_morph <- lda(z_morph[non_im_idx, ], grouping = y_typology, CV = TRUE)
acc_morph <- mean(lda_cv_morph$class == y_typology)

# 只用方向谱
lda_cv_dir <- lda(z_dir[non_im_idx, ], grouping = y_typology, CV = TRUE)
acc_dir <- mean(lda_cv_dir$class == y_typology)

cat(sprintf("形态谱: %.1f%% | 方向谱: %.1f%% | 拼接: 50.0%%\n",
            acc_morph * 100, acc_dir * 100))





library(ade4)

# ==============================================================================
# CoIA 分析
# ==============================================================================

# --- 1. 准备两个矩阵（只用 EXP 标本，排除 IM 和 Multiplatform/Biface）---
X_morph <- z_morph[non_im_idx, ]
X_dir   <- z_dir[non_im_idx, ]

# --- 2. 各自做 PCA（dudi.pca 要求已标准化，我们的 z-score 矩阵满足）---
pca_morph <- dudi.pca(as.data.frame(X_morph), 
                      center = FALSE, scale = FALSE,
                      scannf = FALSE, nf = 3)   # 保留3个PC
pca_dir   <- dudi.pca(as.data.frame(X_dir),   
                      center = FALSE, scale = FALSE,
                      scannf = FALSE, nf = 3)

# 检查各自方差解释
cat("形态谱 PCA 方差贡献：\n")
print(round(pca_morph$eig / sum(pca_morph$eig), 3))
cat("方向谱 PCA 方差贡献：\n")
print(round(pca_dir$eig / sum(pca_dir$eig), 3))

# --- 3. Co-Inertia Analysis ---
coia <- coinertia(pca_morph, pca_dir,
                  scannf = FALSE, nf = 2)

# RV 系数：衡量两组数据的整体相关程度
cat(sprintf("\nRV 系数: %.4f\n", coia$RV))

# --- 4. 置换检验（999次，检验共惯量是否显著）---
set.seed(42)
perm_test <- randtest(coia, nrepet = 999)
print(perm_test)
plot(perm_test)   # 置换分布图

# --- 5. 提取公共空间坐标用于后续分类 ---
# coia$lX：形态谱在公共空间的投影
# coia$lY：方向谱在公共空间的投影
# 两者取平均作为"共同坐标"送入 LDA
coia_scores <- (coia$lX + coia$lY) / 2
colnames(coia_scores) <- paste0("CoIA", 1:ncol(coia_scores))

cat("\nCoIA 公共空间坐标维度：\n")
print(dim(coia_scores))

# --- 6. 用 CoIA 坐标做 LDA ---
lda_coia    <- lda(coia_scores, grouping = y_typology)
lda_coia_cv <- lda(coia_scores, grouping = y_typology, CV = TRUE)

acc_coia <- mean(lda_coia_cv$class == y_typology)
cat(sprintf("\nCoIA + LDA LOO-CV 准确率: %.1f%%\n", acc_coia * 100))

conf_coia <- table(Predicted = lda_coia_cv$class, Actual = y_typology)
print(round(prop.table(conf_coia, margin = 2) * 100, 1))

# --- 7. 可视化：CoIA 双标图 ---
# 两组特征各自的样本投影
coia_plot_df <- data.frame(
  morph1 = coia$lX[, 1], morph2 = coia$lX[, 2],
  dir1   = coia$lY[, 1], dir2   = coia$lY[, 2],
  Typology = y_typology
)

# 用线段连接同一标本在两个空间的投影（越短说明两组特征越一致）
ggplot(coia_plot_df) +
  geom_segment(aes(x = morph1, y = morph2,
                   xend = dir1, yend = dir2,
                   color = Typology),
               alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(x = morph1, y = morph2, color = Typology),
             shape = 16, size = 2.5) +
  geom_point(aes(x = dir1,   y = dir2,   color = Typology),
             shape = 17, size = 2.5) +
  theme_minimal(base_size = 11) +
  labs(title    = "Co-Inertia Analysis: Morphology vs Direction",
       subtitle = "圆点 = 形态谱投影，三角 = 方向谱投影，线段越短两组越一致",
       x = "Axis 1", y = "Axis 2") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- 8. 汇总对比所有方案 ---
cat("\n========== 方法对比汇总 ==========\n")
cat(sprintf("随机基准:       %.1f%%\n", 100 / nlevels(y_typology)))
cat(sprintf("形态谱 alone:   %.1f%%\n", acc_morph * 100))
cat(sprintf("方向谱 alone:   %.1f%%\n", acc_dir   * 100))
cat(sprintf("cbind 拼接:     %.1f%%\n", acc_coia  * 100))  # 你之前的结果
cat(sprintf("CoIA 公共空间:  %.1f%%\n", acc_coia  * 100))





# ==============================================================================
# 5c. 方向统计量（R、E、I）→ LDA 分类
# ==============================================================================

library(readxl)

# --- 辅助函数：计算结果向量长度 R ---
compute_R <- function(ux, uy, uz) {
  n <- length(ux)
  sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2) / n
}

# --- 辅助函数：计算 Fabric 椭球指数 E、I ---
compute_EI <- function(ux, uy, uz) {
  mat <- cbind(ux, uy, uz)
  S   <- t(mat) %*% mat / nrow(mat)          # 方向张量
  eig <- sort(eigen(S, symmetric = TRUE)$values, decreasing = TRUE)
  lambda1 <- eig[1]; lambda2 <- eig[2]; lambda3 <- eig[3]
  E <- 1 - lambda2 / lambda1                 # 椭圆度
  I <- 1 - lambda3 / lambda2                 # 倾斜度（扁率）
  list(E = E, I = I,
       lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
}

# --- 1. 读取片疤方向原始数据 ---
raw <- read_excel(here("analysis/data/raw_data/Scar_orientation_data.xlsx"),
                  sheet = 3)

# --- 2. 计算单位方向向量 ---
raw_dirs <- raw %>%
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

# --- 3. 批量计算 R、E、I ---
results <- raw_dirs %>%
  group_by(ID) %>%
  summarise(
    n_scars = n(),
    R       = compute_R(ux, uy, uz),
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
  dplyr::select(-EI) %>%
  arrange(ID)

cat("R/E/I 计算完成，标本数：", nrow(results), "\n")
print(results %>%
        dplyr::select(ID, n_scars, R, E, I) %>%
        mutate(across(c(R, E, I), \(x) round(x, 4))),
      n = Inf)

# --- 4. 连接 Typology，过滤标本（与 5b 保持一致）---
results_typed <- results %>%
  left_join(SPHARM_direction %>% dplyr::select(ID, Typology), by = "ID") %>%
  filter(str_starts(ID, "EXP"),
         !Typology %in% EXCLUDE_TYPES) 
# %>%
  mutate(
    Typology = case_when(
      Typology %in% LEVALLOIS_MERGE ~ "Levallois",
      TRUE ~ Typology
    ),
    Typology = droplevels(as.factor(Typology))
  ) %>%
  filter(complete.cases(R, E, I))

cat(sprintf("\n可用标本: %d，类别数: %d\n",
            nrow(results_typed), nlevels(results_typed$Typology)))
print(table(results_typed$Typology))

# --- 5. 构建特征矩阵 ---
y_rei     <- results_typed$Typology
X_R       <- results_typed %>% dplyr::select(R)     %>% as.matrix()
X_EI      <- results_typed %>% dplyr::select(E, I)  %>% as.matrix()
X_REI     <- results_typed %>% dplyr::select(R, E, I) %>% scale() %>% as.matrix()

# --- 6. LOO-CV 封装函数 ---
run_lda_cv <- function(X, y) {
  fit_cv <- lda(X, grouping = y, CV = TRUE)
  acc    <- mean(fit_cv$class == y)
  conf   <- table(Predicted = fit_cv$class, Actual = y)
  list(acc = acc, conf = conf, pred = fit_cv$class)
}

cv_R   <- run_lda_cv(X_R,   y_rei)
cv_EI  <- run_lda_cv(X_EI,  y_rei)
cv_REI <- run_lda_cv(X_REI, y_rei)

cat("\n========== R / E+I / R+E+I LOO-CV 准确率 ==========\n")
cat(sprintf("随机基准:  %.1f%%\n", 100 / nlevels(y_rei)))
cat(sprintf("R alone:   %.1f%%\n", cv_R$acc   * 100))
cat(sprintf("E+I alone: %.1f%%\n", cv_EI$acc  * 100))
cat(sprintf("R+E+I:     %.1f%%\n", cv_REI$acc * 100))

cat("\n--- 混淆矩阵（列归一化 %）---\n")
cat("R alone:\n");   print(round(prop.table(cv_R$conf,   margin = 2) * 100, 1))
cat("E+I alone:\n"); print(round(prop.table(cv_EI$conf,  margin = 2) * 100, 1))
cat("R+E+I:\n");     print(round(prop.table(cv_REI$conf, margin = 2) * 100, 1))

# --- 7. LDA 可视化（R+E+I）---
lda_fit_REI <- lda(X_REI, grouping = y_rei)
prop_var_REI <- lda_fit_REI$svd^2 / sum(lda_fit_REI$svd^2)

lda_scores_REI <- predict(lda_fit_REI)$x %>%
  as.data.frame() %>%
  mutate(Typology = y_rei, ID = results_typed$ID)

hull_REI <- lda_scores_REI %>%
  group_by(Typology) %>%
  slice(chull(LD1, LD2)) %>%
  ungroup()

p_rei <- ggplot(lda_scores_REI, aes(x = LD1, y = LD2, color = Typology)) +
  geom_polygon(data = hull_REI,
               aes(fill = Typology, group = Typology),
               alpha = 0.15, color = NA) +
  geom_polygon(data = hull_REI,
               aes(color = Typology, group = Typology),
               fill = NA, linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.85) +
  theme_minimal(base_size = 11) +
  labs(
    title    = "LDA: EXP Specimens by Typology",
    subtitle = sprintf("R + E + I, LOO-CV = %.1f%%", cv_REI$acc * 100),
    x = sprintf("LD1 (%.1f%%)", prop_var_REI[1] * 100),
    y = sprintf("LD2 (%.1f%%)", prop_var_REI[2] * 100)
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p_rei)
ggsave(here("analysis/output/figures/LDA_REI_by_Typology.png"),
       plot = p_rei, width = 7, height = 5.5, dpi = 300)

# --- 8. 全方法汇总对比 ---
cat("\n========== 所有方法汇总对比 ==========\n")
cat(sprintf("随机基准:          %.1f%%\n", 100 / nlevels(y_rei)))
cat(sprintf("SPHARM 形态谱:     %.1f%%\n", acc_morph  * 100))
cat(sprintf("SPHARM 方向谱:     %.1f%%\n", acc_dir    * 100))
cat(sprintf("SPHARM cbind拼接:  %.1f%%\n", accuracy   * 100))
cat(sprintf("SPHARM CoIA:       %.1f%%\n", acc_coia   * 100))
cat(sprintf("R alone:           %.1f%%\n", cv_R$acc   * 100))
cat(sprintf("E+I alone:         %.1f%%\n", cv_EI$acc  * 100))
cat(sprintf("R+E+I:             %.1f%%\n", cv_REI$acc * 100))


















# ==============================================================================
# 6. 保存筛选后数据供下游脚本使用
# ==============================================================================

saveRDS(SPHARM_direction_filter,
        here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
saveRDS(SPHARM_morphology_filter,
        here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

cat("已保存：SPHARM_direction_filter.rds\n")
cat("已保存：SPHARM_morphology_filter.rds\n")