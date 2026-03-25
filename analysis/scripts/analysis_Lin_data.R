library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)
library(htmlwidgets)
library(vegan)
# ============================================================
# 读取数据
# ============================================================
DATA_DIR <- "H:/SPHARM_analysis/analysis/data/derived_data"

df <- read_csv(file.path(DATA_DIR, "SPHARM_direction_lin2024.csv"),
               show_col_types = FALSE)

cat("标本数量：", nrow(df), "\n")
cat("类型分布：\n")
print(df %>% count(Typology))

power_cols <- paste0("power_l", 1:20)
# ============================================================
# 图1：各类型谱熵比较（箱线图）
# ============================================================
p1 <- ggplot(df, aes(x = reorder(Typology, spectral_entropy),
                     y = spectral_entropy,
                     fill = Typology)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  coord_flip() +
  labs(
    title  = "Spectral entropy by typology",
    x      = NULL,
    y      = "Spectral entropy (H)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# ============================================================
# 图2：各类型平均功率谱曲线
# ============================================================
power_cols <- paste0("power_l", 1:20)

df_power <- df %>%
  select(ID, Typology, all_of(power_cols)) %>%
  pivot_longer(
    cols      = all_of(power_cols),
    names_to  = "degree",
    values_to = "power"
  ) %>%
  mutate(degree = as.integer(str_extract(degree, "[0-9]+")))

df_power_mean <- df_power %>%
  group_by(Typology, degree) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    se_power   = sd(power, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p2 <- ggplot(df_power_mean,
             aes(x = degree, y = mean_power,
                 color = Typology, fill = Typology)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(
    aes(ymin = mean_power - se_power,
        ymax = mean_power + se_power),
    alpha = 0.15, color = NA
  ) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(
    title = "Mean power spectrum by typology",
    x     = "Degree (l)",
    y     = "Normalised power",
    color = "Typology",
    fill  = "Typology"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ============================================================
# 图3：功率谱热图（每个标本一行）
# ============================================================
p3 <- ggplot(df_power,
             aes(x = degree, y = ID, fill = power)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "Power") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  facet_grid(Typology ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Power spectrum heatmap",
    x     = "Degree (l)",
    y     = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 7),
    strip.text.y    = element_text(size = 8, angle = 0),
    legend.position = "right"
  )

# ============================================================
# 新增：UMAP 降维分析
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

p4 <- 
  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = transform(df_umap, Typology = NULL), 
             color = "grey90", size = 1) +
  geom_polygon(data = df_hull, aes(fill = Typology), alpha = 0.3) +
  geom_point(aes(color = Typology), size = 1.5) +
  facet_wrap(~Typology, ncol = 4) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.background = element_blank())

p1
p2
p3
p4
# ============================================================
# 导出
# ============================================================
ggsave(file.path(DATA_DIR, "lin2024_spectral_entropy.png"),
       plot = p1, width = 8, height = 6, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_power_spectrum.png"),
       plot = p2, width = 10, height = 6, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_power_heatmap.png"),
       plot = p3, width = 10, height = 12, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_UMAP.png"),
       plot = p4, width = 8, height = 8, dpi = 300)







# ============================================================
# 合并类型
# ============================================================
df <- df %>%
  mutate(Typology_merged = case_when(
    Typology %in% c("Conical bidirectional",
                    "Cylindrical bidirectional")            ~ "bipolar_core",
    Typology %in% c("Conical unidirectional",
                    "Cylindrical unidirectional",
                    "Subconical unidirectional")            ~ "unipolar_core",
    Typology %in% c("Levallois Nubian",
                    "Levallois convergent",
                    "Levallois laminar",
                    "Levallois preferential",
                    "Levallois recurrent")                  ~ "Levallois_core",
    TRUE ~ Typology   # 其余类型保持不变
  ))

# 检查合并结果
cat("合并后类型分布：\n")
print(df %>% count(Typology_merged) %>% arrange(desc(n)))

# ============================================================
# 图1：各类型谱熵比较（箱线图）
# ============================================================
p1 <- ggplot(df, aes(x = reorder(Typology_merged, spectral_entropy),
                     y = spectral_entropy,
                     fill = Typology_merged)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  coord_flip() +
  labs(
    title  = "Spectral entropy by typology",
    x      = NULL,
    y      = "Spectral entropy (H)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# ============================================================
# 图2：各类型平均功率谱曲线
# ============================================================
power_cols <- paste0("power_l", 1:20)

df_power <- df %>%
  select(ID, Typology_merged, all_of(power_cols)) %>%
  pivot_longer(
    cols      = all_of(power_cols),
    names_to  = "degree",
    values_to = "power"
  ) %>%
  mutate(degree = as.integer(str_extract(degree, "[0-9]+")))

df_power_mean <- df_power %>%
  group_by(Typology_merged, degree) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    se_power   = sd(power, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p2 <- ggplot(df_power_mean,
             aes(x = degree, y = mean_power,
                 color = Typology_merged, fill = Typology_merged)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(
    aes(ymin = mean_power - se_power,
        ymax = mean_power + se_power),
    alpha = 0.15, color = NA
  ) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(
    title = "Mean power spectrum by typology",
    x     = "Degree (l)",
    y     = "Normalised power",
    color = "Typology",
    fill  = "Typology"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ============================================================
# 图3：功率谱热图（每个标本一行）
# ============================================================
p3 <- ggplot(df_power,
             aes(x = degree, y = ID, fill = power)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "Power") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  facet_grid(Typology_merged ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Power spectrum heatmap",
    x     = "Degree (l)",
    y     = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 7),
    strip.text.y    = element_text(size = 8, angle = 0),
    legend.position = "right"
  )

# ============================================================
# 新增：UMAP 降维分析
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
  group_by(Typology_merged) %>%
  slice(chull(UMAP1, UMAP2)) %>%
  ungroup()

p4 <- 
  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = transform(df_umap, Typology_merged = NULL), 
             color = "grey90", size = 1) +
  geom_polygon(data = df_hull, aes(fill = Typology_merged), alpha = 0.3) +
  geom_point(aes(color = Typology_merged), size = 1.5) +
  facet_wrap(~Typology_merged, ncol = 3) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.background = element_blank())

p1
p2
p3
p4
# ============================================================
# 导出
# ============================================================
ggsave(file.path(DATA_DIR, "lin2024_spectral_entropy_merged.png"),
       plot = p1, width = 8, height = 6, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_power_spectrum_merged.png"),
       plot = p2, width = 10, height = 6, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_power_heatmap_merged.png"),
       plot = p3, width = 10, height = 12, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_UMAP_merged.png"),
       plot = p4, width = 8, height = 6, dpi = 300)




# 剔除 Biface 并更新数据
df_stat <- df %>% 
  filter(Typology_merged != "Biface") %>%
  mutate(Typology_merged = factor(Typology_merged)) # 重置 factor levels

# 提取特征矩阵
power_matrix_stat <- df_stat %>% select(all_of(power_cols))

# 2. 运行 PERMANOVA 检验
# formula: 矩阵 ~ 分组变量
# permutations: 建议 999 次或更高
permanova_res <- adonis2(power_matrix ~ Typology_merged, 
                         data = df, 
                         permutations = 999, 
                         method = "euclidean")

# 查看全局检验结果
print(permanova_res)


# 计算组间离散度
disper_res <- betadisper(d = dist(power_matrix), group = df$Typology_merged, type = c("centroid"))
permutest(disper_res)



library(pairwiseAdonis)
df_clean <- df %>% 
  filter(Typology_merged != "Biface") %>% 
  mutate(Typology_merged = factor(Typology_merged))

# 提取对应的功率谱矩阵
power_matrix_clean <- df_clean %>% select(all_of(power_cols))

pairwise_perm <- pairwise.adonis(
  x = power_matrix_clean, 
  factors = df_clean$Typology_merged, 
  sim.method = "euclidean", 
  p.adjust.m = "bonferroni" # 严格的校正方法
)

print("--- Pairwise PERMANOVA Results ---")
print(pairwise_perm)



