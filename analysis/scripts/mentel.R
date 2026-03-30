library(tidyverse)
library(vegan)
library(linkET)
library(ggplot2)
library(patchwork)

# ============================================================
# Step 1：统一ID列名并严格对齐标本
# ============================================================
df_morph <- SPHARM_morphology_filter %>% rename(ID = ID)
df_scar  <- SPHARM_direction_filter  %>% rename(ID = ID)

# 找出两个数据集中共有的标本 ID（避免因为个别标本缺失导致矩阵错位）
common_ids <- intersect(df_morph$ID, df_scar$ID)

df_morph <- df_morph %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar  <- df_scar  %>% filter(ID %in% common_ids) %>% arrange(ID)

cat("==== 数据对齐检查 ====\n")
cat("共有标本数量：", length(common_ids), "\n")
cat("ID是否完全匹配：", all(df_morph$ID == df_scar$ID), "\n\n")

# ============================================================
# Step 2：构建功率谱矩阵 (固定为 1-5 阶)
# ============================================================
# 提取形态矩阵 (M1-M5)
morph_power <- df_morph %>%
  select(power_l1:power_l5) %>%
  rename_with(~ paste0("M", 1:5)) %>%
  as.data.frame()

# 提取片疤方向矩阵 (S1-S5)
scar_power <- df_scar %>%
  select(power_l1:power_l5) %>%
  rename_with(~ paste0("S", 1:5)) %>%
  as.data.frame()

rownames(morph_power) <- df_morph$ID
rownames(scar_power)  <- df_scar$ID

# 过滤标准差为零的列（防止缩放时分母为 0 报错）
morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

# ============================================================
# Step 3：构造距离矩阵 (余弦距离)
# ============================================================
normalize_spectra <- function(X) {
  total <- rowSums(X)
  X / ifelse(total > 0, total, 1)
}

cosine_dist <- function(X) {
  sim <- X %*% t(X) / (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  as.dist(1 - sim)
}

# 归一化并计算距离矩阵
X_morph_norm <- morph_power_clean %>% as.matrix() %>% normalize_spectra()
X_scar_norm  <- scar_power_clean  %>% as.matrix() %>% normalize_spectra()

D_morph <- cosine_dist(X_morph_norm)
D_scar  <- cosine_dist(X_scar_norm)

# ============================================================
# Step 4：全局 Mantel 检验
# ============================================================
cat("==== 全局Mantel检验（形态 × 片疤方向）====\n")
mantel_global <- mantel(D_morph, D_scar,
                        method       = "spearman",
                        permutations = 9999)
print(mantel_global)

# ============================================================
# Step 5：逐变量 Mantel 检验 (准备可视化数据)
# ============================================================
# 将两组特征合并，用于计算内部的 Pearson 相关性
spec_df_full <- bind_cols(morph_power_clean, scar_power_clean) %>%
  as.data.frame()
rownames(spec_df_full) <- df_morph$ID

# 循环计算每个阶数与全局矩阵的 Mantel 关系
mantel_rows_full <- map_dfr(colnames(spec_df_full), function(var) {
  x <- spec_df_full[[var]]
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  
  # 对单一变量计算欧式距离
  d_x <- dist(scale(x))
  
  res_m <- mantel(d_x, D_morph, method = "spearman", permutations = 9999)
  res_s <- mantel(d_x, D_scar,  method = "spearman", permutations = 9999)
  
  bind_rows(
    tibble(spec = "Morphology",     env = var,
           r = res_m$statistic, p = res_m$signif),
    tibble(spec = "Scar Direction", env = var,
           r = res_s$statistic, p = res_s$signif)
  )
})

# ============================================================
# Step 6：可视化 (linkET 风格网络图)
# ============================================================

mantel_rows_full <- mantel_rows_full %>%
  mutate(significance = ifelse(p < 0.05, "P≤0.05", "P>0.05"))

p_mantel <- qcorrplot(
  correlate(spec_df_full, method = "spearman"), 
  type = "upper",
  diag = FALSE
) +
  geom_tile(color = "white", linewidth = 0.5) +
  
  geom_couple(
    aes(
      colour = ifelse(p < 0.05, "P≤0.05", "P>0.05"),
      size   = abs(r)
    ),
    data      = mantel_rows_full,
    curvature = 0.15,
    label.params = list(color = "transparent") 
  ) +
  
  scale_fill_viridis_c(
    option = "D",
    limits = c(-1, 1),
    name   = "Spearman's rho",
  ) +
  
  scale_color_manual(
    values = c(
      "P≤0.05"     = "#E6A5A5",
      "P>0.05" = "#BABABA"
    ),
    name = "Mantel test"
  ) +
  
  scale_size_continuous(
    range = c(0.5, 2.5), 
    name  = "Mantel's |r|"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    axis.text       = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    plot.margin     = margin(t = 20, r = 20, b = 20, l = 20),
    axis.title = element_blank()
  ) 

print(p_mantel)

ggsave(
  filename = "analysis/data/derived_data/Mantel_Network_HighRes.png",
  plot     = p_mantel,
  width    = 10,  
  height   = 8,    
  dpi      = 300,   
  bg       = "white" 
)
