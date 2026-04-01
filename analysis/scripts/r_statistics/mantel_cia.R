# ==============================================================================
# mantel_cia.R
# 形态 × 刮痕方向的关联分析：Mantel 检验 + 协惯量分析（CIA）+ PCA
#
# 分析流程：
#   Part A — Mantel 检验
#     1. 对齐标本、构建功率谱矩阵
#     2. 计算余弦距离矩阵
#     3. 全局 Mantel 检验（形态 × 方向）
#     4. 逐变量跨域 Mantel 检验（FDR 校正）
#     5. linkET 网络图可视化
#
#   Part B — 协惯量分析（CIA）
#     1. 零值替换 + ILR 变换
#     2. PCA（dudi.pca）
#     3. Co-inertia 分析 + 置换检验（RV 系数）
#     4. 与 Mantel 检验结果对比
#     5. 可视化（置换分布图 + CIA 双标图 + 碎石图）
#     6. PCA 坐标与载荷图
#
# 前置条件：
#   需先运行 spharm_analysis.R 生成 .rds 文件
#
# 输入：
#   - analysis/data/derived_data/SPHARM_direction_filter.rds
#   - analysis/data/derived_data/SPHARM_morphology_filter.rds
# 输出：
#   - analysis/output/figures/Mantel_Network.png
#   - analysis/output/figures/CIA_ILR_arch.png
#   - analysis/output/figures/CIA_PCA_plots.png
#   - analysis/data/derived_data/CIA_RV_comparison.csv
#   - analysis/data/derived_data/CIA_scores_arch.csv
# ==============================================================================

library(here)
library(tidyverse)
library(vegan)
library(linkET)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(compositions)
library(ade4)


# ==============================================================================
# 读取上游数据（由 spharm_analysis.R 生成）
# ==============================================================================

SPHARM_direction_filter  <- readRDS(here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
SPHARM_morphology_filter <- readRDS(here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))


# ==============================================================================
# Part A：Mantel 检验
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 1：统一 ID 列名并严格对齐标本
# ------------------------------------------------------------------------------

df_morph <- SPHARM_morphology_filter %>% rename(ID = ID)
df_scar  <- SPHARM_direction_filter  %>% rename(ID = ID)

common_ids <- intersect(df_morph$ID, df_scar$ID)

df_morph <- df_morph %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar  <- df_scar  %>% filter(ID %in% common_ids) %>% arrange(ID)

cat("==== 数据对齐检查 ====\n")
cat("共有标本数量：", length(common_ids), "\n")
cat("ID是否完全匹配：", all(df_morph$ID == df_scar$ID), "\n\n")


# ------------------------------------------------------------------------------
# Step 2：构建功率谱矩阵
# ------------------------------------------------------------------------------

morph_power <- df_morph %>%
  select(power_l1:power_l5) %>%
  rename_with(~ paste0("M", 1:5)) %>%
  as.data.frame()

scar_power <- df_scar %>%
  select(power_l1:power_l5) %>%
  rename_with(~ paste0("S", 1:5)) %>%
  as.data.frame()

rownames(morph_power) <- df_morph$ID
rownames(scar_power)  <- df_scar$ID

# 过滤标准差为零的列
morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

cat("==== 矩阵维度检查 ====\n")
cat("形态功率谱：", nrow(morph_power_clean), "标本 ×",
    ncol(morph_power_clean), "阶\n")
cat("片疤功率谱：", nrow(scar_power_clean),  "标本 ×",
    ncol(scar_power_clean),  "阶\n\n")


# ------------------------------------------------------------------------------
# Step 3：构造余弦距离矩阵
# ------------------------------------------------------------------------------

cosine_dist <- function(X) {
  X   <- as.matrix(X)
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  as.dist(1 - sim)
}

D_morph <- cosine_dist(morph_power_clean)
D_scar  <- cosine_dist(scar_power_clean)


# ------------------------------------------------------------------------------
# Step 4：全局 Mantel 检验
# ------------------------------------------------------------------------------

cat("==== 全局Mantel检验（形态 × 片疤方向）====\n")
mantel_global <- mantel(D_morph, D_scar,
                        method = "spearman", permutations = 999)
print(mantel_global)


# ------------------------------------------------------------------------------
# Step 5：逐变量跨域 Mantel 检验（FDR 校正）
# ------------------------------------------------------------------------------

run_cross_mantel <- function(X_single, D_target, from_label, var_label) {
  x <- X_single
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  d_x <- dist(scale(x))
  res <- mantel(d_x, D_target, method = "spearman", permutations = 999)
  tibble(
    from = from_label,
    var  = var_label,
    r    = res$statistic,
    p    = res$signif
  )
}

mantel_M_to_scar <- map_dfr(
  colnames(morph_power_clean),
  ~ run_cross_mantel(morph_power_clean[[.x]], D_scar,
                     from_label = "Scar Direction", var_label = .x)
)

mantel_S_to_morph <- map_dfr(
  colnames(scar_power_clean),
  ~ run_cross_mantel(scar_power_clean[[.x]], D_morph,
                     from_label = "Morphology", var_label = .x)
)

mantel_cross <- bind_rows(mantel_M_to_scar, mantel_S_to_morph) %>%
  mutate(
    p_fdr        = p.adjust(p, method = "fdr"),
    significance = ifelse(p_fdr < 0.05, "P≤0.05", "P>0.05")
  )

cat("\n==== 逐变量跨域Mantel检验结果（FDR校正后）====\n")
mantel_cross %>%
  mutate(across(c(r, p, p_fdr), \(x) round(x, 4))) %>%
  print(n = Inf)

n_before <- sum(mantel_cross$p     < 0.05)
n_after  <- sum(mantel_cross$p_fdr < 0.05)
cat(sprintf("\nFDR校正前显著：%d个，校正后显著：%d个\n", n_before, n_after))


# ------------------------------------------------------------------------------
# Step 6：linkET 网络图可视化
# ------------------------------------------------------------------------------

spec_df_full <- bind_cols(morph_power_clean, scar_power_clean) %>%
  as.data.frame()
rownames(spec_df_full) <- df_morph$ID

p_mantel <- qcorrplot(
  correlate(spec_df_full, method = "spearman"),
  type = "upper", diag = FALSE
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_couple(
    aes(colour = significance, size = abs(r)),
    data         = mantel_cross,
    curvature    = 0.15,
    label.params = list(color = "transparent")
  ) +
  scale_fill_viridis_c(option = "D", limits = c(-1, 1),
                       name = "Spearman's rho") +
  scale_color_manual(
    values = c("P≤0.05" = "#E6A5A5", "P>0.05" = "#BABABA"),
    name   = "Mantel test\n(FDR corrected)"
  ) +
  scale_size_continuous(range = c(0.5, 2.5), name = "Mantel's |r|") +
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    axis.text       = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    plot.margin     = margin(t = 20, r = 20, b = 20, l = 20),
    axis.title      = element_blank()
  )

print(p_mantel)

ggsave(
  here("analysis/output/figures/Mantel_Network.png"),
  plot = p_mantel, width = 10, height = 8, dpi = 300, bg = "white"
)


# ==============================================================================
# Part B：协惯量分析（CIA）— 仅考古标本
# ==============================================================================

morph_arch <- morph_power_clean[
  !str_starts(rownames(morph_power_clean), "IM_"), ]
scar_arch  <- scar_power_clean[
  !str_starts(rownames(scar_power_clean),  "IM_"), ]

stopifnot(all(rownames(morph_arch) == rownames(scar_arch)))
cat("考古标本数量：", nrow(morph_arch), "\n")


# ------------------------------------------------------------------------------
# B-1：零值替换 + ILR 变换
# ------------------------------------------------------------------------------

replace_zeros <- function(X, delta = NULL) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i    <- X[i, ]
    zero_idx <- row_i == 0
    if (!any(zero_idx)) next
    nonzero_min <- min(row_i[!zero_idx])
    d           <- ifelse(is.null(delta), nonzero_min * 0.65, delta)
    n_zero      <- sum(zero_idx)
    row_i[zero_idx]  <- d
    row_i[!zero_idx] <- row_i[!zero_idx] * (1 - n_zero * d)
    X[i, ] <- row_i
  }
  X
}

morph_arch_nozero <- replace_zeros(as.matrix(morph_arch))
scar_arch_nozero  <- replace_zeros(as.matrix(scar_arch))

morph_arch_ilr <- as.data.frame(ilr(morph_arch_nozero))
scar_arch_ilr  <- as.data.frame(ilr(scar_arch_nozero))

rownames(morph_arch_ilr) <- rownames(morph_arch)
rownames(scar_arch_ilr)  <- rownames(scar_arch)

colnames(morph_arch_ilr) <- paste0("M_ilr", seq_len(ncol(morph_arch_ilr)))
colnames(scar_arch_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_arch_ilr)))


# ------------------------------------------------------------------------------
# B-2：PCA（dudi.pca）
# ------------------------------------------------------------------------------

dudi_morph_arch <- dudi.pca(morph_arch_ilr, center = TRUE, scale = TRUE,
                            scannf = FALSE, nf = ncol(morph_arch_ilr))
dudi_scar_arch  <- dudi.pca(scar_arch_ilr,  center = TRUE, scale = TRUE,
                            scannf = FALSE, nf = ncol(scar_arch_ilr))


# ------------------------------------------------------------------------------
# B-3：Co-inertia 分析 + 置换检验
# ------------------------------------------------------------------------------

coin_arch <- coinertia(dudi_morph_arch, dudi_scar_arch,
                       scannf = FALSE, nf = 2)

cat("RV系数（仅考古标本）：", round(coin_arch$RV, 4), "\n")

set.seed(42)
rv_test_arch <- randtest(coin_arch, nrepet = 9999)
cat("\n===== RV置换检验（仅考古标本）=====\n")
print(rv_test_arch)


# ------------------------------------------------------------------------------
# B-4：与 Mantel 检验结果对比
# ------------------------------------------------------------------------------

cat("\n===== 两种方法的全局相关性对比（均仅考古标本）=====\n")
cat(sprintf("Mantel r（余弦距离，Spearman）：%.4f，p = %.3f\n",
            mantel_global$statistic, mantel_global$signif))
cat(sprintf("RV系数（ILR空间，欧氏距离）  ：%.4f，p = %.3f\n",
            coin_arch$RV, rv_test_arch$pvalue))
cat("两种方法结论一致说明结果稳健，不依赖距离度量的选择\n")


# ------------------------------------------------------------------------------
# B-5：CIA 可视化（置换分布图 + 双标图 + 碎石图）
# ------------------------------------------------------------------------------

# --- RV 置换检验分布图 ---
hist_data <- rv_test_arch$plot$hist
rv_null   <- data.frame(mid = hist_data$mids, count = hist_data$counts)

p_rv <- ggplot(rv_null, aes(x = mid, y = count)) +
  geom_col(width = diff(hist_data$breaks)[1] * 0.9,
           fill = "#A1C2E6", color = "white", alpha = 0.8) +
  geom_vline(xintercept = coin_arch$RV,
             color = "#D4619A", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
           x     = coin_arch$RV + diff(range(hist_data$mids)) * 0.03,
           y     = Inf, vjust = 1.5,
           label = sprintf("Observed RV = %.4f\np = %.3f",
                           coin_arch$RV, rv_test_arch$pvalue),
           color = "#D4619A", size = 3.5, hjust = 0) +
  theme_bw(base_size = 10) +
  labs(title = "RV Coefficient Permutation Test (9999 permutations)",
       x = "RV coefficient (null distribution)", y = "Count") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

# --- CIA 双标图 ---
scores_morph <- as.data.frame(coin_arch$lX) %>% rownames_to_column("ID")
scores_scar  <- as.data.frame(coin_arch$lY) %>% rownames_to_column("ID")

scores_combined <- left_join(
  scores_morph %>% select(ID, Axis1_M = AxcX1, Axis2_M = AxcX2),
  scores_scar  %>% select(ID, Axis1_S = AxcY1, Axis2_S = AxcY2),
  by = "ID"
) %>%
  mutate(
    layer = case_when(
      str_detect(ID, "L3") ~ "Layer 3",
      str_detect(ID, "L4") ~ "Layer 4",
      str_detect(ID, "L2") ~ "Layer 2",
      TRUE                  ~ "Other"
    )
  )

cia_inertia <- coin_arch$eig / sum(coin_arch$eig) * 100

p_biplot <- ggplot() +
  geom_segment(
    data = scores_combined,
    aes(x = Axis1_M, y = Axis2_M, xend = Axis1_S, yend = Axis2_S),
    color = "grey70", linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(data = scores_combined,
             aes(x = Axis1_M, y = Axis2_M, color = layer),
             size = 3, alpha = 0.9, shape = 16) +
  geom_point(data = scores_combined,
             aes(x = Axis1_S, y = Axis2_S, color = layer),
             size = 3, alpha = 0.9, shape = 17) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    values = c("Layer 2" = "#E6B89C", "Layer 3" = "#A1C2E6",
               "Layer 4" = "#FFBAE0"),
    name = "Layer"
  ) +
  theme_bw(base_size = 10) +
  labs(
    title = "Co-inertia Analysis Biplot",
    x = sprintf("CIA Axis 1 (%.1f%%)", cia_inertia[1]),
    y = sprintf("CIA Axis 2 (%.1f%%)", cia_inertia[2])
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey50")
  )

# --- CIA 碎石图 ---
eig_df <- data.frame(axis = seq_along(coin_arch$eig), inertia = cia_inertia)

p_scree_cia <- ggplot(eig_df, aes(x = axis, y = inertia)) +
  geom_col(fill = "#B5D5B5", alpha = 0.9, width = 0.6) +
  geom_line(color = "#4A9A4A", linewidth = 0.8) +
  geom_point(color = "#4A9A4A", size = 2.5) +
  scale_x_continuous(breaks = eig_df$axis) +
  theme_bw(base_size = 10) +
  labs(title = "CIA Scree Plot",
       x = "CIA Axis", y = "% Inertia explained") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

# --- 拼合 CIA 三图并保存 ---
p_final_cia <- (p_rv | p_biplot) / p_scree_cia +
  plot_annotation(
    title = "Co-inertia Analysis (ILR-transformed, archaeological specimens only)",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  )

ggsave(
  here("analysis/output/figures/CIA_ILR_arch.png"),
  plot = p_final_cia, width = 12, height = 10, dpi = 300, bg = "white"
)
cat("图已保存：CIA_ILR_arch.png\n")


# ------------------------------------------------------------------------------
# B-6：PCA 坐标与载荷图
# ------------------------------------------------------------------------------

layer_colors <- c("Layer 2" = "#E6B89C",
                  "Layer 3" = "#A1C2E6",
                  "Layer 4" = "#FFBAE0")

# 提取 PCA 坐标与载荷
morph_scores   <- as.data.frame(dudi_morph_arch$li) %>%
  rownames_to_column("ID") %>%
  mutate(layer = case_when(str_detect(ID, "L3") ~ "Layer 3",
                           str_detect(ID, "L4") ~ "Layer 4",
                           str_detect(ID, "L2") ~ "Layer 2",
                           TRUE ~ "Other"))
morph_loadings <- as.data.frame(dudi_morph_arch$co) %>%
  rownames_to_column("variable") %>% mutate(label = variable)
morph_var      <- round(dudi_morph_arch$eig / sum(dudi_morph_arch$eig) * 100, 1)

scar_scores    <- as.data.frame(dudi_scar_arch$li) %>%
  rownames_to_column("ID") %>%
  mutate(layer = case_when(str_detect(ID, "L3") ~ "Layer 3",
                           str_detect(ID, "L4") ~ "Layer 4",
                           str_detect(ID, "L2") ~ "Layer 2",
                           TRUE ~ "Other"))
scar_loadings  <- as.data.frame(dudi_scar_arch$co) %>%
  rownames_to_column("variable") %>% mutate(label = variable)
scar_var       <- round(dudi_scar_arch$eig / sum(dudi_scar_arch$eig) * 100, 1)

cat("===== 形态PCA方差解释 =====\n")
for (i in seq_along(morph_var)) cat(sprintf("  PC%d: %.1f%%\n", i, morph_var[i]))
cat("\n===== 方向PCA方差解释 =====\n")
for (i in seq_along(scar_var))  cat(sprintf("  PC%d: %.1f%%\n", i, scar_var[i]))

# 绘图函数
plot_pca_scores <- function(scores_df, var_pct, title_str,
                            pc_x = 1, pc_y = 2, outlier_n = 5) {
  x_col <- paste0("Axis", pc_x)
  y_col <- paste0("Axis", pc_y)
  scores_df <- scores_df %>%
    mutate(dist       = sqrt(.data[[x_col]]^2 + .data[[y_col]]^2),
           is_outlier = rank(-dist) <= outlier_n)
  p <- ggplot(scores_df, aes(x = .data[[x_col]], y = .data[[y_col]],
                             color = layer)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_point(size = 3, alpha = 0.85) +
    geom_text_repel(data = scores_df %>% filter(is_outlier),
                    aes(label = ID), size = 2.8, max.overlaps = 20,
                    show.legend = FALSE, color = "grey30") +
    scale_color_manual(values = layer_colors, name = "Layer") +
    theme_bw(base_size = 10) +
    labs(title = title_str,
         x = sprintf("PC%d (%.1f%%)", pc_x, var_pct[pc_x]),
         y = sprintf("PC%d (%.1f%%)", pc_y, var_pct[pc_y])) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "right")
  p
}

plot_pca_loadings <- function(loadings_df, var_pct, title_str,
                              pc_x = 1, pc_y = 2) {
  x_col <- paste0("Comp", pc_x)
  y_col <- paste0("Comp", pc_y)
  ggplot(loadings_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_segment(aes(xend = 0, yend = 0),
                 arrow = arrow(length = unit(0.2, "cm"), ends = "first"),
                 color = "#D4619A", linewidth = 0.8, alpha = 0.8) +
    geom_point(color = "#D4619A", size = 2.5) +
    geom_text_repel(aes(label = label), size = 3,
                    color = "grey30", max.overlaps = 20) +
    theme_bw(base_size = 10) +
    labs(title = title_str,
         x = sprintf("PC%d (%.1f%%)", pc_x, var_pct[pc_x]),
         y = sprintf("PC%d (%.1f%%)", pc_y, var_pct[pc_y])) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
}

plot_scree <- function(var_pct, title_str) {
  df <- data.frame(PC = seq_along(var_pct), pct = var_pct,
                   cumsum = cumsum(var_pct))
  ggplot(df, aes(x = PC)) +
    geom_col(aes(y = pct), fill = "#A1C2E6", alpha = 0.8, width = 0.6) +
    geom_line(aes(y = cumsum), color = "#D4619A",
              linewidth = 0.9, linetype = "dashed") +
    geom_point(aes(y = cumsum), color = "#D4619A", size = 2.5) +
    geom_text(aes(y = pct, label = sprintf("%.1f%%", pct)),
              vjust = -0.4, size = 2.8, color = "grey30") +
    scale_x_continuous(breaks = df$PC) +
    scale_y_continuous(limits = c(0, 110),
                       sec.axis = sec_axis(~ ., name = "Cumulative %")) +
    theme_bw(base_size = 10) +
    labs(title = title_str, x = "Principal Component",
         y = "% Variance explained") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
}

# 生成并拼合 PCA 图
p_morph_panel <- (
  plot_scree(morph_var, "Morphology PCA — Scree Plot") |
    plot_pca_scores(morph_scores, morph_var,
                    "Morphology PCA — Sample Scores (ILR space)") |
    plot_pca_loadings(morph_loadings, morph_var,
                      "Morphology PCA — Variable Loadings")
) +
  plot_annotation(
    title = "Morphology PCA (ILR-transformed power spectra, l = 1–5)",
    theme = theme(plot.title = element_text(face = "bold", size = 12,
                                            hjust = 0.5))
  )

p_scar_panel <- (
  plot_scree(scar_var, "Scar Direction PCA — Scree Plot") |
    plot_pca_scores(scar_scores, scar_var,
                    "Scar Direction PCA — Sample Scores (ILR space)") |
    plot_pca_loadings(scar_loadings, scar_var,
                      "Scar Direction PCA — Variable Loadings")
) +
  plot_annotation(
    title = "Scar Direction PCA (ILR-transformed power spectra, l = 1–4)",
    theme = theme(plot.title = element_text(face = "bold", size = 12,
                                            hjust = 0.5))
  )

p_final_pca <- p_morph_panel / p_scar_panel

ggsave(
  here("analysis/output/figures/CIA_PCA_plots.png"),
  plot = p_final_pca, width = 15, height = 10, dpi = 300, bg = "white"
)
cat("图已保存：CIA_PCA_plots.png\n")


# ==============================================================================
# 保存数值结果
# ==============================================================================

tibble(
  method  = c("Mantel (cosine, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic,      coin_arch$RV),
  p_value = c(mantel_global$signif,         rv_test_arch$pvalue)
) %>%
  write_csv(here("analysis/data/derived_data/CIA_RV_comparison.csv"))

scores_combined %>%
  write_csv(here("analysis/data/derived_data/CIA_scores_arch.csv"))

cat("结果已保存：CIA_RV_comparison.csv / CIA_scores_arch.csv\n")