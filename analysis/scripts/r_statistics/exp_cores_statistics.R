# ==============================================================================
# mantel_cia_exp.R
# 形态 × 刮痕方向关联分析——实验标本（EXP）
#
# 分析框架：
#   L1：整体 Mantel + CoIA（基线）
#   L2-A：分组 Mantel（按 Typology）
#   L2-B：CoIA 箭头长度分组差异
#   L2-C：CoIA 箭头方位圆形统计
#   L2-D：spectral_entropy × Typology（方向 & 形态数据）
#
# 输入：
#   analysis/data/derived_data/SPHARM_direction.csv
#   analysis/data/derived_data/SPHARM_morphology.csv
#   analysis/data/raw_data/SDG_core_metric.xlsx
#
# 输出（figures/）：
#   EXP_L1_Mantel_Network.png
#   EXP_L1_CIA_Biplot.png
#   EXP_L2_Mantel_Grouped_dumbbell.png
#   EXP_L2_Mantel_Grouped_heatmap.png
#   EXP_L2_Arrow_Length_typology.png
#   EXP_L2_Arrow_Direction_rose_typology.png
#   EXP_L2D_SE_Direction_Typology_boxplot.png    # 新增
#   EXP_L2D_SE_Morphology_Typology_boxplot.png   # 新增
#
# 输出（derived_data/）：
#   EXP_L1_results.csv
#   EXP_L2_grouped_mantel.csv
#   EXP_L2_arrow_stats.csv
#   EXP_L2_circular_stats.csv
#   EXP_CIA_scores_full.csv
#   EXP_L2D_SE_desc_stats.csv                    # 新增
#   EXP_L2D_SE_dunn_results.csv                  # 新增
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(vegan)
library(linkET)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(compositions)
library(ade4)
library(circular)
library(FSA)       # dunnTest()          ← L2-D 新增依赖


# ==============================================================================
# ---- 全局辅助函数 ----
# ==============================================================================

cosine_dist <- function(X) {
  X   <- as.matrix(X)
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  as.dist(1 - sim)
}

replace_zeros <- function(X, delta = NULL) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i      <- X[i, ]
    zero_idx   <- row_i == 0
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

extract_subdist <- function(D_full, ids) {
  as.dist(as.matrix(D_full)[ids, ids])
}

safe_filter_groups <- function(meta_df, group_col, min_n = 3) {
  counts <- table(meta_df[[group_col]], useNA = "no")
  valid  <- names(counts[counts >= min_n])
  if (length(valid) < 2) {
    cat(sprintf(
      "  [跳过] %s：有效组数不足（需 ≥ 2 组每组 ≥ %d 件）。当前：%s\n",
      group_col, min_n,
      paste(names(counts), counts, sep = "=", collapse = ", ")
    ))
    return(NULL)
  }
  meta_df %>%
    filter(!is.na(.data[[group_col]]),
           .data[[group_col]] %in% valid)
}

circ_stats_one <- function(angles_rad) {
  circ_obj <- circular(angles_rad, type = "angles",
                       units = "radians", modulo = "2pi")
  mean_rad <- as.numeric(mean.circular(circ_obj)) %% (2 * pi)
  list(
    mean_rad = mean_rad,
    mean_deg = mean_rad * 180 / pi,
    rho      = as.numeric(rho.circular(circ_obj))
  )
}

# Watson 两样本置换检验（circular）
watson_perm_test <- function(x1, x2, B = 9999) {
  a1 <- circular(x1, type = "angles", units = "radians", modulo = "2pi")
  a2 <- circular(x2, type = "angles", units = "radians", modulo = "2pi")
  obs_u2 <- as.numeric(watson.two.test(a1, a2)$statistic)
  x_all  <- c(x1, x2)
  n1     <- length(x1)
  n_all  <- length(x_all)
  perm_u2 <- replicate(B, {
    idx <- sample.int(n_all)
    as.numeric(watson.two.test(
      circular(x_all[idx[1:n1]],            type = "angles", units = "radians", modulo = "2pi"),
      circular(x_all[idx[(n1 + 1):n_all]], type = "angles", units = "radians", modulo = "2pi")
    )$statistic)
  })
  list(statistic = obs_u2,
       p.value   = (sum(perm_u2 >= obs_u2) + 1) / (B + 1))
}


# ==============================================================================
# ---- 数据准备 ----
# ==============================================================================

POWER_COLS <- paste0("power_l", 1:5)   # 修改为实际功率谱列名

SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"),
                              show_col_types = FALSE)
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"),
                              show_col_types = FALSE)
metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

# 将 Typology 并入形态数据
SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

# 通用筛选函数
filter_spharm <- function(df, meta = NULL) {
  result <- df %>%
    select(ID, Typology, SHE, spectral_entropy, all_of(POWER_COLS))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}

SPHARM_direction_filter  <- filter_spharm(SPHARM_direction,  metric_data)
SPHARM_morphology_filter <- filter_spharm(SPHARM_morphology, metric_data)

# 划分子集：EXP（含 IM_）
split_by_group <- function(df) {
  list(
    exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
    sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_"))
  )
}
dir_splits   <- split_by_group(SPHARM_direction_filter)
morph_splits <- split_by_group(SPHARM_morphology_filter)

# ---------- 对齐实验标本（排除 IM_）----------
df_scar_all  <- dir_splits$exp_im
df_morph_all <- morph_splits$exp_im

common_ids <- intersect(df_morph_all$ID, df_scar_all$ID)
df_morph_all <- df_morph_all %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar_all  <- df_scar_all  %>% filter(ID %in% common_ids) %>% arrange(ID)
stopifnot(all(df_morph_all$ID == df_scar_all$ID))

cat("==== 数据对齐（EXP + IM）====\n")
cat("共有标本：", length(common_ids),
    "；ID 完全匹配：", all(df_morph_all$ID == df_scar_all$ID), "\n\n")

# 功率谱矩阵（全集，含 IM_）
morph_power_all <- df_morph_all %>%
  select(all_of(POWER_COLS)) %>%
  rename_with(~ paste0("M", seq_along(.))) %>%
  as.data.frame()
scar_power_all <- df_scar_all %>%
  select(all_of(POWER_COLS)) %>%
  rename_with(~ paste0("S", seq_along(.))) %>%
  as.data.frame()
rownames(morph_power_all) <- df_morph_all$ID
rownames(scar_power_all)  <- df_scar_all$ID

morph_power_clean <- morph_power_all[, sapply(morph_power_all, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power_all[,  sapply(scar_power_all,  sd, na.rm = TRUE) > 0]

D_morph_all <- cosine_dist(morph_power_clean)
D_scar_all  <- cosine_dist(scar_power_clean)

# 纯实验标本（排除 IM_）
exp_ids <- rownames(morph_power_clean)[!str_starts(rownames(morph_power_clean), "IM_")]
cat("纯实验标本数量（不含 IM_）：", length(exp_ids), "\n")

morph_exp <- morph_power_clean[exp_ids, ]
scar_exp  <- scar_power_clean[exp_ids, ]

D_morph_exp <- extract_subdist(D_morph_all, exp_ids)
D_scar_exp  <- extract_subdist(D_scar_all,  exp_ids)

# 元数据（实验标本）
meta_exp <- df_morph_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, Typology) %>%
  left_join(metric_data, by = "ID")

# ---------- 合并所有 Levallois 变体 ----------
# 凡 Typology 中含 "Levallois"（不区分大小写）的均归并为 "Levallois"
meta_exp <- meta_exp %>%
  mutate(Typology = if_else(
    str_detect(Typology, regex("levallois", ignore_case = TRUE)),
    "Levallois",
    Typology
  ))

# 同步更新 scores_combined 中的 Typology（CoIA 完成后会重新 left_join，此处预处理即可）
# 也同步到 df_morph_all / df_scar_all（用于后续筛选）
df_morph_all <- df_morph_all %>%
  mutate(Typology = if_else(
    str_detect(Typology, regex("levallois", ignore_case = TRUE)),
    "Levallois",
    Typology
  ))

cat("\n==== 实验标本元数据（Levallois 已合并）====\n")
cat("Typology：\n"); print(table(meta_exp$Typology, useNA = "ifany"))

meta_typology <- safe_filter_groups(meta_exp, "Typology")

# 调色板（按 Typology 实际水平自动生成，或手动覆盖）
typology_levels <- sort(unique(meta_exp$Typology[!is.na(meta_exp$Typology)]))
set.seed(42)
default_pal <- setNames(
  colorRampPalette(c("#7EB8C9", "#E6B89C", "#FFBAE0",
                     "#A1C2E6", "#6271A1", "#C6DEA4",
                     "#90C49A", "#D4A5A3", "#D6D6D6"))(length(typology_levels)),
  typology_levels
)
# 如需手动指定，取消注释并填写：
# typology_pal <- c("Levallois" = "#7EB8C9", "TypeB" = "#E6B89C", ...)
typology_pal <- default_pal


# ==============================================================================
# ========== 第一层：整体 Mantel + CoIA（基线）==========
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  第一层：整体 Mantel + CoIA — 建立基线（EXP）             ##\n")
cat("################################################################\n")


# ------------------------------------------------------------------------------
# L1-1：全局 Mantel + linkET 网络图
# ------------------------------------------------------------------------------

cat("\n==== L1-1：全局 Mantel ====\n")
mantel_global <- mantel(D_morph_exp, D_scar_exp,
                        method = "spearman", permutations = 9999)
print(mantel_global)

run_cross_mantel <- function(X_single, D_target, from_label, var_label,
                             n_perm = 999) {
  if (sd(X_single, na.rm = TRUE) == 0) return(NULL)
  res <- mantel(dist(scale(X_single)), D_target,
                method = "spearman", permutations = n_perm)
  tibble(from = from_label, var = var_label,
         r = res$statistic, p = res$signif)
}

morph_exp_df <- as.data.frame(morph_exp)
scar_exp_df  <- as.data.frame(scar_exp)

mantel_cross_l1 <- bind_rows(
  map_dfr(colnames(morph_exp_df),
          ~ run_cross_mantel(morph_exp_df[[.x]], D_scar_exp,
                             "Scar Direction", .x)),
  map_dfr(colnames(scar_exp_df),
          ~ run_cross_mantel(scar_exp_df[[.x]], D_morph_exp,
                             "Morphology", .x))
) %>%
  mutate(
    p_holm       = p.adjust(p, method = "holm"),
    significance = ifelse(p_holm < 0.05, "P≤0.05", "P>0.05")
  )

spec_exp_full <- bind_cols(morph_exp_df, scar_exp_df)

p_mantel_net <- qcorrplot(
  correlate(spec_exp_full, method = "spearman"),
  type = "upper", diag = FALSE
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_couple(
    aes(colour = significance, size = abs(r)),
    data         = mantel_cross_l1,
    curvature    = 0.15,
    label.params = list(color = "transparent")
  ) +
  scale_fill_viridis_c(option = "D", limits = c(-1, 1),
                       name = "Spearman's rho") +
  scale_color_manual(
    values = c("P≤0.05" = "#E6A5A5", "P>0.05" = "#BABABA"),
    name   = "Mantel test\n(Holm corrected)"
  ) +
  scale_size_continuous(range = c(0.5, 2.5), name = "Mantel's |r|") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(20, 20, 20, 20))

ggsave(here("analysis/output/figures/EXP_L1_Mantel_Network.png"),
       plot = p_mantel_net, width = 10, height = 8, dpi = 300, bg = "white")
cat("图已保存：EXP_L1_Mantel_Network.png\n")


# ------------------------------------------------------------------------------
# L1-2：CoIA + RV 置换检验
# ------------------------------------------------------------------------------

cat("\n==== L1-2：CoIA ====\n")

morph_exp_ilr <- as.data.frame(ilr(replace_zeros(as.matrix(morph_exp))))
scar_exp_ilr  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_exp))))
rownames(morph_exp_ilr) <- exp_ids
rownames(scar_exp_ilr)  <- exp_ids
colnames(morph_exp_ilr) <- paste0("M_ilr", seq_len(ncol(morph_exp_ilr)))
colnames(scar_exp_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_exp_ilr)))

dudi_morph <- dudi.pca(morph_exp_ilr, center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(morph_exp_ilr))
dudi_scar  <- dudi.pca(scar_exp_ilr,  center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(scar_exp_ilr))

coin_exp   <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
cia_inertia <- coin_exp$eig / sum(coin_exp$eig) * 100

cat("RV 系数：", round(coin_exp$RV, 4), "\n")

set.seed(42)
rv_test <- randtest(coin_exp, nrepet = 9999)
cat("\nRV 置换检验：\n"); print(rv_test)

# CIA 坐标 + 箭头
scores_morph <- as.data.frame(coin_exp$lX) %>% rownames_to_column("ID")
scores_scar  <- as.data.frame(coin_exp$lY) %>% rownames_to_column("ID")

scores_combined <- left_join(
  scores_morph %>% select(ID, Axis1_M = AxcX1, Axis2_M = AxcX2),
  scores_scar  %>% select(ID, Axis1_S = AxcY1, Axis2_S = AxcY2),
  by = "ID"
) %>%
  mutate(
    arrow_length = sqrt((Axis1_M - Axis1_S)^2 + (Axis2_M - Axis2_S)^2),
    arrow_angle  = atan2(Axis2_S - Axis2_M, Axis1_S - Axis1_M)
  ) %>%
  left_join(meta_exp %>% select(ID, Typology), by = "ID")

# CIA 双标图
p_cia_biplot <- ggplot() +
  geom_segment(
    data = scores_combined,
    aes(x = Axis1_M, y = Axis2_M,
        xend = Axis1_S, yend = Axis2_S,
        color = arrow_length),
    linewidth = 0.7, alpha = 0.75,
    arrow = arrow(length = unit(0.10, "cm"), type = "closed")
  ) +
  geom_point(data = scores_combined,
             aes(x = Axis1_M, y = Axis2_M, shape = Typology),
             size = 2.8, fill = "white", color = "grey25",
             stroke = 1.1, alpha = 0.9) +
  geom_point(data = scores_combined,
             aes(x = Axis1_S, y = Axis2_S, shape = Typology),
             size = 2.8, alpha = 0.85, color = "grey55") +
  geom_text_repel(
    data = scores_combined %>% slice_max(arrow_length, n = 5),
    aes(x = (Axis1_M + Axis1_S) / 2,
        y = (Axis2_M + Axis2_S) / 2,
        label = ID),
    size = 2.3, color = "grey30", max.overlaps = 15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.35) +
  scale_color_viridis_c(option = "C", name = "Decoupling\n(arrow length)",
                        direction = -1) +
  scale_shape_manual(
    values = setNames(c(21, 22, 24, 23, 25, 20, 15, 17, 16, 18)[
      seq_along(typology_levels)], typology_levels),
    name = "Typology\n(open = morph / grey = scar)"
  ) +
  theme_bw(base_size = 10) +
  labs(
    title    = sprintf("CoIA Biplot (EXP)  |  RV = %.3f, p = %.3f",
                       coin_exp$RV, rv_test$pvalue),
    subtitle = "Arrow: morphology → scar-direction coordinate; length = decoupling",
    x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
    y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
    legend.position = "right"
  )

ggsave(here("analysis/output/figures/EXP_L1_CIA_Biplot.png"),
       plot = p_cia_biplot, width = 10, height = 8, dpi = 300, bg = "white")
cat("图已保存：EXP_L1_CIA_Biplot.png\n")

l1_results <- tibble(
  method  = c("Mantel (cosine, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_exp$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(exp_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/EXP_L1_results.csv"))
cat("\n第一层结论：\n"); print(l1_results)


# ==============================================================================
# ========== 第二层：联合证据 ==========
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  第二层：联合证据（EXP）                                   ##\n")
cat("################################################################\n")


# ==============================================================================
# ---- L2-A：分组 Mantel ----
# ==============================================================================

cat("\n---------- L2-A：各类型独立 Mantel ----------\n")
cat("策略：对每个 Typology 子集分别运行一次独立 Mantel 检验，\n")
cat("      p 值为该子集内置换检验原始值，不做跨组合并校正。\n\n")

# 单组 Mantel：仅用该类型内部的标本
mantel_one_type <- function(type_val, meta_df, D_morph_full, D_scar_full,
                            n_perm = 9999) {
  ids <- meta_df %>% filter(Typology == type_val) %>% pull(ID)
  cat(sprintf("  → %s (n = %d) ...", type_val, length(ids)))
  if (length(ids) < 5) {
    cat(" 跳过（n < 5）\n")
    return(NULL)
  }
  res <- mantel(extract_subdist(D_morph_full, ids),
                extract_subdist(D_scar_full,  ids),
                method = "spearman", permutations = n_perm)
  cat(sprintf(" r = %.4f, p = %.4f\n", res$statistic, res$signif))
  tibble(
    Typology    = type_val,
    n           = length(ids),
    mantel_r    = res$statistic,
    p_value     = res$signif,
    significance = case_when(
      res$signif < 0.001 ~ "***",
      res$signif < 0.01  ~ "**",
      res$signif < 0.05  ~ "*",
      res$signif < 0.10  ~ ".",
      TRUE               ~ "ns"
    )
  )
}

# 所有有效 Typology 水平
all_types <- meta_exp %>%
  filter(!is.na(Typology)) %>%
  count(Typology) %>%
  arrange(desc(n)) %>%
  pull(Typology)

mantel_by_typology <- map_dfr(all_types,
                              ~ mantel_one_type(.x, meta_exp,
                                                D_morph_exp, D_scar_exp))

# 兼容后续可视化字段（保持列名一致）
mantel_by_typology <- mantel_by_typology %>%
  mutate(group_var = "Typology", group = Typology, p_raw = p_value,
         p_fdr = p_value)   # 独立检验，p_fdr 即原始 p

cat("\n==== 各类型 Mantel 结果汇总 ====\n")
print(mantel_by_typology %>%
        select(Typology, n, mantel_r, p_value, significance) %>%
        mutate(across(c(mantel_r, p_value), ~ round(.x, 4))))

# 合并 + 可视化
l2_mantel <- mantel_by_typology %>%
  mutate(
    group_var_label = "Typology",
    group = fct_reorder(as.factor(group), mantel_r)
  )

if (!is.null(l2_mantel) && nrow(l2_mantel) > 0) {
  
  # 哑铃图
  p_l2_dumbbell <- ggplot(l2_mantel,
                          aes(x = mantel_r, y = group,
                              color = significance)) +
    geom_vline(xintercept = mantel_global$statistic,
               linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "grey70", linewidth = 0.4) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = sprintf("n=%d", n)),
              hjust = -0.35, size = 2.5, color = "grey40") +
    scale_color_manual(
      values = c("***" = "#C0392B", "**" = "#E67E22", "*" = "#F1C40F",
                 "."  = "#27AE60", "ns" = "#95A5A6"),
      name = "Significance\n(per-type raw p)"
    ) +
    theme_bw(base_size = 10) +
    labs(
      title    = "EXP L2-A: Per-type Mantel Test (Morphology × Scar Direction)",
      subtitle = "Each type tested independently; dashed line = global baseline; p = permutation p (no pooling)",
      x = "Mantel r (Spearman)", y = "Typology"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      legend.position = "right"
    )
  
  ggsave(here("analysis/output/figures/EXP_L2_Mantel_Grouped_dumbbell.png"),
         plot = p_l2_dumbbell,
         width = 10, height = max(4, nrow(l2_mantel) * 0.7 + 2),
         dpi = 300, bg = "white")
  cat("图已保存：EXP_L2_Mantel_Grouped_dumbbell.png\n")
  
  # 热图
  heatmap_df <- l2_mantel %>%
    mutate(y_int = as.integer(group))
  
  y_labels <- heatmap_df %>%
    distinct(y_int, group) %>%
    arrange(y_int) %>%
    mutate(group_chr = as.character(group))
  
  p_l2_heatmap <- ggplot(heatmap_df,
                         aes(x = group_var_label, y = y_int,
                             fill = mantel_r)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("r = %.3f\n%s", mantel_r, significance)),
              size = 2.8, color = "grey20") +
    scale_fill_gradient2(
      low = "#3B82C4", mid = "white", high = "#C0392B",
      midpoint = 0, limits = c(-0.3, 0.3),
      oob = scales::squish, name = "Mantel r"
    ) +
    scale_y_continuous(
      breaks = y_labels$y_int,
      labels = y_labels$group_chr,
      expand = expansion(add = 0.5)
    ) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    theme_bw(base_size = 10) +
    labs(
      title    = "EXP L2-A: Per-type Mantel r (Morphology × Scar Direction)",
      subtitle = sprintf("Global baseline: r = %.3f (p = %.3f) | p = per-type permutation, no pooling",
                         mantel_global$statistic, mantel_global$signif),
      x = "Typology", y = "Subgroup"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      axis.text.y   = element_text(size = 8)
    )
  
  ggsave(here("analysis/output/figures/EXP_L2_Mantel_Grouped_heatmap.png"),
         plot   = p_l2_heatmap,
         width  = 6,
         height = max(4, nrow(l2_mantel) * 0.6 + 2),
         dpi = 300, bg = "white")
  cat("图已保存：EXP_L2_Mantel_Grouped_heatmap.png\n")
  
  write_csv(l2_mantel, here("analysis/data/derived_data/EXP_L2_grouped_mantel.csv"))
}


# ==============================================================================
# ---- L2-B：CoIA 箭头长度分组差异 ----
# ==============================================================================

cat("\n---------- L2-B：箭头长度（Typology）----------\n")

run_arrow_length_analysis <- function(group_col, group_label, palette) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 3) %>%
    pull(.data[[group_col]]) %>%
    unique()
  
  if (length(valid_groups) < 2) {
    cat(sprintf("  [跳过] %s 箭头长度检验：有效分组不足\n", group_label))
    return(invisible(NULL))
  }
  
  sub_df <- scores_combined %>%
    filter(.data[[group_col]] %in% valid_groups)
  
  cat(sprintf("\n----- %s × 箭头长度 -----\n", group_label))
  kw <- kruskal.test(reformulate(group_col, "arrow_length"), data = sub_df)
  print(kw)
  
  pw <- pairwise.wilcox.test(sub_df$arrow_length, sub_df[[group_col]],
                             p.adjust.method = "holm", exact = FALSE)
  print(pw)
  
  p <- ggplot(sub_df,
              aes(x = .data[[group_col]],
                  y = arrow_length,
                  fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.5) +
    geom_jitter(aes(color = .data[[group_col]]), width = 0.15,
                size = 2.5, alpha = 0.75, show.legend = FALSE) +
    geom_text_repel(
      data = sub_df %>%
        group_by(.data[[group_col]]) %>%
        slice_max(arrow_length, n = 1),
      aes(label = ID), size = 2.4, color = "grey40",
      max.overlaps = 10, show.legend = FALSE
    ) +
    annotate("text",
             x     = 1.5,
             y     = max(sub_df$arrow_length) * 1.02,
             label = sprintf("Kruskal-Wallis: χ² = %.2f, p = %.3f",
                             kw$statistic, kw$p.value),
             size = 3, color = "grey30", hjust = 0.5) +
    scale_fill_manual(values  = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    theme_bw(base_size = 10) +
    labs(
      title    = sprintf("EXP L2-B: CoIA Arrow Length by %s", group_label),
      subtitle = "Longer arrows = greater morphology–technique decoupling",
      x = NULL, y = "Arrow length (CoIA space)"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      axis.text.x   = element_text(angle = 20, hjust = 1)
    )
  
  fname <- sprintf("analysis/output/figures/EXP_L2_Arrow_Length_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  
  list(sub_df = sub_df, kw = kw, pw = pw)
}

res_len_typology <- run_arrow_length_analysis("Typology", "Typology", typology_pal)

cat("\n==== 箭头长度描述统计（Typology）====\n")
scores_combined %>%
  filter(!is.na(Typology)) %>%
  group_by(Typology) %>%
  summarise(
    n      = n(),
    mean   = round(mean(arrow_length),   4),
    median = round(median(arrow_length), 4),
    sd     = round(sd(arrow_length),     4),
    min    = round(min(arrow_length),    4),
    max    = round(max(arrow_length),    4),
    .groups = "drop"
  ) %>%
  print()


# ==============================================================================
# ---- L2-C：CoIA 箭头方位圆形统计 ----
# ==============================================================================

cat("\n---------- L2-C：箭头方位（Typology）----------\n")

run_circular_analysis <- function(group_col, group_label, palette) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 5) %>%
    pull(.data[[group_col]]) %>%
    unique()
  
  if (length(valid_groups) < 2) {
    cat(sprintf("  [跳过] %s 圆形统计：有效组数不足（需 ≥ 2 组每组 ≥ 5 件）\n",
                group_label))
    return(invisible(NULL))
  }
  
  sub_df <- scores_combined %>%
    filter(.data[[group_col]] %in% valid_groups)
  
  # 描述统计
  cat(sprintf("\n----- %s 圆形描述统计 -----\n", group_label))
  circ_desc <- map_dfr(valid_groups, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(
      group_var       = group_col,
      group           = g,
      n               = length(angles),
      mean_dir_deg    = round(cs$mean_deg, 2),
      concentration_r = round(cs$rho,     4)
    )
  })
  print(circ_desc)
  
  # Rayleigh 检验
  cat(sprintf("\n----- %s Rayleigh 检验 -----\n", group_label))
  rayleigh_res <- map_dfr(valid_groups, function(g) {
    angles   <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    circ_obj <- circular(angles, type = "angles", units = "radians", modulo = "2pi")
    rt       <- rayleigh.test(circ_obj)
    cat(sprintf("  %s: U = %.4f, p = %.4f → %s\n",
                g, rt$statistic, rt$p.value,
                ifelse(rt$p.value < 0.05, "方位集中（非随机）", "方位分散（随机）")))
    tibble(
      group_var  = group_col,
      group      = g,
      rayleigh_U = round(rt$statistic, 4),
      rayleigh_p = round(rt$p.value,   4),
      conclusion = ifelse(rt$p.value < 0.05, "concentrated", "uniform")
    )
  })
  
  # Watson 两样本检验
  watson_res <- NULL
  if (length(valid_groups) >= 2) {
    cat(sprintf("\n----- %s Watson 两样本检验（置换 p 值）-----\n", group_label))
    pairs <- combn(valid_groups, 2, simplify = FALSE)
    watson_res <- map_dfr(pairs, function(pair) {
      x1 <- sub_df %>% filter(.data[[group_col]] == pair[1]) %>% pull(arrow_angle)
      x2 <- sub_df %>% filter(.data[[group_col]] == pair[2]) %>% pull(arrow_angle)
      wt <- watson_perm_test(x1, x2, B = 9999)
      cat(sprintf("  %s vs %s: U² = %.4f, p = %.4f → %s\n",
                  pair[1], pair[2], wt$statistic, wt$p.value,
                  ifelse(wt$p.value < 0.05, "分布不同", "差异不显著")))
      tibble(
        group_var    = group_col,
        group1       = pair[1],
        group2       = pair[2],
        U2_statistic = round(wt$statistic, 4),
        p_value      = round(wt$p.value,   4),
        conclusion   = ifelse(wt$p.value < 0.05, "different", "ns")
      )
    })
  }
  
  # 玫瑰图
  mean_dirs <- map_dfr(valid_groups, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(!!group_col := g, mean_deg = cs$mean_deg)
  })
  
  rose_df <- sub_df %>%
    mutate(
      angle_deg = arrow_angle * 180 / pi,
      angle_deg = ifelse(angle_deg < 0, angle_deg + 360, angle_deg)
    )
  
  p_rose <- ggplot(rose_df, aes(x = angle_deg, fill = .data[[group_col]])) +
    geom_histogram(binwidth = 22.5, boundary = 0,
                   color = "white", linewidth = 0.3,
                   alpha = 0.75, position = "identity", closed = "left") +
    geom_vline(data = mean_dirs,
               aes(xintercept = mean_deg, color = .data[[group_col]]),
               linewidth = 0.9, linetype = "dashed", alpha = 0.85) +
    coord_polar(start = -pi / 2, direction = 1) +
    scale_x_continuous(
      limits = c(0, 360),
      breaks = seq(0, 315, by = 45),
      labels = c("0°\n(+CoIA1)", "45°", "90°\n(+CoIA2)",
                 "135°", "180°\n(-CoIA1)", "225°", "270°\n(-CoIA2)", "315°")
    ) +
    scale_fill_manual(values  = palette, name = group_label,
                      aesthetics = c("fill", "color")) +
    facet_wrap(reformulate(group_col),
               ncol = min(length(valid_groups), 3)) +
    theme_bw(base_size = 10) +
    labs(
      title    = sprintf("EXP L2-C: CoIA Arrow Direction by %s", group_label),
      subtitle = "Dashed line = mean direction; direction = axis of decoupling in CoIA space",
      x = NULL, y = "Count"
    ) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle    = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 9)
    )
  
  n_g   <- length(valid_groups)
  fname <- sprintf("analysis/output/figures/EXP_L2_Arrow_Direction_rose_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p_rose,
         width = min(4 + n_g * 3, 16), height = 5, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  
  list(desc = circ_desc, rayleigh = rayleigh_res, watson = watson_res)
}

res_circ_typology <- run_circular_analysis("Typology", "Typology", typology_pal)


# ==============================================================================
# ---- L2-D：spectral_entropy × Typology ----
# 复用：meta_exp（含 Typology，Levallois 已合并）、typology_pal、typology_levels
# 数据源：df_scar_all（方向）、df_morph_all（形态），均已在数据准备阶段过滤为 EXP+IM
# ==============================================================================

cat("\n---------- L2-D：spectral_entropy × Typology ----------\n")
cat("指标：刮痕方向数据（SE_direction）& 形态数据（SE_morphology）\n")
cat("检验：Kruskal-Wallis + Dunn 事后检验（FDR 校正）\n\n")

# ---- 构建 SE 分析数据框（只取纯实验标本，Levallois 已合并）----
se_df <- df_scar_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, SE_direction = spectral_entropy) %>%
  left_join(
    df_morph_all %>%
      filter(ID %in% exp_ids) %>%
      select(ID, SE_morphology = spectral_entropy),
    by = "ID"
  ) %>%
  left_join(
    meta_exp %>% select(ID, Typology),
    by = "ID"
  ) %>%
  filter(!is.na(Typology)) %>%
  # 过滤样本量 < 3 的类型
  group_by(Typology) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(Typology = fct_reorder(Typology, SE_direction, median, na.rm = TRUE))

cat("各 Typology 样本量：\n")
print(count(se_df, Typology))

# ---- 描述统计 ----
se_desc <- se_df %>%
  group_by(Typology) %>%
  summarise(
    n            = n(),
    mean_dir     = round(mean(SE_direction,   na.rm = TRUE), 4),
    median_dir   = round(median(SE_direction, na.rm = TRUE), 4),
    sd_dir       = round(sd(SE_direction,     na.rm = TRUE), 4),
    mean_morph   = round(mean(SE_morphology,   na.rm = TRUE), 4),
    median_morph = round(median(SE_morphology, na.rm = TRUE), 4),
    sd_morph     = round(sd(SE_morphology,     na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("\n==== SE 描述统计 ====\n")
print(se_desc)
write_csv(se_desc,
          here("analysis/data/derived_data/EXP_L2D_SE_desc_stats.csv"))

# ---- Kruskal-Wallis + Dunn 检验（内部辅助，不新增外部函数）----
run_kw_dunn_se <- function(df, y_col, label) {
  cat(sprintf("\n----- %s -----\n", label))
  
  # Kruskal-Wallis
  kw <- kruskal.test(reformulate("Typology", y_col), data = df)
  cat(sprintf("  Kruskal-Wallis: χ² = %.4f, df = %d, p = %.4f → %s\n",
              kw$statistic, kw$parameter, kw$p.value,
              ifelse(kw$p.value < 0.05, "组间有显著差异", "差异不显著")))
  
  # Dunn 事后检验（FSA，FDR 校正）
  dunn_raw <- dunnTest(
    x      = df[[y_col]],
    g      = df[["Typology"]],
    method = "holm"
  )$res                    # 取结果数据框
  
  # 拆分 "GroupA - GroupB" 为两列，添加显著性星号
  dunn <- dunn_raw %>%
    separate(Comparison, into = c("group1", "group2"),
             sep = " - ", remove = FALSE) %>%
    rename(statistic = Z, p = P.unadj, p.adj = P.adj) %>%
    mutate(
      p.adj.signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        p.adj < 0.10  ~ ".",
        TRUE          ~ "ns"
      )
    ) %>%
    select(Comparison, group1, group2, statistic, p, p.adj, p.adj.signif)
  
  cat("  Dunn 事后检验（FSA / Holm）：\n")
  print(dunn)
  list(kw = kw, dunn = dunn)
}

res_se_dir   <- run_kw_dunn_se(se_df %>% filter(!is.na(SE_direction)),
                               "SE_direction",  "SE（方向数据）× Typology")
res_se_morph <- run_kw_dunn_se(se_df %>% filter(!is.na(SE_morphology)),
                               "SE_morphology", "SE（形态数据）× Typology")

# 保存 Dunn 结果
bind_rows(
  res_se_dir$dunn   %>% mutate(source = "direction"),
  res_se_morph$dunn %>% mutate(source = "morphology")
) %>%
  write_csv(here("analysis/data/derived_data/EXP_L2D_SE_dunn_results.csv"))
cat("\n已保存：EXP_L2D_SE_dunn_results.csv\n")

# ---- 箱线图辅助函数（局部，仅供 L2-D 使用）----
make_se_boxplot <- function(df, y_col, y_label, kw_res, dunn_res,
                            title_suffix, fname) {
  
  # 显著对（p.adj < 0.05）
  sig_pairs <- dunn_res %>%
    filter(p.adj < 0.05) %>%
    select(group1, group2, p.adj.signif)
  
  # ---- 构建显著性标注坐标（手动绘制，兼容 FSA 输出格式）----
  sig_pairs <- dunn_res %>% filter(p.adj < 0.05)
  
  y_vals    <- df[[y_col]]
  y_max     <- max(y_vals, na.rm = TRUE)
  y_range   <- diff(range(y_vals, na.rm = TRUE))
  step      <- y_range * 0.10          # 每条连线的垂直间距
  
  # 将 Typology 水平映射到 x 轴数字位置
  type_lvls <- levels(df[["Typology"]])
  
  sig_annot <- if (nrow(sig_pairs) > 0) {
    sig_pairs %>%
      mutate(
        x1    = match(group1, type_lvls),
        x2    = match(group2, type_lvls),
        y_bar = y_max + step * row_number(),   # 各连线高度依次递增
        x_mid = (x1 + x2) / 2
      )
  } else {
    NULL
  }
  
  p <- ggplot(df,
              aes(x = Typology, y = .data[[y_col]],
                  fill = Typology, color = Typology)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha         = 0.55,
      linewidth     = 0.55,
      width         = 0.55
    ) +
    geom_jitter(
      width       = 0.14,
      size        = 2.2,
      alpha       = 0.80,
      stroke      = 0.3,
      shape       = 21,
      color       = "grey20",
      fill        = NA,
      show.legend = FALSE
    ) +
    # 各组样本量标注在箱体下方
    stat_summary(
      fun.data = function(x) {
        ypos <- min(x, na.rm = TRUE) - y_range * 0.04
        data.frame(y = ypos, ymin = ypos, ymax = ypos,
                   label = paste0("n=", length(x)))
      },
      geom        = "text",
      aes(label = after_stat(label)),
      size        = 2.8,
      color       = "grey45",
      vjust       = 1,
      show.legend = FALSE
    ) +
    # 显著性括号（手动绘制）
    {
      if (!is.null(sig_annot) && nrow(sig_annot) > 0) {
        tip <- y_range * 0.012
        list(
          # 水平主线
          geom_segment(
            data = sig_annot,
            aes(x = x1, xend = x2, y = y_bar, yend = y_bar),
            inherit.aes = FALSE, color = "grey30", linewidth = 0.4
          ),
          # 左竖线
          geom_segment(
            data = sig_annot,
            aes(x = x1, xend = x1, y = y_bar - tip, yend = y_bar),
            inherit.aes = FALSE, color = "grey30", linewidth = 0.4
          ),
          # 右竖线
          geom_segment(
            data = sig_annot,
            aes(x = x2, xend = x2, y = y_bar - tip, yend = y_bar),
            inherit.aes = FALSE, color = "grey30", linewidth = 0.4
          ),
          # 显著性标签
          geom_text(
            data = sig_annot,
            aes(x = x_mid, y = y_bar + y_range * 0.015,
                label = p.adj.signif),
            inherit.aes = FALSE, size = 3.2, color = "grey25"
          )
        )
      } else {
        NULL
      }
    } +
    # KW 结果角落标注
    annotate(
      "text",
      x     = Inf, y = Inf,
      label = sprintf("Kruskal-Wallis\nχ² = %.2f, p = %.3f",
                      kw_res$statistic, kw_res$p.value),
      hjust  = 1.05, vjust = 1.15,
      size   = 3.0,
      color  = ifelse(kw_res$p.value < 0.05, "#C0392B", "grey50")
    ) +
    scale_fill_manual(values  = typology_pal, guide = "none") +
    scale_color_manual(values = typology_pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.08))) +
    theme_bw(base_size = 11) +
    labs(
      title    = sprintf("EXP L2-D: Spectral Entropy by Typology — %s", title_suffix),
      subtitle = "Pairwise: Dunn test (Holm corrected); open circles = individual specimens",
      x        = "Stone Core Typology",
      y        = y_label
    ) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle      = element_text(hjust = 0.5, size = 8.5, color = "grey55"),
      axis.text.x        = element_text(angle = 25, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  ggsave(here(fname), plot = p,
         width  = max(6, length(typology_levels) * 1.3 + 2),
         height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  p
}

p_se_dir <- make_se_boxplot(
  df           = se_df %>% filter(!is.na(SE_direction)),
  y_col        = "SE_direction",
  y_label      = "Spectral Entropy (Scar Direction)",
  kw_res       = res_se_dir$kw,
  dunn_res     = res_se_dir$dunn,
  title_suffix = "Scar Direction",
  fname        = "analysis/output/figures/EXP_L2D_SE_Direction_Typology_boxplot.png"
)

p_se_morph <- make_se_boxplot(
  df           = se_df %>% filter(!is.na(SE_morphology)),
  y_col        = "SE_morphology",
  y_label      = "Spectral Entropy (Morphology)",
  kw_res       = res_se_morph$kw,
  dunn_res     = res_se_morph$dunn,
  title_suffix = "Morphology",
  fname        = "analysis/output/figures/EXP_L2D_SE_Morphology_Typology_boxplot.png"
)


# ==============================================================================
# ---- 保存衍生数据 ----
# ==============================================================================

# 圆形统计汇总
if (!is.null(res_circ_typology)) {
  circ_out <- left_join(
    res_circ_typology$desc,
    res_circ_typology$rayleigh,
    by = c("group_var", "group")
  )
  write_csv(circ_out, here("analysis/data/derived_data/EXP_L2_circular_stats.csv"))
  cat("已保存：EXP_L2_circular_stats.csv\n")
}

# 箭头统计汇总
scores_combined %>%
  select(ID, Typology, arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/EXP_L2_arrow_stats.csv"))
cat("已保存：EXP_L2_arrow_stats.csv\n")

# CoIA 坐标全表
scores_combined %>%
  write_csv(here("analysis/data/derived_data/EXP_CIA_scores_full.csv"))
cat("已保存：EXP_CIA_scores_full.csv\n")


# ==============================================================================
# ---- 汇总打印 ----
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  分析结果汇总（EXP）                                       ##\n")
cat("################################################################\n")

cat("\n【第一层：基线】\n")
cat(sprintf("  Mantel r = %.4f, p = %.3f  →  %s\n",
            mantel_global$statistic, mantel_global$signif,
            ifelse(mantel_global$signif < 0.05, "显著相关", "独立（ns）")))
cat(sprintf("  RV       = %.4f, p = %.3f  →  %s\n",
            coin_exp$RV, rv_test$pvalue,
            ifelse(rv_test$pvalue < 0.05, "显著协变", "独立（ns）")))

cat("\n【第二层 A：各类型独立 Mantel（Typology）】\n")
if (!is.null(l2_mantel) && nrow(l2_mantel) > 0) {
  n_sig <- sum(l2_mantel$p_value < 0.05, na.rm = TRUE)
  cat(sprintf("  共检验 %d 个类型（每类型独立置换检验），p < 0.05：%d 个\n",
              nrow(l2_mantel), n_sig))
  if (n_sig == 0) {
    cat("  → 所有类型内部均不显著，形态–技术独立性稳健\n")
  } else {
    cat("  → 显著类型：\n")
    l2_mantel %>%
      filter(p_value < 0.05) %>%
      select(Typology, n, mantel_r, p_value, significance) %>%
      mutate(across(c(mantel_r, p_value), ~ round(.x, 4))) %>%
      print()
  }
}

cat("\n【第二层 B：箭头长度】参见 EXP_L2_Arrow_Length_typology.png\n")
cat("【第二层 C：箭头方位】参见 EXP_L2_Arrow_Direction_rose_typology.png 及 EXP_L2_circular_stats.csv\n")

cat("\n【第二层 D：spectral_entropy × Typology】\n")
cat(sprintf("  SE（方向）Kruskal-Wallis: χ² = %.4f, p = %.4f → %s\n",
            res_se_dir$kw$statistic, res_se_dir$kw$p.value,
            ifelse(res_se_dir$kw$p.value < 0.05, "组间有显著差异", "差异不显著")))
cat(sprintf("  SE（形态）Kruskal-Wallis: χ² = %.4f, p = %.4f → %s\n",
            res_se_morph$kw$statistic, res_se_morph$kw$p.value,
            ifelse(res_se_morph$kw$p.value < 0.05, "组间有显著差异", "差异不显著")))
cat("  详见 EXP_L2D_SE_Direction_Typology_boxplot.png / EXP_L2D_SE_Morphology_Typology_boxplot.png\n")
cat("       EXP_L2D_SE_desc_stats.csv / EXP_L2D_SE_dunn_results.csv\n")

cat("\n\n========== EXP 分析全部完成 ==========\n")
cat("主要输出文件：\n")
cat("  EXP_L1_results.csv\n")
cat("  EXP_L2_grouped_mantel.csv\n")
cat("  EXP_L2_arrow_stats.csv\n")
cat("  EXP_L2_circular_stats.csv\n")
cat("  EXP_CIA_scores_full.csv\n")
cat("  EXP_L2D_SE_desc_stats.csv\n")
cat("  EXP_L2D_SE_dunn_results.csv\n")