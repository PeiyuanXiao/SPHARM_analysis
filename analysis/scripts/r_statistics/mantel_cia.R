# ==============================================================================
# mantel_cia_extended.R
# 形态 × 刮痕方向的关联分析（完整版）
#
# 分析流程：
#   Part A — Mantel 检验（原有）
#   Part B — 协惯量分析 CIA（原有）
#   Part C — PERMANOVA 分组显著性检验（新增）
#     1. 提取分组变量（层位 / 原料）
#     2. PERMANOVA（adonis2）：形态 ~ 层位 / 原料
#     3. PERMDISP：组内离散度齐性检验
#     4. 事后两两 PERMANOVA（FDR 校正）
#     5. 可视化：CAP 约束排序图
#   Part D — CIA 箭头长度解耦分析（新增）
#     1. 计算每个标本的箭头长度
#     2. 按层位 / 原料统计（Kruskal-Wallis + 两两 Wilcoxon）
#     3. 可视化（箱线图 + CIA 散点图箭头着色）
#     4. 保存数值结果
#
# 输入：
#   - analysis/data/derived_data/SPHARM_direction_filter.rds
#   - analysis/data/derived_data/SPHARM_morphology_filter.rds
#   （可选）若原料信息不在 ID 中，需提供含 ID + RawMaterial 列的外部表格
#
# 输出（新增）：
#   - analysis/output/figures/PERMANOVA_CAP_layer.png
#   - analysis/output/figures/PERMANOVA_CAP_rawmat.png
#   - analysis/output/figures/Arrow_Length_layer.png
#   - analysis/output/figures/Arrow_Length_rawmat.png
#   - analysis/output/figures/CIA_Arrow_colored.png
#   - analysis/data/derived_data/PERMANOVA_results.csv
#   - analysis/data/derived_data/Arrow_length_stats.csv
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
library(ggforce)       # stat_ellipse 替代方案，可选


# ==============================================================================
# ---- 读取上游数据（由 spharm_analysis.R 生成）----
# ==============================================================================

SPHARM_direction_filter  <- readRDS(here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
SPHARM_morphology_filter <- readRDS(here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))


# ==============================================================================
# ---- Part A：Mantel 检验（与原脚本相同，保持完整）----
# ==============================================================================

df_morph <- SPHARM_morphology_filter %>% rename(ID = ID)
df_scar  <- SPHARM_direction_filter  %>% rename(ID = ID)

common_ids <- intersect(df_morph$ID, df_scar$ID)
df_morph <- df_morph %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar  <- df_scar  %>% filter(ID %in% common_ids) %>% arrange(ID)

cat("==== 数据对齐检查 ====\n")
cat("共有标本数量：", length(common_ids), "\n")
cat("ID 是否完全匹配：", all(df_morph$ID == df_scar$ID), "\n\n")

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

morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

cosine_dist <- function(X) {
  X   <- as.matrix(X)
  sim <- X %*% t(X) / (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  as.dist(1 - sim)
}

D_morph <- cosine_dist(morph_power_clean)
D_scar  <- cosine_dist(scar_power_clean)

cat("==== 全局 Mantel 检验（形态 × 刮痕方向）====\n")
mantel_global <- mantel(D_morph, D_scar, method = "spearman", permutations = 999)
print(mantel_global)

run_cross_mantel <- function(X_single, D_target, from_label, var_label) {
  x <- X_single
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  d_x <- dist(scale(x))
  res <- mantel(d_x, D_target, method = "spearman", permutations = 999)
  tibble(from = from_label, var = var_label, r = res$statistic, p = res$signif)
}

mantel_M_to_scar <- map_dfr(
  colnames(morph_power_clean),
  ~ run_cross_mantel(morph_power_clean[[.x]], D_scar, "Scar Direction", .x)
)
mantel_S_to_morph <- map_dfr(
  colnames(scar_power_clean),
  ~ run_cross_mantel(scar_power_clean[[.x]], D_morph, "Morphology", .x)
)

mantel_cross <- bind_rows(mantel_M_to_scar, mantel_S_to_morph) %>%
  mutate(
    p_fdr        = p.adjust(p, method = "fdr"),
    significance = ifelse(p_fdr < 0.05, "P≤0.05", "P>0.05")
  )

# linkET 图（略，与原脚本相同）
# ggsave(here("analysis/output/figures/Mantel_Network.png"), ...)


# ==============================================================================
# ---- Part B：CIA（与原脚本相同，保持完整）----
# ==============================================================================

morph_arch <- morph_power_clean[!str_starts(rownames(morph_power_clean), "IM_"), ]
scar_arch  <- scar_power_clean[!str_starts(rownames(scar_power_clean),  "IM_"), ]
stopifnot(all(rownames(morph_arch) == rownames(scar_arch)))
cat("考古标本数量：", nrow(morph_arch), "\n")

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

morph_arch_ilr <- as.data.frame(ilr(replace_zeros(as.matrix(morph_arch))))
scar_arch_ilr  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_arch))))
rownames(morph_arch_ilr) <- rownames(morph_arch)
rownames(scar_arch_ilr)  <- rownames(scar_arch)
colnames(morph_arch_ilr) <- paste0("M_ilr", seq_len(ncol(morph_arch_ilr)))
colnames(scar_arch_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_arch_ilr)))

dudi_morph_arch <- dudi.pca(morph_arch_ilr, center = TRUE, scale = TRUE,
                            scannf = FALSE, nf = ncol(morph_arch_ilr))
dudi_scar_arch  <- dudi.pca(scar_arch_ilr,  center = TRUE, scale = TRUE,
                            scannf = FALSE, nf = ncol(scar_arch_ilr))

coin_arch <- coinertia(dudi_morph_arch, dudi_scar_arch, scannf = FALSE, nf = 2)

cat("RV 系数（仅考古标本）：", round(coin_arch$RV, 4), "\n")
set.seed(42)
rv_test_arch <- randtest(coin_arch, nrepet = 9999)
print(rv_test_arch)

# CIA 坐标（供后续 Part D 使用）
scores_morph <- as.data.frame(coin_arch$lX) %>% rownames_to_column("ID")
scores_scar  <- as.data.frame(coin_arch$lY) %>% rownames_to_column("ID")

cia_inertia <- coin_arch$eig / sum(coin_arch$eig) * 100


# ==============================================================================
# ---- Part C：PERMANOVA 分组显著性检验（新增）----
# ==============================================================================

# ------------------------------------------------------------------------------
# C-0：构建仅含考古标本的距离矩阵
# ------------------------------------------------------------------------------

arch_ids <- rownames(morph_arch)

# 从全局距离矩阵中提取考古标本子矩阵
D_morph_arch <- as.dist(
  as.matrix(D_morph)[arch_ids, arch_ids]
)
D_scar_arch <- as.dist(
  as.matrix(D_scar)[arch_ids, arch_ids]
)


# ------------------------------------------------------------------------------
# C-1：提取分组变量
# ------------------------------------------------------------------------------
# 层位：从 ID 中提取 L2 / L3 / L4
# 原料：
#   方案 A——若 ID 中包含原料编码（如 "FL"=燧石, "QZ"=石英岩），在此修改正则
#   方案 B——若原料信息在外部表中，取消下方注释并载入
# ------------------------------------------------------------------------------

meta_arch <- tibble(ID = arch_ids) %>%
  mutate(
    layer = case_when(
      str_detect(ID, "L2") ~ "Layer 2",
      str_detect(ID, "L3") ~ "Layer 3",
      str_detect(ID, "L4") ~ "Layer 4",
      TRUE                  ~ "Other"
    ),
    raw_material = case_when(
      str_detect(ID, "FL") ~ "Flint",
      str_detect(ID, "QZ") ~ "Quartzite",
      str_detect(ID, "CH") ~ "Chert",
      TRUE                  ~ "Other"
    )
  )

cat("==== 分组诊断 ====\n")
cat("层位分布：\n"); print(table(meta_arch$layer))
cat("原料分布：\n"); print(table(meta_arch$raw_material))

# ------------------------------------------------------------------------------
# 辅助函数：安全筛选——同时满足（1）每组 ≥ 3 件 AND（2）至少 2 个组
# 返回过滤后的 meta 子集，若不满足条件返回 NULL 并打印提示
# ------------------------------------------------------------------------------

safe_filter_groups <- function(meta_df, group_col, min_n = 3) {
  counts <- table(meta_df[[group_col]])
  valid  <- names(counts[counts >= min_n])
  
  if (length(valid) < 2) {
    cat(sprintf(
      "  [跳过] %s：有效组数不足（需 ≥ 2 组，每组 ≥ %d 件）。当前：%s\n",
      group_col, min_n,
      paste(names(counts), counts, sep = "=", collapse = ", ")
    ))
    return(NULL)
  }
  
  meta_df %>% filter(.data[[group_col]] %in% valid)
}

meta_layer  <- safe_filter_groups(meta_arch, "layer")
meta_rawmat <- safe_filter_groups(meta_arch, "raw_material")


# ------------------------------------------------------------------------------
# 提取子距离矩阵（仅在对应 meta 不为 NULL 时执行）
# ------------------------------------------------------------------------------

extract_subdist <- function(D_full, ids) {
  as.dist(as.matrix(D_full)[ids, ids])
}

if (!is.null(meta_layer)) {
  D_morph_layer <- extract_subdist(D_morph_arch, meta_layer$ID)
  D_scar_layer  <- extract_subdist(D_scar_arch,  meta_layer$ID)
}

if (!is.null(meta_rawmat)) {
  D_morph_rawmat <- extract_subdist(D_morph_arch, meta_rawmat$ID)
  D_scar_rawmat  <- extract_subdist(D_scar_arch,  meta_rawmat$ID)
}

# --- 方案 B：若原料在外部表格中 ---
# raw_material_table <- read_csv(here("analysis/data/raw_data/raw_material.csv"))
#   # 表格需含 ID + RawMaterial 两列
# meta_arch <- meta_arch %>%
#   left_join(raw_material_table, by = "ID") %>%
#   rename(raw_material = RawMaterial)

cat("\n==== 分组变量分布 ====\n")
cat("层位：\n"); print(table(meta_arch$layer))
cat("原料：\n"); print(table(meta_arch$raw_material))

# 过滤掉样本量过小的组（< 3 件），避免 PERMANOVA 不稳定
valid_layers  <- names(table(meta_arch$layer)[table(meta_arch$layer) >= 3])
valid_rawmats <- names(table(meta_arch$raw_material)[table(meta_arch$raw_material) >= 3])

meta_layer  <- meta_arch %>% filter(layer %in% valid_layers)
meta_rawmat <- meta_arch %>% filter(raw_material %in% valid_rawmats)

# 对应子距离矩阵
D_morph_layer  <- as.dist(as.matrix(D_morph_arch)[meta_layer$ID,  meta_layer$ID])
D_scar_layer   <- as.dist(as.matrix(D_scar_arch)[meta_layer$ID,   meta_layer$ID])
D_morph_rawmat <- as.dist(as.matrix(D_morph_arch)[meta_rawmat$ID, meta_rawmat$ID])
D_scar_rawmat  <- as.dist(as.matrix(D_scar_arch)[meta_rawmat$ID,  meta_rawmat$ID])


# ------------------------------------------------------------------------------
# C-2：全局 PERMANOVA（adonis2）
# ------------------------------------------------------------------------------

run_permanova <- function(dist_mat, group_vec, group_name, domain_name,
                          n_perm = 9999) {
  
  # 检查分组水平
  group_vec <- as.character(group_vec)
  n_levels  <- length(unique(group_vec))
  if (n_levels < 2) {
    cat(sprintf("  [跳过] %s ~ %s：只有 %d 个分组水平\n",
                domain_name, group_name, n_levels))
    return(NULL)
  }
  
  df_tmp <- data.frame(group = factor(group_vec))
  
  # add = "lingoes"：对非欧距离矩阵的负特征值做 Lingoes 加性校正
  res <- adonis2(dist_mat ~ group, data = df_tmp,
                 permutations = n_perm,
                 add = "lingoes")
  
  cat(sprintf("\n----- %s ~ %s -----\n", domain_name, group_name))
  print(res)
  
  tibble(
    domain   = domain_name,
    grouping = group_name,
    R2       = res$R2[1],
    F_value  = res$F[1],
    p_value  = res$`Pr(>F)`[1],
    df_group = res$Df[1],
    df_resid = res$Df[2]
  )
}

cat("\n\n========== PERMANOVA 结果 ==========\n")

permanova_results <- bind_rows(
  if (!is.null(meta_layer)) {
    bind_rows(
      run_permanova(D_morph_layer, meta_layer$layer, "Layer", "Morphology"),
      run_permanova(D_scar_layer,  meta_layer$layer, "Layer", "Scar Direction")
    )
  },
  if (!is.null(meta_rawmat)) {
    bind_rows(
      run_permanova(D_morph_rawmat, meta_rawmat$raw_material, "Raw Material", "Morphology"),
      run_permanova(D_scar_rawmat,  meta_rawmat$raw_material, "Raw Material", "Scar Direction")
    )
  }
)

cat("\n==== PERMANOVA 汇总 ====\n")
permanova_results %>%
  mutate(across(c(R2, F_value, p_value), ~ round(.x, 4))) %>%
  print()


# ------------------------------------------------------------------------------
# C-3：PERMDISP——组内离散度齐性检验
# ------------------------------------------------------------------------------
# 目的：区分"组间均值差异"（PERMANOVA 检验目标）与"组内离散度差异"
# 若 PERMDISP 显著，需在论文中说明：组间差异部分来源于离散度不均

run_permdisp <- function(dist_mat, group_vec, group_name, domain_name) {
  
  group_vec <- as.character(group_vec)
  if (length(unique(group_vec)) < 2) {
    cat(sprintf("  [跳过] %s ~ %s：分组水平不足\n", domain_name, group_name))
    return(NULL)
  }
  
  # Lingoes 校正后再做 betadisper
  D_corrected <- as.dist(
    as.matrix(dist_mat) +
      abs(min(cmdscale(dist_mat, eig = TRUE)$eig[
        cmdscale(dist_mat, eig = TRUE)$eig < 0
      ], 0))
  )
  
  bd  <- betadisper(D_corrected, factor(group_vec))
  res <- permutest(bd, permutations = 9999)
  
  cat(sprintf("\n----- %s ~ %s -----\n", domain_name, group_name))
  print(res)
  
  tibble(
    domain   = domain_name,
    grouping = group_name,
    F_value  = res$tab$F[1],
    p_value  = res$tab$`Pr(>F)`[1]
  )
}

cat("\n\n========== PERMDISP 结果 ==========\n")

permdisp_results <- bind_rows(
  if (!is.null(meta_layer)) {
    bind_rows(
      run_permdisp(D_morph_layer, meta_layer$layer, "Layer", "Morphology"),
      run_permdisp(D_scar_layer,  meta_layer$layer, "Layer", "Scar Direction")
    )
  },
  if (!is.null(meta_rawmat)) {
    bind_rows(
      run_permdisp(D_morph_rawmat, meta_rawmat$raw_material, "Raw Material", "Morphology"),
      run_permdisp(D_scar_rawmat,  meta_rawmat$raw_material, "Raw Material", "Scar Direction")
    )
  }
)

# ------------------------------------------------------------------------------
# C-4：事后两两 PERMANOVA（FDR 校正）
# ------------------------------------------------------------------------------
# 若全局 PERMANOVA 显著，才有必要做事后两两比较

pairwise_permanova <- function(dist_mat, group_vec, group_name,
                               domain_name, n_perm = 9999) {
  
  group_vec <- as.character(group_vec)
  groups    <- unique(group_vec)
  if (length(groups) < 2) return(NULL)
  
  pairs <- combn(groups, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    idx   <- group_vec %in% pair
    d_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
    g_sub <- factor(group_vec[idx])
    df_tmp <- data.frame(group = g_sub)
    
    res <- adonis2(d_sub ~ group, data = df_tmp,
                   permutations = n_perm,
                   add = "lingoes")
    
    tibble(
      domain  = domain_name,
      grouping = group_name,
      group1  = pair[1],
      group2  = pair[2],
      R2      = res$R2[1],
      F_value = res$F[1],
      p_raw   = res$`Pr(>F)`[1]
    )
  })
}

cat("\n\n========== 事后两两 PERMANOVA ==========\n")

pairwise_results <- bind_rows(
  if (!is.null(meta_layer)) {
    bind_rows(
      pairwise_permanova(D_morph_layer, meta_layer$layer, "Layer", "Morphology"),
      pairwise_permanova(D_scar_layer,  meta_layer$layer, "Layer", "Scar Direction")
    )
  },
  if (!is.null(meta_rawmat)) {
    bind_rows(
      pairwise_permanova(D_morph_rawmat, meta_rawmat$raw_material, "Raw Material", "Morphology"),
      pairwise_permanova(D_scar_rawmat,  meta_rawmat$raw_material, "Raw Material", "Scar Direction")
    )
  }
) %>%
  mutate(
    p_fdr        = p.adjust(p_raw, method = "fdr"),
    significance = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

cat("\n两两 PERMANOVA（FDR 校正）：\n")
pairwise_results %>%
  mutate(across(c(R2, F_value, p_raw, p_fdr), ~ round(.x, 4))) %>%
  print(n = Inf)


# ------------------------------------------------------------------------------
# C-5：CAP 约束排序图
# ------------------------------------------------------------------------------

plot_cap <- function(dist_mat, meta_df, group_col,
                     title_str, palette_vec) {
  
  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  # ---- 分组数据 ----
  df_group <- data.frame(group = meta_df[[group_col]])
  
  # ---- CAP ----
  cap_res  <- capscale(dist_mat ~ group, data = df_group)
  
  # ---- 提取坐标 ----
  all_scores <- scores(cap_res, display = "sites", choices = 1:2)
  
  cap_scores <- as.data.frame(all_scores) %>%
    mutate(ID    = meta_df$ID,
           group = meta_df[[group_col]])
  
  # ---- 动态列名 ----
  x_col <- colnames(cap_scores)[1]
  y_col <- colnames(cap_scores)[2]
  
  # ---- 特征值 ----
  eig_all <- eigenvals(cap_res)
  eig_pos <- eig_all[eig_all > 0]
  
  eig_df <- data.frame(name = names(eig_all),
                       value = as.numeric(eig_all))
  
  pct_df <- eig_df %>%
    filter(name %in% c(x_col, y_col)) %>%
    mutate(pct = round(value / sum(eig_pos) * 100, 1))
  
  x_pct <- pct_df$pct[pct_df$name == x_col]
  y_pct <- pct_df$pct[pct_df$name == y_col]
  
  x_lab <- sprintf("%s (%.1f%%)", x_col, x_pct)
  y_lab <- sprintf("%s (%.1f%%)", y_col, y_pct)
  
  # ---- 作图 ----
  p <- ggplot(cap_scores,
              aes(x = .data[[x_col]],
                  y = .data[[y_col]],
                  color = group)) +
    
    stat_ellipse(aes(fill = group),
                 geom = "polygon",
                 alpha = 0.08,
                 level = 0.9,
                 linetype = "dashed",
                 linewidth = 0.5) +
    
    geom_hline(yintercept = 0,
               linetype = "dotted",
               color = "grey60",
               linewidth = 0.4) +
    
    geom_vline(xintercept = 0,
               linetype = "dotted",
               color = "grey60",
               linewidth = 0.4) +
    
    geom_point(size = 3, alpha = 0.85) +
    
    geom_text_repel(aes(label = ID),
                    size = 2.2,
                    color = "grey40",
                    max.overlaps = 8,
                    show.legend = FALSE) +
    
    scale_color_manual(values = palette_vec,
                       name = group_col,
                       aesthetics = c("color", "fill")) +
    
    theme_bw(base_size = 10) +
    
    labs(title = title_str,
         x = x_lab,
         y = y_lab) +
    
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "right"
    )
  
  return(p)
}

layer_pal <- c("Layer 2" = "#E6B89C",
               "Layer 3" = "#A1C2E6",
               "Layer 4" = "#FFBAE0",
               "Other"   = "#BBBBBB")

rawmat_pal <- c("Flint"     = "#8ECFC9",
                "Quartzite" = "#FFBE7A",
                "Chert"     = "#FA7F6F",
                "Other"     = "#BBBBBB")

cap_morph_layer <- plot_cap(D_morph_layer, meta_layer, "layer",
                            "Morphology — CAP by Layer", layer_pal)

cap_scar_layer  <- plot_cap(D_scar_layer, meta_layer, "layer",
                            "Scar Direction — CAP by Layer", layer_pal)

library(patchwork)

p_cap_layer <- (cap_morph_layer | cap_scar_layer) +
  plot_annotation(
    title = "Constrained Analysis of Principal Coordinates (CAP) — Layer",
    theme = theme(plot.title = element_text(face = "bold",
                                            size = 12,
                                            hjust = 0.5))
  )

ggsave(
  here::here("analysis/output/figures/PERMANOVA_CAP_layer.png"),
  plot = p_cap_layer,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)

cat("图已保存：PERMANOVA_CAP_layer.png\n")

if (!is.null(meta_rawmat) &&
    length(unique(meta_rawmat$raw_material)) > 1) {
  
  cap_morph_rawmat <- plot_cap(D_morph_rawmat, meta_rawmat, "raw_material",
                               "Morphology — CAP by Raw Material", rawmat_pal)
  
  cap_scar_rawmat  <- plot_cap(D_scar_rawmat, meta_rawmat, "raw_material",
                               "Scar Direction — CAP by Raw Material", rawmat_pal)
  
  p_cap_rawmat <- (cap_morph_rawmat | cap_scar_rawmat) +
    plot_annotation(
      title = "Constrained Analysis of Principal Coordinates (CAP) — Raw Material",
      theme = theme(plot.title = element_text(face = "bold",
                                              size = 12,
                                              hjust = 0.5))
    )
  
  ggsave(
    here::here("analysis/output/figures/PERMANOVA_CAP_rawmat.png"),
    plot = p_cap_rawmat,
    width = 14,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("图已保存：PERMANOVA_CAP_rawmat.png\n")
  
} else {
  cat("[跳过] 原料 CAP 图：分组不足或数据为空\n")
}


# ==============================================================================
# ---- Part D：CIA 箭头长度解耦分析（新增）----
# ==============================================================================

# ------------------------------------------------------------------------------
# D-1：计算 CIA 箭头长度并合并元数据
# ------------------------------------------------------------------------------

scores_combined <- left_join(
  scores_morph %>% select(ID, Axis1_M = AxcX1, Axis2_M = AxcX2),
  scores_scar  %>% select(ID, Axis1_S = AxcY1, Axis2_S = AxcY2),
  by = "ID"
) %>%
  mutate(
    # 箭头长度 = 形态坐标与技术坐标之间的欧氏距离
    arrow_length = sqrt((Axis1_M - Axis1_S)^2 + (Axis2_M - Axis2_S)^2),
    
    # 箭头方向（弧度，可用于判断解耦方向）
    arrow_angle = atan2(Axis2_S - Axis2_M, Axis1_S - Axis1_M),
    
    # 分组：层位
    layer = case_when(
      str_detect(ID, "L2") ~ "Layer 2",
      str_detect(ID, "L3") ~ "Layer 3",
      str_detect(ID, "L4") ~ "Layer 4",
      TRUE                  ~ "Other"
    )
  ) %>%
  # 合并原料信息
  left_join(meta_arch %>% select(ID, raw_material), by = "ID")

cat("\n==== 箭头长度描述统计（按层位）====\n")
scores_combined %>%
  group_by(layer) %>%
  summarise(
    n      = n(),
    mean   = round(mean(arrow_length),  4),
    median = round(median(arrow_length), 4),
    sd     = round(sd(arrow_length),     4),
    min    = round(min(arrow_length),    4),
    max    = round(max(arrow_length),    4)
  ) %>%
  print()

cat("\n==== 箭头长度描述统计（按原料）====\n")
scores_combined %>%
  group_by(raw_material) %>%
  summarise(
    n      = n(),
    mean   = round(mean(arrow_length),  4),
    median = round(median(arrow_length), 4),
    sd     = round(sd(arrow_length),     4),
    min    = round(min(arrow_length),    4),
    max    = round(max(arrow_length),    4)
  ) %>%
  print()


# ------------------------------------------------------------------------------
# D-2：Kruskal-Wallis 检验 + 事后两两 Wilcoxon（FDR 校正）
# ------------------------------------------------------------------------------

# 原料组间检验：仅在有有效分组时执行
if (length(valid_arrow_rawmats) >= 2) {
  
  arrow_rawmat_df <- scores_combined %>%
    filter(raw_material %in% valid_arrow_rawmats)
  
  cat("\n----- 原料 × 箭头长度：Kruskal-Wallis -----\n")
  kw_rawmat <- kruskal.test(arrow_length ~ raw_material, data = arrow_rawmat_df)
  print(kw_rawmat)
  
  cat("\n----- 原料：事后两两 Wilcoxon（FDR 校正）-----\n")
  pairwise_rawmat <- pairwise.wilcox.test(
    arrow_rawmat_df$arrow_length,
    arrow_rawmat_df$raw_material,
    p.adjust.method = "fdr",
    exact = FALSE
  )
  print(pairwise_rawmat)
  
} else {
  cat("\n  [跳过] 原料箭头长度检验：有效分组不足（当前所有标本均归为 'Other'）\n")
  cat("  请检查 meta_arch 中 raw_material 的 case_when 规则是否匹配实际 ID 格式\n")
  arrow_rawmat_df <- NULL
  kw_rawmat       <- NULL
  pairwise_rawmat <- NULL
}


# ------------------------------------------------------------------------------
# D-3：箭头长度可视化
# ------------------------------------------------------------------------------

# 从 Wilcoxon 矩阵提取显著性标注的辅助函数
get_sig_label <- function(p_mat, g1, g2) {
  p <- tryCatch(p_mat[g1, g2], error = function(e) p_mat[g2, g1])
  if (is.na(p))  return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# --- 层位箱线图 ---
p_arrow_layer <- ggplot(
  arrow_layer_df,
  aes(x = layer, y = arrow_length, fill = layer)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.5) +
  geom_jitter(aes(color = layer), width = 0.15, size = 2.5,
              alpha = 0.75, show.legend = FALSE) +
  geom_text_repel(
    data = arrow_layer_df %>%
      group_by(layer) %>%
      filter(arrow_length == max(arrow_length)),
    aes(label = ID), size = 2.4, color = "grey40",
    max.overlaps = 10, show.legend = FALSE
  ) +
  annotate(
    "text", x = 1.5, y = max(arrow_layer_df$arrow_length) * 1.02,
    label = sprintf("Kruskal-Wallis: χ² = %.2f, p = %.3f",
                    kw_layer$statistic, kw_layer$p.value),
    size = 3, color = "grey30", hjust = 0.5
  ) +
  scale_fill_manual(values  = layer_pal, guide = "none") +
  scale_color_manual(values = layer_pal, guide = "none") +
  theme_bw(base_size = 10) +
  labs(
    title    = "CIA Arrow Length by Layer",
    subtitle = "Longer arrows = greater morphology–technique decoupling",
    x        = NULL,
    y        = "Arrow length (CIA space)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50")
  )

ggsave(
  here("analysis/output/figures/Arrow_Length_layer.png"),
  plot = p_arrow_layer, width = 7, height = 6, dpi = 300, bg = "white"
)
cat("图已保存：Arrow_Length_layer.png\n")

# --- 原料箱线图 ---
if (!is.null(arrow_rawmat_df)) {
  
  p_arrow_rawmat <- ggplot(
    arrow_rawmat_df,
    aes(x = raw_material, y = arrow_length, fill = raw_material)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.5) +
    geom_jitter(aes(color = raw_material), width = 0.15, size = 2.5,
                alpha = 0.75, show.legend = FALSE) +
    annotate(
      "text", x = 1.5, y = max(arrow_rawmat_df$arrow_length) * 1.02,
      label = sprintf("Kruskal-Wallis: χ² = %.2f, p = %.3f",
                      kw_rawmat$statistic, kw_rawmat$p.value),
      size = 3, color = "grey30", hjust = 0.5
    ) +
    scale_fill_manual(values  = rawmat_pal, guide = "none") +
    scale_color_manual(values = rawmat_pal, guide = "none") +
    theme_bw(base_size = 10) +
    labs(
      title    = "CIA Arrow Length by Raw Material",
      subtitle = "Longer arrows = greater morphology–technique decoupling",
      x = NULL, y = "Arrow length (CIA space)"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50")
    )
  
  ggsave(
    here("analysis/output/figures/Arrow_Length_rawmat.png"),
    plot = p_arrow_rawmat, width = 7, height = 6, dpi = 300, bg = "white"
  )
  cat("图已保存：Arrow_Length_rawmat.png\n")
  
} else {
  cat("  [跳过] 原料箭头长度图：无有效分组\n")
}


# ------------------------------------------------------------------------------
# D-4：CIA 散点图（箭头按长度着色）
# ------------------------------------------------------------------------------
# 在原有 CIA 双标图基础上，将箭头颜色映射为解耦程度，
# 点的形状区分层位，方便同时观察两个维度

p_cia_arrow_colored <- ggplot() +
  # 连接线（颜色 = 箭头长度）
  geom_segment(
    data = scores_combined,
    aes(x = Axis1_M, y = Axis2_M,
        xend = Axis1_S, yend = Axis2_S,
        color = arrow_length),
    linewidth = 0.8, alpha = 0.8,
    arrow = arrow(length = unit(0.12, "cm"), type = "closed")
  ) +
  # 形态坐标（圆点）
  geom_point(data = scores_combined,
             aes(x = Axis1_M, y = Axis2_M, shape = layer),
             size = 3, fill = "white", color = "grey30",
             stroke = 1.2, alpha = 0.9) +
  # 技术坐标（三角）
  geom_point(data = scores_combined,
             aes(x = Axis1_S, y = Axis2_S, shape = layer),
             size = 3, alpha = 0.9, color = "grey50") +
  # 最长箭头（解耦最大）标注 ID
  geom_text_repel(
    data = scores_combined %>%
      slice_max(arrow_length, n = 5),
    aes(x = (Axis1_M + Axis1_S) / 2,
        y = (Axis2_M + Axis2_S) / 2,
        label = ID),
    size = 2.5, color = "grey30", max.overlaps = 15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  scale_color_viridis_c(
    option = "C", name = "Arrow length\n(decoupling)",
    direction = -1
  ) +
  scale_shape_manual(
    values = c("Layer 2" = 21, "Layer 3" = 22,
               "Layer 4" = 24, "Other"   = 20),
    name = "Layer\n(● morph / △ scar)"
  ) +
  theme_bw(base_size = 10) +
  labs(
    title    = "CIA Biplot — Arrow Length as Morphology–Technique Decoupling",
    subtitle = "Warm colors = high decoupling; Cool colors = high coupling",
    x = sprintf("CIA Axis 1 (%.1f%%)", cia_inertia[1]),
    y = sprintf("CIA Axis 2 (%.1f%%)", cia_inertia[2])
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
    legend.position = "right"
  )

ggsave(
  here("analysis/output/figures/CIA_Arrow_colored.png"),
  plot = p_cia_arrow_colored, width = 10, height = 8, dpi = 300, bg = "white"
)
cat("图已保存：CIA_Arrow_colored.png\n")


# ==============================================================================
# ---- 保存全部数值结果 ----
# ==============================================================================

# PERMANOVA 汇总表
bind_rows(
  permanova_results %>%
    mutate(test = "PERMANOVA", statistic = F_value) %>%
    select(test, domain, grouping, R2, statistic, p_value),
  permdisp_results %>%
    mutate(test = "PERMDISP", R2 = NA,
           statistic = F_value, p_value = p_value) %>%
    select(test, domain, grouping, R2, statistic, p_value)
) %>%
  write_csv(here("analysis/data/derived_data/PERMANOVA_results.csv"))

# 箭头长度汇总表
scores_combined %>%
  select(ID, layer, raw_material, arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/Arrow_length_stats.csv"))

cat("\n结果已保存：\n")
cat("  PERMANOVA_results.csv\n")
cat("  Arrow_length_stats.csv\n")

cat("\n\n========== 全部分析完成 ==========\n")