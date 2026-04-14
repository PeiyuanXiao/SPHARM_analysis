# ==============================================================================
# SDG_cores_statistics.R
#
# 第一层：整体 Mantel + CIA
#   形态与片疤方向是否整体独立？
#   附：CIA 箭头长度）+ 箭头方位描述
#
# 第二层：分组 Mantel（层位 × 原料 × 类型）
#   这种独立性是否在所有条件下都成立？
#   联合证据：
#     - 分组 Mantel r
#     - 箭头长度分组差异
#     - 箭头方位分组差异
#
# 第三层：PERMANOVA（层位 × 原料 × 类型）
#   形态和片疤方向各自的分组结构
#
# 输入：
#   - analysis/data/derived_data/SPHARM_direction_filter.rds
#   - analysis/data/derived_data/SPHARM_morphology_filter.rds
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#
# 输出：
#   figures/L1_Mantel_Network.png
#   figures/L1_CIA_Biplot.png
#   figures/L2_Mantel_Grouped_dumbbell.png
#   figures/L2_Mantel_Grouped_heatmap.png
#   figures/L2_Arrow_Length_layer.png
#   figures/L2_Arrow_Length_rawmat.png       （有效分组时）
#   figures/L2_Arrow_Length_coretype.png     （有效分组时）
#   figures/L2_Arrow_Direction_rose_layer.png
#   figures/L2_Arrow_Direction_rose_rawmat.png   （有效分组时）
#   figures/L2_Arrow_Direction_rose_coretype.png （有效分组时）
#   figures/L3_PERMANOVA_R2.png
#   figures/L3_PERMANOVA_CAP_layer.png
#   figures/L3_PERMANOVA_CAP_raw_material.png
#   figures/L3_PERMANOVA_CAP_core_type.png
#   derived_data/L1_results.csv
#   derived_data/L2_grouped_mantel.csv
#   derived_data/L2_arrow_stats.csv
#   derived_data/L2_circular_stats.csv
#   derived_data/L3_permanova.csv
#   derived_data/CIA_scores_full.csv
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
    row_i       <- X[i, ]
    zero_idx    <- row_i == 0
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

# 圆形统计：单组均值方向和集中度
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


# ==============================================================================
# ---- 数据准备（所有层共用）----
# ==============================================================================

SPHARM_direction_filter  <- readRDS(here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
SPHARM_morphology_filter <- readRDS(here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

df_morph <- SPHARM_morphology_filter
df_scar  <- SPHARM_direction_filter

common_ids <- intersect(df_morph$ID, df_scar$ID)
df_morph   <- df_morph %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar    <- df_scar  %>% filter(ID %in% common_ids) %>% arrange(ID)

cat("==== 数据对齐 ====\n")
cat("共有标本：", length(common_ids), "；ID 完全匹配：",
    all(df_morph$ID == df_scar$ID), "\n\n")

morph_power <- df_morph %>%
  select(power_l1:power_l4) %>%
  rename_with(~ paste0("M", 1:4)) %>%
  as.data.frame()
scar_power <- df_scar %>%
  select(power_l1:power_l4) %>%
  rename_with(~ paste0("S", 1:4)) %>%
  as.data.frame()
rownames(morph_power) <- df_morph$ID
rownames(scar_power)  <- df_scar$ID

morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

D_morph <- cosine_dist(morph_power_clean)
D_scar  <- cosine_dist(scar_power_clean)

morph_arch <- morph_power_clean[!str_starts(rownames(morph_power_clean), "IM_"), ]
scar_arch  <- scar_power_clean[!str_starts(rownames(scar_power_clean),  "IM_"), ]
stopifnot(all(rownames(morph_arch) == rownames(scar_arch)))
arch_ids <- rownames(morph_arch)
cat("考古标本数量：", length(arch_ids), "\n")

D_morph_arch <- extract_subdist(D_morph, arch_ids)
D_scar_arch  <- extract_subdist(D_scar,  arch_ids)

# 外部元数据
core_meta_raw <- read_excel(
  "H:/SPHARM_analysis/analysis/data/raw_data/SDG_core_metric.xlsx"
)

cat("\n==== 外部表格诊断 ====\n")
cat("列名："); print(colnames(core_meta_raw))
cat("\n前 3 行：\n"); print(head(core_meta_raw, 3))

core_meta <- core_meta_raw %>%
  select(
    ID           = ID,
    raw_material = Raw_mat,
    core_type    = Core_type_Li_merged
  ) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))

meta_arch <- tibble(ID = arch_ids) %>%
  mutate(
    layer = case_when(
      str_detect(ID, "L2") ~ "Layer 2",
      str_detect(ID, "L3") ~ "Layer 3",
      str_detect(ID, "L4") ~ "Layer 4",
      TRUE                  ~ "Other"
    )
  ) %>%
  left_join(core_meta, by = "ID")

cat("\n==== 元数据合并诊断 ====\n")
cat("层位：\n");     print(table(meta_arch$layer,        useNA = "ifany"))
cat("原料：\n");     print(table(meta_arch$raw_material,  useNA = "ifany"))
cat("石核类型：\n"); print(table(meta_arch$core_type,     useNA = "ifany"))
unmatched <- meta_arch %>% filter(is.na(raw_material)) %>% pull(ID)
if (length(unmatched) > 0) {
  cat("未匹配 ID（请检查格式）：\n"); print(unmatched)
}

meta_layer    <- safe_filter_groups(meta_arch, "layer")
meta_rawmat   <- safe_filter_groups(meta_arch, "raw_material")
meta_coretype <- safe_filter_groups(meta_arch, "core_type")

# 调色板（请按实际类别名称调整）
layer_pal <- c(
  "Layer 2" = "#E6B89C",
  "Layer 3" = "#A1C2E6",
  "Layer 4" = "#FFBAE0"
)
rawmat_pal <- c(
  "chert"     = "#8ECFC9",
  "sandstone" = "#FFBE7A"
)
coretype_pal <- c(
  "Unifacial_unidirection"  = "#7EB8C9",
  "Unifacial_centripetal"   = "#6271A1",
  "Bifacial_adjacent"       = "#C6DEA4",
  "Bifacial_independent"    = "#90C49A",
  "Bifacial_centripetal"    = "#609988",
  "Multifacial"             = "#D4A5A3",
  "Core_on_flake"           = "#D6D6D6",
  "Handaxe"                 = "#EDF0A3",
  "Pick"                    = "#F7E09E"
)


# ==============================================================================
# ========== 第一层：整体 Mantel + CIA（基线）==========
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  第一层：整体 Mantel + CIA — 建立基线                      ##\n")
cat("################################################################\n")


# ------------------------------------------------------------------------------
# L1-1：全局 Mantel 检验 + linkET 网络图
# ------------------------------------------------------------------------------

cat("\n==== L1-1：全局 Mantel 检验 ====\n")
mantel_global <- mantel(D_morph_arch, D_scar_arch,
                        method = "spearman", permutations = 9999)
print(mantel_global)

run_cross_mantel <- function(X_single, D_target, from_label, var_label,
                             n_perm = 999) {
  if (sd(X_single, na.rm = TRUE) == 0) return(NULL)
  d_x <- dist(scale(X_single))
  res <- mantel(d_x, D_target, method = "spearman", permutations = n_perm)
  tibble(from = from_label, var = var_label,
         r = res$statistic, p = res$signif)
}

morph_arch_df <- as.data.frame(morph_arch)
scar_arch_df  <- as.data.frame(scar_arch)

mantel_cross_l1 <- bind_rows(
  map_dfr(colnames(morph_arch_df),
          ~ run_cross_mantel(morph_arch_df[[.x]], D_scar_arch,
                             "Scar Direction", .x)),
  map_dfr(colnames(scar_arch_df),
          ~ run_cross_mantel(scar_arch_df[[.x]], D_morph_arch,
                             "Morphology", .x))
) %>%
  mutate(
    p_fdr        = p.adjust(p, method = "fdr"),
    significance = ifelse(p_fdr < 0.05, "P≤0.05", "P>0.05")
  )

spec_arch_full <- bind_cols(morph_arch_df, scar_arch_df)

p_mantel_net <- qcorrplot(
  correlate(spec_arch_full, method = "spearman"),
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
    name   = "Mantel test\n(FDR corrected)"
  ) +
  scale_size_continuous(range = c(0.5, 2.5), name = "Mantel's |r|") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "right", plot.margin = margin(20, 20, 20, 20))

ggsave(here("analysis/output/figures/L1_Mantel_Network.png"),
       plot = p_mantel_net, width = 10, height = 8, dpi = 300, bg = "white")
cat("图已保存：L1_Mantel_Network.png\n")


# ------------------------------------------------------------------------------
# L1-2：CIA + RV 置换检验
# ------------------------------------------------------------------------------

cat("\n==== L1-2：CIA ====\n")

morph_arch_ilr <- as.data.frame(ilr(replace_zeros(as.matrix(morph_arch))))
scar_arch_ilr  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_arch))))
rownames(morph_arch_ilr) <- arch_ids
rownames(scar_arch_ilr)  <- arch_ids
colnames(morph_arch_ilr) <- paste0("M_ilr", seq_len(ncol(morph_arch_ilr)))
colnames(scar_arch_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_arch_ilr)))

dudi_morph <- dudi.pca(morph_arch_ilr, center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(morph_arch_ilr))
dudi_scar  <- dudi.pca(scar_arch_ilr,  center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(scar_arch_ilr))

coin_arch   <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
cia_inertia <- coin_arch$eig / sum(coin_arch$eig) * 100
cat("RV 系数：", round(coin_arch$RV, 4), "\n")

set.seed(42)
rv_test <- randtest(coin_arch, nrepet = 9999)
cat("\nRV 置换检验：\n"); print(rv_test)

# CIA 坐标 + 箭头长度/方位
scores_morph <- as.data.frame(coin_arch$lX) %>% rownames_to_column("ID")
scores_scar  <- as.data.frame(coin_arch$lY) %>% rownames_to_column("ID")

scores_combined <- left_join(
  scores_morph %>% select(ID, Axis1_M = AxcX1, Axis2_M = AxcX2),
  scores_scar  %>% select(ID, Axis1_S = AxcY1, Axis2_S = AxcY2),
  by = "ID"
) %>%
  mutate(
    arrow_length = sqrt((Axis1_M - Axis1_S)^2 + (Axis2_M - Axis2_S)^2),
    arrow_angle  = atan2(Axis2_S - Axis2_M, Axis1_S - Axis1_M)
  ) %>%
  left_join(meta_arch %>% select(ID, layer, raw_material, core_type), by = "ID")

# CIA 双标图
p_cia_biplot <- ggplot() +
  geom_segment(
    data = scores_combined,
    aes(x = Axis1_M, y = Axis2_M, xend = Axis1_S, yend = Axis2_S,
        color = arrow_length),
    linewidth = 0.7, alpha = 0.75,
    arrow = arrow(length = unit(0.10, "cm"), type = "closed")
  ) +
  geom_point(data = scores_combined,
             aes(x = Axis1_M, y = Axis2_M, shape = layer),
             size = 2.8, fill = "white", color = "grey25",
             stroke = 1.1, alpha = 0.9) +
  geom_point(data = scores_combined,
             aes(x = Axis1_S, y = Axis2_S, shape = layer),
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
    values = c("Layer 2" = 21, "Layer 3" = 22,
               "Layer 4" = 24, "Other"   = 20),
    name = "Layer\n(● morph / △ scar)"
  ) +
  theme_bw(base_size = 10) +
  labs(
    title    = sprintf("CIA Biplot  |  RV = %.3f, p = %.3f",
                       coin_arch$RV, rv_test$pvalue),
    subtitle = "Arrow: morphology → scar-direction coordinate; length = decoupling",
    x = sprintf("CIA Axis 1 (%.1f%%)", cia_inertia[1]),
    y = sprintf("CIA Axis 2 (%.1f%%)", cia_inertia[2])
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
    legend.position = "right"
  )

ggsave(here("analysis/output/figures/L1_CIA_Biplot.png"),
       plot = p_cia_biplot, width = 10, height = 8, dpi = 300, bg = "white")
cat("图已保存：L1_CIA_Biplot.png\n")

l1_results <- tibble(
  method  = c("Mantel (cosine, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_arch$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(arch_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/L1_results.csv"))
cat("\n第一层结论：\n"); print(l1_results)


# ==============================================================================
# ========== 第二层：联合证据 ==========
# ==============================================================================
# 三条证据链：
#   A. 分组 Mantel r（形态–技术相关性在各子组内是否显著）
#   B. CIA 箭头长度分组差异（各子组解耦量是否不同）
#   C. CIA 箭头方位分组差异（各子组解耦方向是否集中/一致）
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  第二层：联合证据                                          ##\n")
cat("################################################################\n")


# ==============================================================================
# ---- 第二层 A：分组 Mantel ----
# ==============================================================================

cat("\n---------- L2-A：分组 Mantel ----------\n")

mantel_within_group <- function(group_val, group_col, D_morph_full,
                                D_scar_full, meta_df, n_perm = 9999) {
  ids <- meta_df %>%
    filter(.data[[group_col]] == group_val) %>%
    pull(ID)
  if (length(ids) < 5) {
    cat(sprintf("  [跳过] %s = %s：样本量 %d < 5\n",
                group_col, group_val, length(ids)))
    return(NULL)
  }
  res <- mantel(extract_subdist(D_morph_full, ids),
                extract_subdist(D_scar_full, ids),
                method = "spearman", permutations = n_perm)
  tibble(group_var = group_col, group = group_val,
         n = length(ids), mantel_r = res$statistic, p_raw = res$signif)
}

run_grouped_mantel <- function(meta_obj, group_col, label) {
  if (is.null(meta_obj)) {
    cat(sprintf("  [跳过] %s 分组 Mantel：分组不足\n", label))
    return(NULL)
  }
  map_dfr(unique(meta_obj[[group_col]]),
          ~ mantel_within_group(.x, group_col,
                                D_morph_arch, D_scar_arch, meta_arch)) %>%
    mutate(
      p_fdr        = p.adjust(p_raw, method = "fdr"),
      significance = case_when(
        p_fdr < 0.001 ~ "***",
        p_fdr < 0.01  ~ "**",
        p_fdr < 0.05  ~ "*",
        p_fdr < 0.10  ~ ".",
        TRUE          ~ "ns"
      )
    )
}

mantel_by_layer    <- run_grouped_mantel(meta_layer,    "layer",        "层位")
mantel_by_rawmat   <- run_grouped_mantel(meta_rawmat,   "raw_material", "原料")
mantel_by_coretype <- run_grouped_mantel(meta_coretype, "core_type",    "石核类型")

cat("\n层位分组 Mantel：\n")
if (!is.null(mantel_by_layer))
  print(mantel_by_layer %>% mutate(across(c(mantel_r, p_raw, p_fdr), ~ round(.x, 4))))

cat("\n原料分组 Mantel：\n")
if (!is.null(mantel_by_rawmat))
  print(mantel_by_rawmat %>% mutate(across(c(mantel_r, p_raw, p_fdr), ~ round(.x, 4))))

cat("\n石核类型分组 Mantel：\n")
if (!is.null(mantel_by_coretype))
  print(mantel_by_coretype %>% mutate(across(c(mantel_r, p_raw, p_fdr), ~ round(.x, 4))))

# 合并 + 可视化
l2_mantel <- bind_rows(
  mantel_by_layer, mantel_by_rawmat, mantel_by_coretype
) %>%
  mutate(
    group_var_label = recode(group_var,
                             "layer"        = "Layer",
                             "raw_material" = "Raw Material",
                             "core_type"    = "Core Type"),
    group = as.factor(group)
  ) %>%
  group_by(group_var_label) %>%
  mutate(group = fct_reorder(group, mantel_r)) %>%
  ungroup()

if (nrow(l2_mantel) > 0) {
  
  # 哑铃图
  p_l2_dumbbell <- ggplot(l2_mantel,
                          aes(x = mantel_r, y = group,
                              color = significance,
                              shape = group_var_label)) +
    geom_vline(xintercept = mantel_global$statistic,
               linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "grey70", linewidth = 0.4) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = sprintf("n=%d", n)),
              hjust = -0.35, size = 2.5, color = "grey40") +
    facet_wrap(~ group_var_label, scales = "free_y", ncol = 1) +
    scale_color_manual(
      values = c("***" = "#C0392B", "**" = "#E67E22", "*" = "#F1C40F",
                 "."  = "#27AE60", "ns" = "#95A5A6"),
      name = "Significance\n(FDR corrected)"
    ) +
    scale_shape_manual(
      values = c("Layer" = 16, "Raw Material" = 17, "Core Type" = 15),
      name   = "Grouping variable"
    ) +
    theme_bw(base_size = 10) +
    labs(
      title    = "L2-A: Within-group Mantel Test (Morphology × Scar Direction)",
      subtitle = "Each point = Mantel r within one subgroup; dashed line = global baseline",
      x = "Mantel r (Spearman, FDR corrected)", y = NULL
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      strip.text    = element_text(face = "bold", size = 9),
      legend.position = "right"
    )
  
  ggsave(here("analysis/output/figures/L2_Mantel_Grouped_dumbbell.png"),
         plot = p_l2_dumbbell,
         width = 10, height = max(4, nrow(l2_mantel) * 0.7 + 2),
         dpi = 300, bg = "white")
  cat("图已保存：L2_Mantel_Grouped_dumbbell.png\n")
  
  # 热图
  heatmap_df <- l2_mantel %>%
    mutate(group = fct_reorder(group, mantel_r, .fun = mean))
  
  y_labels <- heatmap_df %>%
    mutate(y_int = as.integer(group)) %>%
    distinct(y_int, group) %>%
    arrange(y_int) %>%
    mutate(group_chr = as.character(group))
  
  heatmap_df <- heatmap_df %>%
    mutate(y_int = as.integer(group))
  
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
      title    = "L2-A: Within-group Mantel r (Morphology × Scar Direction)",
      subtitle = sprintf("Global baseline: r = %.3f (ns)", mantel_global$statistic),
      x = "Grouping variable", y = "Subgroup"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      axis.text.y   = element_text(size = 8)
    )
  
  n_groups <- length(unique(heatmap_df$group))
  n_vars   <- length(unique(heatmap_df$group_var_label))
  
  ggsave(here("analysis/output/figures/L2_Mantel_Grouped_heatmap.png"),
         plot   = p_l2_heatmap,
         width  = max(6, n_vars * 3 + 2),
         height = max(4, n_groups * 0.6 + 2),
         dpi = 300, bg = "white")
  cat("图已保存：L2_Mantel_Grouped_heatmap.png\n")
  
  write_csv(l2_mantel, here("analysis/data/derived_data/L2_grouped_mantel.csv"))
}


# ==============================================================================
# ---- 第二层 B：CIA 箭头长度分组差异 ----
# ==============================================================================

cat("\n---------- L2-B：箭头长度分组差异 ----------\n")

# Kruskal-Wallis + 两两 Wilcoxon + 箱线图
run_arrow_length_analysis <- function(group_col, group_label, palette) {
  
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
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
                             p.adjust.method = "fdr", exact = FALSE)
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
      title    = sprintf("L2-B: CIA Arrow Length by %s", group_label),
      subtitle = "Longer arrows = greater morphology–technique decoupling",
      x = NULL, y = "Arrow length (CIA space)"
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      axis.text.x   = element_text(angle = 15, hjust = 1)
    )
  
  fname <- sprintf("analysis/output/figures/L2_Arrow_Length_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  
  list(sub_df = sub_df, kw = kw, pw = pw)
}

res_len_layer    <- run_arrow_length_analysis("layer",        "Layer",       layer_pal)
res_len_rawmat   <- run_arrow_length_analysis("raw_material", "Raw_Material",rawmat_pal)
res_len_coretype <- run_arrow_length_analysis("core_type",    "Core_Type",   coretype_pal)

# 箭头长度描述统计汇总
cat("\n==== 箭头长度描述统计 ====\n")
for (gc in c("layer", "raw_material", "core_type")) {
  cat(sprintf("\n--- %s ---\n", gc))
  scores_combined %>%
    filter(!is.na(.data[[gc]]), .data[[gc]] != "Other") %>%
    group_by(.data[[gc]]) %>%
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
}


# ==============================================================================
# ---- 第二层 C：CIA 箭头方位圆形统计 ----
# ==============================================================================

cat("\n---------- L2-C：箭头方位圆形统计 ----------\n")

run_circular_analysis <- function(group_col, group_label, palette) {
  
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
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
  
  # --- 描述统计 ---
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
  
  # --- Rayleigh 检验 ---
  cat(sprintf("\n----- %s Rayleigh 检验 -----\n", group_label))
  rayleigh_res <- map_dfr(valid_groups, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    circ_obj <- circular(angles, type = "angles",
                         units = "radians", modulo = "2pi")
    rt <- rayleigh.test(circ_obj)
    cat(sprintf("  %s: statistic = %.4f, p = %.4f → %s\n",
                g, rt$statistic, rt$p.value,
                ifelse(rt$p.value < 0.05, "方位集中（非随机）", "方位分散（随机）")))
    tibble(
      group_var   = group_col,
      group       = g,
      rayleigh_U  = round(rt$statistic, 4),
      rayleigh_p  = round(rt$p.value,   4),
      conclusion  = ifelse(rt$p.value < 0.05, "concentrated", "uniform")
    )
  })
  
  # --- Watson 两样本检验（置换 p 值）---
  watson_perm_test <- function(x1, x2, B = 9999) {
    # 观测值
    a1 <- circular(x1, type = "angles", units = "radians", modulo = "2pi")
    a2 <- circular(x2, type = "angles", units = "radians", modulo = "2pi")
    wt_obs <- watson.two.test(a1, a2)
    obs_u2 <- as.numeric(wt_obs$statistic)
    
    # 置换分布
    x_all <- c(x1, x2)
    n1 <- length(x1)
    n_all <- length(x_all)
    
    perm_u2 <- replicate(B, {
      idx <- sample.int(n_all)
      p1 <- x_all[idx[1:n1]]
      p2 <- x_all[idx[(n1 + 1):n_all]]
      
      wt_perm <- watson.two.test(
        circular(p1, type = "angles", units = "radians", modulo = "2pi"),
        circular(p2, type = "angles", units = "radians", modulo = "2pi")
      )
      as.numeric(wt_perm$statistic)
    })
    
    p_perm <- (sum(perm_u2 >= obs_u2) + 1) / (B + 1)
    
    list(statistic = obs_u2, p.value = p_perm)
  }
  
  if (length(valid_groups) >= 2) {
    cat(sprintf("\n----- %s Watson 两样本检验（置换 p 值） -----\n", group_label))
    pairs <- combn(valid_groups, 2, simplify = FALSE)
    
    watson_res <- map_dfr(pairs, function(pair) {
      x1 <- sub_df %>%
        filter(.data[[group_col]] == pair[1]) %>%
        pull(arrow_angle)
      
      x2 <- sub_df %>%
        filter(.data[[group_col]] == pair[2]) %>%
        pull(arrow_angle)
      
      wt <- watson_perm_test(x1, x2, B = 9999)
      
      cat(sprintf("  %s vs %s: U² = %.4f, p = %.4f → %s\n",
                  pair[1], pair[2], wt$statistic, wt$p.value,
                  ifelse(wt$p.value < 0.05, "分布不同", "差异不显著")))
      
      tibble(
        group_var   = group_col,
        group1      = pair[1],
        group2      = pair[2],
        U2_statistic = round(wt$statistic, 4),
        p_value      = round(wt$p.value, 4),
        conclusion   = ifelse(wt$p.value < 0.05, "different", "ns")
      )
    })
  } else {
    watson_res <- NULL
  }
  
  # --- 玫瑰图 ---
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
  
  p_rose <- ggplot(rose_df, aes(x = angle_deg,
                                fill = .data[[group_col]])) +
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
      labels = c("0°\n(+CIA1)", "45°", "90°\n(+CIA2)",
                 "135°", "180°\n(-CIA1)", "225°", "270°\n(-CIA2)", "315°")
    ) +
    scale_fill_manual(values  = palette, name = group_label,
                      aesthetics = c("fill", "color")) +
    facet_wrap(reformulate(group_col),
               ncol = min(length(valid_groups), 3)) +
    theme_bw(base_size = 10) +
    labs(
      title    = sprintf("L2-C: CIA Arrow Direction by %s", group_label),
      subtitle = "Dashed line = mean direction; direction = axis of decoupling in CIA space",
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
  
  n_g <- length(valid_groups)
  fname <- sprintf("analysis/output/figures/L2_Arrow_Direction_rose_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p_rose,
         width = min(4 + n_g * 3, 16), height = 5, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  
  list(desc = circ_desc, rayleigh = rayleigh_res, watson = watson_res)
}

res_circ_layer    <- run_circular_analysis("layer",        "Layer",       layer_pal)
res_circ_rawmat   <- run_circular_analysis("raw_material", "Raw_Material",rawmat_pal)
res_circ_coretype <- run_circular_analysis("core_type",    "Core_Type",   coretype_pal)

# 保存圆形统计汇总
circ_desc_all <- bind_rows(
  res_circ_layer$desc,
  res_circ_rawmat$desc,
  res_circ_coretype$desc
)
rayleigh_all <- bind_rows(
  res_circ_layer$rayleigh,
  res_circ_rawmat$rayleigh,
  res_circ_coretype$rayleigh
)

if (nrow(circ_desc_all) > 0 && nrow(rayleigh_all) > 0) {
  left_join(circ_desc_all, rayleigh_all, by = c("group_var", "group")) %>%
    write_csv(here("analysis/data/derived_data/L2_circular_stats.csv"))
  cat("已保存：L2_circular_stats.csv\n")
}

# 保存箭头长度汇总
scores_combined %>%
  select(ID, layer, raw_material, core_type,
         arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/L2_arrow_stats.csv"))
cat("已保存：L2_arrow_stats.csv\n")


# ==============================================================================
# ========== 第三层：PERMANOVA ==========
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  第三层：PERMANOVA — 分组结构分析                          ##\n")
cat("################################################################\n")

run_permanova <- function(dist_mat, group_vec, group_name, domain_name,
                          n_perm = 9999) {
  group_vec <- as.character(group_vec)
  if (length(unique(group_vec)) < 2) {
    cat(sprintf("  [跳过] %s ~ %s：分组水平不足\n", domain_name, group_name))
    return(NULL)
  }
  df_tmp <- data.frame(group = factor(group_vec))
  res    <- adonis2(dist_mat ~ group, data = df_tmp,
                    permutations = n_perm, add = "lingoes")
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

run_permdisp <- function(dist_mat, group_vec, group_name, domain_name) {
  group_vec <- as.character(group_vec)
  if (length(unique(group_vec)) < 2) return(NULL)
  eig_vals    <- cmdscale(dist_mat, eig = TRUE)$eig
  correction  <- abs(min(c(eig_vals[eig_vals < 0], 0)))
  D_corrected <- as.dist(as.matrix(dist_mat) + correction)
  bd  <- betadisper(D_corrected, factor(group_vec))
  res <- permutest(bd, permutations = 9999)
  cat(sprintf("  PERMDISP — %s ~ %s: F = %.3f, p = %.3f\n",
              domain_name, group_name,
              res$tab$F[1], res$tab$`Pr(>F)`[1]))
  tibble(domain = domain_name, grouping = group_name,
         F_disp = res$tab$F[1], p_disp = res$tab$`Pr(>F)`[1])
}

pairwise_permanova <- function(dist_mat, group_vec, group_name,
                               domain_name, n_perm = 9999) {
  group_vec <- as.character(group_vec)
  groups    <- unique(group_vec)
  if (length(groups) < 2) return(NULL)
  map_dfr(combn(groups, 2, simplify = FALSE), function(pair) {
    idx    <- group_vec %in% pair
    d_sub  <- as.dist(as.matrix(dist_mat)[idx, idx])
    df_tmp <- data.frame(group = factor(group_vec[idx]))
    res    <- adonis2(d_sub ~ group, data = df_tmp,
                      permutations = n_perm, add = "lingoes")
    tibble(domain = domain_name, grouping = group_name,
           group1 = pair[1], group2 = pair[2],
           R2 = res$R2[1], F_value = res$F[1],
           p_raw = res$`Pr(>F)`[1])
  })
}

# 构建子距离矩阵
if (!is.null(meta_layer)) {
  D_morph_layer <- extract_subdist(D_morph_arch, meta_layer$ID)
  D_scar_layer  <- extract_subdist(D_scar_arch,  meta_layer$ID)
}
if (!is.null(meta_rawmat)) {
  D_morph_rawmat <- extract_subdist(D_morph_arch, meta_rawmat$ID)
  D_scar_rawmat  <- extract_subdist(D_scar_arch,  meta_rawmat$ID)
}
if (!is.null(meta_coretype)) {
  D_morph_coretype <- extract_subdist(D_morph_arch, meta_coretype$ID)
  D_scar_coretype  <- extract_subdist(D_scar_arch,  meta_coretype$ID)
}

cat("\n==== L3-1：全局 PERMANOVA ====\n")

permanova_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    run_permanova(D_morph_layer,    meta_layer$layer,          "Layer",        "Morphology"),
    run_permanova(D_scar_layer,     meta_layer$layer,          "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    run_permanova(D_morph_rawmat,   meta_rawmat$raw_material,  "Raw Material", "Morphology"),
    run_permanova(D_scar_rawmat,    meta_rawmat$raw_material,  "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    run_permanova(D_morph_coretype, meta_coretype$core_type,   "Core Type",    "Morphology"),
    run_permanova(D_scar_coretype,  meta_coretype$core_type,   "Core Type",    "Scar Direction")
  )
)

cat("\n==== L3-2：PERMDISP ====\n")

permdisp_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    run_permdisp(D_morph_layer,    meta_layer$layer,          "Layer",        "Morphology"),
    run_permdisp(D_scar_layer,     meta_layer$layer,          "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    run_permdisp(D_morph_rawmat,   meta_rawmat$raw_material,  "Raw Material", "Morphology"),
    run_permdisp(D_scar_rawmat,    meta_rawmat$raw_material,  "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    run_permdisp(D_morph_coretype, meta_coretype$core_type,   "Core Type",    "Morphology"),
    run_permdisp(D_scar_coretype,  meta_coretype$core_type,   "Core Type",    "Scar Direction")
  )
)

cat("\n==== L3-3：两两 PERMANOVA ====\n")

pairwise_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    pairwise_permanova(D_morph_layer,    meta_layer$layer,          "Layer",        "Morphology"),
    pairwise_permanova(D_scar_layer,     meta_layer$layer,          "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    pairwise_permanova(D_morph_rawmat,   meta_rawmat$raw_material,  "Raw Material", "Morphology"),
    pairwise_permanova(D_scar_rawmat,    meta_rawmat$raw_material,  "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    pairwise_permanova(D_morph_coretype, meta_coretype$core_type,   "Core Type",    "Morphology"),
    pairwise_permanova(D_scar_coretype,  meta_coretype$core_type,   "Core Type",    "Scar Direction")
  )
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

# R² 比较图
p_r2_compare <- permanova_results %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    domain   = factor(domain,   levels = c("Morphology", "Scar Direction")),
    grouping = factor(grouping, levels = c("Layer", "Raw Material", "Core Type"))
  ) %>%
  ggplot(aes(x = grouping, y = R2, fill = domain)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sig_label, y = R2 + 0.003),
            position = position_dodge(width = 0.7),
            vjust = 0, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(
    values = c("Morphology"     = "#A1C2E6",
               "Scar Direction" = "#FFBAE0"),
    name   = "Domain"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.12))) +
  theme_bw(base_size = 10) +
  labs(
    title    = "L3: PERMANOVA — Proportion of Variance Explained (R²)",
    subtitle = "*** p<0.001  ** p<0.01  * p<0.05  ns = not significant",
    x = "Grouping factor", y = "R² (proportion of variance)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
    legend.position = "right"
  )

ggsave(here("analysis/output/figures/L3_PERMANOVA_R2.png"),
       plot = p_r2_compare, width = 8, height = 5, dpi = 300, bg = "white")
cat("图已保存：L3_PERMANOVA_R2.png\n")

# CAP 图
plot_cap <- function(dist_mat, meta_df, group_col, title_str, palette_vec) {
  df_group   <- data.frame(group = meta_df[[group_col]])
  cap_res    <- capscale(dist_mat ~ group, data = df_group)
  all_scores <- scores(cap_res, display = "sites", choices = 1:2)
  cap_scores <- as.data.frame(all_scores) %>%
    mutate(ID = meta_df$ID, group = meta_df[[group_col]])
  x_col <- colnames(cap_scores)[1]
  y_col <- colnames(cap_scores)[2]
  eig_all <- eigenvals(cap_res)
  eig_pos <- eig_all[eig_all > 0]
  eig_df  <- tibble(name = names(eig_all), value = as.numeric(eig_all))
  pct_df  <- eig_df %>%
    filter(name %in% c(x_col, y_col)) %>%
    mutate(pct = round(value / sum(eig_pos) * 100, 1))
  x_lab <- sprintf("%s (%.1f%%)", x_col, pct_df$pct[pct_df$name == x_col])
  y_lab <- sprintf("%s (%.1f%%)", y_col, pct_df$pct[pct_df$name == y_col])
  ggplot(cap_scores,
         aes(x = .data[[x_col]], y = .data[[y_col]], color = group)) +
    stat_ellipse(aes(fill = group), geom = "polygon",
                 alpha = 0.08, level = 0.9,
                 linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.35) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_text_repel(aes(label = ID), size = 2.0, color = "grey40",
                    max.overlaps = 6, show.legend = FALSE) +
    scale_color_manual(values = palette_vec, name = group_col,
                       aesthetics = c("color", "fill")) +
    theme_bw(base_size = 10) +
    labs(title = title_str, x = x_lab, y = y_lab) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "right")
}

save_cap_pair <- function(meta_obj, D_m, D_s, group_col, palette, tag) {
  if (is.null(meta_obj)) {
    cat(sprintf("  [跳过] %s CAP 图：有效分组不足\n", tag))
    return(invisible(NULL))
  }
  get_stat <- function(dom) {
    permanova_results %>%
      filter(domain == dom, grouping == tag) %>%
      summarise(lab = sprintf("R²=%.3f, p=%.3f", R2, p_value)) %>%
      pull(lab)
  }
  p1 <- plot_cap(D_m, meta_obj, group_col,
                 sprintf("Morphology — CAP by %s", tag), palette)
  p2 <- plot_cap(D_s, meta_obj, group_col,
                 sprintf("Scar Direction — CAP by %s", tag), palette)
  p_combined <- (p1 | p2) +
    plot_annotation(
      title = sprintf(
        "CAP — %s  |  Morphology: %s  |  Scar Direction: %s", tag,
        ifelse(nrow(filter(permanova_results,
                           domain == "Morphology", grouping == tag)) > 0,
               get_stat("Morphology"), "—"),
        ifelse(nrow(filter(permanova_results,
                           domain == "Scar Direction", grouping == tag)) > 0,
               get_stat("Scar Direction"), "—")
      ),
      theme = theme(plot.title = element_text(face = "bold",
                                              size = 10, hjust = 0.5))
    )
  fname <- sprintf("analysis/output/figures/L3_PERMANOVA_CAP_%s.png",
                   tolower(str_replace_all(tag, " ", "_")))
  ggsave(here(fname), plot = p_combined,
         width = 14, height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
}

save_cap_pair(meta_layer,    D_morph_layer,    D_scar_layer,
              "layer",        layer_pal,    "Layer")
save_cap_pair(meta_rawmat,   D_morph_rawmat,   D_scar_rawmat,
              "raw_material", rawmat_pal,   "Raw Material")
save_cap_pair(meta_coretype, D_morph_coretype, D_scar_coretype,
              "core_type",   coretype_pal, "Core Type")


# ==============================================================================
# ---- 三层结果汇总打印 ----
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("##  三层分析结果汇总                                          ##\n")
cat("################################################################\n")

cat("\n【第一层：基线】\n")
cat(sprintf("  Mantel r = %.4f, p = %.3f  →  %s\n",
            mantel_global$statistic, mantel_global$signif,
            ifelse(mantel_global$signif < 0.05, "显著相关", "独立（ns）")))
cat(sprintf("  RV       = %.4f, p = %.3f  →  %s\n",
            coin_arch$RV, rv_test$pvalue,
            ifelse(rv_test$pvalue < 0.05, "显著协变", "独立（ns）")))

cat("\n【第二层 A：分组 Mantel】\n")
if (nrow(l2_mantel) > 0) {
  n_sig <- sum(l2_mantel$p_fdr < 0.05, na.rm = TRUE)
  cat(sprintf("  共检验 %d 个子组，FDR 校正后显著：%d 个\n",
              nrow(l2_mantel), n_sig))
  if (n_sig == 0) {
    cat("  → 所有子组均不显著，形态–技术独立性稳健\n")
  } else {
    cat("  → 以下子组显著（潜在新发现）：\n")
    l2_mantel %>%
      filter(p_fdr < 0.05) %>%
      select(group_var, group, n, mantel_r, p_fdr, significance) %>%
      print()
  }
}

cat("\n【第二层 B：箭头长度】参见 L2_Arrow_Length_*.png 及描述统计\n")
cat("【第二层 C：箭头方位】参见 L2_Arrow_Direction_rose_*.png 及 L2_circular_stats.csv\n")

cat("\n【第三层：PERMANOVA — 各域分组结构】\n")
permanova_results %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(domain, grouping, R2, F_value, p_value, sig) %>%
  mutate(across(c(R2, F_value, p_value), ~ round(.x, 4))) %>%
  arrange(domain, grouping) %>%
  print()


# ==============================================================================
# ---- 保存数值结果 ----
# ==============================================================================

write_csv(l1_results,
          here("analysis/data/derived_data/L1_results.csv"))

if (exists("l2_mantel") && nrow(l2_mantel) > 0)
  write_csv(l2_mantel,
            here("analysis/data/derived_data/L2_grouped_mantel.csv"))

bind_rows(
  permanova_results %>%
    mutate(test = "PERMANOVA", statistic = F_value) %>%
    select(test, domain, grouping, R2, statistic, p_value),
  permdisp_results %>%
    mutate(test = "PERMDISP", R2 = NA,
           statistic = F_disp, p_value = p_disp) %>%
    select(test, domain, grouping, R2, statistic, p_value),
  pairwise_results %>%
    mutate(test = "Pairwise_PERMANOVA",
           statistic = F_value, p_value = p_fdr) %>%
    select(test, domain, grouping, R2, statistic, p_value)
) %>%
  write_csv(here("analysis/data/derived_data/L3_permanova.csv"))

scores_combined %>%
  write_csv(here("analysis/data/derived_data/CIA_scores_full.csv"))

cat("\n\n========== 全部分析完成 ==========\n")
cat("主要输出文件：\n")
cat("  L1_results.csv\n")
cat("  L2_grouped_mantel.csv\n")
cat("  L2_arrow_stats.csv\n")
cat("  L2_circular_stats.csv\n")
cat("  L3_permanova.csv\n")
cat("  CIA_scores_full.csv\n")