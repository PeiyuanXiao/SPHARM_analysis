# ==============================================================================
# SDG_cores_statistics.R
#
# 第一层：整体 Mantel + CoIA
#   形态与片疤方向是否整体独立？
#   附：CoIA 箭头长度 + 箭头方位描述
#   附：CoIA 桑基图（ILR -> PCA -> CoIA 贡献流）
#
# 第二层：分组 Mantel（层位 × 原料 × 类型）
#   这种独立性是否在所有条件下都成立？
#   联合证据：
#     - 分组 Mantel r
#     - 箭头长度分组差异
#     - 箭头方位分组差异
#     - spectral_entropy × Layer / Raw_mat / Core_type
#
# 第三层：PERMANOVA（层位 × 原料 × 类型）
#   形态和片疤方向各自的分组结构
#
# 输入：
#   - analysis/data/derived_data/SPHARM_direction_filter.rds
#   - analysis/data/derived_data/SPHARM_morphology_filter.rds
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#
# 输出（figures/）：
#   L1_Mantel_Network.png
#   L1_CoIA_Biplot_layer.png
#   L1_CoIA_Biplot_raw_material.png
#   L1_CoIA_Biplot_core_type.png
#   L1_CoIA_Diagnostics.png
#   L1_CoIA_Sankey.png
#   L2_Arrow_Length_layer.png
#   L2_Arrow_Length_raw_material.png
#   L2_Arrow_Length_core_type.png
#   L2_Arrow_Direction_rose_layer.png
#   L2_Arrow_Direction_rose_raw_material.png
#   L2_Arrow_Direction_rose_core_type.png
#   L2D_SE_Direction_Layer_boxplot.png
#   L2D_SE_Morphology_Layer_boxplot.png
#   L2D_SE_Direction_RawMat_boxplot.png
#   L2D_SE_Morphology_RawMat_boxplot.png
#   L2D_SE_Direction_CoreType_boxplot.png
#   L2D_SE_Morphology_CoreType_boxplot.png
#   L_CoIA_composite.png
#
# 输出（derived_data/）：
#   L1_results.csv
#   L2_grouped_mantel.csv
#   L2_arrow_stats.csv
#   L2_circular_stats.csv
#   L2D_SE_desc_stats.csv
#   L2D_SE_dunn_results.csv
#   L3_permanova.csv
#   CoIA_scores_full.csv
#   CoIA_coords_full.csv
#   PCA_CoIA_contribution.csv
#   SDG_morph_ILR_scores.csv
#   SDG_scar_ILR_scores.csv
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
library(FSA)


# ==============================================================================
# ---- 全局常量：色盘与顺序 ----
# ==============================================================================

LAYER_ORDER <- c("Layer 2", "Layer 3", "Layer 4")

CORETYPE_ORDER <- c(
  "Unifacial_unidirection",
  "Unifacial_centripetal",
  "Bifacial_adjacent",
  "Bifacial_independent",
  "Bifacial_centripetal",
  "Multifacial",
  "Core_on_flake"
)

# 剔除的石核类型
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")

layer_pal <- c(
  "Layer 2" = "#802520",
  "Layer 3" = "#5C7F71",
  "Layer 4" = "#BA8530"
)

rawmat_pal <- c(
  "chert"     = "#4A6E8A",
  "sandstone" = "#802520"
)

coretype_pal <- c(
  "Unifacial_unidirection" = "#5C7F71",
  "Unifacial_centripetal"  = "#BA8530",
  "Bifacial_adjacent"      = "#802520",
  "Bifacial_independent"   = "#B26538",
  "Bifacial_centripetal"   = "#788C4A",
  "Multifacial"            = "#8A7A68",
  "Core_on_flake"          = "#4A6E8A"
)


# ==============================================================================
# ---- 全局辅助函数 ----
# ==============================================================================

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
      "  [跳过] %s：有效组数不足（需 >= 2 组每组 >= %d 件）。当前：%s\n",
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
# ---- 数据准备（所有层共用）----
# ==============================================================================

POWER_COLS_DIR   <- paste0("power_l", 1:6)   # 方向谱 l=1-6
POWER_COLS_MORPH <- paste0("power_l", 1:8)   # 形态谱 l=1-8

SPHARM_direction_filter  <- readRDS(here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
SPHARM_morphology_filter <- readRDS(here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

df_morph_raw <- SPHARM_morphology_filter
df_scar_raw  <- SPHARM_direction_filter

common_ids <- intersect(df_morph_raw$ID, df_scar_raw$ID)
df_morph_raw <- df_morph_raw %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar_raw  <- df_scar_raw  %>% filter(ID %in% common_ids) %>% arrange(ID)

cat("==== 数据对齐 ====\n")
cat("共有标本：", length(common_ids),
    "；ID 完全匹配：", all(df_morph_raw$ID == df_scar_raw$ID), "\n\n")

# 外部元数据
core_meta_raw <- read_excel(
  here("analysis/data/raw_data/SDG_core_metric.xlsx")
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

# 提取 power 特征矩阵（方向谱和形态谱使用各自的列数）
morph_power <- df_morph_raw %>%
  select(all_of(POWER_COLS_MORPH)) %>%
  rename_with(~ paste0("M", seq_along(.))) %>%
  as.data.frame()
scar_power <- df_scar_raw %>%
  select(all_of(POWER_COLS_DIR)) %>%
  rename_with(~ paste0("S", seq_along(.))) %>%
  as.data.frame()
rownames(morph_power) <- df_morph_raw$ID
rownames(scar_power)  <- df_scar_raw$ID

morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

# ILR 变换 -> 欧氏距离
morph_ilr_all <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power_clean))))
scar_ilr_all  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power_clean))))
rownames(morph_ilr_all) <- rownames(morph_power_clean)
rownames(scar_ilr_all)  <- rownames(scar_power_clean)

D_morph_all <- dist(morph_ilr_all)
D_scar_all  <- dist(scar_ilr_all)

# 分离考古标本（去除 IM_ 参照件）
arch_ids <- rownames(morph_power_clean)[!str_starts(rownames(morph_power_clean), "IM_") &
                                          !str_starts(rownames(morph_power_clean), "EXP")]
cat("考古标本数量（不含 IM_）：", length(arch_ids), "\n")

morph_arch     <- morph_power_clean[arch_ids, ]
scar_arch      <- scar_power_clean[arch_ids, ]
morph_ilr_arch <- morph_ilr_all[arch_ids, ]
scar_ilr_arch  <- scar_ilr_all[arch_ids, ]

D_morph_arch <- extract_subdist(D_morph_all, arch_ids)
D_scar_arch  <- extract_subdist(D_scar_all,  arch_ids)

# 构建元数据
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

# 剔除 Handaxe 和 Pick
meta_arch <- meta_arch %>%
  filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))

# 更新 arch_ids 以反映剔除后的标本集
arch_ids <- meta_arch$ID
morph_arch     <- morph_power_clean[arch_ids, ]
scar_arch      <- scar_power_clean[arch_ids, ]
morph_ilr_arch <- morph_ilr_all[arch_ids, ]
scar_ilr_arch  <- scar_ilr_all[arch_ids, ]
D_morph_arch   <- extract_subdist(D_morph_all, arch_ids)
D_scar_arch    <- extract_subdist(D_scar_all,  arch_ids)

cat("剔除 Handaxe/Pick 后考古标本数量：", length(arch_ids), "\n")

# 合并 spectral_entropy
meta_arch <- meta_arch %>%
  left_join(
    df_scar_raw  %>% select(ID, SE_direction  = spectral_entropy),
    by = "ID"
  ) %>%
  left_join(
    df_morph_raw %>% select(ID, SE_morphology = spectral_entropy),
    by = "ID"
  )

cat("\n==== 元数据合并诊断 ====\n")
cat("层位：\n");     print(table(meta_arch$layer,        useNA = "ifany"))
cat("原料：\n");     print(table(meta_arch$raw_material,  useNA = "ifany"))
cat("石核类型：\n"); print(table(meta_arch$core_type,     useNA = "ifany"))
unmatched <- meta_arch %>% filter(is.na(raw_material)) %>% pull(ID)
if (length(unmatched) > 0) {
  cat("未匹配 ID（请检查格式）：\n"); print(unmatched)
}

meta_layer    <- safe_filter_groups(meta_arch %>% filter(layer != "Other"), "layer")
meta_rawmat   <- safe_filter_groups(meta_arch, "raw_material")
meta_coretype <- safe_filter_groups(meta_arch, "core_type")


# ==============================================================================
# ========== 第一层：整体 Mantel + CoIA ==========
# ==============================================================================

cat("\n\n")
cat("##  第一层：整体 Mantel + CoIA — 建立基线                     ##\n")


# ------------------------------------------------------------------------------
# L1-1：全局 Mantel + linkET 网络图
# ------------------------------------------------------------------------------

cat("\n==== L1-1：全局 Mantel ====\n")
mantel_global <- mantel(D_morph_arch, D_scar_arch,
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
    p_holm       = p.adjust(p, method = "holm"),
    significance = ifelse(p_holm < 0.05, "P\u22640.05", "P>0.05")
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
  scale_fill_gradient2(
    low      = "#802520",
    mid      = "#F5EDDC",
    high     = "#4A6E8A",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman's rho"
  ) +
  scale_color_manual(
    values = c("P\u22640.05" = "#E6A5A5", "P>0.05" = "#BABABA"),
    name   = "Mantel test\n(Holm corrected)"
  ) +
  scale_size_continuous(range = c(0.5, 2.5), name = "Mantel's |r|") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(20, 20, 20, 20))

ggsave(here("analysis/output/figures/L1_Mantel_Network.png"),
       plot = p_mantel_net, width = 10, height = 8, dpi = 300, bg = "white")
cat("图已保存：L1_Mantel_Network.png\n")


# ------------------------------------------------------------------------------
# L1-2：CoIA + RV 置换检验
# ------------------------------------------------------------------------------

cat("\n==== L1-2：CoIA ====\n")

morph_arch_ilr <- morph_ilr_arch
scar_arch_ilr  <- scar_ilr_arch
colnames(morph_arch_ilr) <- paste0("M_ilr", seq_len(ncol(morph_arch_ilr)))
colnames(scar_arch_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_arch_ilr)))

dudi_morph <- dudi.pca(morph_arch_ilr, center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(morph_arch_ilr))
dudi_scar  <- dudi.pca(scar_arch_ilr,  center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(scar_arch_ilr))

# PCA 详细报告
report_pca <- function(dudi_obj, label) {
  n_ax <- length(dudi_obj$eig)
  eig  <- dudi_obj$eig
  pct  <- eig / sum(eig) * 100
  cum  <- cumsum(pct)
  
  cat(sprintf("\n====== %s PCA 报告 ======\n", label))
  
  cat("\n-- 特征值与方差解释 --\n")
  eig_df <- tibble(
    Axis       = paste0("PC", seq_len(n_ax)),
    Eigenvalue = round(eig, 4),
    Pct_var    = round(pct, 2),
    Cumul_pct  = round(cum, 2)
  )
  print(as.data.frame(eig_df))
  
  cat("\n-- 变量载荷（c1：ILR 变量在各 PCA 轴上的载荷）--\n")
  load_df <- as.data.frame(dudi_obj$c1)
  colnames(load_df) <- paste0("PC", seq_len(ncol(load_df)))
  n_ilr <- nrow(load_df)
  load_df$ILR_meaning <- sapply(seq_len(n_ilr), function(k) {
    num_ids <- paste0("l", seq_len(k))
    den_id  <- paste0("l", k + 1)
    sprintf("log(geomean(%s) / %s)",
            paste(num_ids, collapse = "+"), den_id)
  })
  print(load_df)
  
  cat("\n-- 样本得分描述统计（前2轴）--\n")
  score_df <- as.data.frame(dudi_obj$li)[, 1:min(2, n_ax), drop = FALSE]
  colnames(score_df) <- paste0("PC", seq_len(ncol(score_df)))
  score_stats <- score_df %>%
    pivot_longer(everything(), names_to = "Axis", values_to = "Score") %>%
    group_by(Axis) %>%
    summarise(
      mean = round(mean(Score), 4), sd = round(sd(Score), 4),
      min  = round(min(Score),  4), max = round(max(Score), 4),
      .groups = "drop"
    )
  print(as.data.frame(score_stats))
  
  cat("\n-- 各主轴主导变量（|载荷| 最大）--\n")
  load_num <- as.data.frame(dudi_obj$c1)
  for (ax in seq_len(min(2, n_ax))) {
    col <- load_num[[ax]]
    idx <- which.max(abs(col))
    cat(sprintf("  PC%d (%.1f%% var)：主导 ILR%d（载荷 %+.4f）-> %s\n",
                ax, pct[ax], idx, col[idx], load_df$ILR_meaning[idx]))
  }
  invisible(list(eig_df = eig_df, load_df = load_df))
}

pca_report_morph <- report_pca(dudi_morph, "形态谱")
pca_report_scar  <- report_pca(dudi_scar,  "方向谱")

coin_arch   <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
cia_inertia <- coin_arch$eig / sum(coin_arch$eig) * 100

cat("RV 系数：", round(coin_arch$RV, 4), "\n")

set.seed(42)
rv_test <- randtest(coin_arch, nrepet = 9999)
cat("\nRV 置换检验：\n"); print(rv_test)

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
  left_join(meta_arch %>% select(ID, layer, raw_material, core_type),
            by = "ID") %>%
  filter(!str_starts(ID, "EXP"))

# CoIA 坐标输出
cat("\n==== CoIA 样本坐标 ====\n")
cia_coords <- scores_combined %>%
  select(ID, layer, raw_material, core_type,
         Morph_Axis1 = Axis1_M, Morph_Axis2 = Axis2_M,
         Scar_Axis1  = Axis1_S, Scar_Axis2  = Axis2_S,
         arrow_length, arrow_angle)
print(cia_coords %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        as.data.frame())
write_csv(cia_coords,
          here("analysis/data/derived_data/CoIA_coords_full.csv"))
cat("已保存：CoIA_coords_full.csv\n")

# PCA 轴对 CoIA 轴贡献
cat("\n==== 各端 PCA 轴对 CoIA 轴的贡献（weight^2）====\n")

compute_pca_cia_contribution <- function(a_mat, pct_vec, endpoint_label) {
  df <- as.data.frame(a_mat)
  colnames(df) <- paste0("CoIA_Ax", seq_len(ncol(df)))
  df$PC       <- paste0("PC", seq_len(nrow(df)))
  df$var_pct  <- pct_vec[seq_len(nrow(df))]
  df$endpoint <- endpoint_label
  for (ax in colnames(df)[startsWith(colnames(df), "CoIA_Ax")]) {
    df[[paste0(ax, "_w2")]] <- round(df[[ax]]^2, 4)
  }
  for (ax in paste0("CoIA_Ax", seq_len(ncol(a_mat)))) {
    w2col   <- paste0(ax, "_w2")
    rel_col <- paste0(ax, "_contrib_pct")
    df[[rel_col]] <- round(df[[w2col]] / sum(df[[w2col]]) * 100, 1)
  }
  df %>% select(endpoint, PC, var_pct,
                starts_with("CoIA_Ax1"), starts_with("CoIA_Ax2"))
}

morph_contrib <- compute_pca_cia_contribution(
  coin_arch$aX,
  round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1),
  "Morphology"
)
scar_contrib <- compute_pca_cia_contribution(
  coin_arch$aY,
  round(dudi_scar$eig / sum(dudi_scar$eig) * 100, 1),
  "Scar direction"
)

pca_cia_contrib <- bind_rows(morph_contrib, scar_contrib)
cat("\n形态端：\n"); print(morph_contrib %>% select(-endpoint) %>% as.data.frame())
cat("\n方向端：\n"); print(scar_contrib  %>% select(-endpoint) %>% as.data.frame())
write_csv(pca_cia_contrib,
          here("analysis/data/derived_data/PCA_CoIA_contribution.csv"))
cat("已保存：PCA_CoIA_contribution.csv\n")


# ------------------------------------------------------------------------------
# CoIA 辅助可视化：诊断图
# ------------------------------------------------------------------------------

eig_df <- tibble(
  axis       = paste0("Axis ", seq_along(coin_arch$eig)),
  eigenvalue = coin_arch$eig,
  pct        = coin_arch$eig / sum(coin_arch$eig) * 100,
  cum_pct    = cumsum(pct)
)

p_scree <- ggplot(eig_df, aes(x = axis, y = pct)) +
  geom_col(fill = "#5C7F71", alpha = 0.85, width = 0.55) +
  geom_line(aes(y = cum_pct, group = 1), color = "#802520", linewidth = 0.8) +
  geom_point(aes(y = cum_pct), color = "#802520", size = 2.8) +
  geom_text(aes(y = pct + 1.5, label = sprintf("%.1f%%", pct)),
            size = 3, color = "grey30") +
  scale_y_continuous(
    name     = "Explained co-inertia (%)",
    sec.axis = sec_axis(~ ., name = "Cumulative (%)")
  ) +
  theme_bw(base_size = 10) +
  labs(
    title    = "CoIA Scree Plot (SDG)",
    subtitle = sprintf("RV = %.3f, p = %.3f", coin_arch$RV, rv_test$pvalue),
    x = "CoIA Axis"
  ) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"))

morph_pct <- round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1)
scar_pct  <- round(dudi_scar$eig  / sum(dudi_scar$eig)  * 100, 1)

morph_load <- as.data.frame(coin_arch$aX) %>%
  rownames_to_column("variable") %>%
  rename(Axis1 = AxcX1, Axis2 = AxcX2) %>%
  mutate(
    pct   = morph_pct[as.integer(str_extract(variable, "[0-9]+"))],
    variable_label = sprintf("Morph-PCA%s\n(%.1f%% var)",
                             str_extract(variable, "[0-9]+"), pct),
    endpoint = "Morphology"
  )

scar_load <- as.data.frame(coin_arch$aY) %>%
  rownames_to_column("variable") %>%
  rename(Axis1 = AxcY1, Axis2 = AxcY2) %>%
  mutate(
    pct   = scar_pct[as.integer(str_extract(variable, "[0-9]+"))],
    variable_label = sprintf("Dir-PCA%s\n(%.1f%% var)",
                             str_extract(variable, "[0-9]+"), pct),
    endpoint = "Scar Direction"
  )

circle_df <- tibble(angle = seq(0, 2 * pi, length.out = 300),
                    x = cos(angle), y = sin(angle))

make_loading_plot <- function(load_df, title_str, col_fill) {
  ggplot(load_df) +
    geom_path(data = circle_df, aes(x = x, y = y),
              color = "grey80", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_segment(aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
                 arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
                 color = col_fill, linewidth = 0.9, alpha = 0.85) +
    geom_label(aes(x = Axis1 * 1.12, y = Axis2 * 1.12,
                   label = variable_label),
               size = 2.6, color = "grey20",
               label.size = 0.2, fill = "white", alpha = 0.85,
               lineheight = 0.85) +
    coord_fixed(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
    theme_bw(base_size = 10) +
    labs(title = title_str,
         x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
         y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          panel.grid = element_blank())
}

p_load_morph <- make_loading_plot(morph_load,
                                  "Morphology PCA axes on CoIA space", "#5C7F71")
p_load_scar  <- make_loading_plot(scar_load,
                                  "Scar direction PCA axes on CoIA space", "#BA8530")

p_cia_diagnostics <- (p_scree | p_load_morph | p_load_scar) +
  plot_annotation(
    title   = "CoIA Axis Diagnostics (SDG)",
    caption = paste(
      "Left: scree plot. Middle: morphology PCA axis loadings on CoIA axes.",
      "Right: scar-direction PCA axis loadings on CoIA axes.",
      "\nBoth panels use PCA-axis projections (aX / aY) for symmetric interpretation.",
      "Arrow length = contribution to CoIA structure."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0)
    )
  )

ggsave(here("analysis/output/figures/L1_CoIA_Diagnostics.png"),
       plot = p_cia_diagnostics, width = 15, height = 5.5, dpi = 300, bg = "white")
cat("图已保存：L1_CoIA_Diagnostics.png\n")


# ------------------------------------------------------------------------------
# CoIA 双标图通用函数
# show_color_legend：TRUE = 独立图显示颜色图例；FALSE = 组合图用，颜色图例移至方向图
# ------------------------------------------------------------------------------
make_coia_biplot <- function(group_col, group_label, palette,
                             group_order = NULL, fname_tag = NULL,
                             show_color_legend = TRUE) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 1) {
    cat(sprintf("  [跳过] %s CoIA 双标图：无有效分组\n", group_label))
    return(invisible(NULL))
  }
  if (!is.null(group_order)) {
    lvls <- intersect(group_order, valid_groups)
  } else {
    lvls <- sort(valid_groups)
  }
  
  sub_seg <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = lvls))
  
  scores_long <- bind_rows(
    sub_seg %>%
      select(ID, !!group_col, x = Axis1_M, y = Axis2_M, arrow_length) %>%
      mutate(endpoint = "Morphology"),
    sub_seg %>%
      select(ID, !!group_col, x = Axis1_S, y = Axis2_S, arrow_length) %>%
      mutate(endpoint = "Scar direction")
  ) %>%
    mutate(
      endpoint     = factor(endpoint, levels = c("Morphology", "Scar direction")),
      !!group_col := factor(.data[[group_col]], levels = lvls)
    )
  
  endpoint_shapes <- c("Morphology" = 21, "Scar direction" = 24)
  endpoint_sizes  <- c("Morphology" = 3.0, "Scar direction" = 2.6)
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_segment(
      data = sub_seg,
      aes(x = Axis1_M, y = Axis2_M, xend = Axis1_S, yend = Axis2_S,
          color = .data[[group_col]]),
      linewidth = 0.45, alpha = 0.45, lineend = "round"
    ) +
    geom_point(
      data = scores_long,
      aes(x = x, y = y,
          fill  = .data[[group_col]],
          color = .data[[group_col]],
          shape = endpoint,
          size  = endpoint),
      stroke = 0.5, alpha = 0.90
    ) +
    scale_color_manual(values = palette, name = group_label, breaks = lvls) +
    scale_fill_manual(values  = palette, name = group_label, breaks = lvls) +
    scale_shape_manual(values = endpoint_shapes, name = "Endpoint") +
    scale_size_manual(values  = endpoint_sizes,  name = "Endpoint") +
    theme_bw() +
    labs(
      x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
      y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])
    ) +
    guides(
      # 颜色图例：独立图显示，组合图隐藏（移至方向图）
      color = if (show_color_legend) {
        guide_legend(order = 1,
                     override.aes = list(shape = 21, size = 3),
                     title = group_label)
      } else {
        "none"
      },
      fill  = "none",
      # Endpoint shape 图例始终显示
      shape = guide_legend(order = 2,
                           override.aes = list(fill  = "grey60",
                                               color = "grey30",
                                               size  = c(3.0, 2.6)),
                           title = "Endpoint"),
      size  = "none"
    ) +
    theme(
      panel.grid.major.x   = element_blank(),
      panel.grid.major.y   = element_blank(),
      panel.grid.minor     = element_blank(),
      legend.position      = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.box           = "vertical",
      legend.box.just      = "left",
      legend.background    = element_rect(fill  = alpha("white", 0.75),
                                          color = "grey80", linewidth = 0.3),
      legend.key.size      = unit(0.45, "cm"),
      legend.text          = element_text(size = 8),
      legend.margin        = margin(4, 6, 4, 6)
    )
  
  tag   <- if (!is.null(fname_tag)) fname_tag else tolower(str_replace_all(group_label, " ", "_"))
  fname <- sprintf("analysis/output/figures/L1_CoIA_Biplot_%s.png", tag)
  ggsave(here(fname), plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  p
}

# 独立图（保留颜色图例）
make_coia_biplot("layer",        "Layer",        layer_pal,    LAYER_ORDER,    "layer")
make_coia_biplot("raw_material", "Raw Material", rawmat_pal,   NULL,           "raw_material")
make_coia_biplot("core_type",    "Core Type",    coretype_pal, CORETYPE_ORDER, "core_type")

# 组合图用版本（颜色图例隐藏，移至方向图右侧）
p_coia_layer    <- make_coia_biplot("layer",        "Layer",        layer_pal,
                                    LAYER_ORDER,    "layer",        show_color_legend = TRUE)
p_coia_rawmat   <- make_coia_biplot("raw_material", "Raw Material", rawmat_pal,
                                    NULL,           "raw_material", show_color_legend = TRUE)
p_coia_coretype <- make_coia_biplot("core_type",    "Core Type",    coretype_pal,
                                    CORETYPE_ORDER, "core_type",    show_color_legend = TRUE)

l1_results <- tibble(
  method  = c("Mantel (ILR, Euclidean, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_arch$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(arch_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/L1_results.csv"))
cat("\n第一层结论：\n"); print(l1_results)


# ==============================================================================
# ---- L1-3：CoIA 桑基图（ILR -> PCA -> CoIA 贡献流）----
# ==============================================================================

cat("\n==== L1-3：CoIA 桑基图 ====\n")

c1_morph <- as.matrix(dudi_morph$c1)
c1_scar  <- as.matrix(dudi_scar$c1)

colnames(c1_morph) <- paste0("MPC", seq_len(ncol(c1_morph)))
colnames(c1_scar)  <- paste0("DPC", seq_len(ncol(c1_scar)))
rownames(c1_morph) <- paste0("Milr", seq_len(nrow(c1_morph)))
rownames(c1_scar)  <- paste0("Dilr", seq_len(nrow(c1_scar)))

w2_ilr_mpc <- c1_morph^2
w2_ilr_dpc <- c1_scar^2

a_morph <- as.matrix(coin_arch$aX)
a_scar  <- as.matrix(coin_arch$aY)

w2_mpc_cia <- a_morph^2
w2_dpc_cia <- a_scar^2

morph_var_pct_full <- round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1)
scar_var_pct_full  <- round(dudi_scar$eig  / sum(dudi_scar$eig)  * 100, 1)

keep_mpc <- which(cumsum(morph_var_pct_full) <= 95 | seq_along(morph_var_pct_full) == 1)
keep_dpc <- which(cumsum(scar_var_pct_full)  <= 95 | seq_along(scar_var_pct_full)  == 1)
keep_mpc <- keep_mpc[keep_mpc <= ncol(c1_morph)]
keep_dpc <- keep_dpc[keep_dpc <= ncol(c1_scar)]
n_cia_ax <- min(2, ncol(a_morph))

w2_ilr_mpc <- w2_ilr_mpc[, keep_mpc, drop = FALSE]
w2_ilr_dpc <- w2_ilr_dpc[, keep_dpc, drop = FALSE]
w2_mpc_cia <- w2_mpc_cia[keep_mpc, seq_len(n_cia_ax), drop = FALSE]
w2_dpc_cia <- w2_dpc_cia[keep_dpc, seq_len(n_cia_ax), drop = FALSE]

n_ilr_m <- nrow(w2_ilr_mpc)
n_ilr_d <- nrow(w2_ilr_dpc)
n_mpc   <- ncol(w2_ilr_mpc)
n_dpc   <- ncol(w2_ilr_dpc)

cat(sprintf(
  "  形态端：%d ILR -> %d MorphPC；方向端：%d ILR -> %d DirPC；共享 %d CoIA 轴\n",
  n_ilr_m, n_mpc, n_ilr_d, n_dpc, n_cia_ax
))

h_ilr_m_raw <- rowSums(w2_ilr_mpc)
h_ilr_d_raw <- rowSums(w2_ilr_dpc)
h_mpc_raw   <- rowSums(w2_mpc_cia)
h_dpc_raw   <- rowSums(w2_dpc_cia)
h_cia_raw   <- colSums(w2_mpc_cia) + colSums(w2_dpc_cia)

MORPH_HEIGHT <- 220
SCAR_HEIGHT  <- 220
NODE_GAP     <- 5

scale_nodes <- function(h_raw, total_h, n_nodes, gap = NODE_GAP) {
  avail <- total_h - gap * (n_nodes - 1)
  h_raw * (avail / sum(h_raw))
}

h_ilr_m <- scale_nodes(h_ilr_m_raw, MORPH_HEIGHT, n_ilr_m)
h_mpc   <- scale_nodes(h_mpc_raw,   MORPH_HEIGHT, n_mpc)
h_ilr_d <- scale_nodes(h_ilr_d_raw, SCAR_HEIGHT,  n_ilr_d)
h_dpc   <- scale_nodes(h_dpc_raw,   SCAR_HEIGHT,  n_dpc)
h_cia   <- scale_nodes(h_cia_raw,
                       MORPH_HEIGHT + SCAR_HEIGHT + NODE_GAP * (n_cia_ax - 1),
                       n_cia_ax)

y_start <- 52

y_ilr_m <- y_start + cumsum(c(0, head(h_ilr_m + NODE_GAP, -1)))
y_mpc   <- y_start + cumsum(c(0, head(h_mpc   + NODE_GAP, -1)))

y_ilr_d <- y_start + MORPH_HEIGHT + NODE_GAP +
  cumsum(c(0, head(h_ilr_d + NODE_GAP, -1)))
y_dpc   <- y_start + MORPH_HEIGHT + NODE_GAP +
  cumsum(c(0, head(h_dpc   + NODE_GAP, -1)))

cia_total_h <- sum(h_cia) + NODE_GAP * (n_cia_ax - 1)
two_ends_h  <- MORPH_HEIGHT + NODE_GAP + SCAR_HEIGHT
cia_offset  <- (two_ends_h - cia_total_h) / 2
y_cia <- y_start + cia_offset +
  cumsum(c(0, head(h_cia + NODE_GAP, -1)))

SVG_W  <- 680
SVG_H  <- ceiling(y_start + MORPH_HEIGHT + NODE_GAP + SCAR_HEIGHT + 30)
NODE_W <- 100
x_col1 <- 20
x_col2 <- 240
x_col3 <- 470

col_morph_fill   <- "#B5D4F4"; col_morph_stroke <- "#185FA5"
col_morph_text   <- "#0C447C"; col_morph_band   <- "#85B7EB"
col_scar_fill    <- "#F5C4B3"; col_scar_stroke  <- "#993C1D"
col_scar_text    <- "#712B13"; col_scar_band    <- "#F0997B"
col_cia_fill     <- "#CECBF6"; col_cia_stroke   <- "#534AB7"
col_cia_text     <- "#3C3489"
col_band_cia_m   <- "#AFA9EC"
col_band_cia_d   <- "#C9A8E0"

svg_rect <- function(x, y, w, h, fill, stroke, rx = 5, sw = 0.8)
  sprintf(
    '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%d" fill="%s" stroke="%s" stroke-width="%.1f"/>',
    x, y, w, h, rx, fill, stroke, sw
  )

svg_band <- function(x1, y1t, y1b, x2, y2t, y2b, fill, opacity = 0.42) {
  xm <- (x1 + x2) / 2
  sprintf(
    paste0('<path d="M%.1f,%.1f C%.1f,%.1f %.1f,%.1f %.1f,%.1f ',
           'L%.1f,%.1f C%.1f,%.1f %.1f,%.1f %.1f,%.1f Z" ',
           'fill="%s" fill-opacity="%.2f"/>'),
    x1, y1t, xm, y1t, xm, y2t, x2, y2t,
    x2, y2b, xm, y2b, xm, y1b, x1, y1b,
    fill, opacity
  )
}

svg_text <- function(x, y, txt, size = 11, weight = 400,
                     anchor = "middle", fill = "#2C2C2A")
  sprintf(
    '<text x="%.1f" y="%.1f" font-size="%d" font-weight="%d" text-anchor="%s" fill="%s" dominant-baseline="central">%s</text>',
    x, y, size, weight, anchor, fill, txt
  )

lines <- character(0)
push  <- function(...) lines <<- c(lines, ...)

push(sprintf(
  '<svg xmlns="http://www.w3.org/2000/svg" width="%d" viewBox="0 0 %d %d">',
  SVG_W, SVG_W, SVG_H
))
push('<title>ILR to PCA to CoIA contribution flow (SDG)</title>')
push('<desc>Sankey-style flow diagram showing ILR spectral variable contributions through PCA axes to CoIA axes.</desc>')
push(sprintf('<rect width="%d" height="%d" fill="white"/>', SVG_W, SVG_H))

push(svg_text(x_col1 + NODE_W / 2, 20, "ILR variables",    12, 500, fill = "#2C2C2A"))
push(svg_text(x_col1 + NODE_W / 2, 36, "spectral log-contrasts", 10, 400, fill = "#5F5E5A"))
push(svg_text(x_col2 + NODE_W / 2, 20, "PCA axes",          12, 500, fill = "#2C2C2A"))
push(svg_text(x_col2 + NODE_W / 2, 36, "per endpoint",      10, 400, fill = "#5F5E5A"))
push(svg_text(x_col3 + NODE_W / 2, 20, "CoIA axes",         12, 500, fill = "#2C2C2A"))
push(svg_text(x_col3 + NODE_W / 2, 36, "shared structure",  10, 400, fill = "#5F5E5A"))

outlet_ilr_m <- rep(0, n_ilr_m); inlet_mpc  <- rep(0, n_mpc)
outlet_ilr_d <- rep(0, n_ilr_d); inlet_dpc  <- rep(0, n_dpc)
inlet_cia    <- rep(0, n_cia_ax)

for (i in seq_len(n_ilr_m)) for (j in seq_len(n_mpc)) {
  ww <- w2_ilr_mpc[i, j]; if (ww < 0.004) next
  bw_out <- ww * h_ilr_m[i]; bw_in <- ww * h_mpc[j]
  y1t <- y_ilr_m[i] + outlet_ilr_m[i]; y1b <- y1t + bw_out
  y2t <- y_mpc[j]   + inlet_mpc[j];    y2b <- y2t + bw_in
  push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                col_morph_band, 0.28 + 0.32 * ww))
  outlet_ilr_m[i] <- outlet_ilr_m[i] + bw_out
  inlet_mpc[j]    <- inlet_mpc[j]    + bw_in
}

for (i in seq_len(n_ilr_d)) for (j in seq_len(n_dpc)) {
  ww <- w2_ilr_dpc[i, j]; if (ww < 0.004) next
  bw_out <- ww * h_ilr_d[i]; bw_in <- ww * h_dpc[j]
  y1t <- y_ilr_d[i] + outlet_ilr_d[i]; y1b <- y1t + bw_out
  y2t <- y_dpc[j]   + inlet_dpc[j];    y2b <- y2t + bw_in
  push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                col_scar_band, 0.28 + 0.32 * ww))
  outlet_ilr_d[i] <- outlet_ilr_d[i] + bw_out
  inlet_dpc[j]    <- inlet_dpc[j]    + bw_in
}

outlet_mpc <- rep(0, n_mpc)
for (j in seq_len(n_mpc)) for (k in seq_len(n_cia_ax)) {
  ww <- w2_mpc_cia[j, k]; if (ww < 0.004) next
  tot_k  <- colSums(w2_mpc_cia)[k] + colSums(w2_dpc_cia)[k]
  bw_out <- ww * h_mpc[j]
  bw_in  <- ww * h_cia[k] * (colSums(w2_mpc_cia)[k] / tot_k)
  y1t <- y_mpc[j] + outlet_mpc[j]; y1b <- y1t + bw_out
  y2t <- y_cia[k] + inlet_cia[k];  y2b <- y2t + bw_in
  push(svg_band(x_col2 + NODE_W, y1t, y1b, x_col3, y2t, y2b,
                col_band_cia_m, 0.22 + 0.36 * ww))
  outlet_mpc[j] <- outlet_mpc[j] + bw_out
  inlet_cia[k]  <- inlet_cia[k]  + bw_in
}

outlet_dpc <- rep(0, n_dpc)
for (j in seq_len(n_dpc)) for (k in seq_len(n_cia_ax)) {
  ww <- w2_dpc_cia[j, k]; if (ww < 0.004) next
  tot_k  <- colSums(w2_mpc_cia)[k] + colSums(w2_dpc_cia)[k]
  bw_out <- ww * h_dpc[j]
  bw_in  <- ww * h_cia[k] * (colSums(w2_dpc_cia)[k] / tot_k)
  y1t <- y_dpc[j] + outlet_dpc[j]; y1b <- y1t + bw_out
  y2t <- y_cia[k] + inlet_cia[k];  y2b <- y2t + bw_in
  push(svg_band(x_col2 + NODE_W, y1t, y1b, x_col3, y2t, y2b,
                col_band_cia_d, 0.22 + 0.36 * ww))
  outlet_dpc[j] <- outlet_dpc[j] + bw_out
  inlet_cia[k]  <- inlet_cia[k]  + bw_in
}

for (i in seq_len(n_ilr_m)) {
  push(svg_rect(x_col1, y_ilr_m[i], NODE_W, h_ilr_m[i], col_morph_fill, col_morph_stroke))
  push(svg_text(x_col1 + NODE_W / 2, y_ilr_m[i] + h_ilr_m[i] / 2,
                sprintf("ilr%d (morph)", i), 11, 500, fill = col_morph_text))
}
for (i in seq_len(n_ilr_d)) {
  push(svg_rect(x_col1, y_ilr_d[i], NODE_W, h_ilr_d[i], col_scar_fill, col_scar_stroke))
  push(svg_text(x_col1 + NODE_W / 2, y_ilr_d[i] + h_ilr_d[i] / 2,
                sprintf("ilr%d (dir)", i), 11, 500, fill = col_scar_text))
}
for (j in seq_len(n_mpc)) {
  push(svg_rect(x_col2, y_mpc[j], NODE_W, h_mpc[j], col_morph_fill, col_morph_stroke))
  push(svg_text(x_col2 + NODE_W / 2, y_mpc[j] + h_mpc[j] / 2,
                sprintf("Morph-PC%d", j), 11, 500, fill = col_morph_text))
}
for (j in seq_len(n_dpc)) {
  push(svg_rect(x_col2, y_dpc[j], NODE_W, h_dpc[j], col_scar_fill, col_scar_stroke))
  push(svg_text(x_col2 + NODE_W / 2, y_dpc[j] + h_dpc[j] / 2,
                sprintf("Dir-PC%d", j), 11, 500, fill = col_scar_text))
}
for (k in seq_len(n_cia_ax)) {
  push(svg_rect(x_col3, y_cia[k], NODE_W, h_cia[k], col_cia_fill, col_cia_stroke))
  push(svg_text(x_col3 + NODE_W / 2, y_cia[k] + h_cia[k] / 2,
                c("CoIA Axis 1", "CoIA Axis 2")[k], 11, 500, fill = col_cia_text))
}

push("</svg>")

svg_path <- tempfile(fileext = ".svg")
png_path <- here("analysis/output/figures/L1_CoIA_Sankey.png")
writeLines(lines, svg_path, useBytes = FALSE)

if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_png(svg_path, png_path, width = SVG_W * 2)
  cat("PNG 已保存（via rsvg）：L1_CoIA_Sankey.png\n")
} else if (nzchar(Sys.which("rsvg-convert"))) {
  system2("rsvg-convert",
          args = c("-d", "300", "-p", "300", "-o", png_path, svg_path))
  cat("PNG 已保存（via rsvg-convert）：L1_CoIA_Sankey.png\n")
} else if (nzchar(Sys.which("inkscape"))) {
  system2("inkscape",
          args = c("--export-filename", png_path, "--export-dpi", "300", svg_path))
  cat("PNG 已保存（via Inkscape）：L1_CoIA_Sankey.png\n")
} else {
  cat("[警告] 未检测到 rsvg / rsvg-convert / Inkscape，无法生成 PNG。\n")
  cat("       请运行 install.packages('rsvg') 后重新执行此段。\n")
  cat("       SVG 中间文件位于：", svg_path, "\n")
}
unlink(svg_path)

cat("\n==== L1-3 桑基图完成 ====\n")


# ==============================================================================
# ========== 第二层：联合证据 ==========
# ==============================================================================

cat("\n\n")
cat("##  第二层：联合证据                                          ##\n")


# ------------------------------------------------------------------------------
# L2-A：分组 Mantel
# ------------------------------------------------------------------------------

cat("\n---------- L2-A：分组 Mantel ----------\n")

mantel_within_group <- function(group_val, group_col, D_morph_full,
                                D_scar_full, meta_df, n_perm = 9999) {
  ids <- meta_df %>%
    filter(.data[[group_col]] == group_val) %>%
    pull(ID)
  cat(sprintf("  -> %s (n = %d) ...", group_val, length(ids)))
  if (length(ids) < 5) {
    cat(" 跳过（n < 5）\n"); return(NULL)
  }
  res <- mantel(extract_subdist(D_morph_full, ids),
                extract_subdist(D_scar_full,  ids),
                method = "spearman", permutations = n_perm)
  cat(sprintf(" r = %.4f, p = %.4f\n", res$statistic, res$signif))
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
                                D_morph_arch, D_scar_arch, meta_arch))
}

mantel_by_layer    <- run_grouped_mantel(meta_layer,    "layer",        "层位")
mantel_by_rawmat   <- run_grouped_mantel(meta_rawmat,   "raw_material", "原料")
mantel_by_coretype <- run_grouped_mantel(meta_coretype, "core_type",    "石核类型")

cat("\n层位分组 Mantel：\n")
if (!is.null(mantel_by_layer))
  print(mantel_by_layer %>% mutate(across(c(mantel_r, p_raw), ~ round(.x, 4))))
cat("\n原料分组 Mantel：\n")
if (!is.null(mantel_by_rawmat))
  print(mantel_by_rawmat %>% mutate(across(c(mantel_r, p_raw), ~ round(.x, 4))))
cat("\n石核类型分组 Mantel：\n")
if (!is.null(mantel_by_coretype))
  print(mantel_by_coretype %>% mutate(across(c(mantel_r, p_raw), ~ round(.x, 4))))

l2_mantel <- bind_rows(
  mantel_by_layer, mantel_by_rawmat, mantel_by_coretype
) %>%
  mutate(
    p_holm = p.adjust(p_raw, method = "holm"),
    significance = case_when(
      p_holm < 0.001 ~ "***", p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",   p_holm < 0.10  ~ ".",
      TRUE           ~ "ns"
    ),
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
  write_csv(l2_mantel, here("analysis/data/derived_data/L2_grouped_mantel.csv"))
  cat("已保存：L2_grouped_mantel.csv\n")
}


# ------------------------------------------------------------------------------
# L2-B：CoIA 箭头长度
# ------------------------------------------------------------------------------

cat("\n---------- L2-B：箭头长度分组差异 ----------\n")

run_arrow_length_analysis <- function(group_col, group_label, palette,
                                      group_order = NULL) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 3) %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [跳过] %s 箭头长度检验：有效分组不足\n", group_label))
    return(invisible(NULL))
  }
  
  if (!is.null(group_order)) {
    valid_groups_plot <- intersect(group_order, valid_groups)
  } else {
    valid_groups_plot <- valid_groups
  }
  
  sub_df_stat <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups)
  cat(sprintf("\n----- %s x 箭头长度 -----\n", group_label))
  kw <- kruskal.test(reformulate(group_col, "arrow_length"), data = sub_df_stat)
  print(kw)
  pw <- pairwise.wilcox.test(sub_df_stat$arrow_length, sub_df_stat[[group_col]],
                             p.adjust.method = "holm", exact = FALSE)
  print(pw)
  
  sub_df <- scores_combined %>%
    filter(.data[[group_col]] %in% valid_groups_plot) %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_plot))
  
  p <- ggplot(sub_df,
              aes(x = .data[[group_col]], y = arrow_length,
                  fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2.5,
                 alpha = 0.25, linewidth = 0.5) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.7, shape = 16) +
    stat_summary(fun = mean, geom = "point",
                 shape = 16, size = 4, color = "white") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Kruskal-Wallis\nchi\u00b2 = %.2f, P = %.3f",
                             kw$statistic, kw$p.value),
             hjust = 1.05, vjust = 1.2, size = 4, color = "grey40") +
    scale_fill_manual(values  = palette) +
    scale_color_manual(values = palette) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = 9.5),
      legend.position    = "none"
    ) +
    labs(x = NULL, y = "CoIA line length")
  
  fname <- sprintf("analysis/output/figures/L2_Arrow_Length_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  list(sub_df = sub_df, kw = kw, pw = pw, p = p)
}

res_len_layer    <- run_arrow_length_analysis("layer",        "Layer",
                                              layer_pal,    LAYER_ORDER)
res_len_rawmat   <- run_arrow_length_analysis("raw_material", "Raw_Material",
                                              rawmat_pal)
res_len_coretype <- run_arrow_length_analysis("core_type",    "Core_Type",
                                              coretype_pal, CORETYPE_ORDER)

cat("\n==== 箭头长度描述统计 ====\n")
for (gc in c("layer", "raw_material", "core_type")) {
  cat(sprintf("\n--- %s ---\n", gc))
  scores_combined %>%
    filter(!is.na(.data[[gc]]), .data[[gc]] != "Other") %>%
    group_by(.data[[gc]]) %>%
    summarise(n      = n(),
              mean   = round(mean(arrow_length),   4),
              median = round(median(arrow_length), 4),
              sd     = round(sd(arrow_length),     4),
              min    = round(min(arrow_length),    4),
              max    = round(max(arrow_length),    4),
              .groups = "drop") %>%
    print()
}


# ------------------------------------------------------------------------------
# L2-C：CoIA 箭头方位圆形统计（线性展示）
# plot_rose() 支持 show_color_legend 参数：
#   TRUE  = 独立图，颜色图例显示在右侧（默认）
#   FALSE = 组合图内部调用，不显示颜色图例
#   "right_panel" = 组合图最右列，显示颜色图例（legend 放在图右侧）
# ------------------------------------------------------------------------------

cat("\n---------- L2-C：箭头方位圆形统计 ----------\n")

plot_rose <- function(res, palette, show_color_legend = TRUE) {
  group_col <- res$group_col
  kde_df    <- res$kde_df
  mean_dirs <- res$mean_dirs
  
  kde_linear <- kde_df %>%
    mutate(angle_centered = if_else(angle_deg > 180, angle_deg - 360, angle_deg),
           !!group_col := factor(.data[[group_col]],
                                 levels = levels(kde_df[[group_col]])))
  
  mean_linear <- mean_dirs %>%
    mutate(angle_centered = if_else(mean_deg > 180, mean_deg - 360, mean_deg),
           !!group_col := factor(.data[[group_col]],
                                 levels = levels(kde_df[[group_col]])))
  
  rayleigh_labels <- res$rayleigh %>%
    mutate(
      label = case_when(
        rayleigh_p < 0.001 ~ "Rayleigh\nP < 0.001",
        rayleigh_p < 0.05  ~ sprintf("Rayleigh\nP = %.3f", rayleigh_p),
        TRUE               ~ sprintf("Rayleigh\nP = %.3f", rayleigh_p)
      ),
      !!group_col := factor(group, levels = levels(kde_df[[group_col]]))
    )
  
  # 图例位置与显示策略
  if (isTRUE(show_color_legend)) {
    # 独立图：右侧显示颜色图例
    legend_pos    <- "right"
    color_guide   <- guide_legend(title = group_col)
  } else if (identical(show_color_legend, "right_panel")) {
    # 组合图最右列：右侧显示颜色图例
    legend_pos    <- "right"
    color_guide   <- guide_legend(title = group_col)
  } else {
    # 组合图非末列：隐藏颜色图例
    legend_pos    <- "none"
    color_guide   <- "none"
  }
  
  ggplot(kde_linear,
         aes(x     = angle_centered,
             y     = density,
             fill  = .data[[group_col]],
             color = .data[[group_col]])) +
    geom_area(alpha = 0.4, color = NULL) +
    geom_vline(
      data      = mean_linear,
      aes(xintercept = angle_centered, color = .data[[group_col]]),
      linewidth = 0.5, linetype = "dashed", alpha = 0.75
    ) +
    geom_text(
      data = rayleigh_labels,
      aes(label = label),
      x = 170, y = Inf,
      hjust = 0.8, vjust = 1.4,
      size = 3, color = "grey35",
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(-180, 180),
      breaks = seq(-180, 180, by = 45),
      labels = seq(-180, 180, by = 45)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values  = palette, guide = color_guide) +
    scale_color_manual(values = palette, guide = "none") +
    # 组合图内：ncol=1 垂直分面，分面标签隐藏
    facet_wrap(reformulate(group_col), ncol = 1, scales = "free_y") +
    labs(x    = "CoIA line direction (\u00b0)",
         y    = "von Mises KDE") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      panel.grid.major.y  = element_blank(),
      # 分面标签：独立图时显示，组合图时隐藏
      strip.text          = if (isTRUE(show_color_legend)) {
        element_text(face = "bold", size = 9)
      } else {
        element_blank()
      },
      strip.background    = if (isTRUE(show_color_legend)) {
        element_rect(fill = "#EBEBEB", color = "#EBEBEB")
      } else {
        element_blank()
      },
      axis.text.x         = element_text(size = 7.5),
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      legend.position     = legend_pos,
      legend.key.size     = unit(0.45, "cm"),
      legend.text         = element_text(size = 8),
      legend.title        = element_text(size = 8.5, face = "bold")
    )
}

# ---- 分析函数 ----
run_circular_analysis <- function(group_col, group_label, palette,
                                  group_order = NULL) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 5) %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [跳过] %s 圆形统计：有效组数不足\n", group_label))
    return(invisible(NULL))
  }
  
  if (!is.null(group_order)) {
    valid_groups_ordered <- intersect(group_order, valid_groups)
  } else {
    valid_groups_ordered <- sort(valid_groups)
  }
  
  sub_df <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups)
  
  cat(sprintf("\n----- %s 圆形描述统计 -----\n", group_label))
  circ_desc <- map_dfr(valid_groups_ordered, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(group_var = group_col, group = g, n = length(angles),
           mean_dir_deg = round(cs$mean_deg, 2),
           concentration_r = round(cs$rho, 4))
  })
  print(circ_desc)
  
  cat(sprintf("\n----- %s Rayleigh 检验 -----\n", group_label))
  rayleigh_res <- map_dfr(valid_groups_ordered, function(g) {
    angles   <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    circ_obj <- circular(angles, type = "angles", units = "radians", modulo = "2pi")
    rt       <- rayleigh.test(circ_obj)
    cat(sprintf("  %s: U = %.4f, p = %.4f -> %s\n",
                g, rt$statistic, rt$p.value,
                ifelse(rt$p.value < 0.05, "方位集中", "方位分散")))
    tibble(group_var = group_col, group = g,
           rayleigh_U = round(rt$statistic, 4),
           rayleigh_p = round(rt$p.value,   4),
           conclusion = ifelse(rt$p.value < 0.05, "concentrated", "uniform"))
  })
  
  watson_res <- NULL
  if (length(valid_groups) >= 2) {
    cat(sprintf("\n----- %s Watson 两样本检验 -----\n", group_label))
    pairs <- combn(valid_groups, 2, simplify = FALSE)
    watson_res <- map_dfr(pairs, function(pair) {
      x1 <- sub_df %>% filter(.data[[group_col]] == pair[1]) %>% pull(arrow_angle)
      x2 <- sub_df %>% filter(.data[[group_col]] == pair[2]) %>% pull(arrow_angle)
      wt <- watson_perm_test(x1, x2, B = 9999)
      cat(sprintf("  %s vs %s: U2 = %.4f, p = %.4f -> %s\n",
                  pair[1], pair[2], wt$statistic, wt$p.value,
                  ifelse(wt$p.value < 0.05, "分布不同", "差异不显著")))
      tibble(group_var    = group_col,
             group1       = pair[1], group2 = pair[2],
             U2_statistic = round(wt$statistic, 4),
             p_value      = round(wt$p.value,   4),
             conclusion   = ifelse(wt$p.value < 0.05, "different", "ns"))
    })
  }
  
  mean_dirs <- map_dfr(valid_groups_ordered, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(!!group_col := g, mean_deg = cs$mean_deg)
  }) %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_ordered))
  
  compute_circular_kde <- function(angles_deg, bw = 25, n = 360) {
    circ <- circular(angles_deg * pi / 180,
                     type = "angles", units = "radians", modulo = "2pi")
    dens <- density(circ, bw = bw, n = n)
    angle_deg <- as.numeric(dens$x) * 180 / pi %% 360
    density   <- as.numeric(dens$y)
    tibble(angle_deg = c(angle_deg, angle_deg[1]),
           density   = c(density,   density[1]))
  }
  
  kde_df <- sub_df %>%
    mutate(
      angle_deg = arrow_angle * 180 / pi,
      angle_deg = ifelse(angle_deg < 0, angle_deg + 360, angle_deg),
      !!group_col := factor(.data[[group_col]], levels = valid_groups_ordered)
    ) %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~ compute_circular_kde(.x$angle_deg, bw = 25)) %>%
    ungroup() %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_ordered))
  
  # 独立图：用默认 show_color_legend = TRUE（横向分面，有标签）
  p_rose_standalone <- plot_rose(
    list(group_col = group_col, kde_df = kde_df,
         mean_dirs = mean_dirs, rayleigh = rayleigh_res),
    palette,
    show_color_legend = TRUE
  )
  
  n_g   <- length(valid_groups_ordered)
  fname <- sprintf("analysis/output/figures/L2_Arrow_Direction_rose_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p_rose_standalone,
         width = min(n_g * 2.8 + 1, 18), height = 4.5, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  
  # 组合图用：传入原始数据和配置，供组合图阶段按需重建
  list(desc      = circ_desc,
       rayleigh  = rayleigh_res,
       watson    = watson_res,
       kde_df    = kde_df,
       mean_dirs = mean_dirs,
       group_col = group_col,
       palette   = palette,
       p_rose    = p_rose_standalone)
}

res_circ_layer    <- run_circular_analysis("layer",        "Layer",
                                           layer_pal,    LAYER_ORDER)
res_circ_rawmat   <- run_circular_analysis("raw_material", "Raw_Material",
                                           rawmat_pal)
res_circ_coretype <- run_circular_analysis("core_type",    "Core_Type",
                                           coretype_pal, CORETYPE_ORDER)


# ------------------------------------------------------------------------------
# L2-D：spectral_entropy × Layer / Raw_mat / Core_type
# ------------------------------------------------------------------------------

cat("\n---------- L2-D：spectral_entropy 分组分析 ----------\n")

se_df_base <- meta_arch %>%
  filter(ID %in% arch_ids) %>%
  select(ID, layer, raw_material, core_type,
         SE_direction, SE_morphology)

run_kw_dunn_se <- function(df, y_col, group_col, label) {
  sub_df <- df %>%
    filter(!is.na(.data[[y_col]]), !is.na(.data[[group_col]]),
           .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 3) %>% ungroup() %>%
    mutate(!!group_col := fct_reorder(.data[[group_col]],
                                      .data[[y_col]], median, na.rm = TRUE))
  n_groups <- length(unique(sub_df[[group_col]]))
  if (n_groups < 2) {
    cat(sprintf("  [跳过] %s：有效分组不足\n", label))
    return(NULL)
  }
  cat(sprintf("\n----- %s -----\n", label))
  kw <- kruskal.test(reformulate(group_col, y_col), data = sub_df)
  cat(sprintf("  Kruskal-Wallis: chi^2 = %.4f, df = %d, p = %.4f -> %s\n",
              kw$statistic, kw$parameter, kw$p.value,
              ifelse(kw$p.value < 0.05, "组间有显著差异", "差异不显著")))
  
  if (n_groups == 2) {
    cat("  [注] 仅2组，使用 Wilcoxon 检验替代 Dunn 检验\n")
    grp_levels <- levels(sub_df[[group_col]])
    wt <- wilcox.test(
      sub_df[[y_col]][sub_df[[group_col]] == grp_levels[1]],
      sub_df[[y_col]][sub_df[[group_col]] == grp_levels[2]],
      exact = FALSE
    )
    p_raw <- wt$p.value
    p_adj <- p_raw
    dunn <- tibble(
      Comparison   = paste(grp_levels[1], "-", grp_levels[2]),
      group1       = grp_levels[1],
      group2       = grp_levels[2],
      statistic    = as.numeric(wt$statistic),
      p            = p_raw,
      p.signif     = case_when(
        p_raw < 0.001 ~ "***", p_raw < 0.01 ~ "**",
        p_raw < 0.05  ~ "*",   p_raw < 0.10 ~ ".", TRUE ~ "ns"),
      p.adj        = p_adj,
      p.adj.signif = case_when(
        p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
        p_adj < 0.05  ~ "*",   p_adj < 0.10 ~ ".", TRUE ~ "ns")
    )
    cat(sprintf("  Wilcoxon: W = %.1f, p = %.4f (%s)\n",
                wt$statistic, p_raw, dunn$p.signif))
  } else {
    dunn_raw <- dunnTest(x = sub_df[[y_col]], g = sub_df[[group_col]],
                         method = "holm")$res
    dunn <- dunn_raw %>%
      separate(Comparison, into = c("group1", "group2"),
               sep = " - ", remove = FALSE) %>%
      rename(statistic = Z, p = P.unadj, p.adj = P.adj) %>%
      mutate(
        p.signif     = case_when(
          p     < 0.001 ~ "***", p     < 0.01 ~ "**",
          p     < 0.05  ~ "*",   p     < 0.10 ~ ".", TRUE ~ "ns"),
        p.adj.signif = case_when(
          p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
          p.adj < 0.05  ~ "*",   p.adj < 0.10 ~ ".", TRUE ~ "ns")
      ) %>%
      select(Comparison, group1, group2, statistic, p, p.signif, p.adj, p.adj.signif)
    cat("  Dunn 事后检验（Holm）：\n")
    print(dunn)
  }
  list(kw = kw, dunn = dunn, sub_df = sub_df,
       y_col = y_col, group_col = group_col)
}

se_desc_stats <- function(df, y_col, group_col) {
  df %>%
    filter(!is.na(.data[[y_col]]), !is.na(.data[[group_col]]),
           .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>%
    summarise(n          = n(),
              mean_val   = round(mean(.data[[y_col]],   na.rm = TRUE), 4),
              median_val = round(median(.data[[y_col]], na.rm = TRUE), 4),
              sd_val     = round(sd(.data[[y_col]],     na.rm = TRUE), 4),
              .groups = "drop") %>%
    mutate(SE_type    = y_col,
           group_type = group_col)
}

make_se_boxplot <- function(res_obj, y_label, title_suffix, fname, palette) {
  if (is.null(res_obj)) return(invisible(NULL))
  df        <- res_obj$sub_df
  y_col     <- res_obj$y_col
  group_col <- res_obj$group_col
  kw_res    <- res_obj$kw
  dunn_res  <- res_obj$dunn
  type_lvls <- levels(df[[group_col]])
  y_vals    <- df[[y_col]]
  y_max     <- max(y_vals, na.rm = TRUE)
  y_range   <- diff(range(y_vals, na.rm = TRUE))
  step      <- y_range * 0.10
  
  sig_pairs <- dunn_res %>% filter(p.adj < 0.05)
  
  sig_annot <- if (nrow(sig_pairs) > 0) {
    sig_pairs %>%
      mutate(annot_type = "sig",
             group1 = as.character(group1),
             group2 = as.character(group2),
             x1    = match(group1, type_lvls),
             x2    = match(group2, type_lvls),
             y_bar = y_max + step * row_number(),
             x_mid = (x1 + x2) / 2,
             label = p.adj.signif) %>%
      filter(!is.na(x1), !is.na(x2))
  } else NULL
  
  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[y_col]],
                      fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25, linewidth = 0.5) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.7, shape = 16) +
    stat_summary(fun = mean, geom = "point",
                 shape = 16, size = 4, color = "white") +
    {
      if (!is.null(sig_annot) && nrow(sig_annot) > 0) {
        tip <- y_range * 0.012
        list(
          geom_segment(data = sig_annot,
                       aes(x = x1, xend = x2, y = y_bar, yend = y_bar,
                           linetype = annot_type),
                       inherit.aes = FALSE, color = "grey30", linewidth = 0.4),
          geom_segment(data = sig_annot,
                       aes(x = x1, xend = x1, y = y_bar - tip, yend = y_bar),
                       inherit.aes = FALSE, color = "grey30", linewidth = 0.4,
                       linetype = "solid"),
          geom_segment(data = sig_annot,
                       aes(x = x2, xend = x2, y = y_bar - tip, yend = y_bar),
                       inherit.aes = FALSE, color = "grey30", linewidth = 0.4,
                       linetype = "solid"),
          geom_text(data = sig_annot,
                    aes(x = x_mid, y = y_bar + y_range * 0.015, label = label),
                    inherit.aes = FALSE, size = 3.2, color = "grey25"),
          scale_linetype_manual(values = c("sig" = "solid"), guide = "none")
        )
      } else NULL
    } +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Kruskal-Wallis\nchi\u00b2 = %.2f, p = %.3f",
                             kw_res$statistic, kw_res$p.value),
             hjust = 1.05, vjust = 1.2, size = 4,
             color = ifelse(kw_res$p.value < 0.05, "#802520", "grey50")) +
    scale_fill_manual(values  = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    scale_x_discrete(expand = expansion(add = 0.6)) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.15))) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 9.5),
      axis.text.y        = element_text(size = 9.5),
      legend.position    = "none"
    ) +
    labs(x = NULL, y = y_label)
  
  n_grp <- length(unique(df[[group_col]]))
  ggsave(here(fname), plot = p,
         width  = max(6, n_grp * 1.3 + 2),
         height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  p
}

# ---- 按 Layer 分析 ----
cat("\n==== SE x Layer ====\n")
res_se_dir_layer   <- run_kw_dunn_se(se_df_base, "SE_direction",  "layer", "SE（方向）× Layer")
res_se_morph_layer <- run_kw_dunn_se(se_df_base, "SE_morphology", "layer", "SE（形态）× Layer")

p_se_dir_layer_box   <- make_se_boxplot(res_se_dir_layer,   "Spectral Entropy (Scar Direction)", "Layer",
                                        "analysis/output/figures/L2D_SE_Direction_Layer_boxplot.png",   layer_pal)
p_se_morph_layer_box <- make_se_boxplot(res_se_morph_layer, "Spectral Entropy (Morphology)",     "Layer",
                                        "analysis/output/figures/L2D_SE_Morphology_Layer_boxplot.png",  layer_pal)

# ---- 按 Raw Material 分析 ----
cat("\n==== SE x Raw Material ====\n")
res_se_dir_rawmat   <- run_kw_dunn_se(se_df_base, "SE_direction",  "raw_material", "SE（方向）× Raw Material")
res_se_morph_rawmat <- run_kw_dunn_se(se_df_base, "SE_morphology", "raw_material", "SE（形态）× Raw Material")

p_se_dir_rawmat_box   <- make_se_boxplot(res_se_dir_rawmat,   "Spectral Entropy (Scar Direction)", "Raw Material",
                                         "analysis/output/figures/L2D_SE_Direction_RawMat_boxplot.png",   rawmat_pal)
p_se_morph_rawmat_box <- make_se_boxplot(res_se_morph_rawmat, "Spectral Entropy (Morphology)",     "Raw Material",
                                         "analysis/output/figures/L2D_SE_Morphology_RawMat_boxplot.png",  rawmat_pal)

# ---- 按 Core Type 分析 ----
cat("\n==== SE x Core Type ====\n")
res_se_dir_coretype   <- run_kw_dunn_se(se_df_base, "SE_direction",  "core_type", "SE（方向）× Core Type")
res_se_morph_coretype <- run_kw_dunn_se(se_df_base, "SE_morphology", "core_type", "SE（形态）× Core Type")

p_se_dir_coretype_box   <- make_se_boxplot(res_se_dir_coretype,   "Spectral Entropy (Scar Direction)", "Core Type",
                                           "analysis/output/figures/L2D_SE_Direction_CoreType_boxplot.png",   coretype_pal)
p_se_morph_coretype_box <- make_se_boxplot(res_se_morph_coretype, "Spectral Entropy (Morphology)",     "Core Type",
                                           "analysis/output/figures/L2D_SE_Morphology_CoreType_boxplot.png",  coretype_pal)

# ---- 汇总描述统计与 Dunn 结果 ----
se_desc_all <- bind_rows(
  se_desc_stats(se_df_base, "SE_direction",  "layer"),
  se_desc_stats(se_df_base, "SE_morphology", "layer"),
  se_desc_stats(se_df_base, "SE_direction",  "raw_material"),
  se_desc_stats(se_df_base, "SE_morphology", "raw_material"),
  se_desc_stats(se_df_base, "SE_direction",  "core_type"),
  se_desc_stats(se_df_base, "SE_morphology", "core_type")
)
write_csv(se_desc_all, here("analysis/data/derived_data/L2D_SE_desc_stats.csv"))
cat("已保存：L2D_SE_desc_stats.csv\n")

se_dunn_all <- bind_rows(
  if (!is.null(res_se_dir_layer))      res_se_dir_layer$dunn      %>% mutate(source = "SE_direction",  group_var = "layer"),
  if (!is.null(res_se_morph_layer))    res_se_morph_layer$dunn    %>% mutate(source = "SE_morphology", group_var = "layer"),
  if (!is.null(res_se_dir_rawmat))     res_se_dir_rawmat$dunn     %>% mutate(source = "SE_direction",  group_var = "raw_material"),
  if (!is.null(res_se_morph_rawmat))   res_se_morph_rawmat$dunn   %>% mutate(source = "SE_morphology", group_var = "raw_material"),
  if (!is.null(res_se_dir_coretype))   res_se_dir_coretype$dunn   %>% mutate(source = "SE_direction",  group_var = "core_type"),
  if (!is.null(res_se_morph_coretype)) res_se_morph_coretype$dunn %>% mutate(source = "SE_morphology", group_var = "core_type")
)
write_csv(se_dunn_all, here("analysis/data/derived_data/L2D_SE_dunn_results.csv"))
cat("已保存：L2D_SE_dunn_results.csv\n")


# ==============================================================================
# ---- 保存 L2 衍生数据 ----
# ==============================================================================

circ_desc_all <- bind_rows(
  if (!is.null(res_circ_layer))    res_circ_layer$desc,
  if (!is.null(res_circ_rawmat))   res_circ_rawmat$desc,
  if (!is.null(res_circ_coretype)) res_circ_coretype$desc
)
rayleigh_all <- bind_rows(
  if (!is.null(res_circ_layer))    res_circ_layer$rayleigh,
  if (!is.null(res_circ_rawmat))   res_circ_rawmat$rayleigh,
  if (!is.null(res_circ_coretype)) res_circ_coretype$rayleigh
)
if (nrow(circ_desc_all) > 0 && nrow(rayleigh_all) > 0) {
  left_join(circ_desc_all, rayleigh_all, by = c("group_var", "group")) %>%
    write_csv(here("analysis/data/derived_data/L2_circular_stats.csv"))
  cat("已保存：L2_circular_stats.csv\n")
}

scores_combined %>%
  select(ID, layer, raw_material, core_type,
         arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/L2_arrow_stats.csv"))
cat("已保存：L2_arrow_stats.csv\n")

scores_combined %>%
  write_csv(here("analysis/data/derived_data/CoIA_scores_full.csv"))
cat("已保存：CoIA_scores_full.csv\n")

morph_ilr_arch %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/SDG_morph_ILR_scores.csv"))
cat("已保存：SDG_morph_ILR_scores.csv\n")

scar_ilr_arch %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/SDG_scar_ILR_scores.csv"))
cat("已保存：SDG_scar_ILR_scores.csv\n")


# ==============================================================================
# ---- 组合图：CoIA × 箭头长度 × 方向（3行 × 3列）+ 外部图片行 ----
# ==============================================================================
cat("\n==== 组合图：CoIA / 箭头长度 / 方向 ====\n")

make_rose_for_composite <- function(res_circ) {
  if (is.null(res_circ)) return(NULL)
  plot_rose(
    list(group_col = res_circ$group_col,
         kde_df    = res_circ$kde_df,
         mean_dirs = res_circ$mean_dirs,
         rayleigh  = res_circ$rayleigh),
    res_circ$palette,
    show_color_legend = FALSE
  ) + theme(plot.margin = margin(4, 4, 4, 4))
}

strip_margin <- function(p) p + theme(plot.margin = margin(4, 4, 4, 4))

get_len_plot <- function(res_len) {
  if (is.null(res_len)) return(NULL)
  res_len$p + theme(plot.margin = margin(4, 4, 4, 4))
}

make_composite_row <- function(p_coia, p_len, res_circ, row_tag) {
  p_rose <- make_rose_for_composite(res_circ)
  if (is.null(p_coia) || is.null(p_len) || is.null(p_rose)) {
    cat(sprintf("  [跳过] %s 行：子图不完整\n", row_tag))
    return(NULL)
  }
  (strip_margin(p_coia) | get_len_plot(p_len) | p_rose) +
    plot_layout(widths = c(5, 2, 2))
}

row_layer    <- make_composite_row(p_coia_layer,    res_len_layer,    res_circ_layer,    "Layer")
row_rawmat   <- make_composite_row(p_coia_rawmat,   res_len_rawmat,   res_circ_rawmat,   "Raw Material")
row_coretype <- make_composite_row(p_coia_coretype, res_len_coretype, res_circ_coretype, "Core Type")

rows_valid <- Filter(Negate(is.null), list(row_layer, row_rawmat, row_coretype))

if (length(rows_valid) > 0) {
  n_rows     <- length(rows_valid)
  n_subplots <- n_rows * 3  # 每行3列
  
  # ---- 第一步：主体组合图（不加 plot_annotation，手动标签见下） ----
  # 先给每个行内子图手动打 A-I 标签
  all_tags <- LETTERS[seq_len(n_subplots)]
  tag_idx  <- 1L
  
  rows_tagged <- lapply(rows_valid, function(row) {
    # row 是 patchwork 对象，无法直接按格子访问，
    # 所以在 make_composite_row 里重建时打标签更可靠。
    # 此处用 wrap_elements 冻结整行，标签已在子图层打好（见下方改法）
    row
  })
  
  # 更稳妥：在 make_composite_row 返回前，对三个子图分别手动加 tag
  # 重新构建带标签的行
  make_tagged_row <- function(p_coia, p_len, res_circ, tags) {
    p_rose <- make_rose_for_composite(res_circ)
    if (is.null(p_coia) || is.null(p_len) || is.null(p_rose)) return(NULL)
    
    p1 <- strip_margin(p_coia)    + labs(tag = tags[1]) + theme(plot.tag = element_text(size = 11, face = "bold"))
    p2 <- get_len_plot(p_len)     + labs(tag = tags[2]) + theme(plot.tag = element_text(size = 11, face = "bold"))
    p3 <- p_rose                  + labs(tag = tags[3]) + theme(plot.tag = element_text(size = 11, face = "bold"))
    
    (p1 | p2 | p3) + plot_layout(widths = c(5, 2, 2))
  }
  
  # 重建带标签的各行
  inputs <- list(
    list(p_coia_layer,    res_len_layer,    res_circ_layer),
    list(p_coia_rawmat,   res_len_rawmat,   res_circ_rawmat),
    list(p_coia_coretype, res_len_coretype, res_circ_coretype)
  )
  valid_mask <- !sapply(rows_valid, is.null)  # rows_valid 已过滤，长度即 n_rows
  
  tagged_rows <- vector("list", n_rows)
  for (i in seq_len(n_rows)) {
    tags_i <- LETTERS[((i - 1) * 3 + 1):(i * 3)]
    inp    <- inputs[[i]]
    tagged_rows[[i]] <- make_tagged_row(inp[[1]], inp[[2]], inp[[3]], tags_i)
  }
  
  # ---- 第二步：主体拼图（冻结，防止加外部图时布局被破坏）----
  p_main <- Reduce(`/`, tagged_rows) +
    plot_layout(heights = rep(1, n_rows))
  
  # ---- 第三步：外部图片，手动加最后一个标签 ----
  next_tag    <- LETTERS[n_subplots + 1]
  external_img <- png::readPNG(here("asset/Axis_trajectory_SDG.png"))
  grob_img     <- grid::rasterGrob(external_img, interpolate = TRUE)
  p_external   <- wrap_elements(full = grob_img) +
    labs(tag = next_tag) +
    theme(
      plot.tag    = element_text(size = 11, face = "bold"),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # ---- 第四步：最终拼图 ----
  p_final <- wrap_elements(full = p_main) / p_external +
    plot_layout(heights = c(n_rows, 1))
  
  ggsave(
    here("analysis/output/figures/L_CoIA_composite.png"),
    plot   = p_final,
    width  = 12,
    height = n_rows * 5 + 3,  # 额外3英寸给外部图片行
    dpi    = 300,
    bg     = "white"
  )
  cat(sprintf("组合图已保存：L_CoIA_composite.png（%d 行 × 3 列 + 外部图片行）\n", n_rows))
  
} else {
  cat("  [跳过] 所有行均无有效子图，未生成组合图\n")
}
# ==============================================================================
# ========== 第三层：PERMANOVA ==========
# ==============================================================================

cat("\n\n")
cat("##  第三层：PERMANOVA — 分组结构分析                          ##\n")

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
  tibble(domain   = domain_name,
         grouping = group_name,
         R2       = res$R2[1],
         F_value  = res$F[1],
         p_value  = res$`Pr(>F)`[1],
         df_group = res$Df[1],
         df_resid = res$Df[2])
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
    run_permanova(D_morph_layer,    meta_layer$layer,         "Layer",        "Morphology"),
    run_permanova(D_scar_layer,     meta_layer$layer,         "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    run_permanova(D_morph_rawmat,   meta_rawmat$raw_material, "Raw Material", "Morphology"),
    run_permanova(D_scar_rawmat,    meta_rawmat$raw_material, "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    run_permanova(D_morph_coretype, meta_coretype$core_type,  "Core Type",    "Morphology"),
    run_permanova(D_scar_coretype,  meta_coretype$core_type,  "Core Type",    "Scar Direction")
  )
)

cat("\n==== L3-2：PERMDISP ====\n")

permdisp_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    run_permdisp(D_morph_layer,    meta_layer$layer,         "Layer",        "Morphology"),
    run_permdisp(D_scar_layer,     meta_layer$layer,         "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    run_permdisp(D_morph_rawmat,   meta_rawmat$raw_material, "Raw Material", "Morphology"),
    run_permdisp(D_scar_rawmat,    meta_rawmat$raw_material, "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    run_permdisp(D_morph_coretype, meta_coretype$core_type,  "Core Type",    "Morphology"),
    run_permdisp(D_scar_coretype,  meta_coretype$core_type,  "Core Type",    "Scar Direction")
  )
)

cat("\n==== L3-3：两两 PERMANOVA ====\n")

pairwise_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    pairwise_permanova(D_morph_layer,    meta_layer$layer,         "Layer",        "Morphology"),
    pairwise_permanova(D_scar_layer,     meta_layer$layer,         "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    pairwise_permanova(D_morph_rawmat,   meta_rawmat$raw_material, "Raw Material", "Morphology"),
    pairwise_permanova(D_scar_rawmat,    meta_rawmat$raw_material, "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    pairwise_permanova(D_morph_coretype, meta_coretype$core_type,  "Core Type",    "Morphology"),
    pairwise_permanova(D_scar_coretype,  meta_coretype$core_type,  "Core Type",    "Scar Direction")
  )
) %>%
  mutate(
    p_holm = p.adjust(p_raw, method = "holm"),
    significance = case_when(
      p_holm < 0.001 ~ "***", p_holm < 0.01 ~ "**",
      p_holm < 0.05  ~ "*",   TRUE           ~ "ns"
    )
  )

cat("\n两两 PERMANOVA（Holm 校正）：\n")
pairwise_results %>%
  mutate(across(c(R2, F_value, p_raw, p_holm), ~ round(.x, 4))) %>%
  print(n = Inf)


# ==============================================================================
# ---- 保存数值结果 ----
# ==============================================================================

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
           statistic = F_value, p_value = p_holm) %>%
    select(test, domain, grouping, R2, statistic, p_value)
) %>%
  write_csv(here("analysis/data/derived_data/L3_permanova.csv"))
cat("已保存：L3_permanova.csv\n")

# ==============================================================================
# ---- 组合图：SE × 分组（形态谱 / 方向谱）----
# ==============================================================================
cat("\n==== 组合图：SE 分组箱线图 ====\n")

add_tag <- function(p, tag) {
  p + labs(tag = tag) + theme(plot.tag = element_text(size = 11, face = "bold"))
}

row_morph <- (add_tag(p_se_morph_layer_box, "A") |
                add_tag(p_se_morph_rawmat_box, "B") |
                add_tag(p_se_morph_coretype_box, "C")) +
  plot_layout(widths = c(1, 1, 3))

row_dir <- (add_tag(p_se_dir_layer_box, "D") |
              add_tag(p_se_dir_rawmat_box, "E") |
              add_tag(p_se_dir_coretype_box, "F")) +
  plot_layout(widths = c(1, 1, 3))

p_se_composite <- (row_morph / row_dir) +
  plot_layout(heights = c(1, 1))

ggsave(
  here("analysis/output/figures/L2D_SE_composite.png"),
  plot   = p_se_composite,
  width  = 12,
  height = 10,
  dpi    = 300,
  bg     = "white"
)
cat("组合图已保存：L2D_SE_composite.png\n")

# ==============================================================================
# ---- 汇总打印 ----
# ==============================================================================

cat("\n\n")
cat("##  三层分析结果汇总（SDG）                                   ##\n")

cat("\n【第一层：基线】\n")
cat(sprintf("  Mantel r = %.4f, p = %.3f  ->  %s\n",
            mantel_global$statistic, mantel_global$signif,
            ifelse(mantel_global$signif < 0.05, "显著相关", "独立（ns）")))
cat(sprintf("  RV       = %.4f, p = %.3f  ->  %s\n",
            coin_arch$RV, rv_test$pvalue,
            ifelse(rv_test$pvalue < 0.05, "显著协变", "独立（ns）")))

cat("\n【第二层 A：分组 Mantel（Holm 校正）】\n")
if (nrow(l2_mantel) > 0) {
  n_sig <- sum(l2_mantel$p_holm < 0.05, na.rm = TRUE)
  cat(sprintf("  共检验 %d 个子组，Holm 校正后显著：%d 个\n",
              nrow(l2_mantel), n_sig))
  if (n_sig == 0) {
    cat("  -> 所有子组均不显著，形态-技术独立性稳健\n")
  } else {
    cat("  -> 以下子组显著（潜在新发现）：\n")
    l2_mantel %>% filter(p_holm < 0.05) %>%
      select(group_var, group, n, mantel_r, p_holm, significance) %>%
      print()
  }
}

cat("\n【第二层 B】参见 L2_Arrow_Length_*.png 及描述统计\n")
cat("【第二层 C】参见 L2_Arrow_Direction_rose_*.png 及 L2_circular_stats.csv\n")

cat("\n【第二层 D：spectral_entropy 分组分析（Holm 校正）】\n")
for (res_obj in list(res_se_dir_layer, res_se_morph_layer,
                     res_se_dir_rawmat, res_se_morph_rawmat,
                     res_se_dir_coretype, res_se_morph_coretype)) {
  if (!is.null(res_obj)) {
    cat(sprintf("  %s x %s: chi^2 = %.4f, p = %.4f -> %s\n",
                res_obj$y_col, res_obj$group_col,
                res_obj$kw$statistic, res_obj$kw$p.value,
                ifelse(res_obj$kw$p.value < 0.05, "显著", "ns")))
  }
}

cat("\n【第三层：PERMANOVA】\n")
permanova_results %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
      p_value < 0.05  ~ "*",   TRUE            ~ "ns"
    )
  ) %>%
  select(domain, grouping, R2, F_value, p_value, sig) %>%
  mutate(across(c(R2, F_value, p_value), ~ round(.x, 4))) %>%
  arrange(domain, grouping) %>%
  print()

cat("\n【CoIA 双标图】L1_CoIA_Biplot_layer.png / raw_material.png / core_type.png\n")
cat("【组合图】L_CoIA_composite.png（CoIA | 箭头长度 | 方向，3行 × 3列）\n")
cat("【桑基图】L1_CoIA_Sankey.png\n")

cat("\n\n========== SDG 分析全部完成 ==========\n")
cat("主要输出文件：\n")
cat("  L1_results.csv\n")
cat("  L1_CoIA_Sankey.png\n")
cat("  L_CoIA_composite.png\n")
cat("  L2_grouped_mantel.csv\n")
cat("  L2_arrow_stats.csv\n")
cat("  L2_circular_stats.csv\n")
cat("  L2D_SE_desc_stats.csv\n")
cat("  L2D_SE_dunn_results.csv\n")
cat("  L3_permanova.csv\n")
cat("  CoIA_scores_full.csv\n")
cat("  CoIA_coords_full.csv\n")
cat("  PCA_CoIA_contribution.csv\n")
cat("  SDG_morph_ILR_scores.csv\n")
cat("  SDG_scar_ILR_scores.csv\n")