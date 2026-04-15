# ==============================================================================
# SDG_cores_statistics.R
#
# 第一层：整体 Mantel + CIA
#   形态与片疤方向是否整体独立？
#   附：CIA 箭头长度 + 箭头方位描述
#   附：CIA 桑基图（ILR -> PCA -> CoIA 贡献流）
#
# 第二层：分组 Mantel（层位 × 原料 × 类型）
#   这种独立性是否在所有条件下都成立？
#   联合证据：
#     - 分组 Mantel r
#     - 箭头长度分组差异
#     - 箭头方位分组差异
#     - spectral_entropy × Typology / Layer / Raw_mat
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
#   L1_CIA_Biplot.png
#   L1_CIA_Diagnostics.png
#   L1_CIA_Sankey.png
#   L2_Mantel_Grouped_dumbbell.png
#   L2_Mantel_Grouped_heatmap.png
#   L2_Arrow_Length_layer.png
#   L2_Arrow_Length_raw_material.png       （有效分组时）
#   L2_Arrow_Length_core_type.png          （有效分组时）
#   L2_Arrow_Direction_rose_layer.png
#   L2_Arrow_Direction_rose_raw_material.png   （有效分组时）
#   L2_Arrow_Direction_rose_core_type.png      （有效分组时）
#   L2D_SE_Direction_Layer_boxplot.png
#   L2D_SE_Morphology_Layer_boxplot.png
#   L2D_SE_Direction_RawMat_boxplot.png
#   L2D_SE_Morphology_RawMat_boxplot.png
#   L2D_SE_Direction_CoreType_boxplot.png
#   L2D_SE_Morphology_CoreType_boxplot.png
#   L3_PERMANOVA_R2.png
#   L3_PERMANOVA_CAP_layer.png
#   L3_PERMANOVA_CAP_raw_material.png
#   L3_PERMANOVA_CAP_core_type.png
#
# 输出（derived_data/）：
#   L1_results.csv
#   L2_grouped_mantel.csv
#   L2_arrow_stats.csv
#   L2_circular_stats.csv
#   L2D_SE_desc_stats.csv
#   L2D_SE_dunn_results.csv
#   L3_permanova.csv
#   CIA_scores_full.csv
#   CIA_coords_full.csv
#   PCA_CoIA_contribution.csv
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
# ---- 全局辅助函数 ----
# ==============================================================================

# 【修改1】移除 cosine_dist，统一使用 ILR + 欧氏距离

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

POWER_COLS <- paste0("power_l", 1:5)

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

# 提取 power 特征矩阵
morph_power <- df_morph_raw %>%
  select(all_of(POWER_COLS)) %>%
  rename_with(~ paste0("M", seq_along(.))) %>%
  as.data.frame()
scar_power <- df_scar_raw %>%
  select(all_of(POWER_COLS)) %>%
  rename_with(~ paste0("S", seq_along(.))) %>%
  as.data.frame()
rownames(morph_power) <- df_morph_raw$ID
rownames(scar_power)  <- df_scar_raw$ID

morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

# 【修改1】ILR 变换 -> 欧氏距离（与 EXP 脚本对齐，替换原余弦距离）
morph_ilr_all <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power_clean))))
scar_ilr_all  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power_clean))))
rownames(morph_ilr_all) <- rownames(morph_power_clean)
rownames(scar_ilr_all)  <- rownames(scar_power_clean)

D_morph_all <- dist(morph_ilr_all)
D_scar_all  <- dist(scar_ilr_all)

# 分离考古标本（去除 IM_ 参照件）
arch_ids <- rownames(morph_power_clean)[!str_starts(rownames(morph_power_clean), "IM_")]
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

# 调色板
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

# 【修改4】多重检验校正：Holm（原脚本为 FDR）
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
  scale_fill_viridis_c(option = "D", limits = c(-1, 1),
                       name = "Spearman's rho") +
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

# PCA 详细报告（与 EXP 对齐）
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
            by = "ID")

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
          here("analysis/data/derived_data/CIA_coords_full.csv"))
cat("已保存：CIA_coords_full.csv\n")

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
# CoIA 辅助可视化：诊断图（与 EXP 对齐）
# ------------------------------------------------------------------------------

eig_df <- tibble(
  axis       = paste0("Axis ", seq_along(coin_arch$eig)),
  eigenvalue = coin_arch$eig,
  pct        = coin_arch$eig / sum(coin_arch$eig) * 100,
  cum_pct    = cumsum(pct)
)

p_scree <- ggplot(eig_df, aes(x = axis, y = pct)) +
  geom_col(fill = "#7EB8C9", alpha = 0.85, width = 0.55) +
  geom_line(aes(y = cum_pct, group = 1), color = "#6271A1", linewidth = 0.8) +
  geom_point(aes(y = cum_pct), color = "#6271A1", size = 2.8) +
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
                                  "Morphology PCA axes on CoIA space", "#7EB8C9")
p_load_scar  <- make_loading_plot(scar_load,
                                  "Scar direction PCA axes on CoIA space", "#E6B89C")

p_cia_diagnostics <- (p_scree | p_load_morph | p_load_scar) +
  plot_annotation(
    title   = "CoIA Axis Diagnostics (SDG)",
    caption = "Left: scree plot. Middle: morphology PCA axis loadings on CoIA axes. Right: scar-direction PCA axis loadings on CoIA axes.\nBoth panels use PCA-axis projections (aX / aY) for symmetric interpretation.",
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0)
    )
  )

ggsave(here("analysis/output/figures/L1_CIA_Diagnostics.png"),
       plot = p_cia_diagnostics, width = 15, height = 5.5, dpi = 300, bg = "white")
cat("图已保存：L1_CIA_Diagnostics.png\n")


# CoIA 主双标图（按 layer 着色，保留原有逻辑）
p_cia_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_segment(
    data = scores_combined %>% filter(!is.na(layer), layer != "Other"),
    aes(x = Axis1_M, y = Axis2_M, xend = Axis1_S, yend = Axis2_S,
        color = layer),
    linewidth = 0.45, alpha = 0.45, lineend = "round"
  ) +
  geom_point(
    data = scores_combined %>% filter(!is.na(layer), layer != "Other"),
    aes(x = Axis1_M, y = Axis2_M, fill = layer, color = layer),
    shape = 21, size = 3.0, stroke = 0.5, alpha = 0.90
  ) +
  geom_point(
    data = scores_combined %>% filter(!is.na(layer), layer != "Other"),
    aes(x = Axis1_S, y = Axis2_S, fill = layer, color = layer),
    shape = 24, size = 2.6, stroke = 0.5, alpha = 0.90
  ) +
  geom_text_repel(
    data = scores_combined %>%
      filter(!is.na(layer), layer != "Other") %>%
      slice_max(arrow_length, n = 5),
    aes(x = (Axis1_M + Axis1_S) / 2, y = (Axis2_M + Axis2_S) / 2,
        label = ID, color = layer),
    size = 2.2, show.legend = FALSE,
    max.overlaps = 20, segment.color = "grey60", segment.linewidth = 0.3
  ) +
  scale_color_manual(values = layer_pal, name = "Layer") +
  scale_fill_manual(values  = layer_pal, name = "Layer") +
  theme_bw(base_size = 10) +
  labs(
    title = sprintf("CoIA Biplot (SDG)  |  RV = %.3f, p = %.3f",
                    coin_arch$RV, rv_test$pvalue),
    subtitle = "Circle = Morphology endpoint; Triangle = Scar-direction endpoint; Arrow = decoupling",
    x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
    y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])
  ) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle   = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
        legend.position = "right")

ggsave(here("analysis/output/figures/L1_CIA_Biplot.png"),
       plot = p_cia_biplot, width = 10, height = 8, dpi = 300, bg = "white")
cat("图已保存：L1_CIA_Biplot.png\n")

l1_results <- tibble(
  method  = c("Mantel (ILR, Euclidean, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_arch$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(arch_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/L1_results.csv"))
cat("\n第一层结论：\n"); print(l1_results)


# ==============================================================================
# ---- L1-3：CoIA 桑基图（ILR -> PCA -> CoIA 贡献流）【新增，与 EXP 对齐】----
# ==============================================================================

cat("\n==== L1-3：CoIA 桑基图 ====\n")

# ---- 提取权重矩阵 ----
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

# 保留累计方差 <= 95% 的轴（至少保留第1轴）
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

# ---- 节点高度 ----
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

# ---- 画布参数 ----
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

# ---- SVG 辅助函数 ----
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

# ---- 组装 SVG ----
lines <- character(0)
push  <- function(...) lines <<- c(lines, ...)

push(sprintf(
  '<svg xmlns="http://www.w3.org/2000/svg" width="%d" viewBox="0 0 %d %d">',
  SVG_W, SVG_W, SVG_H
))
push('<title>ILR to PCA to CoIA contribution flow (SDG)</title>')
push('<desc>Sankey-style flow diagram showing ILR spectral variable contributions through PCA axes to CoIA axes.</desc>')
push(sprintf('<rect width="%d" height="%d" fill="white"/>', SVG_W, SVG_H))

push(svg_text(x_col1 + NODE_W / 2, 20, "ILR variables",   12, 500, fill = "#2C2C2A"))
push(svg_text(x_col1 + NODE_W / 2, 36, "spectral log-contrasts", 10, 400, fill = "#5F5E5A"))
push(svg_text(x_col2 + NODE_W / 2, 20, "PCA axes",        12, 500, fill = "#2C2C2A"))
push(svg_text(x_col2 + NODE_W / 2, 36, "per endpoint",    10, 400, fill = "#5F5E5A"))
push(svg_text(x_col3 + NODE_W / 2, 20, "CoIA axes",       12, 500, fill = "#2C2C2A"))
push(svg_text(x_col3 + NODE_W / 2, 36, "shared structure",10, 400, fill = "#5F5E5A"))

# 流带
outlet_ilr_m <- rep(0, n_ilr_m)
inlet_mpc    <- rep(0, n_mpc)
outlet_ilr_d <- rep(0, n_ilr_d)
inlet_dpc    <- rep(0, n_dpc)
inlet_cia    <- rep(0, n_cia_ax)

for (i in seq_len(n_ilr_m)) {
  for (j in seq_len(n_mpc)) {
    ww <- w2_ilr_mpc[i, j]
    if (ww < 0.004) next
    bw_out <- ww * h_ilr_m[i]; bw_in <- ww * h_mpc[j]
    y1t <- y_ilr_m[i] + outlet_ilr_m[i]; y1b <- y1t + bw_out
    y2t <- y_mpc[j]   + inlet_mpc[j];    y2b <- y2t + bw_in
    push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                  col_morph_band, 0.28 + 0.32 * ww))
    outlet_ilr_m[i] <- outlet_ilr_m[i] + bw_out
    inlet_mpc[j]    <- inlet_mpc[j]    + bw_in
  }
}

for (i in seq_len(n_ilr_d)) {
  for (j in seq_len(n_dpc)) {
    ww <- w2_ilr_dpc[i, j]
    if (ww < 0.004) next
    bw_out <- ww * h_ilr_d[i]; bw_in <- ww * h_dpc[j]
    y1t <- y_ilr_d[i] + outlet_ilr_d[i]; y1b <- y1t + bw_out
    y2t <- y_dpc[j]   + inlet_dpc[j];    y2b <- y2t + bw_in
    push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                  col_scar_band, 0.28 + 0.32 * ww))
    outlet_ilr_d[i] <- outlet_ilr_d[i] + bw_out
    inlet_dpc[j]    <- inlet_dpc[j]    + bw_in
  }
}

outlet_mpc <- rep(0, n_mpc)
for (j in seq_len(n_mpc)) {
  for (k in seq_len(n_cia_ax)) {
    ww <- w2_mpc_cia[j, k]
    if (ww < 0.004) next
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
}

outlet_dpc <- rep(0, n_dpc)
for (j in seq_len(n_dpc)) {
  for (k in seq_len(n_cia_ax)) {
    ww <- w2_dpc_cia[j, k]
    if (ww < 0.004) next
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
}

# 节点
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
cia_labels <- c("CoIA Axis 1", "CoIA Axis 2")
for (k in seq_len(n_cia_ax)) {
  push(svg_rect(x_col3, y_cia[k], NODE_W, h_cia[k], col_cia_fill, col_cia_stroke))
  push(svg_text(x_col3 + NODE_W / 2, y_cia[k] + h_cia[k] / 2,
                cia_labels[k], 11, 500, fill = col_cia_text))
}

push("</svg>")

# ---- PNG 转换 ----
svg_path <- tempfile(fileext = ".svg")
png_path <- here("analysis/output/figures/L1_CIA_Sankey.png")
writeLines(lines, svg_path, useBytes = FALSE)

if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_png(svg_path, png_path, width = SVG_W * 2)
  cat("PNG 已保存（via rsvg）：L1_CIA_Sankey.png\n")
} else if (nzchar(Sys.which("rsvg-convert"))) {
  system2("rsvg-convert",
          args = c("-d", "300", "-p", "300", "-o", png_path, svg_path))
  cat("PNG 已保存（via rsvg-convert）：L1_CIA_Sankey.png\n")
} else if (nzchar(Sys.which("inkscape"))) {
  system2("inkscape",
          args = c("--export-filename", png_path, "--export-dpi", "300", svg_path))
  cat("PNG 已保存（via Inkscape）：L1_CIA_Sankey.png\n")
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

# 【修改4】Holm 校正（原脚本为 FDR）
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
      name = "Significance\n(Holm corrected)"
    ) +
    scale_shape_manual(
      values = c("Layer" = 16, "Raw Material" = 17, "Core Type" = 15),
      name   = "Grouping variable"
    ) +
    theme_bw(base_size = 10) +
    labs(
      title    = "L2-A: Within-group Mantel Test (Morphology x Scar Direction)",
      subtitle = "Each point = Mantel r within one subgroup; dashed line = global baseline; p = Holm corrected",
      x = "Mantel r (Spearman, Holm corrected)", y = NULL
    ) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.subtitle   = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
          strip.text      = element_text(face = "bold", size = 9),
          legend.position = "right")
  
  ggsave(here("analysis/output/figures/L2_Mantel_Grouped_dumbbell.png"),
         plot = p_l2_dumbbell,
         width = 10, height = max(4, nrow(l2_mantel) * 0.7 + 2),
         dpi = 300, bg = "white")
  cat("图已保存：L2_Mantel_Grouped_dumbbell.png\n")
  
  heatmap_df <- l2_mantel %>%
    mutate(group = fct_reorder(group, mantel_r, .fun = mean),
           y_int = as.integer(group))
  y_labels <- heatmap_df %>% distinct(y_int, group) %>% arrange(y_int) %>%
    mutate(group_chr = as.character(group))
  
  p_l2_heatmap <- ggplot(heatmap_df,
                         aes(x = group_var_label, y = y_int, fill = mantel_r)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("r = %.3f\n%s", mantel_r, significance)),
              size = 2.8, color = "grey20") +
    scale_fill_gradient2(
      low = "#3B82C4", mid = "white", high = "#C0392B",
      midpoint = 0, limits = c(-0.3, 0.3),
      oob = scales::squish, name = "Mantel r"
    ) +
    scale_y_continuous(breaks = y_labels$y_int, labels = y_labels$group_chr,
                       expand = expansion(add = 0.5)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    theme_bw(base_size = 10) +
    labs(
      title    = "L2-A: Within-group Mantel r (Morphology x Scar Direction)",
      subtitle = sprintf("Global baseline: r = %.3f (p = %.3f) | p = Holm corrected",
                         mantel_global$statistic, mantel_global$signif),
      x = "Grouping variable", y = "Subgroup"
    ) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
          axis.text.y   = element_text(size = 8))
  
  ggsave(here("analysis/output/figures/L2_Mantel_Grouped_heatmap.png"),
         plot   = p_l2_heatmap,
         width  = max(6, length(unique(heatmap_df$group_var_label)) * 3 + 2),
         height = max(4, nrow(l2_mantel) * 0.6 + 2),
         dpi = 300, bg = "white")
  cat("图已保存：L2_Mantel_Grouped_heatmap.png\n")
  
  write_csv(l2_mantel, here("analysis/data/derived_data/L2_grouped_mantel.csv"))
}


# ------------------------------------------------------------------------------
# L2-B：CoIA 箭头长度
# ------------------------------------------------------------------------------

cat("\n---------- L2-B：箭头长度分组差异 ----------\n")

run_arrow_length_analysis <- function(group_col, group_label, palette) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 3) %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [跳过] %s 箭头长度检验：有效分组不足\n", group_label))
    return(invisible(NULL))
  }
  sub_df <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups)
  cat(sprintf("\n----- %s x 箭头长度 -----\n", group_label))
  kw <- kruskal.test(reformulate(group_col, "arrow_length"), data = sub_df)
  print(kw)
  # 【修改4】Holm 校正
  pw <- pairwise.wilcox.test(sub_df$arrow_length, sub_df[[group_col]],
                             p.adjust.method = "holm", exact = FALSE)
  print(pw)
  p <- ggplot(sub_df, aes(x = .data[[group_col]], y = arrow_length,
                          fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.5) +
    geom_jitter(aes(color = .data[[group_col]]), width = 0.15,
                size = 2.5, alpha = 0.75, show.legend = FALSE) +
    geom_text_repel(
      data = sub_df %>% group_by(.data[[group_col]]) %>% slice_max(arrow_length, n = 1),
      aes(label = ID), size = 2.4, color = "grey40",
      max.overlaps = 10, show.legend = FALSE
    ) +
    annotate("text", x = 1.5, y = max(sub_df$arrow_length) * 1.02,
             label = sprintf("Kruskal-Wallis: chi^2 = %.2f, p = %.3f",
                             kw$statistic, kw$p.value),
             size = 3, color = "grey30", hjust = 0.5) +
    scale_fill_manual(values  = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    theme_bw(base_size = 10) +
    labs(title    = sprintf("L2-B: CoIA Arrow Length by %s", group_label),
         subtitle = "Longer arrows = greater morphology-technique decoupling",
         x = NULL, y = "Arrow length (CoIA space)") +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
          axis.text.x   = element_text(angle = 15, hjust = 1))
  fname <- sprintf("analysis/output/figures/L2_Arrow_Length_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  list(sub_df = sub_df, kw = kw, pw = pw)
}

res_len_layer    <- run_arrow_length_analysis("layer",        "Layer",       layer_pal)
res_len_rawmat   <- run_arrow_length_analysis("raw_material", "Raw_Material",rawmat_pal)
res_len_coretype <- run_arrow_length_analysis("core_type",    "Core_Type",   coretype_pal)

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
# L2-C：CoIA 箭头方位圆形统计
# ------------------------------------------------------------------------------

cat("\n---------- L2-C：箭头方位圆形统计 ----------\n")

run_circular_analysis <- function(group_col, group_label, palette) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 5) %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [跳过] %s 圆形统计：有效组数不足\n", group_label))
    return(invisible(NULL))
  }
  sub_df <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups)
  
  cat(sprintf("\n----- %s 圆形描述统计 -----\n", group_label))
  circ_desc <- map_dfr(valid_groups, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(group_var = group_col, group = g, n = length(angles),
           mean_dir_deg = round(cs$mean_deg, 2),
           concentration_r = round(cs$rho, 4))
  })
  print(circ_desc)
  
  cat(sprintf("\n----- %s Rayleigh 检验 -----\n", group_label))
  rayleigh_res <- map_dfr(valid_groups, function(g) {
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
  } else {
    watson_res <- NULL
  }
  
  mean_dirs <- map_dfr(valid_groups, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs <- circ_stats_one(angles)
    tibble(!!group_col := g, mean_deg = cs$mean_deg)
  })
  
  rose_df <- sub_df %>%
    mutate(angle_deg = arrow_angle * 180 / pi,
           angle_deg = ifelse(angle_deg < 0, angle_deg + 360, angle_deg))
  
  compute_circular_kde <- function(angles_deg, bw = 25, n = 360) {
    circ <- circular(angles_deg * pi / 180,
                     type = "angles", units = "radians", modulo = "2pi")
    dens <- density(circ, bw = bw, n = n)
    angle_deg <- as.numeric(dens$x) * 180 / pi %% 360
    density   <- as.numeric(dens$y)
    tibble(angle_deg = c(angle_deg, angle_deg[1]),
           density   = c(density,   density[1]))
  }
  
  kde_df <- rose_df %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~ compute_circular_kde(.x$angle_deg, bw = 25)) %>%
    ungroup()
  
  p_rose <- ggplot(kde_df,
                   aes(x = angle_deg, y = density,
                       fill = .data[[group_col]],
                       color = .data[[group_col]])) +
    geom_col(width = 0.1, alpha = 0.4, position = "identity") +
    geom_segment(data = mean_dirs,
                 aes(x = mean_deg, xend = mean_deg,
                     y = 0, yend = Inf,
                     color = .data[[group_col]]),
                 linewidth = 0.9, linetype = "dashed", alpha = 0.85) +
    coord_polar(theta = "x", start = -pi / 2, direction = 1) +
    scale_x_continuous(
      limits = c(0, 360),
      breaks = seq(0, 315, by = 45),
      labels = c("0\n(+CoIA1)", "45", "90\n(+CoIA2)",
                 "135", "180\n(-CoIA1)", "225", "270\n(-CoIA2)", "315")
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values  = palette, name = group_label) +
    scale_color_manual(values = palette, guide = "none") +
    facet_wrap(reformulate(group_col), ncol = min(length(valid_groups), 3)) +
    theme_bw(base_size = 10) +
    labs(title    = sprintf("L2-C: CoIA Direction Density by %s", group_label),
         subtitle = "von Mises kernel density (circular) with mean direction",
         x = NULL, y = NULL) +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.subtitle    = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
          axis.text.y      = element_blank(),
          axis.ticks.y     = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold", size = 9))
  
  n_g   <- length(valid_groups)
  fname <- sprintf("analysis/output/figures/L2_Arrow_Direction_rose_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  ggsave(here(fname), plot = p_rose,
         width = min(4 + n_g * 3, 16), height = 8, dpi = 600, bg = "white")
  cat(sprintf("图已保存：%s\n", basename(fname)))
  list(desc = circ_desc, rayleigh = rayleigh_res, watson = watson_res)
}

res_circ_layer    <- run_circular_analysis("layer",        "Layer",       layer_pal)
res_circ_rawmat   <- run_circular_analysis("raw_material", "Raw_Material",rawmat_pal)
res_circ_coretype <- run_circular_analysis("core_type",    "Core_Type",   coretype_pal)


# ------------------------------------------------------------------------------
# L2-D：spectral_entropy × Typology / Layer / Raw_mat【新增，与 EXP 对齐】
# ------------------------------------------------------------------------------

cat("\n---------- L2-D：spectral_entropy 分组分析 ----------\n")

se_df_base <- meta_arch %>%
  filter(ID %in% arch_ids) %>%
  select(ID, layer, raw_material, core_type,
         SE_direction, SE_morphology)

# KW + 事后检验通用函数（Holm 校正）
# 2组时用 Wilcoxon 替代 dunnTest（dunnTest 在 k=2 时内部排序会报错）
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
    # 只有2组：直接用 Wilcoxon，构造与 Dunn 输出相同结构的 tibble
    cat("  [注] 仅2组，使用 Wilcoxon 检验替代 Dunn 检验\n")
    grp_levels <- levels(sub_df[[group_col]])
    wt <- wilcox.test(
      sub_df[[y_col]][sub_df[[group_col]] == grp_levels[1]],
      sub_df[[y_col]][sub_df[[group_col]] == grp_levels[2]],
      exact = FALSE
    )
    p_raw <- wt$p.value
    p_adj <- p_raw   # 只有1对，Holm 校正 = 原始 p
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
                wt$statistic, p_raw,
                dunn$p.signif))
  } else {
    # 3组及以上：使用 dunnTest（Holm 校正）
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

# SE 描述统计输出函数
se_desc_stats <- function(df, y_col, group_col) {
  df %>%
    filter(!is.na(.data[[y_col]]), !is.na(.data[[group_col]]),
           .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>%
    summarise(n            = n(),
              mean_val     = round(mean(.data[[y_col]],   na.rm = TRUE), 4),
              median_val   = round(median(.data[[y_col]], na.rm = TRUE), 4),
              sd_val       = round(sd(.data[[y_col]],     na.rm = TRUE), 4),
              .groups = "drop") %>%
    mutate(SE_type    = y_col,
           group_type = group_col)
}

# 箱线图绘制函数（与 EXP make_se_boxplot 对齐）
make_se_boxplot <- function(res_obj, y_label, title_suffix, fname,
                            palette) {
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
  
  sig_pairs   <- dunn_res %>% filter(p.adj < 0.05)
  trend_pairs <- dunn_res %>% filter(p < 0.05, p.adj >= 0.05)
  
  all_annot_base <- bind_rows(
    sig_pairs   %>% mutate(annot_type = "sig"),
    trend_pairs %>% mutate(annot_type = "trend")
  )
  sig_annot <- if (nrow(all_annot_base) > 0) {
    all_annot_base %>%
      mutate(x1    = match(group1, type_lvls),
             x2    = match(group2, type_lvls),
             y_bar = y_max + step * row_number(),
             x_mid = (x1 + x2) / 2,
             label = if_else(annot_type == "sig", p.adj.signif, "\u2020"))
  } else NULL
  
  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[y_col]],
                      fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, linewidth = 0.55, width = 0.55) +
    geom_jitter(width = 0.14, size = 2.2, alpha = 0.80, stroke = 0.3,
                shape = 21, color = "grey20", fill = NA, show.legend = FALSE) +
    stat_summary(
      fun.data = function(x) {
        ypos <- min(x, na.rm = TRUE) - y_range * 0.04
        data.frame(y = ypos, ymin = ypos, ymax = ypos,
                   label = paste0("n=", length(x)))
      },
      geom = "text", aes(label = after_stat(label)),
      size = 2.8, color = "grey45", vjust = 1, show.legend = FALSE
    ) +
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
          scale_linetype_manual(values = c("sig" = "solid", "trend" = "dashed"),
                                guide = "none")
        )
      } else NULL
    } +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Kruskal-Wallis\nchi^2 = %.2f, p = %.3f",
                             kw_res$statistic, kw_res$p.value),
             hjust = 1.05, vjust = 1.15, size = 3.0,
             color = ifelse(kw_res$p.value < 0.05, "#C0392B", "grey50")) +
    scale_fill_manual(values  = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.08))) +
    theme_bw(base_size = 11) +
    labs(title    = sprintf("L2-D: Spectral Entropy by %s — %s", title_suffix, y_label),
         subtitle = "Solid bracket: Holm p.adj < 0.05; dashed bracket: raw p < 0.05 trend",
         x = NULL, y = y_label) +
    theme(plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
          plot.subtitle      = element_text(hjust = 0.5, size = 8.5, color = "grey55"),
          axis.text.x        = element_text(angle = 25, hjust = 1, size = 9),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())
  
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

make_se_boxplot(res_se_dir_layer,   "Spectral Entropy (Scar Direction)", "Layer",
                "analysis/output/figures/L2D_SE_Direction_Layer_boxplot.png", layer_pal)
make_se_boxplot(res_se_morph_layer, "Spectral Entropy (Morphology)",     "Layer",
                "analysis/output/figures/L2D_SE_Morphology_Layer_boxplot.png", layer_pal)

# ---- 按 Raw Material 分析 ----
cat("\n==== SE x Raw Material ====\n")
res_se_dir_rawmat   <- run_kw_dunn_se(se_df_base, "SE_direction",  "raw_material", "SE（方向）× Raw Material")
res_se_morph_rawmat <- run_kw_dunn_se(se_df_base, "SE_morphology", "raw_material", "SE（形态）× Raw Material")

make_se_boxplot(res_se_dir_rawmat,   "Spectral Entropy (Scar Direction)", "Raw Material",
                "analysis/output/figures/L2D_SE_Direction_RawMat_boxplot.png", rawmat_pal)
make_se_boxplot(res_se_morph_rawmat, "Spectral Entropy (Morphology)",     "Raw Material",
                "analysis/output/figures/L2D_SE_Morphology_RawMat_boxplot.png", rawmat_pal)

# ---- 按 Core Type 分析 ----
cat("\n==== SE x Core Type ====\n")
res_se_dir_coretype   <- run_kw_dunn_se(se_df_base, "SE_direction",  "core_type", "SE（方向）× Core Type")
res_se_morph_coretype <- run_kw_dunn_se(se_df_base, "SE_morphology", "core_type", "SE（形态）× Core Type")

make_se_boxplot(res_se_dir_coretype,   "Spectral Entropy (Scar Direction)", "Core Type",
                "analysis/output/figures/L2D_SE_Direction_CoreType_boxplot.png", coretype_pal)
make_se_boxplot(res_se_morph_coretype, "Spectral Entropy (Morphology)",     "Core Type",
                "analysis/output/figures/L2D_SE_Morphology_CoreType_boxplot.png", coretype_pal)

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

scores_combined %>%
  select(ID, layer, raw_material, core_type,
         arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/L2_arrow_stats.csv"))
cat("已保存：L2_arrow_stats.csv\n")

scores_combined %>%
  write_csv(here("analysis/data/derived_data/CIA_scores_full.csv"))
cat("已保存：CIA_scores_full.csv\n")


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

# 【修改4】Holm 校正（原脚本为 FDR）
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

# R² 比较图
p_r2_compare <- permanova_results %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
      p_value < 0.05  ~ "*",   TRUE            ~ "ns"
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
    values = c("Morphology" = "#A1C2E6", "Scar Direction" = "#FFBAE0"),
    name   = "Domain"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.12))) +
  theme_bw(base_size = 10) +
  labs(title    = "L3: PERMANOVA — Proportion of Variance Explained (R²)",
       subtitle = "*** p<0.001  ** p<0.01  * p<0.05  ns = not significant",
       x = "Grouping factor", y = "R² (proportion of variance)") +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle   = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
        legend.position = "right")

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

cat("\n【桑基图】L1_CIA_Sankey.png\n")

cat("\n\n========== SDG 分析全部完成 ==========\n")
cat("主要输出文件：\n")
cat("  L1_results.csv\n")
cat("  L1_CIA_Sankey.png\n")
cat("  L2_grouped_mantel.csv\n")
cat("  L2_arrow_stats.csv\n")
cat("  L2_circular_stats.csv\n")
cat("  L2D_SE_desc_stats.csv\n")
cat("  L2D_SE_dunn_results.csv\n")
cat("  L3_permanova.csv\n")
cat("  CIA_scores_full.csv\n")
cat("  CIA_coords_full.csv\n")
cat("  PCA_CoIA_contribution.csv\n")