# ==============================================================================
# exp_cores_statistics_v2.R
#
# 主要改动：
#   1. 对功率谱做 ILR 变换后先 PCA，保留累计方差 >= 95% 的轴（至少1轴）
#   2. 用截断 PCA 坐标计算欧氏距离，替代原始 ILR 距离
#   3. PERMANOVA（adonis2）替代 Mantel 作为 L1 主检验（+ 保留 Mantel 供比较）
#   4. CoIA 使用截断 dudi.pca 对象（nf = k_morph / k_scar）
#   5. 分组 PERMANOVA 替代分组 Mantel（L2-A）
#   6. 所有统计输出保存至 derived_data/，无任何图形输出
#
# 输入：
#   analysis/data/derived_data/SPHARM_direction.csv
#   analysis/data/derived_data/SPHARM_morphology.csv
#   analysis/data/raw_data/SDG_core_metric.xlsx
#
# 输出（derived_data/）：
#   EXP_PCA_truncation_summary.csv      <- 各端截断信息
#   EXP_L1_results.csv                  <- Mantel + RV 基线
#   EXP_L1_PERMANOVA_results.csv        <- 全局 PERMANOVA
#   EXP_L2_grouped_permanova.csv        <- 分组 PERMANOVA（Typology）
#   EXP_L2_grouped_mantel.csv           <- 分组 Mantel（保留供比较）
#   EXP_L2_arrow_stats.csv              <- CoIA 箭头长度/方向
#   EXP_L2_circular_stats.csv           <- 圆形统计
#   EXP_CIA_scores_full.csv             <- CoIA 样本坐标
#   EXP_CIA_coords_full.csv             <- CoIA 坐标（含箭头）
#   EXP_PCA_CoIA_contribution.csv       <- PCA 轴对 CoIA 贡献
#   EXP_L2D_SE_desc_stats.csv           <- SE 描述统计
#   EXP_L2D_SE_dunn_results.csv         <- SE Dunn 检验
#   EXP_morph_ILR_scores.csv
#   EXP_scar_ILR_scores.csv
#   EXP_morph_PCA_scores_trunc.csv      <- 截断后形态 PCA 坐标
#   EXP_scar_PCA_scores_trunc.csv       <- 截断后方向 PCA 坐标
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(vegan)
library(ggplot2)
library(compositions)
library(ade4)
library(circular)
library(FSA)

# ==============================================================================
# ---- 全局常量 ----
# ==============================================================================

TYPOLOGY_COLORS <- c(
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Unidirectional" = "#BA8530",
  "Multiplatform"  = "#8A7A68",
  "Bidirectional"  = "#788C4A"
)

TYPOLOGY_ORDER <- c(
  "Unidirectional",
  "Bidirectional",
  "Levallois",
  "Discoid",
  "Multiplatform"
)

VAR_THRESHOLD <- 0.95   # 截断阈值

POWER_COLS_DIR   <- paste0("power_l", 1:6)
POWER_COLS_MORPH <- paste0("power_l", 1:8)


# ==============================================================================
# ---- 辅助函数 ----
# ==============================================================================

replace_zeros <- function(X, delta = NULL) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i    <- X[i, ]
    zero_idx <- row_i == 0
    if (!any(zero_idx)) next
    nonzero_min          <- min(row_i[!zero_idx])
    d                    <- ifelse(is.null(delta), nonzero_min * 0.65, delta)
    n_zero               <- sum(zero_idx)
    row_i[zero_idx]      <- d
    row_i[!zero_idx]     <- row_i[!zero_idx] * (1 - n_zero * d)
    X[i, ]               <- row_i
  }
  X
}

extract_subdist <- function(D_full, ids) {
  as.dist(as.matrix(D_full)[ids, ids])
}

# PCA 截断：返回保留轴数、截断坐标、截断 dudi 对象（nf 已设为 k）
truncate_pca <- function(dudi_obj, threshold = VAR_THRESHOLD, label = "") {
  eig      <- dudi_obj$eig
  pct      <- eig / sum(eig)
  cum_pct  <- cumsum(pct)
  k        <- max(1L, which(cum_pct >= threshold)[1])  # 至少保留 1 轴
  k        <- min(k, length(eig))                       # 不超过总轴数
  
  cat(sprintf(
    "  [%s] 总轴 %d，保留前 %d 轴，累计方差 %.1f%%\n",
    label, length(eig), k, cum_pct[k] * 100
  ))
  
  # 截断坐标（样本在前 k 个主轴上的得分）
  scores_trunc <- as.data.frame(dudi_obj$li)[, seq_len(k), drop = FALSE]
  colnames(scores_trunc) <- paste0("PC", seq_len(k))
  
  # 重建一个 nf=k 的 dudi 对象供 coinertia 使用
  dudi_trunc        <- dudi_obj
  dudi_trunc$nf     <- k
  dudi_trunc$li     <- dudi_obj$li[,  seq_len(k), drop = FALSE]
  dudi_trunc$l1     <- dudi_obj$l1[,  seq_len(k), drop = FALSE]
  dudi_trunc$co     <- dudi_obj$co[,  seq_len(k), drop = FALSE]
  dudi_trunc$c1     <- dudi_obj$c1[,  seq_len(k), drop = FALSE]
  dudi_trunc$eig    <- dudi_obj$eig[  seq_len(k)]
  
  list(
    k            = k,
    cum_pct      = cum_pct[k],
    scores_trunc = scores_trunc,
    dudi_trunc   = dudi_trunc,
    eig_full     = eig,
    pct_full     = pct,
    cum_pct_full = cum_pct
  )
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
  obs_u2  <- as.numeric(watson.two.test(a1, a2)$statistic)
  x_all   <- c(x1, x2); n1 <- length(x1); n_all <- length(x_all)
  perm_u2 <- replicate(B, {
    idx <- sample.int(n_all)
    as.numeric(watson.two.test(
      circular(x_all[idx[1:n1]],             type = "angles", units = "radians", modulo = "2pi"),
      circular(x_all[idx[(n1 + 1):n_all]],   type = "angles", units = "radians", modulo = "2pi")
    )$statistic)
  })
  list(
    statistic = obs_u2,
    p.value   = (sum(perm_u2 >= obs_u2) + 1) / (B + 1)
  )
}


# ==============================================================================
# ---- 数据读取与预处理 ----
# ==============================================================================

cat("==== 读取数据 ====\n")
SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"),
                              show_col_types = FALSE)
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"),
                              show_col_types = FALSE)
metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

filter_spharm <- function(df, power_cols, meta = NULL) {
  result <- df %>%
    select(ID, Typology, SHE, spectral_entropy, all_of(power_cols))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}

SPHARM_direction_filter  <- filter_spharm(SPHARM_direction,  POWER_COLS_DIR,   metric_data)
SPHARM_morphology_filter <- filter_spharm(SPHARM_morphology, POWER_COLS_MORPH, metric_data)

split_by_group <- function(df) {
  list(
    exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
    sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_"))
  )
}
dir_splits   <- split_by_group(SPHARM_direction_filter)
morph_splits <- split_by_group(SPHARM_morphology_filter)

df_scar_all  <- dir_splits$exp_im
df_morph_all <- morph_splits$exp_im

common_ids   <- intersect(df_morph_all$ID, df_scar_all$ID)
df_morph_all <- df_morph_all %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar_all  <- df_scar_all  %>% filter(ID %in% common_ids) %>% arrange(ID)
stopifnot(all(df_morph_all$ID == df_scar_all$ID))

cat("共有标本：", length(common_ids), "\n")

# 原始功率谱矩阵
morph_power_all <- df_morph_all %>%
  select(all_of(POWER_COLS_MORPH)) %>%
  rename_with(~ paste0("M", seq_along(.))) %>%
  as.data.frame()
scar_power_all  <- df_scar_all %>%
  select(all_of(POWER_COLS_DIR)) %>%
  rename_with(~ paste0("S", seq_along(.))) %>%
  as.data.frame()
rownames(morph_power_all) <- df_morph_all$ID
rownames(scar_power_all)  <- df_scar_all$ID

morph_power_clean <- morph_power_all[, sapply(morph_power_all, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power_all[,  sapply(scar_power_all,  sd, na.rm = TRUE) > 0]

# ILR 变换
morph_ilr_all <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power_clean))))
scar_ilr_all  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power_clean))))
rownames(morph_ilr_all) <- rownames(morph_power_clean)
rownames(scar_ilr_all)  <- rownames(scar_power_clean)

# 筛选纯实验标本
exp_ids <- rownames(morph_power_clean)[
  !str_starts(rownames(morph_power_clean), "IM_") &
    rownames(morph_power_clean) != "EXP43_Biface"]
cat("纯实验标本数量（不含 IM_）：", length(exp_ids), "\n")

morph_ilr_exp <- morph_ilr_all[exp_ids, ]
scar_ilr_exp  <- scar_ilr_all[exp_ids, ]

# 元数据
meta_exp <- df_morph_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, Typology) %>%
  left_join(metric_data, by = "ID") %>%
  mutate(Typology = if_else(
    str_detect(Typology, regex("levallois", ignore_case = TRUE)),
    "Levallois", Typology
  ))

df_morph_all <- df_morph_all %>%
  mutate(Typology = if_else(
    str_detect(Typology, regex("levallois", ignore_case = TRUE)),
    "Levallois", Typology
  ))

cat("\nTypology 分布：\n"); print(table(meta_exp$Typology, useNA = "ifany"))

typology_levels <- TYPOLOGY_ORDER[TYPOLOGY_ORDER %in%
                                    unique(meta_exp$Typology[
                                      meta_exp$Typology != "Biface" &
                                        !is.na(meta_exp$Typology)])]


# ==============================================================================
# ---- PCA + 95% 截断 ----
# ==============================================================================

cat("\n\n## PCA + 95% 方差截断 ##\n")

colnames(morph_ilr_exp) <- paste0("M_ilr", seq_len(ncol(morph_ilr_exp)))
colnames(scar_ilr_exp)  <- paste0("S_ilr", seq_len(ncol(scar_ilr_exp)))

# 全量 PCA（nf = 所有轴）
dudi_morph_full <- dudi.pca(morph_ilr_exp, center = TRUE, scale = TRUE,
                            scannf = FALSE, nf = ncol(morph_ilr_exp))
dudi_scar_full  <- dudi.pca(scar_ilr_exp,  center = TRUE, scale = TRUE,
                            scannf = FALSE, nf = ncol(scar_ilr_exp))

# 截断
trunc_morph <- truncate_pca(dudi_morph_full, VAR_THRESHOLD, "Morphology")
trunc_scar  <- truncate_pca(dudi_scar_full,  VAR_THRESHOLD, "Scar Direction")

# 截断后的 PCA 坐标（用于距离矩阵 & PERMANOVA）
morph_pca_trunc <- trunc_morph$scores_trunc
scar_pca_trunc  <- trunc_scar$scores_trunc
rownames(morph_pca_trunc) <- exp_ids
rownames(scar_pca_trunc)  <- exp_ids

# 截断后欧氏距离矩阵
D_morph_trunc <- dist(morph_pca_trunc)
D_scar_trunc  <- dist(scar_pca_trunc)

# 截断摘要
trunc_summary <- bind_rows(
  tibble(
    endpoint        = "Morphology",
    n_ilr_vars      = ncol(morph_ilr_exp),
    n_pca_full      = length(trunc_morph$eig_full),
    k_retained      = trunc_morph$k,
    cum_var_pct     = round(trunc_morph$cum_pct * 100, 2),
    threshold_pct   = VAR_THRESHOLD * 100,
    eigenvalues_all = paste(round(trunc_morph$eig_full, 4), collapse = "; "),
    pct_var_all     = paste(round(trunc_morph$pct_full * 100, 2), collapse = "; ")
  ),
  tibble(
    endpoint        = "Scar Direction",
    n_ilr_vars      = ncol(scar_ilr_exp),
    n_pca_full      = length(trunc_scar$eig_full),
    k_retained      = trunc_scar$k,
    cum_var_pct     = round(trunc_scar$cum_pct * 100, 2),
    threshold_pct   = VAR_THRESHOLD * 100,
    eigenvalues_all = paste(round(trunc_scar$eig_full, 4), collapse = "; "),
    pct_var_all     = paste(round(trunc_scar$pct_full * 100, 2), collapse = "; ")
  )
)
cat("\n==== PCA 截断摘要 ====\n")
print(trunc_summary %>% select(-eigenvalues_all, -pct_var_all) %>% as.data.frame())
write_csv(trunc_summary, here("analysis/data/derived_data/EXP_PCA_truncation_summary.csv"))
cat("已保存：EXP_PCA_truncation_summary.csv\n")

# 逐轴方差细节
cat("\n形态端各轴方差：\n")
tibble(
  PC      = paste0("PC", seq_along(trunc_morph$pct_full)),
  pct     = round(trunc_morph$pct_full * 100, 2),
  cum_pct = round(trunc_morph$cum_pct_full * 100, 2),
  retained = seq_along(trunc_morph$pct_full) <= trunc_morph$k
) %>% print(n = Inf)

cat("\n方向端各轴方差：\n")
tibble(
  PC       = paste0("PC", seq_along(trunc_scar$pct_full)),
  pct      = round(trunc_scar$pct_full * 100, 2),
  cum_pct  = round(trunc_scar$cum_pct_full * 100, 2),
  retained = seq_along(trunc_scar$pct_full) <= trunc_scar$k
) %>% print(n = Inf)

# 保存截断 PCA 坐标
morph_pca_trunc %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/EXP_morph_PCA_scores_trunc.csv"))
scar_pca_trunc %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/EXP_scar_PCA_scores_trunc.csv"))
cat("已保存：截断 PCA 坐标（morph / scar）\n")


# ==============================================================================
# ---- L1：全局 PERMANOVA + Mantel + CoIA ----
# ==============================================================================

cat("\n\n## 第一层：全局检验（截断 PCA 距离）##\n")

# ---------- L1-1：全局 Mantel ----------
cat("\n---- 全局 Mantel（截断 PCA 欧氏距离，Spearman）----\n")
set.seed(42)
mantel_global <- mantel(D_morph_trunc, D_scar_trunc,
                        method = "spearman", permutations = 9999)
print(mantel_global)

# ---------- L1-2：全局 PERMANOVA ----------
cat("\n---- 全局 PERMANOVA（截断 PCA 坐标，adonis2）----\n")

# 将形态 PCA 坐标与方向 PCA 坐标拼合，检验两端的联合结构
# 策略：把方向 PCA 坐标作为"解释变量矩阵"，以形态距离为响应
# 同时做反向（以方向距离为响应，形态 PCA 为解释）

set.seed(42)
perm_morph_resp <- adonis2(
  D_morph_trunc ~ .,
  data        = scar_pca_trunc,
  method      = "euclidean",
  permutations = 9999,
  by          = "margin"
)
cat("\n[PERMANOVA] 响应 = 形态距离，解释 = 方向 PCA 坐标：\n")
print(perm_morph_resp)

set.seed(42)
perm_scar_resp <- adonis2(
  D_scar_trunc ~ .,
  data        = morph_pca_trunc,
  method      = "euclidean",
  permutations = 9999,
  by          = "margin"
)
cat("\n[PERMANOVA] 响应 = 方向距离，解释 = 形态 PCA 坐标：\n")
print(perm_scar_resp)

# 提取 R² 与 p（总模型行）
extract_perm_global <- function(adonis_obj, response_label, predictor_label) {
  # adonis2 的最后一行是 "Total"，倒数第二行是 "Residual"，之前是各预测变量
  df_res  <- as.data.frame(adonis_obj)
  total_r2 <- sum(df_res$R2[!rownames(df_res) %in% c("Residual", "Total")],
                  na.rm = TRUE)
  min_p    <- min(df_res$`Pr(>F)`, na.rm = TRUE)
  tibble(
    response    = response_label,
    predictor   = predictor_label,
    R2_total    = round(total_r2, 4),
    min_p_value = round(min_p,    4)
  )
}

perm_global_summary <- bind_rows(
  extract_perm_global(perm_morph_resp, "Morphology distance",      "Scar PCA coords"),
  extract_perm_global(perm_scar_resp,  "Scar direction distance",  "Morphology PCA coords")
)
cat("\n==== 全局 PERMANOVA 汇总 ====\n")
print(perm_global_summary)
write_csv(perm_global_summary,
          here("analysis/data/derived_data/EXP_L1_PERMANOVA_results.csv"))
cat("已保存：EXP_L1_PERMANOVA_results.csv\n")

# ---------- L1-3：CoIA（截断 dudi 对象）----------
cat("\n---- CoIA（截断 PCA 轴）----\n")
coin_exp    <- coinertia(trunc_morph$dudi_trunc, trunc_scar$dudi_trunc,
                         scannf = FALSE, nf = 2)
cia_inertia <- coin_exp$eig / sum(coin_exp$eig) * 100

cat("RV 系数：", round(coin_exp$RV, 4), "\n")
set.seed(42)
rv_test <- randtest(coin_exp, nrepet = 9999)
cat("\nRV 置换检验：\n"); print(rv_test)

# 样本坐标
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

cat("\n==== CoIA 样本坐标（前6行）====\n")
print(head(scores_combined %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))

cia_coords <- scores_combined %>%
  select(ID, Typology,
         Morph_Axis1 = Axis1_M, Morph_Axis2 = Axis2_M,
         Scar_Axis1  = Axis1_S, Scar_Axis2  = Axis2_S,
         arrow_length, arrow_angle)
write_csv(cia_coords, here("analysis/data/derived_data/EXP_CIA_coords_full.csv"))
write_csv(scores_combined, here("analysis/data/derived_data/EXP_CIA_scores_full.csv"))
cat("已保存：EXP_CIA_coords_full.csv / EXP_CIA_scores_full.csv\n")

# PCA 轴对 CoIA 轴贡献
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

morph_pct_trunc <- round(trunc_morph$pct_full[seq_len(trunc_morph$k)] * 100, 1)
scar_pct_trunc  <- round(trunc_scar$pct_full[ seq_len(trunc_scar$k)]  * 100, 1)

morph_contrib <- compute_pca_cia_contribution(coin_exp$aX, morph_pct_trunc, "Morphology")
scar_contrib  <- compute_pca_cia_contribution(coin_exp$aY, scar_pct_trunc,  "Scar direction")

pca_cia_contrib <- bind_rows(morph_contrib, scar_contrib)
cat("\n==== PCA 轴对 CoIA 轴贡献 ====\n")
cat("形态端：\n"); print(morph_contrib %>% select(-endpoint) %>% as.data.frame())
cat("方向端：\n"); print(scar_contrib  %>% select(-endpoint) %>% as.data.frame())
write_csv(pca_cia_contrib,
          here("analysis/data/derived_data/EXP_PCA_CoIA_contribution.csv"))
cat("已保存：EXP_PCA_CoIA_contribution.csv\n")

l1_results <- tibble(
  method  = c("Mantel (trunc-PCA, Euclidean, Spearman)", "RV (trunc-PCA, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_exp$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(exp_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/EXP_L1_results.csv"))
cat("\n第一层结论：\n"); print(l1_results)


# ==============================================================================
# ---- L2-A：分组 PERMANOVA（Typology）----
# ==============================================================================

cat("\n\n## L2-A：分组 PERMANOVA（Typology，截断 PCA 距离）##\n")

# 合并截断坐标供单表分析
morph_scar_combined <- cbind(morph_pca_trunc, scar_pca_trunc)

permanova_one_type <- function(type_val, meta_df,
                               D_morph_full, D_scar_full,
                               morph_df, scar_df,
                               n_perm = 9999) {
  ids <- meta_df %>% filter(Typology == type_val) %>% pull(ID)
  cat(sprintf("  -> %s (n = %d) ...", type_val, length(ids)))
  if (length(ids) < 5) {
    cat(" 跳过（n < 5）\n")
    return(NULL)
  }
  # 形态距离 ~ 方向 PCA 坐标
  D_m  <- extract_subdist(D_morph_full, ids)
  sub_scar <- scar_df[ids, , drop = FALSE]
  set.seed(42)
  pm <- adonis2(D_m ~ ., data = as.data.frame(sub_scar),
                method = "euclidean", permutations = n_perm, by = "margin")
  r2_m <- sum(as.data.frame(pm)$R2[
    !rownames(as.data.frame(pm)) %in% c("Residual", "Total")], na.rm = TRUE)
  p_m  <- min(as.data.frame(pm)$`Pr(>F)`, na.rm = TRUE)
  
  # 方向距离 ~ 形态 PCA 坐标
  D_s  <- extract_subdist(D_scar_full, ids)
  sub_morph <- morph_df[ids, , drop = FALSE]
  set.seed(42)
  ps <- adonis2(D_s ~ ., data = as.data.frame(sub_morph),
                method = "euclidean", permutations = n_perm, by = "margin")
  r2_s <- sum(as.data.frame(ps)$R2[
    !rownames(as.data.frame(ps)) %in% c("Residual", "Total")], na.rm = TRUE)
  p_s  <- min(as.data.frame(ps)$`Pr(>F)`, na.rm = TRUE)
  
  cat(sprintf(" morph~scar R²=%.3f p=%.4f | scar~morph R²=%.3f p=%.4f\n",
              r2_m, p_m, r2_s, p_s))
  
  tibble(
    Typology        = type_val,
    n               = length(ids),
    R2_morph_resp   = round(r2_m, 4),
    p_morph_resp    = round(p_m,  4),
    R2_scar_resp    = round(r2_s, 4),
    p_scar_resp     = round(p_s,  4),
    sig_morph       = case_when(
      p_m < 0.001 ~ "***", p_m < 0.01 ~ "**",
      p_m < 0.05  ~ "*",   p_m < 0.10 ~ ".", TRUE ~ "ns"),
    sig_scar        = case_when(
      p_s < 0.001 ~ "***", p_s < 0.01 ~ "**",
      p_s < 0.05  ~ "*",   p_s < 0.10 ~ ".", TRUE ~ "ns")
  )
}

all_types <- meta_exp %>%
  filter(!is.na(Typology)) %>%
  count(Typology) %>%
  arrange(desc(n)) %>%
  pull(Typology)

permanova_by_typology <- map_dfr(
  all_types,
  ~ permanova_one_type(.x, meta_exp,
                       D_morph_trunc, D_scar_trunc,
                       morph_pca_trunc, scar_pca_trunc)
)

cat("\n==== 各类型 PERMANOVA 结果汇总 ====\n")
print(permanova_by_typology %>% as.data.frame())
write_csv(permanova_by_typology,
          here("analysis/data/derived_data/EXP_L2_grouped_permanova.csv"))
cat("已保存：EXP_L2_grouped_permanova.csv\n")

# 保留分组 Mantel 供比较
cat("\n---- 分组 Mantel（对比用）----\n")
mantel_one_type <- function(type_val, meta_df, D_morph_full, D_scar_full,
                            n_perm = 9999) {
  ids <- meta_df %>% filter(Typology == type_val) %>% pull(ID)
  cat(sprintf("  -> %s (n = %d) ...", type_val, length(ids)))
  if (length(ids) < 5) { cat(" 跳过\n"); return(NULL) }
  set.seed(42)
  res <- mantel(extract_subdist(D_morph_full, ids),
                extract_subdist(D_scar_full,  ids),
                method = "spearman", permutations = n_perm)
  cat(sprintf(" r = %.4f, p = %.4f\n", res$statistic, res$signif))
  tibble(Typology = type_val, n = length(ids),
         mantel_r = res$statistic, p_value = res$signif,
         significance = case_when(
           res$signif < 0.001 ~ "***", res$signif < 0.01 ~ "**",
           res$signif < 0.05  ~ "*",   res$signif < 0.10 ~ ".", TRUE ~ "ns"))
}

mantel_by_typology <- map_dfr(all_types,
                              ~ mantel_one_type(.x, meta_exp,
                                                D_morph_trunc, D_scar_trunc))
write_csv(mantel_by_typology,
          here("analysis/data/derived_data/EXP_L2_grouped_mantel.csv"))
cat("已保存：EXP_L2_grouped_mantel.csv\n")


# ==============================================================================
# ---- L2-B：CoIA 箭头长度（Typology）----
# ==============================================================================

cat("\n\n## L2-B：箭头长度 × Typology ##\n")

valid_groups_len <- scores_combined %>%
  filter(!is.na(Typology), Typology != "Biface") %>%
  group_by(Typology) %>% filter(n() >= 3) %>%
  pull(Typology) %>% unique()

sub_len <- scores_combined %>%
  filter(Typology %in% valid_groups_len) %>%
  mutate(Typology = factor(Typology,
                           levels = intersect(TYPOLOGY_ORDER, valid_groups_len)))

kw_len <- kruskal.test(arrow_length ~ Typology, data = sub_len)
cat("Kruskal-Wallis（箭头长度 ~ Typology）：\n"); print(kw_len)

pw_len <- pairwise.wilcox.test(sub_len$arrow_length, sub_len$Typology,
                               p.adjust.method = "holm", exact = FALSE)
cat("Pairwise Wilcoxon（Holm 校正）：\n"); print(pw_len)

desc_len <- sub_len %>%
  group_by(Typology) %>%
  summarise(
    n      = n(),
    mean   = round(mean(arrow_length),   4),
    median = round(median(arrow_length), 4),
    sd     = round(sd(arrow_length),     4),
    min    = round(min(arrow_length),    4),
    max    = round(max(arrow_length),    4),
    .groups = "drop"
  )
cat("\n描述统计：\n"); print(desc_len)

scores_combined %>%
  select(ID, Typology, arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/EXP_L2_arrow_stats.csv"))
cat("已保存：EXP_L2_arrow_stats.csv\n")


# ==============================================================================
# ---- L2-C：箭头方位圆形统计（Typology）----
# ==============================================================================

cat("\n\n## L2-C：箭头方位圆形统计 ##\n")

valid_groups_circ <- scores_combined %>%
  filter(!is.na(Typology), Typology != "Biface") %>%
  group_by(Typology) %>% filter(n() >= 5) %>%
  pull(Typology) %>% unique()
valid_groups_circ <- intersect(TYPOLOGY_ORDER, valid_groups_circ)

circ_desc <- map_dfr(valid_groups_circ, function(g) {
  angles <- scores_combined %>%
    filter(Typology == g) %>% pull(arrow_angle)
  cs <- circ_stats_one(angles)
  tibble(group = g, n = length(angles),
         mean_dir_deg    = round(cs$mean_deg, 2),
         concentration_r = round(cs$rho, 4))
})
cat("圆形描述统计：\n"); print(circ_desc)

rayleigh_res <- map_dfr(valid_groups_circ, function(g) {
  angles   <- scores_combined %>% filter(Typology == g) %>% pull(arrow_angle)
  circ_obj <- circular(angles, type = "angles", units = "radians", modulo = "2pi")
  rt <- rayleigh.test(circ_obj)
  cat(sprintf("  %s: U = %.4f, p = %.4f\n", g, rt$statistic, rt$p.value))
  tibble(group = g,
         rayleigh_U  = round(rt$statistic, 4),
         rayleigh_p  = round(rt$p.value,   4),
         conclusion  = ifelse(rt$p.value < 0.05, "concentrated", "uniform"))
})

cat("\nWatson 两样本检验（各对）：\n")
pairs_circ <- combn(valid_groups_circ, 2, simplify = FALSE)
watson_res <- map_dfr(pairs_circ, function(pair) {
  x1 <- scores_combined %>% filter(Typology == pair[1]) %>% pull(arrow_angle)
  x2 <- scores_combined %>% filter(Typology == pair[2]) %>% pull(arrow_angle)
  wt <- watson_perm_test(x1, x2, B = 9999)
  cat(sprintf("  %s vs %s: U2 = %.4f, p = %.4f\n",
              pair[1], pair[2], wt$statistic, wt$p.value))
  tibble(group1 = pair[1], group2 = pair[2],
         U2_statistic = round(wt$statistic, 4),
         p_value      = round(wt$p.value,   4),
         conclusion   = ifelse(wt$p.value < 0.05, "different", "ns"))
})

circ_out <- left_join(circ_desc, rayleigh_res, by = "group")
write_csv(circ_out, here("analysis/data/derived_data/EXP_L2_circular_stats.csv"))
cat("\n已保存：EXP_L2_circular_stats.csv\n")
cat("\nWatson 结果：\n"); print(watson_res)


# ==============================================================================
# ---- L2-D：spectral_entropy × Typology ----
# ==============================================================================

cat("\n\n## L2-D：spectral_entropy × Typology ##\n")

se_df <- df_scar_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, SE_direction = spectral_entropy) %>%
  left_join(
    df_morph_all %>% filter(ID %in% exp_ids) %>%
      select(ID, SE_morphology = spectral_entropy),
    by = "ID"
  ) %>%
  left_join(meta_exp %>% select(ID, Typology), by = "ID") %>%
  filter(!is.na(Typology), Typology != "Biface") %>%
  group_by(Typology) %>% filter(n() >= 3) %>% ungroup() %>%
  mutate(Typology = factor(Typology,
                           levels = intersect(TYPOLOGY_ORDER, unique(Typology))))

se_desc <- se_df %>%
  group_by(Typology) %>%
  summarise(
    n            = n(),
    mean_dir     = round(mean(SE_direction,    na.rm = TRUE), 4),
    median_dir   = round(median(SE_direction,  na.rm = TRUE), 4),
    sd_dir       = round(sd(SE_direction,      na.rm = TRUE), 4),
    mean_morph   = round(mean(SE_morphology,   na.rm = TRUE), 4),
    median_morph = round(median(SE_morphology, na.rm = TRUE), 4),
    sd_morph     = round(sd(SE_morphology,     na.rm = TRUE), 4),
    .groups = "drop"
  )
cat("SE 描述统计：\n"); print(se_desc)
write_csv(se_desc, here("analysis/data/derived_data/EXP_L2D_SE_desc_stats.csv"))

run_kw_dunn_se <- function(df, y_col, label) {
  cat(sprintf("\n----- %s -----\n", label))
  kw <- kruskal.test(reformulate("Typology", y_col), data = df)
  cat(sprintf("  KW chi²=%.4f df=%d p=%.4f\n",
              kw$statistic, kw$parameter, kw$p.value))
  dunn_raw <- dunnTest(x = df[[y_col]], g = df[["Typology"]],
                       method = "bh")$res
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
  print(dunn)
  list(kw = kw, dunn = dunn)
}

res_se_dir   <- run_kw_dunn_se(se_df %>% filter(!is.na(SE_direction)),
                               "SE_direction",  "SE（方向）× Typology")
res_se_morph <- run_kw_dunn_se(se_df %>% filter(!is.na(SE_morphology)),
                               "SE_morphology", "SE（形态）× Typology")

bind_rows(
  res_se_dir$dunn   %>% mutate(source = "direction"),
  res_se_morph$dunn %>% mutate(source = "morphology")
) %>%
  write_csv(here("analysis/data/derived_data/EXP_L2D_SE_dunn_results.csv"))
cat("已保存：EXP_L2D_SE_dunn_results.csv\n")


# ==============================================================================
# ---- 保存原始 ILR 得分 ----
# ==============================================================================

morph_ilr_exp %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/EXP_morph_ILR_scores.csv"))
scar_ilr_exp %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/EXP_scar_ILR_scores.csv"))
cat("已保存：ILR 得分（morph / scar）\n")


# ==============================================================================
# ---- 最终汇总打印 ----
# ==============================================================================

cat("\n\n========== 分析结果汇总 ==========\n")

cat("\n[PCA 截断]\n")
cat(sprintf("  形态端：%d ILR -> 保留 %d PC（累计方差 %.1f%%）\n",
            ncol(morph_ilr_exp), trunc_morph$k, trunc_morph$cum_pct * 100))
cat(sprintf("  方向端：%d ILR -> 保留 %d PC（累计方差 %.1f%%）\n",
            ncol(scar_ilr_exp),  trunc_scar$k,  trunc_scar$cum_pct  * 100))

cat("\n[L1 全局检验]\n")
cat(sprintf("  Mantel r = %.4f, p = %.4f  -> %s\n",
            mantel_global$statistic, mantel_global$signif,
            ifelse(mantel_global$signif < 0.05, "显著", "ns")))
cat(sprintf("  RV       = %.4f, p = %.4f  -> %s\n",
            coin_exp$RV, rv_test$pvalue,
            ifelse(rv_test$pvalue < 0.05, "显著协变", "独立")))
cat(sprintf("  PERMANOVA (morph~scar): R²=%.4f, p=%.4f\n",
            perm_global_summary$R2_total[1], perm_global_summary$min_p_value[1]))
cat(sprintf("  PERMANOVA (scar~morph): R²=%.4f, p=%.4f\n",
            perm_global_summary$R2_total[2], perm_global_summary$min_p_value[2]))

cat("\n[L2-A 分组 PERMANOVA]\n")
print(permanova_by_typology %>%
        select(Typology, n, R2_morph_resp, sig_morph, R2_scar_resp, sig_scar) %>%
        as.data.frame())

cat("\n[L2-D SE × Typology]\n")
cat(sprintf("  方向 SE KW: p = %.4f  形态 SE KW: p = %.4f\n",
            res_se_dir$kw$p.value, res_se_morph$kw$p.value))

cat("\n所有输出已保存至 analysis/data/derived_data/\n")