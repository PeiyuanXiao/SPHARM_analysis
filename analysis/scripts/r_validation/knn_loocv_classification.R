# ==============================================================================
# mrm_variance_decomposition.R
# 基于 MRM（距离矩阵多元回归）的三方法方差分解
#
# 框架：Lichstein (2007) Oecologia — Multiple Regression on distance Matrices
# 实现：ecodist::MRM()（置换检验替代参数检验处理非独立性）
#
# 分析流程：
#   第一步：三组 MRM，轮换因变量（各方法 → 另两个方法为预测变量）
#           → 每组得到 R²（他人可解释比例）和 1-R²（独有信息比例）
#   第二步：以 SPHARM 为因变量，额外运行两个简单 MRM（单预测变量）
#           → 完整方差分解（SPI 独有 / Fabric 独有 / 共享 / SPHARM 独有）
#   第三步：提取 SPHARM ~ SPI + Fabric 的残差矩阵
#           → 按残差大小排序标本对，识别"SPHARM 独有差异"案例
#
# 注意：MRM 默认使用 Pearson（mrank=FALSE），与 partial Mantel 的 Spearman 不同
#       R² 解释方式与 OLS 一致，但 p 值来自置换检验
#
# 输入：（与 mantel_method_redundancy.R 相同数据源）
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 3）
#   - analysis/data/derived_data/directions_aligned_svd.csv
#   - analysis/data/derived_data/SPHARM_direction.csv
#
# 输出：
#   - analysis/output/figures/mrm_variance_decomp.png
#   - analysis/data/derived_data/mrm_results.csv
#   - analysis/data/derived_data/mrm_spharm_residuals.csv
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(ecodist)
library(patchwork)
library(glue)
library(gridExtra) 
library(grid)

# ==============================================================================
# 参数
# ==============================================================================

NPERM <- 999
set.seed(42)

# ==============================================================================
# 工具函数（与既有脚本保持一致）
# ==============================================================================

compute_R <- function(ux, uy, uz) {
  sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2) /
    sum(sqrt(ux^2 + uy^2 + uz^2))
}

compute_EI <- function(ux, uy, uz) {
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / nrow(U)
  lambda <- sort(eigen(T_mat, symmetric = TRUE)$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_)
  )
}

prepare_scaled_matrix <- function(X) {
  X      <- as.matrix(X)
  col_sd <- apply(X, 2, sd, na.rm = TRUE)
  keep   <- is.finite(col_sd) & col_sd > 1e-10
  X      <- X[, keep, drop = FALSE]
  if (ncol(X) == 0) stop("所有列均为常数列，无法标准化。")
  scale(X)
}

dist_euclidean_scaled <- function(X) {
  as.dist(stats::dist(prepare_scaled_matrix(X), method = "euclidean"))
}

make_sig <- function(p) case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**",
  p < 0.05  ~ "*",   TRUE     ~ "ns"
)

# ==============================================================================
# 读取数据
# ==============================================================================

exp_data <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 3)

directions_aligned <- read_csv(
  here("analysis/data/derived_data/directions_aligned_svd.csv"),
  show_col_types = FALSE)

spharm_all <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE)

# ==============================================================================
# 过滤实验石核，计算 SPI / Fabric
# ==============================================================================

exp_ids        <- unique(exp_data$ID)
directions_exp <- directions_aligned %>% filter(ID %in% exp_ids)

rei_exp <- directions_exp %>%
  group_by(ID) %>%
  summarise(
    SPI = compute_R(ux, uy, uz),
    E   = compute_EI(ux, uy, uz)$E,
    I   = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

# ==============================================================================
# 合并 SPHARM 数据
# ==============================================================================

power_cols <- paste0("power_l", 1:5)

if (!all(power_cols %in% colnames(spharm_all))) {
  stop("SPHARM 文件中缺少 power_l1~power_l5 列。")
}

spharm_exp <- spharm_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, all_of(power_cols))

df <- rei_exp %>%
  left_join(spharm_exp, by = "ID") %>%
  filter(!is.na(SPI), !is.na(E), !is.na(I), !is.na(power_l1)) %>%
  arrange(ID)

n   <- nrow(df)
ids <- df$ID

cat(sprintf("进入分析的标本数：%d\n\n", n))

# ==============================================================================
# 构造距离矩阵
# ==============================================================================

D_spi <- dist_euclidean_scaled(df %>% select(SPI))
D_fab <- dist_euclidean_scaled(df %>% select(E, I))
D_sph <- dist_euclidean_scaled(df %>% select(all_of(power_cols)))

# ==============================================================================
# 第一步：三组 MRM，轮换因变量
# ==============================================================================

cat("======================================================\n")
cat(" 第一步：三组 MRM（轮换因变量）\n")
cat(sprintf(" nperm = %d，Pearson（mrank = FALSE）\n", NPERM))
cat("======================================================\n\n")

set.seed(42); mrm_sph <- MRM(D_sph ~ D_spi + D_fab, nperm = NPERM, mrank = FALSE)
cat("  SPHARM ~ SPI + Fabric  完成\n")
set.seed(42); mrm_spi <- MRM(D_spi ~ D_sph + D_fab, nperm = NPERM, mrank = FALSE)
cat("  SPI ~ SPHARM + Fabric  完成\n")
set.seed(42); mrm_fab <- MRM(D_fab ~ D_spi + D_sph, nperm = NPERM, mrank = FALSE)
cat("  Fabric ~ SPI + SPHARM  完成\n\n")

# 结果提取辅助函数
# ecodist::MRM 输出：
#   $coef     : matrix，行 = 截距 + 各预测变量，列 = coef / pval
#   $r.squared: named vector，R2 + pval
#   $F.test   : named vector，F + F.pval
extract_mrm_summary <- function(mrm_obj, dv_label) {
  r2   <- as.numeric(mrm_obj$r.squared["R2"])
  r2_p <- as.numeric(mrm_obj$r.squared["pval"])
  coef_mat <- mrm_obj$coef
  
  # 用位置索引：第 1 列是系数，第 2 列是 p 值
  # 行顺序：截距(1), 第一预测变量(2), 第二预测变量(3)
  pred_rows <- coef_mat[-1, , drop = FALSE]
  coef_tbl  <- tibble(
    dv        = dv_label,
    predictor = rownames(pred_rows),
    beta      = as.numeric(pred_rows[, 1]),   # 系数
    p_beta    = as.numeric(pred_rows[, 2]),   # p 值
    sig_beta  = make_sig(as.numeric(pred_rows[, 2]))
  )
  
  list(
    r2       = r2,
    unique   = 1 - r2,
    r2_pval  = r2_p,
    r2_sig   = make_sig(r2_p),
    dv       = dv_label,
    coef_tbl = coef_tbl
  )
}

res_sph <- extract_mrm_summary(mrm_sph, "SPHARM power")
res_spi <- extract_mrm_summary(mrm_spi, "SPI")
res_fab <- extract_mrm_summary(mrm_fab, "Fabric")

summary_main <- tibble(
  dv      = c("SPHARM power", "SPI", "Fabric"),
  r2      = c(res_sph$r2,     res_spi$r2,     res_fab$r2),
  unique  = c(res_sph$unique, res_spi$unique, res_fab$unique),
  r2_pval = c(res_sph$r2_pval,res_spi$r2_pval,res_fab$r2_pval),
  r2_sig  = c(res_sph$r2_sig, res_spi$r2_sig, res_fab$r2_sig)
)

coef_main <- bind_rows(
  res_sph$coef_tbl,
  res_spi$coef_tbl,
  res_fab$coef_tbl
)

# 打印主分析结果
cat("R²（另两个方法联合可解释的方差比例）：\n\n")
cat(sprintf("  %-22s  %7s  %7s  %8s  %s\n",
            "因变量", "R²", "1-R²", "p(R²)", "sig"))
cat("  ", paste(rep("-", 58), collapse = ""), "\n", sep = "")
for (i in seq_len(nrow(summary_main))) {
  r <- summary_main[i, ]
  cat(sprintf("  %-22s  %7.4f  %7.4f  %8.4f  %s\n",
              r$dv, r$r2, r$unique, r$r2_pval, r$r2_sig))
}

cat("\n回归系数：\n\n")
cat(sprintf("  %-22s  %-25s  %8s  %8s  %s\n",
            "因变量", "预测变量", "β", "p(β)", "sig"))
cat("  ", paste(rep("-", 76), collapse = ""), "\n", sep = "")
for (i in seq_len(nrow(coef_main))) {
  r <- coef_main[i, ]
  cat(sprintf("  %-22s  %-25s  %8.4f  %8.4f  %s\n",
              r$dv, r$predictor, r$beta, r$p_beta, r$sig_beta))
}

# ==============================================================================
# 第二步：SPHARM 为因变量的完整方差分解
# 额外运行两个简单 MRM（单个预测变量）
# 方差分解公式（2 个预测变量的标准分拆）：
#   [a] SPI 独有    = R²(sph ~ spi+fab) - R²(sph ~ fab only)
#   [b] Fabric 独有 = R²(sph ~ spi+fab) - R²(sph ~ spi only)
#   [c] 共享        = R²(sph~spi) + R²(sph~fab) - R²(sph~spi+fab)
#   [d] 不可解释    = 1 - R²(sph ~ spi+fab)   ← SPHARM 独有信息
# ==============================================================================

cat("\n======================================================\n")
cat(" 第二步：SPHARM 完整方差分解\n")
cat("======================================================\n\n")

set.seed(42); mrm_sph_spi_only <- MRM(D_sph ~ D_spi, nperm = NPERM, mrank = FALSE)
set.seed(42); mrm_sph_fab_only <- MRM(D_sph ~ D_fab, nperm = NPERM, mrank = FALSE)

r2_full     <- as.numeric(mrm_sph$r.squared["R2"])
r2_spi_only <- as.numeric(mrm_sph_spi_only$r.squared["R2"])
r2_fab_only <- as.numeric(mrm_sph_fab_only$r.squared["R2"])

vp_spi    <- r2_full - r2_fab_only                        # [a] SPI 独有
vp_fab    <- r2_full - r2_spi_only                        # [b] Fabric 独有
vp_shared <- r2_spi_only + r2_fab_only - r2_full          # [c] 共享（可能为负）
vp_resid  <- 1 - r2_full                                  # [d] SPHARM 独有

# 共线性估计误差可能导致 shared 或单独项为负，修正到 0（保持 [d] 不变）
vp_shared <- max(vp_shared, 0)
vp_spi    <- max(vp_spi,    0)
vp_fab    <- max(vp_fab,    0)

cat(sprintf("  R²(SPHARM ~ SPI only)   = %.4f\n", r2_spi_only))
cat(sprintf("  R²(SPHARM ~ Fabric only)= %.4f\n", r2_fab_only))
cat(sprintf("  R²(SPHARM ~ SPI+Fabric) = %.4f  [全模型]\n\n", r2_full))

vp_table <- tibble(
  component = c("[a] Unique to SPI",
                "[b] Unique to Fabric",
                "[c] Shared (SPI ∩ Fabric)",
                "[d] Unexplained — SPHARM unique"),
  value     = c(vp_spi, vp_fab, vp_shared, vp_resid),
  pct       = value * 100
)

cat("方差分解：\n\n")
for (i in seq_len(nrow(vp_table))) {
  cat(sprintf("  %-38s : %.4f  (%.1f%%)\n",
              vp_table$component[i],
              vp_table$value[i],
              vp_table$pct[i]))
}

cat(sprintf("\n  总计（应 = 1.0000）：%.4f\n", sum(vp_table$value))  )

# ==============================================================================
# 第三步：残差矩阵提取
# 来源：SPHARM ~ SPI + Fabric（全模型）
# 残差 = SPHARM 观测距离 - MRM 拟合距离
# 正残差：SPHARM 比 SPI+Fabric 预测的"更远"→ SPHARM 捕捉到了额外差异
# 负残差：SPHARM 比 SPI+Fabric 预测的"更近"→ SPHARM 低估了差异
# ==============================================================================

cat("\n======================================================\n")
cat(" 第三步：残差矩阵提取（SPHARM ~ SPI + Fabric）\n")
cat("======================================================\n\n")

# ecodist::MRM $coef 行顺序：截距(1), 第一预测变量(2), 第二预测变量(3)
coef_raw <- mrm_sph$coef[, 1]
b_int    <- coef_raw[1]
b_spi    <- coef_raw[2]
b_fab    <- coef_raw[3]

cat(sprintf("  截距 β₀ = %.4f\n", b_int))
cat(sprintf("  β(SPI)  = %.4f  (p = %.4f  %s)\n",
            b_spi,
            mrm_sph$coef["D_spi", 2],
            make_sig(mrm_sph$coef["D_spi", 2])))
cat(sprintf("  β(Fab)  = %.4f  (p = %.4f  %s)\n\n",
            b_fab,
            mrm_sph$coef["D_fab", 2],
            make_sig(mrm_sph$coef["D_fab", 2])))

# 下三角向量（与 MRM 内部用法一致）
vec_sph    <- as.vector(D_sph)
vec_spi    <- as.vector(D_spi)
vec_fab    <- as.vector(D_fab)

fitted_vec <- b_int + b_spi * vec_spi + b_fab * vec_fab
resid_vec  <- vec_sph - fitted_vec
resid_sd   <- sd(resid_vec)
resid_std_vec <- resid_vec / resid_sd

cat(sprintf("  残差均值：%.6f（理论接近 0）\n", mean(resid_vec)))
cat(sprintf("  残差标准差：%.4f\n", resid_sd))
cat(sprintf("  残差范围：[%.4f, %.4f]\n\n", min(resid_vec), max(resid_vec)))

# 重构为对称矩阵
resid_mat <- matrix(0, n, n)
resid_mat[lower.tri(resid_mat)] <- resid_vec
resid_mat <- resid_mat + t(resid_mat)
rownames(resid_mat) <- colnames(resid_mat) <- ids

# 构建标本对表（按残差绝对值排序）
pairs_idx <- which(upper.tri(resid_mat), arr.ind = TRUE)
D_spi_mat <- as.matrix(D_spi)
D_fab_mat <- as.matrix(D_fab)
D_sph_mat <- as.matrix(D_sph)

resid_pairs <- tibble(
  ID_i      = ids[pairs_idx[, 1]],
  ID_j      = ids[pairs_idx[, 2]],
  resid     = resid_mat[pairs_idx],
  resid_std = resid_std_vec,
  d_spharm  = D_sph_mat[pairs_idx],
  d_spi     = D_spi_mat[pairs_idx],
  d_fabric  = D_fab_mat[pairs_idx],
  fitted    = fitted_vec
) %>%
  arrange(desc(abs(resid_std)))

cat("  残差绝对值最大的标本对（前 10）：\n\n")
cat(sprintf("  %-18s  %-18s  %9s  %8s  %8s  %8s\n",
            "ID_i", "ID_j", "resid_std", "d_SPHARM", "d_SPI", "d_Fabric"))
cat("  ", paste(rep("-", 78), collapse = ""), "\n", sep = "")
resid_pairs %>%
  head(10) %>%
  pwalk(function(ID_i, ID_j, resid_std, d_spharm, d_spi, d_fabric, ...) {
    cat(sprintf("  %-18s  %-18s  %9.3f  %8.4f  %8.4f  %8.4f\n",
                ID_i, ID_j, resid_std, d_spharm, d_spi, d_fabric))
  })

n_pos  <- sum(resid_pairs$resid > 0)
n_neg  <- sum(resid_pairs$resid < 0)
n_out2 <- sum(abs(resid_pairs$resid_std) > 2)

cat(sprintf(
  "\n  正残差（SPHARM 距离 > 预测）：%d 对\n  负残差（SPHARM 距离 < 预测）：%d 对\n  |std.resid| > 2 的异常对：%d 对\n",
  n_pos, n_neg, n_out2
))

# ==============================================================================
# 可视化
# ==============================================================================

dv_order  <- c("SPI", "Fabric", "SPHARM power")
bar_pal   <- c(
  "Explained by other two (R²)" = "#7BAED4",
  "Unique / unexplained (1−R²)" = "#E8C4D8"
)

# ── 图 A：轮换因变量 R² 柱状图 ─────────────────────────────────────────────

label_df <- summary_main %>%
  mutate(
    dv    = factor(dv, levels = dv_order),
    label = glue("R²={round(r2,4)}\n({r2_sig})"),
    y_pos = r2 / 2
  )

p_rotating <- summary_main %>%
  mutate(dv = factor(dv, levels = dv_order)) %>%
  pivot_longer(c(r2, unique), names_to = "part", values_to = "val") %>%
  mutate(part = factor(part,
                       levels = c("r2", "unique"),
                       labels = names(bar_pal))) %>%
  ggplot(aes(x = dv, y = val, fill = part)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4, alpha = 0.92) +
  geom_text(
    data    = label_df,
    mapping = aes(x = dv, y = y_pos, label = label, fill = NULL),
    size = 3.6, fontface = "bold", color = "grey15"
  ) +
  scale_fill_manual(values = bar_pal, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    legend.position  = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title    = "MRM: Information Uniqueness — Rotating Dependent Variable",
    subtitle = glue(
      "R² = variance explained by the other two methods  |  ",
      "1−R² = unique / unexplained  |  nperm = {NPERM}"
    ),
    x = "Dependent Variable", y = "Proportion of Variance"
  )

# ── 图 B：SPHARM 方差分解 Venn（手绘双圆 + 标注）──────────────────────────

# 用 ggplot2 手绘两个部分重叠的圆，面积近似比例于 R²
# 圆 1（SPI）中心在 (-0.18, 0)，圆 2（Fabric）中心在 (0.18, 0)
# 重叠区域面积 ∝ vp_shared，总 R² = 左圆 + 右圆 - 重叠

theta   <- seq(0, 2 * pi, length.out = 300)
r_base  <- 0.32   # 基准半径（视觉用，不严格等比）
sep     <- 0.16   # 圆心间距控制重叠程度（共享越大，sep 越小）

# 根据 vp_shared 动态调整圆心间距（共享越多，圆越重叠）
sep_adj <- r_base * (1 - vp_shared / max(r2_full, 0.01))
sep_adj <- pmin(pmax(sep_adj, 0.05), 0.30)

cx1 <- -sep_adj; cy <- 0
cx2 <-  sep_adj

circ1 <- tibble(
  x   = r_base * cos(theta) + cx1,
  y   = r_base * sin(theta) + cy,
  grp = "SPI"
)
circ2 <- tibble(
  x   = r_base * cos(theta) + cx2,
  y   = r_base * sin(theta) + cy,
  grp = "Fabric"
)

p_venn <- ggplot() +
  geom_polygon(data = circ1,
               aes(x, y), fill = "#A1C2E6", alpha = 0.45, color = "#5B92C4", linewidth = 0.8) +
  geom_polygon(data = circ2,
               aes(x, y), fill = "#FFBAE0", alpha = 0.45, color = "#D472A8", linewidth = 0.8) +
  # 标注区域
  annotate("text", x = cx1 - r_base * 0.4, y = 0,
           label = glue("SPI unique\n[a] {round(vp_spi*100,1)}%"),
           size = 3.5, fontface = "bold", color = "#3A7BC8") +
  annotate("text", x = cx2 + r_base * 0.4, y = 0,
           label = glue("Fabric unique\n[b] {round(vp_fab*100,1)}%"),
           size = 3.5, fontface = "bold", color = "#C4528A") +
  annotate("text", x = 0, y = 0,
           label = glue("Shared\n[c] {round(vp_shared*100,1)}%"),
           size = 3.2, fontface = "bold", color = "grey25") +
  # 圆外标注不可解释部分
  annotate("text", x = 0, y = -(r_base + 0.15),
           label = glue("[d] Unexplained (SPHARM unique) = {round(vp_resid*100,1)}%"),
           size = 3.5, fontface = "italic", color = "grey35") +
  annotate("segment", x = 0, xend = 0,
           y = -(r_base + 0.08), yend = -(r_base - 0.02),
           color = "grey50", linewidth = 0.5,
           arrow = arrow(length = unit(0.2, "cm"))) +
  coord_equal(xlim = c(-0.72, 0.72), ylim = c(-0.68, 0.50)) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Variance Partitioning: SPHARM ~ SPI + Fabric",
    subtitle = glue("R²(full) = {round(r2_full, 4)}  |  R²(SPI only) = {round(r2_spi_only,4)}  |  R²(Fabric only) = {round(r2_fab_only,4)}")
  )

# ── 图 C：残差诊断三联图 ────────────────────────────────────────────────────

p_resid_fit <- ggplot(tibble(fitted = fitted_vec, resid = resid_vec),
                      aes(x = fitted, y = resid)) +
  geom_point(size = 0.9, alpha = 0.28, color = "#6B7280") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#EF4444", linewidth = 0.7) +
  geom_smooth(method = "loess", se = FALSE,
              color = "#374151", linewidth = 0.7, linetype = "dotted") +
  theme_bw(base_size = 9) +
  labs(title = "Residuals vs Fitted",
       x = "Fitted SPHARM distance", y = "Residual") +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))

p_resid_hist <- ggplot(tibble(resid = resid_std_vec), aes(x = resid)) +
  geom_histogram(bins = 40, fill = "#D4619A", alpha = 0.75, color = "white") +
  geom_vline(xintercept = c(-2, 0, 2),
             linetype   = c("dashed", "solid", "dashed"),
             color      = c("grey50", "#374151", "grey50"),
             linewidth  = 0.7) +
  theme_bw(base_size = 9) +
  labs(title = "Standardized Residual Distribution",
       x = "Standardized residual", y = "Count") +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))

# 残差 vs SPI 距离散点图（用残差着色，突出 SPHARM 独有结构）
p_resid_scatter <- ggplot(
  resid_pairs,
  aes(x = d_spi, y = d_spharm, color = resid_std)
) +
  geom_point(size = 1.1, alpha = 0.45) +
  geom_abline(slope     = b_spi,
              intercept = b_int + b_fab * mean(vec_fab),
              color = "grey30", linewidth = 0.7, linetype = "dashed") +
  scale_color_gradient2(
    low = "#3B82F6", mid = "grey85", high = "#EF4444",
    midpoint = 0, limits = c(-3, 3), oob = scales::squish,
    name = "Std.\nresidual"
  ) +
  theme_bw(base_size = 9) +
  labs(
    title    = "SPHARM vs SPI Distance, Colored by Residual",
    subtitle = "Red = SPHARM encodes extra differentiation | Blue = SPHARM underestimates",
    x = "SPI distance", y = "SPHARM distance"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 9, hjust = 0.5),
    plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "grey40")
  )

p_resid_panel <- (p_resid_fit | p_resid_hist | p_resid_scatter) +
  plot_annotation(
    title = "Residual Analysis: SPHARM ~ SPI + Fabric MRM",
    theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  )

# ── 组合 ─────────────────────────────────────────────────────────────────────

p_final <- (p_rotating | p_venn) / p_resid_panel +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(
    title    = "MRM Variance Decomposition: SPI / Fabric / SPHARM",
    subtitle = glue(
      "Multiple Regression on distance Matrices (Lichstein 2007)  |  ",
      "nperm = {NPERM}  |  *** p<0.001  ** p<0.01  * p<0.05  ns p≥0.05"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40")
    )
  )

# ==============================================================================
# 保存输出
# ==============================================================================

ggsave(
  here("analysis/output/figures/mrm_variance_decomp.png"),
  plot  = p_final,
  width = 14, height = 10,
  dpi   = 300, bg = "white"
)
cat("\n图已保存：mrm_variance_decomp.png\n")

# 汇总结果 CSV
bind_rows(
  summary_main %>%
    mutate(analysis = "rotating_dv",
           detail   = NA_character_) %>%
    select(analysis, dv, detail, r2, unique, r2_pval, r2_sig),
  tibble(
    analysis = "variance_partitioning",
    dv       = "SPHARM power",
    detail   = vp_table$component,
    r2       = vp_table$value,
    unique   = NA_real_,
    r2_pval  = NA_real_,
    r2_sig   = NA_character_
  )
) %>%
  write_csv(here("analysis/data/derived_data/mrm_results.csv"))
cat("结果已保存：mrm_results.csv\n")

# 残差矩阵 CSV（标本对级别）
resid_pairs %>%
  select(ID_i, ID_j, resid, resid_std, d_spharm, d_spi, d_fabric, fitted) %>%
  write_csv(here("analysis/data/derived_data/mrm_spharm_residuals.csv"))
cat("残差已保存：mrm_spharm_residuals.csv\n")

# ==============================================================================
# 结论摘要
# ==============================================================================

cat("\n====== 结论 ======\n\n")
cat("【主分析：轮换因变量 R²】\n\n")
for (i in seq_len(nrow(summary_main))) {
  r <- summary_main[i, ]
  cat(sprintf(
    "  %-14s：R² = %.4f（%s），%.1f%% 可被另两个方法解释，%.1f%% 独有\n",
    r$dv, r$r2, r$r2_sig, r$r2 * 100, r$unique * 100
  ))
}

cat(glue(
  "\n【SPHARM 方差分解】\n\n",
  "  全模型 R²         = {round(r2_full, 4)}（SPI + Fabric 联合可解释）\n",
  "    [a] SPI 独有    = {round(vp_spi*100, 1)}%\n",
  "    [b] Fabric 独有 = {round(vp_fab*100, 1)}%\n",
  "    [c] 共享        = {round(vp_shared*100, 1)}%\n",
  "    [d] SPHARM 独有（不可解释）= {round(vp_resid*100, 1)}%\n"
))

cat(glue(
  "\n【残差分析】\n\n",
  "  正残差对（SPHARM 比预测更远，潜在独有差异）：{n_pos} 对\n",
  "  负残差对：{n_neg} 对\n",
  "  |std.resid| > 2 的极端对：{n_out2} 对\n",
  "  完整残差矩阵已保存至 mrm_spharm_residuals.csv\n\n"
))

cat("【方法论说明】\n\n")
cat("  MRM 用置换检验处理距离矩阵非独立性（Lichstein 2007）\n")
cat("  使用 Pearson（mrank=FALSE），R² 解释等同标准 OLS\n")
cat("  若需 Spearman 一致性，可将 mrank=FALSE 改为 mrank=TRUE 重跑\n")
cat("  方差分解 [c] 共享项若为负（估计误差），已自动修正为 0\n")





# ==============================================================================
# 参数
# ==============================================================================

N_TOP_POS   <- 5    # 候选正残差对数（脚本自动筛选）
N_TOP_NEG   <- 3    # 候选负残差对数
N_SHOW_POS  <- 3    # 最终展示正残差案例数
N_SHOW_NEG  <- 1    # 最终展示负残差案例数

KDE_KAPPA   <- 15   # vMF 核函数集中度（越大越平滑）
KDE_GRID_N  <- 100  # 球面投影格点密度（每轴）

POWER_COLS  <- paste0("power_l", 1:5)

# 颜色方案
COL_POS_A  <- "#D4619A"    # 正残差案例：标本 A
COL_POS_B  <- "#7BAED4"    # 正残差案例：标本 B
COL_NEG_A  <- "#E8A838"    # 负残差案例：标本 A
COL_NEG_B  <- "#5BAD8F"    # 负残差案例：标本 B
COL_KDE_HI <- "#B02070"    # KDE 高密度颜色

# ==============================================================================
# 工具函数
# ==============================================================================

compute_R <- function(ux, uy, uz) {
  sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2) /
    sum(sqrt(ux^2 + uy^2 + uz^2))
}

compute_EI <- function(ux, uy, uz) {
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / nrow(U)
  lambda <- sort(eigen(T_mat, symmetric = TRUE)$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_)
  )
}

# ------------------------------------------------------------------------------
# vMF-KDE：Lambert 等面积方位投影上的球面密度估计
# 输入：n×3 单位向量矩阵，kappa = 集中度参数
# 输出：可供 ggplot geom_raster 使用的 tibble（x, y, density）
# ------------------------------------------------------------------------------
vmf_kde_lambert <- function(dirs, kappa = 15, grid_n = 100) {
  dirs <- as.matrix(dirs)
  dirs <- dirs / sqrt(rowSums(dirs^2))   # 确保单位化
  
  # Lambert 等面积方位投影格点（上半球 + 下半球折叠）
  # x, y ∈ [-1, 1]，对应球面上半球
  g    <- seq(-1, 1, length.out = grid_n)
  grid <- expand.grid(lx = g, ly = g) %>%
    filter(lx^2 + ly^2 <= 1)
  
  r2 <- grid$lx^2 + grid$ly^2
  lz <- sqrt(pmax(1 - r2, 0))   # 上半球
  
  # 球面格点 3D 坐标
  gx <- grid$lx * sqrt(2 - r2)
  gy <- grid$ly * sqrt(2 - r2)
  gz <- lz
  
  # vMF 核：density ∝ exp(kappa * u · v)，对所有数据点求和
  # 同时考虑下半球（轴向数据：u 和 -u 等价）
  density <- numeric(nrow(grid))
  for (i in seq_len(nrow(dirs))) {
    dot_pos <- gx * dirs[i, 1] + gy * dirs[i, 2] + gz * dirs[i, 3]
    dot_neg <- -dot_pos   # 对称轴（下半球折叠至上半球）
    density <- density + exp(kappa * pmax(dot_pos, dot_neg))
  }
  
  tibble(
    x       = grid$lx,
    y       = grid$ly,
    density = density / max(density)   # 归一化到 [0,1]
  )
}

# ------------------------------------------------------------------------------
# 单个标本的 KDE 球面密度图
# ------------------------------------------------------------------------------
plot_kde <- function(dirs, specimen_id, fill_hi = COL_KDE_HI,
                     kappa = KDE_KAPPA, grid_n = KDE_GRID_N,
                     subtitle_extra = "") {
  kde_df <- vmf_kde_lambert(dirs, kappa = kappa, grid_n = grid_n)
  
  # 外圆（投影边界）
  theta_circ <- seq(0, 2 * pi, length.out = 300)
  circ <- tibble(x = cos(theta_circ), y = sin(theta_circ))
  
  # 数据点投影（Lambert）
  dirs <- as.matrix(dirs)
  dirs <- dirs / sqrt(rowSums(dirs^2))
  # 折叠下半球到上半球
  dirs[dirs[, 3] < 0, ] <- -dirs[dirs[, 3] < 0, ]
  r2   <- dirs[, 1]^2 + dirs[, 2]^2
  lx   <- dirs[, 1] / sqrt(2 - r2 + 1e-10)
  ly   <- dirs[, 2] / sqrt(2 - r2 + 1e-10)
  pts  <- tibble(x = lx, y = ly)
  
  ggplot() +
    geom_raster(data = kde_df, aes(x, y, fill = density),
                interpolate = TRUE) +
    scale_fill_gradient(
      low = "white", high = fill_hi,
      limits = c(0, 1), guide = "none"
    ) +
    geom_path(data = circ, aes(x, y),
              color = "grey40", linewidth = 0.5) +
    geom_point(data = pts, aes(x, y),
               shape = 3, size = 1.5, color = "grey20", alpha = 0.8) +
    coord_equal(xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08)) +
    theme_void(base_size = 9) +
    labs(
      title    = specimen_id,
      subtitle = subtitle_extra
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 9,
                                   hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5,
                                   color = "grey45", margin = margin(b = 1))
    )
}

# ------------------------------------------------------------------------------
# 单对标本的 SPHARM 能量谱对比图
# ------------------------------------------------------------------------------
plot_power <- function(pow_a, pow_b, id_a, id_b,
                       col_a = COL_POS_A, col_b = COL_POS_B) {
  df <- tibble(
    l      = rep(1:5, 2),
    power  = c(pow_a, pow_b),
    spec   = rep(c(id_a, id_b), each = 5)
  ) %>%
    mutate(spec = factor(spec, levels = c(id_a, id_b)))
  
  ggplot(df, aes(x = factor(l), y = power, fill = spec)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.62, alpha = 0.88) +
    scale_fill_manual(values = c(col_a, col_b),
                      name = NULL) +
    theme_bw(base_size = 9) +
    theme(
      plot.title       = element_text(face = "bold", size = 9, hjust = 0.5),
      legend.position  = "bottom",
      legend.text      = element_text(size = 7.5),
      legend.key.size  = unit(0.45, "cm"),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = "SPHARM Power Spectrum (l1–l5)",
      x     = "Degree (l)",
      y     = "Power"
    )
}

# ------------------------------------------------------------------------------
# 数值对比表（SPI / E / I / SPHARM power）
# ------------------------------------------------------------------------------
make_table_grob <- function(row_a, row_b,
                            pow_a, pow_b,
                            id_a, id_b,
                            d_spi, d_fab, d_sph, resid_std) {
  fmt <- function(x) sprintf("%.4f", x)
  fmt2 <- function(x) sprintf("%.3f", x)
  
  tbl <- data.frame(
    Metric    = c("SPI", "Fabric E", "Fabric I",
                  "SPHARM l1", "SPHARM l2", "SPHARM l3",
                  "SPHARM l4", "SPHARM l5",
                  "─────────", "Δ SPI dist",
                  "Δ Fabric dist", "Δ SPHARM dist",
                  "Std. residual"),
    check.names = FALSE
  )
  tbl[[id_a]] <- c(
    fmt(row_a$SPI), fmt(row_a$E), fmt(row_a$I),
    fmt2(pow_a[1]), fmt2(pow_a[2]), fmt2(pow_a[3]),
    fmt2(pow_a[4]), fmt2(pow_a[5]),
    "─────", "","","",""
  )
  tbl[[id_b]] <- c(
    fmt(row_b$SPI), fmt(row_b$E), fmt(row_b$I),
    fmt2(pow_b[1]), fmt2(pow_b[2]), fmt2(pow_b[3]),
    fmt2(pow_b[4]), fmt2(pow_b[5]),
    "─────",
    fmt(d_spi), fmt(d_fab), fmt(d_sph),
    sprintf("%.3f", resid_std)
  )
  
  tt <- gridExtra::tableGrob(
    tbl, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size  = 7.5,
      core       = list(fg_params = list(cex = 0.82)),
      colhead    = list(fg_params = list(cex = 0.85, fontface = "bold"))
    )
  )
  tt
}

# ==============================================================================
# 读取数据
# ==============================================================================

exp_data <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 3)

directions_aligned <- read_csv(
  here("analysis/data/derived_data/directions_aligned_svd.csv"),
  show_col_types = FALSE)

spharm_all <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE)

residuals_df <- read_csv(
  here("analysis/data/derived_data/mrm_spharm_residuals.csv"),
  show_col_types = FALSE)

exp_ids        <- unique(exp_data$ID)
directions_exp <- directions_aligned %>% filter(ID %in% exp_ids)

rei_exp <- directions_exp %>%
  group_by(ID) %>%
  summarise(
    SPI = compute_R(ux, uy, uz),
    E   = compute_EI(ux, uy, uz)$E,
    I   = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

spharm_exp <- spharm_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, all_of(POWER_COLS))

specimen_meta <- rei_exp %>%
  left_join(spharm_exp, by = "ID") %>%
  filter(!is.na(power_l1))

# ==============================================================================
# 筛选候选标本对
# ==============================================================================

# 正残差候选：std_resid 最大
top_pos <- residuals_df %>%
  arrange(desc(resid_std)) %>%
  head(N_TOP_POS)

# 负残差候选：std_resid 最小（最负）
top_neg <- residuals_df %>%
  arrange(resid_std) %>%
  head(N_TOP_NEG)

cat("====== 正残差候选对（SPHARM-unique differentiation）======\n\n")
top_pos %>%
  select(ID_i, ID_j, resid_std, d_spharm, d_spi, d_fabric) %>%
  print()

cat("\n====== 负残差候选对（Undetected differentiation）======\n\n")
top_neg %>%
  select(ID_i, ID_j, resid_std, d_spharm, d_spi, d_fabric) %>%
  print()

# ==============================================================================
# 手动选定展示标本对
# ==============================================================================
# !! 根据上方打印结果手动填写，优先选择：
#    正残差：不同传统类型（如 Discoid vs Levallois），且 SPI/E/I 值接近
#    负残差：SPI/Fabric 距离大，但 SPHARM 距离小
# !! 若不确定，先跑脚本看候选输出，再填写下方列表
#
# 格式：list(i = "ID_A", j = "ID_B")

pairs_pos <- top_pos %>%
  slice(1:N_SHOW_POS) %>%
  pmap(function(ID_i, ID_j, ...) list(i = ID_i, j = ID_j))

pairs_neg <- top_neg %>%
  slice(1:N_SHOW_NEG) %>%
  pmap(function(ID_i, ID_j, ...) list(i = ID_i, j = ID_j))

# ==============================================================================
# 构造每对标本的面板
# ==============================================================================

make_pair_panel <- function(pair, resid_row,
                            col_a, col_b,
                            panel_title) {
  id_a <- pair$i
  id_b <- pair$j
  
  # 方向向量
  dirs_a <- directions_exp %>% filter(ID == id_a) %>%
    select(ux, uy, uz) %>% as.matrix()
  dirs_b <- directions_exp %>% filter(ID == id_b) %>%
    select(ux, uy, uz) %>% as.matrix()
  
  # 统计量
  row_a <- specimen_meta %>% filter(ID == id_a)
  row_b <- specimen_meta %>% filter(ID == id_b)
  
  pow_a <- as.numeric(row_a[, POWER_COLS])
  pow_b <- as.numeric(row_b[, POWER_COLS])
  
  subtitle_a <- glue("SPI={round(row_a$SPI,3)}  E={round(row_a$E,3)}  I={round(row_a$I,3)}")
  subtitle_b <- glue("SPI={round(row_b$SPI,3)}  E={round(row_b$E,3)}  I={round(row_b$I,3)}")
  
  # (1) KDE 图
  p_kde_a <- plot_kde(dirs_a, id_a, fill_hi = col_a,
                      subtitle_extra = subtitle_a)
  p_kde_b <- plot_kde(dirs_b, id_b, fill_hi = col_b,
                      subtitle_extra = subtitle_b)
  
  # (2) 能量谱图
  p_pow <- plot_power(pow_a, pow_b, id_a, id_b,
                      col_a = col_a, col_b = col_b)
  
  # (3) 数值表
  tbl_grob <- make_table_grob(
    row_a, row_b, pow_a, pow_b, id_a, id_b,
    d_spi      = resid_row$d_spi,
    d_fab      = resid_row$d_fabric,
    d_sph      = resid_row$d_spharm,
    resid_std  = resid_row$resid_std
  )
  
  # 组合面板：左 = 两个 KDE，右 = 能量谱 + 表格
  p_kde_row  <- (p_kde_a | p_kde_b) +
    plot_layout(widths = c(1, 1))
  
  p_right_top <- p_pow
  
  # 将 tableGrob 包装成 ggplot
  p_tbl <- ggplot() +
    annotation_custom(tbl_grob) +
    theme_void()
  
  p_right <- (p_right_top / p_tbl) +
    plot_layout(heights = c(1.2, 1))
  
  p_panel <- (p_kde_row | p_right) +
    plot_layout(widths = c(1.1, 1)) +
    plot_annotation(
      title = panel_title,
      theme = theme(
        plot.title    = element_text(face  = "bold", size = 10,
                                     hjust = 0.5,
                                     color = if (grepl("unique", panel_title, ignore.case = TRUE))
                                       "#8B1A4A" else "#2A5E8A"),
        plot.background = element_rect(
          fill  = if (grepl("unique", panel_title, ignore.case = TRUE))
            "#FFF5FA" else "#F5FAFF",
          color = if (grepl("unique", panel_title, ignore.case = TRUE))
            "#D4619A" else "#7BAED4",
          linewidth = 1
        )
      )
    )
  
  p_panel
}

# ==============================================================================
# 生成所有面板
# ==============================================================================

cat("\n生成正残差面板...\n")

panels_pos <- map2(
  pairs_pos,
  seq_along(pairs_pos),
  function(pair, idx) {
    resid_row <- residuals_df %>%
      filter((ID_i == pair$i & ID_j == pair$j) |
               (ID_i == pair$j & ID_j == pair$i)) %>%
      slice(1)
    
    make_pair_panel(
      pair, resid_row,
      col_a       = COL_POS_A,
      col_b       = COL_POS_B,
      panel_title = glue(
        "SPHARM-unique differentiation #{idx}  |  ",
        "{pair$i} vs {pair$j}  |  ",
        "Std. residual = +{round(resid_row$resid_std, 2)}"
      )
    )
  }
)

cat("生成负残差面板...\n")

panels_neg <- map2(
  pairs_neg,
  seq_along(pairs_neg),
  function(pair, idx) {
    resid_row <- residuals_df %>%
      filter((ID_i == pair$i & ID_j == pair$j) |
               (ID_i == pair$j & ID_j == pair$i)) %>%
      slice(1)
    
    make_pair_panel(
      pair, resid_row,
      col_a       = COL_NEG_A,
      col_b       = COL_NEG_B,
      panel_title = glue(
        "Phase-dependent differentiation undetected by power spectrum #{idx}  |  ",
        "{pair$i} vs {pair$j}  |  ",
        "Std. residual = {round(resid_row$resid_std, 2)}"
      )
    )
  }
)

# ==============================================================================
# 组合全图
# ==============================================================================

all_panels <- c(panels_pos, panels_neg)

p_final <- wrap_plots(all_panels, ncol = 1) +
  plot_annotation(
    title    = "MRM Residual Case Studies: What SPHARM Power Spectrum Captures Beyond SPI and Fabric",
    subtitle = paste(
      glue("Top {N_SHOW_POS} positive-residual pairs (SPHARM sees more differentiation than SPI/Fabric predict)"),
      glue("+ {N_SHOW_NEG} negative-residual pair (control: SPHARM sees less differentiation)"),
      sep = "  |  "
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
      plot.margin   = margin(10, 10, 10, 10)
    )
  )

# ==============================================================================
# 保存
# ==============================================================================

n_panels  <- length(all_panels)
fig_h     <- 4.8 * n_panels + 0.8

ggsave(
  here("analysis/output/figures/mrm_case_study_panels.png"),
  plot   = p_final,
  width  = 14,
  height = fig_h,
  dpi    = 300,
  bg     = "white"
)
cat(glue("\n图已保存：mrm_case_study_panels.png  ({n_panels} 个面板，高度 {round(fig_h,1)} 英寸)\n"))

# ==============================================================================
# 打印候选对完整统计（辅助手动复核）
# ==============================================================================

cat("\n====== 展示对完整统计 ======\n\n")

walk(c(pairs_pos, pairs_neg), function(pair) {
  row_a <- specimen_meta %>% filter(ID == pair$i)
  row_b <- specimen_meta %>% filter(ID == pair$j)
  resid_row <- residuals_df %>%
    filter((ID_i == pair$i & ID_j == pair$j) |
             (ID_i == pair$j & ID_j == pair$i)) %>%
    slice(1)
  
  cat(glue("── {pair$i} vs {pair$j}  (std.resid = {round(resid_row$resid_std, 3)})\n"))
  cat(glue("   {pair$i}:  SPI={round(row_a$SPI,4)}  E={round(row_a$E,4)}  I={round(row_a$I,4)}\n"))
  cat(glue("   {pair$j}:  SPI={round(row_b$SPI,4)}  E={round(row_b$E,4)}  I={round(row_b$I,4)}\n"))
  cat(glue("   d_SPI={round(resid_row$d_spi,4)}  d_Fabric={round(resid_row$d_fabric,4)}  d_SPHARM={round(resid_row$d_spharm,4)}\n\n"))
})

