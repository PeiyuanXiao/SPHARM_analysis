# ==============================================================================
# mrm_analysis_and_visualization.R
# 方法验证：MRM 回归分析 + 案例可视化
#
# 框架：Lichstein (2007) Oecologia — Multiple Regression on distance Matrices
#
# ── 模块一：MRM 回归分析 ──────────────────────────────────────────────────────
#   第一步：三组 MRM，轮换因变量
#           → R²（他人可解释比例）/ 1-R²（独有信息比例）/ 回归系数
#   第二步：SPHARM 为因变量的完整方差分解
#           → [a] SPI 独有 / [b] Fabric 独有 / [c] 共享 / [d] SPHARM 独有
#   第三步：提取 SPHARM ~ SPI + Fabric 的残差矩阵
#           → 按残差排序，识别 SPHARM 独有差异案例
#   可视化：R² 柱状图 + 方差分解 Venn 图 + 残差诊断三联图
#
# ── 模块二：案例可视化 ────────────────────────────────────────────────────────
#   对选定标本对，分别输出：
#   (1) KDE 球面密度图（Lambert 等面积投影，双标本并排）
#   (2) SPI / E / I 哑铃图 + SPHARM 能量谱折线图（拼合）
#   数值统计打印到控制台
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 3）
#   - analysis/data/derived_data/directions_aligned_svd.csv
#   - analysis/data/derived_data/SPHARM_direction.csv
#
# 输出：
#   - analysis/output/figures/mrm_variance_decomp.png
#   - analysis/output/figures/mrm_case_{idx}_..._KDE.png
#   - analysis/output/figures/mrm_case_{idx}_..._metrics.png
#   - analysis/data/derived_data/mrm_results.csv
#   - analysis/data/derived_data/mrm_spharm_residuals.csv
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(ecodist)       # MRM()
library(patchwork)
library(glue)
library(colorspace)    # lighten() / darken()

# ==============================================================================
# 全局参数
# ==============================================================================

NPERM      <- 999
KDE_KAPPA  <- 15
KDE_GRID_N <- 100
POWER_COLS <- paste0("power_l", 1:5)
COL_POS_A  <- "#D4619A"
COL_POS_B  <- "#7BAED4"

set.seed(42)

# ==============================================================================
# 公共工具函数
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

if (!all(POWER_COLS %in% colnames(spharm_all))) {
  stop("SPHARM 文件中缺少 power_l1~power_l5 列。")
}

spharm_exp <- spharm_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, all_of(POWER_COLS))

df <- rei_exp %>%
  left_join(spharm_exp, by = "ID") %>%
  filter(!is.na(SPI), !is.na(E), !is.na(I), !is.na(power_l1)) %>%
  arrange(ID)

specimen_meta <- df   # 供模块二使用

n   <- nrow(df)
ids <- df$ID

cat(sprintf("进入分析的标本数：%d\n\n", n))

# ==============================================================================
# ████████████████████████████████████████████████████████████████████████████
# 模块一：MRM 回归分析
# ████████████████████████████████████████████████████████████████████████████
# ==============================================================================

D_spi <- dist_euclidean_scaled(df %>% select(SPI))
D_fab <- dist_euclidean_scaled(df %>% select(E, I))
D_sph <- dist_euclidean_scaled(df %>% select(all_of(POWER_COLS)))

# ------------------------------------------------------------------------------
# 第一步：三组 MRM，轮换因变量
# ------------------------------------------------------------------------------

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

extract_mrm_summary <- function(mrm_obj, dv_label) {
  r2       <- as.numeric(mrm_obj$r.squared["R2"])
  r2_p     <- as.numeric(mrm_obj$r.squared["pval"])
  coef_mat <- mrm_obj$coef
  pred_rows <- coef_mat[-1, , drop = FALSE]
  coef_tbl  <- tibble(
    dv        = dv_label,
    predictor = rownames(pred_rows),
    beta      = as.numeric(pred_rows[, 1]),
    p_beta    = as.numeric(pred_rows[, 2]),
    sig_beta  = make_sig(as.numeric(pred_rows[, 2]))
  )
  list(
    r2 = r2, unique = 1 - r2,
    r2_pval = r2_p, r2_sig = make_sig(r2_p),
    dv = dv_label, coef_tbl = coef_tbl
  )
}

res_sph <- extract_mrm_summary(mrm_sph, "SPHARM power")
res_spi <- extract_mrm_summary(mrm_spi, "SPI")
res_fab <- extract_mrm_summary(mrm_fab, "Fabric")

summary_main <- tibble(
  dv      = c("SPHARM power", "SPI", "Fabric"),
  r2      = c(res_sph$r2,      res_spi$r2,      res_fab$r2),
  unique  = c(res_sph$unique,  res_spi$unique,  res_fab$unique),
  r2_pval = c(res_sph$r2_pval, res_spi$r2_pval, res_fab$r2_pval),
  r2_sig  = c(res_sph$r2_sig,  res_spi$r2_sig,  res_fab$r2_sig)
)

coef_main <- bind_rows(res_sph$coef_tbl, res_spi$coef_tbl, res_fab$coef_tbl)

cat("R²（另两个方法联合可解释的方差比例）：\n\n")
cat(sprintf("  %-22s  %7s  %7s  %8s  %s\n", "因变量", "R²", "1-R²", "p(R²)", "sig"))
cat("  ", paste(rep("-", 58), collapse = ""), "\n", sep = "")
for (i in seq_len(nrow(summary_main))) {
  r <- summary_main[i, ]
  cat(sprintf("  %-22s  %7.4f  %7.4f  %8.4f  %s\n",
              r$dv, r$r2, r$unique, r$r2_pval, r$r2_sig))
}

cat("\n回归系数：\n\n")
cat(sprintf("  %-22s  %-25s  %8s  %8s  %s\n", "因变量", "预测变量", "β", "p(β)", "sig"))
cat("  ", paste(rep("-", 76), collapse = ""), "\n", sep = "")
for (i in seq_len(nrow(coef_main))) {
  r <- coef_main[i, ]
  cat(sprintf("  %-22s  %-25s  %8.4f  %8.4f  %s\n",
              r$dv, r$predictor, r$beta, r$p_beta, r$sig_beta))
}

# ------------------------------------------------------------------------------
# 第二步：SPHARM 完整方差分解
# ------------------------------------------------------------------------------

cat("\n======================================================\n")
cat(" 第二步：SPHARM 完整方差分解\n")
cat("======================================================\n\n")

set.seed(42); mrm_sph_spi_only <- MRM(D_sph ~ D_spi, nperm = NPERM, mrank = FALSE)
set.seed(42); mrm_sph_fab_only <- MRM(D_sph ~ D_fab, nperm = NPERM, mrank = FALSE)

r2_full     <- as.numeric(mrm_sph$r.squared["R2"])
r2_spi_only <- as.numeric(mrm_sph_spi_only$r.squared["R2"])
r2_fab_only <- as.numeric(mrm_sph_fab_only$r.squared["R2"])

vp_spi    <- max(r2_full - r2_fab_only, 0)
vp_fab    <- max(r2_full - r2_spi_only, 0)
vp_shared <- max(r2_spi_only + r2_fab_only - r2_full, 0)
vp_resid  <- 1 - r2_full

cat(sprintf("  R²(SPHARM ~ SPI only)    = %.4f\n", r2_spi_only))
cat(sprintf("  R²(SPHARM ~ Fabric only) = %.4f\n", r2_fab_only))
cat(sprintf("  R²(SPHARM ~ SPI+Fabric)  = %.4f  [全模型]\n\n", r2_full))

vp_table <- tibble(
  component = c("[a] Unique to SPI", "[b] Unique to Fabric",
                "[c] Shared (SPI ∩ Fabric)", "[d] Unexplained — SPHARM unique"),
  value = c(vp_spi, vp_fab, vp_shared, vp_resid),
  pct   = c(vp_spi, vp_fab, vp_shared, vp_resid) * 100
)

cat("方差分解：\n\n")
for (i in seq_len(nrow(vp_table))) {
  cat(sprintf("  %-38s : %.4f  (%.1f%%)\n",
              vp_table$component[i], vp_table$value[i], vp_table$pct[i]))
}
cat(sprintf("\n  总计（应 = 1.0000）：%.4f\n", sum(vp_table$value)))

# ------------------------------------------------------------------------------
# 第三步：残差矩阵提取
# ------------------------------------------------------------------------------

cat("\n======================================================\n")
cat(" 第三步：残差矩阵提取（SPHARM ~ SPI + Fabric）\n")
cat("======================================================\n\n")

coef_raw <- mrm_sph$coef[, 1]
b_int    <- coef_raw[1]
b_spi    <- coef_raw[2]
b_fab    <- coef_raw[3]

cat(sprintf("  截距 β₀ = %.4f\n", b_int))
cat(sprintf("  β(SPI)  = %.4f  (p = %.4f  %s)\n",
            b_spi, mrm_sph$coef["D_spi", 2], make_sig(mrm_sph$coef["D_spi", 2])))
cat(sprintf("  β(Fab)  = %.4f  (p = %.4f  %s)\n\n",
            b_fab, mrm_sph$coef["D_fab", 2], make_sig(mrm_sph$coef["D_fab", 2])))

vec_sph       <- as.vector(D_sph)
vec_spi       <- as.vector(D_spi)
vec_fab       <- as.vector(D_fab)
fitted_vec    <- b_int + b_spi * vec_spi + b_fab * vec_fab
resid_vec     <- vec_sph - fitted_vec
resid_sd      <- sd(resid_vec)
resid_std_vec <- resid_vec / resid_sd

cat(sprintf("  残差均值：%.6f\n  残差标准差：%.4f\n  残差范围：[%.4f, %.4f]\n\n",
            mean(resid_vec), resid_sd, min(resid_vec), max(resid_vec)))

resid_mat <- matrix(0, n, n)
resid_mat[lower.tri(resid_mat)] <- resid_vec
resid_mat <- resid_mat + t(resid_mat)
rownames(resid_mat) <- colnames(resid_mat) <- ids

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

n_pos  <- sum(resid_pairs$resid > 0)
n_neg  <- sum(resid_pairs$resid < 0)
n_out2 <- sum(abs(resid_pairs$resid_std) > 2)

cat(sprintf(
  "  正残差对：%d  负残差对：%d  |std.resid|>2 的极端对：%d\n\n",
  n_pos, n_neg, n_out2
))

# ------------------------------------------------------------------------------
# 模块一可视化
# ------------------------------------------------------------------------------

dv_order <- c("SPI", "Fabric", "SPHARM power")
bar_pal  <- c(
  "Explained by other two (R²)" = "#7BAED4",
  "Unique / unexplained (1−R²)" = "#E8C4D8"
)

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
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle      = element_text(size = 8, hjust = 0.5, color = "grey40"),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title    = "MRM: Information Uniqueness — Rotating Dependent Variable",
    subtitle = glue("R² = variance explained by the other two methods  |  1−R² = unique / unexplained  |  nperm = {NPERM}"),
    x = "Dependent Variable", y = "Proportion of Variance"
  )

# Venn 图
theta   <- seq(0, 2 * pi, length.out = 300)
r_base  <- 0.32
sep_adj <- pmin(pmax(r_base * (1 - vp_shared / max(r2_full, 0.01)), 0.05), 0.30)
cx1 <- -sep_adj; cx2 <- sep_adj

circ1 <- tibble(x = r_base * cos(theta) + cx1, y = r_base * sin(theta), grp = "SPI")
circ2 <- tibble(x = r_base * cos(theta) + cx2, y = r_base * sin(theta), grp = "Fabric")

p_venn <- ggplot() +
  geom_polygon(data = circ1, aes(x, y),
               fill = "#A1C2E6", alpha = 0.45, color = "#5B92C4", linewidth = 0.8) +
  geom_polygon(data = circ2, aes(x, y),
               fill = "#FFBAE0", alpha = 0.45, color = "#D472A8", linewidth = 0.8) +
  annotate("text", x = cx1 - r_base * 0.4, y = 0,
           label = glue("SPI unique\n[a] {round(vp_spi*100,1)}%"),
           size = 3.5, fontface = "bold", color = "#3A7BC8") +
  annotate("text", x = cx2 + r_base * 0.4, y = 0,
           label = glue("Fabric unique\n[b] {round(vp_fab*100,1)}%"),
           size = 3.5, fontface = "bold", color = "#C4528A") +
  annotate("text", x = 0, y = 0,
           label = glue("Shared\n[c] {round(vp_shared*100,1)}%"),
           size = 3.2, fontface = "bold", color = "grey25") +
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
    subtitle = glue("R²(full)={round(r2_full,4)}  |  R²(SPI only)={round(r2_spi_only,4)}  |  R²(Fabric only)={round(r2_fab_only,4)}")
  )

# 残差诊断
p_resid_fit <- ggplot(tibble(fitted = fitted_vec, resid = resid_vec),
                      aes(x = fitted, y = resid)) +
  geom_point(size = 0.9, alpha = 0.28, color = "#6B7280") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#EF4444", linewidth = 0.7) +
  geom_smooth(method = "loess", se = FALSE,
              color = "#374151", linewidth = 0.7, linetype = "dotted") +
  theme_bw(base_size = 9) +
  labs(title = "Residuals vs Fitted",
       x = "Fitted SPHARM distance", y = "Residual") +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))

p_resid_hist <- ggplot(tibble(resid = resid_std_vec), aes(x = resid)) +
  geom_histogram(bins = 40, fill = "#D4619A", alpha = 0.75, color = "white") +
  geom_vline(xintercept = c(-2, 0, 2),
             linetype = c("dashed", "solid", "dashed"),
             color    = c("grey50", "#374151", "grey50"),
             linewidth = 0.7) +
  theme_bw(base_size = 9) +
  labs(title = "Standardized Residual Distribution",
       x = "Standardized residual", y = "Count") +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))

p_resid_scatter <- ggplot(resid_pairs,
                          aes(x = d_spi, y = d_spharm, color = resid_std)) +
  geom_point(size = 1.1, alpha = 0.45) +
  geom_abline(slope = b_spi,
              intercept = b_int + b_fab * mean(vec_fab),
              color = "grey30", linewidth = 0.7, linetype = "dashed") +
  scale_color_gradient2(low = "#3B82F6", mid = "grey85", high = "#EF4444",
                        midpoint = 0, limits = c(-3, 3), oob = scales::squish,
                        name = "Std.\nresidual") +
  theme_bw(base_size = 9) +
  labs(title = "SPHARM vs SPI Distance, Colored by Residual",
       subtitle = "Red = SPHARM encodes extra differentiation | Blue = SPHARM underestimates",
       x = "SPI distance", y = "SPHARM distance") +
  theme(
    plot.title    = element_text(face = "bold", size = 9, hjust = 0.5),
    plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "grey40")
  )

p_resid_panel <- (p_resid_fit | p_resid_hist | p_resid_scatter) +
  plot_annotation(
    title = "Residual Analysis: SPHARM ~ SPI + Fabric MRM",
    theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  )

p_mrm_final <- (p_rotating | p_venn) / p_resid_panel +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(
    title    = "MRM Variance Decomposition: SPI / Fabric / SPHARM",
    subtitle = glue("Multiple Regression on distance Matrices (Lichstein 2007)  |  nperm = {NPERM}  |  *** p<0.001  ** p<0.01  * p<0.05  ns p≥0.05"),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave(here("analysis/output/figures/mrm_variance_decomp.png"),
       plot = p_mrm_final, width = 14, height = 10, dpi = 300, bg = "white")
cat("图已保存：mrm_variance_decomp.png\n")

# 保存 CSV
bind_rows(
  summary_main %>%
    mutate(analysis = "rotating_dv", detail = NA_character_) %>%
    select(analysis, dv, detail, r2, unique, r2_pval, r2_sig),
  tibble(
    analysis = "variance_partitioning", dv = "SPHARM power",
    detail = vp_table$component, r2 = vp_table$value,
    unique = NA_real_, r2_pval = NA_real_, r2_sig = NA_character_
  )
) %>% write_csv(here("analysis/data/derived_data/mrm_results.csv"))

resid_pairs %>%
  select(ID_i, ID_j, resid, resid_std, d_spharm, d_spi, d_fabric, fitted) %>%
  write_csv(here("analysis/data/derived_data/mrm_spharm_residuals.csv"))

cat("结果已保存：mrm_results.csv\n残差已保存：mrm_spharm_residuals.csv\n\n")

# 结论摘要
cat("====== 结论 ======\n\n")
for (i in seq_len(nrow(summary_main))) {
  r <- summary_main[i, ]
  cat(sprintf("  %-14s：R² = %.4f（%s），%.1f%% 可被另两个方法解释，%.1f%% 独有\n",
              r$dv, r$r2, r$r2_sig, r$r2 * 100, r$unique * 100))
}
cat(glue("\n  全模型 R²={round(r2_full,4)}  [a]SPI独有={round(vp_spi*100,1)}%  [b]Fabric独有={round(vp_fab*100,1)}%  [c]共享={round(vp_shared*100,1)}%  [d]SPHARM独有={round(vp_resid*100,1)}%\n\n"))

# ==============================================================================
# ████████████████████████████████████████████████████████████████████████████
# 模块二：案例可视化
# ████████████████████████████████████████████████████████████████████████████
# ==============================================================================

residuals_df <- read_csv(
  here("analysis/data/derived_data/mrm_spharm_residuals.csv"),
  show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 可视化工具函数
# ------------------------------------------------------------------------------

vmf_kde_lambert <- function(dirs, kappa = 15, grid_n = 100) {
  dirs <- as.matrix(dirs)
  dirs <- dirs / sqrt(rowSums(dirs^2))
  g    <- seq(-1, 1, length.out = grid_n)
  grid <- expand.grid(lx = g, ly = g) %>% filter(lx^2 + ly^2 <= 1)
  r2   <- grid$lx^2 + grid$ly^2
  lz   <- sqrt(pmax(1 - r2, 0))
  gx   <- grid$lx * sqrt(2 - r2)
  gy   <- grid$ly * sqrt(2 - r2)
  gz   <- lz
  density <- numeric(nrow(grid))
  for (i in seq_len(nrow(dirs))) {
    dot_pos <- gx * dirs[i,1] + gy * dirs[i,2] + gz * dirs[i,3]
    density  <- density + exp(kappa * pmax(dot_pos, -dot_pos))
  }
  tibble(x = grid$lx, y = grid$ly, density = density / max(density))
}

plot_kde <- function(dirs, specimen_id,
                     fill_hi = COL_POS_A,
                     kappa   = KDE_KAPPA,
                     grid_n  = KDE_GRID_N) {
  kde_df <- vmf_kde_lambert(dirs, kappa = kappa, grid_n = grid_n)
  circ   <- tibble(x = cos(seq(0, 2*pi, length.out = 300)),
                   y = sin(seq(0, 2*pi, length.out = 300)))
  col_light <- colorspace::lighten(fill_hi, 0.6)
  col_deep  <- colorspace::darken(fill_hi,  0.35)
  
  ggplot() +
    geom_raster(data = kde_df, aes(x, y, fill = density), interpolate = TRUE) +
    scale_fill_gradientn(
      colors = c("white", col_light, fill_hi, col_deep),
      values = c(0, 0.25, 0.65, 1), limits = c(0, 1), guide = "none"
    ) +
    geom_contour(data = kde_df, aes(x, y, z = density),
                 color = "white", alpha = 0.6, linewidth = 0.35,
                 breaks = c(0.2, 0.4, 0.6, 0.8)) +
    geom_path(data = circ, aes(x, y), color = "grey30", linewidth = 0.9) +
    annotate("point", x = 0, y = 0, shape = 3, size = 2.8,
             color = "grey45", stroke = 0.8) +
    coord_equal(xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08)) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5,
                                margin = margin(t = 6, b = 6))
    ) +
    labs(title = specimen_id)
}

# SPI / E / I 哑铃图
plot_sei <- function(row_a, row_b, id_a, id_b,
                     col_a = COL_POS_A, col_b = COL_POS_B) {
  df_sei <- tibble(
    metric = factor(c("SPI", "Fabric E", "Fabric I"),
                    levels = c("Fabric I", "Fabric E", "SPI")),
    val_a  = c(as.numeric(row_a$SPI), as.numeric(row_a$E), as.numeric(row_a$I)),
    val_b  = c(as.numeric(row_b$SPI), as.numeric(row_b$E), as.numeric(row_b$I))
  )
  
  ggplot(df_sei) +
    geom_segment(aes(x = val_a, xend = val_b, y = metric, yend = metric),
                 color = "grey70", linewidth = 1.2) +
    geom_point(aes(x = val_a, y = metric), color = col_a, size = 4.5) +
    geom_point(aes(x = val_b, y = metric), color = col_b, size = 4.5) +
    geom_text(aes(x = val_a, y = metric, label = round(val_a, 3)),
              vjust = -1, size = 3, color = col_a) +
    geom_text(aes(x = val_b, y = metric, label = round(val_b, 3)),
              vjust = -1, size = 3, color = col_b) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.05)) +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.title.y       = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey85")
    ) +
    labs(title = "SPI & Fabric (E / I)", x = "Value")
}

# SPHARM 能量谱折线图
plot_power <- function(pow_a, pow_b, id_a, id_b,
                       col_a = COL_POS_A, col_b = COL_POS_B) {
  tibble(
    l     = rep(1:5, 2),
    power = c(pow_a, pow_b),
    spec  = rep(c(id_a, id_b), each = 5)
  ) %>%
    mutate(spec = factor(spec, levels = c(id_a, id_b))) %>%
    ggplot(aes(x = l, y = power, color = spec, group = spec)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    scale_color_manual(values = c(col_a, col_b), name = NULL) +
    scale_x_continuous(breaks = 1:5, labels = paste0("l", 1:5)) +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position    = "bottom",
      legend.text        = element_text(size = 8),
      legend.key.size    = unit(0.45, "cm"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(title = "SPHARM Power Spectrum (l1–l5)",
         x = "Degree", y = "Power")
}

# 构造单对面板
make_pair_panel <- function(pair, resid_row, col_a, col_b, panel_title) {
  id_a <- pair$i
  id_b <- pair$j
  
  dirs_a <- directions_exp %>% filter(ID == id_a) %>%
    select(ux, uy, uz) %>% as.matrix()
  dirs_b <- directions_exp %>% filter(ID == id_b) %>%
    select(ux, uy, uz) %>% as.matrix()
  
  row_a  <- specimen_meta %>% filter(ID == id_a)
  row_b  <- specimen_meta %>% filter(ID == id_b)
  pow_a  <- as.numeric(row_a[, POWER_COLS])
  pow_b  <- as.numeric(row_b[, POWER_COLS])
  
  # 控制台打印
  sep <- paste(rep("=", 60), collapse = "")
  dsh <- paste(rep("-", 60), collapse = "")
  cat(glue("\n{sep}\n{id_a}  vs  {id_b}\nStd. residual = {round(as.numeric(resid_row$resid_std), 3)}\n{dsh}\n"))
  cat(sprintf("%-14s  %-12s  %-12s\n", "Metric", id_a, id_b))
  cat(sprintf("%-14s  %-12.4f  %-12.4f\n", "SPI",      row_a$SPI, row_b$SPI))
  cat(sprintf("%-14s  %-12.4f  %-12.4f\n", "Fabric E", row_a$E,   row_b$E))
  cat(sprintf("%-14s  %-12.4f  %-12.4f\n", "Fabric I", row_a$I,   row_b$I))
  for (k in 1:5) {
    cat(sprintf("%-14s  %-12.4f  %-12.4f\n", paste0("SPHARM l", k), pow_a[k], pow_b[k]))
  }
  cat(dsh, "\n")
  cat(sprintf("%-14s  %-12.4f\n", "d_SPI",    as.numeric(resid_row$d_spi)))
  cat(sprintf("%-14s  %-12.4f\n", "d_Fabric", as.numeric(resid_row$d_fabric)))
  cat(sprintf("%-14s  %-12.4f\n", "d_SPHARM", as.numeric(resid_row$d_spharm)))
  cat(sep, "\n\n")
  
  list(
    kde_a     = plot_kde(dirs_a, id_a, fill_hi = col_a),
    kde_b     = plot_kde(dirs_b, id_b, fill_hi = col_b),
    sei       = plot_sei(row_a, row_b, id_a, id_b, col_a, col_b),
    power     = plot_power(pow_a, pow_b, id_a, id_b, col_a, col_b),
    title     = panel_title,
    resid_row = resid_row
  )
}

# ------------------------------------------------------------------------------
# 预览候选对
# ------------------------------------------------------------------------------

N_PREVIEW <- 20
cat(glue("====== 正残差候选对（前 {N_PREVIEW} ）======\n\n"))

preview_df <- residuals_df %>% arrange(desc(resid_std)) %>% head(N_PREVIEW)

walk(seq_len(nrow(preview_df)), function(k) {
  r     <- preview_df[k, ]
  row_a <- specimen_meta %>% filter(ID == r$ID_i)
  row_b <- specimen_meta %>% filter(ID == r$ID_j)
  cat(sprintf(
    "#%-2d  %-34s vs %-34s\n     std.resid=%6.3f  d_SPHARM=%.4f  d_SPI=%.4f  d_Fabric=%.4f\n     %-34s  SPI=%.4f  E=%.4f  I=%.4f\n     %-34s  SPI=%.4f  E=%.4f  I=%.4f\n\n",
    k, r$ID_i, r$ID_j,
    r$resid_std, r$d_spharm, r$d_spi, r$d_fabric,
    r$ID_i, row_a$SPI, row_a$E, row_a$I,
    r$ID_j, row_b$SPI, row_b$E, row_b$I
  ))
})

# ------------------------------------------------------------------------------
# 手动填写展示对
# ------------------------------------------------------------------------------

pairs_pos <- list(
  list(i = "EXP18_Discoid",                    j = "EXP29_Levallois preferential"),
  list(i = "EXP19_Levallois preferential",     j = "EXP27_Levallois recurrent"),
  list(i = "EXP11_Cylindrical bidirectional",  j = "EXP20_Cylindrical bidirectional")
)

# ------------------------------------------------------------------------------
# 生成面板并保存
# ------------------------------------------------------------------------------

cat("\n生成面板...\n\n")

panels_pos <- map2(
  pairs_pos,
  seq_along(pairs_pos),
  function(pair, idx) {
    resid_row <- residuals_df %>%
      filter((ID_i == pair$i & ID_j == pair$j) |
               (ID_i == pair$j & ID_j == pair$i)) %>%
      slice(1)
    
    make_pair_panel(
      pair        = pair,
      resid_row   = resid_row,
      col_a       = COL_POS_A,
      col_b       = COL_POS_B,
      panel_title = glue(
        "SPHARM-unique differentiation #{idx}  |  ",
        "{pair$i} vs {pair$j}  |  ",
        "Std. residual = +{round(as.numeric(resid_row$resid_std), 2)}"
      )
    )
  }
)

walk2(panels_pos, seq_along(panels_pos), function(panel, idx) {
  pair <- pairs_pos[[idx]]
  stem <- glue("mrm_case_{idx}_{pair$i}_vs_{pair$j}") %>%
    str_replace_all("[^A-Za-z0-9_\\-]", "_") %>%
    str_replace_all("_+", "_") %>%
    str_trunc(120, ellipsis = "")
  
  # KDE 双图（两标本并排）
  p_kde <- (panel$kde_a | panel$kde_b) +
    plot_layout(widths = c(1, 1))
  
  ggsave(
    here("analysis/output/figures", glue("{stem}_KDE.png")),
    plot = p_kde, width = 8, height = 4.8, dpi = 300, bg = "white"
  )
  
  # SPI/E/I 哑铃图 + 能量谱折线图（拼合）
  p_metrics <- (panel$sei | panel$power) +
    plot_layout(widths = c(1, 1.4))
  
  ggsave(
    here("analysis/output/figures", glue("{stem}_metrics.png")),
    plot = p_metrics, width = 9, height = 4.2, dpi = 300, bg = "white"
  )
  
  cat(glue("  [{idx}] 已保存：{stem}_KDE.png  +  {stem}_metrics.png\n\n"))
})

cat("全部完成。\n")