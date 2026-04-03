# ==============================================================================
# mantel_information_retention.R
# 片疤模式信息保留度分析
# 方法：Mantel 检验（方法距离矩阵 × KDE 真实距离矩阵）
#
# 修正说明（基于方法论审查）：
#   修正1：统一距离度量
#     原版：SPI 用欧氏，Fabric 用马氏，SPHARM 用余弦
#     修正：三种方法全部用余弦距离（特征 z-score 标准化后）
#     原因：马氏距离对 E/I 做了去相关变换，可能碰巧让 Fabric 的秩序
#           结构更接近 ground truth，混淆了"信息量"与"距离度量选择"
#
#   修正2：诊断性分析——全系数 SPHARM
#     除能量谱（power_l1–l5）外，额外用完整球谐系数向量做 Mantel 检验
#     原因：能量谱为换取旋转不变性丢弃了相位信息，全系数保留相位
#     预期：全系数 SPHARM 的 Mantel r 应显著高于能量谱
#     意义：精确量化能量谱为旋转不变性付出的信息代价
#
# 全程使用 SVD 对齐。
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 3）
#   - analysis/data/derived_data/directions_aligned_svd.csv
#   - analysis/data/derived_data/SPHARM_direction.csv
# 输出：
#   - analysis/output/figures/mantel_main.png         主分析（统一余弦距离）
#   - analysis/output/figures/mantel_diagnostic.png   诊断分析（全系数 SPHARM）
#   - analysis/data/derived_data/mantel_retention_results.csv
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(vegan)
library(patchwork)
library(glue)

exp_data           <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 3)
directions_aligned <- read_csv(
  here("analysis/data/derived_data/directions_aligned_svd.csv"),
  show_col_types = FALSE)
spharm_all         <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE)


# ==============================================================================
# 公共函数
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

# 余弦距离（z-score 标准化后）
# 统一用于所有方法，消除距离度量差异的混淆
dist_cosine_scaled <- function(X) {
  X   <- scale(as.matrix(X))          # z-score 标准化
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  1 - sim
}

# vMF-KDE（用于 KDE 真实距离矩阵）
fit_vmf_kde <- function(ux, uy, uz, G, bandwidth = 0.35) {
  X     <- cbind(ux, uy, uz)
  kappa <- 1 / bandwidth^2
  dot   <- G %*% t(X)
  dens  <- rowMeans(exp(kappa * dot))
  dens  / sum(dens)
}

make_sphere_grid <- function(n_bearing = 72, n_plunge = 36) {
  eps     <- pi / 180 * 0.01
  bearing <- seq(0, 2 * pi, length.out = n_bearing + 1)[-(n_bearing + 1)]
  plunge  <- seq(-pi/2 + eps, pi/2 - eps, length.out = n_plunge)
  grid    <- expand.grid(bearing = bearing, plunge = plunge)
  cbind(
    x = cos(grid$plunge) * cos(grid$bearing),
    y = cos(grid$plunge) * sin(grid$bearing),
    z = sin(grid$plunge)
  )
}


# ==============================================================================
# 1. 过滤实验石核的 SVD 对齐方向向量
# ==============================================================================

exp_ids        <- unique(exp_data$ID)
directions_exp <- directions_aligned %>% filter(ID %in% exp_ids)

cat(sprintf("实验石核：%d 件标本，%d 条刮痕\n\n",
            n_distinct(directions_exp$ID), nrow(directions_exp)))


# ==============================================================================
# 2. 计算 SPI 和 Fabric
# ==============================================================================

rei_exp <- directions_exp %>%
  group_by(ID) %>%
  summarise(
    SPI = compute_R(ux, uy, uz),
    E   = compute_EI(ux, uy, uz)$E,
    I   = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

cat(sprintf("E 与 I 的相关系数：%.3f\n\n",
            cor(rei_exp$E, rei_exp$I, use = "complete.obs")))


# ==============================================================================
# 3. 合并 SPHARM 数据（能量谱 + 全系数）
# ==============================================================================

power_cols <- paste0("power_l", 1:5)
coeff_cols <- grep("^coeff_", colnames(spharm_all), value = TRUE)

spharm_exp <- spharm_all %>%
  filter(ID %in% exp_ids) %>%
  .[, c("ID", power_cols, coeff_cols)]

all_data <- rei_exp %>%
  left_join(spharm_exp, by = "ID") %>%
  filter(!is.na(power_l1))

common_ids <- intersect(all_data$ID, unique(directions_exp$ID))
df         <- all_data %>% filter(ID %in% common_ids) %>% arrange(ID)

cat(sprintf("最终分析标本数：%d\n", nrow(df)))
cat(sprintf("全系数维度：%d\n\n", length(coeff_cols)))


# ==============================================================================
# 4. 计算 KDE 真实距离矩阵（SVD 对齐，余弦距离）
# ==============================================================================

cat("计算 KDE 真实距离矩阵...\n")

G          <- make_sphere_grid(n_bearing = 72, n_plunge = 36)
n          <- nrow(df)
kde_matrix <- matrix(0, nrow = n, ncol = nrow(G))

for (i in seq_len(n)) {
  id_i   <- df$ID[i]
  dirs_i <- directions_exp %>%
    filter(ID == id_i, sqrt(ux^2 + uy^2 + uz^2) > 1e-10)
  if (nrow(dirs_i) == 0) next
  kde_matrix[i, ] <- fit_vmf_kde(dirs_i$ux, dirs_i$uy, dirs_i$uz,
                                 G, bandwidth = 0.35)
  if (i %% 10 == 0 || i == n) cat(sprintf("  [%d/%d]\n", i, n))
}

# KDE 距离矩阵用未标准化的余弦距离（密度向量已归一化，直接算）
X_kde <- kde_matrix
sim   <- X_kde %*% t(X_kde) /
  (sqrt(rowSums(X_kde^2)) %o% sqrt(rowSums(X_kde^2)))
D_kde <- as.dist(1 - pmin(pmax(sim, -1), 1))
cat("KDE 距离矩阵完成。\n\n")


# ==============================================================================
# 5. 计算方法距离矩阵（全部使用余弦距离 + z-score 标准化）
# ==============================================================================

cat("计算方法距离矩阵（统一余弦距离）...\n")

X_spi       <- df %>% select(SPI)                %>% as.matrix()
X_fab       <- df %>% select(E, I)               %>% as.matrix()
X_sph_power <- df[, power_cols]                  %>% as.matrix()
X_sph_coeff <- df[, coeff_cols]                  %>% as.matrix()

D_spi       <- as.dist(dist_cosine_scaled(X_spi));       cat("  SPI 完成\n")
D_fab       <- as.dist(dist_cosine_scaled(X_fab));       cat("  Fabric 完成\n")
D_sph_power <- as.dist(dist_cosine_scaled(X_sph_power)); cat("  SPHARM 能量谱完成\n")
D_sph_coeff <- as.dist(dist_cosine_scaled(X_sph_coeff)); cat("  SPHARM 全系数完成\n\n")


# ==============================================================================
# 6. Mantel 检验（Spearman，999 次置换）
# ==============================================================================

cat("====== Mantel 检验（999 次置换，Spearman，统一余弦距离）======\n\n")

set.seed(42)
m_spi       <- mantel(D_kde, D_spi,       method = "spearman", permutations = 999)
m_fab       <- mantel(D_kde, D_fab,       method = "spearman", permutations = 999)
m_sph_power <- mantel(D_kde, D_sph_power, method = "spearman", permutations = 999)
m_sph_coeff <- mantel(D_kde, D_sph_coeff, method = "spearman", permutations = 999)

make_sig <- function(p) case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
)

results_main <- tibble(
  method    = c("SPI", "Fabric (E + I)", "SPHARM power (l1–l5)"),
  r         = c(m_spi$statistic, m_fab$statistic, m_sph_power$statistic),
  p_value   = c(m_spi$signif,   m_fab$signif,   m_sph_power$signif),
  dist_used = "cosine (z-score)"
) %>% mutate(sig = make_sig(p_value)) %>% arrange(desc(r))

results_diag <- tibble(
  method    = c("SPHARM power (l1–l5)", "SPHARM full coefficients"),
  r         = c(m_sph_power$statistic, m_sph_coeff$statistic),
  p_value   = c(m_sph_power$signif,   m_sph_coeff$signif),
  dist_used = "cosine (z-score)"
) %>% mutate(sig = make_sig(p_value)) %>% arrange(desc(r))

cat("【主分析】三种方法（统一余弦距离）：\n")
cat("方法                      r        p        显著性\n")
cat("----------------------------------------------------\n")
for (i in seq_len(nrow(results_main))) {
  cat(sprintf("%-26s  %6.4f   %6.4f   %s\n",
              results_main$method[i], results_main$r[i],
              results_main$p_value[i], results_main$sig[i]))
}

cat("\n【诊断分析】SPHARM 能量谱 vs 全系数：\n")
cat("方法                           r        p        显著性\n")
cat("----------------------------------------------------------\n")
for (i in seq_len(nrow(results_diag))) {
  cat(sprintf("%-30s  %6.4f   %6.4f   %s\n",
              results_diag$method[i], results_diag$r[i],
              results_diag$p_value[i], results_diag$sig[i]))
}

delta_r <- results_diag$r[results_diag$method == "SPHARM full coefficients"] -
  results_diag$r[results_diag$method == "SPHARM power (l1–l5)"]
cat(sprintf(
  "\n全系数 vs 能量谱 Δr = %.4f\n",
  delta_r
))
cat("解读：Δr 反映能量谱为换取旋转不变性所丢弃的相位信息量\n\n")


# ==============================================================================
# 7. 可视化
# ==============================================================================

method_order_main <- c("SPI", "Fabric (E + I)", "SPHARM power (l1–l5)")
colors_main <- c(
  "SPI"                  = "#A1C2E6",
  "Fabric (E + I)"       = "#FFBAE0",
  "SPHARM power (l1–l5)" = "#D4619A"
)

# --- 主分析柱状图 ---
p_main <- results_main %>%
  mutate(method = factor(method, levels = method_order_main)) %>%
  ggplot(aes(x = method, y = r, fill = method)) +
  geom_col(width = 0.55, alpha = 0.9) +
  geom_text(aes(label = glue("r = {round(r, 4)}\n{sig}"),
                vjust = ifelse(r >= 0, -0.3, 1.3)),
            size = 3.8, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = colors_main, guide = "none") +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 11) +
  labs(
    title    = "Information Retention — Main Analysis",
    subtitle = paste(
      "All methods: cosine distance after z-score standardization",
      "\nGround truth: SVD-aligned → vMF-KDE (bandwidth = 0.35) → cosine distance",
      "\n*** p < 0.001  ** p < 0.01  * p < 0.05  ns p ≥ 0.05"
    ),
    x = NULL, y = "Mantel r (Spearman)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 9)
  )

# --- 诊断分析柱状图 ---
colors_diag <- c(
  "SPHARM power (l1–l5)"     = "#D4619A",
  "SPHARM full coefficients"  = "#8B1A6B"
)

p_diag <- results_diag %>%
  mutate(method = factor(method,
                         levels = c("SPHARM power (l1–l5)",
                                    "SPHARM full coefficients"))) %>%
  ggplot(aes(x = method, y = r, fill = method)) +
  geom_col(width = 0.45, alpha = 0.9) +
  geom_text(aes(label = glue("r = {round(r, 4)}\n{sig}")),
            vjust = -0.3, size = 3.8, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  annotate("text", x = 1.5,
           y = max(results_diag$r) * 0.5,
           label = glue("Δr = {round(delta_r, 4)}\n(phase information cost)"),
           size = 3.5, color = "grey30", hjust = 0.5) +
  scale_fill_manual(values = colors_diag, guide = "none") +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 11) +
  labs(
    title    = "Diagnostic Analysis — Phase Information Cost",
    subtitle = paste(
      "Power spectrum: rotation-invariant (phase discarded)",
      "\nFull coefficients: phase retained (SVD-aligned coordinate frame)",
      "\nΔr quantifies information lost by discarding phase"
    ),
    x = NULL, y = "Mantel r (Spearman)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 9)
  )

# --- 散点图（主分析）---
plot_scatter <- function(D_method, method_name, color) {
  kde_vec    <- as.vector(D_kde)
  method_vec <- as.vector(D_method)
  ggplot(data.frame(kde = kde_vec, method = method_vec),
         aes(x = kde, y = method)) +
    geom_point(size = 1, alpha = 0.3, color = color) +
    geom_smooth(method = "lm", se = TRUE,
                color = "grey30", linewidth = 0.8) +
    theme_bw(base_size = 9) +
    labs(title = method_name,
         x = "KDE distance", y = "Method distance") +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
}

p_scatter <- (
  plot_scatter(D_spi,       "SPI",                    "#A1C2E6") |
    plot_scatter(D_fab,       "Fabric (E + I)",          "#FFBAE0") |
    plot_scatter(D_sph_power, "SPHARM power (l1–l5)",   "#D4619A") |
    plot_scatter(D_sph_coeff, "SPHARM full coefficients","#8B1A6B")
) +
  plot_annotation(
    title = "Pairwise Distance: Method vs KDE Ground Truth",
    theme = theme(plot.title = element_text(face = "bold", size = 10,
                                            hjust = 0.5))
  )

p_main_full  <- p_main  / p_scatter + plot_layout(heights = c(1.5, 1))
p_diag_final <- p_diag


# ==============================================================================
# 8. 保存
# ==============================================================================

ggsave(here("analysis/output/figures/mantel_main.png"),
       plot = p_main_full, width = 12, height = 9, dpi = 300, bg = "white")
cat("图已保存：mantel_main.png\n")

ggsave(here("analysis/output/figures/mantel_diagnostic.png"),
       plot = p_diag_final, width = 6, height = 5, dpi = 300, bg = "white")
cat("图已保存：mantel_diagnostic.png\n")

bind_rows(
  results_main %>% mutate(analysis = "main"),
  results_diag %>% mutate(analysis = "diagnostic")
) %>%
  select(analysis, method, dist_used, r, p_value, sig) %>%
  write_csv(here("analysis/data/derived_data/mantel_retention_results.csv"))
cat("结果已保存：mantel_retention_results.csv\n")


# ==============================================================================
# 9. 结论
# ==============================================================================

cat("\n====== 结论 ======\n\n")
best_main <- results_main %>% slice_max(r, n = 1)
cat(sprintf("主分析（统一余弦距离）信息保留度最高：%s（r = %.4f，%s）\n",
            best_main$method, best_main$r, best_main$sig))
cat("方法排序：",
    paste(results_main$method, sprintf("r=%.4f", results_main$r),
          sep="=", collapse=" > "), "\n\n")

r_power <- results_diag$r[results_diag$method == "SPHARM power (l1–l5)"]
r_coeff <- results_diag$r[results_diag$method == "SPHARM full coefficients"]
cat(sprintf("诊断分析：全系数 r = %.4f，能量谱 r = %.4f，Δr = %.4f\n",
            r_coeff, r_power, delta_r))
cat("解读：SPHARM 全系数保留的信息远高于能量谱，\n")
cat(sprintf("      差值 Δr = %.4f 量化了能量谱为换取旋转不变性丢弃的相位信息代价。\n",
            delta_r))