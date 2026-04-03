# ==============================================================================
# mantel_pairwise_angle_retention.R
# 片疤方向信息保留度分析
#
# 核心思想：
#   1) 对每件标本，计算片疤方向向量两两之间的夹角分布
#   2) 将角度分布离散化为直方图，作为 rotation-invariant ground truth
#   3) 比较 SPI / Fabric / SPHARM 的特征距离矩阵与该 ground truth 的对应程度
#
# 主分析：
#   - ground truth: pairwise-angle histogram
#   - 方法距离：统一使用 z-score 后的欧氏距离
#
# 敏感性分析：
#   - SPI: 欧氏距离
#   - Fabric: 马氏距离
#   - SPHARM power / full coeff: 余弦距离
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 3）
#   - analysis/data/derived_data/directions_aligned_svd.csv
#   - analysis/data/derived_data/SPHARM_direction.csv
#
# 输出：
#   - analysis/output/figures/pairwise_angle_main.png
#   - analysis/output/figures/pairwise_angle_sensitivity.png
#   - analysis/data/derived_data/pairwise_angle_retention_results.csv
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(vegan)
library(patchwork)
library(glue)

# ==============================================================================
# 参数
# ==============================================================================

n_angle_bins <- 36
permutations <- 999
set.seed(42)

# ==============================================================================
# 读取数据
# ==============================================================================

exp_data <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"),
  sheet = 3
)

directions_aligned <- read_csv(
  here("analysis/data/derived_data/directions_aligned_svd.csv"),
  show_col_types = FALSE
)

spharm_all <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE
)

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

prepare_scaled_matrix <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 0) stop("输入矩阵没有可用列。")
  col_sd <- apply(X, 2, sd, na.rm = TRUE)
  keep   <- is.finite(col_sd) & col_sd > 1e-10
  X      <- X[, keep, drop = FALSE]
  if (ncol(X) == 0) stop("所有列都是常数列，无法标准化。")
  scale(X)
}

dist_euclidean_scaled <- function(X) {
  X <- prepare_scaled_matrix(X)
  as.dist(stats::dist(X, method = "euclidean"))
}

dist_cosine_scaled <- function(X) {
  X <- prepare_scaled_matrix(X)
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  as.dist(1 - sim)
}

dist_mahalanobis <- function(X, ridge = 1e-6) {
  X <- as.matrix(X)
  if (ncol(X) < 2) {
    # 单变量时，马氏距离退化为欧氏距离
    return(as.dist(stats::dist(scale(X), method = "euclidean")))
  }
  S <- stats::cov(X, use = "pairwise.complete.obs")
  S <- S + diag(ridge, ncol(S))
  Sinv <- solve(S)
  
  n <- nrow(X)
  D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d <- X[i, ] - X[j, ]
      D[i, j] <- sqrt(drop(t(d) %*% Sinv %*% d))
    }
  }
  as.dist(D)
}

# 计算每件标本的 pairwise-angle 直方图（概率分布）
compute_pairwise_angle_hist <- function(ux, uy, uz, breaks) {
  U <- cbind(ux, uy, uz)
  U <- U / sqrt(rowSums(U^2))
  
  if (nrow(U) < 2) return(rep(NA_real_, length(breaks) - 1))
  
  dot <- U %*% t(U)
  dot <- pmin(pmax(dot, -1), 1)
  theta <- acos(dot)
  theta_vec <- theta[upper.tri(theta)]
  
  if (length(theta_vec) == 0) return(rep(NA_real_, length(breaks) - 1))
  
  h <- hist(
    theta_vec,
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )
  
  p <- h$counts / sum(h$counts)
  as.numeric(p)
}

dist_hellinger <- function(P) {
  P <- as.matrix(P)
  rs <- rowSums(P)
  if (any(rs <= 0 | !is.finite(rs))) stop("存在无法构造概率分布的行。")
  P <- P / rs
  
  sp <- sqrt(P)
  n  <- nrow(sp)
  D  <- matrix(0, n, n)
  
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      D[i, j] <- sqrt(sum((sp[i, ] - sp[j, ])^2)) / sqrt(2)
    }
  }
  as.dist(D)
}

make_sig <- function(p) case_when(
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE      ~ "ns"
)

run_mantel_suite <- function(D_truth, method_dists, permutations = 999) {
  out <- map_dfr(names(method_dists), function(nm) {
    m <- mantel(D_truth, method_dists[[nm]],
                method = "spearman",
                permutations = permutations)
    
    tibble(
      method    = nm,
      r         = as.numeric(m$statistic),
      p_value   = as.numeric(m$signif),
      sig       = make_sig(as.numeric(m$signif))
    )
  })
  
  out %>% arrange(desc(r))
}

plot_distance_scatter <- function(D_truth, D_method, title, color) {
  dfp <- tibble(
    truth  = as.vector(D_truth),
    method = as.vector(D_method)
  )
  
  ggplot(dfp, aes(x = truth, y = method)) +
    geom_point(size = 1, alpha = 0.28, color = color) +
    geom_smooth(method = "lm", se = TRUE,
                color = "grey30", linewidth = 0.8) +
    theme_bw(base_size = 9) +
    labs(title = title, x = "Ground truth distance", y = "Method distance") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5)
    )
}

# ==============================================================================
# 1. 过滤实验石核
# ==============================================================================

exp_ids <- unique(exp_data$ID)
directions_exp <- directions_aligned %>% filter(ID %in% exp_ids)

cat(sprintf(
  "实验石核：%d 件标本，%d 条刮痕\n\n",
  n_distinct(directions_exp$ID), nrow(directions_exp)
))

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

cat(sprintf(
  "E 与 I 的相关系数：%.3f\n\n",
  cor(rei_exp$E, rei_exp$I, use = "complete.obs")
))

# ==============================================================================
# 3. 合并 SPHARM 数据
# ==============================================================================

power_cols <- paste0("power_l", 1:5)

if (!all(power_cols %in% colnames(spharm_all))) {
  stop("SPHARM 文件中缺少 power_l1~power_l5 列。")
}

# 尝试自动识别全系数列
coeff_patterns <- c("^coeff_", "^coef_", "^shc_", "^Y_", "^real_", "^imag_", "^re_", "^im_")
coeff_cols <- unique(unlist(lapply(coeff_patterns, function(p) {
  grep(p, colnames(spharm_all), value = TRUE)
})))
coeff_cols <- setdiff(coeff_cols, c("ID", power_cols))

has_coeff <- length(coeff_cols) > 0

if (has_coeff) {
  cat(sprintf("检测到全系数列：%d 个\n\n", length(coeff_cols)))
} else {
  cat("未检测到可用的全系数列；将跳过全系数诊断分析。\n\n")
}

spharm_exp <- spharm_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, all_of(power_cols), any_of(coeff_cols))

all_data <- rei_exp %>%
  left_join(spharm_exp, by = "ID") %>%
  filter(!is.na(power_l1))

common_ids <- intersect(all_data$ID, unique(directions_exp$ID))
df <- all_data %>% filter(ID %in% common_ids) %>% arrange(ID)

cat(sprintf("初始合并后标本数：%d\n", nrow(df)))

# ==============================================================================
# 4. 构造 rotation-invariant ground truth：pairwise-angle 直方图
# ==============================================================================

cat("计算 pairwise-angle ground truth ...\n")

angle_breaks <- seq(0, pi, length.out = n_angle_bins + 1)
angle_breaks[length(angle_breaks)] <- pi + 1e-8

angle_mat <- matrix(NA_real_, nrow = nrow(df), ncol = n_angle_bins)

for (i in seq_len(nrow(df))) {
  id_i <- df$ID[i]
  dirs_i <- directions_exp %>%
    filter(ID == id_i, sqrt(ux^2 + uy^2 + uz^2) > 1e-10)
  
  if (nrow(dirs_i) < 2) next
  
  angle_mat[i, ] <- compute_pairwise_angle_hist(
    dirs_i$ux, dirs_i$uy, dirs_i$uz,
    breaks = angle_breaks
  )
  
  if (i %% 10 == 0 || i == nrow(df)) {
    cat(sprintf("  [%d/%d]\n", i, nrow(df)))
  }
}

keep <- complete.cases(angle_mat)
dropped <- sum(!keep)

df <- df[keep, ] %>% arrange(ID)
angle_mat <- angle_mat[keep, , drop = FALSE]

cat(sprintf("最终进入分析的标本数：%d\n", nrow(df)))
cat(sprintf("因 angle 分布无法构造而剔除的标本数：%d\n\n", dropped))

D_truth <- dist_hellinger(angle_mat)
cat("Rotation-invariant ground truth 距离矩阵完成。\n\n")

# ==============================================================================
# 5. 方法距离矩阵
#    主分析：统一使用 z-score 后的欧氏距离
#    敏感性分析：SPI=欧氏，Fabric=马氏，SPHARM=余弦
# ==============================================================================

X_spi       <- df %>% select(SPI)          %>% as.matrix()
X_fab       <- df %>% select(E, I)         %>% as.matrix()
X_sph_power <- df %>% select(all_of(power_cols)) %>% as.matrix()

D_spi_main       <- dist_euclidean_scaled(X_spi)
D_fab_main       <- dist_euclidean_scaled(X_fab)
D_sph_power_main <- dist_euclidean_scaled(X_sph_power)

method_dists_main <- list(
  "SPI"                  = D_spi_main,
  "Fabric (E + I)"       = D_fab_main,
  "SPHARM power (l1–l5)" = D_sph_power_main
)

if (has_coeff) {
  X_sph_coeff <- df %>% select(all_of(coeff_cols)) %>% as.matrix()
  D_sph_coeff_main <- dist_euclidean_scaled(X_sph_coeff)
  method_dists_main[["SPHARM full coefficients"]] <- D_sph_coeff_main
}

# 敏感性分析：每种方法用其更“自然”的度量
D_spi_nat       <- dist_euclidean_scaled(X_spi)  # 一维下就是标准化绝对差
D_fab_nat       <- dist_mahalanobis(X_fab)
D_sph_power_nat <- dist_cosine_scaled(X_sph_power)

method_dists_nat <- list(
  "SPI"                  = D_spi_nat,
  "Fabric (Mahalanobis)" = D_fab_nat,
  "SPHARM power (cosine)"= D_sph_power_nat
)

if (has_coeff) {
  X_sph_coeff <- df %>% select(all_of(coeff_cols)) %>% as.matrix()
  D_sph_coeff_nat <- dist_cosine_scaled(X_sph_coeff)
  method_dists_nat[["SPHARM full coefficients"]] <- D_sph_coeff_nat
}

# ==============================================================================
# 6. Mantel 检验
# ==============================================================================

cat("====== Mantel 检验：主分析（统一 z-score 欧氏距离）======\n\n")
results_main <- run_mantel_suite(D_truth, method_dists_main, permutations = permutations)
results_main <- results_main %>%
  mutate(
    analysis  = "main",
    dist_used = "ground truth: Hellinger(pairwise-angle hist); methods: euclidean(z-score)"
  )

cat("方法                           r        p        显著性\n")
cat("---------------------------------------------------------------\n")
for (i in seq_len(nrow(results_main))) {
  cat(sprintf(
    "%-30s  %6.4f   %6.4f   %s\n",
    results_main$method[i], results_main$r[i],
    results_main$p_value[i], results_main$sig[i]
  ))
}

cat("\n====== Mantel 检验：敏感性分析（自然距离）======\n\n")
results_nat <- run_mantel_suite(D_truth, method_dists_nat, permutations = permutations)
results_nat <- results_nat %>%
  mutate(
    analysis  = "sensitivity",
    dist_used = "ground truth: Hellinger(pairwise-angle hist); methods: mixed natural metrics"
  )

cat("方法                           r        p        显著性\n")
cat("---------------------------------------------------------------\n")
for (i in seq_len(nrow(results_nat))) {
  cat(sprintf(
    "%-30s  %6.4f   %6.4f   %s\n",
    results_nat$method[i], results_nat$r[i],
    results_nat$p_value[i], results_nat$sig[i]
  ))
}

# 如果有 full coefficients，给一个诊断性 Δr
if (has_coeff && "SPHARM full coefficients" %in% results_main$method) {
  r_power_main <- results_main$r[results_main$method == "SPHARM power (l1–l5)"]
  r_coeff_main <- results_main$r[results_main$method == "SPHARM full coefficients"]
  delta_main <- r_coeff_main - r_power_main
  
  r_power_nat <- results_nat$r[results_nat$method == "SPHARM power (cosine)"]
  r_coeff_nat <- results_nat$r[results_nat$method == "SPHARM full coefficients"]
  delta_nat <- r_coeff_nat - r_power_nat
  
  cat(sprintf("\n主分析：全系数 vs 能量谱 Δr = %.4f\n", delta_main))
  cat(sprintf("敏感性分析：全系数 vs 能量谱 Δr = %.4f\n\n", delta_nat))
}

# ==============================================================================
# 7. 可视化
# ==============================================================================

method_order_main <- names(method_dists_main)
palette_main <- c(
  "SPI"                  = "#A1C2E6",
  "Fabric (E + I)"       = "#FFBAE0",
  "SPHARM power (l1–l5)" = "#D4619A",
  "SPHARM full coefficients" = "#8B1A6B"
)

p_main <- results_main %>%
  mutate(method = factor(method, levels = method_order_main)) %>%
  ggplot(aes(x = method, y = r, fill = method)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(
    aes(label = glue("r = {round(r, 4)}\n{sig}")),
    vjust = -0.35, size = 3.6, fontface = "bold"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = palette_main, guide = "none") +
  scale_y_continuous(limits = c(0, max(results_main$r, na.rm = TRUE) + 0.15)) +
  theme_bw(base_size = 11) +
  labs(
    title = "Information Retention — Pairwise-Angle Ground Truth",
    subtitle = paste(
      "Ground truth = Hellinger distance between pairwise-angle histograms",
      "\nMethods = z-score standardized Euclidean distance",
      "\n*** p < 0.001  ** p < 0.01  * p < 0.05  ns p ≥ 0.05"
    ),
    x = NULL,
    y = "Mantel r (Spearman)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 9)
  )

p_nat <- results_nat %>%
  mutate(method = factor(method, levels = names(method_dists_nat))) %>%
  ggplot(aes(x = method, y = r, fill = method)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(
    aes(label = glue("r = {round(r, 4)}\n{sig}")),
    vjust = -0.35, size = 3.6, fontface = "bold"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = palette_main, guide = "none") +
  scale_y_continuous(limits = c(0, max(results_nat$r, na.rm = TRUE) + 0.15)) +
  theme_bw(base_size = 11) +
  labs(
    title = "Sensitivity Analysis — Method-Specific Distances",
    subtitle = paste(
      "SPI: Euclidean   Fabric: Mahalanobis   SPHARM: Cosine",
      "\nGround truth still = Hellinger(pairwise-angle hist)"
    ),
    x = NULL,
    y = "Mantel r (Spearman)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 9)
  )

# 散点图：主分析
p_scatter_main <- (
  plot_distance_scatter(D_truth, D_spi_main,       "SPI",                  "#A1C2E6") |
    plot_distance_scatter(D_truth, D_fab_main,       "Fabric (E + I)",       "#FFBAE0") |
    plot_distance_scatter(D_truth, D_sph_power_main, "SPHARM power (l1–l5)", "#D4619A")
)

if (has_coeff) {
  p_scatter_main <- p_scatter_main |
    plot_distance_scatter(D_truth, D_sph_coeff_main, "SPHARM full coefficients", "#8B1A6B")
}

p_scatter_main <- p_scatter_main +
  plot_annotation(
    title = "Pairwise Distances: Methods vs Rotation-Invariant Ground Truth",
    theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  )

p_main_full <- p_main / p_scatter_main + plot_layout(heights = c(1.2, 1))
p_nat_full  <- p_nat

# ==============================================================================
# 8. 保存
# ==============================================================================

ggsave(
  here("analysis/output/figures/pairwise_angle_main.png"),
  plot = p_main_full, width = 12, height = 9, dpi = 300, bg = "white"
)
cat("图已保存：pairwise_angle_main.png\n")

ggsave(
  here("analysis/output/figures/pairwise_angle_sensitivity.png"),
  plot = p_nat_full, width = 7, height = 5, dpi = 300, bg = "white"
)
cat("图已保存：pairwise_angle_sensitivity.png\n")

bind_rows(results_main, results_nat) %>%
  select(analysis, method, dist_used, r, p_value, sig) %>%
  write_csv(here("analysis/data/derived_data/pairwise_angle_retention_results.csv"))

cat("结果已保存：pairwise_angle_retention_results.csv\n")

# ==============================================================================
# 9. 结论
# ==============================================================================

cat("\n====== 结论 ======\n\n")

best_main <- results_main %>% slice_max(r, n = 1)
cat(sprintf(
  "主分析（rotation-invariant ground truth）信息保留度最高：%s（r = %.4f，%s）\n",
  best_main$method, best_main$r, best_main$sig
))
cat("主分析排序：",
    paste(results_main$method, sprintf("r=%.4f", results_main$r),
          sep = "=", collapse = " > "), "\n\n")

best_nat <- results_nat %>% slice_max(r, n = 1)
cat(sprintf(
  "敏感性分析信息保留度最高：%s（r = %.4f，%s）\n",
  best_nat$method, best_nat$r, best_nat$sig
))
cat("敏感性分析排序：",
    paste(results_nat$method, sprintf("r=%.4f", results_nat$r),
          sep = "=", collapse = " > "), "\n\n")

cat("解读：pairwise-angle 直方图提供了一个完全旋转不变、且不依赖任何摘要方法的参考空间。\n")
cat("主分析回答的是：谁最能保留原始方向分布的内部角结构。\n")
cat("敏感性分析回答的是：当使用各方法更自然的距离度量时，结论是否稳定。\n")






df_analysis <- df %>%
  select(ID, SPI, E, I, starts_with("power_l")) %>%
  drop_na()

# SPHARM PCA
X_sph <- df_analysis %>%
  select(starts_with("power_l")) %>%
  scale()

pca <- prcomp(X_sph, center = FALSE, scale. = FALSE)

# 取前2个主成分（解释大部分方差）
df_analysis <- df_analysis %>%
  mutate(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2]
  )

m_spi_1 <- lm(PC1 ~ SPI, data = df_analysis)
m_spi_2 <- lm(PC2 ~ SPI, data = df_analysis)

m_fab_1 <- lm(PC1 ~ E + I, data = df_analysis)
m_fab_2 <- lm(PC2 ~ E + I, data = df_analysis)

m_all_1 <- lm(PC1 ~ SPI + E + I, data = df_analysis)
m_all_2 <- lm(PC2 ~ SPI + E + I, data = df_analysis)

summary(m_spi_1)$r.squared
summary(m_fab_1)$r.squared
summary(m_all_1)$r.squared


df_analysis <- df_analysis %>%
  mutate(
    res_PC1 = resid(m_all_1),
    res_PC2 = resid(m_all_2)
  )

ggplot(df_analysis, aes(x = res_PC1, y = res_PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Residual structure (SPHARM not explained by SPI + Fabric)",
    x = "Residual PC1",
    y = "Residual PC2"
  )



