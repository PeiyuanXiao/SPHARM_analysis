# ==============================================================================
# mantel_method_redundancy.R
# 方法间冗余性分析：SPI / Fabric / SPHARM 的 partial Mantel 检验
#
# 核心问题：
#   三个方法的距离矩阵之间，在控制第三方后是否仍有显著关联？
#   即各方法编码了多少独有的结构信息？
#
# 分析矩阵：
#   D_spi   : SPI（z-score 欧氏距离）
#   D_fab   : Fabric E + I（z-score 欧氏距离）
#   D_sph   : SPHARM power l1–l5（z-score 欧氏距离）
#
# 输出：
#   analysis/output/figures/method_redundancy_mantel.png
#   analysis/data/derived_data/method_redundancy_results.csv
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

permutations <- 999
set.seed(42)

# ==============================================================================
# 公共函数（与主脚本保持一致）
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
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE      ~ "ns"
)

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
# 过滤实验石核，计算 SPI 和 Fabric
# ==============================================================================

exp_ids         <- unique(exp_data$ID)
directions_exp  <- directions_aligned %>% filter(ID %in% exp_ids)

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
  filter(
    !is.na(SPI),
    !is.na(E), !is.na(I),
    !is.na(power_l1)
  ) %>%
  arrange(ID)

cat(sprintf("进入分析的标本数：%d\n\n", nrow(df)))

# ==============================================================================
# 构造三个方法的距离矩阵
# ==============================================================================

X_spi <- df %>% select(SPI)               %>% as.matrix()
X_fab <- df %>% select(E, I)              %>% as.matrix()
X_sph <- df %>% select(all_of(power_cols)) %>% as.matrix()

D_spi <- dist_euclidean_scaled(X_spi)
D_fab <- dist_euclidean_scaled(X_fab)
D_sph <- dist_euclidean_scaled(X_sph)

# ==============================================================================
# 简单 Mantel（两两之间，无控制变量）
# ==============================================================================

cat("======================================================\n")
cat(" 简单 Mantel 检验（两两之间）\n")
cat("======================================================\n\n")

run_simple <- function(Da, Db, label_a, label_b) {
  m <- mantel(Da, Db, method = "spearman", permutations = permutations)
  cat(sprintf(
    "%-30s vs %-30s  r = %6.4f  p = %6.4f  %s\n",
    label_a, label_b,
    as.numeric(m$statistic), as.numeric(m$signif),
    make_sig(as.numeric(m$signif))
  ))
  tibble(
    type    = "simple",
    A       = label_a,
    B       = label_b,
    control = NA_character_,
    r       = as.numeric(m$statistic),
    p_value = as.numeric(m$signif),
    sig     = make_sig(as.numeric(m$signif))
  )
}

r_simple <- bind_rows(
  run_simple(D_spi, D_fab, "SPI",    "Fabric"),
  run_simple(D_spi, D_sph, "SPI",    "SPHARM power"),
  run_simple(D_fab, D_sph, "Fabric", "SPHARM power")
)

# ==============================================================================
# Partial Mantel（控制第三个矩阵后的两两关系）
# ==============================================================================

cat("\n======================================================\n")
cat(" Partial Mantel 检验（控制第三个方法后）\n")
cat("======================================================\n\n")

run_partial <- function(Da, Db, Dc, label_a, label_b, label_c) {
  m <- mantel.partial(
    Da, Db, Dc,
    method = "spearman",
    permutations = permutations
  )
  cat(sprintf(
    "%-30s vs %-30s  | %-18s  r = %6.4f  p = %6.4f  %s\n",
    label_a, label_b, label_c,
    as.numeric(m$statistic), as.numeric(m$signif),
    make_sig(as.numeric(m$signif))
  ))
  tibble(
    type    = "partial",
    A       = label_a,
    B       = label_b,
    control = label_c,
    r       = as.numeric(m$statistic),
    p_value = as.numeric(m$signif),
    sig     = make_sig(as.numeric(m$signif))
  )
}

r_partial <- bind_rows(
  # SPI vs Fabric，控制 SPHARM
  run_partial(D_spi, D_fab, D_sph,
              "SPI", "Fabric", "SPHARM power"),
  
  # SPI vs SPHARM，控制 Fabric
  run_partial(D_spi, D_sph, D_fab,
              "SPI", "SPHARM power", "Fabric"),
  
  # Fabric vs SPHARM，控制 SPI
  run_partial(D_fab, D_sph, D_spi,
              "Fabric", "SPHARM power", "SPI")
)

# ==============================================================================
# 汇总结果并打印解读
# ==============================================================================

results_all <- bind_rows(r_simple, r_partial)

cat("\n======================================================\n")
cat(" 解读\n")
cat("======================================================\n\n")

for (i in seq_len(nrow(r_partial))) {
  row <- r_partial[i, ]
  direction <- if (row$r > 0) "正相关" else "负相关"
  conclusion <- if (row$p_value < 0.05) {
    glue("控制 {row$control} 后，{row$A} 与 {row$B} 仍有显著{direction}（r={round(row$r,4)}, {row$sig}），说明两者不完全冗余。")
  } else {
    glue("控制 {row$control} 后，{row$A} 与 {row$B} 关联消失（r={round(row$r,4)}, {row$sig}），说明两者的共同结构可能主要由 {row$control} 解释。")
  }
  cat(conclusion, "\n\n")
}

# ==============================================================================
# 可视化
# ==============================================================================

# --- 图1：简单 Mantel 气泡图（r 值矩阵形式）---

bubble_df <- r_simple %>%
  mutate(
    label = glue("r = {round(r, 4)}\n{sig}")
  ) %>%
  bind_rows(
    mutate(., A_tmp = B, B = A, A = A_tmp) %>% select(-A_tmp)
  )

# 构造对角线（自身相关 = 1）
methods_order <- c("SPI", "Fabric", "SPHARM power")
diag_df <- tibble(
  A = methods_order,
  B = methods_order,
  r = 1,
  label = "r = 1.000"
)

bubble_full <- bind_rows(bubble_df, diag_df) %>%
  mutate(
    A = factor(A, levels = methods_order),
    B = factor(B, levels = rev(methods_order))
  )

p_bubble <- ggplot(bubble_full, aes(x = A, y = B)) +
  geom_point(aes(size = abs(r), fill = r),
             shape = 21, color = "grey30", stroke = 0.5) +
  geom_text(aes(label = label),
            size = 3.2, vjust = 2.2, color = "grey20") +
  scale_size_continuous(range = c(4, 18), guide = "none") +
  scale_fill_gradient2(
    low  = "#3B82F6", mid = "white", high = "#EF4444",
    midpoint = 0, limits = c(-1, 1),
    name = "Mantel r"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Simple Mantel: Pairwise Method Correlations",
    subtitle = "Spearman correlation between z-score Euclidean distance matrices"
  )

# --- 图2：Partial Mantel 前后 r 值对比（dumbbell 图）---

compare_df <- r_simple %>%
  mutate(pair = glue("{A} vs {B}")) %>%
  select(pair, r_simple = r, sig_simple = sig) %>%
  left_join(
    r_partial %>%
      mutate(pair = glue("{A} vs {B}")) %>%
      select(pair, r_partial = r, sig_partial = sig, control),
    by = "pair"
  ) %>%
  filter(!is.na(r_partial)) %>%
  mutate(
    pair_label = glue("{pair}\n(控制 {control})"),
    delta_r    = r_partial - r_simple
  )

p_dumbbell <- compare_df %>%
  ggplot(aes(y = reorder(pair_label, r_simple))) +
  geom_segment(
    aes(x = r_simple, xend = r_partial,
        yend = reorder(pair_label, r_simple)),
    color = "grey60", linewidth = 1.2,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  geom_point(aes(x = r_simple),  color = "#3B82F6", size = 4) +
  geom_point(aes(x = r_partial), color = "#EF4444", size = 4) +
  geom_text(
    aes(x = r_simple,
        label = glue("r={round(r_simple,4)} {sig_simple}")),
    hjust = 1.15, size = 3, color = "#3B82F6"
  ) +
  geom_text(
    aes(x = r_partial,
        label = glue("r={round(r_partial,4)} {sig_partial}")),
    hjust = -0.15, size = 3, color = "#EF4444"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  scale_x_continuous(
    limits = function(x) c(min(x) - 0.15, max(x) + 0.15),
    name   = "Mantel r (Spearman)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.title.y  = element_blank(),
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = "Partial Mantel: Change in r After Controlling for Third Method",
    subtitle = "Blue dot = simple Mantel r  |  Red dot = partial Mantel r  |  Arrow shows direction of change"
  )

# --- 组合输出 ---

p_combined <- p_bubble / p_dumbbell +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title    = "Method Redundancy Analysis: SPI / Fabric / SPHARM power",
    subtitle = glue(
      "Spearman Mantel, {permutations} permutations  |  ",
      "*** p<0.001  ** p<0.01  * p<0.05  ns p≥0.05"
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
  here("analysis/output/figures/method_redundancy_mantel.png"),
  plot  = p_combined,
  width = 10, height = 11,
  dpi   = 300, bg = "white"
)
cat("\n图已保存：method_redundancy_mantel.png\n")

results_all %>%
  arrange(type, A, B) %>%
  write_csv(here("analysis/data/derived_data/method_redundancy_results.csv"))
cat("结果已保存：method_redundancy_results.csv\n")