# ==============================================================================
# power_order_selection.R
# 功率谱阶数筛选脚本——形态谱 & 方向谱（全 20 阶）
#
# 分析内容：
#   1. 各阶描述性指标汇总（均值、标准差、方差、变异系数、衰减率、SNR）
#   2. 四联可视化图：方差图、变异系数图、衰减率图、累积能量图
#
# 输入：
#   analysis/data/derived_data/SPHARM_direction.csv
#   analysis/data/derived_data/SPHARM_morphology.csv
#
# 输出：
#   analysis/output/figures/OrderSelection_Diagnostics.png
#   analysis/data/derived_data/OrderSelection_stats_direction.csv
#   analysis/data/derived_data/OrderSelection_stats_morphology.csv
# ==============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)


# ==============================================================================
# ---- 参数设置 ----
# ==============================================================================

POWER_COLS_ALL <- paste0("power_l", 1:20)   # 全 20 阶
EXP_PREFIX     <- "EXP"
CANDIDATE_CUTOFFS <- c(4, 5, 6)             # 候选截断点（图上标注竖线）


# ==============================================================================
# ---- 读取数据 ----
# ==============================================================================

SPHARM_direction  <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE
)
SPHARM_morphology <- read_csv(
  here("analysis/data/derived_data/SPHARM_morphology.csv"),
  show_col_types = FALSE
)

dir_exp   <- SPHARM_direction  %>% filter(str_starts(ID, EXP_PREFIX))
morph_exp <- SPHARM_morphology %>% filter(str_starts(ID, EXP_PREFIX))

cat(sprintf("EXP 标本数：方向谱 n=%d，形态谱 n=%d\n",
            nrow(dir_exp), nrow(morph_exp)))

check_cols <- function(df, label) {
  missing   <- setdiff(POWER_COLS_ALL, colnames(df))
  available <- intersect(POWER_COLS_ALL, colnames(df))
  if (length(missing) > 0) {
    warning(sprintf("%s 缺少列：%s", label, paste(missing, collapse = ", ")))
    cat(sprintf("  %s：实际可用 %d 阶（%s ~ %s）\n",
                label, length(available),
                available[1], available[length(available)]))
  } else {
    cat(sprintf("  %s：全部 20 阶均可用\n", label))
  }
  available
}

dir_cols   <- check_cols(dir_exp,   "方向谱")
morph_cols <- check_cols(morph_exp, "形态谱")


# ==============================================================================
# ---- 计算各阶描述性指标 ----
# ==============================================================================

compute_order_stats <- function(df, cols, label) {
  
  mat      <- df %>% select(all_of(cols)) %>% as.matrix()
  n_orders <- length(cols)
  
  col_means  <- colMeans(mat, na.rm = TRUE)
  col_sds    <- apply(mat, 2, sd,  na.rm = TRUE)
  col_vars   <- apply(mat, 2, var, na.rm = TRUE)
  col_cvs    <- col_sds / col_means * 100
  col_snr    <- col_means / col_sds
  row_sums   <- rowSums(mat, na.rm = TRUE)
  total_mean <- mean(row_sums)
  cumul_pct  <- cumsum(col_means) / total_mean * 100
  pct_each   <- col_means / total_mean * 100
  decay_rate <- c(NA, col_means[-1] / col_means[-n_orders])
  
  stats_df <- tibble(
    source      = label,
    order       = seq_len(n_orders),
    order_label = paste0("l=", seq_len(n_orders)),
    mean        = round(col_means,  6),
    sd          = round(col_sds,    6),
    variance    = round(col_vars,   6),
    cv_pct      = round(col_cvs,    2),
    snr         = round(col_snr,    4),
    pct_energy  = round(pct_each,   3),
    cumul_pct   = round(cumul_pct,  3),
    decay_rate  = round(decay_rate, 4)
  )
  
  cat(sprintf("\n======== %s 各阶描述性统计 ========\n", label))
  print(stats_df %>%
          select(order, mean, sd, variance, cv_pct, snr,
                 pct_energy, cumul_pct, decay_rate) %>%
          as.data.frame())
  
  cat(sprintf("\n  行和：min=%.4f, max=%.4f, mean=%.4f, sd=%.4f\n",
              min(row_sums), max(row_sums), mean(row_sums), sd(row_sums)))
  
  snr_drop <- which(diff(col_snr) < -0.3)
  if (length(snr_drop) > 0)
    cat(sprintf("  SNR 明显下降（前一阶）：l=%s\n",
                paste(snr_drop, collapse = ", ")))
  
  for (thr in c(90, 95, 99, 99.5)) {
    k <- which(cumul_pct >= thr)[1]
    cat(sprintf("  累积能量 >= %5.1f%%：前 %d 阶\n", thr, k))
  }
  
  stats_df
}

stats_dir   <- compute_order_stats(dir_exp,   dir_cols,   "方向谱")
stats_morph <- compute_order_stats(morph_exp, morph_cols, "形态谱")

write_csv(stats_dir,
          here("analysis/data/derived_data/OrderSelection_stats_direction.csv"))
write_csv(stats_morph,
          here("analysis/data/derived_data/OrderSelection_stats_morphology.csv"))
cat("\n已保存：OrderSelection_stats_direction.csv\n")
cat("已保存：OrderSelection_stats_morphology.csv\n")


# ==============================================================================
# ---- 可视化 ----
# ==============================================================================

stats_all <- bind_rows(stats_dir, stats_morph) %>%
  mutate(source = factor(source, levels = c("形态谱", "方向谱")))

pal_source <- c("形态谱" = "#7EB8C9", "方向谱" = "#E6B89C")
cutoff_df  <- tibble(order = CANDIDATE_CUTOFFS)

theme_diag <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.subtitle    = element_text(hjust = 0.5, size = 8, color = "grey50"),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  )

x_scale <- scale_x_continuous(
  breaks = 1:max(stats_all$order),
  labels = paste0("l=", 1:max(stats_all$order)),
  expand = expansion(add = 0.4)
)

cutoff_lines <- geom_vline(
  data = cutoff_df, aes(xintercept = order),
  linetype = "dashed", color = "grey40", linewidth = 0.45, alpha = 0.75
)

cutoff_labels <- geom_text(
  data = cutoff_df,
  aes(x = order, y = Inf, label = paste0("l=", order)),
  vjust = 1.4, hjust = -0.15, size = 2.5, color = "grey35"
)

# 图1：方差
p_var <- ggplot(stats_all, aes(x = order, y = variance, color = source)) +
  cutoff_lines +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 2.2, alpha = 0.9) +
  cutoff_labels +
  scale_color_manual(values = pal_source) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  theme_diag +
  labs(title    = "Variance per order",
       subtitle = "Absolute between-specimen variability",
       x = NULL, y = "Variance")

# 图2：变异系数
p_cv <- ggplot(stats_all, aes(x = order, y = cv_pct, color = source)) +
  cutoff_lines +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 2.2, alpha = 0.9) +
  cutoff_labels +
  geom_hline(yintercept = 100, linetype = "dotted",
             color = "#C0392B", linewidth = 0.5, alpha = 0.7) +
  annotate("text", x = max(stats_all$order) * 0.82, y = 103,
           label = "CV = 100% (SNR = 1)", size = 2.3,
           color = "#C0392B", hjust = 0) +
  scale_color_manual(values = pal_source) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_diag +
  labs(title    = "Coefficient of variation per order (%)",
       subtitle = "CV = SD / mean x 100; above 100% = noise-dominated",
       x = NULL, y = "CV (%)")

# 图3：衰减率
p_decay <- ggplot(stats_all %>% filter(!is.na(decay_rate)),
                  aes(x = order, y = decay_rate, color = source)) +
  cutoff_lines +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 2.2, alpha = 0.9) +
  cutoff_labels +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "grey55", linewidth = 0.45) +
  annotate("text", x = max(stats_all$order) * 0.82, y = 0.52,
           label = "ratio = 0.5", size = 2.3, color = "grey45", hjust = 0) +
  scale_color_manual(values = pal_source) +
  x_scale +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.02, 0.12))) +
  theme_diag +
  labs(title    = "Energy decay rate per order",
       subtitle = "Ratio = mean(l) / mean(l-1); elbow = signal-noise boundary",
       x = NULL, y = "Decay rate")

# 图4：累积能量
ref_lines <- tibble(yval = c(90, 95, 99), label = c("90%", "95%", "99%"))

p_cumul <- ggplot(stats_all, aes(x = order, y = cumul_pct, color = source)) +
  geom_hline(data = ref_lines, aes(yintercept = yval),
             linetype = "dotted", color = "grey55", linewidth = 0.45) +
  geom_text(data = ref_lines,
            aes(x = max(stats_all$order) * 0.9, y = yval + 0.8, label = label),
            size = 2.3, color = "grey45", hjust = 0, inherit.aes = FALSE) +
  cutoff_lines +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 2.2, alpha = 0.9) +
  cutoff_labels +
  scale_color_manual(values = pal_source, name = NULL) +
  x_scale +
  scale_y_continuous(limits = c(0, 102), breaks = seq(0, 100, by = 20),
                     expand = expansion(mult = c(0.01, 0.04))) +
  theme_diag +
  theme(legend.position = "right", legend.text = element_text(size = 9)) +
  labs(title    = "Cumulative energy (%)",
       subtitle = "Proportion of total power spectrum energy retained",
       x = NULL, y = "Cumulative energy (%)")

# 组合
p_combined <- (p_var | p_cv) / (p_decay | p_cumul) +
  plot_annotation(
    title    = "Power Spectrum Order Selection Diagnostics",
    subtitle = sprintf(
      "EXP specimens: morphology n=%d, direction n=%d | Dashed lines = candidate cutoffs (l=%s)",
      nrow(morph_exp), nrow(dir_exp),
      paste(CANDIDATE_CUTOFFS, collapse = ", ")
    ),
    caption = paste0(
      "Variance: absolute between-specimen variability. ",
      "CV: SD/mean x 100%; CV > 100% = noise-dominated. ",
      "Decay rate: mean(l)/mean(l-1); rapid drop marks signal-noise boundary. ",
      "Cumulative energy: proportion of total mean power retained."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
      plot.caption  = element_text(size = 7, color = "grey55",
                                   hjust = 0, lineheight = 1.3)
    )
  )

ggsave(here("analysis/output/figures/OrderSelection_Diagnostics.png"),
       plot = p_combined, width = 14, height = 10, dpi = 300, bg = "white")
cat("\n图已保存：OrderSelection_Diagnostics.png\n")


# ==============================================================================
# ---- 辅助：各阶指标对照打印 ----
# ==============================================================================

cat("\n======== 各阶指标对照 ========\n")
bind_rows(stats_dir, stats_morph) %>%
  select(source, order, mean, sd, snr, cv_pct, cumul_pct, decay_rate) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = 50)

cat("\n========== 阶数筛选脚本执行完成 ==========\n")
cat("输出文件：\n")
cat("  OrderSelection_Diagnostics.png\n")
cat("  OrderSelection_stats_direction.csv\n")
cat("  OrderSelection_stats_morphology.csv\n")