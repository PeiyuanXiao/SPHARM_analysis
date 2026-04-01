# ==============================================================================
# validate_rotation_all.R
# 旋转不变性验证：R / E / I + SPHARM 功率谱 + 谱熵
# 方法：Bland-Altman 分析（逐对比较三种对齐方式）
#
# 前置条件：
#   1. 已运行 align_svd.R     → directions_raw.csv + directions_aligned_svd.csv
#   2. 已运行 align_lin2024.R → directions_aligned_lin2024.csv
#   3. 已运行 python kde_to_spharm_main.py --source all
#      → validation/raw|svd|lin2024/SPHARM_direction.csv
#
# 输出：
#   控制台：所有指标的 bias / LoA 汇总表
#   analysis/output/figures/validation_ba_REI.png
#   analysis/output/figures/validation_ba_spharm.png
#   analysis/output/figures/validation_ba_entropy.png
# ==============================================================================

library(here)
library(tidyverse)
library(patchwork)
library(glue)


# ==============================================================================
# 公共函数
# ==============================================================================

# --- 方向统计量计算 ---
compute_R <- function(ux, uy, uz) {
  resultant <- sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2)
  total_len <- sum(sqrt(ux^2 + uy^2 + uz^2))
  resultant / total_len
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

# --- Bland-Altman 核心计算 ---
# 返回 bias、SD、95% LoA 及绘图所需数据框
bland_altman_calc <- function(x, y) {
  mean_xy   <- (x + y) / 2
  diff_xy   <- x - y
  bias      <- mean(diff_xy,  na.rm = TRUE)
  sd_diff   <- sd(diff_xy,    na.rm = TRUE)
  loa_upper <- bias + 1.96 * sd_diff
  loa_lower <- bias - 1.96 * sd_diff
  list(
    df        = data.frame(mean = mean_xy, diff = diff_xy),
    bias      = bias,
    sd_diff   = sd_diff,
    loa_upper = loa_upper,
    loa_lower = loa_lower
  )
}

# --- 单张 Bland-Altman 图 ---
plot_ba <- function(ba, title_str, x_label = "Mean", y_label = "Difference") {
  ggplot(ba$df, aes(x = mean, y = diff)) +
    geom_hline(yintercept = ba$bias,      color = "#D4619A",
               linewidth = 0.9, linetype = "dashed") +
    geom_hline(yintercept = ba$loa_upper, color = "#A1C2E6",
               linewidth = 0.7, linetype = "dotted") +
    geom_hline(yintercept = ba$loa_lower, color = "#A1C2E6",
               linewidth = 0.7, linetype = "dotted") +
    geom_point(size = 2, alpha = 0.75, color = "grey40") +
    annotate("text",
             x = -Inf, y = ba$bias,      hjust = -0.1, vjust = -0.4,
             label = sprintf("Bias = %.2e", ba$bias),
             color = "#D4619A", size = 2.8) +
    annotate("text",
             x = -Inf, y = ba$loa_upper, hjust = -0.1, vjust = -0.4,
             label = sprintf("+1.96 SD = %.2e", ba$loa_upper),
             color = "#A1C2E6", size = 2.5) +
    annotate("text",
             x = -Inf, y = ba$loa_lower, hjust = -0.1, vjust =  1.4,
             label = sprintf("-1.96 SD = %.2e", ba$loa_lower),
             color = "#A1C2E6", size = 2.5) +
    theme_bw(base_size = 9) +
    labs(title = title_str, x = x_label, y = y_label) +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
}

# --- 汇总表辅助 ---
summary_row <- function(metric, pair, ba) {
  tibble(
    metric    = metric,
    pair      = pair,
    bias      = ba$bias,
    sd_diff   = ba$sd_diff,
    loa_lower = ba$loa_lower,
    loa_upper = ba$loa_upper
  )
}


# ==============================================================================
# Part A：R / E / I 验证
# ==============================================================================

cat("====== Part A: R / E / I ======\n\n")

# --- 读取方向向量 CSV ---
read_directions <- function(source) {
  path <- here(glue("analysis/data/derived_data/directions_{source}.csv"))
  if (!file.exists(path))
    stop(glue("找不到：{path}"))
  read_csv(path, show_col_types = FALSE) %>% mutate(source = source)
}

dirs <- bind_rows(
  read_directions("raw"),
  read_directions("aligned_svd"),
  read_directions("aligned_lin2024")
)

# 统一 source 标签
dirs <- dirs %>%
  mutate(source = case_when(
    source == "raw"              ~ "raw",
    source == "aligned_svd"     ~ "svd",
    source == "aligned_lin2024" ~ "lin2024"
  ))

# --- 计算 R / E / I ---
rei <- dirs %>%
  group_by(ID, source) %>%
  summarise(
    SPI = compute_R(ux, uy, uz),
    Elongation = compute_EI(ux, uy, uz)$E,
    Isotropy = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

# 宽格式（每行一个标本，三列分别对应三个 source）
rei_wide <- rei %>%
  pivot_wider(names_from = source,
              values_from = c(SPI, Elongation, Isotropy),
              names_glue = "{.value}_{source}")

common_ids_rei <- rei_wide %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  pull(ID)

rei_wide <- rei_wide %>% filter(ID %in% common_ids_rei)
cat(sprintf("R/E/I 验证：%d 个标本\n\n", nrow(rei_wide)))

# --- 三对 × 三指标 = 9 张 B-A 图 ---
pairs_label <- list(
  c("raw", "svd"),
  c("raw", "lin2024"),
  c("svd", "lin2024")
)

ba_plots_rei <- list()
summary_table <- tibble()

for (metric in c("SPI", "Elongation", "Isotropy")) {
  for (pair in pairs_label) {
    src_a <- pair[1]; src_b <- pair[2]
    col_a <- glue("{metric}_{src_a}")
    col_b <- glue("{metric}_{src_b}")
    pair_label <- glue("{src_a} vs {src_b}")
    
    ba  <- bland_altman_calc(rei_wide[[col_a]], rei_wide[[col_b]])
    plt <- plot_ba(ba,
                   title_str = glue("{metric}: {pair_label}"),
                   x_label   = glue("Mean {metric}"),
                   y_label   = glue("Δ{metric} ({src_a} − {src_b})"))
    
    ba_plots_rei[[glue("{metric}_{pair_label}")]] <- plt
    summary_table <- bind_rows(summary_table,
                               summary_row(metric, pair_label, ba))
  }
}

# 拼图：3 行（指标）× 3 列（对比组）
p_rei <- wrap_plots(ba_plots_rei, ncol = 3) +
  plot_annotation(
    title    = "Bland-Altman test: SPI / Elongation / Isotropy across alignment methods",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave(here("analysis/output/figures/validation_ba_REI.png"),
       plot = p_rei, width = 12, height = 10, dpi = 300, bg = "white")
cat("图已保存：validation_ba_REI.png\n\n")


# ==============================================================================
# Part B：SPHARM 功率谱验证（各阶 Bland-Altman，分面展示）
# ==============================================================================

cat("====== Part B: SPHARM 功率谱 ======\n\n")

# --- 读取三份功率谱 CSV ---
read_spharm <- function(source) {
  path <- here(glue("analysis/data/derived_data/validation/{source}/SPHARM_direction.csv"))
  if (!file.exists(path))
    stop(glue("找不到：{path}\n请先运行：python kde_to_spharm_main.py --source all"))
  read_csv(path, show_col_types = FALSE) %>% mutate(source = source)
}

spharm_all <- bind_rows(
  read_spharm("raw"),
  read_spharm("svd"),
  read_spharm("lin2024")
)

power_cols <- spharm_all %>%
  select(starts_with("power_l")) %>%
  colnames()

# 宽格式
spharm_wide <- spharm_all %>%
  select(ID, source, all_of(power_cols), spectral_entropy) %>%
  pivot_wider(names_from  = source,
              values_from = c(all_of(power_cols), spectral_entropy),
              names_glue  = "{.value}__{source}")   # 双下划线分隔，避免列名歧义

common_ids_spharm <- spharm_wide %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  pull(ID)
spharm_wide <- spharm_wide %>% filter(ID %in% common_ids_spharm)
cat(sprintf("SPHARM 验证：%d 个标本，%d 阶\n\n",
            nrow(spharm_wide), length(power_cols)))

# --- 各阶 Bland-Altman，结果收进长表 ---
ba_power_summary <- tibble()

for (pair in pairs_label) {
  src_a <- pair[1]; src_b <- pair[2]
  pair_label <- glue("{src_a} vs {src_b}")
  
  for (pcol in power_cols) {
    col_a <- glue("{pcol}__{src_a}")
    col_b <- glue("{pcol}__{src_b}")
    ba    <- bland_altman_calc(spharm_wide[[col_a]], spharm_wide[[col_b]])
    
    degree <- as.integer(str_remove(pcol, "power_l"))
    ba_power_summary <- bind_rows(ba_power_summary, tibble(
      degree    = degree,
      pair      = pair_label,
      bias      = ba$bias,
      loa_upper = ba$loa_upper,
      loa_lower = ba$loa_lower,
      sd_diff   = ba$sd_diff
    ))
    
    summary_table <- bind_rows(summary_table,
                               summary_row(pcol, pair_label, ba))
  }
}

# 绘图：bias + LoA 随阶次变化（三对叠在同一图，分面）
p_spharm <- ggplot(ba_power_summary, aes(x = degree)) +
  geom_ribbon(aes(ymin = loa_lower, ymax = loa_upper, fill = pair),
              alpha = 0.15) +
  geom_line(aes(y = bias,      color = pair), linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    values = c("raw vs svd"      = "#D4619A",
               "raw vs lin2024"  = "#4A9A4A",
               "svd vs lin2024"  = "#3B8BD4"),
    name = "Comparison"
  ) +
  scale_fill_manual(
    values = c("raw vs svd"      = "#D4619A",
               "raw vs lin2024"  = "#4A9A4A",
               "svd vs lin2024"  = "#3B8BD4"),
    name = "Comparison"
  ) +
  scale_x_continuous(breaks = seq(0, max(ba_power_summary$degree), by = 2)) +
  theme_bw(base_size = 10) +
  labs(
    title    = "Bland-Altman Summary: SPHARM Power per Degree",
    subtitle = "Line = bias，Ribbon = 95% LoA",
    x        = "Spherical Harmonic Degree (l)",
    y        = "Difference (Method A − Method B)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

ggsave(here("analysis/output/figures/validation_ba_spharm.png"),
       plot = p_spharm, width = 11, height = 6, dpi = 300, bg = "white")
cat("图已保存：validation_ba_spharm.png\n\n")


# ==============================================================================
# Part C：谱熵 Bland-Altman
# ==============================================================================

cat("====== Part C: 谱熵 ======\n\n")

ba_plots_entropy <- list()

for (pair in pairs_label) {
  src_a <- pair[1]; src_b <- pair[2]
  pair_label <- glue("{src_a} vs {src_b}")
  col_a <- glue("spectral_entropy__{src_a}")
  col_b <- glue("spectral_entropy__{src_b}")
  
  ba  <- bland_altman_calc(spharm_wide[[col_a]], spharm_wide[[col_b]])
  plt <- plot_ba(ba,
                 title_str = pair_label,
                 x_label   = "Mean spectral entropy",
                 y_label   = glue("Δ entropy ({src_a} − {src_b})"))
  
  ba_plots_entropy[[pair_label]] <- plt
  summary_table <- bind_rows(summary_table,
                             summary_row("spectral_entropy", pair_label, ba))
}

p_entropy <- wrap_plots(ba_plots_entropy, ncol = 3) +
  plot_annotation(
    title    = "Bland-Altman: Spectral Entropy Across Alignment Methods",
    subtitle = "Pink dashed = bias，Blue dotted = 95% LoA",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave(here("analysis/output/figures/validation_ba_entropy.png"),
       plot = p_entropy, width = 12, height = 4, dpi = 300, bg = "white")
cat("图已保存：validation_ba_entropy.png\n\n")


# ==============================================================================
# 汇总表
# ==============================================================================

cat("====== 所有指标 Bland-Altman 汇总 ======\n\n")

summary_table %>%
  mutate(across(c(bias, sd_diff, loa_lower, loa_upper),
                \(x) formatC(x, format = "e", digits = 3))) %>%
  arrange(metric, pair) %>%
  print(n = Inf)

write_csv(
  summary_table %>% arrange(metric, pair),
  here("analysis/data/derived_data/validation_ba_summary.csv")
)
cat("\n汇总表已保存：validation_ba_summary.csv\n")

# ==============================================================================
# 具体数值汇总：三种方法并排输出
# ==============================================================================

cat("\n====== 具体数值汇总：R / E / I ======\n\n")

# --- R / E / I 三方法并排 ---
rei_compare <- rei_wide %>%
  select(ID,
         R_raw, R_svd, R_lin2024,
         E_raw, E_svd, E_lin2024,
         I_raw, I_svd, I_lin2024) %>%
  arrange(ID)

rei_compare %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)

write_csv(rei_compare,
          here("analysis/data/derived_data/validation_values_REI.csv"))
cat("已保存：validation_values_REI.csv\n")


cat("\n====== 具体数值汇总：谱熵 ======\n\n")

# --- 谱熵三方法并排 ---
entropy_compare <- spharm_wide %>%
  select(ID,
         entropy_raw     = `spectral_entropy__raw`,
         entropy_svd     = `spectral_entropy__svd`,
         entropy_lin2024 = `spectral_entropy__lin2024`) %>%
  arrange(ID)

entropy_compare %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)

write_csv(entropy_compare,
          here("analysis/data/derived_data/validation_values_entropy.csv"))
cat("已保存：validation_values_entropy.csv\n")


cat("\n====== 具体数值汇总：SPHARM 功率谱（l=1–5）======\n\n")

# --- 功率谱三方法并排（仅展示 l=1-5，聚焦分析用阶次）---
# 输出格式：每行一个标本，每个阶次显示三列（raw / svd / lin2024）
power_compare <- spharm_wide %>%
  select(ID, matches("^power_l[1-5]__")) %>%
  arrange(ID)

# 整理列名为可读格式：power_l1_raw, power_l1_svd, ...
power_compare <- power_compare %>%
  rename_with(
    ~ str_replace(., "__", "_"),
    matches("^power_l[1-5]__")
  )

power_compare %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)

write_csv(power_compare,
          here("analysis/data/derived_data/validation_values_power_l1_5.csv"))
cat("已保存：validation_values_power_l1_5.csv\n")


cat("\n====== 全阶功率谱数值已保存（不在控制台打印）======\n")

# --- 全阶功率谱（l=0–20）保存为 CSV，供完整参考 ---
power_all_compare <- spharm_wide %>%
  select(ID, matches("^power_l[0-9]+__")) %>%
  rename_with(~ str_replace(., "__", "_"),
              matches("^power_l[0-9]+__")) %>%
  arrange(ID)

write_csv(power_all_compare,
          here("analysis/data/derived_data/validation_values_power_all.csv"))
cat("已保存：validation_values_power_all.csv\n")

