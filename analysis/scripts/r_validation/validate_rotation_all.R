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
#   analysis/output/figures/validation_ba_svd_rotated_spharm.png
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
    geom_hline(yintercept = ba$bias,      color = "#802520",
               linewidth = 0.9, linetype = "dashed") +
    geom_hline(yintercept = ba$loa_upper, color = "#5C7F71",
               linewidth = 0.7, linetype = "dotted") +
    geom_hline(yintercept = ba$loa_lower, color = "#5C7F71",
               linewidth = 0.7, linetype = "dotted") +
    geom_point(size = 2.5, alpha = 0.75, color = "#B8B8B8") +
    scale_y_continuous(labels = \(x) formatC(x, format = "e", digits = 1)) +
    annotate("text",
             x = -Inf, y = ba$bias,      hjust = -0.1, vjust = -0.4,
             label = sprintf("Bias = %.2e", ba$bias),
             color = "#802520", size = 2.8) +
    annotate("text",
             x = -Inf, y = ba$loa_upper, hjust = -0.1, vjust = -0.4,
             label = sprintf("+1.96 SD = %.2e", ba$loa_upper),
             color = "#5C7F71", size = 2.5) +
    annotate("text",
             x = -Inf, y = ba$loa_lower, hjust = -0.1, vjust =  1.4,
             label = sprintf("-1.96 SD = %.2e", ba$loa_lower),
             color = "#5C7F71", size = 2.5) +
    theme_bw(base_size = 9) +
    labs(title = title_str, x = x_label, y = y_label) +
    theme(panel.grid = element_blank(), 
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
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
) %>%
  filter(!str_starts(ID, "IM_"))

# 统一 source 标签
dirs <- dirs %>%
  mutate(source = case_when(
    source == "raw"              ~ "raw",
    source == "aligned_svd"     ~ "svd",
    source == "aligned_lin2024" ~ "lin2024"
  ))

# --- 计算 R / E / I（修复：避免重复特征值分解）---
rei <- dirs %>%
  group_by(ID, source) %>%
  summarise(
    R  = compute_R(ux, uy, uz),
    ei = list(compute_EI(ux, uy, uz)),
    .groups = "drop"
  ) %>%
  mutate(
    E = map_dbl(ei, "E"),
    I = map_dbl(ei, "I")
  ) %>%
  select(-ei)

# 宽格式（每行一个标本，三列分别对应三个 source）
rei_wide <- rei %>%
  pivot_wider(names_from = source,
              values_from = c(R, E, I),
              names_glue = "{.value}_{source}")

common_ids_rei <- rei_wide %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  pull(ID)

rei_wide <- rei_wide %>% filter(ID %in% common_ids_rei)
cat(sprintf("R/E/I 验证：%d 个标本\n\n", nrow(rei_wide)))

# --- 对比组定义（语义命名）---
# pair 列表：每项为 c(source_a, source_b, 显示标签)
pairs_label <- list(
  c("raw",    "svd",     "none-aligned vs techno-aligned"),
  c("raw",    "lin2024", "none-aligned vs morph-aligned"),
  c("svd",    "lin2024", "techno-aligned vs morph-aligned")
)

ba_plots_rei <- list()
summary_table <- tibble()

metric_labels <- c(
  R = "SPI",
  E = "Elongation",
  I = "Isotropy"
)

for (metric in c("R", "E", "I")) {
  for (pair in pairs_label) {
    src_a      <- pair[1]
    src_b      <- pair[2]
    pair_label <- pair[3]   # 使用语义名称
    col_a      <- glue("{metric}_{src_a}")
    col_b      <- glue("{metric}_{src_b}")
    
    metric_name <- metric_labels[[metric]]
    
    ba  <- bland_altman_calc(rei_wide[[col_a]], rei_wide[[col_b]])
    plt <- plot_ba(ba,
                   title_str = glue("{metric_name}: {pair_label}"),
                   x_label   = glue("Mean of measures"),
                   y_label   = glue("Difference of measures"))
    
    ba_plots_rei[[glue("{metric}_{pair_label}")]] <- plt
    summary_table <- bind_rows(summary_table,
                               summary_row(metric, pair_label, ba))
  }
}

# 拼图：3 行（指标）× 3 列（对比组）
p_rei <- wrap_plots(ba_plots_rei, ncol = 3)
p_rei

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
) %>%
  filter(!str_starts(ID, "IM_"))

# 单独保留 svd 数据框供 Part D 使用
df_svd <- spharm_all %>% filter(source == "svd")

power_cols <- spharm_all %>%
  select(starts_with("power_l")) %>%
  colnames()

# 宽格式
spharm_wide <- spharm_all %>%
  select(ID, source, all_of(power_cols), spectral_entropy) %>%
  pivot_wider(names_from  = source,
              values_from = c(all_of(power_cols), spectral_entropy),
              names_glue  = "{.value}__{source}")

common_ids_spharm <- spharm_wide %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  pull(ID)
spharm_wide <- spharm_wide %>% filter(ID %in% common_ids_spharm)
cat(sprintf("SPHARM 验证：%d 个标本，%d 阶\n\n",
            nrow(spharm_wide), length(power_cols)))

# --- 各阶 Bland-Altman，结果收进长表 ---
ba_power_summary <- tibble()

for (pair in pairs_label) {
  src_a      <- pair[1]
  src_b      <- pair[2]
  pair_label <- pair[3]   # 语义名称
  
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

power_compare <- spharm_wide %>%
  select(ID, matches("^power_l[1-5]__")) %>%
  arrange(ID)

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

power_all_compare <- spharm_wide %>%
  select(ID, matches("^power_l[0-9]+__")) %>%
  rename_with(~ str_replace(., "__", "_"),
              matches("^power_l[0-9]+__")) %>%
  arrange(ID)

write_csv(power_all_compare,
          here("analysis/data/derived_data/validation_values_power_all.csv"))
cat("已保存：validation_values_power_all.csv\n")


# ==============================================================================
# Part C：实证验证 — SVD 对齐 vs SVD 对齐 + 随机 Z 轴旋转
# ==============================================================================

cat("\n====== Part D: SVD vs SVD + 随机 Z 轴旋转 ======\n\n")

# --- 读取 svd_rotated 功率谱 ---
path_rotated <- here("analysis/data/derived_data/validation/svd_rotated/SPHARM_direction.csv")

if (!file.exists(path_rotated)) {
  stop(glue(
    "找不到：{path_rotated}\n",
    "请先运行：\n",
    "  1. python rotate_svd_directions.py\n",
    "  2. python kde_to_spharm_main.py --source svd_rotated"
  ))
}

df_svd_rotated <- read_csv(path_rotated, show_col_types = FALSE) %>%
  mutate(source = "svd_rotated") %>%
  filter(!str_starts(ID, "IM_"))

# 取与 svd 共有的标本
common_ids_rot <- intersect(
  df_svd %>% pull(ID),
  df_svd_rotated %>% pull(ID)
)
cat(sprintf("svd vs svd_rotated 验证：%d 个标本\n\n", length(common_ids_rot)))

# 宽格式合并
spharm_rot_wide <- bind_rows(
  df_svd         %>% filter(ID %in% common_ids_rot),
  df_svd_rotated %>% filter(ID %in% common_ids_rot)
) %>%
  select(ID, source, all_of(power_cols), spectral_entropy) %>%
  pivot_wider(
    names_from  = source,
    values_from = c(all_of(power_cols), spectral_entropy),
    names_glue  = "{.value}__{source}"
  )


# ------------------------------------------------------------------------------
# C-1：各阶功率谱 Bland-Altman
# ------------------------------------------------------------------------------

ba_rot_summary <- tibble()

for (pcol in power_cols) {
  col_svd <- glue("{pcol}__svd")
  col_rot <- glue("{pcol}__svd_rotated")
  ba      <- bland_altman_calc(spharm_rot_wide[[col_svd]],
                               spharm_rot_wide[[col_rot]])
  degree  <- as.integer(str_remove(pcol, "power_l"))
  
  ba_rot_summary <- bind_rows(ba_rot_summary, tibble(
    degree    = degree,
    bias      = ba$bias,
    loa_upper = ba$loa_upper,
    loa_lower = ba$loa_lower,
    sd_diff   = ba$sd_diff
  ))
  
  summary_table <- bind_rows(summary_table,
                             summary_row(pcol, "perturbed vs unperturbed", ba))
}

ba_rot_plot_df <- ba_rot_summary %>%
  mutate(pair = "perturbed vs unperturbed")

ba_combined <- bind_rows(
  ba_power_summary,
  ba_rot_plot_df
) %>%
  mutate(pair = factor(pair, levels = c(
    "none-aligned vs techno-aligned",
    "none-aligned vs morph-aligned",
    "techno-aligned vs morph-aligned",
    "perturbed vs unperturbed"
  )))

pair_colors <- c(
  "none-aligned vs techno-aligned" = "#802520",
  "none-aligned vs morph-aligned"  = "#BA8530",
  "techno-aligned vs morph-aligned"= "#5C7F71",
  "perturbed vs unperturbed"       = "#F5EDDC"
)

p_rot_spharm <- ggplot(ba_combined, aes(x = degree)) +
  geom_ribbon(aes(ymin = loa_lower, ymax = loa_upper,
                  fill = pair),
              alpha = 0.15) +
  geom_line(aes(y = bias, color = pair),
            linewidth = 1.25) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.2) +
  scale_color_manual(values = pair_colors) +
  scale_fill_manual(values  = pair_colors) +
  scale_x_continuous(breaks = seq(0, max(ba_combined$degree), by = 1)) +
  theme_bw(base_size = 10) +
  labs(
    x        = "SP-SPHARM power spectra degree (l)",
    y        = "Difference (bias ± 95% LoA)",
    color    = NULL,
    fill     = NULL
  ) +
  theme(
    panel.grid.major.x = element_line(color = "grey50", linewidth = 0.3, linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle   = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
    legend.position      = c(0.9, 0.2),
    legend.justification = c(0.5, 0.5),
    legend.background    = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )
p_rot_spharm

# ------------------------------------------------------------------------------
# C-3：数值汇总与结论
# ------------------------------------------------------------------------------

cat("==== Part D 数值汇总 ====\n")

loa_width_intermethod <- ba_power_summary %>%
  group_by(pair) %>%
  summarise(max_loa_width = max(loa_upper - loa_lower), .groups = "drop")

loa_width_intramethod <- ba_rot_summary %>%
  summarise(max_loa_width = max(loa_upper - loa_lower)) %>%
  mutate(pair = "perturbed vs unperturbed")

cat("\n各对比组最大 LoA 宽度（越小说明一致性越高）：\n")
bind_rows(loa_width_intermethod, loa_width_intramethod) %>%
  mutate(max_loa_width = formatC(max_loa_width, format = "e", digits = 3)) %>%
  print()

max_intra <- max(ba_rot_summary$loa_upper - ba_rot_summary$loa_lower)
max_inter <- max(ba_power_summary$loa_upper - ba_power_summary$loa_lower)

cat(sprintf(
  "\n坐标系内扰动 LoA 宽度（%.2e）%s 方法间差异 LoA 宽度（%.2e）\n",
  max_intra,
  ifelse(max_intra < max_inter, "<", "≥"),
  max_inter
))

if (max_intra < max_inter * 0.1) {
  cat("结论：坐标系内随机扰动远小于方法间系统差异（< 10%），\n")
  cat("      固定 SVD 对齐坐标系后，分析结果不受坐标系内随机误差影响。✓\n")
} else if (max_intra < max_inter) {
  cat("结论：坐标系内随机扰动小于方法间系统差异，\n")
  cat("      固定 SVD 对齐坐标系后，分析结果基本不受影响。\n")
} else {
  cat("结论：坐标系内随机扰动不可忽略，建议检查对齐质量或 KDE 参数。\n")
}

# 更新汇总表
write_csv(
  summary_table %>% arrange(metric, pair),
  here("analysis/data/derived_data/validation_ba_summary.csv")
)
cat("\n完整汇总表已更新：validation_ba_summary.csv\n")

# ==============================================================================
# Final combined plot
# ==============================================================================

p_rotational_invariance_validity <- p_rei / p_rot_spharm +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold")
    )
  )
p_rotational_invariance_validity

ggsave(
  here("analysis/output/figures/validation_combined.png"),
  plot = p_rotational_invariance_validity,
  width = 12,
  height = 14,
  dpi = 600,
  bg = "white"
)
