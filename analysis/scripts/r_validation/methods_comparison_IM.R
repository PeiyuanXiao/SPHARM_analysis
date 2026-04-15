# ==============================================================================
# methods_comparison_IM.R
# 理想模型（IM）方法有效性验证
#
# 分析流程：
#   Part A — 方向统计量计算
#     1. 计算所有标本的 R（Mean Resultant Length）、E（Elongation）、I（Isotropy）
#
#   Part B — SPHARM 功率谱分析（仅 IM 标本）
#     2. 各阶方差折线图
#     3. 功率谱分面折线图（每类理想模型一个面板）
#
#   Part C — 三种方法的判别能力比较
#     4. 合并 R/E/I 与 SPHARM 功率谱数据
#     5. 成对欧氏距离热图：R vs E+I vs Power l1–l4
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 1：IM 标本）
#   - analysis/data/derived_data/SPHARM_direction.csv
# 输出：
#   - analysis/output/figures/IM_Direction_Variance.png
#   - analysis/output/figures/IM_Direction_PowerSpectrum_Facet.png
#   - analysis/output/figures/IM_3_methods_comparison.png
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(ggrepel)
library(patchwork)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(base::`%*%`)
conflicted::conflicts_prefer(var::stats)

# ==============================================================================
# 公共函数（整个脚本共用，避免重复定义）
# ==============================================================================

# compute_R()
# 计算 Mean Resultant Length R（Clarkson 法）
# 分子：合力向量的模；分母：所有原始刮痕向量的模长之和
compute_R <- function(dx, dy, dz) {
  resultant_magnitude <- sqrt(sum(dx)^2 + sum(dy)^2 + sum(dz)^2)
  total_length        <- sum(sqrt(dx^2 + dy^2 + dz^2))
  resultant_magnitude / total_length
}

# compute_EI()
# 计算 Elongation（E）和 Isotropy（I）
# 基于方向张量的特征值：λ1 ≥ λ2 ≥ λ3
# 返回：E、I 及三个特征值（供诊断用）
compute_EI <- function(ux, uy, uz) {
  n      <- length(ux)
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / n
  eig    <- eigen(T_mat, symmetric = TRUE)
  lambda <- sort(eig$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)   # 防止数值误差导致负特征值
  list(
    E       = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I       = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_),
    lambda1 = lambda[1],
    lambda2 = lambda[2],
    lambda3 = lambda[3]
  )
}


# ==============================================================================
# 1. 读取数据
# ==============================================================================

raw          <- read_excel(here("analysis/data/raw_data/Scar_orientation_data.xlsx"),
                           sheet = 1)
SPHARM_direction <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"))


# ==============================================================================
# Part A：方向统计量（R、E、I）— 全部标本
# ==============================================================================

# 计算单位方向向量（保留原始 raw 不变，结果存入新对象）
raw_dirs <- raw %>%
  mutate(
    dx     = End_X - Start_X,
    dy     = End_Y - Start_Y,
    dz     = End_Z - Start_Z,
    length = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(length > 1e-10) %>%
  mutate(
    ux = dx / length,
    uy = dy / length,
    uz = dz / length
  )

# 批量计算 R、E、I
results <- raw_dirs %>%
  group_by(ID) %>%
  summarise(
    n_scars = n(),
    R       = compute_R(ux, uy, uz),
    EI      = list(compute_EI(ux, uy, uz)),
    .groups = "drop"
  ) %>%
  mutate(
    E       = map_dbl(EI, "E"),
    I       = map_dbl(EI, "I"),
    lambda1 = map_dbl(EI, "lambda1"),
    lambda2 = map_dbl(EI, "lambda2"),
    lambda3 = map_dbl(EI, "lambda3")
  ) %>%
  select(-EI) %>%
  arrange(ID)

results %>%
  select(ID, n_scars, R, E, I) %>%
  mutate(across(c(R, E, I), \(x) round(x, 4))) %>%
  print(n = Inf)


# ==============================================================================
# Part B：SPHARM 功率谱分析 — 仅 IM 标本
# ==============================================================================

SPHARM_direction_IM <- SPHARM_direction %>%
  filter(str_starts(ID, "IM_"))

cat("理想模型标本数：", nrow(SPHARM_direction_IM), "\n")
cat("包含类型：\n")
print(SPHARM_direction_IM$ID)

# --- B-1：各阶方差折线图 ---
power_cols_IM <- SPHARM_direction_IM %>%
  select(starts_with("power_l")) %>%
  colnames()

variance_IM <- SPHARM_direction_IM %>%
  select(all_of(power_cols_IM)) %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(),
               names_to = "degree_label", values_to = "variance") %>%
  mutate(degree = as.integer(str_remove(degree_label, "power_l"))) %>%
  arrange(degree) %>%
  select(degree, variance) %>%
  mutate(
    var_pct    = variance / sum(variance) * 100,
    var_cumsum = cumsum(var_pct)
  )

variance_IM %>%
  select(degree, variance, var_pct, var_cumsum) %>%
  mutate(across(c(var_pct, var_cumsum), \(x) round(x, 2))) %>%
  print(n = Inf)

p_variance_IM <- ggplot(variance_IM, aes(x = degree, y = variance)) +
  geom_line(color = "#FFBAE0", linewidth = 1, alpha = 0.9) +
  geom_point(color = "#FFBAE0", size = 3, alpha = 0.9) +
  scale_x_continuous(
    breaks = seq(min(variance_IM$degree), max(variance_IM$degree), by = 1)
  ) +
  scale_y_continuous(labels = scales::scientific) +
  theme_bw() +
  labs(
    title = "Direction SPHARM Variance per Degree\n(Ideal Models only)",
    x     = "Spherical Harmonic Degree (l)",
    y     = "Variance"
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 10, hjust = 0.5),
    panel.grid.minor.x = element_blank()
  )

print(p_variance_IM)
ggsave(here("analysis/output/figures/IM_Direction_Variance.png"),
       plot = p_variance_IM, width = 8, height = 5, dpi = 300)

# --- B-2：功率谱分面折线图 ---
SPHARM_direction_IM_filter <- SPHARM_direction_IM %>%
  select(ID, Typology, SHE, spectral_entropy, power_l1:power_l4)

print(SPHARM_direction_IM_filter, n = Inf)

id_order <- c(
  "Cylindrical unipolar cortical", "Cylindrical unipolar scarred",
  "Cylindrical bipolar",
  "Conical unipolar cortical",     "Conical unipolar scarred",
  "Discoid",                       "Discoid unifacial",
  "Levallois preferential",        "Levallois convergent",
  "Levallois laminar",
  "Biface",                        "Multiplatform"
)

df_long <- SPHARM_direction_IM_filter %>%
  select(ID, Typology, power_l1:power_l4) %>%
  pivot_longer(cols = power_l1:power_l4,
               names_to = "degree_label", values_to = "power") %>%
  mutate(
    degree   = as.integer(str_remove(degree_label, "power_l")),
    ID_clean = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence() %>%
      factor(levels = id_order)
  )

p_facet <- ggplot(df_long, aes(x = degree, y = power)) +
  geom_area(fill = "#D4619A", alpha = 0.1) +
  geom_line(color = "#D4619A", linewidth = 0.2) +
  geom_point(fill = "#D4619A", color = "white",
             size = 2, shape = 21, stroke = 0.3) +
  facet_wrap(~ ID_clean, ncol = 3) +
  scale_x_continuous(breaks = 1:4, labels = paste0("l=", 1:4)) +
  scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
  theme_bw(base_size = 10) +
  labs(title = "Direction SPHARM Power Spectrum in Ideal Models",
       y = "Normalised Power") +
  theme(
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle      = element_text(size = 8, hjust = 0.5,
                                      color = "grey50", margin = margin(b = 8)),
    strip.text         = element_text(face = "bold", size = 8.5),
    strip.background   = element_rect(fill = "grey80", color = "grey80"),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 8),
    axis.text.y        = element_text(size = 7)
  )

print(p_facet)
ggsave(here("analysis/output/figures/IM_Direction_PowerSpectrum_Facet.png"),
       plot = p_facet, width = 8, height = 6, dpi = 600)

df_long %>%
  select(ID_clean, degree, power) %>%
  pivot_wider(names_from = degree, values_from = power,
              names_prefix = "l=") %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)


# ==============================================================================
# Part C：三种方法的判别能力比较
# ==============================================================================

# --- C-1：计算 IM 标本的 R / E / I ---
fabric_IM <- raw %>%
  filter(str_starts(ID, "IM_")) %>%
  mutate(
    dx  = End_X - Start_X,
    dy  = End_Y - Start_Y,
    dz  = End_Z - Start_Z,
    len = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(len > 1e-10) %>%
  mutate(ux = dx / len, uy = dy / len, uz = dz / len) %>%
  group_by(ID) %>%
  summarise(
    R = compute_R(ux, uy, uz),
    E = compute_EI(ux, uy, uz)$E,
    I = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

# --- C-2：合并数据集 ---
df <- SPHARM_direction_IM_filter %>%
  select(ID, Typology, spectral_entropy, power_l1:power_l4) %>%
  left_join(fabric_IM, by = "ID") %>%
  mutate(
    label = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

cat("===== 合并后数据集（n =", nrow(df), "）=====\n")
print(df %>% select(label, R, E, I, power_l1:power_l4), n = Inf)

# --- C-3：成对欧氏距离热图 ---
label_order <- id_order   # 复用 Part B 中定义的排列顺序

make_dist_df <- function(X, labels, method_name) {
  X_scaled <- scale(X)
  d        <- as.matrix(dist(X_scaled, method = "euclidean"))
  rownames(d) <- colnames(d) <- labels
  as.data.frame(d) %>%
    rownames_to_column("From") %>%
    pivot_longer(-From, names_to = "To", values_to = "distance") %>%
    mutate(method = method_name)
}

labs        <- df$label
X_R_mat     <- df %>% select(R)                 %>% as.matrix()
X_EI_mat    <- df %>% select(E, I)              %>% as.matrix()
X_power_mat <- df %>% select(power_l1:power_l4) %>% as.matrix()

dist_all <- bind_rows(
  make_dist_df(X_R_mat,     labs, "R"),
  make_dist_df(X_EI_mat,    labs, "E + I"),
  make_dist_df(X_power_mat, labs, "Power l1–l4")
) %>%
  mutate(
    From   = factor(From,   levels = label_order),
    To     = factor(To,     levels = label_order),
    method = factor(method, levels = c("R", "E + I", "Power l1–l4"))
  )

avg_dist <- dist_all %>%
  filter(From != To) %>%
  group_by(method) %>%
  summarise(mean_dist = round(mean(distance), 3), .groups = "drop")

cat("\n===== 各方法的标准化平均类间距离 =====\n")
cat("（越大 = 整体判别能力越强）\n")
print(avg_dist)

dist_all_upper <- dist_all %>%
  dplyr::filter(as.numeric(From) < as.numeric(To))

p_heatmap <- ggplot(dist_all_upper, aes(x = To, y = From, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(From != To, sprintf("%.1f", distance), "")),
    size = 2.2
  ) +
  facet_wrap(~ method, ncol = 3) +
  scale_fill_gradient2(
    low  = "#5C7F71",  
    mid  = "white", 
    high = "#802520",
    midpoint = 2,
    name = "Euclidean\ndistance\n(standardised)"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 8) +
  labs(x = NULL, y = NULL) +
  theme(
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#EBEBEB", color = "#EBEBEB"),
    axis.text        = element_text(size = 6.5)
  )

print(p_heatmap)
ggsave(here("analysis/output/figures/IM_3_methods_comparison.png"),
       plot = p_heatmap, width = 12, height = 5, dpi = 300)
