library(tidyverse)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
raw <- read_excel("analysis/data/raw_data/Scar_orientation_data.xlsx", sheet = 1)

# 计算单位方向向量==============================================================
raw <- raw |>
  mutate(
    dx     = End_X - Start_X,
    dy     = End_Y - Start_Y,
    dz     = End_Z - Start_Z,
    length = sqrt(dx^2 + dy^2 + dz^2)
  ) |>
  filter(length > 1e-10) |>  
  mutate(
    ux = dx / length,
    uy = dy / length,
    uz = dz / length
  )

# Mean Resultant Length R ======================================================
compute_R <- function(ux, uy, uz) {
  n  <- length(ux)
  Rx <- sum(ux)
  Ry <- sum(uy)
  Rz <- sum(uz)
  sqrt(Rx^2 + Ry^2 + Rz^2) / n
}

# Elongation (E) & Isotropy (I) ================================================
compute_EI <- function(ux, uy, uz) {
  n <- length(ux)
  U <- cbind(ux, uy, uz)             # n × 3 矩阵
  
  # 方向张量（外积的均值）
  T_mat <- (t(U) %*% U) / n         # 3 × 3
  
  # 特征值（降序）
  eig    <- eigen(T_mat, symmetric = TRUE)
  lambda <- sort(eig$values, decreasing = TRUE)   # λ1 ≥ λ2 ≥ λ3
  
  # 防止数值误差导致负特征值
  lambda <- pmax(lambda, 0)
  
  E <- ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_)
  I <- ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_)
  
  list(
    E       = E,
    I       = I,
    lambda1 = lambda[1],
    lambda2 = lambda[2],
    lambda3 = lambda[3]
  )
}

# 批量计算
results <- raw |>
  group_by(ID) |>
  summarise(
    n_scars = n(),
    
    # R
    R = compute_R(ux, uy, uz),
    
    # E 和 I（从列表中展开）
    EI      = list(compute_EI(ux, uy, uz)),
    .groups = "drop"
  ) |>
  mutate(
    E       = map_dbl(EI, "E"),
    I       = map_dbl(EI, "I"),
    lambda1 = map_dbl(EI, "lambda1"),
    lambda2 = map_dbl(EI, "lambda2"),
    lambda3 = map_dbl(EI, "lambda3")
  ) |>
  select(-EI) |>
  arrange(ID)

results |>
  select(ID, n_scars, R, E, I) |>
  mutate(across(c(R, E, I), \(x) round(x, 4))) |>
  print(n = Inf)

# 理想模型片疤方向SPHARM========================================================
SPHARM_direction <- read_csv("analysis/data/derived_data/SPHARM_direction.csv")

SPHARM_direction_IM <- SPHARM_direction %>%
  filter(str_starts(ID, "IM_"))

cat("理想模型标本数：", nrow(SPHARM_direction_IM), "\n")
cat("包含类型：\n")
print(SPHARM_direction_IM$ID)

# 基于理想模型计算每阶方差
power_cols_IM <- SPHARM_direction_IM %>%
  select(starts_with("power_l")) %>%
  colnames()

variance_IM <- SPHARM_direction_IM %>%
  select(all_of(power_cols_IM)) %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "degree_label",
    values_to = "variance"
  ) %>%
  mutate(
    degree = as.integer(str_remove(degree_label, "power_l"))
  ) %>%
  arrange(degree) %>%
  select(degree, variance)

print(variance_IM, n = Inf)

# 方差折线图
variance_IM <- variance_IM %>%
  mutate(
    var_pct  = variance / sum(variance) * 100,
    var_cumsum = cumsum(var_pct)
  )

variance_IM %>%
  select(degree, variance, var_pct, var_cumsum) %>%
  mutate(across(c(var_pct, var_cumsum), \(x) round(x, 2))) %>%
  print(n = Inf)

# 绘图
p_variance_IM <- ggplot(variance_IM, aes(x = degree, y = variance)) +
  geom_line(color  = "#FFBAE0", linewidth = 1, alpha = 0.9) +
  geom_point(color = "#FFBAE0", size = 3,      alpha = 0.9) +
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
    plot.title       = element_text(face = "bold", size = 10, hjust = 0.5),
    panel.grid.minor.x = element_blank()
  )

print(p_variance_IM)

ggsave(
  "analysis/data/derived_data/IM_Direction_Variance.png",
  plot   = p_variance_IM,
  width  = 8,
  height = 5,
  dpi    = 300
)

# 筛选阶数
SPHARM_direction_IM_filter <- SPHARM_direction_IM %>%
  select(
    ID,
    Typology,
    SHE,
    spectral_entropy,
    power_l1:power_l4
  )

print(SPHARM_direction_IM_filter, n = Inf)


df_long <- SPHARM_direction_IM_filter %>%
  select(ID, Typology, power_l1:power_l4) %>%
  pivot_longer(
    cols         = power_l1:power_l4,
    names_to     = "degree_label",
    values_to    = "power"
  ) %>%
  mutate(
    degree = as.integer(str_remove(degree_label, "power_l")),
    # 清理标签：去掉 IM_ 前缀，下划线换空格，方便图中显示
    ID_clean = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

# 固定分面顺序（按类型逻辑分组排列）
id_order <- c(
  "Cylindrical unipolar cortical",
  "Cylindrical unipolar scarred",
  "Cylindrical bipolar",
  "Conical unipolar cortical",
  "Conical unipolar scarred",
  "Discoid",
  "Levallois preferential",
  "Levallois convergent",
  "Levallois laminar",
  "Biface",
  "Multiplatform"
)

df_long <- df_long %>%
  mutate(ID_clean = factor(ID_clean, levels = id_order))

# ==============================================================================
# 2. 分面折线图
# ==============================================================================

p <- ggplot(df_long, aes(x = degree, y = power)) +
  
  # 填充面积，强调功率分布形状
  geom_line(color = "#D4619A", linewidth = 0.9) +
  geom_point(color = "#D4619A", size = 2.2, shape = 16) +
  
  # 每个标本一个分面
  facet_wrap(
    ~ ID_clean,
    ncol   = 3,
    scales = "free_y"     
  ) +
  
  scale_x_continuous(
    breaks = 1:4,
    labels = paste0("l=", 1:4)
  ) +
  scale_y_continuous(
    labels = scales::label_scientific(digits = 2)
  ) +
  
  theme_bw(base_size = 10) +
  labs(
    title    = "Direction SPHARM Power Spectrum in Ideal Models",
    y        = "Normalised Power"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey50", margin = margin(b = 8)),
    strip.text       = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey80", color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 7)
  )

print(p)

ggsave(
  "analysis/data/derived_data/IM_Direction_PowerSpectrum_Facet.png",
  plot   = p,
  width  = 10,
  height = 10,
  dpi    = 300
)

df_long %>%
  select(ID_clean, degree, power) %>%
  pivot_wider(names_from = degree, values_from = power,
              names_prefix = "l=") %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)



# 留一法模型验证=======================================================
install.packages(c("MASS", "caret", "patchwork", "ggrepel"))
library(MASS)     
library(caret)     
library(ggplot2)
library(patchwork)

# ==============================================================================
# 0. 读取数据
# ==============================================================================

raw <- read_excel(
  "H:/SPHARM_analysis/analysis/data/raw_data/Scar_orientation_data.xlsx",
  sheet = 1
)


# ==============================================================================
# 1. 计算 R / E / I
# ==============================================================================

compute_R <- function(ux, uy, uz)
  sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2) / length(ux)

compute_EI <- function(ux, uy, uz) {
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / nrow(U)
  lambda <- sort(eigen(T_mat, symmetric = TRUE)$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E = ifelse(lambda[1] > 1e-10, 1 - lambda[2]/lambda[1], NA_real_),
    I = ifelse(lambda[1] > 1e-10,     lambda[3]/lambda[1], NA_real_)
  )
}

fabric_IM <- raw %>%
  filter(str_starts(ID, "IM_")) %>%
  mutate(
    dx = End_X - Start_X, dy = End_Y - Start_Y, dz = End_Z - Start_Z,
    len = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(len > 1e-10) %>%
  mutate(ux = dx/len, uy = dy/len, uz = dz/len) %>%
  group_by(ID) %>%
  summarise(
    R = compute_R(ux, uy, uz),
    E = compute_EI(ux, uy, uz)$E,
    I = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

# ==============================================================================
# 2. 合并数据集
# ==============================================================================

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

# ==============================================================================
# 4. 成对欧氏距离热图
# 对每套特征分别做 z-score 标准化后计算成对距离
# 距离越大 = 这两个类型在该方法下越容易区分
# ==============================================================================

make_dist_df <- function(X, labels, method_name) {
  X_scaled <- scale(X)
  d        <- as.matrix(dist(X_scaled, method = "euclidean"))
  rownames(d) <- colnames(d) <- labels
  as.data.frame(d) %>%
    rownames_to_column("From") %>%
    pivot_longer(-From, names_to = "To", values_to = "distance") %>%
    mutate(method = method_name)
}

# 标签顺序（按类型逻辑排列）
label_order <- c(
  "Cylindrical unipolar cortical", "Cylindrical unipolar scarred",
  "Cylindrical bipolar",
  "Conical unipolar cortical", "Conical unipolar scarred",
  "Discoid",
  "Levallois preferential", "Levallois convergent", "Levallois laminar",
  "Biface", "Multiplatform"
)

labs <- df$label

X_R_mat     <- df %>% select(R)               %>% as.matrix()
X_EI_mat    <- df %>% select(E, I)            %>% as.matrix()
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

# 每种方法的平均类间距离
avg_dist <- dist_all %>%
  filter(From != To) %>%
  group_by(method) %>%
  summarise(mean_dist = round(mean(distance), 3), .groups = "drop")

cat("\n===== 各方法的标准化平均类间距离 =====\n")
cat("（越大 = 整体判别能力越强）\n")
print(avg_dist)

# 热图
p_heatmap <- ggplot(dist_all,
                    aes(x = To, y = From, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(From != To, sprintf("%.1f", distance), "")),
            size = 2.2) +
  facet_wrap(~ method, ncol = 3) +
  scale_fill_gradient(low = "white", high = "#E84D89",
                      name = "Euclidean\ndistance\n(standardised)") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 8) +
  labs(
    x = NULL, y = NULL
  ) +
  theme(
    plot.title   = element_text(face = "bold", size = 10, hjust = 0.5),
    strip.text   = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey80", color = "grey80"),
    axis.text    = element_text(size = 6.5)
  )

ggsave(
  "analysis/data/derived_data/IM_3_methods_comparison.png",
  plot = p_heatmap, width = 12, height = 5, dpi = 300
)


