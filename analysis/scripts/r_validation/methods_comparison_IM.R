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
#     3. 功率谱折线图（每类理想模型单独输出一张图）
#
#   Part C — 三种方法的判别能力比较
#     4. 合并 R/E/I 与 SPHARM 功率谱数据
#     5. 成对欧氏距离热图：SPI vs Fabric vs SPHARM
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 1：IM 标本）
#   - analysis/data/derived_data/SPHARM_direction.csv
# 输出：
#   - analysis/output/figures/IM_Direction_Variance.png
#   - analysis/output/figures/PowerSpectrum_<type>.png（每类一张）
#   - analysis/output/figures/IM_3_methods_comparison.png
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(ggrepel)
library(patchwork)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(base::`%*%`)

# ==============================================================================
# 公共函数
# ==============================================================================

compute_R <- function(dx, dy, dz) {
  resultant_magnitude <- sqrt(sum(dx)^2 + sum(dy)^2 + sum(dz)^2)
  total_length        <- sum(sqrt(dx^2 + dy^2 + dz^2))
  resultant_magnitude / total_length
}

compute_EI <- function(ux, uy, uz) {
  n      <- length(ux)
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / n
  eig    <- eigen(T_mat, symmetric = TRUE)
  lambda <- sort(eig$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_)
  )
}

# ==============================================================================
# 1. 读取数据
# ==============================================================================

raw <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 1)

SPHARM_direction <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE)

# ==============================================================================
# Part A：方向统计量（R、E、I）— 全部标本
# ==============================================================================

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

results <- raw_dirs %>%
  group_by(ID) %>%
  summarise(
    n_scars = n(),
    R       = compute_R(ux, uy, uz),
    E       = compute_EI(ux, uy, uz)$E,
    I       = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  ) %>%
  arrange(ID)

results %>%
  mutate(across(c(R, E, I), \(x) round(x, 4))) %>%
  print(n = Inf)

# ==============================================================================
# Part B：SPHARM 功率谱分析 — 仅 IM 标本
# ==============================================================================

SPHARM_IM <- SPHARM_direction %>%
  filter(str_starts(ID, "IM_")) %>%
  select(ID, Typology, spectral_entropy, power_l1:power_l20)

cat("理想模型标本数：", nrow(SPHARM_IM), "\n")
print(SPHARM_IM$ID)

# --- B-1：各阶方差折线图 ---
variance_IM <- SPHARM_IM %>%
  select(starts_with("power_l")) %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(),
               names_to = "degree_label", values_to = "variance") %>%
  mutate(degree = as.integer(str_remove(degree_label, "power_l"))) %>%
  arrange(degree) %>%
  mutate(
    var_pct    = variance / sum(variance) * 100,
    var_cumsum = cumsum(var_pct)
  )

variance_IM %>%
  mutate(across(c(var_pct, var_cumsum), \(x) round(x, 2))) %>%
  print(n = Inf)

p_variance_IM <- ggplot(variance_IM, aes(x = degree, y = variance)) +
  geom_line(color = "#FFBAE0", linewidth = 1, alpha = 0.9) +
  geom_point(color = "#FFBAE0", size = 3, alpha = 0.9) +
  scale_x_continuous(breaks = seq(min(variance_IM$degree),
                                  max(variance_IM$degree), by = 1)) +
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

# --- B-2：功率谱折线图（每类单独输出）---
df_long <- SPHARM_IM %>%
  pivot_longer(
    cols      = starts_with("power_l"),
    names_to  = "degree_label",
    values_to = "power"
  ) %>%
  mutate(
    degree   = as.integer(str_remove(degree_label, "power_l")),
    ID_label = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

y_min <- min(df_long$power, na.rm = TRUE)
y_max <- max(df_long$power, na.rm = TRUE)

df_long %>%
  group_by(ID_label) %>%
  group_split() %>%
  walk(function(df_sub) {
    
    type_name <- as.character(df_sub$ID_label[1])
    
    p <- ggplot(df_sub, aes(x = degree, y = power)) +
      geom_line(color = "#B26538", linewidth = 0.5, alpha = 0.7) +
      geom_point(color = "#B26538", size = 2.5, shape = 16, alpha = 0.7) +
      scale_x_continuous(breaks = 1:20) +
      scale_y_continuous(
        labels = scales::label_scientific(digits = 2),
        limits = c(y_min, y_max)
      ) +
      theme_classic() +
      labs(
        title = type_name,
        x     = "Spherical Harmonic Degree (l)",
        y     = "Normalised Power"
      ) +
      theme(
        plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_text(size = 9),
        axis.text.y        = element_text(size = 9)
      )
    
    file_name <- type_name %>%
      str_replace_all(" ", "_") %>%
      str_to_lower()
    
    ggsave(
      here(paste0("analysis/output/figures/PowerSpectrum_", file_name, ".png")),
      plot = p, width = 6, height = 3, dpi = 600
    )
    
    cat("已保存：", type_name, "\n")
  })

# ==============================================================================
# Part C：三种方法的判别能力比较
# ==============================================================================

# --- C-1：IM 标本 Fabric（R/E/I）---
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
id_order <- c(
  "Cylindrical unipolar cortical", "Cylindrical unipolar scarred",
  "Cylindrical bipolar",
  "Conical unipolar cortical",     "Conical unipolar scarred",
  "Discoid",                       "Discoid unifacial",
  "Levallois preferential",        "Levallois convergent",
  "Levallois laminar",
  "Biface",                        "Multiplatform"
)

df_im <- SPHARM_IM %>%
  left_join(fabric_IM, by = "ID") %>%
  mutate(
    label = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

cat("===== 合并后数据集（n =", nrow(df_im), "）=====\n")
print(df_im %>% select(label, R, E, I, power_l1:power_l4), n = Inf)

# --- C-3：成对欧氏距离热图 ---
make_dist_df <- function(X, labels, method_name) {
  d <- as.matrix(dist(scale(X), method = "euclidean"))
  rownames(d) <- colnames(d) <- labels
  as.data.frame(d) %>%
    rownames_to_column("From") %>%
    pivot_longer(-From, names_to = "To", values_to = "distance") %>%
    mutate(method = method_name)
}

labs <- df_im$label

dist_all <- bind_rows(
  make_dist_df(df_im %>% select(R)                 %>% as.matrix(), labs, "SPI"),
  make_dist_df(df_im %>% select(E, I)              %>% as.matrix(), labs, "Fabric"),
  make_dist_df(df_im %>% select(power_l1:power_l4) %>% as.matrix(), labs, "SPHARM")
) %>%
  mutate(
    From   = factor(From,   levels = id_order),
    To     = factor(To,     levels = id_order),
    method = factor(method, levels = c("SPI", "Fabric", "SPHARM"))
  )

avg_dist <- dist_all %>%
  filter(From != To) %>%
  group_by(method) %>%
  summarise(mean_dist = round(mean(distance), 2), .groups = "drop")

cat("\n===== 各方法的标准化平均类间距离 =====\n")
cat("（越大 = 整体判别能力越强）\n")
print(avg_dist)

dist_all_upper <- dist_all %>%
  dplyr::filter(as.numeric(From) < as.numeric(To))

p_heatmap <- ggplot(dist_all_upper, aes(x = To, y = From, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", distance)), size = 2.2) +
  facet_wrap(~ method, ncol = 3) +
  scale_fill_gradient2(
    low      = "#5C7F71",
    mid      = "#F5EDDC",
    high     = "#802520",
    midpoint = 2,
    name     = "Euclidean\ndistance\n(standardised)"
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