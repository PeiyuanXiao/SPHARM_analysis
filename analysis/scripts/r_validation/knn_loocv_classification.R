# ==============================================================================
# silhouette_analysis.R
# 实验石核片疤模式分析方法有效性对比
# 方法：轮廓系数（Silhouette Score）
#
# 研究问题：
#   三种方法（SPI、Fabric、SPHARM）能否使同类石核的片疤模式
#   在特征空间中彼此接近、与异类石核保持距离？
#   → 轮廓系数直接量化类内凝聚度与类间分离度的比值
#
# 与理想模型分析的关系：
#   理想模型（每类 n=1）→ 成对类间距离热图（无法计算类内距离）
#   实验标本（每类 n≥2）→ 轮廓系数（同时量化类内和类间距离）
#   两者互补，共同构成完整的方法有效性论证
#
# 距离度量：
#   SPI           → 欧氏距离（标准化后）
#   Fabric（E+I） → 马氏距离（处理 E/I 强负相关 r = −0.731）
#   SPHARM        → 余弦距离
#
# 输入：
#   - analysis/data/raw_data/Scar_orientation_data.xlsx（sheet 3：实验石核）
#   - analysis/data/derived_data/SPHARM_direction_lin2024.csv
# 输出：
#   - analysis/output/figures/silhouette_overall.png   总体轮廓系数柱状图
#   - analysis/output/figures/silhouette_bytype.png    各类型轮廓系数热图
#   - analysis/output/figures/silhouette_violin.png    各标本轮廓系数分布图
#   - analysis/data/derived_data/silhouette_results.csv
# ==============================================================================

library(here)
library(tidyverse)
library(readxl)
library(patchwork)
library(glue)
library(MASS)   # ginv()，马氏距离协方差矩阵求逆备用

exp_data             <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 3)
exp_SPHARM_direction <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction_lin2024.csv"),
  show_col_types = FALSE)


# ==============================================================================
# 公共函数：方向统计量计算
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


# ==============================================================================
# 距离矩阵计算函数
# ==============================================================================

# 欧氏距离（z-score 标准化后）
dist_euclidean <- function(X) {
  X_scaled <- scale(X)
  as.matrix(dist(X_scaled, method = "euclidean"))
}

# 余弦距离：1 - cos(x, y)
dist_cosine <- function(X) {
  X   <- as.matrix(X)
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  1 - sim
}

# 马氏距离（用全局协方差矩阵）
dist_mahalanobis <- function(X) {
  X     <- as.matrix(X)
  S     <- cov(X)
  S_inv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  n     <- nrow(X)
  D     <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      d        <- X[i, ] - X[j, ]
      D[i, j]  <- sqrt(max(0, t(d) %*% S_inv %*% d))
    }
  }
  D
}


# ==============================================================================
# 轮廓系数计算函数
# ==============================================================================

# 对每件标本计算轮廓系数
# s(i) = (b(i) - a(i)) / max(a(i), b(i))
#   a(i)：与同类其他标本的平均距离（类内凝聚度）
#   b(i)：与最近异类所有标本的平均距离（类间分离度）

compute_silhouette <- function(D, labels) {
  n      <- nrow(D)
  types  <- unique(labels)
  s      <- numeric(n)
  a      <- numeric(n)
  b      <- numeric(n)
  
  for (i in seq_len(n)) {
    same_class  <- which(labels == labels[i] & seq_len(n) != i)
    other_types <- types[types != labels[i]]
    
    # a(i)：类内平均距离
    if (length(same_class) == 0) {
      # 单件类型无法计算，设为 NA
      s[i] <- NA_real_
      a[i] <- NA_real_
      b[i] <- NA_real_
      next
    }
    a[i] <- mean(D[i, same_class])
    
    # b(i)：最近异类的平均距离
    b_candidates <- sapply(other_types, function(t) {
      idx <- which(labels == t)
      mean(D[i, idx])
    })
    b[i] <- min(b_candidates)
    
    # 轮廓系数
    s[i] <- (b[i] - a[i]) / max(a[i], b[i])
  }
  
  data.frame(
    label = labels,
    a     = a,
    b     = b,
    s     = s
  )
}


# ==============================================================================
# 1. 计算 SPI 和 Fabric
# ==============================================================================

raw_dirs <- exp_data %>%
  mutate(
    dx  = End_X - Start_X,
    dy  = End_Y - Start_Y,
    dz  = End_Z - Start_Z,
    len = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(len > 1e-10) %>%
  mutate(ux = dx / len, uy = dy / len, uz = dz / len)

rei_exp <- raw_dirs %>%
  group_by(ID, Typology) %>%
  summarise(
    SPI = compute_R(ux, uy, uz),
    E   = compute_EI(ux, uy, uz)$E,
    I   = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

cat(sprintf("E 与 I 的相关系数：%.3f（Fabric 使用马氏距离）\n\n",
            cor(rei_exp$E, rei_exp$I, use = "complete.obs")))


# ==============================================================================
# 2. 合并 SPHARM 数据
# ==============================================================================

spharm_exp <- exp_SPHARM_direction %>%
  select(ID, power_l1:power_l5)

all_data <- rei_exp %>%
  left_join(spharm_exp, by = "ID") %>%
  filter(!is.na(power_l1))


# ==============================================================================
# 3. 剔除单件类型（轮廓系数需要每类至少 2 件）
# ==============================================================================

single_types <- all_data %>% count(Typology) %>%
  filter(n < 2) %>% pull(Typology)

cat("剔除单件类型：", paste(single_types, collapse = ", "), "\n\n")

df <- all_data %>% filter(!Typology %in% single_types)

cat(sprintf("最终：%d 件标本，%d 种类型\n\n",
            nrow(df), n_distinct(df$Typology)))
df %>% count(Typology, sort = TRUE) %>% print()


# ==============================================================================
# 4. 计算三种方法的距离矩阵
# ==============================================================================

labels <- df$Typology

X_spi <- df %>% select(SPI)           %>% as.matrix()
X_fab <- df %>% select(E, I)          %>% as.matrix()
X_sph <- df %>% select(power_l1:power_l5) %>% as.matrix()

cat("\n计算距离矩阵...\n")
D_spi <- dist_euclidean(X_spi)
cat("  SPI（欧氏）：完成\n")
D_fab <- dist_mahalanobis(X_fab)
cat("  Fabric（马氏）：完成\n")
D_sph <- dist_cosine(X_sph)
cat("  SPHARM（余弦）：完成\n\n")


# ==============================================================================
# 5. 计算轮廓系数
# ==============================================================================

sil_spi <- compute_silhouette(D_spi, labels) %>%
  mutate(method = "SPI",            ID = df$ID)
sil_fab <- compute_silhouette(D_fab, labels) %>%
  mutate(method = "Fabric (E + I)", ID = df$ID)
sil_sph <- compute_silhouette(D_sph, labels) %>%
  mutate(method = "SPHARM (l1–l5)", ID = df$ID)

sil_all <- bind_rows(sil_spi, sil_fab, sil_sph)


# ==============================================================================
# 6. 数值汇总
# ==============================================================================

# 总体轮廓系数（各标本的均值，忽略单件类型 NA）
overall_sil <- sil_all %>%
  group_by(method) %>%
  summarise(
    mean_s   = round(mean(s, na.rm = TRUE), 4),
    median_s = round(median(s, na.rm = TRUE), 4),
    n_pos    = sum(s > 0, na.rm = TRUE),
    n_total  = sum(!is.na(s)),
    pct_pos  = round(mean(s > 0, na.rm = TRUE) * 100, 1),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean_s))

cat("====== 总体轮廓系数汇总 ======\n\n")
overall_sil %>%
  mutate(
    dist_used = case_when(
      method == "SPI"            ~ "euclidean",
      method == "Fabric (E + I)" ~ "mahalanobis",
      TRUE                       ~ "cosine"
    )
  ) %>%
  select(method, dist_used, mean_s, median_s, n_pos, n_total, pct_pos) %>%
  print()

cat("\n说明：\n")
cat("  mean_s   = 平均轮廓系数（越高越好，>0 表示类内紧凑度优于类间分离度）\n")
cat("  median_s = 中位轮廓系数（对异常值更稳健）\n")
cat("  pct_pos  = 轮廓系数 > 0 的标本比例（越高越好）\n\n")

# 各类型平均轮廓系数
cat("====== 各类型平均轮廓系数 ======\n\n")
sil_all %>%
  group_by(method, label) %>%
  summarise(mean_s = round(mean(s, na.rm = TRUE), 3),
            n = n(), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_s,
              names_sort = TRUE) %>%
  arrange(label) %>%
  print(n = Inf)


# ==============================================================================
# 7. 可视化 A：总体轮廓系数柱状图
# ==============================================================================

method_order  <- c("SPI", "Fabric (E + I)", "SPHARM (l1–l5)")
method_colors <- c(
  "SPI"            = "#A1C2E6",
  "Fabric (E + I)" = "#FFBAE0",
  "SPHARM (l1–l5)" = "#D4619A"
)

p_overall <- overall_sil %>%
  mutate(method = factor(method, levels = method_order)) %>%
  ggplot(aes(x = method, y = mean_s, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_col(width = 0.55, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.4f", mean_s),
                vjust = ifelse(mean_s >= 0, -0.4, 1.4)),
            size = 3.8, fontface = "bold") +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(breaks = seq(-1, 1, 0.05)) +
  theme_bw(base_size = 11) +
  labs(
    title    = "Mean Silhouette Score by Method",
    subtitle = paste("SPI: Euclidean | Fabric: Mahalanobis | SPHARM: Cosine",
                     "\nPositive = intra-class distance < inter-class distance"),
    x        = NULL,
    y        = "Mean silhouette score"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 10)
  )


# ==============================================================================
# 8. 可视化 B：各类型轮廓系数热图
# ==============================================================================

type_sil <- sil_all %>%
  group_by(method, label) %>%
  summarise(mean_s = round(mean(s, na.rm = TRUE), 3),
            .groups = "drop") %>%
  mutate(method = factor(method, levels = method_order))

p_bytype <- ggplot(type_sil,
                   aes(x = method, y = fct_rev(label), fill = mean_s)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label  = sprintf("%.3f", mean_s),
                color  = ifelse(abs(mean_s) > 0.3, "white", "grey20")),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(
    low      = "#3B8BD4",
    mid      = "white",
    high     = "#D4619A",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Silhouette\nscore"
  ) +
  scale_color_identity() +
  theme_bw(base_size = 10) +
  labs(
    title    = "Silhouette Score by Type and Method",
    subtitle = "Pink = intra-class compact (good) | Blue = inter-class closer (poor)",
    x        = NULL,
    y        = NULL
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 10, face = "bold"),
    axis.text.y   = element_text(size = 9)
  )


# ==============================================================================
# 9. 可视化 C：各标本轮廓系数分布（小提琴图 + 抖动点）
# ==============================================================================

p_violin <- sil_all %>%
  filter(!is.na(s)) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  ggplot(aes(x = method, y = s, fill = method, color = method)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_violin(alpha = 0.3, linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.4, linewidth = 0.8,
               color = "grey20") +
  scale_fill_manual(values  = method_colors, guide = "none") +
  scale_color_manual(values = method_colors, guide = "none") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  theme_bw(base_size = 11) +
  labs(
    title    = "Distribution of Silhouette Scores (per specimen)",
    subtitle = "Crossbar = mean | Points = individual specimens",
    x        = NULL,
    y        = "Silhouette score"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(size = 10)
  )


# ==============================================================================
# 10. 保存图片
# ==============================================================================

ggsave(here("analysis/output/figures/silhouette_overall.png"),
       plot = p_overall, width = 7, height = 5, dpi = 300, bg = "white")
cat("图已保存：silhouette_overall.png\n")

ggsave(here("analysis/output/figures/silhouette_bytype.png"),
       plot = p_bytype, width = 8, height = 6, dpi = 300, bg = "white")
cat("图已保存：silhouette_bytype.png\n")

ggsave(here("analysis/output/figures/silhouette_violin.png"),
       plot = p_violin, width = 7, height = 5, dpi = 300, bg = "white")
cat("图已保存：silhouette_violin.png\n")


# ==============================================================================
# 11. 保存详细结果
# ==============================================================================

sil_all %>%
  select(method, ID, label, a, b, s) %>%
  arrange(method, label, ID) %>%
  write_csv(here("analysis/data/derived_data/silhouette_results.csv"))
cat("详细结果已保存：silhouette_results.csv\n")


# ==============================================================================
# 12. 结论
# ==============================================================================

cat("\n====== 结论 ======\n\n")
best <- overall_sil %>% slice_max(mean_s, n = 1)
worst <- overall_sil %>% slice_min(mean_s, n = 1)

cat(sprintf(
  "平均轮廓系数最高：%s（%.4f）\n",
  best$method, best$mean_s
))
cat(sprintf(
  "平均轮廓系数最低：%s（%.4f）\n",
  worst$method, worst$mean_s
))
cat("\n方法排序（均值轮廓系数）：\n")
cat(paste(overall_sil$method,
          sprintf("%.4f", overall_sil$mean_s),
          sep = " = ", collapse = " > "), "\n\n")

# 各方法轮廓系数 > 0 的标本比例
cat("轮廓系数 > 0 的标本比例（类内距离 < 类间距离）：\n")
overall_sil %>%
  select(method, pct_pos, n_pos, n_total) %>%
  mutate(result = glue("{pct_pos}% ({n_pos}/{n_total})")) %>%
  select(method, result) %>%
  print()