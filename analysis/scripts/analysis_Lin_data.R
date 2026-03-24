library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)
library(htmlwidgets)
# ============================================================
# 读取数据
# ============================================================
DATA_DIR <- "H:/SPHARM_analysis/analysis/data/derived_data"

df <- read_csv(file.path(DATA_DIR, "SPHARM_direction_lin2024.csv"),
               show_col_types = FALSE)

cat("标本数量：", nrow(df), "\n")
cat("类型分布：\n")
print(df %>% count(Typology))

power_cols <- paste0("power_l", 1:20)
# ============================================================
# 图1：各类型谱熵比较（箱线图）
# ============================================================
p1 <- ggplot(df, aes(x = reorder(Typology, spectral_entropy),
                     y = spectral_entropy,
                     fill = Typology)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  coord_flip() +
  labs(
    title  = "Spectral entropy by typology",
    x      = NULL,
    y      = "Spectral entropy (H)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# ============================================================
# 图2：各类型平均功率谱曲线
# ============================================================
power_cols <- paste0("power_l", 1:20)

df_power <- df %>%
  select(ID, Typology, all_of(power_cols)) %>%
  pivot_longer(
    cols      = all_of(power_cols),
    names_to  = "degree",
    values_to = "power"
  ) %>%
  mutate(degree = as.integer(str_extract(degree, "[0-9]+")))

df_power_mean <- df_power %>%
  group_by(Typology, degree) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    se_power   = sd(power, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p2 <- ggplot(df_power_mean,
             aes(x = degree, y = mean_power,
                 color = Typology, fill = Typology)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(
    aes(ymin = mean_power - se_power,
        ymax = mean_power + se_power),
    alpha = 0.15, color = NA
  ) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(
    title = "Mean power spectrum by typology",
    x     = "Degree (l)",
    y     = "Normalised power",
    color = "Typology",
    fill  = "Typology"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ============================================================
# 图3：功率谱热图（每个标本一行）
# ============================================================
p3 <- ggplot(df_power,
             aes(x = degree, y = ID, fill = power)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "Power") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  facet_grid(Typology ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Power spectrum heatmap",
    x     = "Degree (l)",
    y     = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 7),
    strip.text.y    = element_text(size = 8, angle = 0),
    legend.position = "right"
  )

# ============================================================
# 新增：UMAP 降维分析
# ============================================================
# 1. 准备数据：提取数值矩阵并进行简单预处理
umap_input <- df %>% 
  select(all_of(power_cols)) %>% 
  as.matrix()

# 2. 运行 UMAP (设置随机种子保证结果可重复)
set.seed(42) 
umap_config <- umap.defaults
umap_results <- umap(umap_input, config = umap_config)

# 3. 将结果合并回原始数据框
df_umap <- df %>%
  mutate(
    UMAP1 = umap_results$layout[, 1],
    UMAP2 = umap_results$layout[, 2]
  )

# 4. 绘制 UMAP 散点图
df_hull <- df_umap %>%
  group_by(Typology) %>%
  slice(chull(UMAP1, UMAP2)) %>%
  ungroup()

p4 <- 
  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = transform(df_umap, Typology = NULL), 
             color = "grey90", size = 1) +
  geom_polygon(data = df_hull, aes(fill = Typology), alpha = 0.3) +
  geom_point(aes(color = Typology), size = 1.5) +
  facet_wrap(~Typology, ncol = 4) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.background = element_blank())

p1
p2
p3
p4
# ============================================================
# 导出
# ============================================================
ggsave(file.path(DATA_DIR, "lin2024_spectral_entropy.png"),
       plot = p1, width = 8, height = 6, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_power_spectrum.png"),
       plot = p2, width = 10, height = 6, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_power_heatmap.png"),
       plot = p3, width = 10, height = 12, dpi = 300)

ggsave(file.path(DATA_DIR, "lin2024_UMAP.png"),
       plot = p4, width = 8, height = 8, dpi = 300)

























library(plotly)

# ============================================================
# 1. 读取数据
# ============================================================

DATA_DIR <- "H:/SPHARM_analysis/analysis/data/raw_data/Lin_2024_scar_data"

df <- read_csv(
  file.path(DATA_DIR, "Scar_Vectors_Lin2024.csv"),
  show_col_types = FALSE
)


# ============================================================
# 2. 构建向量 + 单位向量
# ============================================================
df_vec <- df %>%
  mutate(
    dx = X2 - X1,
    dy = Y2 - Y1,
    dz = Z2 - Z1,
    norm = sqrt(dx^2 + dy^2 + dz^2),
    
    # 防止除0
    ux = ifelse(norm == 0, 0, dx / norm),
    uy = ifelse(norm == 0, 0, dy / norm),
    uz = ifelse(norm == 0, 0, dz / norm)
  )

# ============================================================
# 3. 获取标本列表
# ============================================================
specimens <- unique(df_vec$ Specimen_ID)

# ============================================================
# 4. 初始化 plotly
# ============================================================
fig <- plot_ly()

# ============================================================
# 5. 循环添加每个标本
# ============================================================
scale_factor <- 20   # 控制线长度
cone_size <- 10      # 控制箭头大小

for (i in seq_along(specimens)) {
  
  sub <- df_vec %>% filter(Specimen_ID == specimens[i])
  
  fig <- fig %>%
    
    # ---- 点（起点）
    add_trace(
      data = sub,
      x = ~X1, y = ~Y1, z = ~Z1,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 2),
      visible = (i == 1),
      showlegend = FALSE
    ) %>%
    
    # ---- 向量线
    add_trace(
      data = sub,
      x = c(rbind(sub$X1, sub$X1 + sub$ux * scale_factor, NA)),
      y = c(rbind(sub$Y1, sub$Y1 + sub$uy * scale_factor, NA)),
      z = c(rbind(sub$Z1, sub$Z1 + sub$uz * scale_factor, NA)),
      type = "scatter3d",
      mode = "lines",
      line = list(width = 3),
      visible = (i == 1),
      showlegend = FALSE
    ) %>%
    
    # ---- 箭头（cone）
    add_trace(
      data = sub,
      type = "cone",
      
      # 🔥 放在向量终点
      x = ~ (X1 + ux * scale_factor),
      y = ~ (Y1 + uy * scale_factor),
      z = ~ (Z1 + uz * scale_factor),
      
      # 🔥 方向
      u = ~ux,
      v = ~uy,
      w = ~uz,
      
      sizemode = "absolute",
      sizeref = 0.1,        # 👈 小很多！！
      
      anchor = "tip",     # 🔥 关键：箭头对齐到尖端
      
      showscale = FALSE,
      visible = (i == 1))
}

# ============================================================
# 6. dropdown 控制可见性
# ============================================================
n_traces_per_specimen <- 3

visibility_matrix <- lapply(seq_along(specimens), function(i) {
  v <- rep(FALSE, length(specimens) * n_traces_per_specimen)
  idx_start <- (i - 1) * n_traces_per_specimen + 1
  v[idx_start:(idx_start + 2)] <- TRUE
  v
})

fig <- fig %>%
  layout(
    updatemenus = list(
      list(
        buttons = lapply(seq_along(specimens), function(i) {
          list(
            method = "update",
            args = list(list(visible = visibility_matrix[[i]])),
            label = specimens[i]
          )
        }),
        direction = "down",
        x = 0.05,
        y = 0.95
      )
    ),
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z"),
      aspectmode = "data"
    )
  )

# ============================================================
# 7. 导出 HTML
# ============================================================
saveWidget(fig, "scar_vectors_3D.html", selfcontained = TRUE)







