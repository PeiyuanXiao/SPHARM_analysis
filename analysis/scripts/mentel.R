library(tidyverse)
library(vegan)
library(linkET)
library(ggplot2)
library(patchwork)

# ============================================================
# 参数配置
# ============================================================
DATA_DIR   <- "H:/SDG_Lithic_Analysis/analysis/data/drived_data"
LMAX_MORPH <- 20
LMAX_SCAR  <- 20
N_PERM     <- 9999

# ============================================================
# Step 1：读取数据
# ============================================================
df_morph <- read_csv(file.path(DATA_DIR, "SPHARM_results.csv"),
                     show_col_types = FALSE)
df_scar  <- read_csv(file.path(DATA_DIR, "spharm_direction.csv"),
                     show_col_types = FALSE)
df_sphy  <- read_csv(file.path(DATA_DIR, "sphericity_iso.csv"),
                     show_col_types = FALSE)
df_curv  <- read_csv(file.path(DATA_DIR, "curvature.csv"),
                     show_col_types = FALSE)

# 统一顺序
df_scar <- df_scar %>%
  mutate(ID = as.character(ID)) %>%
  slice(match(as.character(df_morph$specimen_id), ID))
df_sphy <- df_sphy %>%
  mutate(ID = as.character(ID)) %>%
  slice(match(as.character(df_morph$specimen_id), ID))
df_curv <- df_curv %>%
  mutate(filename = as.character(filename)) %>%
  slice(match(as.character(df_morph$specimen_id), filename))

cat("ID匹配检查（形态×方向）：",
    all(as.character(df_morph$specimen_id) == df_scar$ID), "\n")
cat("ID匹配检查（形态×球度）：",
    all(as.character(df_morph$specimen_id) == df_sphy$ID), "\n")
cat("ID匹配检查（形态×曲率）：",
    all(as.character(df_morph$specimen_id) == df_curv$filename), "\n\n")

# ============================================================
# Step 2：构建功率谱矩阵
# ============================================================
morph_power <- df_morph %>%
  select(paste0("power_degree_", 1:LMAX_MORPH)) %>%
  rename_with(~ paste0("M", 1:LMAX_MORPH)) %>%
  as.data.frame()

scar_power <- df_scar %>%
  select(paste0("power_l", 1:LMAX_SCAR)) %>%
  rename_with(~ paste0("S", 1:LMAX_SCAR)) %>%
  as.data.frame()

rownames(morph_power) <- as.character(df_morph$specimen_id)
rownames(scar_power)  <- as.character(df_morph$specimen_id)

# 过滤标准差为零的列
sd_morph <- sapply(morph_power, sd, na.rm = TRUE)
sd_scar  <- sapply(scar_power,  sd, na.rm = TRUE)

morph_power_clean <- morph_power[, sd_morph > 0]
scar_power_clean  <- scar_power[,  sd_scar  > 0]

cat("形态保留阶数：", ncol(morph_power_clean), "\n")
cat("片疤保留阶数：", ncol(scar_power_clean),  "\n\n")

# 合并为完整功率谱矩阵（右侧热图）
spec_df_full <- bind_cols(morph_power_clean, scar_power_clean) %>%
  as.data.frame()
rownames(spec_df_full) <- as.character(df_morph$specimen_id)

# ============================================================
# Step 3：构造距离矩阵
# ============================================================
normalize_spectra <- function(X) {
  total <- rowSums(X)
  X / ifelse(total > 0, total, 1)
}

cosine_dist <- function(X) {
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  as.dist(1 - sim)
}

X_morph_norm <- morph_power_clean %>% as.matrix() %>% normalize_spectra()
X_scar_norm  <- scar_power_clean  %>% as.matrix() %>% normalize_spectra()

D_morph <- cosine_dist(X_morph_norm)
D_scar  <- cosine_dist(X_scar_norm)

# ============================================================
# Step 4：全局Mantel检验
# ============================================================
cat("==== 全局Mantel检验（形态 × 片疤方向）====\n")
mantel_global <- mantel(D_morph, D_scar,
                        method       = "pearson",
                        permutations = N_PERM)
print(mantel_global)

# ============================================================
# Step 5：逐变量Mantel检验
# 每个功率谱阶 vs D_morph 和 D_scar
# ============================================================
mantel_rows_full <- map_dfr(colnames(spec_df_full), function(var) {
  x <- spec_df_full[[var]]
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  d_x <- dist(scale(x))
  
  res_m <- mantel(d_x, D_morph,
                  method = "pearson", permutations = N_PERM)
  res_s <- mantel(d_x, D_scar,
                  method = "pearson", permutations = N_PERM)
  
  bind_rows(
    tibble(spec = "Morphology",     env = var,
           r = res_m$statistic, p = res_m$signif),
    tibble(spec = "Scar Direction", env = var,
           r = res_s$statistic, p = res_s$signif)
  )
})

# ============================================================
# Step 6：绘图
# ============================================================
p_mantel_full <- qcorrplot(
  correlate(spec_df_full),
  type = "upper",
  diag = FALSE
) +
  geom_square() +
  geom_couple(
    aes(colour = p, size = abs(r)),
    data      = mantel_rows_full,
    curvature = 0.05
  ) +
  scale_fill_viridis_c(
    option = "D",        
    limits = c(-1, 1),
    name   = "Pearson's R"
  ) +
  scale_color_viridis_c(
    option    = "plasma", 
    direction = -1,       
    limits    = c(0, 1),
    name      = "Mantel's p"
  ) +
  scale_size_continuous(
    range = c(0.3, 2.0),
    name  = "Mantel's |r|"
  ) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    axis.text       = element_text(size = 6),
    plot.title      = element_text(face = "bold", size = 12),
    plot.margin      = margin(t = 10, r = 10, b = 30, l = 80)
  )

ggsave(
  file.path(DATA_DIR, "mantel_linkET_full.png"),
  plot   = p_mantel_full,
  width  = 18,
  height = 16,
  dpi    = 300
)
